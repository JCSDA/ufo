! (C) Copyright 2017 UCAR
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 

!> Fortran module to handle radiance observations

module ufo_radiance_tlad_mod
  
  use ioda_obsdb_mod
  use ioda_obs_vectors
  use ufo_vars_mod
  use ioda_locs_mod
  use ufo_geovals_mod
  use kinds  
  use ufo_basis_tlad_mod, only: ufo_basis_tlad
  
  use crtm_module
  
  implicit none
  
  public :: ufo_radiance_tlad
  private
  integer, parameter :: max_string=800
  
  !> Fortran derived type for radiance trajectory
  type, extends(ufo_basis_tlad) :: ufo_radiance_tlad
     
   contains
     procedure :: delete  => ufo_radiance_tlad_delete
     procedure :: settraj => ufo_radiance_tlad_settraj 
     procedure :: eqv_tl  => ufo_radiance_tlad_eqv_tl
     procedure :: eqv_ad  => ufo_radiance_tlad_eqv_ad
  end type ufo_radiance_tlad

contains

  ! ------------------------------------------------------------------------------
  
  subroutine ufo_radiance_tlad_delete(self)
    implicit none
    class(ufo_radiance_tlad), intent(inout)  :: self
    
    self%ltraj = .false.
    
  end subroutine ufo_radiance_tlad_delete
  
  ! ------------------------------------------------------------------------------
  
  subroutine ufo_radiance_tlad_settraj(self, geovals, obss)
    implicit none
    
    class(ufo_radiance_tlad), intent(inout) :: self
    type(ufo_geovals),        intent(in)    :: geovals
    !  type(obs_vector),         intent(inout) :: hofx
    type(ioda_obsdb), target, intent(in)    :: obss
    
    character(len=*), parameter :: myname_="ufo_radiance_tlad_settraj"
    character(max_string) :: err_msg
    
    type(ioda_obsdb), pointer    :: Radiance => NULL()
    
    !*************************************************************************************
    !******* Begin CRTM block ************************************************************
    !*************************************************************************************
    
    ! --------------------------
    ! Some non-CRTM-y Parameters
    ! --------------------------
    character(*), parameter :: PROGRAM_NAME = 'ufo_radiance_tlad_mod.F90'
    
    ! ============================================================================
    ! STEP 2. **** SET UP SOME parameterS FOR THE CRTM RUN ****
    !
    
    ! Directory location of coefficients
    character(*), parameter :: ENDIAN_TYPE='little_endian'
    character(*), parameter :: COEFFICIENT_PATH='Data/'  
    
!!$    ! Profile dimensions
!!$    !** UFO to provide n_Layers, n_Absorbers, n_Clouds, n_Aerosols
    
    integer, parameter :: n_Absorbers = 3  !** UFO
    integer, parameter :: n_Clouds    = 2  !** UFO
    integer, parameter :: n_Aerosols  = 0  !** UFO
    
    integer :: n_Profiles ! = 806  !** required because of the rank of the atm and sfc structures
    integer :: n_Layers   ! = 71 !64 !** UFO  !** need a way to populate this... 
    
    ! Sensor information
    integer     , parameter :: n_Sensors = 1  !** each call to CRTM will be for a single sensor type (zenith/scan angle will be different)
    !  character(*), parameter :: SENSOR_ID(n_Sensors) = (/'cris399_npp','atms_npp   '/)  !** example of how to list multiple sensors
    character(*), parameter :: SENSOR_ID(n_Sensors) = (/'amsua_n19'/)  !** UFO to provide sensor name
    
    !** these remaining items are still missing from UFO -> CRTM, likely available from locs.
!!$    REAL(fp), parameter :: Latitude  = 46.3369_fp
!!$    REAL(fp), parameter :: Longitude = 354.4514_fp
!!$    REAL(fp), parameter :: Elevation = 161_fp
!!$    REAL(fp), parameter :: Obs_Time  = -1.83777777777778_fp
    
    !** NOTE: From CRTM_Parameters.f90, the maximum zenith angle is fixed at:
    !** REAL(fp), PUBLIC, parameter :: MAX_TRANS_ZENITH_ANGLE = 63.6122_fp !corresponding to amass 2.25
    
    ! ============================================================================
    
    ! ---------
    ! Local Variables
    ! ---------
    character(256) :: message, version
    integer        :: err_stat, alloc_stat
    integer        :: n_Channels
    integer        :: l, m, n, nc, i
    real(fp)       :: cf
    
    ! ============================================================================
    ! STEP 3. **** DEFINE THE CRTM INTERFACE STRUCTURES ****
    !
    ! 3a. Define the "non-demoninational" arguments
    ! ---------------------------------------------
    type(CRTM_ChannelInfo_type)             :: chinfo(n_Sensors)
    type(CRTM_Geometry_type),   allocatable :: geo(:)
    
    
    ! 3b. Define the FORWARD variables
    ! --------------------------------
    type(CRTM_Atmosphere_type), allocatable :: atm(:)
    type(CRTM_Surface_type),    allocatable :: sfc(:)
    type(CRTM_RTSolution_type), allocatable :: rts(:,:)
    
    
    ! 3c. Define the K-MATRIX variables
    ! ---------------------------------
    type(CRTM_Atmosphere_type), allocatable :: atm_K(:,:)
    type(CRTM_Surface_type)   , allocatable :: sfc_K(:,:)
    type(CRTM_RTSolution_type), allocatable :: rts_K(:,:)
    ! ============================================================================
    
    type(ufo_geoval), pointer :: geoval
    character(MAXVARLEN)      :: varname
    logical                   :: lfound
    integer                   :: ivar, nobs, nlocs
    real(fp)                  :: rmse
    real(fp), allocatable     :: diff(:,:)
    
    !** Refer Radiance to obss 
    Radiance => obss
    
    !** get number of profiles and number of layers from geovals input
    n_Profiles = geovals%nobs
    n_Layers   = geovals%geovals(1)%nval
    
    ! Program header
    ! --------------
    
    call CRTM_Version( Version )
    call Program_Message( PROGRAM_NAME, &
         'Check/example program for the CRTM Forward and K-Matrix functions using '//&
         ENDIAN_type//' coefficient datafiles', &
         'CRTM Version: '//TRIM(Version) )
    
    ! ============================================================================
    ! STEP 4. **** INITIALIZE THE CRTM ****
    !
    ! 4a. Initialise all the sensors at once
    ! --------------------------------------
    !** NOTE: CRTM_Init points to the various binary files needed for CRTM.  See the
    !**       CRTM_Lifecycle.f90 for more details. 
    write( *,'(/5x,"Initializing the CRTM...")' )
    err_stat = CRTM_Init( SENSOR_ID, &
         chinfo, &
         File_Path=COEFFICIENT_PATH, &
         Quiet=.TRUE.)
    if ( err_stat /= SUCCESS ) THEN
       message = 'Error initializing CRTM'
       call Display_Message( PROGRAM_NAME, message, FAILURE )
       stop
    end if
    
    ! ============================================================================
    ! Begin loop over sensors
    !** UFO: this loop isn't necessary if we're calling CRTM for each sensor -- it's
    !        not clear to me whether it's more efficient to call all sensors at once
    !        or do each one individually.  I'm leaving this capability intact. BTJ 
    ! 
    ! ----------------------------------------------------------------------------
    Sensor_Loop:do n = 1, n_Sensors
       
       ! ==========================================================================
       ! STEP 5. **** ALLOCATE STRUCTURE ARRAYS ****
       !
       ! 5a. Determine the number of channels
       !     for the current sensor
       ! ------------------------------------
       n_Channels = CRTM_ChannelInfo_n_Channels(chinfo(n))
       
       ! 5b. Allocate the ARRAYS
       ! -----------------------
       allocate( geo( n_Profiles ),          &
            atm( n_Profiles ),               &
            sfc( n_Profiles ),               &
            rts( n_Channels, n_Profiles ),   &
            atm_K( n_Channels, n_Profiles ), &
            sfc_K( n_Channels, n_Profiles ), &
            rts_K( n_Channels, n_Profiles ), &
            STAT = alloc_stat )
       if ( alloc_stat /= 0 ) THEN
          message = 'Error allocating structure arrays'
          call Display_Message( PROGRAM_NAME, message, FAILURE )
          stop
       end if
       
       ! 5c. Allocate the STRUCTURE INTERNALS
       !     NOTE: Only the Atmosphere structures
       !           are allocated in this example
       ! ----------------------------------------
       ! The input FORWARD structure
       call CRTM_Atmosphere_Create( atm, n_Layers, n_Absorbers, n_Clouds, n_Aerosols )
       if ( ANY(.NOT. CRTM_Atmosphere_Associated(atm)) ) THEN
          message = 'Error allocating CRTM Forward Atmosphere structure'
          call Display_Message( PROGRAM_NAME, message, FAILURE )
          stop
       end if
       
       call CRTM_Surface_Create(sfc, n_Channels)
       if ( ANY(.NOT. CRTM_Surface_Associated(sfc)) ) THEN
          message = 'Error allocating CRTM Surface structure'
          call Display_Message( PROGRAM_NAME, message, FAILURE )
          stop
       end if
       
       ! The output K-MATRIX structure
       call CRTM_Atmosphere_Create( atm_K, n_Layers, n_Absorbers, n_Clouds, n_Aerosols )
       if ( ANY(.NOT. CRTM_Atmosphere_Associated(atm_K)) ) THEN
          message = 'Error allocating CRTM K-matrix Atmosphere structure'
          call Display_Message( PROGRAM_NAME, message, FAILURE )
          stop
       end if
       
       call CRTM_Surface_Create(sfc_K, n_Channels)
       if ( ANY(.NOT. CRTM_Surface_Associated(sfc_K)) ) THEN
          message = 'Error allocating CRTM K-matrix Surface structure'
          call Display_Message( PROGRAM_NAME, message, FAILURE )
          stop
       end if
       
       ! ==========================================================================
       
       ! ==========================================================================
       ! STEP 6. **** ASSIGN INPUT DATA ****
       !
       ! 6a. Atmosphere and Surface input
       !     NOTE: that this is the hard part (in my opinion :o). The mechanism by
       !     by which the atmosphere and surface data are loaded in to their
       !     respective structures below was done purely to keep the step-by-step
       !     instructions in this program relatively "clean".
       ! ------------------------------------------------------------------------
       !** UFO NOTE: this is where input data from UFO/OOPS will be loaded
       !**           subroutines not necessary, but helps cleanly separate atmos
       !**           and surface data. 
       
       call Load_Atm_Data()   !** NOTE: could be moved out of sensor loop
       
       !** NOTE:  need to add in aerosol data to read routine
       
       call Load_Sfc_Data()   !** NOTE: could be moved out of sensor loop
       
       ! 6b. Geometry input
       ! ------------------
       ! All profiles are given the same value
       !  The Sensor_Scan_Angle is optional.
       call Load_Geom_Data()
       
       ! ==========================================================================
       
       ! ==========================================================================
       ! STEP 7. **** INITIALIZE THE K-MATRIX ARGUMENTS ****
       !
       ! 7a. Zero the K-matrix OUTPUT structures
       ! ---------------------------------------
       !** UFO: these structures will be used in the adjoint, so will need to be
       !**      passed back out. 
       call CRTM_Atmosphere_Zero( atm_K )
       call CRTM_Surface_Zero( sfc_K )
       
       ! 7b. Inintialize the K-matrix INPUT so
       !     that the results are dTb/dx
       ! -------------------------------------
       rts_K%Radiance               = ZERO
       rts_K%Brightness_Temperature = ONE
       ! ==========================================================================
       
       ! ==========================================================================
       ! STEP 8. **** call THE CRTM FUNCTIONS FOR THE CURRENT SENSOR ****
       !
       call CRTM_Atmosphere_Inspect(atm(12))
       call CRTM_Surface_Inspect(sfc(12))
       call CRTM_Geometry_Inspect(geo(12))
       call CRTM_ChannelInfo_Inspect(chinfo(1))
       
       !       write( *, '( /5x, "Calling the CRTM functions for ",a,"..." )' ) TRIM(SENSOR_ID(n))
       
       ! 8a. The forward model call for each sensor
       ! -----------------------------------------------
       err_stat = CRTM_Forward( atm, &  ! Input
            sfc                    , &  ! Input
            geo                    , &  ! Input
            chinfo(n:n)            , &  ! Input
            rts          )              ! Output
       if ( err_stat /= SUCCESS ) THEN
          message = 'Error calling CRTM Forward Model for '//TRIM(SENSOR_ID(n))
          call Display_Message( PROGRAM_NAME, message, FAILURE )
          stop
       end if
       
       
       ! 8b. The K-matrix model
       ! ----------------------
       err_stat = CRTM_K_Matrix( atm, &  ! FORWARD  Input
            sfc                     , &  ! FORWARD  Input
            rts_K                   , &  ! K-MATRIX Input
            geo                     , &  ! Input
            chinfo(n:n)             , &  ! Input
            atm_K                   , &  ! K-MATRIX Output
            sfc_K                   , &  ! K-MATRIX Output
            rts          )               ! FORWARD  Output
       if ( err_stat /= SUCCESS ) THEN
          message = 'Error calling CRTM K-Matrix Model for '//TRIM(SENSOR_ID(n))
          call Display_Message( PROGRAM_NAME, message, FAILURE )
          stop
       end if
       ! ==========================================================================

       ! ==========================================================================
       ! STEP 9. **** CLEAN UP FOR NEXT SENSOR ****
       !
       ! 9a. Deallocate the structures
       ! -----------------------------
       call CRTM_Geometry_Destroy(geo)
       call CRTM_Atmosphere_Destroy(atm_K)
       call CRTM_Atmosphere_Destroy(atm)
       call CRTM_RTSolution_Destroy(rts_K)
       call CRTM_RTSolution_Destroy(rts)
       call CRTM_Surface_Destroy(sfc)
       call CRTM_Surface_Destroy(sfc_K)
       
       !** NOTE: Not 100% clear if any of the RTS structures need to be destroyed. 
       
       ! 9b. Deallocate the arrays !** NOTE: this is required
       ! -------------------------
       deallocate(geo, atm, sfc, rts, rts_K, sfc_K, atm_K, STAT = alloc_stat)
       if ( alloc_stat /= 0 ) THEN
          message = 'Error deallocating structure arrays'
          call Display_Message( PROGRAM_NAME, message, FAILURE )
          stop
       end if
       ! ==========================================================================
       
    end do Sensor_Loop
    
    ! ==========================================================================
    ! 10. **** DESTROY THE CRTM ****
    !
    write( *, '( /5x, "Destroying the CRTM..." )' )
    err_stat = CRTM_Destroy( chinfo )
    if ( err_stat /= SUCCESS ) THEN
       message = 'Error destroying CRTM'
       call Display_Message( PROGRAM_NAME, message, FAILURE )
       stop
    end if
    ! ==========================================================================
    
    ! ==========================================================================
    !                Below are some internal procedures that load the
    !                necessary input structures with some pretend data
    ! ==========================================================================
    
    !
    ! Internal subprogam to load some test profile data
    !
    subroutine Load_Atm_Data()
      implicit none
      ! Local variables
      integer :: nc, NL
      integer :: k1, k2
      
      ! 4a.1 Profile #1
      ! ---------------
      ! ...Profile and absorber definitions
      
      do k1 = 1,geovals%nvar
         varname = geovals%variables%fldnames(k1)
         print *, k1, varname
      end do
      
      !** populate the atmosphere structures for CRTM (atm(k1), for the k1-th profile)
      do k1 = 1,n_Profiles
         lfound = ufo_geovals_get_var(geovals, var_tv, geoval)
         atm(k1)%Temperature(1:n_Layers) = geoval%vals(:,k1) 
         
         lfound = ufo_geovals_get_var(geovals, var_prs, geoval)
         atm(k1)%Pressure(1:n_Layers) = geoval%vals(:,k1) 
         
         lfound = ufo_geovals_get_var(geovals, var_prsi, geoval)
         atm(k1)%Level_Pressure(0:n_Layers) = geoval%vals(:,k1)
         
         atm(k1)%Climatology            = US_STANDARD_ATMOSPHERE
         
         atm(k1)%Absorber_Id(1:1)       = (/ H2O_ID /)
         atm(k1)%Absorber_Units(1:1)    = (/ MASS_MIXING_RATIO_UNITS /)
         lfound = ufo_geovals_get_var(geovals, var_mixr, geoval)
         atm(k1)%Absorber(1:n_Layers,1) = geoval%vals(:,k1) 
         
         atm(k1)%Absorber_Id(2:2)       = (/ O3_ID /)
         atm(k1)%Absorber_Units(2:2)    = (/ VOLUME_MIXING_RATIO_UNITS /)
         lfound = ufo_geovals_get_var(geovals, var_oz, geoval)
         atm(k1)%Absorber(1:n_Layers,2) = geoval%vals(:,k1) 
         
         atm(k1)%Absorber_Id(3:3)       = (/ CO2_ID /)
         atm(k1)%Absorber_Units(3:3)    = (/ VOLUME_MIXING_RATIO_UNITS /)
         lfound = ufo_geovals_get_var(geovals, var_co2, geoval)
         atm(k1)%Absorber(1:n_Layers,3) = geoval%vals(:,k1)
         
         atm(k1)%Cloud(1)%Type = WATER_CLOUD
         lfound = ufo_geovals_get_var(geovals, var_clw, geoval)
         atm(k1)%Cloud(1)%Water_Content    = geoval%vals(:,k1)
         lfound = ufo_geovals_get_var(geovals, var_clwefr, geoval)
         atm(k1)%Cloud(1)%Effective_Radius = geoval%vals(:,k1)
         
         atm(k1)%Cloud(2)%Type = ICE_CLOUD
         lfound = ufo_geovals_get_var(geovals, var_cli, geoval)
         atm(k1)%Cloud(2)%Water_Content    = geoval%vals(:,k1)
         lfound = ufo_geovals_get_var(geovals, var_cliefr, geoval)
         atm(k1)%Cloud(2)%Effective_Radius = geoval%vals(:,k1)
      end do
      
    end subroutine Load_Atm_Data
    
    
    !
    ! Internal subprogam to load some test surface data
    !
    subroutine Load_Sfc_Data()
      implicit none
      integer  :: k1
      real(fp) :: sfrac
      
      ! 4a.0 Surface type definitions for default SfcOptics definitions
      !      For IR and VIS, this is the NPOESS reflectivities.
      ! ---------------------------------------------------------------
      integer, parameter :: TUNDRA_SURFACE_type         = 10  ! NPOESS Land surface type for IR/VIS Land SfcOptics
      integer, parameter :: SCRUB_SURFACE_type          =  7  ! NPOESS Land surface type for IR/VIS Land SfcOptics
      integer, parameter :: COARSE_SOIL_type            =  1  ! Soil type                for MW land SfcOptics
      integer, parameter :: GROUNDCOVER_VEGETATION_type =  7  ! Vegetation type          for MW Land SfcOptics
      integer, parameter :: BARE_SOIL_VEGETATION_type   = 11  ! Vegetation type          for MW Land SfcOptics
      integer, parameter :: SEA_WATER_type              =  1  ! Water type               for all SfcOptics
      integer, parameter :: FRESH_SNOW_type             =  2  ! NPOESS Snow type         for IR/VIS SfcOptics
      integer, parameter :: FRESH_ICE_type              =  1  ! NPOESS Ice type          for IR/VIS SfcOptics
      
      real(kind_real), allocatable :: Radiance_Tbobs(:,:)
      type(obs_vector)             :: TmpOvec
      integer                      :: ch
      
      ! 4a.1 Surface Characteristics
      ! ---------------
      ! ...Land surface characteristics
      
      allocate(Radiance_Tbobs(n_Channels, n_Profiles))
      call ioda_obsvec_setup(TmpOvec, Radiance%nobs)
      call ioda_obsdb_var_to_ovec(Radiance, TmpOvec, "Observation")
      Radiance_Tbobs = reshape(TmpOvec%values, (/n_Channels, n_Profiles/))
      
      do k1 = 1,n_Profiles
         sfc(k1)%sensordata%sensor_id        = chinfo(1)%sensor_id
         sfc(k1)%sensordata%wmo_sensor_id    = chinfo(1)%wmo_sensor_id
         sfc(k1)%sensordata%wmo_satellite_id = chinfo(1)%wmo_satellite_id
         sfc(k1)%sensordata%sensor_channel   = chinfo(1)%sensor_channel
         
         sfc(k1)%Water_Type         = SEA_WATER_type    !** NOTE: need to check how to determine fresh vs sea water types (salinity???)
         lfound                     = ufo_geovals_get_var(geovals, var_sfc_wspeed, geoval)
         sfc(k1)%Wind_Speed         = geoval%vals(1,k1) 
         lfound                     = ufo_geovals_get_var(geovals, var_sfc_wdir, geoval)
         sfc(k1)%Wind_Direction     = geoval%vals(1,k1) 
         lfound                     = ufo_geovals_get_var(geovals, var_sfc_wfrac, geoval)
         sfc(k1)%Water_Coverage     = geoval%vals(1,k1) 
         lfound                     = ufo_geovals_get_var(geovals, var_sfc_wtmp, geoval)
         sfc(k1)%Water_Temperature  = geoval%vals(1,k1) 
         lfound                     = ufo_geovals_get_var(geovals, var_sfc_ifrac, geoval)
         sfc(k1)%Ice_Coverage       = geoval%vals(1,k1) 
         lfound                     = ufo_geovals_get_var(geovals, var_sfc_itmp, geoval)
         sfc(k1)%Ice_Temperature    = geoval%vals(1,k1) 
         lfound                     = ufo_geovals_get_var(geovals, var_sfc_sfrac, geoval)
         sfc(k1)%Snow_Coverage      = geoval%vals(1,k1) 
         lfound                     = ufo_geovals_get_var(geovals, var_sfc_stmp, geoval)
         sfc(k1)%Snow_Temperature   = geoval%vals(1,k1) 
         lfound                     = ufo_geovals_get_var(geovals, var_sfc_sdepth, geoval)
         sfc(k1)%Snow_Depth         = geoval%vals(1,k1)
         lfound                     = ufo_geovals_get_var(geovals, var_sfc_landtyp, geoval)
         sfc(k1)%Land_Type          = geoval%vals(1,k1)    !** NOTE:  is this Land_Type same as CRTM's land type??
         lfound                     = ufo_geovals_get_var(geovals, var_sfc_lfrac, geoval)
         sfc(k1)%Land_Coverage      = geoval%vals(1,k1) 
         lfound                     = ufo_geovals_get_var(geovals, var_sfc_ltmp, geoval)
         sfc(k1)%Land_Temperature   = geoval%vals(1,k1) 
         lfound                     = ufo_geovals_get_var(geovals, var_sfc_lai, geoval)
         sfc(k1)%Lai                = geoval%vals(1,k1) 
         lfound                     = ufo_geovals_get_var(geovals, var_sfc_vegfrac, geoval)
         sfc(k1)%Vegetation_Fraction = geoval%vals(1,k1) 
         lfound                     = ufo_geovals_get_var(geovals, var_sfc_vegtyp, geoval)
         sfc(k1)%Vegetation_Type    = geoval%vals(1,k1) 
         lfound                     = ufo_geovals_get_var(geovals, var_sfc_soiltyp, geoval)
         sfc(k1)%Soil_Type          = geoval%vals(1,k1) 
         lfound                     = ufo_geovals_get_var(geovals, var_sfc_soilm, geoval)
         sfc(k1)%Soil_Moisture_Content = geoval%vals(1,k1) 
         lfound                     = ufo_geovals_get_var(geovals, var_sfc_soilt, geoval)
         sfc(k1)%Soil_Temperature   = geoval%vals(1,k1) 
         do ch = 1, n_Channels
            sfc(k1)%sensordata%tb(ch) = Radiance_TbObs(ch, k1)  !** required to match GSI simulated TBs over snow and ice surfaces
         end do
      end do
      deallocate(Radiance_Tbobs)
      call ioda_obsvec_delete(TmpOvec)
      
    end subroutine Load_Sfc_Data
    
    !
    ! Internal subprogam to load some test geometry data
    !
    subroutine Load_Geom_Data()
      implicit none
      
      type(obs_vector) :: TmpOvec
      
      call ioda_obsvec_setup(TmpOvec, Radiance%nobs)
      
      call ioda_obsdb_var_to_ovec(Radiance, TmpOvec, "Sat_Zenith_Angle")
      geo(:)%Sensor_Zenith_Angle = TmpOvec%values(::n_Channels)
      call ioda_obsdb_var_to_ovec(Radiance, TmpOvec, "Sol_Zenith_Angle")
      geo(:)%Source_Zenith_Angle = TmpOvec%values(::n_Channels)
      call ioda_obsdb_var_to_ovec(Radiance, TmpOvec, "Sat_Azimuth_Angle")
      geo(:)%Sensor_Azimuth_Angle = TmpOvec%values(::n_Channels)
      call ioda_obsdb_var_to_ovec(Radiance, TmpOvec, "Sol_Azimuth_Angle")
      geo(:)%Source_Azimuth_Angle = TmpOvec%values(::n_Channels)
      call ioda_obsdb_var_to_ovec(Radiance, TmpOvec, "Scan_Position")
      geo(:)%Ifov = TmpOvec%values(::n_Channels) 
      call ioda_obsdb_var_to_ovec(Radiance, TmpOvec, "Scan_Angle")
      geo(:)%Sensor_Scan_Angle = TmpOvec%values(::n_Channels)
      
      call ioda_obsvec_delete(TmpOvec)
      
    end subroutine Load_Geom_Data
    
    self%ltraj = .true.

  end subroutine ufo_radiance_tlad_settraj

  ! ------------------------------------------------------------------------------
  
  subroutine ufo_radiance_tlad_eqv_tl(self, geovals, hofx, obss)
    implicit none
    class(ufo_radiance_tlad), intent(in)     :: self
    type(ufo_geovals),     intent(in)     :: geovals
    type(obs_vector),      intent(inout)  :: hofx
    type(ioda_obsdb),      intent(in)     :: obss
    
    character(len=*), parameter :: myname_="ufo_radiance_tlad_eqv_tl"
    character(max_string) :: err_msg
    
    ! Nothing here yet
    
  end subroutine ufo_radiance_tlad_eqv_tl
  
  ! ------------------------------------------------------------------------------
  
  subroutine ufo_radiance_tlad_eqv_ad(self, geovals, hofx, obss)
    implicit none
    class(ufo_radiance_tlad), intent(in) :: self
    type(ufo_geovals),     intent(inout)    :: geovals
    type(obs_vector),      intent(in)       :: hofx
    type(ioda_obsdb),      intent(in)       :: obss
    
    character(len=*), parameter :: myname_="ufo_radiance_tlad_eqv_ad"
    character(max_string) :: err_msg
    
    ! Nothing here yet
    
  end subroutine ufo_radiance_tlad_eqv_ad
  
! ------------------------------------------------------------------------------
  
end module ufo_radiance_tlad_mod



  
    
    
