! (C) Copyright 2017 UCAR
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 

!> Fortran module to handle radiance observations

module ufo_radiance_mod
  
  use iso_c_binding
  use config_mod
  use duration_mod
  use ufo_obs_data
  use ufo_obs_vectors
  use ufo_vars_mod
  use ufo_locs_mod
  use ufo_geovals_mod
  use kinds
  
  use crtm_module
  implicit none
  private
  
  ! ------------------------------------------------------------------------------
  
  !> Fortran derived type for stream function observations for the QG model
  type :: ufo_obsoper
     integer :: nothing_here_yet
  end type ufo_obsoper
  
#define LISTED_TYPE ufo_obsoper
  
  !> Linked list interface - defines registry_t type
#include "linkedList_i.f"
  
  !> Global registry
  type(registry_t) :: ufo_radiance_registry
  
  ! ------------------------------------------------------------------------------
contains
  ! ------------------------------------------------------------------------------
  !> Linked list implementation
#include "linkedList_c.f"
  
  ! ------------------------------------------------------------------------------
  
  subroutine c_ufo_radiance_setup(c_key_self, c_conf) bind(c,name='ufo_radiance_setup_f90')
    implicit none
    integer(c_int), intent(inout) :: c_key_self
    type(c_ptr), intent(in)    :: c_conf
    
    type(ufo_obsoper), pointer :: self
    
    call ufo_radiance_registry%init()
    call ufo_radiance_registry%add(c_key_self)
    call ufo_radiance_registry%get(c_key_self, self)
    
  end subroutine c_ufo_radiance_setup
  
  ! ------------------------------------------------------------------------------
  
  subroutine c_ufo_radiance_delete(c_key_self) bind(c,name='ufo_radiance_delete_f90')
    implicit none
    integer(c_int), intent(inout) :: c_key_self
    
    type(ufo_obsoper), pointer :: self
    
    call ufo_radiance_registry%get(c_key_self, self)
    call ufo_radiance_registry%remove(c_key_self)
    
  end subroutine c_ufo_radiance_delete
  
  ! ------------------------------------------------------------------------------
  subroutine ufo_radiance_noeqv(c_key_geovals, c_key_hofx, c_bias)
    implicit none
    integer(c_int), intent(in) :: c_key_geovals
    integer(c_int), intent(in) :: c_key_hofx
    real(c_double), intent(in) :: c_bias
    type(obs_vector), pointer :: hofx
    call ufo_obs_vect_registry%get(c_key_hofx,hofx)
    hofx%values(:) = 1.0
  end subroutine ufo_radiance_noeqv
  
  subroutine ufo_radiance_eqv(c_key_geovals, c_key_hofx, c_bias) bind(c,name='ufo_radiance_eqv_f90')
    implicit none
    integer(c_int), intent(in) :: c_key_geovals
    integer(c_int), intent(in) :: c_key_hofx
    real(c_double), intent(in) :: c_bias
    type(ufo_geovals), pointer  :: geovals
    type(obs_vector), pointer :: hofx

    !*************************************************************************************
    !******* Begin CRTM block ************************************************************
    !*************************************************************************************

    ! --------------------------
    ! Some non-CRTM-y Parameters
    ! --------------------------
    CHARACTER(*), PARAMETER :: PROGRAM_NAME   = 'ufo_radiance_mod.F90'
    
    
    ! ============================================================================
    ! STEP 2. **** SET UP SOME PARAMETERS FOR THE CRTM RUN ****
    !

    ! Directory location of coefficients
    !** temporary local path for storing coefficient files (2 files per sensor), also several non-sensor specific binary files needed for other things
    !** NOTE: for some strange reason, this compiled as little endian, even though BIG_ENDIAN was specified on the compiler flags --BTJ
    CHARACTER(*), PARAMETER :: ENDIAN_TYPE='little_endian'
    CHARACTER(*), PARAMETER :: COEFFICIENT_PATH='Data/'

    ! Profile dimensions
    !** UFO to provide N_LAYERS, N_ABSORBERS, N_CLOUDS, N_AEROSOLS
    INTEGER, PARAMETER :: N_PROFILES  = 1  !** required because of the rank of the atm and sfc structures
    INTEGER, PARAMETER :: N_LAYERS    = 64 !** UFO  !** need a way to populate this... 
    INTEGER, PARAMETER :: N_ABSORBERS = 2  !** UFO
    INTEGER, PARAMETER :: N_CLOUDS    = 0  !** UFO
    INTEGER, PARAMETER :: N_AEROSOLS  = 0  !** UFO
    
    ! Sensor information
    INTEGER     , PARAMETER :: N_SENSORS = 1  !** each call to CRTM will be for a single sensor type (zenith/scan angle will be different)
    !  CHARACTER(*), PARAMETER :: SENSOR_ID(N_SENSORS) = (/'cris399_npp','atms_npp   '/)  !** example of how to list multiple sensors
    CHARACTER(*), PARAMETER :: SENSOR_ID(N_SENSORS) = (/'amsua_n19'/)  !** UFO to provide sensor name
    
    ! Some pretend geometry angles. The scan angle is based
    ! on the default Re (earth radius) and h (satellite height)
    REAL(fp), PARAMETER :: ZENITH_ANGLE = 30.0_fp   !** UFO to provide (however, I would not be against creating a geometry database...--BTJ)
    REAL(fp), PARAMETER :: SCAN_ANGLE   = 26.37293341421_fp  !** UFO to provide 
    !** NOTE: From CRTM_Parameters.f90, the maximum zenith angle is fixed at:
    !** REAL(fp), PUBLIC, PARAMETER :: MAX_TRANS_ZENITH_ANGLE = 63.6122_fp !corresponding to amass 2.25
    !**   I will try to figure out why this is the maximum. --BTJ
    ! ============================================================================
    
    ! ---------
    ! Local Variables
    ! ---------
    CHARACTER(256) :: message, version
    INTEGER        :: err_stat, alloc_stat
    INTEGER        :: n_channels
    INTEGER        :: l, m, n, nc, i
    real(fp)       :: cf
    
    
    ! ============================================================================
    ! STEP 3. **** DEFINE THE CRTM INTERFACE STRUCTURES ****
    !
    ! 3a. Define the "non-demoninational" arguments
    ! ---------------------------------------------
    TYPE(CRTM_ChannelInfo_type)             :: chinfo(N_SENSORS)
    TYPE(CRTM_Geometry_type)                :: geo(N_PROFILES)
    
    
    ! 3b. Define the FORWARD variables
    ! --------------------------------
    TYPE(CRTM_Atmosphere_type)              :: atm(N_PROFILES)
    TYPE(CRTM_Surface_type)                 :: sfc(N_PROFILES)
    TYPE(CRTM_RTSolution_type), ALLOCATABLE :: rts(:,:)
    
    
    ! 3c. Define the K-MATRIX variables
    ! ---------------------------------
    TYPE(CRTM_Atmosphere_type), ALLOCATABLE :: atm_K(:,:)
    TYPE(CRTM_Surface_type)   , ALLOCATABLE :: sfc_K(:,:)
    TYPE(CRTM_RTSolution_type), ALLOCATABLE :: rts_K(:,:)
    ! ============================================================================

    type(ufo_geoval)     :: geoval
    character(MAXVARLEN) :: varname
    logical              :: lfound
    integer              :: ivar

    ! Program header
    ! --------------

    ! Get pointers to geovals and hofx
    call ufo_geovals_registry%get(c_key_geovals,geovals)
    call ufo_obs_vect_registry%get(c_key_hofx,hofx)

!** geovals index and variable names:
!!$ 1   Temperature
!!$ 2   Water vapor
!!$ 3   Pressure
!!$ 4   Level pressure
!!$ 5   Ozone
!!$ 6   Water_Fraction
!!$ 7   Land_Fraction
!!$ 8   Ice_Fraction
!!$ 9   Snow_Fraction
!!$ 10  Water_Temperature
!!$ 11  Land_Temperature
!!$ 12  Ice_Temperature
!!$ 13  Snow_Temperature
!!$ 14  Vegetation_Fraction
!!$ 15  Land_Type_Index
    
    CALL CRTM_Version( Version )
    CALL Program_Message( PROGRAM_NAME, &
         'Check/example program for the CRTM Forward and K-Matrix functions using '//&
         ENDIAN_TYPE//' coefficient datafiles', &
         'CRTM Version: '//TRIM(Version) )
    
    ! ============================================================================
    ! STEP 4. **** INITIALIZE THE CRTM ****
    !
    ! 4a. Initialise all the sensors at once
    ! --------------------------------------
    !** NOTE: CRTM_Init points to the various binary files needed for CRTM.  See the
    !**       CRTM_Lifecycle.f90 for more details. 
    WRITE( *,'(/5x,"Initializing the CRTM...")' )
    err_stat = CRTM_Init( SENSOR_ID, &
         chinfo, &
         File_Path=COEFFICIENT_PATH, &
         Quiet=.TRUE.)
    IF ( err_stat /= SUCCESS ) THEN
       message = 'Error initializing CRTM'
       CALL Display_Message( PROGRAM_NAME, message, FAILURE )
       STOP
    END IF
    
    ! 4b. Output some channel information
    ! -----------------------------------
    n_channels = SUM(CRTM_ChannelInfo_n_Channels(chinfo))
    WRITE( *,'(/5x,"Processing a total of ",i0," channels...", i0, " layers..")' ) n_channels, N_LAYERS
    DO n = 1, N_SENSORS
       WRITE( *,'(7x,i0," from ",a)' ) &
            CRTM_ChannelInfo_n_Channels(chinfo(n)), TRIM(SENSOR_ID(n))
    END DO
    ! ============================================================================
    ! Begin loop over sensors
    !** UFO: this loop isn't necessary if we're calling CRTM for each sensor -- it's
    !        not clear to me whether it's more efficient to call all sensors at once
    !        or do each one individually.  I'm leaving this capability intact.  
    ! 
    ! ----------------------------------------------------------------------------
    Sensor_Loop: DO n = 1, N_SENSORS
       
       ! ==========================================================================
       ! STEP 5. **** ALLOCATE STRUCTURE ARRAYS ****
       !
       ! 5a. Determine the number of channels
       !     for the current sensor
       ! ------------------------------------
       n_channels = CRTM_ChannelInfo_n_Channels(chinfo(n))
       
       ! 5b. Allocate the ARRAYS
       ! -----------------------
       ALLOCATE( rts( n_channels, N_PROFILES ), &
            atm_K( n_channels, N_PROFILES ), &
            sfc_K( n_channels, N_PROFILES ), &
            rts_K( n_channels, N_PROFILES ), &
            STAT = alloc_stat )
       IF ( alloc_stat /= 0 ) THEN
          message = 'Error allocating structure arrays'
          CALL Display_Message( PROGRAM_NAME, message, FAILURE )
          STOP
       END IF
       
       ! 5c. Allocate the STRUCTURE INTERNALS
       !     NOTE: Only the Atmosphere structures
       !           are allocated in this example
       ! ----------------------------------------
       ! The input FORWARD structure
       CALL CRTM_Atmosphere_Create( atm, N_LAYERS, N_ABSORBERS, N_CLOUDS, N_AEROSOLS )
       IF ( ANY(.NOT. CRTM_Atmosphere_Associated(atm)) ) THEN
          message = 'Error allocating CRTM Forward Atmosphere structure'
          CALL Display_Message( PROGRAM_NAME, message, FAILURE )
          STOP
       END IF
       
       ! The output K-MATRIX structure
       CALL CRTM_Atmosphere_Create( atm_K, N_LAYERS, N_ABSORBERS, N_CLOUDS, N_AEROSOLS )
       IF ( ANY(.NOT. CRTM_Atmosphere_Associated(atm_K)) ) THEN
          message = 'Error allocating CRTM K-matrix Atmosphere structure'
          CALL Display_Message( PROGRAM_NAME, message, FAILURE )
          STOP
       END IF
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

       CALL Load_Atm_Data()   !** NOTE: could be moved out of sensor loop
       
       !** NOTE:  need to add in cloud and aerosol data to read routine
       
       CALL Load_Sfc_Data()   !** NOTE: could be moved out of sensor loop
       
       ! 6b. Geometry input
       ! ------------------
       ! All profiles are given the same value
       !  The Sensor_Scan_Angle is optional.
       ! ** UFO NOTE: sensor geometry information will need to be provided by calling
       !          routines -- we can't use hardcoded values.  
       CALL CRTM_Geometry_SetValue( geo, &
            Sensor_Zenith_Angle = ZENITH_ANGLE, &
            Sensor_Scan_Angle   = SCAN_ANGLE )
       ! ==========================================================================
       
       ! ==========================================================================
       ! STEP 7. **** INITIALIZE THE K-MATRIX ARGUMENTS ****
       !
       ! 7a. Zero the K-matrix OUTPUT structures
       ! ---------------------------------------
       !** UFO: these structures will be used in the adjoint, so will need to be
       !**      passed back out. 
       CALL CRTM_Atmosphere_Zero( atm_K )
       CALL CRTM_Surface_Zero( sfc_K )
       
       ! 7b. Inintialize the K-matrix INPUT so
       !     that the results are dTb/dx
       ! -------------------------------------
       rts_K%Radiance               = ZERO
       rts_K%Brightness_Temperature = ONE
       ! ==========================================================================

       ! ==========================================================================
       ! STEP 8. **** CALL THE CRTM FUNCTIONS FOR THE CURRENT SENSOR ****
       !
       call CRTM_Atmosphere_Inspect(atm)
       call CRTM_Surface_Inspect(sfc(1))

       WRITE( *, '( /5x, "Calling the CRTM functions for ",a,"..." )' ) TRIM(SENSOR_ID(n))
       
       ! 8a. The forward model call for each sensor
       ! -----------------------------------------------
       err_stat = CRTM_Forward( atm, &  ! Input
            sfc                    , &  ! Input
            geo                    , &  ! Input
            chinfo(n:n)            , &  ! Input
            rts          )              ! Output
       IF ( err_stat /= SUCCESS ) THEN
          message = 'Error calling CRTM Forward Model for '//TRIM(SENSOR_ID(n))
          CALL Display_Message( PROGRAM_NAME, message, FAILURE )
          STOP
       END IF
       
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
       IF ( err_stat /= SUCCESS ) THEN
          message = 'Error calling CRTM K-Matrix Model for '//TRIM(SENSOR_ID(n))
          CALL Display_Message( PROGRAM_NAME, message, FAILURE )
          STOP
       END IF
       ! ==========================================================================
       
       ! ============================================================================
       ! 8c. **** OUTPUT THE RESULTS TO SCREEN (optional) ****
       !
       ! User should read the user guide or the source code of the routine
       ! CRTM_RTSolution_Inspect in the file CRTM_RTSolution_Define.f90 to
       ! select the needed variables for outputs.  These variables are contained
       ! in the structure RTSolution.
       DO m = 1, N_PROFILES
          WRITE( *,'(//7x,"Profile ",i0," output for ",a )') m, TRIM(Sensor_Id(n))
          DO l = 1, n_Channels
!             WRITE( *, '(/5x,"Channel ",i0," results")') chinfo(n)%Sensor_Channel(l)
             !CALL CRTM_RTSolution_Inspect(rts(l,m))
             print '(A,I4,A2,F12.3)', '[Ch] TB: [', chinfo(n)%Sensor_Channel(l), '] ', rts(l,m)%Brightness_Temperature
          END DO
       END DO

       ! output to hofx structure   
       i = 1
       do m = 1, N_PROFILES
         do l = 1, n_Channels
           hofx%values(i) = rts(l,m)%Brightness_Temperature !AS: I'm guessing here.
           i = i + 1
         enddo
       enddo
       
       ! ==========================================================================
       ! STEP 9. **** CLEAN UP FOR NEXT SENSOR ****
       !
       ! 9a. Deallocate the structures
       ! -----------------------------
       CALL CRTM_Atmosphere_Destroy(atm_K)
       CALL CRTM_Atmosphere_Destroy(atm)
       CALL CRTM_RTSolution_Destroy(rts_K)
       CALL CRTM_RTSolution_Destroy(rts)
       CALL CRTM_Surface_Destroy(sfc)
       CALL CRTM_Surface_Destroy(sfc_K)
       
       !** NOTE: Not 100% clear if any of the RTS structures need to be destroyed. 
       
       ! 9b. Deallocate the arrays !** NOTE: this is required
       ! -------------------------
       DEALLOCATE(rts, rts_K, sfc_K, atm_K, STAT = alloc_stat)
       IF ( alloc_stat /= 0 ) THEN
          message = 'Error allocating structure arrays'
          CALL Display_Message( PROGRAM_NAME, message, FAILURE )
          STOP
       END IF
       ! ==========================================================================
       
    END DO Sensor_Loop
    
    ! ==========================================================================
    ! 10. **** DESTROY THE CRTM ****
    !
    WRITE( *, '( /5x, "Destroying the CRTM..." )' )
    err_stat = CRTM_Destroy( chinfo )
    IF ( err_stat /= SUCCESS ) THEN
       message = 'Error destroying CRTM'
       CALL Display_Message( PROGRAM_NAME, message, FAILURE )
       STOP
    END IF
    ! ==========================================================================
    
  CONTAINS
    
    ! ==========================================================================
    !                Below are some internal procedures that load the
    !                necessary input structures with some pretend data
    ! ==========================================================================
    
    !
    ! Internal subprogam to load some test profile data
    !
    SUBROUTINE Load_Atm_Data()
      ! Local variables
      INTEGER :: nc, NL
      INTEGER :: k1, k2
      
      ! 4a.1 Profile #1
      ! ---------------
      ! ...Profile and absorber definitions (fake/placeholder()

!!$ 1   Temperature
!!$ 2   Water vapor
!!$ 3   Pressure
!!$ 4   Level pressure
!!$ 5   Ozone
!!$ 6   Water_Fraction
!!$ 7   Land_Fraction
!!$ 8   Ice_Fraction
!!$ 9   Snow_Fraction
!!$ 10  Water_Temperature
!!$ 11  Land_Temperature
!!$ 12  Ice_Temperature
!!$ 13  Snow_Temperature
!!$ 14  Vegetation_Fraction
!!$ 15  Land_Type_Index

!!$      !       varname = geovals%variables%fldnames(1)
       lfound = ufo_geovals_get_var(geovals,'Temperature             ', geoval)
       atm(1)%Temperature(N_LAYERS:1:-1) = geoval%vals(N_LAYERS:1:-1,1) !** 1 == iobs, hardcoding for testing
       print *, 'Temperature:', atm(1)%Temperature(1:2), geoval%vals(1:2,1)
       lfound = ufo_geovals_get_var(geovals,'Pressure                ', geoval)
       atm(1)%Pressure(N_LAYERS:1:-1) = geoval%vals(N_LAYERS:1:-1,1) !** 1 == iobs, hardcoding for testing
       print *, 'Pressure:', atm(1)%Pressure(1:2), geoval%vals(1:2,1)
       lfound = ufo_geovals_get_var(geovals,'Level pressure          ', geoval)
       atm(1)%Level_Pressure(0:N_LAYERS) = geoval%vals(N_LAYERS+1:1:-1,1) !** 1 == iobs, hardcoding for testing
       print *, 'level_pressure:', atm(1)%Level_Pressure(0:1), geoval%vals(1:2,1)
       atm(1)%Climatology         = US_STANDARD_ATMOSPHERE
       atm(1)%Absorber_Id(1:1)    = (/ H2O_ID /)
       atm(1)%Absorber_Units(1:1) = (/ MASS_MIXING_RATIO_UNITS /)
       lfound = ufo_geovals_get_var(geovals,'Water vapor             ', geoval)
       atm(1)%Absorber(:,1)       = geoval%vals(N_LAYERS:1:-1,1) !** 1 == iobs, hardcoding for testing
       print *, 'water vapor:', atm(1)%Absorber(1:2,1), geoval%vals(1:2,1)
       atm(1)%Absorber_Id(2:2)    = (/ O3_ID /)
       atm(1)%Absorber_Units(2:2) = (/ VOLUME_MIXING_RATIO_UNITS /)
       lfound = ufo_geovals_get_var(geovals,'Ozone                   ', geoval)
       atm(1)%Absorber(:,2)       = geoval%vals(N_LAYERS:1:-1,1) !** 1 == iobs, hardcoding for testing
       print *, 'Ozone:', atm(1)%Absorber(1:2,2), geoval%vals(1:2,1)


      
!!$      atm(1)%Climatology       = US_STANDARD_ATMOSPHERE
!!$      atm(1)%Absorber_Id(1:1)       = (/ H2O_ID /)
!!$      atm(1)%Absorber_Units(1:1)    = (/ MASS_MIXING_RATIO_UNITS /)
!!$      atm(1)%Absorber(:,1)     = 1.0E+00_fp  !** broadcast fill with some value

!!$      atm(1)%Absorber_Id(2:2)       = (/ O3_ID /)
!!$      atm(1)%Absorber_Units(2:2)    = (/ VOLUME_MIXING_RATIO_UNITS /)
!!$      atm(1)%Absorber(:,2)     = 1.0E-03_fp !** broadcast fill with some fake value
!!$
!!$      ! ...make fake Profile data -- this is where UFO/OOPS will populate the various atm objects needed. 
!!$      do NL = 0,N_LAYERS
!!$         !** generate level pressures according to number of layers (faked)
!!$         atm(1)%Level_Pressure(NL) = 0.714_fp*exp(dble(NL)/(dble(N_LAYERS)/7.34_fp))  !** fake pressure generator
!!$         if (NL >= 1) then
!!$            atm(1)%Pressure(NL) = 0.5_fp * (atm(1)%Level_Pressure(NL) - atm(1)%Level_Pressure(NL-1)) + atm(1)%Level_Pressure(NL-1) !** fake pressure profile
!!$            atm(1)%Temperature(NL) = 40*sin(-4.0_fp*ATAN(1.0_fp)+0.1_fp+dble(NL)/20.0_fp)+267.15_fp  !** fake temperature profile generator
!!$            !          print *, NL, atm(1)%Pressure(NL), atm(1)%Temperature(NL)
!!$         end if
!!$      end do
!!$
!!$      ! Cloud data
!!$      IF ( atm(1)%n_Clouds > 0 ) THEN
!!$         k1 = 75
!!$         k2 = 79
!!$         DO nc = 1, atm(1)%n_Clouds
!!$            atm(1)%Cloud(nc)%Type = SNOW_CLOUD
!!$            atm(1)%Cloud(nc)%Effective_Radius(k1:k2) = 500.0_fp ! microns
!!$            atm(1)%Cloud(nc)%Water_Content(k1:k2)    = 10.0_fp  ! kg/m^2
!!$!            atm(1)%Cloud_Fraction = 0.25 !*** when 2.3.0 is released, uncomment this line.  
!!$         END DO
!!$      END IF

      !*** example of loading aerosol data
!!$      Load_Aerosol_Data_1: IF ( atm(1)%n_Aerosols > 0 ) THEN
!!$         atm(1)%Aerosol(1)%Type = DUST_AEROSOL
!!$         atm(1)%Aerosol(1)%Effective_Radius = (/ ... /) ! microns
!!$         atm(1)%Aerosol(1)%Concentration = (/ ... /)
!!$      end IF Load_Aerosol_Data_1
         
      
    END SUBROUTINE Load_Atm_Data
    
    
    !
    ! Internal subprogam to load some test surface data
    !
    SUBROUTINE Load_Sfc_Data()
      implicit none
      integer :: k1
      character(24) :: sfc_types(4)
      real(fp) :: sfrac
      
      ! 4a.0 Surface type definitions for default SfcOptics definitions
      !      For IR and VIS, this is the NPOESS reflectivities.
      ! ---------------------------------------------------------------
      INTEGER, PARAMETER :: TUNDRA_SURFACE_TYPE         = 10  ! NPOESS Land surface type for IR/VIS Land SfcOptics
      INTEGER, PARAMETER :: SCRUB_SURFACE_TYPE          =  7  ! NPOESS Land surface type for IR/VIS Land SfcOptics
      INTEGER, PARAMETER :: COARSE_SOIL_TYPE            =  1  ! Soil type                for MW land SfcOptics
      INTEGER, PARAMETER :: GROUNDCOVER_VEGETATION_TYPE =  7  ! Vegetation type          for MW Land SfcOptics
      INTEGER, PARAMETER :: BARE_SOIL_VEGETATION_TYPE   = 11  ! Vegetation type          for MW Land SfcOptics
      INTEGER, PARAMETER :: SEA_WATER_TYPE              =  1  ! Water type               for all SfcOptics
      INTEGER, PARAMETER :: FRESH_SNOW_TYPE             =  2  ! NPOESS Snow type         for IR/VIS SfcOptics
      INTEGER, PARAMETER :: FRESH_ICE_TYPE              =  1  ! NPOESS Ice type          for IR/VIS SfcOptics
      
      
      
      ! 4a.1 Profile #1  !** UFO: to be provided by UFO
      ! ---------------
      ! ...Land surface characteristics
!!$      sfc(1)%Land_Coverage     = 0.1_fp
!!$      sfc(1)%Land_Type         = TUNDRA_SURFACE_TYPE
!!$      sfc(1)%Land_Temperature  = 272.0_fp
!!$      sfc(1)%Lai               = 0.17_fp
!!$      sfc(1)%Soil_Type         = COARSE_SOIL_TYPE
!!$      sfc(1)%Vegetation_Type   = GROUNDCOVER_VEGETATION_TYPE
!!$      ! ...Water surface characteristics
!!$      sfc(1)%Water_Coverage    = 0.5_fp
!!$      sfc(1)%Water_Type        = SEA_WATER_TYPE
!!$      sfc(1)%Water_Temperature = 275.0_fp
!!$      ! ...Snow coverage characteristics
!!$      sfc(1)%Snow_Coverage    = 0.25_fp
!!$      sfc(1)%Snow_Type        = FRESH_SNOW_TYPE
!!$      sfc(1)%Snow_Temperature = 265.0_fp
!!$      ! ...Ice surface characteristics
!!$      sfc(1)%Ice_Coverage    = 0.15_fp
!!$      sfc(1)%Ice_Type        = FRESH_ICE_TYPE
!!$      sfc(1)%Ice_Temperature = 269.0_fp

!!$ 6   Water_Fraction
!!$ 7   Land_Fraction
!!$ 8   Ice_Fraction
!!$ 9   Snow_Fraction
!!$ 10  Water_Temperature
!!$ 11  Land_Temperature
!!$ 12  Ice_Temperature
!!$ 13  Snow_Temperature
!!$ 14  Vegetation_Fraction
!!$ 15  Land_Type_Index

      !       varname = geovals%variables%fldnames(1)
       !******                               123456789012345678901234'

      !** loop over all surface fractions (need a way to generalize this to avoid changes in indices)
      sfc_types(1:4) = (/'Water_Fraction          ','Land_Fraction           ', 'Ice_Fraction            ', &
           'Snow_Fraction           '/)
      
      lfound = ufo_geovals_get_var(geovals,sfc_types(1), geoval)
      if (geoval%vals(1,1) > 0.0_fp) then
         sfc(1)%Water_Type        = SEA_WATER_TYPE    !** need to check how to determine fresh vs sea water types (salinity???)
         lfound = ufo_geovals_get_var(geovals,'Water_Fraction          ', geoval)
         sfc(1)%Water_Coverage    = geoval%vals(1,1) !** 1 == iobs, hardcoding for testing
         print '(A,2F12.3)', 'Water Coverage:', sfc(1)%Water_Coverage, geoval%vals(1,1)
         lfound = ufo_geovals_get_var(geovals,'Water_Temperature       ', geoval)
         sfc(1)%Water_Temperature = geoval%vals(1,1) !** 1 == iobs, hardcoding for testing
         print '(A,2F12.3)', 'Water Temperature:', sfc(1)%Water_Temperature, geoval%vals(1,1)
      end if
      lfound = ufo_geovals_get_var(geovals,sfc_types(3), geoval)
      if (geoval%vals(1,1) > 0.0_fp) then
         lfound = ufo_geovals_get_var(geovals,'Ice_Fraction            ', geoval)
         sfc(1)%Ice_Coverage    = geoval%vals(1,1) !** 1 == iobs, hardcoding for testing
         print '(A,2F12.3)', 'Ice Coverage:', sfc(1)%Ice_Coverage, geoval%vals(1,1)
         lfound = ufo_geovals_get_var(geovals,'Ice_Temperature         ', geoval)
         sfc(1)%Ice_Temperature = geoval%vals(1,1) !** 1 == iobs, hardcoding for testing
         print '(A,2F12.3)', 'Ice Temperature:', sfc(1)%Ice_Temperature, geoval%vals(1,1)
      end if
      lfound = ufo_geovals_get_var(geovals,sfc_types(4), geoval)
      if (geoval%vals(1,1) > 0.0_fp) then
         lfound = ufo_geovals_get_var(geovals,'Snow_Fraction           ', geoval)
         sfc(1)%Snow_Coverage    = geoval%vals(1,1) !** 1 == iobs, hardcoding for testing
         print '(A,2F12.3)', 'Snow Coverage:', sfc(1)%Snow_Coverage, geoval%vals(1,1)
         lfound = ufo_geovals_get_var(geovals,'Snow_Temperature        ', geoval)
         sfc(1)%Snow_Temperature = geoval%vals(1,1) !** 1 == iobs, hardcoding for testing
         print '(A,2F12.3)', 'Snow Temperature:', sfc(1)%Snow_Temperature, geoval%vals(1,1)
      end if
      lfound = ufo_geovals_get_var(geovals,sfc_types(2), geoval)
      if (geoval%vals(1,1) > 0.0_fp) then
         lfound = ufo_geovals_get_var(geovals,'Land_Type_Index         ', geoval)
         sfc(1)%Land_Type        = geoval%vals(1,1)    !** is this land_type same as CRTM's land type??
         print '(A,2F12.3)', 'Land Type:', sfc(1)%Land_Type, geoval%vals(1,1)
         lfound = ufo_geovals_get_var(geovals,'Land_Fraction           ', geoval)
         sfc(1)%Land_Coverage    = geoval%vals(1,1) !** 1 == iobs, hardcoding for testing
         print '(A,2F12.3)', 'Land Temperature:', sfc(1)%Land_Coverage, geoval%vals(1,1)
         lfound = ufo_geovals_get_var(geovals,'Land_Temperature        ', geoval)
         sfc(1)%Land_Temperature = geoval%vals(1,1) !** 1 == iobs, hardcoding for testing
         print '(A,2F12.3)', 'Land Temperature:', sfc(1)%Land_Temperature, geoval%vals(1,1)
         lfound = ufo_geovals_get_var(geovals,'Vegetation_Fraction     ', geoval)
         sfc(1)%Lai  = geoval%vals(1,1) !** 1 == iobs, hardcoding for testing
         print '(A,2F12.3)', 'Vegetation Fraction:', sfc(1)%Lai, geoval%vals(1,1)

         !** this wasn't provide by the netcdf file, guessing.  
         sfc(1)%Soil_Type         = COARSE_SOIL_TYPE
         sfc(1)%Vegetation_Type   = GROUNDCOVER_VEGETATION_TYPE
      end if
      
    END SUBROUTINE Load_Sfc_Data
    
  end subroutine ufo_radiance_eqv
  
  ! ------------------------------------------------------------------------------
  
  subroutine c_ufo_radiance_inputs(c_key_self, c_key_vars) bind(c,name='ufo_radiance_inputs_f90')
    implicit none
    integer(c_int), intent(in)    :: c_key_self
    integer(c_int), intent(inout) :: c_key_vars
    
    type(ufo_obsoper), pointer :: self
    type(ufo_vars), pointer :: vars
    
    call ufo_radiance_registry%get(c_key_self, self)
    call ufo_vars_registry%init()
    call ufo_vars_registry%add(c_key_vars)
    call ufo_vars_registry%get(c_key_vars, vars)
    
  end subroutine c_ufo_radiance_inputs
  
  ! ------------------------------------------------------------------------------
  subroutine ufo_radiance_equiv_tl(c_key_geovals, c_key_hofx, c_key_traj, c_bias) &
       & bind(c,name='ufo_radiance_equiv_tl_f90')
    implicit none
    integer(c_int), intent(in) :: c_key_geovals
    integer(c_int), intent(in) :: c_key_hofx
    integer(c_int), intent(in) :: c_key_traj
    real(c_double), intent(in) :: c_bias
  end subroutine ufo_radiance_equiv_tl
  ! ------------------------------------------------------------------------------
  subroutine ufo_radiance_equiv_ad(c_key_gom, c_key_hofx, c_key_traj, c_bias) &
       & bind(c,name='ufo_radiance_equiv_ad_f90')
    implicit none
    integer(c_int), intent(in) :: c_key_gom
    integer(c_int), intent(in) :: c_key_hofx
    integer(c_int), intent(in) :: c_key_traj
    real(c_double), intent(inout) :: c_bias
  end subroutine ufo_radiance_equiv_ad
  ! ------------------------------------------------------------------------------
  
end module ufo_radiance_mod
