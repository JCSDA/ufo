! (C) Copyright 2017 UCAR
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 

!> Fortran module to handle radiance observations

module ufo_radiance_mod
  
  use ufo_obs_data
  use ufo_obs_data_mod
  use ufo_obs_vectors
  use ufo_locs_mod
  use ufo_geovals_mod
  use ufo_geovals_mod_c, only: ufo_geovals_registry
  use ufo_vars_mod
  use kinds
  
  use crtm_module
  implicit none

  public :: ufo_radiance_eqv

  private
  
  
  ! ------------------------------------------------------------------------------
contains
  
  ! ------------------------------------------------------------------------------

  subroutine ufo_radiance_eqv(geovals, obss, hofx) 
    use radDiag_mod, only: RadDiag
    use ufo_obs_data_mod, only: Radiance
    implicit none
    type(ufo_geovals), intent(in)    :: geovals
    type(obs_data),    intent(inout) :: obss
    type(obs_vector),  intent(inout) :: hofx

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
    INTEGER, PARAMETER :: N_PROFILES  = 806  !** required because of the rank of the atm and sfc structures
    INTEGER, PARAMETER :: N_LAYERS    = 71 !64 !** UFO  !** need a way to populate this... 
    INTEGER, PARAMETER :: N_ABSORBERS = 3  !** UFO
    INTEGER, PARAMETER :: N_CLOUDS    = 2  !** UFO
    INTEGER, PARAMETER :: N_AEROSOLS  = 0  !** UFO
    
    ! Sensor information
    INTEGER     , PARAMETER :: N_SENSORS = 1  !** each call to CRTM will be for a single sensor type (zenith/scan angle will be different)
    !  CHARACTER(*), PARAMETER :: SENSOR_ID(N_SENSORS) = (/'cris399_npp','atms_npp   '/)  !** example of how to list multiple sensors
    CHARACTER(*), PARAMETER :: SENSOR_ID(N_SENSORS) = (/'amsua_n19'/)  !** UFO to provide sensor name
    
    ! Some pretend geometry angles. The scan angle is based
    ! on the default Re (earth radius) and h (satellite height)
    REAL(fp), PARAMETER :: ZENITH_ANGLE      = -44.65_fp   !** UFO to provide (however, I would not be against creating a geometry database...--BTJ)
    REAL(fp), PARAMETER :: SCAN_ANGLE        = -35.0_fp  !** UFO to provide
    REAL(fp), PARAMETER :: Latitude          = 46.3369
    REAL(fp), PARAMETER :: Longitude         = 354.4514
    REAL(fp), PARAMETER :: Elevation         = 161
    REAL(fp), PARAMETER :: Obs_Time          = -1.83777777777778
    REAL(fp), PARAMETER :: Scan_Position     = 4
    REAL(fp), PARAMETER :: Sat_Zenith_Angle  = -44.65
    REAL(fp), PARAMETER :: Sat_Azimuth_Angle = 290.23
    REAL(fp), PARAMETER :: Sol_Zenith_Angle  = 118.88
    REAL(fp), PARAMETER :: Sol_Azimuth_Angle = 66.63
                
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

    integer              :: nobs
    integer              :: nlocs
    real(fp)              :: rmse
    real(fp), allocatable :: diff(:,:)

    ! Program header
    ! --------------

    nobs=obss%nobs; nlocs=obss%nlocs
    obss%Obspoint => Radiance

!** geovals index and variable names:
!!$ 1   Temperature
!!$ 2   Water vapor
!!$ 3   Pressure
!!$ 4   Level pressure
!!$ 5   Ozone
!!$ 6   Cloud liquid
!!$ 7   Cloud ice
!!$ 8   Water_Fraction
!!$ 9   Land_Fraction
!!$ 10  Ice_Fraction
!!$ 11  Snow_Fraction
!!$ 12  Water_Temperature
!!$ 13  Land_Temperature
!!$ 14  Ice_Temperature
!!$ 15  Snow_Temperature
!!$ 16  Vegetation_Fraction
!!$ 17  Sfc_Wind_Speed
!!$ 18  Sfc_Wind_Direction
!!$ 19  Lai
!!$ 20  Soil_Moisture
!!$ 21  Soil_Temperature
!!$ 22  Land_Type_Index
!!$ 23  Vegetation_Type
!!$ 24  Soil_Type

    
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
    !WRITE( *,'(/5x,"Processing a total of ",i0," channels...", i0, " layers..")' ) n_channels, N_LAYERS
    DO n = 1, N_SENSORS
       !WRITE( *,'(7x,i0," from ",a)' ) &
       !     CRTM_ChannelInfo_n_Channels(chinfo(n)), TRIM(SENSOR_ID(n))
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
       CALL Load_Geom_Data()

!!$       REAL(fp), PARAMETER :: Latitude          = 46.3369
!!$       REAL(fp), PARAMETER :: Longitude         = 354.4514
!!$       REAL(fp), PARAMETER :: Elevation         = 161
!!$       REAL(fp), PARAMETER :: Obs_Time          = -1.83777777777778
!!$       REAL(fp), PARAMETER :: Scan_Position     = 4  !** at 3.333 degrees per scan position, starting at 48.333, so this is -48.333+4*3.333 = 
!!$       REAL(fp), PARAMETER :: Sat_Zenith_Angle  = -44.65
!!$       REAL(fp), PARAMETER :: Sat_Azimuth_Angle = 290.23
!!$       REAL(fp), PARAMETER :: Sol_Zenith_Angle  = 118.88
!!$       REAL(fp), PARAMETER :: Sol_Azimuth_Angle = 66.63
       
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
       call CRTM_Atmosphere_Inspect(atm(12))
       call CRTM_Surface_Inspect(sfc(12))
       call CRTM_Geometry_Inspect(geo(12))
       call CRTM_ChannelInfo_Inspect(chinfo(1))

!       WRITE( *, '( /5x, "Calling the CRTM functions for ",a,"..." )' ) TRIM(SENSOR_ID(n))
       
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

       allocate(diff(n_channels,n_profiles))
       rmse = 0
       DO m = 1, N_PROFILES
          DO l = 1, n_Channels
             diff(l,m) = rts(l,m)%Brightness_Temperature - (Radiance%datachan(m,l)%tbobs - Radiance%datachan(m,l)%omgnbc)
!             print *, rts(l,m)%Brightness_Temperature, Radiance%datachan(m,l)%tbobs - Radiance%datachan(m,l)%omgnbc
             rmse = rmse + (Radiance%datachan(m,l)%tbobs - Radiance%datachan(m,l)%omgnbc) * (Radiance%datachan(m,l)%tbobs - Radiance%datachan(m,l)%omgnbc)
          END DO
          WRITE( *,'(//7x,"Profile ",i0," output for ",a, " difference:",f12.6 )') m, TRIM(Sensor_Id(n)), maxval(abs(diff(:,m)))
       END DO
       print *, 'Max difference: ', maxval(abs(diff))
       deallocate(diff)

       rmse = sqrt(rmse / (n_profiles * n_channels))
       print *, 'rmse: ', rmse

       ! output to hofx structure   
       hofx%values(:) = 0.0
       i = 1
       do m = 1, N_PROFILES
         do l = 1, n_Channels
           hofx%values(i) = rts(l,m)%Brightness_Temperature
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
      implicit none
      ! Local variables
      INTEGER :: nc, NL
      INTEGER :: k1, k2
      
      ! 4a.1 Profile #1
      ! ---------------
      ! ...Profile and absorber definitions (fake/placeholder()

!!$ 1     Temperature
!!$ 2     Water vapor
!!$ 3     Pressure
!!$ 4     Level pressure
!!$ 5     Ozone
!!$ 6     Cloud liquid
!!$ 7     Cloud ice
!!$ 8     Water_Fraction
!!$ 9     Land_Fraction
!!$ 10    Ice_Fraction
!!$ 11    Snow_Fraction
!!$ 12    Water_Temperature
!!$ 13    Land_Temperature
!!$ 14    Ice_Temperature
!!$ 15    Snow_Temperature
!!$ 16    Vegetation_Fraction
!!$ 17    Sfc_Wind_Speed
!!$ 18    Sfc_Wind_Direction
!!$ 19    Land_Type_Index
      
      do k1 = 1,geovals%nvar
         varname = geovals%variables%fldnames(k1)
         print *, k1, varname
      end do

      !** populate the atmosphere structures for CRTM (atm(k1), for the k1-th profile)
      do k1 = 1,N_PROFILES
         lfound = ufo_geovals_get_var(geovals, var_tv, geoval)
         atm(k1)%Temperature(1:N_LAYERS) = geoval%vals(:,k1) 
         !print *, 'Temperature:', atm(k1)%Temperature(1:2), geoval%vals(1:2,k1)
         lfound = ufo_geovals_get_var(geovals, var_prs, geoval)
         atm(k1)%Pressure(1:N_LAYERS) = geoval%vals(:,k1) 
         !print *, 'Pressure:', atm(k1)%Pressure(1:2), geoval%vals(1:2,k1)
         lfound = ufo_geovals_get_var(geovals, var_prsi, geoval)
         atm(k1)%Level_Pressure(0:N_LAYERS) = geoval%vals(:,k1)
         !print *, 'level_pressure:', atm(k1)%Level_Pressure(0:1), geoval%vals(1:2,k1)
         atm(k1)%Climatology         = US_STANDARD_ATMOSPHERE
         atm(k1)%Absorber_Id(1:1)    = (/ H2O_ID /)
         atm(k1)%Absorber_Units(1:1) = (/ MASS_MIXING_RATIO_UNITS /)
         lfound = ufo_geovals_get_var(geovals, var_mixr, geoval)
         atm(k1)%Absorber(1:N_LAYERS,1)       = geoval%vals(:,k1) 
         !print *, 'water vapor:', atm(k1)%Absorber(1:2,1), geoval%vals(1:2,k1)
         atm(k1)%Absorber_Id(2:2)    = (/ O3_ID /)
         atm(k1)%Absorber_Units(2:2) = (/ VOLUME_MIXING_RATIO_UNITS /)
         lfound = ufo_geovals_get_var(geovals, var_oz, geoval)
         atm(k1)%Absorber(1:N_LAYERS,2)       = geoval%vals(:,k1) 
         !print *, 'Ozone:', atm(k1)%Absorber(1:2,2), geoval%vals(1:2,k1)

         atm(k1)%Absorber_Id(3:3)    = (/ CO2_ID /)
         atm(k1)%Absorber_Units(3:3) = (/ VOLUME_MIXING_RATIO_UNITS /)
         lfound = ufo_geovals_get_var(geovals, var_co2, geoval)
         atm(k1)%Absorber(1:N_LAYERS,3)       = geoval%vals(:,k1)

         atm(k1)%Cloud(1)%Type = WATER_CLOUD
         lfound = ufo_geovals_get_var(geovals, var_clw, geoval)
         atm(k1)%Cloud(1)%Water_Content = geoval%vals(:,k1)
         lfound = ufo_geovals_get_var(geovals, var_clwefr, geoval)
         atm(k1)%Cloud(1)%Effective_Radius = geoval%vals(:,k1)

         atm(k1)%Cloud(2)%Type = ICE_CLOUD
         lfound = ufo_geovals_get_var(geovals, var_cli, geoval)
         atm(k1)%Cloud(2)%Water_Content = geoval%vals(:,k1)
         lfound = ufo_geovals_get_var(geovals, var_cliefr, geoval)
         atm(k1)%Cloud(2)%Effective_Radius = geoval%vals(:,k1)


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

         !print *, 'Ozone:', atm(k1)%Absorber(1:2,2), geoval%vals(1:2,k1)
      enddo

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


      !       varname = geovals%variables%fldnames(1)
       !******                               123456789012345678901234'

      
      do k1 = 1,N_PROFILES
         sfc(k1)%Water_Type         = SEA_WATER_TYPE    !** NOTE: need to check how to determine fresh vs sea water types (salinity???)
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
      end do

    END SUBROUTINE Load_Sfc_Data

    !
    ! Internal subprogam to load some test geometry data
    !
    SUBROUTINE Load_Geom_Data()
      implicit none
      integer :: k1
      do k1 = 1,N_PROFILES
         geo(k1)%Sensor_Zenith_Angle = Radiance%datafix(k1)%satzen_ang
!YT ???         geo(k1)%Sensor_Scan_Angle   = Radiance%datafix(k1)%senscn_ang
         geo(k1)%Source_Zenith_Angle = Radiance%datafix(k1)%solzen_ang
         geo(k1)%Sensor_Azimuth_Angle = Radiance%datafix(k1)%satazm_ang
         geo(k1)%Source_Azimuth_Angle = Radiance%datafix(k1)%solazm_ang
         geo(k1)%Ifov = Radiance%datafix(k1)%senscn_pos
      enddo

    END SUBROUTINE Load_Geom_Data
    
  end subroutine ufo_radiance_eqv

  ! ------------------------------------------------------------------------------
  
end module ufo_radiance_mod
