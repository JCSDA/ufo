! (C) Copyright 2017 UCAR
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 

!> Fortran module to handle aod observations

MODULE ufo_aod_mod
  
  use ioda_obsdb_mod
  use ioda_obs_vectors
  use ufo_vars_mod
  use ioda_locs_mod
  use ufo_geovals_mod
  use kinds
  USE ufo_aod_misc
  use crtm_module
  USE ufo_basis_mod, only: ufo_basis

  implicit none

  public :: ufo_aod

  private
  integer, parameter :: max_string=800  
  LOGICAL, PARAMETER :: ice4qsat=.TRUE.

  REAL(fp),PARAMETER:: &
       &ttp = 2.7316e+2_fp, &
       &psat = 6.1078e+2_fp,&
       &rd = 2.8705e+2_fp,&
       &rv = 4.6150e+2_fp,&
       &cv = 7.1760e+2_fp,&
       &cliq = 4.1855e+3_fp,&
       &csol = 2.1060e+3_fp,&
       &cvap = 1.8460e+3_fp,&
       &hvap = 2.5000e+6_fp,&
       &hfus = 3.3358e+5_fp,&
       &grav = 9.81_fp
  
  REAL(fp),PARAMETER ::  &
       &tmix = ttp-20._fp,&
       &hsub = hvap+hfus,&
       &eps = rd/rv,&
       &omeps=one-eps,&
       &dldt =cvap-cliq,&
       &dldti = cvap-csol,&
       &xa = -(dldt/rv),&
       &xai = -(dldti/rv),&
       &xb = xa+hvap/(rv*ttp),&
       &xbi = xai+hsub/(rv*ttp)
  
  LOGICAL :: flip_vertical

  
  !> Fortran derived type for aod trajectory
  type, extends(ufo_basis) :: ufo_aod
  contains
    procedure :: eqv => ufo_aod_eqv
  end type ufo_aod

contains
  
! ------------------------------------------------------------------------------

  SUBROUTINE ufo_aod_eqv(self, geovals, hofx, obss) 
    implicit none
    class(ufo_aod),    intent(in)    :: self
    type(ufo_geovals), intent(in)    :: geovals
    type(obs_vector),  intent(inout) :: hofx
    type(ioda_obsdb), target, intent(in)  :: obss

    type(obs_vector) :: TmpOvec
    real(kind_real), allocatable :: Aod_Obs(:,:)
    real(kind_real), allocatable :: Omg_Aod(:,:)

    !*************************************************************************************
    !******* Begin CRTM block ************************************************************
    !*************************************************************************************

    ! --------------------------
    ! Some non-CRTM-y Parameters
    ! --------------------------
    CHARACTER(*), PARAMETER :: PROGRAM_NAME   = 'ufo_aod_mod.F90'
    
    
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


    INTEGER, PARAMETER :: N_ABSORBERS = 2  !** UFO
    INTEGER, PARAMETER :: N_CLOUDS    = 0  !** UFO
    INTEGER, PARAMETER :: n_aerosols_gocart_crtm = 14,N_AEROSOLS=n_aerosols_gocart_crtm


    INTEGER :: N_PROFILES 
    INTEGER :: N_LAYERS   

    
    ! Sensor information
    INTEGER     , PARAMETER :: N_SENSORS = 1  
    CHARACTER(*), PARAMETER :: SENSOR_ID(N_SENSORS) = (/'v.viirs-m_npp'/)  !** UFO to provide sensor name

    CHARACTER(max_string) :: err_msg

  
    ! ---------
    ! Local Variables
    ! ---------
    CHARACTER(256) :: message, version
    INTEGER        :: err_stat, alloc_stat
    INTEGER        :: n_channels
    INTEGER        :: l, m, n, nc, i,k
    real(fp)       :: cf

   
    
    ! ============================================================================
    ! STEP 3. **** DEFINE THE CRTM INTERFACE STRUCTURES ****
    !
    ! 3a. Define the "non-demoninational" arguments
    ! ---------------------------------------------
    TYPE(CRTM_ChannelInfo_type)             :: chinfo(N_SENSORS)
    
    
    ! 3b. Define the FORWARD variables
    ! --------------------------------
    TYPE(CRTM_Atmosphere_type), ALLOCATABLE :: atm(:)
    TYPE(CRTM_RTSolution_type), ALLOCATABLE :: rts(:,:)
    
    
    ! 3c. Define the K-MATRIX variables
    ! ---------------------------------
    TYPE(CRTM_Atmosphere_type), ALLOCATABLE :: atm_K(:,:)
    TYPE(CRTM_RTSolution_type), ALLOCATABLE :: rts_K(:,:)
    ! ============================================================================

    TYPE(CRTM_Aerosol_type), ALLOCATABLE :: aerosols(:)

    TYPE(ufo_geoval), pointer :: geoval
    character(MAXVARLEN) :: varname
    integer              :: ivar
    integer              :: ierr

    integer              :: nobs
    integer              :: nlocs
    REAL(fp), allocatable :: rmse(:)
    real(fp), allocatable :: diff(:,:)


    ! Program header
    ! --------------

    nobs=obss%nobs; nlocs=obss%nlocs

    n_profiles=geovals%nobs

    varname=var_aerosols(1)

    call ufo_geovals_get_var(geovals,varname, geoval,status=ierr)
    IF (ierr==0) THEN
       n_layers=SIZE(geoval%vals,1)
    ELSE
       err_msg=TRIM(varname)//' not found - Stopping'
       CALL abor1_ftn(err_msg)
    ENDIF
    

    ALLOCATE(atm(n_profiles),aerosols(n_aerosols))

    
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

       CALL CRTM_RTSolution_Create(rts, n_Layers )
       CALL CRTM_RTSolution_Create(rts_k, n_Layers )

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
       
       IF (n_Aerosols > 0)  CALL Load_Aerosol_Data()

       CALL CRTM_Atmosphere_Zero( atm_K )
       
       ! 7b. Inintialize the K-matrix INPUT so
       !     that the results are dTb/dx
       ! -------------------------------------
       rts_K%Radiance               = ZERO
       rts_K%Brightness_Temperature = ZERO

       DO m = 1, N_PROFILES
          DO l = 1, n_Channels
             rts_K(l,m)%layer_optical_depth = ONE
          ENDDO
       ENDDO

       ! ==========================================================================

       ! ==========================================================================
       ! STEP 8. **** CALL THE CRTM FUNCTIONS FOR THE CURRENT SENSOR ****
       !
       call CRTM_Atmosphere_Inspect(atm)
       call CRTM_ChannelInfo_Inspect(chinfo(1))

! 8b.1 The K-matrix model for AOD
! ----------------------
       err_stat = CRTM_AOD_K( atm,    &  ! FORWARD  Input
            rts_K                   , &  ! K-MATRIX Input
            chinfo(n:n)             , &  ! Input
            rts                     , &  ! FORWARD  Output
            atm_k        )               ! K-MATRIX Output
       IF ( err_stat /= SUCCESS ) THEN
          message = 'Error calling CRTM K-Matrix Model for '//TRIM(SENSOR_ID(n))
          CALL Display_Message( PROGRAM_NAME, message, FAILURE )
          STOP
       END IF

       ! ============================================================================
       ! 8c. **** OUTPUT THE RESULTS TO SCREEN (optional) ****
       !
       ! User should read the user guide or the source code of the routine
       ! CRTM_RTSolution_Inspect in the file CRTM_RTSolution_Define.f90 to
       ! select the needed variables for outputs.  These variables are contained
       ! in the structure RTSolution.

       ALLOCATE(diff(n_channels,n_profiles),rmse(n_channels))

       allocate(Aod_Obs(n_channels, n_profiles))
       allocate(Omg_Aod(n_channels, n_profiles))
       call ioda_obsvec_setup(TmpOvec, obss%nobs)
       call ioda_obsdb_var_to_ovec(obss, TmpOvec, "Observation")
       Aod_Obs = reshape(TmpOvec%values, (/n_channels, n_profiles/))
       call ioda_obsdb_var_to_ovec(obss, TmpOvec, "Obs_Minus_Forecast_unadjusted")
       Omg_Aod = reshape(TmpOvec%values, (/n_channels, n_profiles/))

       rmse = 0
       DO m = 1, N_PROFILES
          DO l = 1, n_Channels
             diff(l,m) = SUM(rts(l,m)%layer_optical_depth(:)) - (Aod_Obs(l,m) - Omg_Aod(l,m))
             rmse(l) = rmse(l) + diff(l,m)**2
          END DO
       ENDDO

       rmse=SQRT(rmse/n_profiles)

       PRINT *,'N_profiles', N_PROFILES
       DO l = 1, n_Channels
          PRINT *, 'Channel: ',l
          PRINT *, 'Max difference: ', MAXVAL(ABS(diff(l,:)))
          PRINT *, 'RMSE: ', rmse(l)
       ENDDO

       DEALLOCATE(diff,rmse)

       deallocate(Aod_Obs)
       deallocate(Omg_Aod)
       call ioda_obsvec_delete(TmpOvec)

       ! output to hofx structure   
       hofx%values(:) = 0.0
       i = 1
       do m = 1, N_PROFILES
         do l = 1, n_Channels
            hofx%values(i) = SUM(rts(l,m)%layer_optical_depth)
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
       
       !** NOTE: Not 100% clear if any of the RTS structures need to be destroyed. 
       
       ! 9b. Deallocate the arrays !** NOTE: this is required
       ! -------------------------
       DEALLOCATE(rts, rts_K, atm_K, STAT = alloc_stat)
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


!!$ 1   Temperature
!!$ 2   Water vapor
!!$ 3   Pressure
!!$ 4   Level pressure

      !** populate the atmosphere structures for CRTM (atm(k1), for the k1-th profile)


      k1=1

      varname=var_prs
      call ufo_geovals_get_var(geovals, varname, geoval)

      IF (geoval%vals(1,k1) > geoval%vals(N_LAYERS,k1)) THEN 
         flip_vertical=.TRUE.
      ELSE
         flip_vertical=.FALSE.
      ENDIF

      DO k1 = 1,N_PROFILES

         varname=var_t
         call ufo_geovals_get_var(geovals, varname, geoval)
         IF (flip_vertical) THEN
            atm(k1)%Temperature(1:N_LAYERS) = geoval%vals(N_LAYERS:1:-1,k1)
         ELSE             
            atm(k1)%Temperature(1:N_LAYERS) = geoval%vals(:,k1) 
         ENDIF
!         print *, 'Temperature:', atm(k1)%Temperature(1:2), geoval%vals(1:2,k1)

         varname=var_prs
         call ufo_geovals_get_var(geovals, varname, geoval)
         IF (flip_vertical) THEN
            atm(k1)%Pressure(1:N_LAYERS) = geoval%vals(N_LAYERS:1:-1,k1)
         ELSE
            atm(k1)%Pressure(1:N_LAYERS) = geoval%vals(:,k1) 
         ENDIF
!         print *, 'Pressure:', atm(k1)%Pressure(1:2), geoval%vals(1:2,k1)


         varname=var_prsi
         call ufo_geovals_get_var(geovals, varname, geoval)
         IF (flip_vertical) THEN
            atm(k1)%Level_Pressure(0:N_LAYERS) = geoval%vals(N_LAYERS+1:1:-1,k1)
         ELSE
            atm(k1)%Level_Pressure(0:N_LAYERS) = geoval%vals(:,k1)
         ENDIF
!         print *, 'level_pressure:', atm(k1)%Level_Pressure(0:1), geoval%vals(1:2,k1)

         atm(k1)%Climatology         = US_STANDARD_ATMOSPHERE

         atm(k1)%Absorber_Id(1:1)    = (/ H2O_ID /)
         atm(k1)%Absorber_Units(1:1) = (/ MASS_MIXING_RATIO_UNITS /)
         varname=var_mixr
         call ufo_geovals_get_var(geovals, varname, geoval)
         IF (flip_vertical) THEN
            atm(k1)%Absorber(1:N_LAYERS,1)       = geoval%vals(N_LAYERS:1:-1,k1)
         ELSE
            atm(k1)%Absorber(1:N_LAYERS,1)       = geoval%vals(:,k1) 
         ENDIF
!         print *, 'water vapor:', atm(k1)%Absorber(1:2,1), geoval%vals(1:2,k1)

         atm(k1)%Absorber_Id(2:2)    = (/ O3_ID /)
         atm(k1)%Absorber_Units(2:2) = (/ VOLUME_MIXING_RATIO_UNITS /)
         atm(k1)%Absorber(1:N_LAYERS,2)=1.e-10

      ENDDO

      
    END SUBROUTINE Load_Atm_Data
    
    SUBROUTINE Load_Aerosol_Data()

      USE crtm_aerosolcoeff, ONLY: AeroC

      REAL(fp), DIMENSION(N_LAYERS) :: ugkg_kgm2,qsat,rh,prsl,tsen,p25

! Local variables
      INTEGER :: nc, NL
      INTEGER :: k1, k2

! 4a.1 Profile #1
! ---------------
! ...Profile and absorber definitions (fake/placeholder()

      CHARACTER(len=MAXVARLEN) :: varname

      INTEGER :: n_aerosols_all

      INTEGER :: i,ii,k,m

      DO m=1,N_PROFILES

!ug2kg && hPa2Pa
         DO k=1,N_LAYERS
            ugkg_kgm2(k)=1.0e-9_fp*(atm(m)%Level_Pressure(k)-&
                 &atm(m)%Level_Pressure(k-1))*100_fp/grav
            prsl(k)=atm(m)%Pressure(N_LAYERS-k+1)*0.1_fp ! must be in cb for genqsat
            tsen(k)=atm(m)%Temperature(N_LAYERS-k+1)
         ENDDO

         CALL genqsat(qsat,tsen,prsl,n_layers,ice4qsat)

         DO k=1,N_LAYERS
            rh(k)=atm(m)%Absorber(k,1)*1.e-3_fp/&
                 &qsat(N_LAYERS-k+1)
         ENDDO

         n_aerosols_all=SIZE(var_aerosols)
         
         IF (n_aerosols_all == naerosols_gocart_esrl) THEN
            
            DO i=1,n_aerosols_all
               varname=var_aerosols(i)
               IF (TRIM(varname) == 'p25') THEN
                  call ufo_geovals_get_var(geovals,varname, geoval)
                  
                  IF (flip_vertical) THEN
                     p25(1:N_LAYERS)=geoval%vals(N_LAYERS:1:-1,m)
                  ELSE
                     p25(1:N_LAYERS)=geoval%vals(1:N_LAYERS,m)
                  ENDIF

                  p25=MAX(p25*ugkg_kgm2,zero)
                  
                  EXIT
               ENDIF
            ENDDO
        
         ELSE
            p25=0_fp
         ENDIF

         ii=0

         DO i=1,n_aerosols_all
            varname=var_aerosols(i)
            IF (TRIM(varname) == 'p25') CYCLE
            ii=ii+1
            call ufo_geovals_get_var(geovals,varname, geoval)

            IF (flip_vertical) THEN
               atm(m)%aerosol(ii)%Concentration(1:N_LAYERS)=geoval%vals(N_LAYERS:1:-1,m)
            ELSE
               atm(m)%aerosol(ii)%Concentration(1:N_LAYERS)=geoval%vals(1:N_LAYERS,m)
            ENDIF

            atm(m)%aerosol(ii)%Concentration=MAX(atm(m)%aerosol(ii)%Concentration*ugkg_kgm2,&
                 &small_value)

            SELECT CASE ( TRIM(varname))
            CASE ('sulf')
               atm(m)%aerosol(ii)%type  = SULFATE_AEROSOL
!rh needs to be from top to bottom
               DO k=1,N_LAYERS
                  atm(m)%Aerosol(ii)%Effective_Radius(k)=&
                       &GOCART_Aerosol_size(atm(m)%aerosol(ii)%type, &
                       &rh(k))
               ENDDO

            CASE ('bc1')
               atm(m)%aerosol(ii)%type  = BLACK_CARBON_AEROSOL
               atm(m)%Aerosol(ii)%Effective_Radius(:)=&
                    &AeroC%Reff(1,atm(m)%aerosol(ii)%type)
            CASE ('bc2')
               atm(m)%aerosol(ii)%type  = BLACK_CARBON_AEROSOL
               DO k=1,N_LAYERS
                  atm(m)%Aerosol(ii)%Effective_Radius(k)=&
                       &GOCART_Aerosol_size(atm(m)%aerosol(ii)%type, &
                       &rh(k))
               ENDDO

            CASE ('oc1')
               atm(m)%aerosol(ii)%type  = ORGANIC_CARBON_AEROSOL
               atm(m)%Aerosol(ii)%Effective_Radius(:)=&
                    &AeroC%Reff(1,atm(m)%aerosol(ii)%type)
            CASE ('oc2')
               atm(m)%aerosol(ii)%type  = ORGANIC_CARBON_AEROSOL
               DO k=1,N_LAYERS
                  atm(m)%Aerosol(ii)%Effective_Radius(k)=&
                       &GOCART_Aerosol_size(atm(m)%aerosol(ii)%type, &
                       &rh(k))
               ENDDO

            CASE ('dust1')
               atm(m)%aerosol(ii)%type  = DUST_AEROSOL
               atm(m)%Aerosol(ii)%Effective_Radius(:)=0.55_fp
               atm(m)%aerosol(ii)%Concentration=&
                    &atm(m)%aerosol(ii)%Concentration+0.78_fp*p25
            CASE ('dust2')
               atm(m)%aerosol(ii)%type  = DUST_AEROSOL
               atm(m)%Aerosol(ii)%Effective_Radius(:)=1.4_fp
               atm(m)%aerosol(ii)%Concentration=&
                    &atm(m)%aerosol(ii)%Concentration+0.22_fp*p25
            CASE ('dust3')
               atm(m)%aerosol(ii)%type  = DUST_AEROSOL
               atm(m)%Aerosol(ii)%Effective_Radius(:)=2.4_fp
            CASE ('dust4')
               atm(m)%aerosol(ii)%type  = DUST_AEROSOL
               atm(m)%Aerosol(ii)%Effective_Radius(:)=4.5_fp
            CASE ('dust5')
               atm(m)%aerosol(ii)%type  = DUST_AEROSOL
               atm(m)%Aerosol(ii)%Effective_Radius(:)=8.0_fp

            CASE ('seas1')
               atm(m)%aerosol(ii)%type  = SEASALT_SSAM_AEROSOL
               DO k=1,N_LAYERS
                  atm(m)%Aerosol(ii)%Effective_Radius(k)=&
                       &GOCART_Aerosol_size(atm(m)%aerosol(ii)%type, &
                       &rh(k))
               ENDDO
            CASE ('seas2')
               atm(m)%aerosol(ii)%type  = SEASALT_SSCM1_AEROSOL
               DO k=1,N_LAYERS
                  atm(m)%Aerosol(ii)%Effective_Radius(k)=&
                       &GOCART_Aerosol_size(atm(m)%aerosol(ii)%type, &
                       &rh(k))
               ENDDO
            CASE ('seas3')
               atm(m)%aerosol(ii)%type  = SEASALT_SSCM2_AEROSOL
               DO k=1,N_LAYERS
                  atm(m)%Aerosol(ii)%Effective_Radius(k)=&
                       &GOCART_Aerosol_size(atm(m)%aerosol(ii)%type, &
                       &rh(k))
               ENDDO
            CASE ('seas4')
               atm(m)%aerosol(ii)%type  = SEASALT_SSCM3_AEROSOL
               DO k=1,N_LAYERS
                  atm(m)%Aerosol(ii)%Effective_Radius(k)=&
                       &GOCART_Aerosol_size(atm(m)%aerosol(ii)%type, &
                       &rh(k))
               ENDDO
            END SELECT

         ENDDO

      ENDDO
      
    END SUBROUTINE Load_Aerosol_Data
    
    FUNCTION GOCART_Aerosol_size( itype,  & ! Input
         eh ) & ! Input in 0-1
                           result( R_eff  )   ! in micrometer
      USE crtm_aerosolcoeff, ONLY: AeroC
      IMPLICIT NONE
      
!
!   modified from a function provided by Quanhua Liu
!
      INTEGER ,INTENT(in) :: itype
      REAL(fp)    ,INTENT(in) :: eh
      
      INTEGER :: j1,j2,k
      REAL(fp)    :: h1
      REAL(fp)    :: R_eff

      j2 = 0
      IF ( eh <= AeroC%RH(1) ) THEN
         j1 = 1
      ELSE IF ( eh >= AeroC%RH(AeroC%n_RH) ) THEN
         j1 = AeroC%n_RH
      ELSE
         DO k = 1, AeroC%n_RH-1
            IF ( eh < AeroC%RH(k+1) .AND. eh > AeroC%RH(k) ) THEN
               j1 = k
               j2 = k+1
               h1 = (eh-AeroC%RH(k))/(AeroC%RH(k+1)-AeroC%RH(k))
               EXIT
            ENDIF
         ENDDO
      ENDIF
      
      IF ( j2 == 0 ) THEN
         R_eff = AeroC%Reff(j1,itype )
      ELSE
         R_eff = (1.0_fp-h1)*AeroC%Reff(j1,itype ) + h1*AeroC%Reff(j2,itype )
      ENDIF
      
      RETURN
      
    END FUNCTION GOCART_Aerosol_size

    SUBROUTINE genqsat(qsat,tsen,prsl,nsig,ice)

!   input argument list:
!     tsen      - input sensibile temperature field (nlat,nlon,nsig)
!     prsl      - input layer mean pressure field (nlat,nlon,nsig)
!     nsig      - number of levels                              
!     ice       - logical flag:  T=include ice and ice-water effects,
!                 depending on t, in qsat calcuations.
!                 otherwise, compute qsat with respect to water surface
!
!   output argument list:
!     qsat      - saturation specific humidity (output)
!
! remarks: see modules used
!
! attributes:
!   language: f90
!   machine:  ibm RS/6000 SP
!

      IMPLICIT NONE

      LOGICAL                               ,INTENT(in   ) :: ice
      REAL(fp),DIMENSION(nsig), INTENT(  out) :: qsat
      REAL(fp),DIMENSION(nsig),INTENT(in   ) :: tsen,prsl
      INTEGER                       ,INTENT(in   ) :: nsig


      INTEGER k
      REAL(fp) pw,tdry,tr,es,es2
      REAL(fp) w,onep3,esmax
      REAL(fp) desidt,deswdt,dwdt,desdt,esi,esw
      REAL(fp) :: mint,estmax
      INTEGER :: lmint

      onep3 = 1.e3_fp

      mint=340._fp
      lmint=1

      DO k=1,nsig
         IF((prsl(k) < 30._fp .AND.  &
              prsl(k) > 2._fp) .AND.  &
              tsen(k) < mint)THEN
            lmint=k
            mint=tsen(k)
         END IF
      END DO

      tdry = mint
      tr = ttp/tdry

      IF (tdry >= ttp .OR. .NOT. ice) THEN
         estmax = psat * (tr**xa) * EXP(xb*(one-tr))
      ELSEIF (tdry < tmix) THEN
         estmax = psat * (tr**xai) * EXP(xbi*(one-tr))
      ELSE
         w  = (tdry - tmix) / (ttp - tmix)
         estmax =  w * psat * (tr**xa) * EXP(xb*(one-tr)) &
              + (one-w) * psat * (tr**xai) * EXP(xbi*(one-tr))
      ENDIF

      DO k = 1,nsig

         tdry = tsen(k)
         tr = ttp/tdry
         IF (tdry >= ttp .OR. .NOT. ice) THEN
            es = psat * (tr**xa) * EXP(xb*(one-tr))
         ELSEIF (tdry < tmix) THEN
            es = psat * (tr**xai) * EXP(xbi*(one-tr))
         ELSE
            esw = psat * (tr**xa) * EXP(xb*(one-tr)) 
            esi = psat * (tr**xai) * EXP(xbi*(one-tr)) 
            w  = (tdry - tmix) / (ttp - tmix)
            es =  w * psat * (tr**xa) * EXP(xb*(one-tr)) &
                 + (one-w) * psat * (tr**xai) * EXP(xbi*(one-tr))
         ENDIF

         pw = onep3*prsl(k)
         esmax = es
         IF(lmint < k)THEN
            esmax=0.1_fp*pw
            esmax=MIN(esmax,estmax)
         END IF
         es2=MIN(es,esmax)
         qsat(k) = eps * es2 / (pw - omeps * es2)

      END DO

      RETURN

    END SUBROUTINE genqsat

  END SUBROUTINE ufo_aod_eqv

! ------------------------------------------------------------------------------


END MODULE ufo_aod_mod
