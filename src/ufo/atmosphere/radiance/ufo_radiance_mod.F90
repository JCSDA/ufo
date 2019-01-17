! (C) Copyright 2017-2018 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

!> Fortran module to handle radiance observations

module ufo_radiance_mod

 use iso_c_binding
 use config_mod
 use kinds

 use ufo_geovals_mod, only: ufo_geovals, ufo_geoval, ufo_geovals_get_var
 use ufo_basis_mod, only: ufo_basis
 use ufo_vars_mod
 use ufo_radiance_utils_mod
! use obsspace_mod
 use crtm_module
 USE rttov_types
 USE rttov_const, ONLY : errorstatus_success, deg2rad
 USE rttov_unix_env

 implicit none
 private

 !> Fortran derived type for radiance trajectory

 type, extends(ufo_basis), public :: ufo_radiance
 private
  type(rad_conf) :: rc
 contains
   procedure :: setup  => ufo_radiance_setup
   procedure :: delete => ufo_radiance_delete
   procedure :: simobs => ufo_radiance_simobs
 end type ufo_radiance

contains

! ------------------------------------------------------------------------------

subroutine ufo_radiance_setup(self, c_conf)

implicit none
class(ufo_radiance), intent(inout) :: self
type(c_ptr),         intent(in)    :: c_conf

 call rad_conf_setup(self%rc,c_conf)

end subroutine ufo_radiance_setup

! ------------------------------------------------------------------------------

subroutine ufo_radiance_delete(self)

implicit none
class(ufo_radiance), intent(inout) :: self

 call rad_conf_delete(self%rc)

end subroutine ufo_radiance_delete

! ------------------------------------------------------------------------------

SUBROUTINE ufo_radiance_simobs(self, geovals, hofx, obss)
  CLASS(ufo_radiance),      INTENT(in) :: self
  TYPE(ufo_geovals),        INTENT(in) :: geovals
  REAL(c_double),        INTENT(inout) :: hofx(:)
  TYPE(c_ptr), VALUE,       INTENT(in) :: obss
  
  IF(TRIM(self%rc%rtmodel) == 'CRTM') THEN
    CALL ufo_radiance_simobs_crtm(self, geovals, hofx, obss)
  ELSEIF(TRIM(self%rc%rtmodel) == 'RTTOV') THEN
    CALL ufo_radiance_simobs_rttov(self, geovals, hofx, obss)
  ENDIF
  
END SUBROUTINE ufo_radiance_simobs

SUBROUTINE ufo_radiance_simobs_rttov(self, geovals, hofx, obss)

  USE ufo_radiance_utils_mod , ONLY : rttov_config

  IMPLICIT NONE

  CLASS(ufo_radiance),      INTENT(in) :: self
  TYPE(ufo_geovals),        INTENT(in) :: geovals
  REAL(c_double),        INTENT(inout) :: hofx(:)
  TYPE(c_ptr), VALUE,       INTENT(in) :: obss
 
  ! Local Variables
  CHARACTER(*), PARAMETER :: PROGRAM_NAME = 'ufo_radiance_mod.F90'
  CHARACTER(255) :: message, version
  INTEGER        :: err_stat, alloc_stat
  INTEGER        :: l, m, n, i, s, ierr
  TYPE(ufo_geoval), POINTER :: temp

  INTEGER :: n_Profiles
  INTEGER :: n_Layers
  INTEGER :: n_Channels

  ! ============================================================================
  ! STEP 3. **** DEFINE THE RTTOV INTERFACE STRUCTURES ****
  !

  TYPE(rttov_chanprof),    POINTER :: chanprof(:)    => NULL() ! Input channel/profile list
  LOGICAL(KIND=jplm),      POINTER :: calcemis(:)    => NULL() ! Flag to indicate calculation of emissivity within RTTOV
  TYPE(rttov_emissivity),  POINTER :: emissivity(:)  => NULL() ! Input/output surface emissivity
  TYPE(rttov_profile),     POINTER :: profiles(:)    => NULL() ! Input profiles
  TYPE(rttov_transmission)         :: transmission             ! Output transmittances
  TYPE(rttov_radiance)             :: radiance                 ! Output radiances

  INTEGER(KIND=jpim)               :: errorstatus              ! Return error status of RTTOV subroutine calls

  INTEGER :: i_inst,j , jch,n_levels, nch, nchannels, nchannels_total, nchans_inst, asw

  INCLUDE 'rttov_direct.interface'
  INCLUDE 'rttov_alloc_direct.interface'

 ! Get number of profile and layers from geovals
 ! ---------------------------------------------
 n_Profiles = geovals%nobs
 call ufo_geovals_get_var(geovals, var_tv, temp)
 n_Layers = temp%nval
 N_LEVELS = N_LAYERS + 1
 nullify(temp)

 hofx(:) = 0.0_kind_real
 errorstatus = 0_jpim
 nchannels_total = 0
 
 asw = 1

 IF( .NOT. rttov_config%rttov_is_setup) THEN
   CALL rttov_config%setup(self%rc, asw)
 ENDIF

 Sensor_Loop:DO i_inst = 1, self%rc%n_Sensors

  nchans_inst = rttov_config%rttov_coef_array(i_inst)%coef%fmv_chn
    
  ! Ensure the options and coefficients are consistent
  CALL rttov_user_options_checkinput(errorstatus, rttov_config%opts, rttov_config%rttov_coef_array(i_inst))
  IF (errorstatus /= errorstatus_success) THEN
    WRITE(*,*) 'error in rttov options'
    CALL rttov_exit(errorstatus)
  ENDIF

  ! --------------------------------------------------------------------------
  ! 3. Allocate RTTOV input and output structures
  ! --------------------------------------------------------------------------

  ! Determine the total number of radiances to simulate (nchanprof).
  ! In this example we simulate all specified channels for each profile, but
  ! in general one can simulate a different number of channels for each profile.
  nchannels = nchans_inst * N_PROFILES

  ! Allocate structures for rttov_direct
  CALL rttov_alloc_direct( &
        errorstatus,             &
        1_jpim,                  &  ! 1 => allocate
        N_PROFILES,              &
        nchannels,               &
        N_LEVELS,                &
        chanprof,                &
        rttov_config%opts,              &
        profiles,                &
        rttov_config%rttov_coef_array(i_inst),&
        transmission,            &
        radiance,                &
        calcemis=calcemis,       &
        emissivity=emissivity,   &
        init=.TRUE._jplm)

  IF (errorstatus /= errorstatus_success) THEN
    WRITE(*,*) 'allocation error for rttov_direct structures'
    CALL rttov_exit(errorstatus)
  ENDIF

  ! --------------------------------------------------------------------------
  ! 4. Build the list of profile/channel indices in chanprof
  ! --------------------------------------------------------------------------

  nch = 0_jpim
  DO j = 1, N_PROFILES
    DO jch = 1, nchans_inst
      nch = nch + 1_jpim
      chanprof(nch)%prof = j
      chanprof(nch)%chan = jch ! only all channels for now. Look at OPS for better implementation.
    ENDDO
  ENDDO

   !Assign the data from the GeoVaLs
   !--------------------------------
!   call Load_Atm_Data(n_Profiles,n_Layers,geovals,atm)
!   call Load_Sfc_Data(n_Profiles,n_Layers,n_Channels,geovals,sfc,chinfo,obss)
!   call Load_Geom_Data(obss,geo)

   !Assign the data from the GeoVaLs
   !--------------------------------

   CALL rttov_Load_Atm_Data(n_Profiles,n_Layers,geovals,obss,profiles)

   call rttov_Load_Geom_Data(obss,profiles)

  ! --------------------------------------------------------------------------
  ! 6. Specify surface emissivity and reflectance
  ! --------------------------------------------------------------------------

  ! In this example we have no values for input emissivities
  emissivity(:) % emis_in = 0._jprb

  ! Calculate emissivity within RTTOV where the input emissivity value is
  ! zero or less (all channels in this case)
  calcemis(:) = (emissivity(:) % emis_in <= 0._jprb)

  ! ! In this example we have no values for input reflectances
  ! reflectance(:) % refl_in = 0._jprb

  ! ! Calculate BRDF within RTTOV where the input BRDF value is zero or less
  ! ! (all channels in this case)
  ! calcrefl(:) = (reflectance(:) % refl_in <= 0._jprb)

   ! --------------------------------------------------------------------------
   ! 7. Call RTTOV forward model
   ! --------------------------------------------------------------------------
   CALL rttov_direct(                &
     errorstatus,              &! out   error flag
     chanprof,                 &! in    channel and profile index structure
     rttov_config%opts,               &! in    options structure
     profiles,                 &! in    profile array
     rttov_config%rttov_coef_array(i_inst), &! in    coefficients structure
     transmission,             &! inout computed transmittances
     radiance,                 &! inout computed radiances
     calcemis    = calcemis,   &! in    flag for internal emissivity calcs
     emissivity  = emissivity)!, &! inout input/output emissivities per channel
   ! calcrefl    = calcrefl,   &! in    flag for internal BRDF calcs
   ! reflectance = reflectance) ! inout input/output BRDFs per channel
       
   IF ( errorstatus /= errorstatus_success ) THEN
     message = 'Error calling RTTOV Forward Model for amsua'!//TRIM(SENSOR_ID(n))
     WRITE(*,*) message
     STOP
   END IF

   ! Put simulated brightness temperature into hofx
   ! ----------------------------------------------
   
   hofx(nchannels_total+1:nchannels_total+nchannels) = radiance%bt(1:nchannels)

   nchannels_total = nchannels_total + nchannels

  ! Allocate structures for rttov_direct
  CALL rttov_alloc_direct( &
        errorstatus,             &
        0_jpim,                  &  ! 0 => deallocate
        N_PROFILES,              &
        nchannels,               &
        N_LEVELS,                &
        chanprof,                &
        rttov_config%opts,              &
        profiles,                &
        rttov_config%rttov_coef_array(i_inst),&
        transmission,            &
        radiance,                &
        calcemis=calcemis,       &
        emissivity=emissivity)

  IF (errorstatus /= errorstatus_success) THEN
    WRITE(*,*) 'deallocation error for rttov_direct structures'
    CALL rttov_exit(errorstatus)
  ENDIF

 end do Sensor_Loop


END SUBROUTINE ufo_radiance_simobs_rttov

SUBROUTINE ufo_radiance_simobs_crtm(self, geovals, hofx, obss)

implicit none
class(ufo_radiance),      intent(in) :: self
type(ufo_geovals),        intent(in) :: geovals
real(c_double),        intent(inout) :: hofx(:)
type(c_ptr), value,       intent(in) :: obss

! Local Variables
character(*), parameter :: PROGRAM_NAME = 'ufo_radiance_mod.F90'
character(255) :: message, version
integer        :: err_stat, alloc_stat
integer        :: l, m, n, i, s
type(ufo_geoval), pointer :: temp

integer :: n_Profiles
integer :: n_Layers
integer :: n_Channels

! Define the "non-demoninational" arguments
type(CRTM_ChannelInfo_type)             :: chinfo(self%rc%n_Sensors)
type(CRTM_Geometry_type),   allocatable :: geo(:)

! Define the FORWARD variables
type(CRTM_Atmosphere_type), allocatable :: atm(:)
type(CRTM_Surface_type),    allocatable :: sfc(:)
type(CRTM_RTSolution_type), allocatable :: rts(:,:)


 ! Get number of profile and layers from geovals
 ! ---------------------------------------------
 n_Profiles = geovals%nobs
 call ufo_geovals_get_var(geovals, var_tv, temp)
 n_Layers = temp%nval
 nullify(temp)


 ! Program header
 ! --------------
 call CRTM_Version( Version )
 call Program_Message( PROGRAM_NAME, &
                       'Check/example program for the CRTM Forward and K-Matrix functions using '//&
                       trim(self%rc%ENDIAN_type)//' coefficient datafiles', &
                       'CRTM Version: '//TRIM(Version) )


 ! Initialise all the sensors at once
 ! ----------------------------------
 !** NOTE: CRTM_Init points to the various binary files needed for CRTM.  See the
 !**       CRTM_Lifecycle.f90 for more details.

 write( *,'(/5x,"Initializing the CRTM...")' )
 err_stat = CRTM_Init( self%rc%SENSOR_ID, &
            chinfo, &
            File_Path=trim(self%rc%COEFFICIENT_PATH), &
            Quiet=.TRUE.)
 if ( err_stat /= SUCCESS ) THEN
   message = 'Error initializing CRTM'
   call Display_Message( PROGRAM_NAME, message, FAILURE )
   stop
 end if


 ! Loop over all sensors. Not necessary if we're calling CRTM for each sensor
 ! ----------------------------------------------------------------------------
 Sensor_Loop:do n = 1, self%rc%n_Sensors


   ! Determine the number of channels for the current sensor
   ! -------------------------------------------------------
   N_Channels = CRTM_ChannelInfo_n_Channels(chinfo(n))


   ! Allocate the ARRAYS
   ! -------------------
   allocate( geo( n_Profiles ),               &
             atm( n_Profiles ),               &
             sfc( n_Profiles ),               &
             rts( N_Channels, n_Profiles ),   &
             STAT = alloc_stat )
   if ( alloc_stat /= 0 ) THEN
      message = 'Error allocating structure arrays'
      call Display_Message( PROGRAM_NAME, message, FAILURE )
      stop
   end if


   ! Create the input FORWARD structure (atm)
   ! ----------------------------------------
   call CRTM_Atmosphere_Create( atm, n_Layers, self%rc%n_Absorbers, self%rc%n_Clouds, self%rc%n_Aerosols )
   if ( ANY(.NOT. CRTM_Atmosphere_Associated(atm)) ) THEN
      message = 'Error allocating CRTM Forward Atmosphere structure'
      CALL Display_Message( PROGRAM_NAME, message, FAILURE )
      STOP
   END IF


   ! Create the input FORWARD structure (sfc)
   ! ----------------------------------------
   call CRTM_Surface_Create(sfc, N_Channels)
   IF ( ANY(.NOT. CRTM_Surface_Associated(sfc)) ) THEN
      message = 'Error allocating CRTM Surface structure'
      CALL Display_Message( PROGRAM_NAME, message, FAILURE )
      STOP
   END IF


   !Assign the data from the GeoVaLs
   !--------------------------------
   call Load_Atm_Data(n_Profiles,n_Layers,geovals,atm)
   call Load_Sfc_Data(n_Profiles,n_Layers,n_Channels,geovals,sfc,chinfo,obss)
   call Load_Geom_Data(obss,geo)


   ! Call THE CRTM inspection
   ! ------------------------
   call CRTM_Atmosphere_Inspect(atm(12))
   call CRTM_Surface_Inspect(sfc(12))
   call CRTM_Geometry_Inspect(geo(12))
   call CRTM_ChannelInfo_Inspect(chinfo(1))


   ! Call the forward model call for each sensor
   ! -------------------------------------------
   err_stat = CRTM_Forward( atm        , &  ! Input
                            sfc        , &  ! Input
                            geo        , &  ! Input
                            chinfo(n:n), &  ! Input
                            rts          )  ! Output
   if ( err_stat /= SUCCESS ) THEN
      message = 'Error calling CRTM Forward Model for '//TRIM(self%rc%SENSOR_ID(n))
      call Display_Message( PROGRAM_NAME, message, FAILURE )
      stop
   end if


   ! Put simulated brightness temperature into hofx
   ! ----------------------------------------------

   !Set to zero and initializ counter
   hofx(:) = 0.0_kind_real
   i = 1

   do m = 1, n_Profiles
     do l = 1, N_Channels

       hofx(i) = rts(l,m)%Brightness_Temperature
       i = i + 1

     end do
   end do


   ! Deallocate the structures
   ! -------------------------
   call CRTM_Geometry_Destroy(geo)
   call CRTM_Atmosphere_Destroy(atm)
   call CRTM_RTSolution_Destroy(rts)
   call CRTM_Surface_Destroy(sfc)


   ! Deallocate all arrays
   ! ---------------------
   deallocate(geo, atm, sfc, rts, STAT = alloc_stat)
   if ( alloc_stat /= 0 ) THEN
      message = 'Error deallocating structure arrays'
      call Display_Message( PROGRAM_NAME, message, FAILURE )
      stop
   end if

 end do Sensor_Loop


 ! Destroy CRTM instance
 ! ---------------------
 write( *, '( /5x, "Destroying the CRTM..." )' )
 err_stat = CRTM_Destroy( chinfo )
 if ( err_stat /= SUCCESS ) THEN
    message = 'Error destroying CRTM'
    call Display_Message( PROGRAM_NAME, message, FAILURE )
    stop
 end if

END SUBROUTINE ufo_radiance_simobs_crtm

! ------------------------------------------------------------------------------

end module ufo_radiance_mod
