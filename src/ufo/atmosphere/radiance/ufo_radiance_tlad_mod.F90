! (C) Copyright 2017-2018 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

!> Fortran module to handle tl/ad for radiance observations

module ufo_radiance_tlad_mod

 use iso_c_binding
 use config_mod
 use kinds

 use ufo_geovals_mod, only: ufo_geovals, ufo_geoval, ufo_geovals_get_var
 use ufo_basis_tlad_mod, only: ufo_basis_tlad
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
 type, extends(ufo_basis_tlad), public :: ufo_radiance_tlad
 private
  type(rad_conf) :: rc
  integer :: n_Profiles
  integer :: n_Layers
  integer :: n_Channels
  type(CRTM_Atmosphere_type), allocatable :: atm_K(:,:)
  type(CRTM_Surface_type), allocatable :: sfc_K(:,:)
  TYPE(rttov_profile), POINTER :: profiles_k(:) => NULL()
  TYPE(rttov_chanprof), POINTER :: chanprof(:) => NULL()
 contains
  procedure :: setup  => ufo_radiance_tlad_setup
  procedure :: delete  => ufo_radiance_tlad_delete
  procedure :: settraj => ufo_radiance_tlad_settraj
  procedure :: simobs_tl  => ufo_radiance_simobs_tl
  procedure :: simobs_ad  => ufo_radiance_simobs_ad
 end type ufo_radiance_tlad

contains

! ------------------------------------------------------------------------------

subroutine ufo_radiance_tlad_setup(self, c_conf)

implicit none
class(ufo_radiance_tlad), intent(inout) :: self
type(c_ptr),              intent(in)    :: c_conf

 call rad_conf_setup(self%rc,c_conf)

end subroutine ufo_radiance_tlad_setup

! ------------------------------------------------------------------------------

subroutine ufo_radiance_tlad_delete(self)

implicit none
class(ufo_radiance_tlad), intent(inout) :: self

 self%ltraj = .false.
 call rad_conf_delete(self%rc)

 if (allocated(self%atm_k)) then
   call CRTM_Atmosphere_Destroy(self%atm_K)
   deallocate(self%atm_k)
 endif

 if (allocated(self%sfc_k)) then
   call CRTM_Surface_Destroy(self%sfc_K)
   deallocate(self%sfc_k)
 endif

end subroutine ufo_radiance_tlad_delete

! ------------------------------------------------------------------------------

SUBROUTINE ufo_radiance_tlad_settraj(self, geovals, obss)
  CLASS(ufo_radiance_tlad), INTENT(inout) :: self
  TYPE(ufo_geovals),        INTENT(in)    :: geovals
  TYPE(c_ptr), VALUE,       INTENT(in)    :: obss
  
  IF(TRIM(self%rc%rtmodel) == 'CRTM') THEN
    CALL ufo_radiance_tlad_settraj_crtm(self, geovals, obss)
  ELSEIF(TRIM(self%rc%rtmodel) == 'RTTOV') THEN
    CALL ufo_radiance_tlad_settraj_rttov(self, geovals, obss)
  ENDIF
  
END SUBROUTINE ufo_radiance_tlad_settraj

! ------------------------------------------------------------------------------

SUBROUTINE ufo_radiance_tlad_settraj_crtm(self, geovals, obss)

implicit none

class(ufo_radiance_tlad), intent(inout) :: self
type(ufo_geovals),        intent(in)    :: geovals
type(c_ptr), value,       intent(in)    :: obss

! Local Variables
character(*), parameter :: PROGRAM_NAME = 'ufo_radiance_mod.F90'
character(255) :: message, version
integer        :: err_stat, alloc_stat
integer        :: n, k1
type(ufo_geoval), pointer :: temp

! Define the "non-demoninational" arguments
type(CRTM_ChannelInfo_type)             :: chinfo(self%rc%n_Sensors)
type(CRTM_Geometry_type),   allocatable :: geo(:)

! Define the FORWARD variables
type(CRTM_Atmosphere_type), allocatable :: atm(:)
type(CRTM_Surface_type),    allocatable :: sfc(:)
type(CRTM_RTSolution_type), allocatable :: rts(:,:)

! Define the K-MATRIX variables
type(CRTM_RTSolution_type), allocatable :: rts_K(:,:)


 ! Get number of profile and layers from geovals
 ! ---------------------------------------------
 self%n_Profiles = geovals%nobs
 call ufo_geovals_get_var(geovals, var_ts, temp)
 self%n_Layers = temp%nval
 nullify(temp)

 ! Program header
 ! --------------
 call CRTM_Version( Version )
 call Program_Message( PROGRAM_NAME, &
                       'Check/example program for the CRTM Forward and K-Matrix (setTraj) functions using '//&
                       trim(self%rc%ENDIAN_type)//' coefficient datafiles', &
                       'CRTM Version: '//TRIM(Version) )


 ! Initialise all the sensors at once
 ! ----------------------------------
 !** NOTE: CRTM_Init points to the various binary files needed for CRTM.  See the
 !**       CRTM_Lifecycle.f90 for more details.

 write( *,'(/5x,"Initializing the CRTM (setTraj) ...")' )
 err_stat = CRTM_Init( self%rc%SENSOR_ID, &
            chinfo, &
            File_Path=trim(self%rc%COEFFICIENT_PATH), &
            Quiet=.TRUE.)
 if ( err_stat /= SUCCESS ) THEN
   message = 'Error initializing CRTM (setTraj)'
   call Display_Message( PROGRAM_NAME, message, FAILURE )
   stop
 end if

 ! Loop over all sensors. Not necessary if we're calling CRTM for each sensor
 ! ----------------------------------------------------------------------------
 Sensor_Loop:do n = 1, self%rc%n_Sensors


   ! Determine the number of channels for the current sensor
   ! -------------------------------------------------------
   self%N_Channels = CRTM_ChannelInfo_n_Channels(chinfo(n))


   ! Allocate the ARRAYS
   ! -------------------
   allocate( geo( self%n_Profiles )                         , &
             atm( self%n_Profiles )                         , &
             sfc( self%n_Profiles )                         , &
             rts( self%N_Channels, self%n_Profiles )        , &
             self%atm_K( self%N_Channels, self%n_Profiles ) , &
             self%sfc_K( self%N_Channels, self%n_Profiles ) , &
             rts_K( self%N_Channels, self%n_Profiles )      , &
             STAT = alloc_stat                                )
   if ( alloc_stat /= 0 ) THEN
      message = 'Error allocating structure arrays (setTraj)'
      call Display_Message( PROGRAM_NAME, message, FAILURE )
      stop
   end if


   ! Create the input FORWARD structure (atm)
   ! ----------------------------------------
   call CRTM_Atmosphere_Create( atm, self%n_Layers, self%rc%n_Absorbers, self%rc%n_Clouds, self%rc%n_Aerosols )
   if ( ANY(.NOT. CRTM_Atmosphere_Associated(atm)) ) THEN
      message = 'Error allocating CRTM Forward Atmosphere structure (setTraj)'
      CALL Display_Message( PROGRAM_NAME, message, FAILURE )
      STOP
   END IF


   ! Create the input FORWARD structure (sfc)
   ! ----------------------------------------
   call CRTM_Surface_Create(sfc, self%N_Channels)
   IF ( ANY(.NOT. CRTM_Surface_Associated(sfc)) ) THEN
      message = 'Error allocating CRTM Surface structure (setTraj)'
      CALL Display_Message( PROGRAM_NAME, message, FAILURE )
      STOP
   END IF


   ! Create output K-MATRIX structure (atm)
   ! --------------------------------------
   call CRTM_Atmosphere_Create( self%atm_K, self%n_Layers, self%rc%n_Absorbers, self%rc%n_Clouds, self%rc%n_Aerosols )
   if ( ANY(.NOT. CRTM_Atmosphere_Associated(self%atm_K)) ) THEN
      message = 'Error allocating CRTM K-matrix Atmosphere structure (setTraj)'
      CALL Display_Message( PROGRAM_NAME, message, FAILURE )
      STOP
   END IF


   ! Create output K-MATRIX structure (sfc)
   ! --------------------------------------
   call CRTM_Surface_Create(self%sfc_K, self%N_Channels)
   IF ( ANY(.NOT. CRTM_Surface_Associated(self%sfc_K)) ) THEN
      message = 'Error allocating CRTM K-matrix Surface structure (setTraj)'
      CALL Display_Message( PROGRAM_NAME, message, FAILURE )
      STOP
   END IF


   !Assign the data from the GeoVaLs
   !--------------------------------
   call Load_Atm_Data(self%N_PROFILES,self%N_LAYERS,geovals,atm)
   call Load_Sfc_Data(self%N_PROFILES,self%N_LAYERS,self%N_Channels,geovals,sfc,chinfo,obss)
   call Load_Geom_Data(obss,geo)


   ! Zero the K-matrix OUTPUT structures
   ! -----------------------------------
   call CRTM_Atmosphere_Zero( self%atm_K )
   call CRTM_Surface_Zero( self%sfc_K )


   ! Inintialize the K-matrix INPUT so that the results are dTb/dx
   ! -------------------------------------------------------------
   rts_K%Radiance               = ZERO
   rts_K%Brightness_Temperature = ONE


   ! Call the K-matrix model
   ! -----------------------
   err_stat = CRTM_K_Matrix( atm         , &  ! FORWARD  Input
                             sfc         , &  ! FORWARD  Input
                             rts_K       , &  ! K-MATRIX Input
                             geo         , &  ! Input
                             chinfo(n:n) , &  ! Input
                             self%atm_K  , &  ! K-MATRIX Output
                             self%sfc_K  , &  ! K-MATRIX Output
                             rts           )  ! FORWARD  Output
   if ( err_stat /= SUCCESS ) THEN
      message = 'Error calling CRTM (setTraj) K-Matrix Model for '//TRIM(self%rc%SENSOR_ID(n))
      call Display_Message( PROGRAM_NAME, message, FAILURE )
      stop
   end if


   ! Deallocate the structures
   ! -------------------------
   call CRTM_Geometry_Destroy(geo)
   call CRTM_Atmosphere_Destroy(atm)
   call CRTM_RTSolution_Destroy(rts_K)
   call CRTM_RTSolution_Destroy(rts)
   call CRTM_Surface_Destroy(sfc)


   ! Deallocate all arrays
   ! ---------------------
   deallocate(geo, atm, sfc, rts, rts_K, STAT = alloc_stat)
   if ( alloc_stat /= 0 ) THEN
      message = 'Error deallocating structure arrays (setTraj)'
      call Display_Message( PROGRAM_NAME, message, FAILURE )
      stop
   end if

 end do Sensor_Loop


 ! Destroy CRTM instance
 ! ---------------------
 write( *, '( /5x, "Destroying the CRTM (setTraj)..." )' )
 err_stat = CRTM_Destroy( chinfo )
 if ( err_stat /= SUCCESS ) THEN
    message = 'Error destroying CRTM (setTraj)'
    call Display_Message( PROGRAM_NAME, message, FAILURE )
    stop
 end if


 ! Set flag that the tracectory was set
 ! ------------------------------------
 self%ltraj = .true.

END SUBROUTINE ufo_radiance_tlad_settraj_crtm

! ------------------------------------------------------------------------------

SUBROUTINE ufo_radiance_tlad_settraj_rttov(self, geovals, obss)

USE ufo_radiance_utils_mod , ONLY : rttov_config

implicit none

class(ufo_radiance_tlad), intent(inout) :: self
type(ufo_geovals),        intent(in)    :: geovals
type(c_ptr), value,       intent(in)    :: obss

! Local Variables
character(*), parameter :: PROGRAM_NAME = 'ufo_radiance_mod.F90'
character(255) :: message, version
integer        :: err_stat, alloc_stat
INTEGER :: i_inst,j , jch,n_levels, nch, nchannels, nchans_total, nchans_inst, asw, ierr
type(ufo_geoval), pointer :: temp

! ============================================================================
! STEP 3. **** DEFINE THE RTTOV INTERFACE STRUCTURES ****
!

LOGICAL(KIND=jplm),      POINTER :: calcemis(:)    => NULL() ! Flag to indicate calculation of emissivity within RTTOV

TYPE(rttov_emissivity),  POINTER :: emissivity(:)  => NULL() ! Input/output surface emissivity
TYPE(rttov_profile),     POINTER :: profiles(:)    => NULL() ! Input profiles
TYPE(rttov_transmission)         :: transmission             ! Output transmittances
TYPE(rttov_radiance)             :: radiance                 ! Output radiances

TYPE(rttov_emissivity),  POINTER :: emissivity_k(:)  => NULL() ! Input/output surface emissivity
TYPE(rttov_transmission)         :: transmission_k             ! Output transmittances
TYPE(rttov_radiance)             :: radiance_k                 ! Output radiances

INTEGER(KIND=jpim)               :: errorstatus              ! Return error status of RTTOV subroutine calls

INCLUDE 'rttov_k.interface'
INCLUDE 'rttov_alloc_k.interface'


 ! Get number of profile and layers from geovals
 ! ---------------------------------------------
 self%n_Profiles = geovals%nobs
 call ufo_geovals_get_var(geovals, var_tv, temp)
 self%n_Layers = temp%nval

 n_levels = self%n_layers + 1
 nullify(temp)

 asw = 1

 IF( .NOT. rttov_config%rttov_is_setup) THEN
   CALL rttov_config % setup(self%rc, asw)
 ENDIF

 Sensor_Loop:DO i_inst = 1, self%rc%n_Sensors

  nchans_inst = rttov_config%rttov_coef_array(i_inst)%coef%fmv_chn
  self%n_channels = nchans_inst
    
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
  nchannels = nchans_inst * self%N_PROFILES

  ! Allocate structures for rttov_direct
  CALL rttov_alloc_k( &
        errorstatus,             &
        1_jpim,                  &  ! 1 => allocate
        self%N_PROFILES,              &
        nchannels,               &
        N_LEVELS,                &
        self%chanprof,           &
        rttov_config%opts,       &
        profiles,                &
        self%profiles_k,         &
        rttov_config%rttov_coef_array(i_inst),&
        transmission,            &
        transmission_k,          &
        radiance,                &
        radiance_k,              &
        calcemis=calcemis,       &
        emissivity=emissivity,   &
        emissivity_k=emissivity_k,   &
        init=.TRUE._jplm)

  IF (errorstatus /= errorstatus_success) THEN
    WRITE(*,*) 'allocation error for rttov_k structures'
    CALL rttov_exit(errorstatus)
  ENDIF

  emissivity_k % emis_out = 0
  emissivity_k % emis_in = 0
  emissivity % emis_out = 0
  
  ! Inintialize the K-matrix INPUT so that the results are dTb/dx
  ! -------------------------------------------------------------

  radiance_k % bt(:) = 1
  radiance_k % total(:) = 1

  ! --------------------------------------------------------------------------
  ! 4. Build the list of profile/channel indices in chanprof
  ! --------------------------------------------------------------------------

  nch = 0_jpim
  DO j = 1, self%N_PROFILES
    DO jch = 1, nchans_inst
      nch = nch + 1_jpim
      self%chanprof(nch)%prof = j
      self%chanprof(nch)%chan = jch ! only all channels for now. Look at OPS for better implementation.
    ENDDO
  ENDDO

   !Assign the data from the GeoVaLs
   !--------------------------------

   CALL rttov_Load_Atm_Data(self%n_Profiles,self%n_Layers,geovals,obss,profiles)

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
   CALL rttov_k(                &
     errorstatus,              &! out   error flag
     self%chanprof(nchans_total + 1:nchans_total + nchans_inst),  &! in    channel and profile index structure
     rttov_config%opts,               &! in    options structure
     profiles,                 &! in    profile array
     self%profiles_k(nchans_total + 1:nchans_total + nchans_inst), &! in    profile array
     rttov_config%rttov_coef_array(i_inst), &! in    coefficients structure
     transmission,             &! inout computed transmittances
     transmission_k,           &! inout computed transmittances
     radiance,                 &! inout computed radiances
     radiance_k,               &! inout computed radiances
     calcemis    = calcemis(nchans_total + 1:nchans_total + nchans_inst),   &! in    flag for internal emissivity calcs
     emissivity  = emissivity(nchans_total + 1:nchans_total + nchans_inst),  &!, &! inout input/output emissivities per channel
     emissivity_k = emissivity_k(nchans_total + 1:nchans_total + nchans_inst))!, &! inout input/output emissivities per channel
   ! calcrefl    = calcrefl,   &! in    flag for internal BRDF calcs
   ! reflectance = reflectance) ! inout input/output BRDFs per channel
       
   IF ( errorstatus /= errorstatus_success ) THEN
     message = 'Error calling RTTOV K Model for amsua'!//TRIM(SENSOR_ID(n))
     WRITE(*,*) message
     STOP
   END IF

   nchans_total = nchans_total + nchannels

 end do Sensor_Loop

 ! Set flag that the tracectory was set
 ! ------------------------------------
 self%ltraj = .true.

END SUBROUTINE ufo_radiance_tlad_settraj_rttov

! ------------------------------------------------------------------------------

SUBROUTINE ufo_radiance_simobs_tl(self, geovals, hofx, obss)

implicit none
class(ufo_radiance_tlad), intent(in) :: self
type(ufo_geovals),        intent(in) :: geovals
real(c_double),        intent(inout) :: hofx(:)
type(c_ptr), value,    intent(in)    :: obss

character(len=*), parameter :: myname_="ufo_radiance_simobs_tl"
character(max_string) :: err_msg
INTEGER :: job, jprofile, jchannel, jlevel, ierr, ichan, prof
type(ufo_geoval), pointer :: tv_d


 ! Initial checks
 ! --------------

 ! Check if trajectory was set
 if (.not. self%ltraj) then
   write(err_msg,*) myname_, ' trajectory wasnt set!'
   call abor1_ftn(err_msg)
 endif

 ! Check if nobs is consistent in geovals & hofx
 if (geovals%nobs /= self%n_Profiles) then
   write(err_msg,*) myname_, ' error: nobs inconsistent!'
   call abor1_ftn(err_msg)
 endif

 ! Initialize hofx
 ! ---------------
 hofx(:) = 0.0_kind_real


 ! Temperature
 ! -----------

 ! Get t from geovals
 call ufo_geovals_get_var(geovals, var_ts, tv_d)

 ! Check model levels is consistent in geovals & crtm
 if (tv_d%nval /= self%n_Layers) then
   write(err_msg,*) myname_, ' error: layers inconsistent!'
   call abor1_ftn(err_msg)
 endif

 IF(self%rc%rtmodel == 'RTTOV') THEN
   DO ichan = 1, self%n_profiles * self%n_Channels
     prof = self%chanprof(ichan)%prof
     hofx(ichan) = hofx(ichan) + &
       SUM(self%profiles_k(ichan)%t(2:) * tv_d%vals(:,prof))
   ENDDO
 ELSE
   ! Multiply by Jacobian and add to hofx
   job = 0
   DO jprofile = 1, self%n_Profiles
     DO jchannel = 1, self%n_Channels
       job = job + 1
       DO jlevel = 1, tv_d%nval
         hofx(job) = hofx(job) + &
           self%atm_K(jchannel,jprofile)%Temperature(jlevel) * tv_d%vals(jlevel,jprofile)
       ENDDO
     ENDDO
   ENDDO
 ENDIF

end subroutine ufo_radiance_simobs_tl

! ------------------------------------------------------------------------------

subroutine ufo_radiance_simobs_ad(self, geovals, hofx, obss)

implicit none
class(ufo_radiance_tlad), intent(in) :: self
type(ufo_geovals),     intent(inout) :: geovals
real(c_double),           intent(in) :: hofx(:)
type(c_ptr), value,    intent(in)    :: obss

character(len=*), parameter :: myname_="ufo_radiance_simobs_ad"
character(max_string) :: err_msg
INTEGER :: job, jprofile, jchannel, jlevel, ierr, ichan, prof
type(ufo_geoval), pointer :: tv_d


 ! Initial checks
 ! --------------

 ! Check if trajectory was set
 if (.not. self%ltraj) then
   write(err_msg,*) myname_, ' trajectory wasnt set!'
   call abor1_ftn(err_msg)
 endif

 ! Check if nobs is consistent in geovals & hofx
 if (geovals%nobs /= self%n_Profiles) then
   write(err_msg,*) myname_, ' error: nobs inconsistent!'
   call abor1_ftn(err_msg)
 endif


 ! Temperature
 ! -----------

 ! Get t from geovals
 call ufo_geovals_get_var(geovals, var_ts, tv_d)

 ! allocate if not yet allocated
 if (.not. allocated(tv_d%vals)) then
    tv_d%nobs = self%n_Profiles
    tv_d%nval = self%n_Layers
    allocate(tv_d%vals(tv_d%nval,tv_d%nobs))
    tv_d%vals = 0.0_kind_real
 endif


 IF(self%rc%rtmodel == 'RTTOV') THEN
   DO ichan = 1, self%n_profiles * self%n_Channels
     prof = self%chanprof(ichan)%prof
     tv_d%vals(:,prof) = tv_d%vals(:,prof) + &
       self%profiles_k(ichan)%t(2:) * hofx(ichan)
   ENDDO
 ELSE
   ! Multiply by Jacobian and add to hofx (adjoint)
   job = 0
   DO jprofile = 1, self%n_Profiles
     DO jchannel = 1, self%n_Channels
       job = job + 1
       DO jlevel = 1, tv_d%nval
         tv_d%vals(jlevel,jprofile) = tv_d%vals(jlevel,jprofile) + &
           self%atm_K(jchannel,jprofile)%Temperature(jlevel) * hofx(job)
       ENDDO
     ENDDO
   ENDDO
 ENDIF




 ! Once all geovals set replace flag
 ! ---------------------------------
 if (.not. geovals%linit ) geovals%linit=.true.


end subroutine ufo_radiance_simobs_ad

! ------------------------------------------------------------------------------

end module ufo_radiance_tlad_mod
