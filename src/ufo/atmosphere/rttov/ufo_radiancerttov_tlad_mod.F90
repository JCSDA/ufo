! (C) Copyright 2017-2018 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

!> Fortran module for radiancerttov tl/ad observation operator

module ufo_radiancerttov_tlad_mod

 use iso_c_binding
 use config_mod
 use kinds

 use ufo_geovals_mod, only: ufo_geovals, ufo_geoval, ufo_geovals_get_var
 use ufo_basis_tlad_mod, only: ufo_basis_tlad
 use ufo_vars_mod
 use ufo_radiancerttov_utils_mod

 USE rttov_types
 USE rttov_const, ONLY : errorstatus_success, deg2rad
 USE rttov_unix_env

 implicit none
 private

 !> Fortran derived type for the tl/ad observation operator
 TYPE, EXTENDS(ufo_basis_tlad), PUBLIC :: ufo_radiancerttov_tlad
 private
  integer :: nprofiles
  integer :: nlayers
  integer :: nchannels

  type(rad_conf) :: rc
  TYPE(rttov_profile), POINTER :: profiles_k(:) => NULL()
  TYPE(rttov_chanprof), POINTER :: chanprof(:) => NULL()

 contains
  procedure :: setup  => ufo_radiancerttov_tlad_setup
  procedure :: delete  => ufo_radiancerttov_tlad_delete
  procedure :: settraj => ufo_radiancerttov_tlad_settraj
  procedure :: simobs_tl  => ufo_radiancerttov_simobs_tl
  procedure :: simobs_ad  => ufo_radiancerttov_simobs_ad
 end type ufo_radiancerttov_tlad

contains

! ------------------------------------------------------------------------------
subroutine ufo_radiancerttov_tlad_setup(self, c_conf)
implicit none
class(ufo_radiancerttov_tlad), intent(inout) :: self
type(c_ptr),              intent(in)    :: c_conf

CALL rad_conf_setup(self % rc,c_conf)

end subroutine ufo_radiancerttov_tlad_setup

! ------------------------------------------------------------------------------
subroutine ufo_radiancerttov_tlad_delete(self)
implicit none
class(ufo_radiancerttov_tlad), intent(inout) :: self

 self % ltraj = .false.
 call rad_conf_delete(self % rc)

end subroutine ufo_radiancerttov_tlad_delete

! ------------------------------------------------------------------------------
subroutine ufo_radiancerttov_tlad_settraj(self, geovals, obss)

USE ufo_radiancerttov_utils_mod , ONLY : config_rttov

implicit none

class(ufo_radiancerttov_tlad), intent(inout) :: self
type(ufo_geovals),       intent(in)          :: geovals
type(c_ptr), value,      intent(in)          :: obss

! Local Variables
character(*), parameter                      :: PROGRAM_NAME = 'ufo_radiancerttov_mod.F90'
character(255)                               :: message, version
integer                                      :: err_stat, alloc_stat
INTEGER                                      :: i_inst,j , jch,nlevels, nch, nchannels, nchans_total, nchans_inst, asw, ierr
type(ufo_geoval), pointer                    :: temp

LOGICAL(KIND=jplm),      POINTER             :: calcemis(:)    => NULL()   ! Flag to indicate calculation of emissivity within RTTOV

TYPE(rttov_emissivity),  POINTER             :: emissivity(:)  => NULL()   ! Input/output surface emissivity
TYPE(rttov_profile),     POINTER             :: profiles(:)    => NULL()   ! Input profiles
TYPE(rttov_transmission)                     :: transmission               ! Output transmittances
TYPE(rttov_radiance)                         :: radiance                   ! Output radiances

TYPE(rttov_emissivity),  POINTER             :: emissivity_k(:)  => NULL() ! Input/output surface emissivity
TYPE(rttov_transmission)                     :: transmission_k             ! Output transmittances
TYPE(rttov_radiance)                         :: radiance_k                 ! Output radiances

INTEGER(KIND=jpim)                           :: errorstatus                ! Return error status of RTTOV subroutine calls

INCLUDE 'rttov_k.interface'
INCLUDE 'rttov_alloc_k.interface'

 ! Get number of profile and layers from geovals
 ! ---------------------------------------------
 self % nprofiles = geovals % nobs
 CALL ufo_geovals_get_var(geovals, var_ts, temp)
 self % nlayers = temp % nval

 nlevels = self % nlayers + 1
 nullify(temp)

 errorstatus = 0_jpim
 nchans_total = 0

 asw = 1

 IF( .NOT. config_rttov % rttov_is_setup) THEN
   CALL config_rttov % setup(self % rc, asw)
 ENDIF

 Sensor_Loop:DO i_inst = 1, self % rc % nSensors

  nchans_inst = config_rttov % rttov_coef_array(i_inst) % coef % fmv_chn
  self % nchannels = nchans_inst
    
  ! Ensure the options and coefficients are consistent
  CALL rttov_user_options_checkinput(errorstatus, config_rttov % opts, config_rttov % rttov_coef_array(i_inst))
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
  nchannels = nchans_inst * self % NPROFILES

  ! Allocate structures for rttov_direct
  CALL rttov_alloc_k(                            &
        errorstatus,                             &
        1_jpim,                                  &  ! 1 => allocate
        self % NPROFILES,                        &
        nchannels,                               &
        NLEVELS,                                 &
        self % chanprof,                         &
        config_rttov % opts,                     &
        profiles,                                &
        self % profiles_k,                       &
        config_rttov % rttov_coef_array(i_inst), &
        transmission,                            &
        transmission_k,                          &
        radiance,                                &
        radiance_k,                              &
        calcemis=calcemis,                       &
        emissivity=emissivity,                   &
        emissivity_k=emissivity_k,               &
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
  DO j = 1, self % NPROFILES
    DO jch = 1, nchans_inst
      nch = nch + 1_jpim
      self % chanprof(nch) % prof = j
      self % chanprof(nch) % chan = jch ! only all channels for now. Look at OPS for better implementation.
    ENDDO
  ENDDO

   !Assign the data from the GeoVaLs
   !--------------------------------

   CALL load_atm_data_rttov(self % nprofiles,self % nlayers,geovals,obss,profiles)

   call load_geom_data_rttov(obss,profiles)

  ! --------------------------------------------------------------------------
  ! 6. Specify surface emissivity and reflectance
  ! --------------------------------------------------------------------------

  ! In this example we have no values for input emissivities
  emissivity(:) % emis_in = 0._jprb

  ! Calculate emissivity within RTTOV where the input emissivity value is
  ! zero or less (all channels in this case)
  calcemis(:) = (emissivity(:) % emis_in <= 0._jprb)

   ! --------------------------------------------------------------------------
   ! 7. Call RTTOV forward model
   ! --------------------------------------------------------------------------
   CALL rttov_k(                                                                 &
     errorstatus,                                                                &! out   error flag
     self % chanprof(nchans_total + 1:nchans_total + nchans_inst),               &! in    channel and profile index structure
     config_rttov % opts,                                                        &! in    options structure
     profiles,                                                                   &! in    profile array
     self % profiles_k(nchans_total + 1:nchans_total + nchans_inst),             &! in    profile array
     config_rttov % rttov_coef_array(i_inst),                                    &! in    coefficients structure
     transmission,                                                               &! inout computed transmittances
     transmission_k,                                                             &! inout computed transmittances
     radiance,                                                                   &! inout computed radiances
     radiance_k,                                                                 &! inout computed radiances
     calcemis    = calcemis(nchans_total + 1:nchans_total + nchans_inst),        &! in    flag for internal emissivity calcs
     emissivity  = emissivity(nchans_total + 1:nchans_total + nchans_inst),      &!, &! inout input/output emissivities per channel
     emissivity_k = emissivity_k(nchans_total + 1:nchans_total + nchans_inst))!, &! inout input/output emissivities per channel
       
   IF ( errorstatus /= errorstatus_success ) THEN
     message = 'Error calling RTTOV K Model for amsua'!//TRIM(SENSOR_ID(n))
     WRITE(*,*) message
     STOP
   END IF

   nchans_total = nchans_total + nchannels

 end do Sensor_Loop

 ! Set flag that the tracectory was set
 ! ------------------------------------
 self % ltraj = .true.

end subroutine ufo_radiancerttov_tlad_settraj

! ------------------------------------------------------------------------------
subroutine ufo_radiancerttov_simobs_tl(self, geovals, hofx, obss)
implicit none
class(ufo_radiancerttov_tlad), intent(in)    :: self
type(ufo_geovals),       intent(in)    :: geovals
real(c_double),          intent(inout) :: hofx(:)
type(c_ptr), value,      intent(in)    :: obss

character(len=*), parameter :: myname_="ufo_radiancerttov_simobs_tl"
character(max_string) :: err_msg
INTEGER :: job, jprofile, jchannel, jlevel, ierr, ichan, prof
type(ufo_geoval), pointer :: tv_d


 ! Initial checks
 ! --------------

 ! Check if trajectory was set
 if (.not. self % ltraj) then
   write(err_msg,*) myname_, ' trajectory wasnt set!'
   call abor1_ftn(err_msg)
 endif

 ! Check if nobs is consistent in geovals & hofx
 if (geovals % nobs /= self % nprofiles) then
   write(err_msg,*) myname_, ' error: nobs inconsistent!'
   call abor1_ftn(err_msg)
 endif

 ! Initialize hofx
 ! ---------------
 hofx(:) = 0.0_kind_real

 ! Temperature
 ! -----------

 ! Get tv from geovals
 call ufo_geovals_get_var(geovals, var_ts, tv_d)

 ! Check model levels is consistent in geovals
 if (tv_d % nval /= self % nlayers) then
   write(err_msg,*) myname_, ' error: layers inconsistent!'
   call abor1_ftn(err_msg)
 endif

 DO ichan = 1, self % nprofiles * self % nchannels
   prof = self % chanprof(ichan) % prof
   hofx(ichan) = hofx(ichan) + &
     SUM(self % profiles_k(ichan) % t(2:) * tv_d % vals(:,prof))
 ENDDO

end subroutine ufo_radiancerttov_simobs_tl

! ------------------------------------------------------------------------------
subroutine ufo_radiancerttov_simobs_ad(self, geovals, hofx, obss)
implicit none
class(ufo_radiancerttov_tlad), intent(in)    :: self
type(ufo_geovals),       intent(inout) :: geovals
real(c_double),          intent(in)    :: hofx(:)
type(c_ptr), value,      intent(in)    :: obss


character(len=*), parameter :: myname_="ufo_radiancerttov_simobs_ad"
character(max_string) :: err_msg
INTEGER :: job, jprofile, jchannel, jlevel, ierr, ichan, prof
type(ufo_geoval), pointer :: tv_d


 ! Initial checks
 ! --------------

 ! Check if trajectory was set
 if (.not. self % ltraj) then
   write(err_msg,*) myname_, ' trajectory wasnt set!'
   call abor1_ftn(err_msg)
 endif

 ! Check if nobs is consistent in geovals & hofx
 if (geovals % nobs /= self % nprofiles) then
   write(err_msg,*) myname_, ' error: nobs inconsistent!'
   call abor1_ftn(err_msg)
 endif


 ! Temperature
 ! -----------

 ! Get tv from geovals
 call ufo_geovals_get_var(geovals, var_ts, tv_d)

 ! allocate if not yet allocated
 if (.not. allocated(tv_d % vals)) then
    tv_d % nobs = self % nprofiles
    tv_d % nval = self % nlayers
    allocate(tv_d % vals(tv_d % nval,tv_d % nobs))
    tv_d % vals = 0.0_kind_real
 endif

 DO ichan = 1, self % nprofiles * self % nchannels
   prof = self % chanprof(ichan) % prof
   tv_d % vals(:,prof) = tv_d % vals(:,prof) + &
     self % profiles_k(ichan) % t(2:) * hofx(ichan)
 ENDDO

 ! Once all geovals set replace flag
 ! ---------------------------------
 if (.not. geovals % linit ) geovals % linit=.true.

end subroutine ufo_radiancerttov_simobs_ad

! ------------------------------------------------------------------------------

end module ufo_radiancerttov_tlad_mod
