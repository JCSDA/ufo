! (C) Copyright 2017-2018 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

!> Fortran module for radiancerttov observation operator

module ufo_radiancerttov_mod

 use iso_c_binding
 use config_mod
 use kinds

 use ufo_geovals_mod, only: ufo_geovals, ufo_geoval, ufo_geovals_get_var
 use ufo_basis_mod, only: ufo_basis
 use ufo_vars_mod
 use obsspace_mod
 use ufo_radiancerttov_utils_mod

 USE rttov_types
 USE rttov_const, ONLY : errorstatus_success, deg2rad
 USE rttov_unix_env

 implicit none
 private

!> Fortran derived type for the observation type
! TODO: fill in if needed
 TYPE, EXTENDS(ufo_basis), PUBLIC :: ufo_radiancerttov
   PRIVATE
   TYPE(rad_conf) :: rc
 CONTAINS
   PROCEDURE :: setup  => ufo_radiancerttov_setup
   PROCEDURE :: delete => ufo_radiancerttov_delete
   PROCEDURE :: simobs => ufo_radiancerttov_simobs
 END TYPE ufo_radiancerttov

contains

! ------------------------------------------------------------------------------
! TODO: add setup of your observation operator (optional)
subroutine ufo_radiancerttov_setup(self, c_conf)
implicit none
class(ufo_radiancerttov), intent(inout) :: self
type(c_ptr),        intent(in)    :: c_conf

CALL rad_conf_setup(self % rc,c_conf)

end subroutine ufo_radiancerttov_setup

! ------------------------------------------------------------------------------
! TODO: add cleanup of your observation operator (optional)
subroutine ufo_radiancerttov_delete(self)
implicit none
class(ufo_radiancerttov), intent(inout) :: self

CALL rad_conf_delete(self % rc)

end subroutine ufo_radiancerttov_delete

! ------------------------------------------------------------------------------
! TODO: put code for your nonlinear observation operator in this routine
! Code in this routine is for radiancerttov only, please remove and replace
SUBROUTINE ufo_radiancerttov_simobs(self, geovals, hofx, obss)

USE ufo_radiancerttov_utils_mod , ONLY : config_rttov

implicit none
class(ufo_radiancerttov), intent(in) :: self
type(ufo_geovals),  intent(in)       :: geovals
real(c_double),     intent(inout)    :: hofx(:)
type(c_ptr), value, intent(in)       :: obss

  ! Local Variables
  CHARACTER(*), PARAMETER          :: PROGRAM_NAME = 'ufo_radiancerttov_mod.F90'
  CHARACTER(255)                   :: message, version
  INTEGER                          :: err_stat, alloc_stat
  INTEGER                          :: l, m, n, i, s, ierr
  TYPE(ufo_geoval), POINTER        :: temp

  INTEGER                          :: nprofiles
  INTEGER                          :: nlayers

  TYPE(rttov_chanprof),    POINTER :: chanprof(:)    => NULL() ! Input channel/profile list
  LOGICAL(KIND=jplm),      POINTER :: calcemis(:)    => NULL() ! Flag to indicate calculation of emissivity within RTTOV
  TYPE(rttov_emissivity),  POINTER :: emissivity(:)  => NULL() ! Input/output surface emissivity
  TYPE(rttov_profile),     POINTER :: profiles(:)    => NULL() ! Input profiles
  TYPE(rttov_transmission)         :: transmission             ! Output transmittances
  TYPE(rttov_radiance)             :: radiance                 ! Output radiances

  INTEGER(KIND=jpim)               :: errorstatus              ! Return error status of RTTOV subroutine calls

  INTEGER                          :: i_inst,j , jch, nlevels, nch, nchannels, nchans_total, nchans_inst, asw

  INCLUDE 'rttov_direct.interface'
  INCLUDE 'rttov_alloc_direct.interface'

 ! Get number of profile and layers from geovals
 ! ---------------------------------------------
 nprofiles = geovals % nobs
 call ufo_geovals_get_var(geovals, var_ts, temp)
 nlayers = temp % nval
 nlevels = nlayers + 1
 nullify(temp)

 hofx(:) = 0.0_kind_real
 errorstatus = 0_jpim
 nchans_total = 0
 
 asw = 1

 IF( .NOT. config_rttov % rttov_is_setup) THEN
   CALL config_rttov % setup(self % rc, asw)
 ENDIF

 Sensor_Loop:DO i_inst = 1, self % rc % nSensors

  nchans_inst = config_rttov % rttov_coef_array(i_inst) % coef % fmv_chn
    
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
  nchannels = nchans_inst * nprofiles

  ! Allocate structures for rttov_direct
  CALL rttov_alloc_direct( &
        errorstatus,             &
        1_jpim,                  &  ! 1 => allocate
        nprofiles,              &
        nchannels,               &
        nlevels,                &
        chanprof,                &
        config_rttov % opts,              &
        profiles,                &
        config_rttov % rttov_coef_array(i_inst),&
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
  DO j = 1, nprofiles
    DO jch = 1, nchans_inst
      nch = nch + 1_jpim
      chanprof(nch) % prof = j
      chanprof(nch) % chan = jch ! only all channels for now. Look at OPS for better implementation.
    ENDDO
  ENDDO

   !Assign the data from the GeoVaLs
   !--------------------------------

   CALL load_atm_data_rttov(nprofiles,nlayers,geovals,obss,profiles)

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
   CALL rttov_direct(                         &
     errorstatus,                             &! out   error flag
     chanprof,                                &! in    channel and profile index structure
     config_rttov % opts,                     &! in    options structure
     profiles,                                &! in    profile array
     config_rttov % rttov_coef_array(i_inst), &! in    coefficients structure
     transmission,                            &! inout computed transmittances
     radiance,                                &! inout computed radiances
     calcemis    = calcemis,                  &! in    flag for internal emissivity calcs
     emissivity  = emissivity)!,              &! inout input/output emissivities per channel
       
   IF ( errorstatus /= errorstatus_success ) THEN
     message = 'Error calling RTTOV Forward Model for amsua'!//TRIM(SENSOR_ID(n))
     WRITE(*,*) message
     STOP
   END IF

   ! Put simulated brightness temperature into hofx
   ! ----------------------------------------------
   
   hofx(nchans_total+1:nchans_total+nchannels) = radiance % bt(1:nchannels)

   nchans_total = nchans_total + nchannels

  ! Allocate structures for rttov_direct
  CALL rttov_alloc_direct(                       &
        errorstatus,                             &
        0_jpim,                                  &  ! 0 => deallocate
        nprofiles,                              &
        nchannels,                               &
        nlevels,                                &
        chanprof,                                &
        config_rttov % opts,                     &
        profiles,                                &
        config_rttov % rttov_coef_array(i_inst), &
        transmission,                            &
        radiance,                                &
        calcemis = calcemis,                     &
        emissivity = emissivity)

  IF (errorstatus /= errorstatus_success) THEN
    WRITE(*,*) 'deallocation error for rttov_direct structures'
    CALL rttov_exit(errorstatus)
  ENDIF

END DO Sensor_Loop

END SUBROUTINE ufo_radiancerttov_simobs

! ------------------------------------------------------------------------------

end module ufo_radiancerttov_mod
