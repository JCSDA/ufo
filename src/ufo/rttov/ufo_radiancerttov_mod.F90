! (C) Copyright 2017-2018 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

!> Fortran module for radiancerttov observation operator

module ufo_radiancerttov_mod

 use fckit_configuration_module, only: fckit_configuration
 use iso_c_binding
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
 type, public :: ufo_radiancerttov
   private
   type(rttov_conf) :: conf
 CONTAINS
   procedure :: setup  => ufo_radiancerttov_setup
   procedure :: delete => ufo_radiancerttov_delete
   procedure :: simobs => ufo_radiancerttov_simobs
 end type ufo_radiancerttov

contains

! ------------------------------------------------------------------------------
subroutine ufo_radiancerttov_setup(self, f_conf)
implicit none
class(ufo_radiancerttov), intent(inout) :: self
type(fckit_configuration), intent(in)   :: f_conf

! Local variables
type(fckit_configuration) :: f_confOpts

! Setup conf
call f_conf % get_or_die("ObsOptions", f_confOpts)
call rttov_conf_setup(self % conf, f_confOpts)

end subroutine ufo_radiancerttov_setup

! ------------------------------------------------------------------------------
subroutine ufo_radiancerttov_delete(self)
implicit none
class(ufo_radiancerttov), intent(inout) :: self

 call rttov_conf_delete(self%conf)

end subroutine ufo_radiancerttov_delete

! ------------------------------------------------------------------------------
SUBROUTINE ufo_radiancerttov_simobs(self, geovals, hofx, obss, channels, ob_info)

use fckit_log_module, only : fckit_log
use ufo_rttovonedvarcheck_ob_mod

implicit none
class(ufo_radiancerttov),    intent(in) :: self
type(ufo_geovals),           intent(in) :: geovals
real(c_double),              intent(inout) :: hofx(:)
type(c_ptr), value,          intent(in) :: obss
integer(c_int),              intent(in) :: channels(:)  ! List of channels to use
type(ufo_rttovonedvarcheck_ob), optional, intent(in) :: ob_info  ! Used for onedvarcheck 

! Local Variables
character(*), parameter          :: PROGRAM_NAME = 'ufo_radiancerttov_mod.F90'
character(255)                   :: message, version
integer                          :: err_stat, alloc_stat
integer                          :: l, m, n, i, s, ierr
type(ufo_geoval), pointer        :: temp

integer                          :: nprofiles

type(rttov_chanprof),    pointer :: chanprof(:)    => NULL() ! Input channel/profile list
logical(kind=jplm),      pointer :: calcemis(:)    => NULL() ! Flag to indicate calculation of emissivity within RTTOV
type(rttov_emissivity),  pointer :: emissivity(:)  => NULL() ! Input/output surface emissivity
type(rttov_profile),     pointer :: profiles(:)    => NULL() ! Input profiles
type(rttov_transmission)         :: transmission             ! Output transmittances
type(rttov_radiance)             :: radiance                 ! Output radiances

integer(kind=jpim)               :: errorstatus              ! Return error status of RTTOV subroutine calls

integer                          :: i_inst,j , jch, nlevels, nch, nchans_total, nchans_inst, asw

integer                          :: nprof_sim, nchan_sim, nchan_max_sim, nprof_max_sim
integer                          :: prof_start, prof_end

character(MAXVARLEN) :: varname

INCLUDE 'rttov_direct.interface'
INCLUDE 'rttov_print_profile.interface'
INCLUDE 'rttov_alloc_direct.interface'
INCLUDE 'rttov_user_profile_checkinput.interface'

! Get number of profile and layers from geovals
! ---------------------------------------------

nprofiles = geovals % nlocs
varname = 'air_pressure' !var_prsi
call ufo_geovals_get_var(geovals, varname, temp)
nlevels = temp % nval ! lfric passing nlevels
nullify(temp)

hofx(:) = 0.0_kind_real
errorstatus = 0_jpim
nchans_total = 0

nchan_max_sim = 40000 ! Maximum number of channels to pass to RTTOV to simulate

if( .NOT. config_rttov % rttov_is_setup) then
  asw = 1
  call config_rttov % setup(self % conf, asw)
else
  write(*,*) "Config rttov already setup"
end if

Sensor_Loop:do i_inst = 1, self % conf % nSensors

!  nchans_inst = config_rttov % rttov_coef_array(i_inst) % coef % fmv_chn
  nchans_inst = SIZE(channels)

  ! Ensure the options and coefficients are consistent
  call rttov_user_options_checkinput(errorstatus, config_rttov % opts, &
                                     config_rttov % rttov_coef_array(i_inst))

  if (errorstatus /= errorstatus_success) then
    write(message,'(A, I6)') 'after rttov_user_options_checkinput: error = ',&
      errorstatus
    call fckit_log%info(message)
    call rttov_exit(errorstatus)
  end if
  ! --------------------------------------------------------------------------
  ! 3. Allocate RTTOV input and output structures
  ! --------------------------------------------------------------------------

  ! Determine the total number of radiances to simulate (nchanprof).
  ! In this example we simulate all specified channels for each profile, but
  ! in general one can simulate a different number of channels for each profile.

  nprof_max_sim = nchan_max_sim / nchans_inst
  nprof_sim = MIN(nprof_max_sim, nprofiles)
  nchan_sim = nprof_sim * nchans_inst

  ! Allocate structures for rttov_direct
  call rttov_alloc_direct( &
        errorstatus,            &
        1_jpim,                 &  ! 1 => allocate
        nprof_sim,              &
        nchan_sim,              &
        nlevels,                &
        chanprof,               &
        config_rttov % opts,    &
        profiles,               &
        config_rttov % rttov_coef_array(i_inst),&
        transmission,            &
        radiance,                &
        calcemis=calcemis,       &
        emissivity=emissivity,   &
        init=.TRUE._jplm)

  if (errorstatus /= errorstatus_success) then
    write(message,'(A, I6)') 'after rttov_alloc_direct error = ', errorstatus
    call fckit_log%info(message)
    call rttov_exit(errorstatus)
  end if 
  ! --------------------------------------------------------------------------
  ! 4. Build the list of profile/channel indices in chanprof
  ! --------------------------------------------------------------------------
  prof_start = 1
  prof_end = nprofiles

  do while ( prof_start <= prof_end)
    
    nch = 0_jpim
    do j = 1, MIN(nprof_sim, prof_end - prof_start + 1)
      do jch = 1, nchans_inst
        
        nch = nch + 1_jpim
        chanprof(nch) % prof = j
        chanprof(nch) % chan = channels(jch)
      end do
    end do

    ! --------------------------------------------------------------------------
    ! 6. Specify surface emissivity and reflectance
    ! --------------------------------------------------------------------------

    ! In this example we have no values for input emissivities
    emissivity(:) % emis_in = 0._jprb

    ! Calculate emissivity within RTTOV where the input emissivity value is
    ! zero or less (all channels in this case)
    calcemis(:) = (emissivity(:) % emis_in <= 0._jprb)

    !Assign the data from the GeoVaLs
    !--------------------------------
    if (present(ob_info)) then
      call load_atm_data_rttov(geovals,obss,profiles,prof_start,ob_info=ob_info)
      call load_geom_data_rttov(obss,profiles,prof_start,ob_info=ob_info)
      emissivity(:) % emis_in = ob_info % emiss(:)
      calcemis(:) = ob_info % calc_emiss(:)
    else
      call load_atm_data_rttov(geovals,obss,profiles,prof_start)
      call load_geom_data_rttov(obss,profiles,prof_start)
    end if

    call rttov_user_profile_checkinput(rttov_errorstatus, &
      config_rttov % opts, &
      config_rttov % rttov_coef_array(i_inst), &
      profiles(1))

    call rttov_print_profile(profiles(1))

    ! --------------------------------------------------------------------------
    ! 7. Call RTTOV forward model
    ! --------------------------------------------------------------------------
    call rttov_direct(                         &
      errorstatus,                             &! out   error flag
      chanprof,                                &! in    channel and profile index structure
      config_rttov % opts,                     &! in    options structure
      profiles,                                &! in    profile array
      config_rttov % rttov_coef_array(i_inst), &! in    coefficients structure
      transmission,                            &! inout computed transmittances
      radiance,                                &! inout computed radiances
      calcemis    = calcemis,                  &! in    flag for internal emissivity calcs
      emissivity  = emissivity)!,              &! inout input/output emissivities per channel

    if ( errorstatus /= errorstatus_success ) then
      write(message,'(A, 2I6)') 'after rttov_direct: error ', errorstatus, i_inst
      call fckit_log%info(message)
    end if

    ! Put simulated brightness temperature into hofx
    ! ----------------------------------------------

    hofx(nchans_total+1:nchans_total+nchan_sim) = radiance % bt(1:nchan_sim)

    nchans_total = nchans_total + nchan_sim

    prof_start = prof_start + nprof_sim
  end do

  ! Allocate structures for rttov_direct
  call rttov_alloc_direct(                       &
        errorstatus,                             &
        0_jpim,                                  &  ! 0 => deallocate
        nprof_sim,                              &
        nchan_sim,                               &
        nlevels,                                &
        chanprof,                                &
        config_rttov % opts,                     &
        profiles,                                &
        config_rttov % rttov_coef_array(i_inst), &
        transmission,                            &
        radiance,                                &
        calcemis = calcemis,                     &
        emissivity = emissivity)

  if (errorstatus /= errorstatus_success) then
    write(message,'(A, 2I6)') &
      'after rttov_alloc_direct (deallocation): errorstatus, i_inst =', &
      errorstatus, i_inst
    call fckit_log%info(message)
    call rttov_exit(errorstatus)
  end if

end do Sensor_Loop

end subroutine ufo_radiancerttov_simobs

! ------------------------------------------------------------------------------

end module ufo_radiancerttov_mod
