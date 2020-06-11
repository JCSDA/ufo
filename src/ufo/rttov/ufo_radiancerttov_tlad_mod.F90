! (C) Copyright 2017-2018 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

!> Fortran module for radiancerttov tl/ad observation operator

module ufo_radiancerttov_tlad_mod

use fckit_configuration_module, only: fckit_configuration
use iso_c_binding
use kinds

use ufo_geovals_mod, only: ufo_geovals, ufo_geoval, ufo_geovals_get_var
use ufo_basis_tlad_mod, only: ufo_basis_tlad
use ufo_vars_mod
use ufo_radiancerttov_utils_mod

use rttov_types
use rttov_const, only : errorstatus_success, deg2rad
use rttov_unix_env

implicit none
private

!> Fortran derived type for radiancerttov trajectory
type, public :: ufo_radiancerttov_tlad
private
integer :: nprofiles
! integer :: nlayers
integer :: nlevels
integer :: nchannels

type(rttov_conf) :: conf
type(rttov_profile), pointer, public :: profiles_k(:) => NULL()
type(rttov_chanprof), pointer :: chanprof(:) => NULL()
logical :: ltraj
contains
  procedure :: setup  => ufo_radiancerttov_tlad_setup
  procedure :: delete  => ufo_radiancerttov_tlad_delete
  procedure :: settraj => ufo_radiancerttov_tlad_settraj
  procedure :: simobs_tl  => ufo_radiancerttov_simobs_tl
  procedure :: simobs_ad  => ufo_radiancerttov_simobs_ad
end type ufo_radiancerttov_tlad

contains

! ------------------------------------------------------------------------------
subroutine ufo_radiancerttov_tlad_setup(self, f_conf)
implicit none
class(ufo_radiancerttov_tlad), intent(inout) :: self
type(fckit_configuration),     intent(in)    :: f_conf

! Local variables
type(fckit_configuration) :: f_confOpts

! Setup conf
call f_conf % get_or_die("ObsOptions", f_confOpts)
call rttov_conf_setup(self % conf, f_confOpts)

end subroutine ufo_radiancerttov_tlad_setup

! ------------------------------------------------------------------------------
subroutine ufo_radiancerttov_tlad_delete(self)
  use ufo_radiancerttov_utils_mod , ONLY : config_rttov

  implicit none
  class(ufo_radiancerttov_tlad), intent(inout) :: self

  integer(kind=jpim) :: errorstatus

  include 'rttov_alloc_prof.interface'

  if (ASSOCIATED(self % profiles_k)) then
    call rttov_alloc_prof(errorstatus, & 
                               SIZE(self % profiles_k), &
                               self % profiles_k, &
                               -1, & ! doesn't matter what nlevels is
                               config_rttov % opts, & 
                               asw = 0)!, &
    deallocate(self % profiles_k)
    deallocate(self % chanprof)
    WRITE(*,*) 'Killing profiles_k and chanprof'
    nullify( self % profiles_k, self % chanprof)
  end if

  self % ltraj = .FALSE.

  call rttov_conf_delete(self % conf)

end subroutine ufo_radiancerttov_tlad_delete

! ------------------------------------------------------------------------------
subroutine ufo_radiancerttov_tlad_settraj(self, geovals, obss, channels, ob_info, BT)

use ufo_radiancerttov_utils_mod , ONLY : config_rttov
use ufo_rttovonedvarcheck_obinfo_mod, ONLY : ObInfo_type

implicit none

class(ufo_radiancerttov_tlad), intent(inout) :: self
type(ufo_geovals),             intent(in)    :: geovals
type(c_ptr), value,            intent(in)    :: obss
integer(c_int),                intent(in)    :: channels(:) ! List of channels to use
type(ObInfo_type), optional,   intent(in)    :: ob_info    ! Used for rttovonedvarcheck
real(kind_real), optional,     intent(out)   :: BT(:)       ! Used for rttovonedvarcheck

! Local Variables
character(*), parameter                      :: PROGRAM_NAME = 'ufo_radiancerttov_tlad_settraj'
character(255)                               :: message, version
integer                                      :: err_stat, alloc_stat
integer                                      :: i_inst,j , jch,nlevels, nch, nchannels, nchans_total, nchans_inst, asw, ierr
type(ufo_geoval), pointer                    :: temp

logical(kind=jplm),      pointer             :: calcemis(:)    => NULL()   ! Flag to indicate calculation of emissivity within RTTOV

type(rttov_emissivity),  pointer             :: emissivity(:)  => NULL()   ! Input/output surface emissivity
type(rttov_profile),     pointer             :: profiles(:)    => NULL()   ! Input profiles
type(rttov_profile),     pointer             :: profiles_k(:)  => NULL()   ! Input profiles
type(rttov_chanprof),    pointer             :: chanprof(:)    => NULL()   ! Input profiles
type(rttov_transmission)                     :: transmission               ! Output transmittances
type(rttov_radiance)                         :: radiance                   ! Output radiances

type(rttov_emissivity),  pointer             :: emissivity_k(:)  => NULL() ! Input/output surface emissivity
type(rttov_transmission)                     :: transmission_k             ! Output transmittances
type(rttov_radiance)                         :: radiance_k                 ! Output radiances

integer(kind=jpim)                           :: errorstatus                ! Return error status of RTTOV subroutine calls

integer                                      :: nprof_sim, nchan_sim, nchan_max_sim, nprof_max_sim, nchan_count
integer                                      :: prof_start, prof_end

character(MAXVARLEN)                         :: varname

integer :: idbg = 1

character(MAXVARLEN) :: label

include 'rttov_k.interface'
include 'rttov_alloc_k.interface'
include 'rttov_alloc_prof.interface'
include 'rttov_print_profile.interface'
include 'rttov_user_profile_checkinput.interface'

! Get number of profile and layers from geovals
! ---------------------------------------------
self % nprofiles = geovals % nlocs
varname = 'air_pressure' ! var_prsi
call ufo_geovals_get_var(geovals, varname, temp)
self % nlevels = temp % nval
! nlayers = self % nlevels - 1

nullify(temp)

errorstatus = 0_jpim
nchan_count = 0

asw = 1

nchan_max_sim = 2500 ! Maximum number of channels to pass to RTTOV to simulate

if( .NOT. config_rttov % rttov_is_setup) then
  call config_rttov % setup(self % conf, asw)
else
  write(*,*) "Config rttov TLAD already setup"
end if

Sensor_Loop:do i_inst = 1, self % conf % nSensors

!  nchans_inst = config_rttov % rttov_coef_array(i_inst) % coef % fmv_chn
  nchans_inst = SIZE(channels)
  self % nchannels = nchans_inst
  nchans_total = self % nchannels * self % nprofiles

  if( .NOT. ASSOCIATED(self % profiles_k)) then
    ! one channel will be ~ 3-15 fields * 8 bytes * ~100 levels = approx 10kB
    ! so 15 million channels * 1e4 bytes / channel will be 1.5e11 bytes (150 GB - too big!!)
    allocate(self % chanprof(nchans_total))
    allocate(self % profiles_k(nchans_total))

    call rttov_alloc_prof(errorstatus, &
                             nchans_total, &
                             self % profiles_k, &
                             self % nlevels, &
                             config_rttov % opts, & 
                             asw,&!, &
!                            coefs = config_rttov % rttov_coef_array(i_inst) ! cld/aer only
                             init = .TRUE.)
  end if

  ! Ensure the options and coefficients are consistent
  call rttov_user_options_checkinput(errorstatus, &
                                     config_rttov % opts, &
                                     config_rttov % rttov_coef_array(i_inst))
  if (errorstatus /= errorstatus_success) then
    WRITE(*,*) 'error in rttov options'
    call rttov_exit(errorstatus)
  end if

  ! --------------------------------------------------------------------------
  ! 3. Allocate RTTOV input and output structures
  ! --------------------------------------------------------------------------

  ! Determine the total number of radiances to simulate (nchanprof).
  ! In this example we simulate all specified channels for each profile, but
  ! in general one can simulate a different number of channels for each profile.

  nprof_max_sim = nchan_max_sim / nchans_inst
  nprof_sim = MIN(nprof_max_sim, self % nprofiles)
  nchan_sim = nprof_sim * nchans_inst

  idbg = idbg + 1

  ! Allocate temporary structures for rttov_k
  call rttov_alloc_k(                            &
        errorstatus,                             &
        1_jpim,                                  &  ! 1 => allocate
        nprof_sim,                               &
        nchan_sim,                               &
        self % NLEVELS,                          &
        chanprof,                                &
        config_rttov % opts,                     &
        profiles,                                &
        profiles_k,                              &
        config_rttov % rttov_coef_array(i_inst), &
        transmission,                            &
        transmission_k,                          &
        radiance,                                &
        radiance_k,                              &
        calcemis=calcemis,                       &
        emissivity=emissivity,                   &
        emissivity_k=emissivity_k,               &
        init=.TRUE._jplm)

  if (errorstatus /= errorstatus_success) then
    WRITE(*,*) 'allocation error for rttov_k structures'
    call rttov_exit(errorstatus)
  end if

  idbg = idbg + 1

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

  prof_start = 1
  prof_end = self % nprofiles

  idbg = idbg + 1

  do while ( prof_start <= prof_end)
    nch = 0_jpim
    do j = 1, MIN(nprof_sim, prof_end - prof_start + 1)
      do jch = 1, nchans_inst
        nch = nch + 1_jpim
        self % chanprof(nchan_count + nch) % prof = j
        self % chanprof(nchan_count + nch) % chan = channels(jch) 
        ! only all channels for now. Look at OPS for better implementation.
        !local
        chanprof(nch) % prof = j
        chanprof(nch) % chan = channels(jch)
      end do
    end do

    idbg = idbg + 1

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

    idbg = idbg + 1

    ! --------------------------------------------------------------------------
    ! 7. Call RTTOV forward model
    ! --------------------------------------------------------------------------
    call rttov_k(                              &
      errorstatus,                             &! out   error flag
      self % chanprof(nchan_count + 1:nchan_count + nch), &! in LOCAL channel and profile index structure
      config_rttov % opts,                     &! in    options structure
      profiles,                                &! in    profile array
      self % profiles_k(nchan_count + 1 : nchan_count + nch), &! inout    profile array
      config_rttov % rttov_coef_array(i_inst), &! in    coefficients structure
      transmission,                            &! inout computed transmittances
      transmission_k,                          &! inout computed transmittances
      radiance,                                &! inout computed radiances
      radiance_k,                              &! inout computed radiances
      calcemis    = calcemis,                  &! in    flag for internal emissivity calcs
      emissivity  = emissivity,                &!, &! inout input/output emissivities per channel
      emissivity_k = emissivity_k)!,           &! inout input/output emissivities per channel

    idbg = idbg + 1

    if ( errorstatus /= errorstatus_success ) then
      message = 'Error calling RTTOV K Model for amsua'!//TRIM(SENSOR_ID(n))
      WRITE(*,*) message
!     STOP
    end if

    prof_start = prof_start + nprof_sim
    nchan_count = nchan_count + nch
  end do
  
  if (present(BT)) BT(:) = radiance % bt(:)

  ! Allocate structures for rttov_k
  call rttov_alloc_k(                        &
    errorstatus,                             &
    0_jpim,                                  &  ! 1 => allocate, 0=> deallocate
    nprof_sim,                               &
    nchan_sim,                               &
    self % NLEVELS,                          &
    chanprof,                                &
    config_rttov % opts,                     &
    profiles,                                &
    profiles_k,                              &
    config_rttov % rttov_coef_array(i_inst), &
    transmission,                            &
    transmission_k,                          &
    radiance,                                &
    radiance_k,                              &
    calcemis=calcemis,                       &
    emissivity=emissivity,                   &
    emissivity_k=emissivity_k,               &
    init=.TRUE._jplm)

    idbg = idbg + 1

end do Sensor_Loop

! Set flag that the tracectory was set
! ------------------------------------
self % ltraj = .true.

end subroutine ufo_radiancerttov_tlad_settraj

! ------------------------------------------------------------------------------
subroutine ufo_radiancerttov_simobs_tl(self, geovals, obss, hofx, channels)
implicit none
class(ufo_radiancerttov_tlad), intent(in) :: self
type(ufo_geovals), intent(in) :: geovals
type(c_ptr), value, intent(in) :: obss
real(c_double), intent(inout) :: hofx(:)
integer(c_int), intent(in) :: channels(:)  !List of channels to use

character(len=*), parameter :: myname_="ufo_radiancerttov_simobs_tl"
character(max_string) :: err_msg
integer :: job, jprofile, jchannel, jlevel, ierr, ichan, prof, nlevels, lev
type(ufo_geoval), pointer :: tv_d

character(MAXVARLEN) :: varname

! Initial checks
! --------------

! Check if trajectory was set
if (.not. self % ltraj) then
  write(err_msg,*) myname_, ' trajectory wasnt set!'
  call abor1_ftn(err_msg)
end if

! Check if nlocs is consistent in geovals & hofx
if (geovals % nlocs /= self % nprofiles) then
  write(err_msg,*) myname_, ' error: nlocs inconsistent!'
  call abor1_ftn(err_msg)
end if

! Initialize hofx
! ---------------
hofx(:) = 0.0_kind_real

! Temperature
! -----------
call ufo_geovals_get_var(geovals, var_ts, tv_d) ! var_ts = air_temperature

! Check model levels is consistent in geovals
if (tv_d % nval /= self % nlevels) then
  write(err_msg,*) myname_, ' error: layers inconsistent!'
  call abor1_ftn(err_msg)
end if

nlevels = SIZE(self % profiles_k(1) % t)

do ichan = 1, self % nprofiles * self % nchannels
  prof = self % chanprof(ichan) % prof

  hofx(ichan) = hofx(ichan) + &
    SUM(self % profiles_k(ichan) % t(nlevels:1:-1) * tv_d % vals(:,prof))

end do

end subroutine ufo_radiancerttov_simobs_tl

! ------------------------------------------------------------------------------
subroutine ufo_radiancerttov_simobs_ad(self, geovals, obss, hofx, channels)
implicit none
class(ufo_radiancerttov_tlad), intent(in)    :: self
type(ufo_geovals),       intent(inout) :: geovals
type(c_ptr), value,      intent(in)    :: obss
real(c_double),          intent(in)    :: hofx(:)
integer(c_int),           intent(in)    :: channels(:)  !List of channels to use

character(len=*), parameter :: myname_="ufo_radiancerttov_simobs_ad"
character(max_string) :: err_msg
integer :: job, jprofile, jchannel, jlevel, ierr, ichan, prof, nlevels, lev
type(ufo_geoval), pointer :: tv_d

character(MAXVARLEN) :: varname

! Initial checks
! --------------

! Check if trajectory was set
if (.not. self % ltraj) then
  write(err_msg,*) myname_, ' trajectory wasnt set!'
  call abor1_ftn(err_msg)
end if

! Check if nlocs is consistent in geovals & hofx
if (geovals % nlocs /= self % nprofiles) then
  write(err_msg,*) myname_, ' error: nlocs inconsistent!'
  call abor1_ftn(err_msg)
end if


! Temperature
! -----------
call ufo_geovals_get_var(geovals, var_ts, tv_d) ! var_ts = air_temperature


! allocate if not yet allocated
if (.not. allocated(tv_d % vals)) then
  tv_d % nlocs = self % nprofiles
  tv_d % nval = self % nlevels
  allocate(tv_d % vals(tv_d % nval,tv_d % nlocs))
  tv_d % vals = 0.0_kind_real
end if

nlevels = SIZE(self % profiles_k(1) % t)


do ichan = 1, self % nprofiles * self % nchannels
  prof = self % chanprof(ichan) % prof

  tv_d % vals(:,prof) = tv_d % vals(:,prof) + &
    self % profiles_k(ichan) % t(nlevels:1:-1) * hofx(ichan)
end do

! Once all geovals set replace flag
! ---------------------------------
if (.not. geovals % linit ) geovals % linit=.true.

end subroutine ufo_radiancerttov_simobs_ad

! ------------------------------------------------------------------------------

end module ufo_radiancerttov_tlad_mod
