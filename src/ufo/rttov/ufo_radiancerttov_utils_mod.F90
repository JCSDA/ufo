! (C) Copyright 2018 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

!> Fortran module to provide code shared between nonlinear and tlm/adm radiance calculations

module ufo_radiancerttov_utils_mod

use fckit_configuration_module, only: fckit_configuration
use iso_c_binding
use kinds

use rttov_types, only : rttov_options, rttov_profile, rttov_coefs, &
                        rttov_radiance, rttov_transmission, rttov_emissivity, &
                        rttov_chanprof
use rttov_const, only : errorstatus_success, deg2rad

use ufo_vars_mod
use ufo_geovals_mod, only: ufo_geovals, ufo_geoval, ufo_geovals_get_var
use ufo_basis_mod, only: ufo_basis

implicit none
private

public rttov_conf
public rttov_conf_setup
public rttov_conf_delete
public load_atm_data_rttov
public load_geom_data_rttov

integer, parameter, public :: max_string=800

integer, public  :: rttov_errorstatus

!Type for general config
type rttov_conf
 integer                         :: nsensors
 character(len=255), allocatable :: SENSOR_ID(:)
 character(len=255)              :: COEFFICIENT_PATH
 logical                         :: mw_clw = .FALSE.  ! flag for mw clw
end type rttov_conf

type conf_type_rttov
  type(rttov_coefs), allocatable :: rttov_coef_array(:)
  type(rttov_options)            :: opts
  logical                        :: rttov_is_setup = .FALSE.
  contains
    PROCEDURE :: set_opts => set_options_rttov
    PROCEDURE :: setup =>    setup_rttov
end type conf_type_rttov

type(conf_type_rttov), PUBLIC :: config_rttov

contains

! ------------------------------------------------------------------------------

subroutine rttov_conf_setup(conf, f_conf)

implicit none
type(rttov_conf), intent(inout)       :: conf
type(fckit_configuration), intent(in) :: f_conf

character(len=:),allocatable :: str

!Number of sensors, each call to RTTOV will be for a single sensor
!type (zenith/scan angle will be different)
conf % nSensors = 1

! Allocate SENSOR_ID
allocate(conf % SENSOR_ID(conf % nSensors))

! Get sensor ID from config
call f_conf % get_or_die("Sensor_ID",str)
conf % SENSOR_ID(conf%nSensors) = str

! Path to coefficient files
call f_conf % get_or_die("CoefficientPath",str)
conf % COEFFICIENT_PATH = str

! Flag to turn on mw cloud liquid water (clw)
if (f_conf % has("mw_clw")) then
  call f_conf % get_or_die("mw_clw", conf % mw_clw)
end if

end subroutine rttov_conf_setup

! -----------------------------------------------------------------------------

subroutine rttov_conf_delete(conf)

implicit none
type(rttov_conf), intent(inout) :: conf

deallocate(conf%SENSOR_ID)

end subroutine rttov_conf_delete

! -----------------------------------------------------------------------------

subroutine load_atm_data_rttov(geovals,obss,profiles,prof_start1,obs_info)

use fckit_log_module, only : fckit_log
use obsspace_mod, only : obsspace_get_db, obsspace_get_nlocs, obsspace_has
use ufo_rttovonedvarcheck_utils_mod, only: ObInfo_type

implicit none

type(ufo_geovals), intent(in) :: geovals
type(c_ptr), VALUE, intent(in) :: obss
type(rttov_profile), intent(inout) :: profiles(:)
integer, OPTIONAL, intent(IN) :: prof_start1
type(ObInfo_type), optional, intent(in) :: obs_info  ! Used for rttovonedvarcheck

! Local variables
integer :: k1, nlocs_total, iprof
integer :: nlevels
integer :: nprofiles
integer :: prof_start

type(ufo_geoval), pointer :: geoval
character(MAXVARLEN) :: varname
character(max_string) :: err_msg

real :: ifrac, sfrac, lfrac
real :: itmp, stmp, ltmp
real :: windsp

real(kind_real), allocatable :: TmpVar(:)
real, parameter :: q_mixratio_to_ppmv  = 1.60771704e+3 ! g/kg -> ppmv
character(255) :: message
logical :: variable_present

if(PRESENT(obs_info)) then
  nlocs_total = 1
else
  nlocs_total = obsspace_get_nlocs(obss)
end if

if(PRESENT(prof_start1)) then
  prof_start = prof_start1
else
  prof_start = 1
end if

nprofiles = SIZE(profiles)
 
nlevels = SIZE(profiles(1)%p)

write(message,'(A, 2I6)') ' load_atm_data_rttov  nprofiles nlevels, ', nprofiles, nlevels
call fckit_log%info(message)

!  prof_end = MIN(prof_start + nprofiles, nlocs_total)

do k1 = 1, geovals%nvar
  varname = geovals%variables(k1)
  write(message,'(I6, A)') k1, varname
  call fckit_log%info(message)
end do

!gas_units are ppmv (moist?)

! gas_units = 1 is mixing_ratio (moist)
! gas_units = 2 is ppmv (moist)
! gas_units is per-profile and cannot be set for individual instruments.
profiles(1:nprofiles)%gas_units = 1

varname = "air_pressure" !var_prsi
call ufo_geovals_get_var(geovals, varname, geoval) ! lfric

do iprof = 1, nprofiles
  profiles(iprof)%p(nlevels:1:-1) = geoval%vals(:,prof_start + iprof - 1) / 100! hPa
!  profiles(iprof)%p(1) = 0.5 * profiles(iprof)%p(2)
end do

varname = "air_temperature" !var_ts
call ufo_geovals_get_var(geovals, varname, geoval) !lfric

do iprof = 1, nprofiles
  profiles(iprof)%t(nlevels:1:-1) = geoval%vals(:,prof_start + iprof - 1) ! K
!    profiles(iprof)%t(1) = profiles(iprof)%t(2)
end do

varname = "specific_humidity" !var_mixr
call ufo_geovals_get_var(geovals, varname, geoval)

do iprof = 1, nprofiles
  profiles(iprof)%q(nlevels:1:-1) = geoval%vals(:,prof_start + iprof - 1) ! kg/kg
!    profiles(iprof)%q(1) = profiles(iprof)%q(2)
end do

if (ASSOCIATED(profiles(1)%o3)) then
  varname = var_oz
  call ufo_geovals_get_var(geovals, varname, geoval)
  do iprof = 1, nprofiles
    profiles(iprof)%o3(nlevels:1:-1) = geoval%vals(:,prof_start + iprof - 1) ! kg/kg
!      profiles(iprof)%o3(1) = profiles(iprof)%o3(2)
  end do
end if

if(ASSOCIATED(profiles(1)%co2)) then
  varname = var_co2
  call ufo_geovals_get_var(geovals, varname, geoval)
  do iprof = 1, nprofiles
    profiles(iprof)%co2(nlevels:1:-1) = geoval%vals(:,prof_start + iprof - 1) ! kg/kg
!      profiles(iprof)%co2(1) = profiles(iprof)%co2(2)
  end do
end if

if (ASSOCIATED(profiles(1)%clw)) then
  call ufo_geovals_get_var(geovals, var_clw, geoval)
  do iprof = 1, nprofiles
    profiles(iprof)%clw(nlevels:1:-1) = geoval%vals(:,prof_start + iprof - 1) ! kg/kg
!      profiles(iprof)%clw(1) = profiles(iprof)%clw(2)
  end do
end if

! Near surface

varname = "air_pressure_at_two_meters_above_surface"
call ufo_geovals_get_var(geovals, varname, geoval) ! lfric
profiles(1:nprofiles)%s2m%p = geoval%vals(1,prof_start:prof_start + nprofiles - 1)/ 100

varname = "air_temperature_at_two_meters_above_surface"
call ufo_geovals_get_var(geovals, varname, geoval) ! lfric
profiles(1:nprofiles)%s2m%t = geoval%vals(1,prof_start:prof_start + nprofiles - 1)

varname = "specific_humidity_at_two_meters_above_surface"
call ufo_geovals_get_var(geovals, varname, geoval) ! lfric
profiles(1:nprofiles)%s2m%q = geoval%vals(1,prof_start:prof_start + nprofiles - 1)

!DAR: O3_2m unused
!        profiles(k1)%s2m%o = profiles(k1)%o3(nlevels)

! ! 10m windspeed - I wonder if this will ultimately be U10 and V10 so won't need to be converted - should expect either.
! call ufo_geovals_get_var(geovals, var_sfc_wspeed, geoval)
! windsp = geoval%vals(1,k1)

! call ufo_geovals_get_var(geovals, var_sfc_wdir, geoval)       
! profiles(k1)%s2m%u             = windsp * COS(geoval%vals(1,k1) * deg2rad)
! profiles(k1)%s2m%v             = windsp * SIN(geoval%vals(1,k1) * deg2rad)

varname = "eastward_wind"
call ufo_geovals_get_var(geovals, varname, geoval)
profiles(1:nprofiles)%s2m%u = geoval%vals(1,prof_start:prof_start + nprofiles - 1)

varname = "northward_wind"
call ufo_geovals_get_var(geovals, varname, geoval)
profiles(1:nprofiles)%s2m%v = geoval%vals(1,prof_start:prof_start + nprofiles - 1)

!Skin
varname = "skin_temperature"
call ufo_geovals_get_var(geovals, varname, geoval)
profiles(1:nprofiles)%skin%t = geoval%vals(1,prof_start:prof_start + nprofiles - 1)

varname = "water_type"
call ufo_geovals_get_var(geovals, varname, geoval)
profiles(1:nprofiles)%skin%watertype = geoval%vals(1,prof_start:prof_start + nprofiles - 1)

varname = "surface_type"
call ufo_geovals_get_var(geovals, varname, geoval)
profiles(1:nprofiles)%skin%surftype = geoval%vals(1,prof_start:prof_start + nprofiles - 1)

!DAR: Salinity fixed for now too
profiles(1:nprofiles)%skin%salinity = 35.0

!DAR: Default fastem parameters. We are not using FASTEM over land so these are unused
do k1 = 1,nprofiles
  profiles(k1)%skin%fastem            = [3.0, 5.0, 15.0, 0.1, 0.3]
end do

! !Land point or sea point
! call ufo_geovals_get_var(geovals, var_sfc_wfrac, geoval)
! if (geoval%vals(1,k1) > 0.5) then
!   profiles(k1)%skin%surftype   = 1 ! sea         

!   call ufo_geovals_get_var(geovals, var_sfc_wtmp, geoval)
!   profiles(k1)%skin%t   = geoval%vals(1,k1)

! else ! land
!   profiles(k1)%skin%surftype   = 0 ! land

!   !determine land, snow and ice fractions and temperatures to determine average temperature

!   call ufo_geovals_get_var(geovals, var_sfc_lfrac, geoval)
!   lfrac   = geoval%vals(1,k1)
!   call ufo_geovals_get_var(geovals, var_sfc_sfrac, geoval)
!   sfrac   = geoval%vals(1,k1)
!   call ufo_geovals_get_var(geovals, var_sfc_ifrac, geoval)
!   ifrac   = geoval%vals(1,k1)

!   call ufo_geovals_get_var(geovals, var_sfc_ltmp, geoval)
!   ltmp   = geoval%vals(1,k1)
!   call ufo_geovals_get_var(geovals, var_sfc_stmp, geoval)
!   stmp   = geoval%vals(1,k1)
!   call ufo_geovals_get_var(geovals, var_sfc_itmp, geoval)
!   itmp   = geoval%vals(1,k1)

!   !Skin temperature is a combination of (i)ce temp, (l)and temp and (s)now temp
!   profiles(k1)%skin%t   = (lfrac * ltmp + sfrac * stmp + ifrac * itmp) / (lfrac + sfrac + ifrac)

! end if

!DAR: Could/should get emissivity here?
! call rttov_get_emissivity()

if(PRESENT(obs_info)) then
  profiles(1)%elevation = obs_info%elevation / 1000.0 ! m -> km
  profiles(1)%latitude = obs_info%latitude
  profiles(1)%longitude = obs_info%longitude

else

  allocate(TmpVar(nprofiles))

  variable_present = obsspace_has(obss, "MetaData", "elevation")
  if (variable_present) then
    call obsspace_get_db(obss, "MetaData", "elevation", TmpVar(prof_start:prof_start + nprofiles - 1) )
    profiles(1:nprofiles)%elevation = TmpVar(prof_start:prof_start + nprofiles - 1) / 1000.0 !m -> km for RTTOV
  else
    write(message,'(A)') &
      'MetaData elevation not in database: check implicit filtering'
    call fckit_log%info(message)
  end if

  variable_present = obsspace_has(obss, "MetaData", "latitude")
  if (variable_present) then
    call obsspace_get_db(obss, "MetaData", "latitude", TmpVar(prof_start:prof_start + nprofiles - 1) )
    profiles(1:nprofiles)%latitude = TmpVar(prof_start:prof_start + nprofiles - 1)
  else
    write(message,'(A)') &
      'MetaData latitude not in database: check implicit filtering'
    call fckit_log%info(message)
  end if

  variable_present = obsspace_has(obss, "MetaData", "longitude")
  if (variable_present) then
    call obsspace_get_db(obss, "MetaData", "longitude", TmpVar(prof_start:prof_start + nprofiles - 1) )
    profiles(1:nprofiles)%longitude = TmpVar(prof_start:prof_start + nprofiles - 1)
  else
    write(message,'(A)') &
    'MetaData longitude not in database: check implicit filtering'
    call fckit_log%info(message)
  end if

end if

end subroutine load_atm_data_rttov

! ------------------------------------------------------------------------------

! Internal subprogam to load some test geometry data
subroutine load_geom_data_rttov(obss,profiles,prof_start1,obs_info)

! Satellite viewing geometry
! DAR: check it's all within limits
use obsspace_mod, only :  obsspace_get_nlocs, obsspace_get_db
use ufo_rttovonedvarcheck_utils_mod, only: ObInfo_type

implicit none

type(c_ptr), VALUE,       intent(in)    :: obss
type(rttov_profile), intent(inout) :: profiles(:)
integer, OPTIONAL, intent(IN) :: prof_start1
type(ObInfo_type), optional, intent(in) :: obs_info  ! Used in rttovonedvarcheck

real(kind_real), allocatable :: TmpVar(:)

integer :: prof_start
integer :: nlocs_total, nprofiles
integer :: nlevels

if(PRESENT(prof_start1)) then
  prof_start = prof_start1
else
  prof_start = 1
end if

if(PRESENT(obs_info)) then
  nlocs_total = 1
  nprofiles = 1
  nlevels = SIZE(profiles(1)%p)
  
  profiles(1)%zenangle    = obs_info%sensor_zenith_angle
  profiles(1)%azangle     = obs_info%sensor_azimuth_angle
  profiles(1)%sunzenangle = obs_info%solar_zenith_angle
  profiles(1)%sunazangle  = obs_info%solar_azimuth_angle
  
else

  nlocs_total = obsspace_get_nlocs(obss)
  nprofiles = SIZE(profiles)

  allocate(TmpVar(nprofiles))

  nlevels = SIZE(profiles(1)%p)

  call obsspace_get_db(obss, "MetaData", "sensor_zenith_angle", TmpVar(prof_start:prof_start + nprofiles - 1))
  profiles(1:nprofiles)%zenangle = TmpVar(prof_start:prof_start + nprofiles - 1)

  call obsspace_get_db(obss, "MetaData", "sensor_azimuth_angle", TmpVar(prof_start:prof_start + nprofiles - 1))
  profiles(1:nprofiles)%azangle = TmpVar(prof_start:prof_start + nprofiles - 1)

  call obsspace_get_db(obss, "MetaData", "solar_zenith_angle", TmpVar(prof_start:prof_start + nprofiles - 1))
  profiles(1:nprofiles)%sunzenangle = TmpVar(prof_start:prof_start + nprofiles - 1)

  call obsspace_get_db(obss, "MetaData", "solar_azimuth_angle", TmpVar(prof_start:prof_start + nprofiles - 1))
  profiles(1:nprofiles)%sunazangle = TmpVar(prof_start:prof_start + nprofiles - 1)

  deallocate(TmpVar)

end if

end subroutine load_geom_data_rttov

! ------------------------------------------------------------------------------

subroutine set_options_rttov(self, conf)
implicit none
class(conf_type_rttov), intent(INOUT) :: self
type(rttov_conf), intent(in) :: conf

self % opts % rt_ir % addsolar            = .FALSE. ! Do not include solar radiation
self % opts % interpolation % addinterp   = .TRUE.  ! Allow interpolation of input profile
self % opts % interpolation % interp_mode = 1       ! Set interpolation method
self % opts % interpolation % reg_limit_extrap = .TRUE. ! Set interpolation methodreg_limit_extrap
self % opts % rt_all % addrefrac          = .TRUE.  ! Include refraction in path calc
self % opts % rt_all % switchrad          = .TRUE.  ! Include refraction in path calc
self % opts % rt_ir % addclouds           = .FALSE. ! Don't include cloud effects
self % opts % rt_ir % addaerosl           = .FALSE. ! Don't include aerosol effects

self % opts % rt_ir % ozone_data          = .FALSE. ! Set the relevant flag to .TRUE.
self % opts % rt_ir % co2_data            = .FALSE. !   when supplying a profile of the
self % opts % rt_ir % n2o_data            = .FALSE. !   given trace gas (ensure the
self % opts % rt_ir % ch4_data            = .FALSE. !   coef file supports the gas)
self % opts % rt_ir % co_data             = .FALSE. !
self % opts % rt_ir % so2_data            = .FALSE. !
self % opts % rt_mw % clw_data            = .FALSE. !

self % opts % config % verbose            = .TRUE.  ! Enable printing of warnings
self % opts % config % apply_reg_limits   = .TRUE.
self % opts % config % do_checkinput      = .TRUE.

! Update based on rttov conf input
if (conf % mw_clw) self % opts % rt_mw % clw_data = .TRUE.

end subroutine set_options_rttov

! ------------------------------------------------------------------------------

subroutine setup_rttov(self, conf, asw)
class(conf_type_rttov) :: self
type(rttov_conf), intent(in) :: conf
integer, intent(IN) :: asw !allocate switch

character(len=255) :: coef_filename
integer :: i_inst

INCLUDE 'rttov_read_coefs.interface'

rttov_errorstatus = 0

if(asw == 1) then

  ! --------------------------------------------------------------------------
  ! 1. Setup rttov options
  ! --------------------------------------------------------------------------
  call self % set_opts(conf)

  ! --------------------------------------------------------------------------
  ! 2. Read coefficients
  ! --------------------------------------------------------------------------
  allocate(self%rttov_coef_array(conf%nSensors))

  do i_inst = 1, conf%nSensors
    coef_filename = &
      TRIM(conf%COEFFICIENT_PATH)//'rtcoef_'//TRIM(conf%SENSOR_ID(i_inst))//'.dat'
    call rttov_read_coefs(rttov_errorstatus, &             !out
                          self%rttov_coef_array(i_inst), & !inout
                          self%opts, &                     !in
                          file_coef=coef_filename)         !in

    if (rttov_errorstatus /= errorstatus_success) then
      WRITE(*,*) 'fatal error reading coefficients'
      !          call rttov_exit(errorstatus)
    else
      WRITE(*,*) 'successfully read' // coef_filename
    end if

  end do

  self%rttov_is_setup =.TRUE.

else !asw == 0
  !Quick and dirty for now
  deallocate(self%rttov_coef_array)
  self%rttov_is_setup =.FALSE.
end if

end subroutine setup_rttov

! ------------------------------------------------------------------------------

subroutine get_var_name(varname_tmplate,n,varname)

character(len=*), intent(in) :: varname_tmplate
integer, intent(in) :: n
character(len=*), intent(out) :: varname

character(len=3) :: chan

! pass in varname_tmplate = "brightness_temperature"
WRITE(chan, '(I0)') n
varname = TRIM(varname_tmplate) // '_' // TRIM(chan) // '_'

end subroutine get_var_name

! -----------------------------------------------------------------------------

subroutine get_var_name_new(varname_tmplate,n,varname)

character(len=*), intent(in) :: varname_tmplate
integer, intent(in) :: n
character(len=*), intent(out) :: varname

character(len=3) :: chan

 ! pass in varname_tmplate = "brigtness_temperature"
 write(chan, '(I0)') n
 varname = trim(varname_tmplate) // '_' // trim(chan)

end subroutine get_var_name_new

! ------------------------------------------------------------------------------

end module ufo_radiancerttov_utils_mod
