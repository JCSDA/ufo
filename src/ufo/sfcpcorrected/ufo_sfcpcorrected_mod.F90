! (C) Copyright 2017-2018 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

!> Fortran module for sfcpcorrected observation operator

module ufo_sfcpcorrected_mod

 use oops_variables_mod
 use ufo_vars_mod
 use missing_values_mod
 use iso_c_binding
 use kinds
 use ufo_constants_mod, only : grav, rd, Lclr, t2tv

 implicit none
 private
 integer, parameter :: max_string = 800

!> Fortran derived type for the observation type
 type, public :: ufo_sfcpcorrected
 private
   type(oops_variables), public :: obsvars
   type(oops_variables), public :: geovars
   character(len=MAXVARLEN)     :: da_psfc_scheme
 contains
   procedure :: setup  => ufo_sfcpcorrected_setup
   procedure :: simobs => ufo_sfcpcorrected_simobs
 end type ufo_sfcpcorrected

 character(len=MAXVARLEN), dimension(5), parameter :: geovars_default = (/ var_ps, var_geomz, var_sfc_geomz, var_tv, var_prs /)

contains

! ------------------------------------------------------------------------------
subroutine ufo_sfcpcorrected_setup(self, f_conf)
use fckit_configuration_module, only: fckit_configuration
implicit none
class(ufo_sfcpcorrected), intent(inout)     :: self
type(fckit_configuration), intent(in) :: f_conf
character(len=:), allocatable         :: str

call self%geovars%push_back(geovars_default)

self%da_psfc_scheme = "UKMO"
if (f_conf%has("da_psfc_scheme")) then
   call f_conf%get_or_die("da_psfc_scheme",str)
   self%da_psfc_scheme = str
end if

end subroutine ufo_sfcpcorrected_setup

! ------------------------------------------------------------------------------
subroutine ufo_sfcpcorrected_simobs(self, geovals, obss, nvars, nlocs, hofx)
use ufo_geovals_mod, only: ufo_geovals, ufo_geoval, ufo_geovals_get_var
use obsspace_mod
use vert_interp_mod
use fckit_log_module,  only : fckit_log
implicit none
class(ufo_sfcpcorrected), intent(in)    :: self
integer, intent(in)               :: nvars, nlocs
type(ufo_geovals),  intent(in)    :: geovals
real(c_double),     intent(inout) :: hofx(nvars, nlocs)
type(c_ptr), value, intent(in)    :: obss

! Local variables
real(c_double)                    :: missing
real(kind_real)                   :: H2000 = 2000.0
integer                           :: nobs, iobs, ivar
real(kind_real),    allocatable   :: cor_psfc(:)
type(ufo_geoval),   pointer       :: model_ps, model_p, model_sfc_geomz, model_tv, model_geomz
character(len=*), parameter       :: myname_="ufo_sfcpcorrected_simobs"
character(max_string)             :: err_msg
real(kind_real)                   :: wf
integer                           :: wi
logical                           :: variable_present
real(kind_real), dimension(:), allocatable :: obs_height, obs_t, obs_q, obs_psfc
real(kind_real), dimension(:), allocatable :: model_tvs, model_zs, model_level1, model_p_2000, model_tv_2000, model_psfc

missing = missing_value(missing)
nobs    = obsspace_get_nlocs(obss)

! check if nobs is consistent in geovals & nlocs
if (geovals%nlocs /= nobs) then
   write(err_msg,*) myname_, ' error: nlocs of model and obs is inconsistent!'
   call abor1_ftn(err_msg)
endif

! cor_psfc: observed surface pressure at model surface height, corresponding to P_o2m in da_intpsfc_prs* subroutines
allocate(cor_psfc(nobs))

! get obs variables
allocate(obs_height(nobs))
allocate(obs_psfc(nobs))
call obsspace_get_db(obss, "MetaData",  "station_elevation",obs_height)
call obsspace_get_db(obss, "ObsValue",  "surface_pressure", obs_psfc)

! get model variables
call ufo_geovals_get_var(geovals, var_ps, model_ps)
call ufo_geovals_get_var(geovals, var_geomz, model_geomz)
call ufo_geovals_get_var(geovals, var_sfc_geomz, model_sfc_geomz)
call ufo_geovals_get_var(geovals, var_tv, model_tv)
call ufo_geovals_get_var(geovals, var_prs, model_p)

if (model_geomz%vals(1,1) .gt. model_geomz%vals(model_geomz%nval,1) ) then
   write(err_msg,'(a)') '  ufo_sfcpcorrected:'//new_line('a')//                   &
                        '  Model vertical height profile is from top to bottom'
   call fckit_log%info(err_msg)
end if

allocate(model_zs(nobs))
allocate(model_level1(nobs))
allocate(model_psfc(nobs))

model_zs = model_sfc_geomz%vals(1,:)
model_level1 = model_geomz%vals(model_geomz%nval,:)   !reverse
model_psfc = model_ps%vals(1,:)

! do terrain height correction, two optional schemes
select case (trim(self%da_psfc_scheme))

case ("WRFDA")
   ! get extra obs values
   variable_present = obsspace_has(obss, "ObsValue", "air_temperature")
   if (variable_present) then
      allocate(obs_t(nobs))
      call obsspace_get_db(obss, "ObsValue", "air_temperature", obs_t)
   end if
   variable_present = obsspace_has(obss, "ObsValue", "specific_humidity")
   if (variable_present) then
      allocate(obs_q(nobs))
      call obsspace_get_db(obss, "ObsValue", "specific_humidity", obs_q)
   end if

   ! get extra model values
   allocate(model_tvs(nobs))
   model_tvs = model_tv%vals(model_tv%nval,:) + Lclr * ( model_level1 - model_zs )  !Lclr = 0.0065 K/m

   ! correction
   call da_intpsfc_prs(nobs, missing, cor_psfc, obs_height, obs_psfc, model_zs, model_tvs, obs_t, obs_q)

   deallocate(obs_t)
   deallocate(obs_q)
   deallocate(model_tvs)

case ("UKMO")

   allocate(model_p_2000(nobs))
   allocate(model_tv_2000(nobs))
   do iobs = 1, nobs
      ! vertical interpolation for getting model P and tv at 2000 m
      call vert_interp_weights(model_geomz%nval, H2000, model_geomz%vals(:,iobs), wi, wf)
      call vert_interp_apply(model_p%nval, model_p%vals(:,iobs), model_p_2000(iobs), wi, wf)
      call vert_interp_apply(model_tv%nval, model_tv%vals(:,iobs), model_tv_2000(iobs), wi, wf)
   end do

   ! correction
   call da_intpsfc_prs_ukmo(nobs, missing, cor_psfc, obs_height, obs_psfc, model_zs, model_psfc, model_tv_2000, model_p_2000)

   deallocate(model_p_2000)
   deallocate(model_tv_2000)

case default
   write(err_msg,*) "ufo_sfcpcorrected_mod.F90: da_psfc_scheme must be WRFDA or UKMO"
   call fckit_log%info(err_msg)
   call abor1_ftn(err_msg)
end select

! update the obs surface pressure
do ivar = 1, nvars
   do iobs = 1, nlocs
      if ( cor_psfc(iobs) /= missing) then
         hofx(ivar,iobs) = obs_psfc(iobs) - cor_psfc(iobs) + model_psfc(iobs)
      else
         hofx(ivar,iobs) = model_psfc(iobs)
      end if
   enddo
enddo

deallocate(obs_height)
deallocate(obs_psfc)

deallocate(model_zs)
deallocate(model_level1)
deallocate(model_psfc)

end subroutine ufo_sfcpcorrected_simobs

! ------------------------------------------------------------------------------
!> \Conduct terrain height correction for surface pressure
!!
!! \This subroutine is based on a subroutine from WRFDA da_intpsfc_prs.inc file
!!  corresponding to sfc_assi_options = 1 in WRFDA's namelist
!!
!! \Date: June 2019: Created
!! \Method: hydrosatic equation
!!
!!  P_o2m = P_o * exp [-grav/rd * (H_m-H_o) / (TV_m + TV_o)/2)
!!
!!  Where:
!!  H_m    = model surface height
!!  H_o    = station height
!!  TV_m   = virtual temperature at model surface height
!!  TV_o   = virtual temperature at station height
!!  P_o2m  = pressure interpolated from station height to model surface height
!!  P_o    = pressure at station height
!!  grav   = gravitational acceleration
!!  rd     = gas constant per mole
!!

subroutine da_intpsfc_prs (nobs, missing, P_o2m, H_o, P_o, H_m, TV_m, T_o, Q_o)
implicit none
integer,                          intent (in)           :: nobs !<total observation number
real(c_double),                   intent (in)           :: missing
real(kind_real), dimension(nobs), intent (out)          :: P_o2m !<observed PS at model sfc height
real(kind_real), dimension(nobs), intent (in)           :: H_o, P_o !<observed Height and PS
real(kind_real), dimension(nobs), intent (in)           :: H_m, TV_m !<model sfc height and TV
real(kind_real), dimension(nobs), intent (in), optional :: T_o, Q_o !<obserbed T and Q
real(kind_real), dimension(nobs)                        :: TV_o, TV

! 1.  model and observation virtual temperature
! ---------------------------------------------

if (present(T_o) .and. present(Q_o)) then
   TV_o = T_o  * (1.0 + t2tv * Q_o)
else if (present(T_o) .and. .not.(present(Q_o))) then
      TV_o = T_o
else
      TV_o = TV_m
end if

TV  = 0.5 * (TV_m + TV_o)

! 2.  extrapolate pressure from station height to model surface height
! --------------------------------------------------------------------

where ( H_o /= missing .and. P_o /= missing )
   P_o2m = P_o * exp ( - (H_m - H_o) * grav / (rd * TV) )
elsewhere
   P_o2m = P_o
end where

end subroutine da_intpsfc_prs

! ------------------------------------------------------------------------------
!> \Conduct terrain height correction for surface pressure
!!
!! \Reference: Ingleby,2013. UKMO Technical Report No: 582. Appendix 1.
!!
!! \Method: integrate the hydrosatic equation dp/dz=-rho*g/RT to get P_m2o first, equation:
!!
!!  (P_m2o/P_m)=(T_m2o/T_m)** (grav/rd*L)
!!
!!  Where:
!!  P_m2o = model surface pressure at station height
!!  P_m   = model surface pressure
!!  T_m   = temperature at model surface height; derived from TV_2000
!!  T_m2o = model surface temperature at station height
!!  grav  = gravitational acceleration
!!  rd    = gas constant per mole
!!  Lclr  = constant lapse rate (0.0065 K/m)
!!
!!  To avoid dirunal/local variations, use TV_2000 (2000 m above the model surface height) instead of direct T_m
!!
!!  T_m = TV_2000 * (P_o / P_2000) ** (rd*L/grav)
!!
!! Where:
!!  P_2000  = background pressure at 2000 m
!!  TV_2000 = background virtual temperature at 2000 m
!!  P_o     = pressure at station height
!!
!!  Finally, in practice, adjust P_o to the model surface height using
!!
!!  P_o2m = P_o * (P_m / P_m2o)
!!

subroutine da_intpsfc_prs_ukmo (nobs, missing, P_o2m, H_o, P_o, H_m, P_m, TV_2000, P_2000)
implicit none
integer,                          intent (in)           :: nobs !<total observation number
real(c_double),                   intent (in)           :: missing
real(kind_real), dimension(nobs), intent (out)          :: P_o2m !<observed PS at model sfc height
real(kind_real), dimension(nobs), intent (in)           :: H_o, P_o !<observed Height and PS
real(kind_real), dimension(nobs), intent (in)           :: H_m, P_m, TV_2000, P_2000 !<model Height, PS, TV at 2000 m, and P at 2000 m
real(kind_real)                                         :: ind !<local variable
real(kind_real), dimension(nobs)                        :: P_m2o, T_m, T_m2o !<local variables: model PS at observed height, model sfc T,
                                                                             !!             and model T  at observed height

! define the constant power exponent
ind = rd * Lclr / grav

where ( H_o /= missing .and. P_o /= missing )
   ! calculate T_m   -- background temperature at model surface height
   !           T_m2o -- background temperature at station height
   T_m = TV_2000 * (P_m / P_2000) ** ind
   T_m2o = T_m + Lclr * ( H_m - H_o)
   P_m2o = P_m * (T_m2o / T_m) ** (1.0 / ind)

   ! In practice, P_o is adjusted to the model surface height
   P_o2m = P_o * P_m / P_m2o
elsewhere
   P_o2m = P_o
end where

end subroutine da_intpsfc_prs_ukmo

! ------------------------------------------------------------------------------

end module ufo_sfcpcorrected_mod
