! (C) Copyright 2017-2018 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

!> Fortran stericheight module for observation operator

module ufo_stericheight_mod

 use iso_c_binding
 use config_mod
 use kinds

 use ufo_geovals_mod, only: ufo_geovals, ufo_geoval, ufo_geovals_get_var
 use ufo_basis_mod, only: ufo_basis
 use ufo_vars_mod
 use obsspace_mod

 implicit none
 private

 integer, parameter :: max_string=800

!> Fortran derived type for the observation type
 type, extends(ufo_basis), public :: ufo_stericheight
 private
 contains
   procedure :: setup  => ufo_stericheight_setup
   procedure :: delete => ufo_stericheight_delete
   procedure :: simobs => ufo_stericheight_simobs
 end type ufo_stericheight

contains

! ------------------------------------------------------------------------------
subroutine ufo_stericheight_setup(self, c_conf)
implicit none
class(ufo_stericheight), intent(inout) :: self
type(c_ptr),        intent(in)    :: c_conf

end subroutine ufo_stericheight_setup

! ------------------------------------------------------------------------------
subroutine ufo_stericheight_delete(self)
implicit none
class(ufo_stericheight), intent(inout) :: self

end subroutine ufo_stericheight_delete

! ------------------------------------------------------------------------------
subroutine ufo_stericheight_simobs(self, geovals, hofx, obss)
implicit none
class(ufo_stericheight), intent(in)    :: self
type(ufo_geovals),  intent(in)    :: geovals
real(c_double),     intent(inout) :: hofx(:)
type(c_ptr), value, intent(in)    :: obss

character(len=*), parameter :: myname_="ufo_stericheight_simobs"
character(max_string) :: err_msg

integer :: iobs
type(ufo_geoval), pointer :: geoval_temp, geoval_salt, geoval_adt

print *,myname_

!dh - this cant be done here *** hofx%nobs = geovals%nobs
print *, myname_, ' nobs: ', geovals%nobs, size(hofx,1)

! check if nobs is consistent in geovals & hofx
if (geovals%nobs /= size(hofx,1)) then
  write(err_msg,*) myname_, ' error: nobs inconsistent!'
  call abor1_ftn(err_msg)
endif

! Sea surface height above geoid (absolute dynamic topography)
call ufo_geovals_get_var(geovals, var_abs_topo, geoval_adt)

! Potential Temperature
call ufo_geovals_get_var(geovals, var_ocn_pot_temp, geoval_temp)

! Salinity
call ufo_geovals_get_var(geovals, var_ocn_salt, geoval_salt)

! Steric height obs operator
do iobs = 1, size(hofx,1)
   hofx(iobs) = geoval_adt%vals(1,iobs)
   write(102,*)hofx(iobs)
enddo

end subroutine ufo_stericheight_simobs


end module ufo_stericheight_mod
