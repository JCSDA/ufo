! (C) Copyright 2017 UCAR
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 

!> Fortran module to handle steric height operator

module ufo_stericheight_mod
  
use ufo_obs_stericheight_mod
use ufo_obs_vectors
use ufo_vars_mod
use ufo_locs_mod
use ufo_geovals_mod
use kinds
  
implicit none
public :: ufo_stericheight
public :: ufo_stericheight_eqv
private
integer, parameter :: max_string=800

!> Fortran derived type for steric height observation operator
type :: ufo_stericheight
end type ufo_stericheight


! ------------------------------------------------------------------------------

contains
 
! ------------------------------------------------------------------------------

subroutine ufo_stericheight_eqv(self, geovals, hofx)
implicit none
type(ufo_stericheight), intent(in) :: self
type(ufo_geovals), intent(in)    :: geovals
type(obs_vector),  intent(inout) :: hofx

character(len=*), parameter :: myname_="ufo_stericheight_eqv"
character(max_string) :: err_msg

integer :: iobs
type(ufo_geoval), pointer :: geoval_temp, geoval_salt, geoval_adt

print *,myname_

hofx%nobs = geovals%nobs
print *, myname_, ' nobs: ', geovals%nobs, hofx%nobs

! check if nobs is consistent in geovals & hofx
if (geovals%nobs /= hofx%nobs) then
  write(err_msg,*) myname_, ' error: nobs inconsistent!'
  call abor1_ftn(err_msg)
endif

! Sea surface height above geoid (absolute dynamic topography)
if (.not. ufo_geovals_get_var(geovals, var_abs_topo, geoval_adt)) then
  write(err_msg,*) myname_, trim(var_abs_topo), ' doesnt exist'
  call abor1_ftn(err_msg)
endif

! Potential Temperature
if (.not. ufo_geovals_get_var(geovals, var_ocn_pot_temp, geoval_temp)) then
  write(err_msg,*) myname_, trim(var_ocn_pot_temp), ' doesnt exist'
  call abor1_ftn(err_msg)
endif

! Salinity
if (.not. ufo_geovals_get_var(geovals, var_ocn_salt, geoval_salt)) then
  write(err_msg,*) myname_, trim(var_ocn_salt), ' doesnt exist'
  call abor1_ftn(err_msg)
endif

! Steric height obs operator
do iobs = 1, hofx%nobs
   hofx%values(iobs) = geoval_adt%vals(1,iobs)
   write(102,*)hofx%values(iobs)
enddo

end subroutine ufo_stericheight_eqv


end module ufo_stericheight_mod
