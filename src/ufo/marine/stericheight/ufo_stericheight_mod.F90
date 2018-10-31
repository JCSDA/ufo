! (C) Copyright 2017-2018 UCAR
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 

!> Fortran module to handle steric height operator

module ufo_stericheight_mod

use iso_c_binding
use ufo_vars_mod
use ufo_geovals_mod
use kinds
  
implicit none
public :: ufo_stericheight
public :: ufo_stericheight_simobs
private
integer, parameter :: max_string=800

!> Fortran derived type for steric height observation operator
type :: ufo_stericheight
end type ufo_stericheight


! ------------------------------------------------------------------------------

contains
 
! ------------------------------------------------------------------------------

subroutine ufo_stericheight_simobs(self, geovals, hofx)
implicit none
type(ufo_stericheight), intent(in) :: self
type(ufo_geovals), intent(in)    :: geovals
real(c_double),  intent(inout) :: hofx(:)

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
