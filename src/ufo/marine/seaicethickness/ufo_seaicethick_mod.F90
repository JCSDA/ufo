! (C) Copyright 2017-2018 UCAR
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 

!> Fortran module to handle ice concentration observations

module ufo_seaicethick_mod
  
use ioda_obs_seaicethick_mod
use ioda_obs_vectors
use ufo_vars_mod
use ioda_locs_mod
use ufo_geovals_mod
use kinds
  
implicit none
public :: ufo_seaicethick
public :: ufo_seaicethick_eqv
private
integer, parameter :: max_string=800

!> Fortran derived type for sea ice fraction observation operator
type :: ufo_seaicethick
end type ufo_seaicethick


! ------------------------------------------------------------------------------

contains
 
! ------------------------------------------------------------------------------

subroutine ufo_seaicethick_eqv(self, geovals, hofx)
implicit none
type(ufo_seaicethick), intent(in) :: self
type(ufo_geovals), intent(in)    :: geovals
type(obs_vector),  intent(inout) :: hofx

character(len=*), parameter :: myname_="ufo_seaicethick_eqv"
character(max_string) :: err_msg

integer :: iobs, icat, ncat
type(ufo_geoval), pointer :: icethick, icefrac

print *, myname_, ' nobs: ', geovals%nobs, hofx%nobs

! check if nobs is consistent in geovals & hofx
if (geovals%nobs /= hofx%nobs) then
  write(err_msg,*) myname_, ' error: nobs inconsistent!'
  call abor1_ftn(err_msg)
endif

! check if sea ice fraction variable is in geovals and get it
call ufo_geovals_get_var(geovals, var_seaicefrac, icefrac)
! check if sea ice thickness variable is in geovals and get it
call ufo_geovals_get_var(geovals, var_seaicethick, icethick)

ncat = icefrac%nval
hofx%values = 0.0
! total sea ice fraction obs operator
do iobs = 1, hofx%nobs
   do icat = 1, ncat
     hofx%values(iobs) = hofx%values(iobs) + icefrac%vals(icat,iobs) * icethick%vals(icat,iobs) / 905.0
   enddo
   write(102,*)hofx%values(iobs)
enddo

end subroutine ufo_seaicethick_eqv

end module ufo_seaicethick_mod
