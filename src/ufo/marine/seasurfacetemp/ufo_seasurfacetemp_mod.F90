! (C) Copyright 2017 UCAR
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 

!> Fortran module to handle ice concentration observations

module ufo_seasurfacetemp_mod
  
use ioda_obs_seasurfacetemp_mod
use ioda_obs_vectors
use ufo_vars_mod
use ioda_locs_mod
use ufo_geovals_mod
use kinds
  
implicit none
public :: ufo_seasurfacetemp
public :: ufo_seasurfacetemp_eqv
private
integer, parameter :: max_string=800

!> Fortran derived type for sea surface temperature observation operator
type :: ufo_seasurfacetemp
end type ufo_seasurfacetemp


! ------------------------------------------------------------------------------

contains
 
! ------------------------------------------------------------------------------

subroutine ufo_seasurfacetemp_eqv(self, geovals, hofx)
implicit none
type(ufo_seasurfacetemp), intent(in) :: self
type(ufo_geovals), intent(in)    :: geovals
type(obs_vector),  intent(inout) :: hofx

character(len=*), parameter :: myname_="ufo_seasurfacetemp_eqv"
character(max_string) :: err_msg

integer :: iobs
type(ufo_geoval), pointer :: geoval

print *, myname_, ' nobs: ', geovals%nobs, hofx%nobs

! check if nobs is consistent in geovals & hofx
if (geovals%nobs /= hofx%nobs) then
  write(err_msg,*) myname_, ' error: nobs inconsistent!'
  call abor1_ftn(err_msg)
endif

! check if sst variables is in geovals and get it
if (.not. ufo_geovals_get_var(geovals, var_ocn_pot_temp, geoval)) then
  write(err_msg,*) myname_, trim(var_ocn_pot_temp), ' doesnt exist'
  call abor1_ftn(err_msg)
endif

! sst obs operator
do iobs = 1, hofx%nobs
   hofx%values(iobs) = geoval%vals(1,iobs)
   write(102,*)hofx%values(iobs)
enddo

end subroutine ufo_seasurfacetemp_eqv

end module ufo_seasurfacetemp_mod
