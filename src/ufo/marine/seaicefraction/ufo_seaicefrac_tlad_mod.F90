! (C) Copyright 2017 UCAR
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 

!> Fortran module to handle ice concentration observations

module ufo_seaicefrac_tlad_mod
  
use ioda_obs_seaicefrac_mod
use ioda_obs_vectors
use ufo_vars_mod
use ioda_locs_mod
use ufo_geovals_mod
use kinds
  
implicit none
public :: ufo_seaicefrac_tlad
public :: ufo_seaicefrac_tlad_settraj
public :: ufo_seaicefrac_tlad_eqv_tl
public :: ufo_seaicefrac_tlad_eqv_ad
private
integer, parameter :: max_string=800

!> Fortran derived type for sea ice fraction observation operator
type :: ufo_seaicefrac_tlad
   integer :: ncat = -1      !< number of ice categories
end type ufo_seaicefrac_tlad


! ------------------------------------------------------------------------------

contains
 

! ------------------------------------------------------------------------------

subroutine ufo_seaicefrac_tlad_settraj(self, geovals)
implicit none
type(ufo_seaicefrac_tlad), intent(inout) :: self
type(ufo_geovals), intent(in)       :: geovals

character(len=*), parameter :: myname_="ufo_seaicefrac_tlad_settraj"
character(max_string) :: err_msg

type(ufo_geoval), pointer :: geoval

! since observation operator is linear, only need to save the number
! of ice categories here, don't care about trajectory itself

! check if sea ice fraction variables is in geovals and get it
call ufo_geovals_get_var(geovals, var_seaicefrac, geoval)

self%ncat = geoval%nval

end subroutine ufo_seaicefrac_tlad_settraj


! ------------------------------------------------------------------------------

subroutine ufo_seaicefrac_tlad_eqv_tl(self, geovals, hofx)
implicit none
type(ufo_seaicefrac_tlad), intent(in) :: self
type(ufo_geovals), intent(in)    :: geovals
type(obs_vector),  intent(inout) :: hofx

character(len=*), parameter :: myname_="ufo_seaicefrac_tlad_eqv_tl"
character(max_string) :: err_msg

integer :: iobs
type(ufo_geoval), pointer :: geoval

print *, myname_, ' nobs: ', geovals%nobs, hofx%nobs

! check if nobs is consistent in geovals & hofx
if (geovals%nobs /= hofx%nobs) then
  write(err_msg,*) myname_, ' error: nobs inconsistent!'
  call abor1_ftn(err_msg)
endif

! check if sea ice fraction variables is in geovals and get it
call ufo_geovals_get_var(geovals, var_seaicefrac, geoval)

! total sea ice fraction obs operator
do iobs = 1, hofx%nobs
   hofx%values(iobs) = sum(geoval%vals(:,iobs))
enddo

end subroutine ufo_seaicefrac_tlad_eqv_tl

! ------------------------------------------------------------------------------

subroutine ufo_seaicefrac_tlad_eqv_ad(self, geovals, hofx)
implicit none
type(ufo_seaicefrac_tlad), intent(in) :: self
type(ufo_geovals), intent(inout) :: geovals
type(obs_vector),  intent(in)    :: hofx

character(len=*), parameter :: myname_="ufo_seaicefrac_tlad_eqv_ad"
character(max_string) :: err_msg

integer :: iobs
type(ufo_geoval), pointer :: geoval

! check if nobs is consistent in geovals & hofx
if (geovals%nobs /= hofx%nobs) then
  write(err_msg,*) myname_, ' error: nobs inconsistent!'
  call abor1_ftn(err_msg)
endif

if (.not. geovals%linit ) geovals%linit=.true.

! check if sea ice fraction variables is in geovals and get it
call ufo_geovals_get_var(geovals, var_seaicefrac, geoval)

if (.not.(allocated(geoval%vals))) then
   if (self%ncat < 1) then
     write(err_msg,*) myname_, ' unknown number of categories'
     call abor1_ftn(err_msg)
   endif
   allocate(geoval%vals(self%ncat,hofx%nobs))
end if

! backward sea ice fraction obs operator
geoval%vals=0.0
do iobs = 1, hofx%nobs
   geoval%vals(:,iobs) = geoval%vals(:,iobs) + hofx%values(iobs)
enddo

end subroutine ufo_seaicefrac_tlad_eqv_ad

end module ufo_seaicefrac_tlad_mod
