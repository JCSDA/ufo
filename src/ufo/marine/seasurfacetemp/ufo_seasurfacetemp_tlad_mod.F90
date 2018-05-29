! (C) Copyright 2017 UCAR
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 

!> Fortran module to handle ice concentration observations

module ufo_seasurfacetemp_tlad_mod
  
use ioda_obs_seasurfacetemp_mod
use ioda_obs_vectors
use ufo_vars_mod
use ioda_locs_mod
use ufo_geovals_mod
use kinds
  
implicit none
public :: ufo_seasurfacetemp_tlad
public :: ufo_seasurfacetemp_tlad_settraj
public :: ufo_seasurfacetemp_tlad_eqv_tl
public :: ufo_seasurfacetemp_tlad_eqv_ad
private
integer, parameter :: max_string=800

!> Fortran derived type for sea ice fraction observation operator
type :: ufo_seasurfacetemp_tlad
end type ufo_seasurfacetemp_tlad


! ------------------------------------------------------------------------------

contains
 

! ------------------------------------------------------------------------------

subroutine ufo_seasurfacetemp_tlad_settraj(self, geovals)
implicit none
type(ufo_seasurfacetemp_tlad), intent(inout) :: self
type(ufo_geovals), intent(in)       :: geovals

character(len=*), parameter :: myname_="ufo_seasurfacetemp_tlad_settraj"
character(max_string) :: err_msg

type(ufo_geoval), pointer :: geoval

! since observation operator is linear, don't care about trajectory itself

end subroutine ufo_seasurfacetemp_tlad_settraj


! ------------------------------------------------------------------------------

subroutine ufo_seasurfacetemp_tlad_eqv_tl(self, geovals, hofx)
implicit none
type(ufo_seasurfacetemp_tlad), intent(in) :: self
type(ufo_geovals), intent(in)    :: geovals
type(obs_vector),  intent(inout) :: hofx

character(len=*), parameter :: myname_="ufo_seasurfacetemp_tlad_eqv_tl"
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
if (.not. ufo_geovals_get_var(geovals, var_ocn_sst, geoval)) then
  write(err_msg,*) myname_, trim(var_ocn_sst), ' doesnt exist'
  call abor1_ftn(err_msg)
endif

! sst obs operator
do iobs = 1, hofx%nobs
   hofx%values(iobs) = geoval%vals(1,iobs)
enddo

end subroutine ufo_seasurfacetemp_tlad_eqv_tl

! ------------------------------------------------------------------------------

subroutine ufo_seasurfacetemp_tlad_eqv_ad(self, geovals, hofx)
implicit none
type(ufo_seasurfacetemp_tlad), intent(in) :: self
type(ufo_geovals), intent(inout) :: geovals
type(obs_vector),  intent(in)    :: hofx

character(len=*), parameter :: myname_="ufo_seasurfacetemp_tlad_eqv_ad"
character(max_string) :: err_msg

integer :: iobs
type(ufo_geoval), pointer :: geoval

! check if nobs is consistent in geovals & hofx
if (geovals%nobs /= hofx%nobs) then
  write(err_msg,*) myname_, ' error: nobs inconsistent!'
  call abor1_ftn(err_msg)
endif

! check if sst variables is in geovals and get it
if (.not. ufo_geovals_get_var(geovals, var_ocn_sst, geoval)) then
  write(err_msg,*) myname_, trim(var_ocn_sst), ' doesnt exist'
  call abor1_ftn(err_msg)
endif

if (.not.(allocated(geoval%vals))) then
   geoval%nval=1
   allocate(geoval%vals(1,hofx%nobs))
end if
!print *,'in ad:',geoval%nval
!read(*,*)
! backward sst obs operator
geoval%vals=0.0
do iobs = 1, hofx%nobs
   geoval%vals(1,iobs) = geoval%vals(1,iobs) + hofx%values(iobs)
enddo

end subroutine ufo_seasurfacetemp_tlad_eqv_ad

end module ufo_seasurfacetemp_tlad_mod
