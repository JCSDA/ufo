! (C) Copyright 2017-2018 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

!> Fortran module for seaicefraction tl/ad observation operator

module ufo_seaicefraction_tlad_mod

 use fckit_configuration_module, only: fckit_configuration 
 use iso_c_binding
 use kinds

 use ufo_geovals_mod, only: ufo_geovals, ufo_geoval, ufo_geovals_get_var
 use ufo_basis_tlad_mod, only: ufo_basis_tlad
 use ufo_vars_mod
 use obsspace_mod
 use missing_values_mod

 implicit none
 private

 integer, parameter :: max_string=800

 !> Fortran derived type for the tl/ad observation operator
 type, extends(ufo_basis_tlad), public :: ufo_seaicefraction_tlad
 private
  integer :: ncat = -1      !< number of ice categories
 contains
  procedure :: setup  => ufo_seaicefraction_tlad_setup
  procedure :: delete  => ufo_seaicefraction_tlad_delete
  procedure :: settraj => ufo_seaicefraction_tlad_settraj
  procedure :: simobs_tl  => ufo_seaicefraction_simobs_tl
  procedure :: simobs_ad  => ufo_seaicefraction_simobs_ad
 end type ufo_seaicefraction_tlad

contains

! ------------------------------------------------------------------------------
subroutine ufo_seaicefraction_tlad_setup(self, f_conf)
implicit none
class(ufo_seaicefraction_tlad), intent(inout) :: self
type(fckit_configuration), intent(in)         :: f_conf

end subroutine ufo_seaicefraction_tlad_setup

! ------------------------------------------------------------------------------
subroutine ufo_seaicefraction_tlad_delete(self)
implicit none
class(ufo_seaicefraction_tlad), intent(inout) :: self

end subroutine ufo_seaicefraction_tlad_delete

! ------------------------------------------------------------------------------
subroutine ufo_seaicefraction_tlad_settraj(self, geovals, obss)
implicit none
class(ufo_seaicefraction_tlad), intent(inout) :: self
type(ufo_geovals),       intent(in)    :: geovals
type(c_ptr), value,      intent(in)    :: obss

type(ufo_geoval), pointer :: geoval

! since observation operator is linear, only need to save the number
! of ice categories here, don't care about trajectory itself

! get sea ice fraction variables from geovals
call ufo_geovals_get_var(geovals, var_seaicefrac, geoval)

self%ncat = geoval%nval

end subroutine ufo_seaicefraction_tlad_settraj

! ------------------------------------------------------------------------------
subroutine ufo_seaicefraction_simobs_tl(self, geovals, hofx, obss)
implicit none
class(ufo_seaicefraction_tlad), intent(in)    :: self
type(ufo_geovals),       intent(in)    :: geovals
real(c_double),          intent(inout) :: hofx(:)
type(c_ptr), value,      intent(in)    :: obss

character(len=*), parameter :: myname_="ufo_seaicefrac_simobs_tl"
character(max_string) :: err_msg

integer :: iobs
type(ufo_geoval), pointer :: geoval

print *, myname_, ' nlocs: ', geovals%nlocs, size(hofx,1)

! check if nlocs is consistent in geovals & hofx
if (geovals%nlocs /= size(hofx,1)) then
  write(err_msg,*) myname_, ' error: nlocs inconsistent!'
  call abor1_ftn(err_msg)
endif

! check if sea ice fraction variables is in geovals and get it
call ufo_geovals_get_var(geovals, var_seaicefrac, geoval)

! total sea ice fraction obs operator
do iobs = 1, size(hofx,1)
   hofx(iobs) = sum(geoval%vals(:,iobs))
enddo

end subroutine ufo_seaicefraction_simobs_tl

! ------------------------------------------------------------------------------
subroutine ufo_seaicefraction_simobs_ad(self, geovals, hofx, obss)
implicit none
class(ufo_seaicefraction_tlad), intent(in)    :: self
type(ufo_geovals),       intent(inout) :: geovals
real(c_double),          intent(in)    :: hofx(:)
type(c_ptr), value,      intent(in)    :: obss

character(len=*), parameter :: myname_="ufo_seaicefrac_simobs_ad"
character(max_string) :: err_msg

integer :: iobs
type(ufo_geoval), pointer :: geoval
real(c_double) :: missing

!> Set missing value
missing = missing_value(missing)

! check if nlocs is consistent in geovals & hofx
if (geovals%nlocs /= size(hofx,1)) then
  write(err_msg,*) myname_, ' error: nlocs inconsistent!'
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
   allocate(geoval%vals(self%ncat,size(hofx,1)))
end if

! backward sea ice fraction obs operator
geoval%vals=0.0
do iobs = 1, size(hofx,1)
 if (hofx(iobs) /= missing) then
   geoval%vals(:,iobs) = geoval%vals(:,iobs) + hofx(iobs)
 end if
enddo

end subroutine ufo_seaicefraction_simobs_ad

end module ufo_seaicefraction_tlad_mod
