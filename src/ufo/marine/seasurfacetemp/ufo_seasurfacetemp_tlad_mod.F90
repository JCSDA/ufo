! (C) Copyright 2017-2018 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

!> Fortran module for seasurfacetemp tl/ad observation operator

module ufo_seasurfacetemp_tlad_mod

 use iso_c_binding
 use config_mod
 use kinds

 use ufo_geovals_mod, only: ufo_geovals, ufo_geoval, ufo_geovals_get_var
 use ufo_basis_tlad_mod, only: ufo_basis_tlad
 use ufo_vars_mod
 use obsspace_mod

 implicit none
 private

 integer, parameter :: max_string=800

 !> Fortran derived type for the tl/ad observation operator
 type, extends(ufo_basis_tlad), public :: ufo_seasurfacetemp_tlad
 private
 contains
  procedure :: setup  => ufo_seasurfacetemp_tlad_setup
  procedure :: delete  => ufo_seasurfacetemp_tlad_delete
  procedure :: settraj => ufo_seasurfacetemp_tlad_settraj
  procedure :: simobs_tl  => ufo_seasurfacetemp_simobs_tl
  procedure :: simobs_ad  => ufo_seasurfacetemp_simobs_ad
 end type ufo_seasurfacetemp_tlad

contains

! ------------------------------------------------------------------------------
subroutine ufo_seasurfacetemp_tlad_setup(self, c_conf)
implicit none
class(ufo_seasurfacetemp_tlad), intent(inout) :: self
type(c_ptr),              intent(in)    :: c_conf

end subroutine ufo_seasurfacetemp_tlad_setup

! ------------------------------------------------------------------------------
subroutine ufo_seasurfacetemp_tlad_delete(self)
implicit none
class(ufo_seasurfacetemp_tlad), intent(inout) :: self

end subroutine ufo_seasurfacetemp_tlad_delete

! ------------------------------------------------------------------------------
subroutine ufo_seasurfacetemp_tlad_settraj(self, geovals, obss)
implicit none
class(ufo_seasurfacetemp_tlad), intent(inout) :: self
type(ufo_geovals),       intent(in)    :: geovals
type(c_ptr), value,      intent(in)    :: obss

! since observation operator is linear, don't care about trajectory itself

end subroutine ufo_seasurfacetemp_tlad_settraj

! ------------------------------------------------------------------------------
subroutine ufo_seasurfacetemp_simobs_tl(self, geovals, hofx, obss)
implicit none
class(ufo_seasurfacetemp_tlad), intent(in)    :: self
type(ufo_geovals),       intent(in)    :: geovals
real(c_double),          intent(inout) :: hofx(:)
type(c_ptr), value,      intent(in)    :: obss

character(len=*), parameter :: myname_="ufo_seasurfacetemp_simobs_tl"
character(max_string) :: err_msg

integer :: iobs
type(ufo_geoval), pointer :: geoval_sst

print *, myname_, ' nobs: ', geovals%nobs, size(hofx,1)

! check if nobs is consistent in geovals & hofx
if (geovals%nobs /= size(hofx,1)) then
  write(err_msg,*) myname_, ' error: nobs inconsistent!'
  call abor1_ftn(err_msg)
endif

! check if sst variables is in geovals and get it
call ufo_geovals_get_var(geovals, var_ocn_sst, geoval_sst)

! sst obs operator
do iobs = 1, size(hofx,1)
   hofx(iobs) = geoval_sst%vals(1,iobs)
enddo

end subroutine ufo_seasurfacetemp_simobs_tl

! ------------------------------------------------------------------------------
subroutine ufo_seasurfacetemp_simobs_ad(self, geovals, hofx, obss)
implicit none
class(ufo_seasurfacetemp_tlad), intent(in)    :: self
type(ufo_geovals),       intent(inout) :: geovals
real(c_double),          intent(in)    :: hofx(:)
type(c_ptr), value,      intent(in)    :: obss

character(len=*), parameter :: myname_="ufo_seasurfacetemp_simobs_ad"
character(max_string) :: err_msg

integer :: iobs
type(ufo_geoval), pointer :: geoval_sst

! check if nobs is consistent in geovals & hofx
if (geovals%nobs /= size(hofx,1)) then
  write(err_msg,*) myname_, ' error: nobs inconsistent!'
  call abor1_ftn(err_msg)
endif

if (.not. geovals%linit ) geovals%linit=.true.

! check if sst variables is in geovals and get it
call ufo_geovals_get_var(geovals, var_ocn_sst, geoval_sst)

if (.not.(allocated(geoval_sst%vals))) then
   geoval_sst%nval=1
   allocate(geoval_sst%vals(1,size(hofx,1)))
end if

! backward sst obs operator
geoval_sst%vals=0.0
do iobs = 1, size(hofx,1)
   geoval_sst%vals(1,iobs) = geoval_sst%vals(1,iobs) + hofx(iobs)
enddo

end subroutine ufo_seasurfacetemp_simobs_ad

! ------------------------------------------------------------------------------

end module ufo_seasurfacetemp_tlad_mod
