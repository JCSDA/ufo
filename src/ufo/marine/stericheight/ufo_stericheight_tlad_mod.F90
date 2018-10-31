! (C) Copyright 2017-2018 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

!> Fortran module to handle steric height operator

module ufo_stericheight_tlad_mod

use iso_c_binding
use ufo_vars_mod
use ufo_geovals_mod
use kinds

implicit none
public :: ufo_stericheight_tlad
public :: ufo_stericheight_tlad_settraj
public :: ufo_stericheight_simobs_tl
public :: ufo_stericheight_simobs_ad
private
integer, parameter :: max_string=800

!> Fortran derived type for steric height observation operator
type :: ufo_stericheight_tlad
   integer :: nl = -1      !< number of levels for T & S
end type ufo_stericheight_tlad


! ------------------------------------------------------------------------------

contains

! ------------------------------------------------------------------------------

subroutine ufo_stericheight_tlad_settraj(self, geovals)
implicit none
type(ufo_stericheight_tlad), intent(inout) :: self
type(ufo_geovals), intent(in)       :: geovals

character(len=*), parameter :: myname_="ufo_stericheight_tlad_settraj"
character(max_string) :: err_msg

type(ufo_geoval), pointer :: geoval

print *, myname_, ' nobs: ', geovals%nobs


call ufo_geovals_get_var(geovals, var_abs_topo, geoval)

!call ufo_geovals_get_var(geovals, var_stericheight, geoval)
print *,'==========================================='
self%nl = geoval%nval
print *, myname_, ' nval: ', geoval%nval

end subroutine ufo_stericheight_tlad_settraj


! ------------------------------------------------------------------------------

subroutine ufo_stericheight_simobs_tl(self, geovals, hofx)!, traj)
implicit none
type(ufo_stericheight_tlad), intent(in) :: self
type(ufo_geovals), intent(in)    :: geovals
real(c_double),  intent(inout) :: hofx(:)
!type(ufo_geovals), intent(in)    :: traj

character(len=*), parameter :: myname_="ufo_stericheight_simobs_tl"
character(max_string) :: err_msg

integer :: iobs
type(ufo_geoval), pointer :: geoval

print *, myname_, ' nobs: ', geovals%nobs, size(hofx,1)


! check if nobs is consistent in geovals & hofx
if (geovals%nobs /= size(hofx,1)) then
  write(err_msg,*) myname_, ' error: nobs inconsistent!'
  call abor1_ftn(err_msg)
endif

! check if sea ice fraction variables is in geovals and get it
!call ufo_geovals_get_var(geovals, var_stericheight, geoval)
call ufo_geovals_get_var(geovals, var_abs_topo, geoval)

! total sea ice fraction obs operator
do iobs = 1, size(hofx,1)
   hofx(iobs) = geoval%vals(1,iobs)
enddo

end subroutine ufo_stericheight_simobs_tl

! ------------------------------------------------------------------------------

subroutine ufo_stericheight_simobs_ad(self, geovals, hofx)
implicit none
type(ufo_stericheight_tlad), intent(in) :: self
type(ufo_geovals), intent(inout) :: geovals
real(c_double),  intent(inout) :: hofx(:)

character(len=*), parameter :: myname_="ufo_stericheight_simobs_ad"
character(max_string) :: err_msg

integer :: iobs
type(ufo_geoval), pointer :: geoval

print *,'&&&&&&&&&&&&7 in adjoint'
read(*,*)

! check if nobs is consistent in geovals & hofx
if (geovals%nobs /= size(hofx,1)) then
  write(err_msg,*) myname_, ' error: nobs inconsistent!'
  call abor1_ftn(err_msg)
endif

call ufo_geovals_get_var(geovals, var_abs_topo, geoval)

! check if sea ice fraction variables is in geovals and get it
!call ufo_geovals_get_var(geovals, var_stericheight, geoval)

if (.not.(allocated(geoval%vals))) then
   if (self%nl < 1) then
      !write(err_msg,*) myname_, ' unknown number of categories'
      !call abor1_ftn(err_msg)
   endif
   !allocate(geoval%vals(self%ncat,size(hofx,1)))
   allocate(geoval%vals(1,size(hofx,1)))
end if

if (.not. geovals%linit ) geovals%linit=.true.

! backward steric height obs operator
geoval%vals=0.0
do iobs = 1, size(hofx,1)
   geoval%vals(1,iobs) = geoval%vals(1,iobs) + size(hofx,1)
enddo

end subroutine ufo_stericheight_simobs_ad

end module ufo_stericheight_tlad_mod
