! (C) Copyright 2017 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

!> Fortran module to handle adt observations

module ufo_adt_tlad_mod

use ufo_vars_mod
use ioda_locs_mod
use ufo_geovals_mod
use kinds
use iso_c_binding
use obsspace_mod

implicit none
public :: ufo_adt_tlad
public :: ufo_adt_tlad_delete
public :: ufo_adt_tlad_settraj
public :: ufo_adt_simobs_tl
public :: ufo_adt_simobs_ad
private
integer, parameter :: max_string=800

!> Fortran derived type for adt observation operator
type :: ufo_adt_tlad
   type(ufo_geoval) :: geoval_adt !< adt (traj)
   logical :: ltraj = .false.   !< trajectory set?
end type ufo_adt_tlad


! ------------------------------------------------------------------------------

contains

! ------------------------------------------------------------------------------

subroutine ufo_adt_tlad_delete(self)
implicit none
type(ufo_adt_tlad), intent(inout) :: self

self%ltraj = .false.

end subroutine ufo_adt_tlad_delete

! ------------------------------------------------------------------------------

subroutine ufo_adt_tlad_settraj(self, geovals)
implicit none
type(ufo_adt_tlad), intent(inout) :: self
type(ufo_geovals), intent(in)       :: geovals

character(len=*), parameter :: myname_="ufo_adt_tlad_settraj"
character(max_string) :: err_msg

type(ufo_geoval), pointer :: geoval_adt

! check if adt variables is in geovals and get it
call ufo_geovals_get_var(geovals, var_abs_topo, geoval_adt)

self%geoval_adt = geoval_adt
self%ltraj    = .true.

end subroutine ufo_adt_tlad_settraj


! ------------------------------------------------------------------------------

subroutine ufo_adt_simobs_tl(self, geovals, hofx)
implicit none
type(ufo_adt_tlad), intent(in) :: self
type(ufo_geovals),  intent(in) :: geovals
real(c_double),  intent(inout) :: hofx(:)

character(len=*), parameter :: myname_="ufo_adt_simobs_tl"
character(max_string) :: err_msg

integer :: iobs, nobs
type(ufo_geoval), pointer :: geoval_adt
real(kind_real) :: offset_obs, offset_hofx

! check if trajectory was set
if (.not. self%ltraj) then
  write(err_msg,*) myname_, ' trajectory wasnt set!'
  call abor1_ftn(err_msg)
endif

! check if nobs is consistent in geovals & hofx
nobs = size(hofx,1)
if (geovals%nobs /= nobs) then
  write(err_msg,*) myname_, ' error: nobs inconsistent!'
  call abor1_ftn(err_msg)
endif

! check if adt variable is in geovals and get it
call ufo_geovals_get_var(geovals, var_abs_topo, geoval_adt)

! Compute offset
offset_hofx=sum(geoval_adt%vals(1,:))/nobs

! adt obs operator
hofx = 0.0
do iobs = 1, nobs
     hofx(iobs) = geoval_adt%vals(1,iobs) - offset_hofx
enddo

end subroutine ufo_adt_simobs_tl

! ------------------------------------------------------------------------------

subroutine ufo_adt_simobs_ad(self, geovals, hofx)
implicit none
type(ufo_adt_tlad),    intent(in) :: self
type(ufo_geovals),  intent(inout) :: geovals
real(c_double),     intent(inout) :: hofx(:)

character(len=*), parameter :: myname_="ufo_adt_simobs_ad"
character(max_string) :: err_msg

integer :: iobs, nobs
type(ufo_geoval), pointer :: geoval_adt
real(kind_real) :: offset_hofx

! check if trajectory was set
if (.not. self%ltraj) then
  write(err_msg,*) myname_, ' trajectory wasnt set!'
  call abor1_ftn(err_msg)
endif

! check if nobs is consistent in geovals & hofx
nobs = size(hofx,1)
if (geovals%nobs /= nobs) then
  write(err_msg,*) myname_, ' error: nobs inconsistent!'
  call abor1_ftn(err_msg)
endif

if (.not. geovals%linit ) geovals%linit=.true.

! check if adt variable is in geovals and get it
call ufo_geovals_get_var(geovals, var_abs_topo, geoval_adt)

! backward adt obs operator

! Compute offset
offset_hofx=sum(hofx)/nobs

if (.not. allocated(geoval_adt%vals))  allocate(geoval_adt%vals(1,nobs))
geoval_adt%vals = 0.0
do iobs = 1, nobs
      geoval_adt%vals(1,iobs) = geoval_adt%vals(1,iobs) + hofx(iobs) - offset_hofx 
enddo

end subroutine ufo_adt_simobs_ad

end module ufo_adt_tlad_mod
