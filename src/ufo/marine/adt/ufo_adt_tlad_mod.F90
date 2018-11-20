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
use fckit_mpi_module, only: fckit_mpi_comm, fckit_mpi_sum
   
implicit none
public :: ufo_adt_tlad
public :: ufo_adt_tlad_delete
public :: ufo_adt_tlad_settraj
public :: ufo_adt_simobs_tl
public :: ufo_adt_simobs_ad
private
integer, parameter :: max_string=800

!> Fortran derived type for linear and adjoint adt observation operator
type :: ufo_adt_tlad
   integer          :: nlocs           !< Local number of obs
   real(c_double)   :: r_miss_val      !< Missing value flag  
   type(ufo_geoval) :: geoval_adt      !< adt (traj)
   logical          :: ltraj = .false. !< trajectory set?
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

subroutine ufo_adt_tlad_settraj(self, geovals, obss)
implicit none
type(ufo_adt_tlad), intent(inout) :: self
type(ufo_geovals),     intent(in) :: geovals
type(c_ptr),    value, intent(in) :: obss

character(len=*), parameter :: myname_="ufo_adt_tlad_settraj"
type(ufo_geoval), pointer :: geoval_adt
real(c_double) :: missing_value

self%nlocs = obsspace_get_nlocs(obss)

! check if adt variables is in geovals and get it
call ufo_geovals_get_var(geovals, var_abs_topo, geoval_adt)

self%geoval_adt = geoval_adt
self%ltraj    = .true.

! Set missing flag
self%r_miss_val = abs(obspace_missing_value())

end subroutine ufo_adt_tlad_settraj


! ------------------------------------------------------------------------------

subroutine ufo_adt_simobs_tl(self, geovals, hofx)
implicit none
type(ufo_adt_tlad), intent(in) :: self
type(ufo_geovals),  intent(in) :: geovals
real(c_double),  intent(inout) :: hofx(:)

character(len=*), parameter :: myname_="ufo_adt_simobs_tl"
character(max_string) :: err_msg
integer :: iobs, nobs, cnt, cnt_glb
type(ufo_geoval), pointer :: geoval_adt
real(kind_real) :: offset_hofx, pe_offset_hofx
type(fckit_mpi_comm) :: f_comm

f_comm = fckit_mpi_comm()

! check if trajectory was set
if (.not. self%ltraj) then
  write(err_msg,*) myname_, ' trajectory wasnt set!'
  call abor1_ftn(err_msg)
endif

! check if nobs is consistent in geovals & hofx
nobs = self%nlocs

if (geovals%nobs /= nobs) then
  write(err_msg,*) myname_, ' error: nobs inconsistent!'
  call abor1_ftn(err_msg)
endif

! check if adt variable is in geovals and get it
call ufo_geovals_get_var(geovals, var_abs_topo, geoval_adt)

! Local offset
pe_offset_hofx = 0.0
cnt = 0
do iobs = 1, self%nlocs
   if (abs(hofx(iobs)).lt.self%r_miss_val) then
      pe_offset_hofx = pe_offset_hofx + geoval_adt%vals(1,iobs)
      cnt = cnt + 1
   end if
end do

! Global offset
call f_comm%allreduce(pe_offset_hofx, offset_hofx, fckit_mpi_sum()) 
call f_comm%allreduce(cnt, cnt_glb, fckit_mpi_sum()) 
offset_hofx = offset_hofx/cnt_glb

! adt obs operator
hofx = 0.0
do iobs = 1, self%nlocs
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

integer :: iobs, nobs, cnt, cnt_glb
type(ufo_geoval), pointer :: geoval_adt
real(kind_real) :: offset_hofx, pe_offset_hofx
type(fckit_mpi_comm) :: f_comm

! check if trajectory was set
if (.not. self%ltraj) then
  write(err_msg,*) myname_, ' trajectory wasnt set!'
  call abor1_ftn(err_msg)
endif

f_comm = fckit_mpi_comm()

! check if nobs is consistent in geovals & hofx
nobs = self%nlocs
if (geovals%nobs /= nobs) then
  write(err_msg,*) myname_, ' error: nobs inconsistent!'
  call abor1_ftn(err_msg)
endif

if (.not. geovals%linit ) geovals%linit=.true.

! check if adt variable is in geovals and get it
call ufo_geovals_get_var(geovals, var_abs_topo, geoval_adt)

! Local offset
pe_offset_hofx = 0.0
cnt = 0
do iobs = 1, self%nlocs
   if (abs(hofx(iobs)).lt.self%r_miss_val) then
      pe_offset_hofx = pe_offset_hofx + hofx(iobs)
      cnt = cnt + 1
   end if
end do

! Global offset
call f_comm%allreduce(pe_offset_hofx, offset_hofx, fckit_mpi_sum()) 
call f_comm%allreduce(cnt, cnt_glb, fckit_mpi_sum()) 
offset_hofx = offset_hofx/cnt_glb

if (.not. allocated(geoval_adt%vals))  allocate(geoval_adt%vals(1,nobs))
geoval_adt%vals = 0.0
do iobs = 1, nobs
   if (abs(hofx(iobs)).lt.self%r_miss_val) then
      geoval_adt%vals(1,iobs) = geoval_adt%vals(1,iobs) + hofx(iobs) - offset_hofx
   end if
enddo

end subroutine ufo_adt_simobs_ad

end module ufo_adt_tlad_mod
