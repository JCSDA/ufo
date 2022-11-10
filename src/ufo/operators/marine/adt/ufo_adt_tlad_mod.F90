! (C) Copyright 2017-2018 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

!> Fortran adt module for tl/ad observation operator

module ufo_adt_tlad_mod

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
 type, extends(ufo_basis_tlad), public :: ufo_adt_tlad
 private
  integer          :: nlocs           !< Local number of obs
  real(c_double)   :: r_miss_val      !< Missing value flag
  type(ufo_geoval) :: geoval_adt      !< adt (traj)
 contains
  procedure :: setup  => ufo_adt_tlad_setup
  procedure :: delete  => ufo_adt_tlad_delete
  procedure :: settraj => ufo_adt_tlad_settraj
  procedure :: simobs_tl  => ufo_adt_simobs_tl
  procedure :: simobs_ad  => ufo_adt_simobs_ad
 end type ufo_adt_tlad

contains

! ------------------------------------------------------------------------------
subroutine ufo_adt_tlad_setup(self, f_conf)
implicit none
class(ufo_adt_tlad), intent(inout)    :: self
type(fckit_configuration), intent(in) :: f_conf

end subroutine ufo_adt_tlad_setup

! ------------------------------------------------------------------------------
subroutine ufo_adt_tlad_delete(self)
implicit none
class(ufo_adt_tlad), intent(inout) :: self

self%ltraj = .false.

end subroutine ufo_adt_tlad_delete

! ------------------------------------------------------------------------------
subroutine ufo_adt_tlad_settraj(self, geovals, obss)
implicit none
class(ufo_adt_tlad), intent(inout) :: self
type(ufo_geovals),       intent(in)    :: geovals
type(c_ptr), value,      intent(in)    :: obss

character(len=*), parameter :: myname_="ufo_adt_tlad_settraj"
type(ufo_geoval), pointer :: geoval_adt

self%nlocs = obsspace_get_nlocs(obss)

! check if adt variables is in geovals and get it
call ufo_geovals_get_var(geovals, var_abs_topo, geoval_adt)

self%geoval_adt = geoval_adt
self%ltraj    = .true.



end subroutine ufo_adt_tlad_settraj

! ------------------------------------------------------------------------------
subroutine ufo_adt_simobs_tl(self, geovals, hofx, obss)
use fckit_mpi_module, only: fckit_mpi_comm, fckit_mpi_sum
implicit none
class(ufo_adt_tlad), intent(in)    :: self
type(ufo_geovals),       intent(in)    :: geovals
real(c_double),          intent(inout) :: hofx(:)
type(c_ptr), value,      intent(in)    :: obss

character(len=*), parameter :: myname_="ufo_adt_simobs_tl"
character(max_string) :: err_msg
integer :: iobs, nlocs, cnt, cnt_glb
type(ufo_geoval), pointer :: geoval_adt
real(kind_real) :: offset, pe_offset
type(fckit_mpi_comm) :: f_comm
real(c_double) :: missing

call obsspace_get_comm(obss, f_comm)

! Set missing flag
missing = missing_value(missing)

! check if trajectory was set
if (.not. self%ltraj) then
  write(err_msg,*) myname_, ' trajectory wasnt set!'
  call abor1_ftn(err_msg)
endif

! check if nlocs is consistent in geovals & hofx
nlocs = self%nlocs

if (geovals%nlocs /= nlocs) then
  write(err_msg,*) myname_, ' error: nlocs inconsistent!'
  call abor1_ftn(err_msg)
endif

! check if adt variable is in geovals and get it
call ufo_geovals_get_var(geovals, var_abs_topo, geoval_adt)

! Local offset
pe_offset = 0.0
cnt = 0
do iobs = 1, self%nlocs
   if (geoval_adt%vals(1,iobs) /= missing) then
      pe_offset = pe_offset - geoval_adt%vals(1,iobs)
      cnt = cnt + 1
   end if
end do

! Global offset
call f_comm%allreduce(pe_offset, offset, fckit_mpi_sum())
call f_comm%allreduce(cnt, cnt_glb, fckit_mpi_sum())
if (cnt_glb > 0) then
   offset = offset/cnt_glb
end if

! adt obs operator
hofx = 0.0
do iobs = 1, self%nlocs
   hofx(iobs) = geoval_adt%vals(1,iobs)
   if (hofx(iobs) /= missing) then
      hofx(iobs) = hofx(iobs) + offset
   end if
enddo

end subroutine ufo_adt_simobs_tl

! ------------------------------------------------------------------------------
subroutine ufo_adt_simobs_ad(self, geovals, hofx, obss)
use fckit_mpi_module, only: fckit_mpi_comm, fckit_mpi_sum
implicit none
class(ufo_adt_tlad), intent(in)    :: self
type(ufo_geovals),       intent(inout) :: geovals
real(c_double),          intent(in)    :: hofx(:)
type(c_ptr), value,      intent(in)    :: obss

character(len=*), parameter :: myname_="ufo_adt_simobs_ad"
character(max_string) :: err_msg

integer :: iobs, nlocs, cnt, cnt_glb
type(ufo_geoval), pointer :: geoval_adt
real(kind_real) :: offset, pe_offset
type(fckit_mpi_comm) :: f_comm
real(c_double) :: missing

call obsspace_get_comm(obss, f_comm)

! Set missing flag
missing = missing_value(missing)

! check if trajectory was set
if (.not. self%ltraj) then
  write(err_msg,*) myname_, ' trajectory wasnt set!'
  call abor1_ftn(err_msg)
endif

! check if nlocs is consistent in geovals & hofx
nlocs = self%nlocs
if (geovals%nlocs /= nlocs) then
  write(err_msg,*) myname_, ' error: nlocs inconsistent!'
  call abor1_ftn(err_msg)
endif

if (.not. geovals%linit ) geovals%linit=.true.

! check if adt variable is in geovals and get it
call ufo_geovals_get_var(geovals, var_abs_topo, geoval_adt)

! Local offset
pe_offset = 0.0
cnt = 0
do iobs = 1, self%nlocs
   if (hofx(iobs) /= missing) then
      pe_offset = pe_offset + hofx(iobs)
      cnt = cnt + 1
   end if
end do

! Global offset
call f_comm%allreduce(pe_offset, offset, fckit_mpi_sum())
call f_comm%allreduce(cnt, cnt_glb, fckit_mpi_sum())
if (cnt_glb > 0) then
   offset = offset/cnt_glb
end if

do iobs = 1, nlocs
   if (geoval_adt%vals(1,iobs) /= missing .and. hofx(iobs) /= missing) then
      geoval_adt%vals(1,iobs) = geoval_adt%vals(1,iobs) + hofx(iobs) - offset
   end if
enddo

end subroutine ufo_adt_simobs_ad

end module ufo_adt_tlad_mod
