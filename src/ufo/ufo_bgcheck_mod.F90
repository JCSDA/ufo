! (C) Copyright 2017-2018 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

!> Fortran module to implement background check

module ufo_bgcheck_mod

use iso_c_binding
use kinds
use ufo_geovals_mod
use obsspace_mod
use config_mod

implicit none
public :: ufo_bgcheck, ufo_bgcheck_create, ufo_bgcheck_delete, ufo_bgcheck_prior, ufo_bgcheck_post
private

integer, parameter :: max_string_length=99 ! Yuk!

! ------------------------------------------------------------------------------

type :: ufo_bgcheck
  real(kind_real) :: threshold
  character(len=max_string_length) :: variable
  type(c_ptr) :: obsdb
end type ufo_bgcheck

! ------------------------------------------------------------------------------
contains
! ------------------------------------------------------------------------------

subroutine ufo_bgcheck_create(self, obspace, conf)
implicit none
type(ufo_bgcheck), intent(inout) :: self
type(c_ptr), value, intent(in)   :: obspace
type(c_ptr), intent(in)          :: conf

self%variable = config_get_string(conf, max_string_length, "variable")
self%threshold = config_get_real(conf, "threshold")
if (self%threshold<=0.0_kind_real) call abor1_ftn("ufo_bgcheck_create: Error threshold")
self%obsdb = obspace

end subroutine ufo_bgcheck_create

! ------------------------------------------------------------------------------

subroutine ufo_bgcheck_delete(self)
implicit none
type(ufo_bgcheck), intent(inout) :: self
end subroutine ufo_bgcheck_delete

! ------------------------------------------------------------------------------

subroutine ufo_bgcheck_prior(self, geovals)
implicit none
type(ufo_bgcheck), intent(in) :: self
type(ufo_geovals), intent(in) :: geovals
end subroutine ufo_bgcheck_prior

! ------------------------------------------------------------------------------

subroutine ufo_bgcheck_post(self, hofx)
implicit none
type(ufo_bgcheck), intent(in) :: self
real(c_double),  intent(in)   :: hofx(:)
integer :: iobs, jobs, ireject
real(kind_real) :: zmax
real(kind_real), allocatable :: yobs(:), yerr(:)
integer(c_int32_t), allocatable :: flags(:)
real(kind_real) :: missing
character(len=max_string_length) :: cerr

missing = obspace_missing_value()
iobs = size(hofx)
allocate(yobs(iobs))
allocate(yerr(iobs))
allocate(flags(iobs))
flags(:) = 0

call obsspace_get_db(self%obsdb, "ObsValue", trim(self%variable), yobs)
cerr = trim(self%variable)//'_err'
call obsspace_get_db(self%obsdb, "ObsError", trim(cerr), yerr)

zmax = 0.0
ireject = 0
do jobs = 1, iobs
  if (hofx(jobs)/=missing .and. yobs(jobs)/=missing .and. yerr(jobs)/=missing) then
    zmax = max(zmax, abs(hofx(jobs)-yobs(jobs))/yerr(jobs))
    if (abs(hofx(jobs)-yobs(jobs)) > yerr(jobs)*self%threshold) then
      flags(jobs) = 2
      ireject = ireject + 1
    endif
  else
    flags(jobs) = 1
  endif
enddo
write(0,*)'ufo_bgcheck_post reject = ',iobs,ireject,zmax

cerr = trim(self%variable)//'_qg'
call obsspace_put_db(self%obsdb, "QC", trim(cerr), flags)

deallocate(yobs)
deallocate(yerr)
deallocate(flags)

end subroutine ufo_bgcheck_post

! ------------------------------------------------------------------------------

end module ufo_bgcheck_mod
