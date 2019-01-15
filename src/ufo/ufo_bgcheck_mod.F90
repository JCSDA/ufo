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
public :: ufo_bgcheck, ufo_bgcheck_create, ufo_bgcheck_delete, ufo_bgcheck_post, &
        & max_string_length
private

integer, parameter :: max_string_length=99 ! Yuk!

! ------------------------------------------------------------------------------

type :: ufo_bgcheck
  integer :: nvars
  character(len=max_string_length), allocatable :: variables(:)
  real(kind_real) :: threshold
  type(c_ptr) :: obsdb
  character(len=max_string_length) :: qcname
end type ufo_bgcheck

! ------------------------------------------------------------------------------
contains
! ------------------------------------------------------------------------------

subroutine ufo_bgcheck_create(self, obspace, conf)
implicit none
type(ufo_bgcheck), intent(inout) :: self
type(c_ptr), value, intent(in)   :: obspace
type(c_ptr), intent(in)          :: conf

self%nvars = size(config_get_string_vector(conf, max_string_length, "variables"))
allocate(self%variables(self%nvars))
self%variables = config_get_string_vector(conf, max_string_length, "variables")
self%threshold = config_get_real(conf, "threshold")
if (self%threshold<=0.0_kind_real) call abor1_ftn("ufo_bgcheck_create: Error threshold")
self%obsdb = obspace
self%qcname = config_get_string(conf, max_string_length, "QCname")

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

subroutine ufo_bgcheck_post(self, hofx, hofxvars)
use fckit_log_module, only : fckit_log
use missing_values_mod
implicit none
type(ufo_bgcheck), intent(in) :: self
real(c_double),  intent(in)   :: hofx(:,:)
character(len=max_string_length) :: hofxvars(:)

integer :: iloc, jvar, ireject, icount, jobs, ivar
real(kind_real), allocatable :: yobs(:), yerr(:)
integer(c_int32_t), allocatable :: flags(:)
real(kind_real) :: missing
character(len=max_string_length) :: var
character(len=250)  :: buf

missing = missing_value(missing)
iloc = obsspace_get_nlocs(self%obsdb)
allocate(yobs(iloc))
allocate(yerr(iloc))
allocate(flags(iloc))

do jvar = 1, self%nvars
  var = trim(self%variables(jvar))
  ivar = find_index(hofxvars, var)

  call obsspace_get_db(self%obsdb, "ObsValue", var, yobs)
  call obsspace_get_db(self%obsdb, "ObsError", var, yerr)
  call obsspace_get_db(self%obsdb, self%qcname, var,flags )

  ireject = 0
  icount = 0
  do jobs = 1, iloc
    if (flags(jobs) == 0) then
      icount = icount + 1
      if (hofx(ivar, jobs)/=missing .and. yobs(jobs)/=missing .and. yerr(jobs)/=missing) then
        if (abs(hofx(ivar, jobs)-yobs(jobs)) > yerr(jobs)*self%threshold) then
          flags(jobs) = 10
          ireject = ireject + 1
        endif
      else
        flags(jobs) = 2
        ireject = ireject + 1
      endif
    endif
  enddo
  write(buf,*)'UFO Background Check: ',ireject,trim(var),' rejected out of ',icount,' (',iloc,' total)'
  call fckit_log%info(buf)

  call obsspace_put_db(self%obsdb, self%qcname, var, flags)
enddo

deallocate(yobs)
deallocate(yerr)
deallocate(flags)

end subroutine ufo_bgcheck_post

! ------------------------------------------------------------------------------

! This function should be moved to a variables_mod to handle variables in fortran
function find_index(vars, var)
character(len=max_string_length) :: vars(:)
character(len=max_string_length) :: var
integer :: find_index, jv
character(len=250)  :: buf

find_index = -1
do jv = 1, size(vars)
  if (trim(vars(jv)) == trim(var)) find_index = jv
enddo
write(buf,*)'UFO Background Check: unknown variable: ',trim(var)
if (find_index < 1) call abor1_ftn(buf)

end function find_index

! ------------------------------------------------------------------------------

end module ufo_bgcheck_mod
