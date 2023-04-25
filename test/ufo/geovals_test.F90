!
! (C) Copyright 2019 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
!

module test_ufo_geovals

use, intrinsic :: iso_c_binding
use kinds
use ufo_geovals_mod
use oops_variables_mod

implicit none
private

!-------------------------------------------------------------------------------
contains
!-------------------------------------------------------------------------------

! All these functions are expected to cause aborts.

subroutine test_ufo_geovals_setup_with_mismatched_nvars_c() &
  bind(c,name='test_ufo_geovals_setup_with_mismatched_nvars_f90')
implicit none

type(ufo_geovals) :: geovals
integer, parameter :: nlocs = 4
type(oops_variables) :: vars
integer, parameter :: nvars = 3  ! mismatch with vars%nvars()
integer(c_size_t) :: nvals(nvars)
integer, parameter :: nsampling_methods = 2
integer(c_size_t) :: npaths_by_method(nsampling_methods)
integer(c_size_t) :: sampling_method_by_var(nvars)

vars = oops_variables()
call vars%push_back("a")
call vars%push_back("b")

nvals(:) = 5
npaths_by_method(:) = 6
sampling_method_by_var(:) = 1

call ufo_geovals_setup(geovals, nlocs, vars, nvars, nvals, &
                       nsampling_methods, npaths_by_method, sampling_method_by_var)

end subroutine test_ufo_geovals_setup_with_mismatched_nvars_c

!-------------------------------------------------------------------------------

subroutine test_ufo_geovals_setup_with_sampling_method_set_to_0_c() &
  bind(c,name='test_ufo_geovals_setup_with_sampling_method_set_to_0_f90')
implicit none

type(ufo_geovals) :: geovals
integer, parameter :: nlocs = 4
type(oops_variables) :: vars
integer, parameter :: nvars = 3
integer(c_size_t) :: nvals(nvars)
integer, parameter :: nsampling_methods = 2
integer(c_size_t) :: npaths_by_method(nsampling_methods)
integer(c_size_t) :: sampling_method_by_var(nvars)

vars = oops_variables()
call vars%push_back("a")
call vars%push_back("b")
call vars%push_back("c")

nvals(:) = 5
npaths_by_method(:) = 6
sampling_method_by_var(:) = 1
sampling_method_by_var(2) = 0  ! out-of-bounds value

call ufo_geovals_setup(geovals, nlocs, vars, nvars, nvals, &
                       nsampling_methods, npaths_by_method, sampling_method_by_var)

end subroutine test_ufo_geovals_setup_with_sampling_method_set_to_0_c

!-------------------------------------------------------------------------------

subroutine test_ufo_geovals_setup_with_sampling_method_set_to_3_c() &
  bind(c,name='test_ufo_geovals_setup_with_sampling_method_set_to_3_f90')
implicit none

type(ufo_geovals) :: geovals
integer, parameter :: nlocs = 4
type(oops_variables) :: vars
integer, parameter :: nvars = 3
integer(c_size_t) :: nvals(nvars)
integer, parameter :: nsampling_methods = 2
integer(c_size_t) :: npaths_by_method(nsampling_methods)
integer(c_size_t) :: sampling_method_by_var(nvars)

vars = oops_variables()
call vars%push_back("a")
call vars%push_back("b")
call vars%push_back("c")

nvals(:) = 5
npaths_by_method(:) = 6
sampling_method_by_var(:) = 1
sampling_method_by_var(2) = 3  ! out-of-bounds value

call ufo_geovals_setup(geovals, nlocs, vars, nvars, nvals, &
                       nsampling_methods, npaths_by_method, sampling_method_by_var)

end subroutine test_ufo_geovals_setup_with_sampling_method_set_to_3_c

!-------------------------------------------------------------------------------

subroutine test_ufo_geovals_partial_setup_with_mismatched_nvars_c() &
  bind(c,name='test_ufo_geovals_partial_setup_with_mismatched_nvars_f90')
implicit none

type(ufo_geovals) :: geovals
integer, parameter :: nlocs = 4
type(oops_variables) :: vars
integer, parameter :: nvars = 3  ! mismatch with vars%nvars()
integer, parameter :: nsampling_methods = 2
integer(c_size_t) :: sampling_method_by_var(nvars)

vars = oops_variables()
call vars%push_back("a")
call vars%push_back("b")

sampling_method_by_var(:) = 1

call ufo_geovals_partial_setup(geovals, nlocs, vars, nvars,  &
                               nsampling_methods, sampling_method_by_var)

end subroutine test_ufo_geovals_partial_setup_with_mismatched_nvars_c

!-------------------------------------------------------------------------------

subroutine test_ufo_geovals_partial_setup_with_sampling_method_set_to_0_c() &
  bind(c,name='test_ufo_geovals_partial_setup_with_sampling_method_set_to_0_f90')
implicit none

type(ufo_geovals) :: geovals
integer, parameter :: nlocs = 4
type(oops_variables) :: vars
integer, parameter :: nvars = 3
integer, parameter :: nsampling_methods = 2
integer(c_size_t) :: sampling_method_by_var(nvars)

vars = oops_variables()
call vars%push_back("a")
call vars%push_back("b")
call vars%push_back("c")

sampling_method_by_var(:) = 1
sampling_method_by_var(2) = 0  ! out-of-bounds value

call ufo_geovals_partial_setup(geovals, nlocs, vars, nvars, &
                               nsampling_methods, sampling_method_by_var)

end subroutine test_ufo_geovals_partial_setup_with_sampling_method_set_to_0_c

!-------------------------------------------------------------------------------

subroutine test_ufo_geovals_partial_setup_with_sampling_method_set_to_3_c() &
  bind(c,name='test_ufo_geovals_partial_setup_with_sampling_method_set_to_3_f90')
implicit none

type(ufo_geovals) :: geovals
integer, parameter :: nlocs = 4
type(oops_variables) :: vars
integer, parameter :: nvars = 3
integer, parameter :: nsampling_methods = 2
integer(c_size_t) :: sampling_method_by_var(nvars)

vars = oops_variables()
call vars%push_back("a")
call vars%push_back("b")
call vars%push_back("c")

sampling_method_by_var(:) = 1
sampling_method_by_var(2) = 3  ! out-of-bounds value

call ufo_geovals_partial_setup(geovals, nlocs, vars, nvars, &
                               nsampling_methods, sampling_method_by_var)

end subroutine test_ufo_geovals_partial_setup_with_sampling_method_set_to_3_c

!-------------------------------------------------------------------------------

subroutine test_ufo_geovals_setup_sampling_method_0_c() &
  bind(c,name='test_ufo_geovals_setup_sampling_method_0_f90')
implicit none

type(ufo_geovals) :: geovals
integer, parameter :: nlocs = 4
type(oops_variables) :: vars
integer, parameter :: nvars = 3
integer, parameter :: nsampling_methods = 2
integer(c_size_t) :: sampling_method_by_var(nvars)
integer, parameter :: sampling_method = 0  ! out-of-bounds value
integer, parameter :: npaths = nlocs
type(ufo_index_range) :: paths_by_loc(nlocs)
integer :: i

vars = oops_variables()
call vars%push_back("a")
call vars%push_back("b")
call vars%push_back("c")

sampling_method_by_var(:) = 1
sampling_method_by_var(2) = 2

call ufo_geovals_partial_setup(geovals, nlocs, vars, nvars, &
                               nsampling_methods, sampling_method_by_var)

do i = 1, nlocs
  paths_by_loc(i)%begin = i
  paths_by_loc(i)%end = i + 1
enddo
call ufo_geovals_setup_sampling_method(geovals, sampling_method, npaths, nlocs, paths_by_loc)

end subroutine test_ufo_geovals_setup_sampling_method_0_c

!-------------------------------------------------------------------------------

subroutine test_ufo_geovals_setup_sampling_method_3_c() &
  bind(c,name='test_ufo_geovals_setup_sampling_method_3_f90')
implicit none

type(ufo_geovals) :: geovals
integer, parameter :: nlocs = 4
type(oops_variables) :: vars
integer, parameter :: nvars = 3
integer, parameter :: nsampling_methods = 2
integer(c_size_t) :: sampling_method_by_var(nvars)
integer, parameter :: sampling_method = 3  ! out-of-bounds value
integer, parameter :: npaths = nlocs
type(ufo_index_range) :: paths_by_loc(nlocs)
integer :: i

vars = oops_variables()
call vars%push_back("a")
call vars%push_back("b")
call vars%push_back("c")

sampling_method_by_var(:) = 1
sampling_method_by_var(2) = 2

call ufo_geovals_partial_setup(geovals, nlocs, vars, nvars, &
                               nsampling_methods, sampling_method_by_var)

do i = 1, nlocs
  paths_by_loc(i)%begin = i
  paths_by_loc(i)%end = i + 1
enddo
call ufo_geovals_setup_sampling_method(geovals, sampling_method, npaths, nlocs, paths_by_loc)

end subroutine test_ufo_geovals_setup_sampling_method_3_c

!-------------------------------------------------------------------------------

subroutine test_ufo_geovals_setup_sampling_method_with_mismatched_nlocs_c() &
  bind(c,name='test_ufo_geovals_setup_sampling_method_with_mismatched_nlocs_f90')
implicit none

type(ufo_geovals) :: geovals
integer, parameter :: nlocs = 4
type(oops_variables) :: vars
integer, parameter :: nvars = 3
integer, parameter :: nsampling_methods = 2
integer(c_size_t) :: sampling_method_by_var(nvars)
integer, parameter :: sampling_method = 1
integer, parameter :: npaths = nlocs
integer, parameter :: mismatched_nlocs = 2 * nlocs
type(ufo_index_range) :: paths_by_loc(mismatched_nlocs)
integer :: i

vars = oops_variables()
call vars%push_back("a")
call vars%push_back("b")
call vars%push_back("c")

sampling_method_by_var(:) = 1
sampling_method_by_var(2) = 2

call ufo_geovals_partial_setup(geovals, nlocs, vars, nvars, &
                               nsampling_methods, sampling_method_by_var)

do i = 1, nlocs
  paths_by_loc(i)%begin = i
  paths_by_loc(i)%end = i + 1
enddo
call ufo_geovals_setup_sampling_method(geovals, sampling_method, npaths, &
                                       mismatched_nlocs, paths_by_loc)

end subroutine test_ufo_geovals_setup_sampling_method_with_mismatched_nlocs_c

!-------------------------------------------------------------------------------

subroutine test_ufo_geovals_setup_trivial_sampling_method_0_c() &
  bind(c,name='test_ufo_geovals_setup_trivial_sampling_method_0_f90')
implicit none

type(ufo_geovals) :: geovals
integer, parameter :: nlocs = 4
type(oops_variables) :: vars
integer, parameter :: nvars = 3
integer, parameter :: nsampling_methods = 2
integer(c_size_t) :: sampling_method_by_var(nvars)
integer, parameter :: sampling_method = 0  ! out-of-bounds value

vars = oops_variables()
call vars%push_back("a")
call vars%push_back("b")
call vars%push_back("c")

sampling_method_by_var(:) = 1
sampling_method_by_var(2) = 2

call ufo_geovals_partial_setup(geovals, nlocs, vars, nvars, &
                               nsampling_methods, sampling_method_by_var)

call ufo_geovals_setup_trivial_sampling_method(geovals, sampling_method)

end subroutine test_ufo_geovals_setup_trivial_sampling_method_0_c

!-------------------------------------------------------------------------------

subroutine test_ufo_geovals_setup_trivial_sampling_method_3_c() &
  bind(c,name='test_ufo_geovals_setup_trivial_sampling_method_3_f90')
implicit none

type(ufo_geovals) :: geovals
integer, parameter :: nlocs = 4
type(oops_variables) :: vars
integer, parameter :: nvars = 3
integer, parameter :: nsampling_methods = 2
integer(c_size_t) :: sampling_method_by_var(nvars)
integer, parameter :: sampling_method = 3  ! out-of-bounds value

vars = oops_variables()
call vars%push_back("a")
call vars%push_back("b")
call vars%push_back("c")

sampling_method_by_var(:) = 1
sampling_method_by_var(2) = 2

call ufo_geovals_partial_setup(geovals, nlocs, vars, nvars, &
                               nsampling_methods, sampling_method_by_var)

call ufo_geovals_setup_trivial_sampling_method(geovals, sampling_method)

end subroutine test_ufo_geovals_setup_trivial_sampling_method_3_c

!-------------------------------------------------------------------------------

end module test_ufo_geovals
