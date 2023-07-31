!
! (C) Copyright 2017-2018 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
!
module ufo_geovals_mod_c

use fckit_configuration_module, only: fckit_configuration
use iso_c_binding
use ufo_geovals_mod
use kinds

implicit none

public :: ufo_geovals_registry

private
integer, parameter :: max_string=800

#define LISTED_TYPE ufo_geovals

!> Linked list interface - defines registry_t type
#include "oops/util/linkedList_i.f"

!> Global registry
type(registry_t) :: ufo_geovals_registry

! ------------------------------------------------------------------------------
contains
! ------------------------------------------------------------------------------
!> Linked list implementation
#include "oops/util/linkedList_c.f"
! ------------------------------------------------------------------------------
!> Setup GeoVaLs (don't store anything; don't do allocation yet)
subroutine ufo_geovals_default_constr_c(c_key_self) bind(c,name='ufo_geovals_default_constr_f90')
implicit none
integer(c_int), intent(inout)  :: c_key_self
type(ufo_geovals), pointer :: self

call ufo_geovals_registry%init()
call ufo_geovals_registry%add(c_key_self)
call ufo_geovals_registry%get(c_key_self, self)
call ufo_geovals_default_constr(self)

end subroutine ufo_geovals_default_constr_c

! ------------------------------------------------------------------------------
!> Setup and allocate GeoVaLs
subroutine ufo_geovals_setup_c(c_key_self, c_nlocs, c_vars, c_nvars, c_sizes, c_nsampling_methods, &
                               c_npaths_by_method, c_sampling_method_by_var, &
                               c_reduced_vars, c_nreduced_vars, c_reduced_sizes, &
                               c_is_sampling_method_trivial) &
                               bind(c,name='ufo_geovals_setup_f90')
use oops_variables_mod
implicit none
integer(c_int), intent(inout)  :: c_key_self
integer(c_int), intent(in)     :: c_nlocs, c_nvars, c_nreduced_vars
type(c_ptr), value, intent(in) :: c_vars, c_reduced_vars
integer(c_size_t), intent(in)  :: c_sizes(c_nvars), c_reduced_sizes(c_nreduced_vars)
integer(c_size_t), intent(in)  :: c_nsampling_methods
integer(c_size_t), intent(in)  :: c_npaths_by_method(c_nsampling_methods)
integer(c_size_t), intent(in)  :: c_sampling_method_by_var(c_nvars)
logical(c_bool), intent(in)    :: c_is_sampling_method_trivial(c_nsampling_methods)

type(ufo_geovals), pointer :: self
type(oops_variables) :: vars, reduced_vars
integer              :: nsampling_methods
integer(c_size_t)    :: sampling_method_by_var(c_nvars)

call ufo_geovals_registry%init()
call ufo_geovals_registry%add(c_key_self)
call ufo_geovals_registry%get(c_key_self, self)

vars = oops_variables(c_vars)
nsampling_methods = c_nsampling_methods
sampling_method_by_var(:) = 1 + c_sampling_method_by_var(:)  ! convert 0- to 1-based indexing
reduced_vars = oops_variables(c_reduced_vars)
call ufo_geovals_setup(self, c_nlocs, vars, c_nvars, c_sizes, &
                       nsampling_methods, c_npaths_by_method, sampling_method_by_var, &
                       reduced_vars, c_nreduced_vars, c_reduced_sizes, c_is_sampling_method_trivial)

end subroutine ufo_geovals_setup_c

! ------------------------------------------------------------------------------
!> Setup GeoVaLs (store nlocs, variables; don't do allocation yet)
subroutine ufo_geovals_partial_setup_c(c_key_self, c_nlocs, c_vars, c_nvars, c_nsampling_methods, &
                                       c_sampling_method_by_var, &
                                       c_reduced_vars, c_is_sampling_method_trivial) &
                                       bind(c,name='ufo_geovals_partial_setup_f90')
use oops_variables_mod
implicit none
integer(c_int), intent(inout)  :: c_key_self
integer(c_int), intent(in)     :: c_nlocs, c_nvars
type(c_ptr), value, intent(in) :: c_vars
integer(c_size_t), intent(in)  :: c_nsampling_methods
integer(c_size_t), intent(in)  :: c_sampling_method_by_var(c_nvars)
type(c_ptr), value, intent(in) :: c_reduced_vars
logical(c_bool), intent(in)    :: c_is_sampling_method_trivial(c_nsampling_methods)

type(ufo_geovals), pointer :: self
type(oops_variables) :: vars, reduced_vars
integer              :: nsampling_methods
integer(c_size_t)    :: sampling_method_by_var(c_nvars)

call ufo_geovals_registry%init()
call ufo_geovals_registry%add(c_key_self)
call ufo_geovals_registry%get(c_key_self, self)

vars = oops_variables(c_vars)
nsampling_methods = c_nsampling_methods
sampling_method_by_var(:) = 1 + c_sampling_method_by_var(:)  ! convert 0- to 1-based indexing
reduced_vars = oops_variables(c_reduced_vars)
call ufo_geovals_partial_setup(self, c_nlocs, vars, c_nvars, nsampling_methods, &
                               sampling_method_by_var, reduced_vars, c_is_sampling_method_trivial)

end subroutine ufo_geovals_partial_setup_c

!> Allocate GeoVaLs
subroutine ufo_geovals_allocate_c(c_key_self, c_nlevels, c_vars) bind(c,name='ufo_geovals_allocate_f90')
use oops_variables_mod
implicit none
integer(c_int), intent(inout)  :: c_key_self
integer(c_int), intent(in)     :: c_nlevels
type(c_ptr), value, intent(in) :: c_vars

type(ufo_geovals), pointer :: self
type(oops_variables) :: vars

call ufo_geovals_registry%get(c_key_self, self)

vars = oops_variables(c_vars)
call ufo_geovals_allocate(self, vars, c_nlevels)

end subroutine ufo_geovals_allocate_c

! ------------------------------------------------------------------------------
!> Specify which interpolation paths produced by a given method sample which observation locations.
subroutine ufo_geovals_setup_sampling_method_c(c_key_self, c_sampling_method, &
                                               c_npaths, c_nlocs, c_paths_by_loc) &
                                               bind(c,name='ufo_geovals_setup_sampling_method_f90')
implicit none
integer(c_int), intent(inout)     :: c_key_self
integer(c_size_t), intent(in)     :: c_sampling_method, c_npaths, c_nlocs
type(ufo_index_range), intent(in) :: c_paths_by_loc(c_nlocs)

type(ufo_geovals), pointer :: self
integer                    :: npaths, nlocs, sampling_method
type(ufo_index_range)      :: paths_by_loc(c_nlocs)
integer(c_size_t)          :: i

call ufo_geovals_registry%get(c_key_self, self)

npaths = c_npaths
nlocs = c_nlocs
! Convert from 0- to 1-based indexing
sampling_method = 1 + c_sampling_method
do i = 1, c_nlocs
  paths_by_loc(i)%begin = 1 + c_paths_by_loc(i)%begin
  paths_by_loc(i)%end = 1 + c_paths_by_loc(i)%end
enddo

call ufo_geovals_setup_sampling_method(self, sampling_method, npaths, nlocs, paths_by_loc)

end subroutine ufo_geovals_setup_sampling_method_c

! ------------------------------------------------------------------------------
subroutine ufo_geovals_add_reduced_vars_c(c_key_self, c_vars, c_nvars, c_sizes) &
  bind(c,name='ufo_geovals_add_reduced_vars_f90')
use oops_variables_mod
implicit none
integer(c_int), intent(inout)  :: c_key_self
type(c_ptr), value, intent(in) :: c_vars
integer(c_int), intent(in)     :: c_nvars
integer(c_size_t), intent(in)  :: c_sizes(c_nvars)

type(ufo_geovals), pointer :: self
type(oops_variables) :: vars

call ufo_geovals_registry%get(c_key_self, self)

vars = oops_variables(c_vars)
call ufo_geovals_add_reduced_vars(self, vars, c_nvars, c_sizes)

end subroutine ufo_geovals_add_reduced_vars_c

! ------------------------------------------------------------------------------
subroutine ufo_geovals_get_vars_c(c_key_self, c_vars, c_format) &
  bind(c,name='ufo_geovals_get_vars_f90')
use oops_variables_mod
implicit none
integer(c_int), intent(in)     :: c_key_self
type(c_ptr), value, intent(in) :: c_vars
integer(c_int), intent(in)     :: c_format

type(ufo_geovals), pointer :: self
type(oops_variables)       :: vars

call ufo_geovals_registry%get(c_key_self, self)

vars = oops_variables(c_vars)
call ufo_geovals_get_vars(self, vars, c_format)

end subroutine ufo_geovals_get_vars_c

! ------------------------------------------------------------------------------
!> Designate an observation location sampling method as "trivial", i.e. one producing a set of
!> interpolation paths such that each location is sampled solely by the path with the same index.
subroutine ufo_geovals_setup_trivial_sampling_method_c(c_key_self, c_sampling_method) &
  bind(c,name='ufo_geovals_setup_trivial_sampling_method_f90')
implicit none
integer(c_int), intent(inout)     :: c_key_self
integer(c_size_t), intent(in)     :: c_sampling_method

type(ufo_geovals), pointer :: self
integer                    :: sampling_method

call ufo_geovals_registry%get(c_key_self, self)

! Convert from 0- to 1-based indexing
sampling_method = 1 + c_sampling_method

call ufo_geovals_setup_trivial_sampling_method(self, sampling_method)

end subroutine ufo_geovals_setup_trivial_sampling_method_c

! ------------------------------------------------------------------------------
!> Copy one GeoVaLs object into another

subroutine ufo_geovals_copy_c(c_key_self, c_key_other) bind(c,name='ufo_geovals_copy_f90')
implicit none
integer(c_int), intent(in)    :: c_key_self
integer(c_int), intent(inout) :: c_key_other
type(ufo_geovals), pointer    :: self
type(ufo_geovals), pointer    :: other

call ufo_geovals_registry%init()
call ufo_geovals_registry%add(c_key_other)
call ufo_geovals_registry%get(c_key_self, self)
call ufo_geovals_registry%get(c_key_other, other)

call ufo_geovals_copy(self, other)

end subroutine ufo_geovals_copy_c

! ------------------------------------------------------------------------------
!> Copy one GeoVaLs location into another object

subroutine ufo_geovals_copy_one_c(c_key_self, c_key_other, c_ind) bind(c,name='ufo_geovals_copy_one_f90')
implicit none
integer(c_int), intent(inout) :: c_key_self
integer(c_int), intent(in)    :: c_key_other
integer(c_int), intent(in)    :: c_ind
type(ufo_geovals), pointer    :: self
type(ufo_geovals), pointer    :: other
integer :: ind

! Convert location index from the C++ to the Fortran convention.
ind = c_ind + 1

call ufo_geovals_registry%init()
call ufo_geovals_registry%add(c_key_self)
call ufo_geovals_registry%get(c_key_self, self)
call ufo_geovals_registry%get(c_key_other, other)

call ufo_geovals_copy_one(self, other, ind)

end subroutine ufo_geovals_copy_one_c

! ------------------------------------------------------------------------------
!> Analytic init

subroutine ufo_geovals_analytic_init_c(c_key_self, c_locs, c_conf) bind(c,name='ufo_geovals_analytic_init_f90')
use ufo_sampled_locations_mod
implicit none
integer(c_int), intent(in) :: c_key_self
type(c_ptr), value, intent(in) :: c_locs
type(c_ptr), value, intent(in) :: c_conf

type(ufo_geovals), pointer :: self
type(ufo_sampled_locations) ::locs
character(len=30) :: ic
character(len=:), allocatable :: str
type(fckit_configuration) :: f_conf

call ufo_geovals_registry%get(c_key_self, self)

f_conf = fckit_configuration(c_conf)
call f_conf%get_or_die("method",str)
ic = str
locs = ufo_sampled_locations(c_locs)
call ufo_geovals_analytic_init(self,locs,ic)

end subroutine ufo_geovals_analytic_init_c

! ------------------------------------------------------------------------------

subroutine ufo_geovals_delete_c(c_key_self) bind(c,name='ufo_geovals_delete_f90')

implicit none
integer(c_int), intent(inout) :: c_key_self

type(ufo_geovals), pointer :: self

call ufo_geovals_registry%get(c_key_self, self)

call ufo_geovals_delete(self)

call ufo_geovals_registry%remove(c_key_self)

end subroutine ufo_geovals_delete_c

! ------------------------------------------------------------------------------

integer(c_int) function ufo_geovals_get_default_format_c(c_key_self) &
  bind(c,name='ufo_geovals_get_default_format_f90')
implicit none
integer(c_int), intent(in) :: c_key_self

type(ufo_geovals), pointer :: self

call ufo_geovals_registry%get(c_key_self, self)

ufo_geovals_get_default_format_c = ufo_geovals_get_default_format(self)

end function ufo_geovals_get_default_format_c

! ------------------------------------------------------------------------------

subroutine ufo_geovals_set_default_format_c(c_key_self, c_format) &
  bind(c,name='ufo_geovals_set_default_format_f90')
implicit none
integer(c_int), intent(inout) :: c_key_self
integer(c_int), intent(in) :: c_format

type(ufo_geovals), pointer :: self

call ufo_geovals_registry%get(c_key_self, self)

call ufo_geovals_set_default_format(self, c_format)

end subroutine ufo_geovals_set_default_format_c

! ------------------------------------------------------------------------------

logical(c_bool) function ufo_geovals_are_reduced_and_sampled_formats_aliased_c(c_key_self, lvar, c_var) &
  bind(c,name='ufo_geovals_are_reduced_and_sampled_formats_aliased_f90')
use ufo_vars_mod, only: MAXVARLEN
use string_f_c_mod
implicit none
integer(c_int), intent(in) :: c_key_self
integer(c_int), intent(in) :: lvar
character(kind=c_char, len=1), intent(in) :: c_var(lvar+1)

type(ufo_geovals), pointer :: self
character(len=MAXVARLEN) :: varname

call ufo_geovals_registry%get(c_key_self, self)
call c_f_string(c_var, varname)

ufo_geovals_are_reduced_and_sampled_formats_aliased_c = &
  ufo_geovals_are_reduced_and_sampled_formats_aliased(self, varname)

end function ufo_geovals_are_reduced_and_sampled_formats_aliased_c

! ------------------------------------------------------------------------------

subroutine ufo_geovals_zero_c(c_key_self) bind(c,name='ufo_geovals_zero_f90')
implicit none
integer(c_int), intent(in) :: c_key_self
type(ufo_geovals), pointer :: self

call ufo_geovals_registry%get(c_key_self, self)

call ufo_geovals_zero(self)

end subroutine ufo_geovals_zero_c

! ------------------------------------------------------------------------------

subroutine ufo_geovals_reorderzdir_c(c_key_self, lvar, c_var, lvar1, c_var1) bind(c,name='ufo_geovals_reorderzdir_f90')
use ufo_vars_mod, only: MAXVARLEN
use string_f_c_mod
implicit none
integer(c_int), intent(in) :: c_key_self
integer(c_int), intent(in) :: lvar
character(kind=c_char, len=1), intent(in) :: c_var(lvar+1)
integer(c_int), intent(in) :: lvar1
character(kind=c_char, len=1), intent(in) :: c_var1(lvar1+1)
character(len=MAXVARLEN) :: varname
character(len=MAXVARLEN) :: vardir
type(ufo_geovals), pointer :: self

call c_f_string(c_var, varname)
call c_f_string(c_var1, vardir)
call ufo_geovals_registry%get(c_key_self, self)

call ufo_geovals_reorderzdir(self, varname, vardir)

end subroutine ufo_geovals_reorderzdir_c

! ------------------------------------------------------------------------------

subroutine ufo_geovals_abs_c(c_key_self) bind(c,name='ufo_geovals_abs_f90')
implicit none
integer(c_int), intent(in) :: c_key_self
type(ufo_geovals), pointer :: self

call ufo_geovals_registry%get(c_key_self, self)

call ufo_geovals_abs(self)

end subroutine ufo_geovals_abs_c

! ------------------------------------------------------------------------------

subroutine ufo_geovals_rms_c(c_key_self,vrms) bind(c,name='ufo_geovals_rms_f90')
implicit none
integer(c_int), intent(in) :: c_key_self
real(c_double), intent(inout) :: vrms
type(ufo_geovals), pointer :: self

call ufo_geovals_registry%get(c_key_self, self)

call ufo_geovals_rms(self,vrms)

end subroutine ufo_geovals_rms_c

! ------------------------------------------------------------------------------

subroutine ufo_geovals_random_c(c_key_self) bind(c,name='ufo_geovals_random_f90')
implicit none
integer(c_int), intent(in) :: c_key_self
type(ufo_geovals), pointer :: self

call ufo_geovals_registry%get(c_key_self, self)

call ufo_geovals_random(self)

end subroutine ufo_geovals_random_c

! ------------------------------------------------------------------------------

subroutine ufo_geovals_scalmult_c(c_key_self, zz) bind(c,name='ufo_geovals_scalmult_f90')
implicit none
integer(c_int), intent(in) :: c_key_self
real(c_double), intent(in) :: zz
type(ufo_geovals), pointer :: self

call ufo_geovals_registry%get(c_key_self, self)

call ufo_geovals_scalmult(self, zz)

end subroutine ufo_geovals_scalmult_c

! ------------------------------------------------------------------------------

subroutine ufo_geovals_profmult_c(c_key_self, nlocs, values) bind(c,name='ufo_geovals_profmult_f90')
implicit none
integer(c_int), intent(in) :: c_key_self
integer(c_int), intent(in) :: nlocs
real(c_float), intent(in) :: values(nlocs)
type(ufo_geovals), pointer :: self

call ufo_geovals_registry%get(c_key_self, self)

call ufo_geovals_profmult(self, nlocs, values)

end subroutine ufo_geovals_profmult_c

! ------------------------------------------------------------------------------

subroutine ufo_geovals_assign_c(c_key_self, c_key_rhs) bind(c,name='ufo_geovals_assign_f90')
implicit none
integer(c_int), intent(in) :: c_key_self
integer(c_int), intent(in) :: c_key_rhs
type(ufo_geovals), pointer :: self
type(ufo_geovals), pointer :: rhs

call ufo_geovals_registry%get(c_key_self, self)
call ufo_geovals_registry%get(c_key_rhs, rhs)

call ufo_geovals_assign(self, rhs)

end subroutine ufo_geovals_assign_c

! ------------------------------------------------------------------------------

subroutine ufo_geovals_add_c(c_key_self, c_key_other) bind(c,name='ufo_geovals_add_f90')
implicit none
integer(c_int), intent(in) :: c_key_self
integer(c_int), intent(in) :: c_key_other
type(ufo_geovals), pointer :: self
type(ufo_geovals), pointer :: other

call ufo_geovals_registry%get(c_key_self, self)
call ufo_geovals_registry%get(c_key_other, other)

call ufo_geovals_add(self, other)

end subroutine ufo_geovals_add_c

! ------------------------------------------------------------------------------

subroutine ufo_geovals_diff_c(c_key_self, c_key_other) bind(c,name='ufo_geovals_diff_f90')
implicit none
integer(c_int), intent(in) :: c_key_self
integer(c_int), intent(in) :: c_key_other
type(ufo_geovals), pointer :: self
type(ufo_geovals), pointer :: other

call ufo_geovals_registry%get(c_key_self, self)
call ufo_geovals_registry%get(c_key_other, other)

call ufo_geovals_diff(self, other)

end subroutine ufo_geovals_diff_c

! ------------------------------------------------------------------------------

subroutine ufo_geovals_schurmult_c(c_key_self, c_key_other) bind(c,name='ufo_geovals_schurmult_f90')
implicit none
integer(c_int), intent(in) :: c_key_self
integer(c_int), intent(in) :: c_key_other
type(ufo_geovals), pointer :: self
type(ufo_geovals), pointer :: other

call ufo_geovals_registry%get(c_key_self, self)
call ufo_geovals_registry%get(c_key_other, other)

call ufo_geovals_schurmult(self, other)

end subroutine ufo_geovals_schurmult_c

! ------------------------------------------------------------------------------

subroutine ufo_geovals_normalize_c(c_key_self, c_key_other) bind(c,name='ufo_geovals_normalize_f90')
implicit none
integer(c_int), intent(in) :: c_key_self
integer(c_int), intent(in) :: c_key_other
type(ufo_geovals), pointer :: self
type(ufo_geovals), pointer :: other

call ufo_geovals_registry%get(c_key_self, self)
call ufo_geovals_registry%get(c_key_other, other)

call ufo_geovals_normalize(self, other)

end subroutine ufo_geovals_normalize_c

! ------------------------------------------------------------------------------

subroutine ufo_geovals_minmaxavg_c(c_key_self, kobs, kvar, pmin, pmax, prms) bind(c,name='ufo_geovals_minmaxavg_f90')
implicit none
integer(c_int), intent(in) :: c_key_self
integer(c_int), intent(inout) :: kobs
integer(c_int), intent(in) :: kvar
real(c_double), intent(inout) :: pmin, pmax, prms
type(ufo_geovals), pointer :: self

call ufo_geovals_registry%get(c_key_self, self)

call ufo_geovals_minmaxavg(self, kobs, kvar, pmin, pmax, prms)

end subroutine ufo_geovals_minmaxavg_c

! ------------------------------------------------------------------------------

subroutine ufo_geovals_nlocs_c(c_key_self, nlocs) bind(c, name='ufo_geovals_nlocs_f90')
use ufo_vars_mod, only: MAXVARLEN
use string_f_c_mod
implicit none
integer(c_int), intent(in) :: c_key_self
integer(c_size_t), intent(out) :: nlocs

type(ufo_geovals), pointer :: self

call ufo_geovals_registry%get(c_key_self, self)

nlocs = self%nlocs

end subroutine ufo_geovals_nlocs_c

! ------------------------------------------------------------------------------

subroutine ufo_geovals_nprofiles_c(c_key_self, lvar, c_var, c_format, nprofiles) &
  bind(c, name='ufo_geovals_nprofiles_f90')
use ufo_vars_mod, only: MAXVARLEN
use string_f_c_mod
implicit none
integer(c_int), intent(in) :: c_key_self
integer(c_int), intent(in) :: lvar
character(kind=c_char, len=1), intent(in) :: c_var(lvar+1)
integer(c_int), intent(in) :: c_format
integer(c_size_t), intent(out) :: nprofiles

type(ufo_geoval), pointer :: geoval
character(len=MAXVARLEN) :: varname
type(ufo_geovals), pointer :: self

call c_f_string(c_var, varname)
call ufo_geovals_registry%get(c_key_self, self)

call ufo_geovals_get_var(self, varname, geoval, c_format)

nprofiles = geoval%nprofiles

end subroutine ufo_geovals_nprofiles_c

! ------------------------------------------------------------------------------

subroutine ufo_geovals_nlevs_c(c_key_self, lvar, c_var, c_format, nlevs) &
  bind(c, name='ufo_geovals_nlevs_f90')
use ufo_vars_mod, only: MAXVARLEN
use string_f_c_mod
implicit none
integer(c_int), intent(in) :: c_key_self
integer(c_int), intent(in) :: lvar
character(kind=c_char, len=1), intent(in) :: c_var(lvar+1)
integer(c_int), intent(in) :: c_format
integer(c_int), intent(out) :: nlevs

type(ufo_geoval), pointer :: geoval
character(len=MAXVARLEN) :: varname
type(ufo_geovals), pointer :: self

call c_f_string(c_var, varname)
call ufo_geovals_registry%get(c_key_self, self)

call ufo_geovals_get_var(self, varname, geoval, c_format)

nlevs = geoval%nval

end subroutine ufo_geovals_nlevs_c

! ------------------------------------------------------------------------------

subroutine ufo_geovals_get2d_c(c_key_self, lvar, c_var, c_format, nprofiles, values) &
  bind(c, name='ufo_geovals_get2d_f90')
use ufo_vars_mod, only: MAXVARLEN
use string_f_c_mod
implicit none
integer(c_int),    intent(in) :: c_key_self
integer(c_int),    intent(in) :: lvar
character(kind=c_char, len=1), intent(in) :: c_var(lvar+1)
integer(c_int),    intent(in) :: c_format
integer(c_int),    intent(in) :: nprofiles
real(c_double), intent(inout) :: values(nprofiles)

character(max_string) :: err_msg
type(ufo_geoval), pointer :: geoval
character(len=MAXVARLEN) :: varname
type(ufo_geovals), pointer :: self

call c_f_string(c_var, varname)
call ufo_geovals_registry%get(c_key_self, self)

call ufo_geovals_get_var(self, varname, geoval, c_format)

if (size(geoval%vals,1) /= 1) then
  write(err_msg,*)'ufo_geovals_get2d_f90',trim(varname),'is not a 2D var:',size(geoval%vals,1), ' levels'
  call abor1_ftn(err_msg)
endif
if (nprofiles /= size(geoval%vals,2)) then
  write(err_msg,*)'ufo_geovals_get2d_f90',trim(varname),'error profiles number:',nprofiles,size(geoval%vals,2)
  call abor1_ftn(err_msg)
endif

values(:) = geoval%vals(1,:)

end subroutine ufo_geovals_get2d_c

! ------------------------------------------------------------------------------

subroutine ufo_geovals_get_c(c_key_self, lvar, c_var, c_format, c_lev, nprofiles, values) &
  bind(c, name='ufo_geovals_get_f90')
use ufo_vars_mod, only: MAXVARLEN
use string_f_c_mod
implicit none
integer(c_int), intent(in) :: c_key_self
integer(c_int), intent(in) :: lvar
character(kind=c_char, len=1), intent(in) :: c_var(lvar+1)
integer(c_int), intent(in) :: c_format
integer(c_int), intent(in) :: c_lev
integer(c_int), intent(in) :: nprofiles
real(c_float), intent(inout) :: values(nprofiles)

character(max_string) :: err_msg
type(ufo_geoval), pointer :: geoval
character(len=MAXVARLEN) :: varname
type(ufo_geovals), pointer :: self
integer(c_int) :: lev

! Convert level index from the C++ to the Fortran convention.
lev = c_lev + 1

call c_f_string(c_var, varname)
call ufo_geovals_registry%get(c_key_self, self)

call ufo_geovals_get_var(self, varname, geoval, c_format)

if (lev<1 .or. lev>size(geoval%vals,1)) then
  write(err_msg,*)'ufo_geovals_get_f90 "',trim(varname),'" level out of range: 1~', &
                  size(geoval%vals,1), ', lev=', lev
  call abor1_ftn(err_msg)
endif
if (nprofiles /= size(geoval%vals,2)) then
  write(err_msg,*)'ufo_geovals_get_f90 "',trim(varname),'" error profiles number:',nprofiles,&
                  ' /= ',size(geoval%vals,2)
  call abor1_ftn(err_msg)
endif

values(:) = geoval%vals(lev,:)

end subroutine ufo_geovals_get_c

! ------------------------------------------------------------------------------

subroutine ufo_geovals_get_profile_c(c_key_self, lvar, c_var, c_format, c_profile, nlevs, values) &
  bind(c, name='ufo_geovals_get_profile_f90')
use ufo_vars_mod, only: MAXVARLEN
use string_f_c_mod
implicit none
integer(c_int), intent(in) :: c_key_self
integer(c_int), intent(in) :: lvar
character(kind=c_char, len=1), intent(in) :: c_var(lvar+1)
integer(c_int), intent(in) :: c_format
integer(c_int), intent(in) :: c_profile
integer(c_int), intent(in) :: nlevs
real(c_double), intent(inout) :: values(nlevs)

character(max_string) :: err_msg
type(ufo_geoval), pointer :: geoval
character(len=MAXVARLEN) :: varname
type(ufo_geovals), pointer :: self
integer(c_int) :: profile

call c_f_string(c_var, varname)
call ufo_geovals_registry%get(c_key_self, self)

call ufo_geovals_get_var(self, varname, geoval, c_format)

! Convert profile index from the C++ to the Fortran convention.
profile = c_profile + 1

if (profile<1 .or. profile>size(geoval%vals,2)) then
  write(err_msg,*)'ufo_geovals_get_profile_f90',trim(varname),'profile index out of range:',profile,size(geoval%vals,2)
  call abor1_ftn(err_msg)
endif
if (nlevs /= size(geoval%vals,1)) then
  write(err_msg,*)'ufo_geovals_get_profile_f90',trim(varname),'incorrect number of levels:',nlevs,size(geoval%vals,1)
  call abor1_ftn(err_msg)
endif

values(:) = geoval%vals(:,profile)

end subroutine ufo_geovals_get_profile_c

! ------------------------------------------------------------------------------

subroutine ufo_geovals_getdouble_c(c_key_self, lvar, c_var, c_format, c_lev, nprofiles, values) &
  bind(c, name='ufo_geovals_getdouble_f90')
use ufo_vars_mod, only: MAXVARLEN
use string_f_c_mod
implicit none
integer(c_int), intent(in) :: c_key_self
integer(c_int), intent(in) :: lvar
character(kind=c_char, len=1), intent(in) :: c_var(lvar+1)
integer(c_int), intent(in) :: c_format
integer(c_int), intent(in) :: c_lev
integer(c_int), intent(in) :: nprofiles
real(c_double), intent(inout) :: values(nprofiles)

type(ufo_geoval), pointer :: geoval
character(len=MAXVARLEN) :: varname
type(ufo_geovals), pointer :: self
integer(c_int) :: lev

! Convert level index from the C++ to the Fortran convention.
lev = c_lev + 1

call c_f_string(c_var, varname)
call ufo_geovals_registry%get(c_key_self, self)
call ufo_geovals_get_var(self, varname, geoval, c_format)
values(:) = geoval%vals(lev,:)

end subroutine ufo_geovals_getdouble_c

! ------------------------------------------------------------------------------

subroutine ufo_geovals_putdouble_c(c_key_self, lvar, c_var, c_format, c_lev, nprofiles, values) &
  bind(c, name='ufo_geovals_putdouble_f90')
use ufo_vars_mod, only: MAXVARLEN
use string_f_c_mod
integer(c_int), intent(in) :: c_key_self
integer(c_int), intent(in) :: lvar
character(kind=c_char, len=1), intent(in) :: c_var(lvar+1)
integer(c_int), intent(in) :: c_format
integer(c_int), intent(in) :: c_lev
integer(c_int), intent(in) :: nprofiles
real(c_double), intent(in) :: values(nprofiles)

type(ufo_geoval), pointer  :: geoval
character(len=MAXVARLEN)   :: varname
type(ufo_geovals), pointer :: self
integer(c_int) :: lev

! Convert level index from the C++ to the Fortran convention.
lev = c_lev + 1

call c_f_string(c_var, varname)
call ufo_geovals_registry%get(c_key_self, self)
call ufo_geovals_get_var(self, varname, geoval, c_format)
geoval%vals(lev,:) = values(:)

end subroutine ufo_geovals_putdouble_c

! ------------------------------------------------------------------------------
subroutine ufo_geovals_put_profile_c(c_key_self, lvar, c_var, c_format, c_profile, nlevs, values) &
  bind(c, name='ufo_geovals_put_profile_f90')
use ufo_vars_mod, only: MAXVARLEN
use string_f_c_mod

implicit none
integer(c_int), intent(in) :: c_key_self
integer(c_int), intent(in) :: lvar
character(kind=c_char, len=1), intent(in) :: c_var(lvar+1)
integer(c_int), intent(in) :: c_format
integer(c_int), intent(in) :: c_profile
integer(c_int), intent(in) :: nlevs
real(c_double), intent(in) :: values(nlevs)

character(max_string) :: err_msg
type(ufo_geoval), pointer :: geoval
character(len=MAXVARLEN) :: varname
type(ufo_geovals), pointer :: self
integer(c_int) :: profile

call c_f_string(c_var, varname)
call ufo_geovals_registry%get(c_key_self, self)
call ufo_geovals_get_var(self, varname, geoval, c_format)

! Convert profile index from the C++ to the Fortran convention.
profile = c_profile + 1

if (profile<1 .or. profile>size(geoval%vals,2)) then
  write(err_msg,*)'ufo_geovals_put_profile_f90',trim(varname),'profile out of range:',profile,size(geoval%vals,2)
  call abor1_ftn(err_msg)
endif
if (nlevs /= size(geoval%vals,1)) then
  write(err_msg,*)'ufo_geovals_put_profile_f90',trim(varname),'incorrect number of levels:',nlevs,size(geoval%vals,1)
  call abor1_ftn(err_msg)
endif

geoval%vals(:,profile) = values(:)

end subroutine ufo_geovals_put_profile_c

! ------------------------------------------------------------------------------
subroutine ufo_geovals_get_profile_indices_grouped_by_loc_c( &
  c_key_self, lvar, c_var, c_format, c_nlocs, c_profile_indices_grouped_by_loc) &
  bind(c, name='ufo_geovals_get_profile_indices_grouped_by_loc_f90')
use ufo_vars_mod, only: ufo_vars_getindex, MAXVARLEN
use string_f_c_mod

implicit none
integer(c_int), intent(in) :: c_key_self
integer(c_int), intent(in) :: lvar
character(kind=c_char, len=1), intent(in) :: c_var(lvar+1)
integer(c_int), intent(in) :: c_format
integer(c_size_t), intent(in) :: c_nlocs
type(ufo_index_range), intent(inout) :: c_profile_indices_grouped_by_loc(c_nlocs)

character(max_string) :: err_msg
character(len=MAXVARLEN) :: varname
type(ufo_geovals), pointer :: self
character(len=MAXVARLEN), pointer :: variables(:)  !< either the sampled or reduced variable list
integer :: ivar, iloc, jv

call c_f_string(c_var, varname)
call ufo_geovals_registry%get(c_key_self, self)

if (c_nlocs /= self%nlocs) then
  write(err_msg,*) 'ufo_geovals_get_profile_indices_grouped_by_loc_f90: nlocs mismatch: received ', &
    c_nlocs, ', expected ', self%nlocs
  call abor1_ftn(err_msg)
endif

if (c_format == ufo_geoval_sampled) then
  variables => self%variables
else
  variables => self%reduced_variables
endif

ivar = ufo_vars_getindex(variables, varname)

if (ivar < 0) then
  write(err_msg,*) 'ufo_geovals_get_profile_indices_grouped_by_loc_f90: ', trim(varname), ' doesn''t exist'
  call abor1_ftn(err_msg)
endif

c_profile_indices_grouped_by_loc(:) = &
  self%sampling_methods(self%sampling_method_by_var(ivar))%paths_by_loc(:)

! Convert profile indices from the Fortran to the C++ convention.
do iloc = 1, self%nlocs
  c_profile_indices_grouped_by_loc(iloc)%begin = c_profile_indices_grouped_by_loc(iloc)%begin - 1
  c_profile_indices_grouped_by_loc(iloc)%end = c_profile_indices_grouped_by_loc(iloc)%end - 1
enddo

end subroutine ufo_geovals_get_profile_indices_grouped_by_loc_c

! ------------------------------------------------------------------------------

subroutine ufo_geovals_maxloc_c(c_key_self, mxval, iloc, ivar) bind(c,name='ufo_geovals_maxloc_f90')
implicit none
integer(c_int), intent(in) :: c_key_self
real(c_double), intent(inout) :: mxval
integer(c_int), intent(inout) :: iloc, ivar
type(ufo_geovals), pointer :: self

call ufo_geovals_registry%get(c_key_self, self)

call ufo_geovals_maxloc(self, mxval, iloc, ivar)

end subroutine ufo_geovals_maxloc_c

! ------------------------------------------------------------------------------

subroutine ufo_geovals_fill_c(c_key, lvar, c_var, c_nloc, c_indx, c_nlev, c_vals, c_levelsTopDown) &
  bind(c, name="ufo_geovals_fill_f90")
use ufo_vars_mod, only: MAXVARLEN
use string_f_c_mod, only: c_f_string
implicit none
integer(c_int), intent(in) :: c_key
integer(c_int), intent(in) :: lvar
character(kind=c_char, len=1), intent(in) :: c_var(lvar+1)
integer(c_int), intent(in) :: c_nloc
integer(c_int), intent(in) :: c_indx(c_nloc)
integer(c_int), intent(in) :: c_nlev
real(c_double), intent(in) :: c_vals(c_nloc, c_nlev)
logical(c_bool), intent(in) :: c_levelsTopDown

type(ufo_geovals), pointer :: geovals
character(len=MAXVARLEN) :: varname

call ufo_geovals_registry%get(c_key, geovals)
call c_f_string(c_var, varname)

call ufo_geovals_fill(geovals, varname, c_nloc, c_indx, c_nlev, c_vals, c_levelsTopDown)

end subroutine ufo_geovals_fill_c

! ------------------------------------------------------------------------------

subroutine ufo_geovals_fillad_c(c_key, lvar, c_var, c_nloc, c_indx, c_nlev, c_vals, c_levelsTopDown) &
  bind(c, name="ufo_geovals_fillad_f90")
use ufo_vars_mod, only: MAXVARLEN
use string_f_c_mod, only: c_f_string
implicit none
integer(c_int), intent(in) :: c_key
integer(c_int), intent(in) :: lvar
character(kind=c_char, len=1), intent(in) :: c_var(lvar+1)
integer(c_int), intent(in) :: c_nloc
integer(c_int), intent(in) :: c_indx(c_nloc)
integer(c_int), intent(in) :: c_nlev
real(c_double), intent(inout) :: c_vals(c_nloc, c_nlev)
logical(c_bool), intent(in) :: c_levelsTopDown

type(ufo_geovals), pointer :: geovals
character(len=MAXVARLEN) :: varname
call c_f_string(c_var, varname)

call ufo_geovals_registry%get(c_key, geovals)

call ufo_geovals_fillad(geovals, varname, c_nloc, c_indx, c_nlev, c_vals, c_levelsTopDown)

end subroutine ufo_geovals_fillad_c

! ------------------------------------------------------------------------------

subroutine ufo_geovals_read_file_c(c_key_self, c_conf, c_obspace, c_vars) &
  bind(c,name='ufo_geovals_read_file_f90')
use oops_variables_mod

implicit none
integer(c_int), intent(inout)  :: c_key_self
type(c_ptr), value, intent(in) :: c_conf
type(c_ptr), value, intent(in) :: c_obspace
type(c_ptr), value, intent(in) :: c_vars

type(ufo_geovals), pointer :: self
character(max_string)      :: filename
integer :: loc_multiplier
logical :: levels_are_top_down
character(len=:), allocatable :: str
type(fckit_configuration) :: f_conf
type(oops_variables)      :: vars

call ufo_geovals_registry%init()
call ufo_geovals_registry%add(c_key_self)
call ufo_geovals_registry%get(c_key_self, self)

! read filename for config
f_conf = fckit_configuration(c_conf)
call f_conf%get_or_die("filename",str)
filename = str

if (f_conf%has("loc_multiplier")) then
  call f_conf%get_or_die("loc_multiplier", loc_multiplier)
else
  loc_multiplier = 1
endif

if (f_conf%has("levels_are_top_down")) then
  call f_conf%get_or_die("levels_are_top_down", levels_are_top_down)
else
  levels_are_top_down = .True.
endif

vars = oops_variables(c_vars)
! read geovals
call ufo_geovals_read_netcdf(self, filename, loc_multiplier, levels_are_top_down, c_obspace, vars)

end subroutine ufo_geovals_read_file_c

! ------------------------------------------------------------------------------

subroutine ufo_geovals_write_file_c(c_key_self, c_conf, c_rank) bind(c,name='ufo_geovals_write_file_f90')
implicit none
integer(c_int), intent(in)     :: c_key_self
type(c_ptr), value, intent(in) :: c_conf
integer(c_size_t), intent(in)  :: c_rank   ! mpi rank (to be added to filename)

type(ufo_geovals), pointer :: self
character(max_string)      :: fout, filename

character(len=10)             :: cproc
integer                       :: ppos
character(len=:), allocatable :: str
type(fckit_configuration)     :: f_conf

! read filename for config
f_conf = fckit_configuration(c_conf)
call f_conf%get_or_die("filename",str)
filename = str

write(cproc,fmt='(i4.4)') c_rank

! Find the left-most dot in the file name, and use that to pick off the file name
! and file extension.
ppos = scan(trim(filename), '.', BACK=.true.)
if (ppos > 0) then
 ! found a file extension
 fout = filename(1:ppos-1) // '_' // trim(adjustl(cproc)) // trim(filename(ppos:))
else
 ! no file extension
 fout = trim(filename) // '_' // trim(adjustl(cproc))
endif

call ufo_geovals_registry%get(c_key_self, self)
call ufo_geovals_write_netcdf(self, fout)

end subroutine ufo_geovals_write_file_c

! ------------------------------------------------------------------------------

end module ufo_geovals_mod_c
