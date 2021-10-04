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

subroutine ufo_geovals_setup_c(c_key_self, c_nlocs, c_vars, c_nvars, c_sizes) bind(c,name='ufo_geovals_setup_f90')
use oops_variables_mod
implicit none
integer(c_int), intent(inout)  :: c_key_self
integer(c_int), intent(in)     :: c_nlocs, c_nvars
type(c_ptr), value, intent(in) :: c_vars
integer(c_size_t), intent(in)     :: c_sizes(c_nvars)

type(ufo_geovals), pointer :: self
type(oops_variables) :: vars

call ufo_geovals_registry%init()
call ufo_geovals_registry%add(c_key_self)
call ufo_geovals_registry%get(c_key_self, self)

vars = oops_variables(c_vars)
call ufo_geovals_setup(self, vars, c_nlocs, c_nvars, c_sizes)

end subroutine ufo_geovals_setup_c

!> Setup GeoVaLs (store nlocs, variables; don't do allocation yet)
subroutine ufo_geovals_partial_setup_c(c_key_self, c_nlocs, c_vars) bind(c,name='ufo_geovals_partial_setup_f90')
use oops_variables_mod
implicit none
integer(c_int), intent(inout)  :: c_key_self
integer(c_int), intent(in)     :: c_nlocs
type(c_ptr), value, intent(in) :: c_vars

type(ufo_geovals), pointer :: self
type(oops_variables) :: vars

call ufo_geovals_registry%init()
call ufo_geovals_registry%add(c_key_self)
call ufo_geovals_registry%get(c_key_self, self)

vars = oops_variables(c_vars)
call ufo_geovals_partial_setup(self, vars, c_nlocs)

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
use ufo_locations_mod
implicit none
integer(c_int), intent(in) :: c_key_self
type(c_ptr), value, intent(in) :: c_locs
type(c_ptr), value, intent(in) :: c_conf

type(ufo_geovals), pointer :: self
type(ufo_locations) ::locs
character(len=30) :: ic
character(len=:), allocatable :: str
type(fckit_configuration) :: f_conf

call ufo_geovals_registry%get(c_key_self, self)

f_conf = fckit_configuration(c_conf)
call f_conf%get_or_die("method",str)
ic = str
locs = ufo_locations(c_locs)
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

subroutine ufo_geovals_split_c(c_key_self, c_key_other1, c_key_other2) bind(c,name='ufo_geovals_split_f90')
implicit none
integer(c_int), intent(in) :: c_key_self, c_key_other1, c_key_other2
type(ufo_geovals), pointer :: self, other1, other2

call ufo_geovals_registry%get(c_key_self, self)
call ufo_geovals_registry%get(c_key_other1, other1)
call ufo_geovals_registry%get(c_key_other2, other2)

call ufo_geovals_split(self, other1, other2)

end subroutine ufo_geovals_split_c

! ------------------------------------------------------------------------------

subroutine ufo_geovals_merge_c(c_key_self, c_key_other1, c_key_other2) bind(c,name='ufo_geovals_merge_f90')
implicit none
integer(c_int), intent(in) :: c_key_self, c_key_other1, c_key_other2
type(ufo_geovals), pointer :: self, other1, other2

call ufo_geovals_registry%get(c_key_self, self)
call ufo_geovals_registry%get(c_key_other1, other1)
call ufo_geovals_registry%get(c_key_other2, other2)

call ufo_geovals_merge(self, other1, other2)

end subroutine ufo_geovals_merge_c

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

subroutine ufo_geovals_nlocs_c(c_key_self, kobs) bind(c, name='ufo_geovals_nlocs_f90')
implicit none
integer(c_int),    intent(in) :: c_key_self
integer(c_size_t), intent(inout) :: kobs
type(ufo_geovals), pointer :: self

call ufo_geovals_registry%get(c_key_self, self)
kobs = self%nlocs

end subroutine ufo_geovals_nlocs_c

! ------------------------------------------------------------------------------

subroutine ufo_geovals_nlevs_c(c_key_self, lvar, c_var, nlevs) bind(c, name='ufo_geovals_nlevs_f90')
use ufo_vars_mod, only: MAXVARLEN
use string_f_c_mod
implicit none
integer(c_int), intent(in) :: c_key_self
integer(c_int), intent(in) :: lvar
character(kind=c_char, len=1), intent(in) :: c_var(lvar+1)
integer(c_int), intent(out) :: nlevs

type(ufo_geoval), pointer :: geoval
character(len=MAXVARLEN) :: varname
type(ufo_geovals), pointer :: self

call c_f_string(c_var, varname)
call ufo_geovals_registry%get(c_key_self, self)

call ufo_geovals_get_var(self, varname, geoval)

nlevs = geoval%nval

end subroutine ufo_geovals_nlevs_c

! ------------------------------------------------------------------------------

subroutine ufo_geovals_get2d_c(c_key_self, lvar, c_var, nlocs, values) bind(c, name='ufo_geovals_get2d_f90')
use ufo_vars_mod, only: MAXVARLEN
use string_f_c_mod
implicit none
integer(c_int),    intent(in) :: c_key_self
integer(c_int),    intent(in) :: lvar
character(kind=c_char, len=1), intent(in) :: c_var(lvar+1)
integer(c_int),    intent(in) :: nlocs
real(c_double), intent(inout) :: values(nlocs)

character(max_string) :: err_msg
type(ufo_geoval), pointer :: geoval
character(len=MAXVARLEN) :: varname
type(ufo_geovals), pointer :: self

call c_f_string(c_var, varname)
call ufo_geovals_registry%get(c_key_self, self)

call ufo_geovals_get_var(self, varname, geoval)

if (size(geoval%vals,1) /= 1) then
  write(err_msg,*)'ufo_geovals_get2d_f90',trim(varname),'is not a 2D var:',size(geoval%vals,1), ' levels'
  call abor1_ftn(err_msg)
endif
if (nlocs /= size(geoval%vals,2)) then
  write(err_msg,*)'ufo_geovals_get2d_f90',trim(varname),'error locs number:',nlocs,size(geoval%vals,2)
  call abor1_ftn(err_msg)
endif

values(:) = geoval%vals(1,:)

end subroutine ufo_geovals_get2d_c

! ------------------------------------------------------------------------------

subroutine ufo_geovals_get_c(c_key_self, lvar, c_var, c_lev, nlocs, values) bind(c, name='ufo_geovals_get_f90')
use ufo_vars_mod, only: MAXVARLEN
use string_f_c_mod
implicit none
integer(c_int), intent(in) :: c_key_self
integer(c_int), intent(in) :: lvar
character(kind=c_char, len=1), intent(in) :: c_var(lvar+1)
integer(c_int), intent(in) :: c_lev
integer(c_int), intent(in) :: nlocs
real(c_float), intent(inout) :: values(nlocs)

character(max_string) :: err_msg
type(ufo_geoval), pointer :: geoval
character(len=MAXVARLEN) :: varname
type(ufo_geovals), pointer :: self
integer(c_int) :: lev

! Convert level index from the C++ to the Fortran convention.
lev = c_lev + 1

call c_f_string(c_var, varname)
call ufo_geovals_registry%get(c_key_self, self)

call ufo_geovals_get_var(self, varname, geoval)

if (lev<1 .or. lev>size(geoval%vals,1)) then
  write(err_msg,*)'ufo_geovals_get_f90 "',trim(varname),'" level out of range: 1~', &
                  size(geoval%vals,1), ', lev=', lev
  call abor1_ftn(err_msg)
endif
if (nlocs /= size(geoval%vals,2)) then
  write(err_msg,*)'ufo_geovals_get_f90 "',trim(varname),'" error locs number:',nlocs,&
                  ' /= ',size(geoval%vals,2)
  call abor1_ftn(err_msg)
endif

values(:) = geoval%vals(lev,:)

end subroutine ufo_geovals_get_c

! ------------------------------------------------------------------------------

subroutine ufo_geovals_get_loc_c(c_key_self, lvar, c_var, c_loc, nlevs, values) bind(c, name='ufo_geovals_get_loc_f90')
use ufo_vars_mod, only: MAXVARLEN
use string_f_c_mod
implicit none
integer(c_int), intent(in) :: c_key_self
integer(c_int), intent(in) :: lvar
character(kind=c_char, len=1), intent(in) :: c_var(lvar+1)
integer(c_int), intent(in) :: c_loc
integer(c_int), intent(in) :: nlevs
real(c_double), intent(inout) :: values(nlevs)

character(max_string) :: err_msg
type(ufo_geoval), pointer :: geoval
character(len=MAXVARLEN) :: varname
type(ufo_geovals), pointer :: self
integer(c_int) :: loc

call c_f_string(c_var, varname)
call ufo_geovals_registry%get(c_key_self, self)

call ufo_geovals_get_var(self, varname, geoval)

! Convert location index from the C++ to the Fortran convention.
loc = c_loc + 1

if (loc<1 .or. loc>size(geoval%vals,2)) then
  write(err_msg,*)'ufo_geovals_get_loc_f90',trim(varname),'location out of range:',loc,size(geoval%vals,2)
  call abor1_ftn(err_msg)
endif
if (nlevs /= size(geoval%vals,1)) then
  write(err_msg,*)'ufo_geovals_get_loc_f90',trim(varname),'incorrect number of levels:',nlevs,size(geoval%vals,1)
  call abor1_ftn(err_msg)
endif

values(:) = geoval%vals(:,loc)

end subroutine ufo_geovals_get_loc_c

! ------------------------------------------------------------------------------

subroutine ufo_geovals_getdouble_c(c_key_self, lvar, c_var, c_lev, nlocs, values)&
  bind(c, name='ufo_geovals_getdouble_f90')
use ufo_vars_mod, only: MAXVARLEN
use string_f_c_mod
implicit none
integer(c_int), intent(in) :: c_key_self
integer(c_int), intent(in) :: lvar
character(kind=c_char, len=1), intent(in) :: c_var(lvar+1)
integer(c_int), intent(in) :: c_lev
integer(c_int), intent(in) :: nlocs
real(c_double), intent(inout) :: values(nlocs)

type(ufo_geoval), pointer :: geoval
character(len=MAXVARLEN) :: varname
type(ufo_geovals), pointer :: self
integer(c_int) :: lev

! Convert level index from the C++ to the Fortran convention.
lev = c_lev + 1

call c_f_string(c_var, varname)
call ufo_geovals_registry%get(c_key_self, self)
call ufo_geovals_get_var(self, varname, geoval)
values(:) = geoval%vals(lev,:)

end subroutine ufo_geovals_getdouble_c

! ------------------------------------------------------------------------------

subroutine ufo_geovals_putdouble_c(c_key_self, lvar, c_var, c_lev, nlocs, values) bind(c, name='ufo_geovals_putdouble_f90')
use ufo_vars_mod, only: MAXVARLEN
use string_f_c_mod
integer(c_int), intent(in) :: c_key_self
integer(c_int), intent(in) :: lvar
character(kind=c_char, len=1), intent(in) :: c_var(lvar+1)
integer(c_int), intent(in) :: c_lev
integer(c_int), intent(in) :: nlocs
real(c_double), intent(in) :: values(nlocs)

type(ufo_geoval), pointer  :: geoval
character(len=MAXVARLEN)   :: varname
type(ufo_geovals), pointer :: self
integer(c_int) :: lev

! Convert level index from the C++ to the Fortran convention.
lev = c_lev + 1

call c_f_string(c_var, varname)
call ufo_geovals_registry%get(c_key_self, self)
call ufo_geovals_get_var(self, varname, geoval)
geoval%vals(lev,:) = values(:)

end subroutine ufo_geovals_putdouble_c

! ------------------------------------------------------------------------------
subroutine ufo_geovals_put_loc_c(c_key_self, lvar, c_var, c_loc, nlevs, values) bind(c, name='ufo_geovals_put_loc_f90')
use ufo_vars_mod, only: MAXVARLEN
use string_f_c_mod

implicit none
integer(c_int), intent(in) :: c_key_self
integer(c_int), intent(in) :: lvar
character(kind=c_char, len=1), intent(in) :: c_var(lvar+1)
integer(c_int), intent(in) :: c_loc
integer(c_int), intent(in) :: nlevs
real(c_double), intent(in) :: values(nlevs)

character(max_string) :: err_msg
type(ufo_geoval), pointer :: geoval
character(len=MAXVARLEN) :: varname
type(ufo_geovals), pointer :: self
integer(c_int) :: loc

call c_f_string(c_var, varname)
call ufo_geovals_registry%get(c_key_self, self)
call ufo_geovals_get_var(self, varname, geoval)

! Convert location index from the C++ to the Fortran convention.
loc = c_loc + 1

if (loc<1 .or. loc>size(geoval%vals,2)) then
  write(err_msg,*)'ufo_geovals_put_loc_f90',trim(varname),'location out of range:',loc,size(geoval%vals,2)
  call abor1_ftn(err_msg)
endif
if (nlevs /= size(geoval%vals,1)) then
  write(err_msg,*)'ufo_geovals_put_loc_f90',trim(varname),'incorrect number of levels:',nlevs,size(geoval%vals,1)
  call abor1_ftn(err_msg)
endif

geoval%vals(:,loc) = values(:)

end subroutine ufo_geovals_put_loc_c

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

subroutine ufo_geovals_read_file_c(c_key_self, c_conf, c_obspace, c_vars) bind(c,name='ufo_geovals_read_file_f90')
use oops_variables_mod
use datetime_mod

implicit none
integer(c_int), intent(inout)  :: c_key_self
type(c_ptr), value, intent(in) :: c_conf
type(c_ptr), value, intent(in) :: c_obspace
type(c_ptr), value, intent(in) :: c_vars

type(ufo_geovals), pointer :: self
character(max_string)      :: filename
integer :: loc_multiplier
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

vars = oops_variables(c_vars)
! read geovals
call ufo_geovals_read_netcdf(self, filename, loc_multiplier, c_obspace, vars)

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
