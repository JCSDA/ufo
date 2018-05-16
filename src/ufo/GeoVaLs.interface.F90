!
! (C) Copyright 2017 UCAR
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
!
module ufo_geovals_mod_c

use iso_c_binding
use ufo_geovals_mod
use ufo_locs_mod
use ufo_locs_mod_c, only : ufo_locs_registry
use ufo_vars_mod
use kinds

implicit none

public :: ufo_geovals_registry

private
integer, parameter :: max_string=800

#define LISTED_TYPE ufo_geovals

!> Linked list interface - defines registry_t type
#include "linkedList_i.f"

!> Global registry
type(registry_t) :: ufo_geovals_registry

! ------------------------------------------------------------------------------
contains
! ------------------------------------------------------------------------------
!> Linked list implementation
#include "linkedList_c.f"
! ------------------------------------------------------------------------------

subroutine ufo_geovals_setup_c(c_key_self, c_key_locs, c_vars) bind(c,name='ufo_geovals_setup_f90')
use config_mod
implicit none
integer(c_int), intent(in) :: c_key_self
integer(c_int), intent(in) :: c_key_locs
type(c_ptr), intent(in)    :: c_vars

type(ufo_geovals), pointer :: self
type(ufo_locs), pointer :: locs
type(ufo_vars) :: vars

call ufo_geovals_registry%init()
call ufo_geovals_registry%add(c_key_self)
call ufo_geovals_registry%get(c_key_self, self)

call ufo_locs_registry%get(c_key_locs,locs)

call ufo_vars_setup(vars, c_vars)

call ufo_geovals_init(self)
call ufo_geovals_setup(self, vars, locs%nlocs)

end subroutine ufo_geovals_setup_c

! ------------------------------------------------------------------------------

subroutine ufo_geovals_create_c(c_key_self) bind(c,name='ufo_geovals_create_f90')

implicit none
integer(c_int), intent(inout) :: c_key_self

type(ufo_geovals), pointer :: self

call ufo_geovals_registry%init()
call ufo_geovals_registry%add(c_key_self)
call ufo_geovals_registry%get(c_key_self, self)

call ufo_geovals_init(self)

end subroutine ufo_geovals_create_c

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

subroutine ufo_geovals_abs_c(c_key_self) bind(c,name='ufo_geovals_abs_f90')
implicit none
integer(c_int), intent(in) :: c_key_self
type(ufo_geovals), pointer :: self

call ufo_geovals_registry%get(c_key_self, self)

call ufo_geovals_abs(self)

end subroutine ufo_geovals_abs_c

! ------------------------------------------------------------------------------

subroutine ufo_geovals_setup_random_c(c_key_self, c_conf, c_vars) bind(c,name='ufo_geovals_setup_random_f90')
use config_mod
implicit none
integer(c_int), intent(in) :: c_key_self
type(c_ptr), intent(in)    :: c_conf
type(c_ptr), intent(in)    :: c_vars

type(ufo_geovals), pointer :: self
type(ufo_vars) :: vars
integer :: nobs

call ufo_geovals_registry%init()
call ufo_geovals_registry%add(c_key_self)
call ufo_geovals_registry%get(c_key_self, self)

!> read variables
call ufo_vars_setup(vars, c_vars)

! randomize
nobs = config_get_int(c_conf, "nobs")
call ufo_geovals_init(self)
call ufo_geovals_setup(self, vars, nobs)
call ufo_geovals_random(self)

end subroutine ufo_geovals_setup_random_c

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

subroutine ufo_geovals_dotprod_c(c_key_self, c_key_other, prod) bind(c,name='ufo_geovals_dotprod_f90')
implicit none
integer(c_int), intent(in) :: c_key_self, c_key_other
real(c_double), intent(inout) :: prod
type(ufo_geovals), pointer :: self, other

call ufo_geovals_registry%get(c_key_self, self)
call ufo_geovals_registry%get(c_key_other, other)

call ufo_geovals_dotprod(self, other, prod)

end subroutine ufo_geovals_dotprod_c

! ------------------------------------------------------------------------------

subroutine ufo_geovals_minmaxavg_c(c_key_self, kobs, pmin, pmax, prms) bind(c,name='ufo_geovals_minmaxavg_f90')
implicit none
integer(c_int), intent(in) :: c_key_self
integer(c_int), intent(inout) :: kobs
real(c_double), intent(inout) :: pmin, pmax, prms
type(ufo_geovals), pointer :: self

call ufo_geovals_registry%get(c_key_self, self)

call ufo_geovals_minmaxavg(self, kobs, pmin, pmax, prms)

end subroutine ufo_geovals_minmaxavg_c

! ------------------------------------------------------------------------------

subroutine ufo_geovals_read_file_c(c_key_self, c_conf, c_vars) bind(c,name='ufo_geovals_read_file_f90')
use config_mod

implicit none
integer(c_int), intent(in) :: c_key_self
type(c_ptr), intent(in)    :: c_conf
type(c_ptr), intent(in)    :: c_vars

type(ufo_geovals), pointer :: self
type(ufo_vars) :: vars
character(max_string) :: filename

call ufo_geovals_registry%init()
call ufo_geovals_registry%add(c_key_self)
call ufo_geovals_registry%get(c_key_self, self)

!> read variables
call ufo_vars_setup(vars, c_vars)

! read filename for config
filename = config_get_string(c_conf,len(filename),"filename")

! read geovals
call ufo_geovals_read_netcdf(self, filename, vars)

end subroutine ufo_geovals_read_file_c

! ------------------------------------------------------------------------------

subroutine ufo_geovals_write_file_c(c_key_self, c_conf) bind(c,name='ufo_geovals_write_file_f90')
use config_mod
implicit none
integer(c_int), intent(in) :: c_key_self
type(c_ptr), intent(in) :: c_conf
type(ufo_geovals), pointer :: self

call ufo_geovals_registry%get(c_key_self, self)

end subroutine ufo_geovals_write_file_c

! ------------------------------------------------------------------------------

end module ufo_geovals_mod_c
