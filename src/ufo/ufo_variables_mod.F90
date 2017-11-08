!
!  (C) Copyright 2017 UCAR
!  
!  This software is licensed under the terms of the Apache Licence Version 2.0
!  which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
!

module ufo_vars_mod

use iso_c_binding
use config_mod

implicit none
private
public :: ufo_vars, ufo_vars_setup, ufo_vars_clone, ufo_vars_delete
public :: ufo_vars_registry

! ------------------------------------------------------------------------------

!> Fortran derived type to represent QG model variables
type :: ufo_vars
  integer :: nv
  character(len=1), allocatable :: fldnames(:) !< Variable identifiers
end type ufo_vars

#define LISTED_TYPE ufo_vars

!> Linked list interface - defines registry_t type
#include "linkedList_i.f"

!> Global registry
type(registry_t) :: ufo_vars_registry

! ------------------------------------------------------------------------------
contains
! ------------------------------------------------------------------------------
!> Linked list implementation
#include "linkedList_c.f"

! ------------------------------------------------------------------------------

subroutine ufo_vars_setup(self, cvars)
implicit none
type(ufo_vars), intent(inout) :: self
character(len=1), intent(in) :: cvars(:)

end subroutine ufo_vars_setup

! ------------------------------------------------------------------------------

subroutine ufo_vars_clone(self, other)
implicit none
type(ufo_vars), intent(inout) :: self
type(ufo_vars), intent(in) :: other
end subroutine ufo_vars_clone

! ------------------------------------------------------------------------------

subroutine ufo_vars_delete(self)
implicit none
type(ufo_vars), intent(inout) :: self
end subroutine ufo_vars_delete

! ------------------------------------------------------------------------------

subroutine c_ufo_vars_create(c_key_self, c_conf) bind(c,name='ufo_var_create_f90')
implicit none
integer(c_int), intent(inout) :: c_key_self
type(c_ptr), intent(in)    :: c_conf

c_key_self=1
return
end subroutine c_ufo_vars_create

! ------------------------------------------------------------------------------

subroutine c_ufo_vars_clone(c_key_self, c_key_other) bind(c,name='ufo_var_clone_f90')
implicit none
integer(c_int), intent(in)    :: c_key_self
integer(c_int), intent(inout) :: c_key_other

c_key_other=c_key_self
end subroutine c_ufo_vars_clone

! ------------------------------------------------------------------------------

subroutine c_ufo_vars_delete(c_key_self) bind(c,name='ufo_var_delete_f90')

implicit none
integer(c_int), intent(inout) :: c_key_self

c_key_self=0
return
end subroutine c_ufo_vars_delete

! ------------------------------------------------------------------------------

end module ufo_vars_mod
