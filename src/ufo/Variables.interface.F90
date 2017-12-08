!
!  (C) Copyright 2017 UCAR
!  
!  This software is licensed under the terms of the Apache Licence Version 2.0
!  which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
!

module ufo_vars_mod_c

use iso_c_binding
use config_mod
use ufo_vars_mod

implicit none

public :: ufo_vars_registry

private

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

subroutine ufo_vars_create_c(c_key_self, c_conf) bind(c,name='ufo_var_create_f90')
implicit none
integer(c_int), intent(inout) :: c_key_self
type(c_ptr), intent(in)    :: c_conf

type(ufo_vars), pointer :: self
character(len=512) :: svars
integer :: nvar
character(len=MAXVARLEN), allocatable :: cvars(:)

call ufo_vars_registry%init()
call ufo_vars_registry%add(c_key_self)
call ufo_vars_registry%get(c_key_self, self)

nvar = config_get_int(c_conf, "nvars")
allocate(cvars(nvar))

svars = config_get_string(c_conf,len(svars),"variables")
read(svars,*) cvars

call ufo_vars_setup(self, cvars)

deallocate(cvars)

end subroutine ufo_vars_create_c

! ------------------------------------------------------------------------------

subroutine ufo_vars_clone_c(c_key_self, c_key_other) bind(c,name='ufo_var_clone_f90')
implicit none
integer(c_int), intent(in)    :: c_key_self
integer(c_int), intent(inout) :: c_key_other

type(ufo_vars), pointer :: self, other

call ufo_vars_registry%get(c_key_self, self)
call ufo_vars_registry%get(c_key_other, other)

call ufo_vars_clone(self, other)

end subroutine ufo_vars_clone_c

! ------------------------------------------------------------------------------

subroutine ufo_vars_delete_c(c_key_self) bind(c,name='ufo_var_delete_f90')

implicit none
integer(c_int), intent(inout) :: c_key_self

type(ufo_vars), pointer :: self

call ufo_vars_registry%get(c_key_self, self)

call ufo_vars_delete(self)

end subroutine ufo_vars_delete_c

! ------------------------------------------------------------------------------

!subroutine ufo_vars_info_c(c_key_self, c_nv, lline, c_line) bind(c,name='ufo_var_info_f90')
subroutine ufo_vars_info_c(c_key_self, c_nv) bind(c,name='ufo_var_info_f90')
implicit none
integer(c_int), intent(in)    :: c_key_self
integer(c_int), intent(inout) :: c_nv
!integer(c_int), intent(in)    :: lline
!character(kind=c_char,len=1), intent(inout) :: c_line(lline+1)

type(ufo_vars), pointer :: self

call ufo_vars_registry%get(c_key_self, self)

c_nv = ufo_vars_nvars(self)

end subroutine ufo_vars_info_c

end module ufo_vars_mod_c
