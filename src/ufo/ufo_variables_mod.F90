
!> Fortran module to handle variables for UFO
module ufo_vars_mod

use iso_c_binding
use config_mod

implicit none
private
public :: ufo_vars, ufo_vars_setup, ufo_vars_clone
public :: ufo_vars_registry

! ------------------------------------------------------------------------------

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
integer :: jj

self%nv = size(cvars)
allocate(self%fldnames(self%nv))
self%fldnames(:)=cvars(:)

end subroutine ufo_vars_setup

! ------------------------------------------------------------------------------

subroutine c_ufo_vars_create(c_key_self, c_conf) bind(c,name='ufo_var_create_f90')
implicit none
integer(c_int), intent(inout) :: c_key_self
type(c_ptr), intent(in)    :: c_conf

type(ufo_vars), pointer :: self
character(len=2) :: svar

call ufo_vars_registry%init()
call ufo_vars_registry%add(c_key_self)
call ufo_vars_registry%get(c_key_self, self)

self%nv = 0

return
end subroutine c_ufo_vars_create

! ------------------------------------------------------------------------------

subroutine c_ufo_vars_clone(c_key_self, c_key_other) bind(c,name='ufo_var_clone_f90')
implicit none
integer(c_int), intent(in)    :: c_key_self
integer(c_int), intent(inout) :: c_key_other

type(ufo_vars), pointer :: self, other

call ufo_vars_registry%get(c_key_self, self)
call ufo_vars_registry%add(c_key_other)
call ufo_vars_registry%get(c_key_other, other)

call ufo_vars_clone(self, other)

end subroutine c_ufo_vars_clone

! ------------------------------------------------------------------------------

subroutine ufo_vars_clone(self, other)
implicit none
type(ufo_vars), intent(in)    :: self
type(ufo_vars), intent(inout) :: other

other%nv = self%nv

end subroutine ufo_vars_clone

! ------------------------------------------------------------------------------

subroutine c_ufo_vars_delete(c_key_self) bind(c,name='ufo_var_delete_f90')

implicit none
integer(c_int), intent(inout) :: c_key_self
type(ufo_vars), pointer :: self

call ufo_vars_registry%get(c_key_self, self)
if (allocated(self%fldnames)) deallocate(self%fldnames)
call ufo_vars_registry%remove(c_key_self)

return
end subroutine c_ufo_vars_delete

! ------------------------------------------------------------------------------

subroutine c_ufo_vars_info(c_key_self, c_nv, c_nl) bind(c,name='ufo_var_info_f90')
implicit none
integer(c_int), intent(in)    :: c_key_self
integer(c_int), intent(inout) :: c_nv
integer(c_int), intent(inout) :: c_nl
type(ufo_vars), pointer :: self

call ufo_vars_registry%get(c_key_self, self)

c_nv = self%nv

return
end subroutine c_ufo_vars_info

! ------------------------------------------------------------------------------

end module ufo_vars_mod
