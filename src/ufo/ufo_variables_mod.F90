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
public :: ufo_vars_registry, ufo_vars_readconfig
public :: ufo_vars_getindex
public :: MAXVARLEN

integer, parameter :: MAXVARLEN=8
! ------------------------------------------------------------------------------

!> Fortran derived type to represent QG model variables
type :: ufo_vars
  integer :: nv
  character(len=MAXVARLEN), allocatable :: fldnames(:) !< Variable identifiers
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
character(len=MAXVARLEN), intent(in) :: cvars(:)
self%nv = size(cvars)

! TODO: a check on whether this var is in the list of defined vars
allocate(self%fldnames(self%nv))
self%fldnames(:) = cvars(:)

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

subroutine ufo_vars_readconfig(self, c_conf) 
implicit none
type(ufo_vars), intent(inout) :: self
type(c_ptr), intent(in)    :: c_conf

character(len=MAXVARLEN) :: svar

svar = config_get_string(c_conf,len(svar),"variables")
self%nv = 1
allocate(self%fldnames(self%nv))
self%fldnames(1) = svar
return
end subroutine ufo_vars_readconfig

! ------------------------------------------------------------------------------
integer function ufo_vars_getindex(self, varname)
implicit none
type(ufo_vars), intent(in)       :: self
character(MAXVARLEN), intent(in) :: varname

integer :: ivar

do ivar = 1, self%nv
  if (self%fldnames(ivar) == varname) then
    exit
  endif
enddo

if (ivar <= self%nv) then
  ufo_vars_getindex = ivar
else
  ufo_vars_getindex = -1
endif

end function ufo_vars_getindex

! ------------------------------------------------------------------------------

subroutine ufo_vars_create_c(c_key_self, c_conf) bind(c,name='ufo_var_create_f90')
implicit none
integer(c_int), intent(inout) :: c_key_self
type(c_ptr), intent(in)    :: c_conf

type(ufo_vars), pointer :: self
character(len=MAXVARLEN) :: svar

call ufo_vars_registry%init()
call ufo_vars_registry%add(c_key_self)
call ufo_vars_registry%get(c_key_self, self)

call ufo_vars_readconfig(self, c_conf)

return
end subroutine ufo_vars_create_c

! ------------------------------------------------------------------------------

subroutine ufo_vars_clone_c(c_key_self, c_key_other) bind(c,name='ufo_var_clone_f90')
implicit none
integer(c_int), intent(in)    :: c_key_self
integer(c_int), intent(inout) :: c_key_other

c_key_other=c_key_self
end subroutine ufo_vars_clone_c

! ------------------------------------------------------------------------------

subroutine ufo_vars_delete_c(c_key_self) bind(c,name='ufo_var_delete_f90')

implicit none
integer(c_int), intent(inout) :: c_key_self

c_key_self=0
return
end subroutine ufo_vars_delete_c

! ------------------------------------------------------------------------------

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

!c_line = self%fldnames(1)
c_nv = self%nv

return
end subroutine ufo_vars_info_c

end module ufo_vars_mod
