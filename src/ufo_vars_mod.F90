! (C) Copyright 2009-2016 ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
! In applying this licence, ECMWF does not waive the privileges and immunities 
! granted to it by virtue of its status as an intergovernmental organisation nor
! does it submit to any jurisdiction.

!> Fortran module to handle variables for the QG model
module qg_vars_mod

use iso_c_binding
use config_mod

implicit none
private
public :: qg_vars, qg_vars_setup, qg_vars_clone
public :: qg_vars_registry

! ------------------------------------------------------------------------------

!> Fortran derived type to represent QG model variables
type :: qg_vars
  integer :: nv
  character(len=1), allocatable :: fldnames(:) !< Variable identifiers
  logical :: lbc
end type qg_vars

#define LISTED_TYPE qg_vars

!> Linked list interface - defines registry_t type
#include "linkedList_i.f"

!> Global registry
type(registry_t) :: qg_vars_registry

! ------------------------------------------------------------------------------
contains
! ------------------------------------------------------------------------------
!> Linked list implementation
#include "linkedList_c.f"

! ------------------------------------------------------------------------------

subroutine qg_vars_setup(self, cvars)
implicit none
type(qg_vars), intent(inout) :: self
character(len=1), intent(in) :: cvars(:)
integer :: jj

self%nv = size(cvars)
self%lbc = .false.

do jj=1,self%nv
  if (cvars(jj)/="x" .and. cvars(jj)/="q" .and. cvars(jj)/="u" .and. cvars(jj)/="v") &
     & call abor1_ftn ("qg_vars_setup: unknown field")
enddo
allocate(self%fldnames(self%nv))
self%fldnames(:)=cvars(:)

end subroutine qg_vars_setup

! ------------------------------------------------------------------------------

subroutine c_qg_vars_create(c_key_self, c_conf) bind(c,name='qg_var_create_f90')
implicit none
integer(c_int), intent(inout) :: c_key_self
type(c_ptr), intent(in)    :: c_conf

type(qg_vars), pointer :: self
character(len=2) :: svar

call qg_vars_registry%init()
call qg_vars_registry%add(c_key_self)
call qg_vars_registry%get(c_key_self, self)

svar = config_get_string(c_conf,len(svar),"variables")
select case (svar)
case ("nl")
  self%nv = 4
  self%lbc = .true.
  allocate(self%fldnames(4))
  self%fldnames(:) = (/"x","q","u","v"/)
case ("tl")
  self%nv = 4
  self%lbc = .false.
  allocate(self%fldnames(4))
  self%fldnames(:) = (/"x","q","u","v"/)
case ("cv")
  self%nv = 1
  self%lbc = .true.
  allocate(self%fldnames(1))
  self%fldnames(1) = "x"
case ("ci")
  self%nv = 1
  self%lbc = .false.
  allocate(self%fldnames(1))
  self%fldnames(1) = "x"
case ("x")
  self%nv = 1
  self%lbc = .false.
  allocate(self%fldnames(1))
  self%fldnames(1) = "x"
case default
  call abor1_ftn("c_qg_vars_create: undefined variables")
end select

return
end subroutine c_qg_vars_create

! ------------------------------------------------------------------------------

subroutine c_qg_vars_clone(c_key_self, c_key_other) bind(c,name='qg_var_clone_f90')
implicit none
integer(c_int), intent(in)    :: c_key_self
integer(c_int), intent(inout) :: c_key_other

type(qg_vars), pointer :: self, other

call qg_vars_registry%get(c_key_self, self)
call qg_vars_registry%add(c_key_other)
call qg_vars_registry%get(c_key_other, other)

call qg_vars_clone(self, other)

end subroutine c_qg_vars_clone

! ------------------------------------------------------------------------------

subroutine qg_vars_clone(self, other)
implicit none
type(qg_vars), intent(in)    :: self
type(qg_vars), intent(inout) :: other

other%nv = self%nv
other%lbc = self%lbc

allocate(other%fldnames(other%nv))
other%fldnames(:)=self%fldnames(:)

end subroutine qg_vars_clone

! ------------------------------------------------------------------------------

subroutine c_qg_vars_delete(c_key_self) bind(c,name='qg_var_delete_f90')

implicit none
integer(c_int), intent(inout) :: c_key_self

type(qg_vars), pointer :: self
call qg_vars_registry%get(c_key_self, self)
deallocate(self%fldnames)
call qg_vars_registry%remove(c_key_self)

return
end subroutine c_qg_vars_delete

! ------------------------------------------------------------------------------

subroutine c_qg_vars_info(c_key_self, c_nv, c_nl) bind(c,name='qg_var_info_f90')
implicit none
integer(c_int), intent(in)    :: c_key_self
integer(c_int), intent(inout) :: c_nv
integer(c_int), intent(inout) :: c_nl
type(qg_vars), pointer :: self

call qg_vars_registry%get(c_key_self, self)

c_nv = self%nv
c_nl = 0
if (self%lbc) c_nl = 1

!self%fldnames(:)

return
end subroutine c_qg_vars_info

! ------------------------------------------------------------------------------

end module qg_vars_mod
