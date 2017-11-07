! (C) Copyright 2009-2016 ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
! In applying this licence, ECMWF does not waive the privileges and immunities 
! granted to it by virtue of its status as an intergovernmental organisation nor
! does it submit to any jurisdiction.

!> Fortran module for streamfunction observations for the QG model
module qg_obsoper_mod

use iso_c_binding
use config_mod
use qg_vars_mod
use kinds

implicit none
private
public :: qg_obsoper, qg_oper_setup, qg_oper_delete
public :: qg_obsoper_registry

! ------------------------------------------------------------------------------

!> Fortran derived type for stream function observations for the QG model
type :: qg_obsoper
  character(len=30) :: request
  type(qg_vars) :: varin
  integer :: ncol
end type qg_obsoper

#define LISTED_TYPE qg_obsoper

!> Linked list interface - defines registry_t type
#include "linkedList_i.f"

!> Global registry
type(registry_t) :: qg_obsoper_registry

! ------------------------------------------------------------------------------
contains
! ------------------------------------------------------------------------------
!> Linked list implementation
#include "linkedList_c.f"

! ------------------------------------------------------------------------------

subroutine qg_oper_setup(self, c_conf, svars, ncol)
implicit none
type(qg_obsoper), intent(inout) :: self
type(c_ptr), intent(in)    :: c_conf
character(len=*), intent(in) :: svars(:)
integer :: ncol

self%request = config_get_string(c_conf, len(self%request), "ObsType")
call qg_vars_setup(self%varin, svars)
self%ncol = ncol

end subroutine qg_oper_setup

! ------------------------------------------------------------------------------

subroutine qg_oper_delete(self)
implicit none
type(qg_obsoper), intent(inout) :: self

deallocate(self%varin%fldnames)

end subroutine qg_oper_delete

! ------------------------------------------------------------------------------

subroutine c_qg_obsoper_inputs(c_key_self, c_key_vars) bind(c,name='qg_obsoper_inputs_f90')
implicit none
integer(c_int), intent(in)    :: c_key_self
integer(c_int), intent(inout) :: c_key_vars

type(qg_obsoper), pointer :: self
type(qg_vars), pointer :: vars

call qg_obsoper_registry%get(c_key_self, self)
call qg_vars_registry%init()
call qg_vars_registry%add(c_key_vars)
call qg_vars_registry%get(c_key_vars, vars)
call qg_vars_clone(self%varin, vars)

end subroutine c_qg_obsoper_inputs

! ------------------------------------------------------------------------------

end module qg_obsoper_mod
