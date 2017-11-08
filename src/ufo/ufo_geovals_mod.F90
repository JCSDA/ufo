!
! (C) Copyright 2017 UCAR
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
!
module ufo_geovals_mod

use iso_c_binding
use ufo_vars_mod
use kinds

implicit none
private
public :: ufo_geovals, geovals_setup
public :: ufo_geovals_registry

! ------------------------------------------------------------------------------

!> Fortran derived type to hold interpolated fields required by the obs operators
type :: ufo_geovals
  integer :: nobs
  integer :: nvar
  integer :: used
  integer, allocatable :: indx(:)
  real(kind=kind_real), allocatable :: values(:,:)
  character(len=1), allocatable :: variables(:)
  logical :: lalloc
end type ufo_geovals

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

subroutine c_ufo_geovals_create(c_key_self) bind(c,name='ufo_geovals_create_f90')

implicit none
integer(c_int), intent(inout) :: c_key_self

type(ufo_geovals), pointer :: self

call ufo_geovals_registry%init()
call ufo_geovals_registry%add(c_key_self)
call ufo_geovals_registry%get(c_key_self, self)

self%lalloc = .false.

end subroutine c_ufo_geovals_create

! ------------------------------------------------------------------------------

subroutine geovals_setup(self, vars, kobs)
implicit none
type(ufo_geovals), intent(inout) :: self
type(ufo_vars), intent(in) :: vars
integer, intent(in) :: kobs(:)

self%nobs=size(kobs)

!self%nvar=vars%nv
self%nvar=1
self%used=0

allocate(self%indx(self%nobs))
!self%indx(:)=kobs(:)
self%indx(1)=1

allocate(self%variables(self%nvar))
!self%variables(:)=vars%fldnames(:)
self%variables(1)="z"

allocate(self%values(self%nvar,self%nobs))

self%lalloc = .true.

end subroutine geovals_setup

! ------------------------------------------------------------------------------

subroutine c_ufo_geovals_delete(c_key_self) bind(c,name='ufo_geovals_delete_f90')

implicit none
integer(c_int), intent(inout) :: c_key_self

type(ufo_geovals), pointer :: self

call ufo_geovals_registry%get(c_key_self, self)
if (self%lalloc) then
  deallocate(self%values)
  deallocate(self%indx)
  deallocate(self%variables)
endif
call ufo_geovals_registry%remove(c_key_self)

end subroutine c_ufo_geovals_delete

! ------------------------------------------------------------------------------

subroutine c_ufo_geovals_zero(c_key_self) bind(c,name='ufo_geovals_zero_f90')
implicit none
integer(c_int), intent(in) :: c_key_self
type(ufo_geovals), pointer :: self
call ufo_geovals_registry%get(c_key_self, self)
self%values(:,:)=0.0_kind_real
end subroutine c_ufo_geovals_zero

! ------------------------------------------------------------------------------

subroutine c_ufo_geovals_random(c_key_self) bind(c,name='ufo_geovals_random_f90')
use random_vectors_mod
implicit none
integer(c_int), intent(in) :: c_key_self
type(ufo_geovals), pointer :: self
call ufo_geovals_registry%get(c_key_self, self)
!call random_vector(self%values(:,:))
self%values(:,:)=1.0_kind_real
end subroutine c_ufo_geovals_random

! ------------------------------------------------------------------------------

subroutine c_ufo_geovals_dotprod(c_key_self, c_key_other, prod) bind(c,name='ufo_geovals_dotprod_f90')
implicit none
integer(c_int), intent(in) :: c_key_self, c_key_other
real(c_double), intent(inout) :: prod
type(ufo_geovals), pointer :: self, other
integer :: jo, jv

call ufo_geovals_registry%get(c_key_self, self)
call ufo_geovals_registry%get(c_key_other, other)
prod=0.0_kind_real
do jo=1,self%nobs
  do jv=1,self%nvar
    prod=prod+self%values(jv,jo)*other%values(jv,jo)
  enddo
enddo

end subroutine c_ufo_geovals_dotprod

! ------------------------------------------------------------------------------

subroutine c_ufo_geovals_minmaxavg(c_key_self, kobs, pmin, pmax, prms) bind(c,name='ufo_geovals_minmaxavg_f90')
implicit none
integer(c_int), intent(in) :: c_key_self
integer(c_int), intent(inout) :: kobs
real(c_double), intent(inout) :: pmin, pmax, prms
type(ufo_geovals), pointer :: self

call ufo_geovals_registry%get(c_key_self, self)

kobs = self%nobs
pmin=minval(self%values(:,:))
pmax=maxval(self%values(:,:))
prms=sqrt(sum(self%values(:,:)**2)/real(self%nobs*self%nvar,kind_real))

end subroutine c_ufo_geovals_minmaxavg

! ------------------------------------------------------------------------------

subroutine ufo_geovals_read_file_c(c_key_self, c_conf) bind(c,name='ufo_geovals_read_file_f90')
use config_mod
use fckit_log_module, only : fckit_log
implicit none
integer(c_int), intent(in) :: c_key_self
type(c_ptr), intent(in)    :: c_conf
type(ufo_geovals), pointer :: self

call ufo_geovals_registry%get(c_key_self, self)

self%nobs=1
self%nvar=1
self%used=0

allocate(self%indx(self%nobs))
self%indx(1)=1

allocate(self%variables(self%nvar))
self%variables(1)="z"

allocate(self%values(self%nvar,self%nobs))

self%values(:,:)=1.0_kind_real

self%lalloc = .true.

end subroutine ufo_geovals_read_file_c

! ------------------------------------------------------------------------------

subroutine ufo_geovals_write_file_c(c_key_self, c_conf) bind(c,name='ufo_geovals_write_file_f90')
use config_mod
use fckit_log_module, only : fckit_log
implicit none
integer(c_int), intent(in) :: c_key_self
type(c_ptr), intent(in) :: c_conf
type(ufo_geovals), pointer :: self

call ufo_geovals_registry%get(c_key_self, self)

end subroutine ufo_geovals_write_file_c

! ------------------------------------------------------------------------------

end module ufo_geovals_mod
