!
! (C) Copyright 2017 UCAR
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
!
module ufo_obscheck_mod

use iso_c_binding
use ufo_vars_mod
use ufo_obs_vectors
use ufo_locs_mod
use ufo_geovals_mod
use kinds

implicit none
private
public :: ufo_obscheck_registry

! ------------------------------------------------------------------------------

!> Fortran derived type to hold interpolated fields required by the obs operators
type :: ufo_obscheck
  integer :: nobs
  integer :: nvar
  integer :: used
  integer, allocatable :: indx(:)
  real(kind=kind_real), allocatable :: values(:,:)
  character(len=1), allocatable :: variables(:)
  logical :: lalloc
end type ufo_obscheck

#define LISTED_TYPE ufo_obscheck

!> Linked list interface - defines registry_t type
#include "linkedList_i.f"

!> Global registry
type(registry_t) :: ufo_obscheck_registry

! ------------------------------------------------------------------------------
contains
! ------------------------------------------------------------------------------
!> Linked list implementation
#include "linkedList_c.f"

! ------------------------------------------------------------------------------

subroutine c_ufo_obscheck_create(c_key_self) bind(c,name='ufo_obscheck_create_f90')

implicit none
integer(c_int), intent(inout) :: c_key_self

type(ufo_obscheck), pointer :: self

call ufo_obscheck_registry%init()
call ufo_obscheck_registry%add(c_key_self)
call ufo_obscheck_registry%get(c_key_self, self)

write(*,*) 'Here in c_ufo_obscheck_create ========='
self%lalloc = .false.

end subroutine c_ufo_obscheck_create

! ------------------------------------------------------------------------------

subroutine ufo_obscheck_read_file_c(c_key_self, c_conf) bind(c,name='ufo_obscheck_read_file_f90')
use config_mod
use fckit_log_module, only : fckit_log
implicit none
integer(c_int), intent(in) :: c_key_self
type(c_ptr), intent(in)    :: c_conf
type(ufo_obscheck), pointer :: self

call ufo_obscheck_registry%get(c_key_self, self)

write(*,*) 'Here in ufo_obscheck_read_file_c========='
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

end subroutine ufo_obscheck_read_file_c

! ------------------------------------------------------------------------------

subroutine c_ufo_postFilter_f90(c_key_geovals, c_key_hofx) bind(c,name='ufo_postFilter_f90')

    implicit none
    integer(c_int), intent(in) :: c_key_geovals
    integer(c_int), intent(in) :: c_key_hofx
    type(ufo_geovals), pointer  :: geovals
    type(obs_vector), pointer :: hofx


write(*,*) 'ufo_postFilter_f90_c========='

end subroutine c_ufo_postFilter_f90

! ------------------------------------------------------------------------------

subroutine c_ufo_priorFilter_f90(c_key_geovals, c_key_hofx) bind(c,name='ufo_priorFilter_f90')

    implicit none
    integer(c_int), intent(in) :: c_key_geovals
    integer(c_int), intent(in) :: c_key_hofx
    type(ufo_geovals), pointer  :: geovals
    type(obs_vector), pointer :: hofx


write(*,*) 'ufo_priorFilter_f90_c========='

end subroutine c_ufo_priorFilter_f90

! ------------------------------------------------------------------------------

end module ufo_obscheck_mod
