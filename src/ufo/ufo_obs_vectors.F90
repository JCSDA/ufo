! (C) Copyright 2009-2016 ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
! In applying this licence, ECMWF does not waive the privileges and immunities 
! granted to it by virtue of its status as an intergovernmental organisation nor
! does it submit to any jurisdiction.

!> Fortran module handling observation vectors

module ufo_obs_vectors

use iso_c_binding
use random_vectors_mod
use kinds

implicit none
private
public :: obs_vector, obsvec_setup
public :: ufo_obs_vect_registry

! ------------------------------------------------------------------------------

!> Fortran derived type to represent an observation vector
type obs_vector
  integer :: nobs=0
  real(kind=kind_real), allocatable :: values(:)
end type obs_vector

#define LISTED_TYPE obs_vector

!> Linked list interface - defines registry_t type
#include "linkedList_i.f"

!> Global registry
type(registry_t) :: ufo_obs_vect_registry

! ------------------------------------------------------------------------------
contains
! ------------------------------------------------------------------------------
!> Linked list implementation
#include "linkedList_c.f"

! ------------------------------------------------------------------------------

subroutine ufo_obsvec_setup_c(c_key_self, c_nobs) bind(c,name='ufo_obsvec_setup_f90')
implicit none
integer(c_int), intent(inout) :: c_key_self
integer(c_int), intent(in) :: c_nobs
type(obs_vector), pointer :: self
integer :: iobs

call ufo_obs_vect_registry%init()
call ufo_obs_vect_registry%add(c_key_self)
call ufo_obs_vect_registry%get(c_key_self,self)
iobs = c_nobs
call obsvec_setup(self, iobs)

end subroutine ufo_obsvec_setup_c

! ------------------------------------------------------------------------------

subroutine obsvec_setup(self, kobs)
implicit none
type(obs_vector), intent(inout) :: self
integer, intent(in) :: kobs

self%nobs=kobs
if (allocated(self%values)) deallocate(self%values)
allocate(self%values(self%nobs))

end subroutine obsvec_setup

! ------------------------------------------------------------------------------

subroutine ufo_obsvec_clone_c(c_key_self, c_key_other) bind(c,name='ufo_obsvec_clone_f90')
implicit none
integer(c_int), intent(in)    :: c_key_self
integer(c_int), intent(inout) :: c_key_other
type(obs_vector), pointer :: self, other

call ufo_obs_vect_registry%get(c_key_self,self)
call ufo_obs_vect_registry%init()
call ufo_obs_vect_registry%add(c_key_other)
call ufo_obs_vect_registry%get(c_key_other,other)
other%nobs=self%nobs
allocate(other%values(other%nobs))

end subroutine ufo_obsvec_clone_c

! ------------------------------------------------------------------------------

subroutine ufo_obsvec_delete_c(c_key_self) bind(c,name='ufo_obsvec_delete_f90')
implicit none
integer(c_int), intent(inout) :: c_key_self
type(obs_vector), pointer :: self

call ufo_obs_vect_registry%get(c_key_self,self)
deallocate(self%values)
call ufo_obs_vect_registry%remove(c_key_self)

end subroutine ufo_obsvec_delete_c

! ------------------------------------------------------------------------------
subroutine ufo_obsvec_assign_c(c_key_self, c_key_rhs) bind(c,name='ufo_obsvec_assign_f90')
implicit none
integer(c_int), intent(in) :: c_key_self
integer(c_int), intent(in) :: c_key_rhs
type(obs_vector), pointer :: self, rhs

call ufo_obs_vect_registry%get(c_key_self,self)
call ufo_obs_vect_registry%get(c_key_rhs,rhs)
if (rhs%nobs/=self%nobs) then
  deallocate(self%values)
  self%nobs=rhs%nobs
  allocate(self%values(self%nobs))
endif
self%values(:)=rhs%values(:)

end subroutine ufo_obsvec_assign_c
! ------------------------------------------------------------------------------
subroutine ufo_obsvec_zero_c(c_key_self) bind(c,name='ufo_obsvec_zero_f90')
implicit none
integer(c_int), intent(in) :: c_key_self
type(obs_vector), pointer :: self

call ufo_obs_vect_registry%get(c_key_self,self)
self%values(:)=0.0_kind_real

end subroutine ufo_obsvec_zero_c
! ------------------------------------------------------------------------------
subroutine ufo_obsvec_mul_scal_c(c_key_self, zz) bind(c,name='ufo_obsvec_mul_scal_f90')
implicit none
integer(c_int), intent(in) :: c_key_self
real(c_double), intent(in) :: zz
type(obs_vector), pointer :: self

call ufo_obs_vect_registry%get(c_key_self,self)
self%values(:)=zz*self%values(:)

end subroutine ufo_obsvec_mul_scal_c
! ------------------------------------------------------------------------------
subroutine ufo_obsvec_add_c(c_key_self, c_key_other) bind(c,name='ufo_obsvec_add_f90')
implicit none
integer(c_int), intent(in) :: c_key_self
integer(c_int), intent(in) :: c_key_other
type(obs_vector), pointer :: self, other

call ufo_obs_vect_registry%get(c_key_self,self)
call ufo_obs_vect_registry%get(c_key_other,other)
self%values(:)=self%values(:)+other%values(:)

end subroutine ufo_obsvec_add_c
! ------------------------------------------------------------------------------
subroutine ufo_obsvec_sub_c(c_key_self, c_key_other) bind(c,name='ufo_obsvec_sub_f90')
implicit none
integer(c_int), intent(in) :: c_key_self
integer(c_int), intent(in) :: c_key_other
type(obs_vector), pointer :: self, other

call ufo_obs_vect_registry%get(c_key_self,self)
call ufo_obs_vect_registry%get(c_key_other,other)
self%values(:)=self%values(:)-other%values(:)

end subroutine ufo_obsvec_sub_c
! ------------------------------------------------------------------------------
subroutine ufo_obsvec_mul_c(c_key_self, c_key_other) bind(c,name='ufo_obsvec_mul_f90')
implicit none
integer(c_int), intent(in) :: c_key_self
integer(c_int), intent(in) :: c_key_other
type(obs_vector), pointer :: self, other

call ufo_obs_vect_registry%get(c_key_self,self)
call ufo_obs_vect_registry%get(c_key_other,other)
self%values(:)=self%values(:)*other%values(:)

end subroutine ufo_obsvec_mul_c
! ------------------------------------------------------------------------------
subroutine ufo_obsvec_div_c(c_key_self, c_key_other) bind(c,name='ufo_obsvec_div_f90')
implicit none
integer(c_int), intent(in) :: c_key_self
integer(c_int), intent(in) :: c_key_other
type(obs_vector), pointer :: self, other

call ufo_obs_vect_registry%get(c_key_self,self)
call ufo_obs_vect_registry%get(c_key_other,other)
self%values(:)=self%values(:)/other%values(:)

end subroutine ufo_obsvec_div_c
! ------------------------------------------------------------------------------
subroutine ufo_obsvec_axpy_c(c_key_self, zz, c_key_other) bind(c,name='ufo_obsvec_axpy_f90')
implicit none
integer(c_int), intent(in) :: c_key_self
real(c_double), intent(in) :: zz
integer(c_int), intent(in) :: c_key_other
type(obs_vector), pointer :: self, other

call ufo_obs_vect_registry%get(c_key_self,self)
call ufo_obs_vect_registry%get(c_key_other,other)
self%values(:)=self%values(:)+zz*other%values(:)

end subroutine ufo_obsvec_axpy_c
! ------------------------------------------------------------------------------
subroutine ufo_obsvec_invert_c(c_key_self) bind(c,name='ufo_obsvec_invert_f90')
implicit none
integer(c_int), intent(in) :: c_key_self
type(obs_vector), pointer :: self

call ufo_obs_vect_registry%get(c_key_self,self)
self%values(:)=1.0_kind_real/self%values(:)

end subroutine ufo_obsvec_invert_c
! ------------------------------------------------------------------------------
subroutine ufo_obsvec_random_c(c_key_self) bind(c,name='ufo_obsvec_random_f90')
implicit none
integer(c_int), intent(in) :: c_key_self
type(obs_vector), pointer :: self

call ufo_obs_vect_registry%get(c_key_self,self)
call random_vector(self%values)

end subroutine ufo_obsvec_random_c
! ------------------------------------------------------------------------------
subroutine ufo_obsvec_dotprod_c(c_key_self, c_key_other, zz) bind(c,name='ufo_obsvec_dotprod_f90')
implicit none
integer(c_int), intent(in) :: c_key_self
integer(c_int), intent(in) :: c_key_other
real(c_double), intent(inout) :: zz
type(obs_vector), pointer :: self, other
integer :: jo

zz=0.0_kind_real
call ufo_obs_vect_registry%get(c_key_self,self)
call ufo_obs_vect_registry%get(c_key_other,other)
do jo=1,self%nobs
  zz = zz + self%values(jo)*other%values(jo)
enddo

end subroutine ufo_obsvec_dotprod_c
! ------------------------------------------------------------------------------
subroutine ufo_obsvec_minmaxavg_c(c_key_self, zmin, zmax, zavg) bind(c,name='ufo_obsvec_minmaxavg_f90')
implicit none
integer(c_int), intent(in)    :: c_key_self
real(c_double), intent(inout) :: zmin, zmax, zavg
type(obs_vector), pointer :: self

call ufo_obs_vect_registry%get(c_key_self,self)
if (self%nobs>0) then
  if (.not.allocated(self%values)) call abor1_ftn("obsvec_minmax: obs vector not allocated")
  zmin = minval(self%values)
  zmax = maxval(self%values)
  zavg = sum(self%values)/real(self%nobs)
else
  zmin=0.0_kind_real
  zmax=0.0_kind_real
  zavg=0.0_kind_real
endif
end subroutine ufo_obsvec_minmaxavg_c
! ------------------------------------------------------------------------------
subroutine ufo_obsvec_nobs_c(c_key_self, kobs) bind(c,name='ufo_obsvec_nobs_f90')
implicit none
integer(c_int), intent(in)    :: c_key_self
integer(c_int), intent(inout) :: kobs
type(obs_vector), pointer :: self

call ufo_obs_vect_registry%get(c_key_self,self)
kobs=self%nobs

end subroutine ufo_obsvec_nobs_c
! ------------------------------------------------------------------------------

end module ufo_obs_vectors
