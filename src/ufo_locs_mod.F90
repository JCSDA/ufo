! (C) Copyright 2009-2016 ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
! In applying this licence, ECMWF does not waive the privileges and immunities 
! granted to it by virtue of its status as an intergovernmental organisation nor
! does it submit to any jurisdiction.

!> Fortran module handling observation locations

module qg_locs_mod

use iso_c_binding
use qg_obs_vectors
use kinds

implicit none
private
public :: qg_locs, qg_loc_setup
public :: qg_locs_registry

! ------------------------------------------------------------------------------

!> Fortran derived type to hold observation locations
type :: qg_locs
  integer :: nloc
  real(kind=kind_real), allocatable :: xyz(:,:)
end type qg_locs

#define LISTED_TYPE qg_locs

!> Linked list interface - defines registry_t type
#include "linkedList_i.f"

!> Global registry
type(registry_t) :: qg_locs_registry

! ------------------------------------------------------------------------------
contains
! ------------------------------------------------------------------------------
!> Linked list implementation
#include "linkedList_c.f"

! ------------------------------------------------------------------------------

subroutine qg_loc_setup(self, lvec)
implicit none
type(qg_locs), intent(inout) :: self
type(obs_vect), intent(in) :: lvec
integer :: jc, jo

self%nloc=lvec%nobs
allocate(self%xyz(3,self%nloc))
do jo=1,self%nloc
  do jc=1,3
    self%xyz(jc,jo)=lvec%values(jc,jo)
  enddo
enddo

end subroutine qg_loc_setup

! ------------------------------------------------------------------------------

subroutine c_qg_loc_delete(key) bind(c,name='qg_loc_delete_f90')

implicit none
integer(c_int), intent(inout) :: key
type(qg_locs), pointer :: self

call qg_locs_registry%get(key,self)
deallocate(self%xyz)
call qg_locs_registry%remove(key)

end subroutine c_qg_loc_delete

! ------------------------------------------------------------------------------

subroutine c_qg_loc_nobs(key, kobs) bind(c,name='qg_loc_nobs_f90')

implicit none
integer(c_int), intent(in) :: key
integer(c_int), intent(inout) :: kobs
type(qg_locs), pointer :: self

call qg_locs_registry%get(key,self)
kobs = self%nloc

end subroutine c_qg_loc_nobs

! ------------------------------------------------------------------------------

end module qg_locs_mod
