! (C) Copyright 2009-2016 ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
! In applying this licence, ECMWF does not waive the privileges and immunities 
! granted to it by virtue of its status as an intergovernmental organisation nor
! does it submit to any jurisdiction.

!> Fortran module handling observation locations

module ufo_locs_mod

use iso_c_binding
use ufo_obs_vectors
use kinds

implicit none
private
public :: ufo_locs, ufo_loc_setup
public :: ufo_locs_registry

! ------------------------------------------------------------------------------

!> Fortran derived type to hold observation locations
type :: ufo_locs
  integer :: nloc
  real(kind=kind_real), allocatable :: xyz(:,:)
end type ufo_locs

#define LISTED_TYPE ufo_locs

!> Linked list interface - defines registry_t type
#include "linkedList_i.f"

!> Global registry
type(registry_t) :: ufo_locs_registry

! ------------------------------------------------------------------------------
contains
! ------------------------------------------------------------------------------
!> Linked list implementation
#include "linkedList_c.f"

! ------------------------------------------------------------------------------

subroutine ufo_loc_setup(self, lvec)
implicit none
type(ufo_locs), intent(inout) :: self
type(obs_vector), intent(in) :: lvec
integer :: jc, jo

self%nloc=lvec%nobs
allocate(self%xyz(3,self%nloc))
do jo=1,self%nloc
  do jc=1,3
    self%xyz(jc,jo)=lvec%values(jo)
  enddo
enddo

end subroutine ufo_loc_setup

! ------------------------------------------------------------------------------

subroutine ufo_loc_delete_c(key) bind(c,name='ufo_loc_delete_f90')

implicit none
integer(c_int), intent(inout) :: key
type(ufo_locs), pointer :: self

call ufo_locs_registry%get(key,self)
deallocate(self%xyz)
call ufo_locs_registry%remove(key)

end subroutine ufo_loc_delete_c

! ------------------------------------------------------------------------------

subroutine ufo_loc_nobs_c(key, kobs) bind(c,name='ufo_loc_nobs_f90')

implicit none
integer(c_int), intent(in) :: key
integer(c_int), intent(inout) :: kobs
type(ufo_locs), pointer :: self

call ufo_locs_registry%get(key,self)
kobs = self%nloc

end subroutine ufo_loc_nobs_c

! ------------------------------------------------------------------------------

end module ufo_locs_mod
