!
!  (C) Copyright 2017 UCAR
!
!  This software is licensed under the terms of the Apache Licence Version 2.0
!  which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
!

module ufo_locs_mod_c

use iso_c_binding
use ufo_locs_mod
use kinds

implicit none

public :: ufo_locs_registry

private

! ------------------------------------------------------------------------------

#define LISTED_TYPE ufo_locs

!> Linked list interface - defines registry_t type
#include "oops/util/linkedList_i.f"

!> Global registry
type(registry_t) :: ufo_locs_registry

! ------------------------------------------------------------------------------
contains
! ------------------------------------------------------------------------------
!> Linked list implementation
#include "oops/util/linkedList_c.f"

! ------------------------------------------------------------------------------

subroutine ufo_locs_create_c(key, klocs, c_obsspace, klats, klons) bind(c,name='ufo_locs_create_f90')

implicit none
integer(c_int), intent(inout) :: key
integer(c_int), intent(in) :: klocs
type(c_ptr), value, intent(in) :: c_obsspace
real(c_double), intent(in) :: klats(klocs)
real(c_double), intent(in) :: klons(klocs)

type(ufo_locs), pointer :: self
real(kind_real) :: lats(klocs)
real(kind_real) :: lons(klocs)

call ufo_locs_registry%setup(key, self)

lats(:) = klats(:)
lons(:) = klons(:)

call ufo_locs_create(self, c_obsspace, klocs, lats, lons)

end subroutine ufo_locs_create_c

! ------------------------------------------------------------------------------

subroutine ufo_locs_setup_c(key, nlocs) bind(c,name='ufo_locs_setup_f90')

implicit none
integer(c_int), intent(inout) :: key
integer(c_int), intent(in) :: nlocs
type(ufo_locs), pointer :: self

call ufo_locs_registry%setup(key, self)

call ufo_locs_setup(self, nlocs)

end subroutine ufo_locs_setup_c

!------------------------------------------------------------------------------

subroutine ufo_locs_copy_c(key, key2) bind(c,name='ufo_locs_copy_f90')

implicit none
integer(c_int), intent(inout) :: key
integer(c_int), intent(in)    :: key2

type(ufo_locs), pointer :: self
type(ufo_locs), pointer :: other

call ufo_locs_registry%setup(key, self)
call ufo_locs_registry%get(key2, other)

call ufo_locs_copy(self, other)

end subroutine ufo_locs_copy_c

! ------------------------------------------------------------------------------

subroutine ufo_locs_delete_c(key) bind(c,name='ufo_locs_delete_f90')

implicit none
integer(c_int), intent(inout) :: key
type(ufo_locs), pointer :: self

call ufo_locs_registry%get(key,self)
call ufo_locs_delete(self)
call ufo_locs_registry%remove(key)

end subroutine ufo_locs_delete_c

! ------------------------------------------------------------------------------

subroutine ufo_locs_nobs_c(key, kobs) bind(c,name='ufo_locs_nobs_f90')

implicit none
integer(c_int), intent(in) :: key
integer(c_int), intent(inout) :: kobs
type(ufo_locs), pointer :: self

call ufo_locs_registry%get(key,self)
kobs = self%nlocs

end subroutine ufo_locs_nobs_c
! ------------------------------------------------------------------------------

subroutine ufo_locs_coords_c(key, idx, mylat, mylon) bind(c,name='ufo_locs_coords_f90')

implicit none
integer(c_int), intent(in) :: key
integer(c_int), intent(in) :: idx
real(c_double), intent(inout) :: mylat,mylon

type(ufo_locs), pointer :: self

call ufo_locs_registry%get(key,self)
mylat = self%lat(idx+1)
mylon = self%lon(idx+1)

end subroutine ufo_locs_coords_c
!---------------------------------------------------------------------------------

subroutine ufo_locs_indx_c(key, idx, indx, max_indx) bind(c,name='ufo_locs_indx_f90')

implicit none
integer(c_int), intent(in) :: key
integer(c_int), intent(in) :: idx
integer(c_int), intent(inout) :: indx
integer(c_int), intent(inout) :: max_indx

type(ufo_locs), pointer :: self

call ufo_locs_registry%get(key, self)
max_indx = self%max_indx
if (max_indx > 0) indx = self%indx(idx+1) - 1 ! the minus to take account of C++ starting from 0


end subroutine ufo_locs_indx_c
!------------------------------------------------------------------------------
subroutine ufo_locs_concatenate_c(key, key2) bind(c,name='ufo_locs_concatenate_f90')

implicit none
integer(c_int), intent(in) :: key
integer(c_int), intent(in) :: key2

type(ufo_locs), pointer :: self
type(ufo_locs), pointer :: other

call ufo_locs_registry%get(key, self)
call ufo_locs_registry%get(key2, other)
call ufo_locs_concatenate(self, other)

end subroutine ufo_locs_concatenate_c
! ------------------------------------------------------------------------------

subroutine ufo_locs_init_c(c_key_self, c_obsspace) bind(c,name='ufo_locs_init_f90')
implicit none
integer(c_int), intent(inout)  :: c_key_self
type(c_ptr), value, intent(in) :: c_obsspace

type(ufo_locs), pointer :: self

call ufo_locs_registry%setup(c_key_self, self)

call ufo_locs_registry%get(c_key_self, self)
call ufo_locs_init(self, c_obsspace)

end subroutine ufo_locs_init_c

end module ufo_locs_mod_c
