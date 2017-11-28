!
!  (C) Copyright 2017 UCAR
!  
!  This software is licensed under the terms of the Apache Licence Version 2.0
!  which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
!

module ufo_locs_mod_c

use iso_c_binding
use ufo_locs_mod

implicit none

public :: ufo_locs_registry

private

! ------------------------------------------------------------------------------

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

end module ufo_locs_mod_c
