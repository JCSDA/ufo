!
! (C) Copyright 2017-2018 UCAR
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
!
module ufo_bgcheck_mod_c

use iso_c_binding
use ufo_bgcheck_mod
use ufo_geovals_mod
use ufo_geovals_mod_c,   only: ufo_geovals_registry

implicit none
private

#define LISTED_TYPE ufo_bgcheck

!> Linked list interface - defines registry_t type
#include "linkedList_i.f"

!> Global registry
type(registry_t) :: ufo_bgcheck_registry

! ------------------------------------------------------------------------------
contains
! ------------------------------------------------------------------------------
!> Linked list implementation
#include "linkedList_c.f"
! ------------------------------------------------------------------------------

subroutine ufo_bgcheck_create_c(c_self, c_obspace, c_conf) bind(c,name='ufo_bgcheck_create_f90')
implicit none
integer(c_int), intent(inout)  :: c_self
type(c_ptr), value, intent(in) :: c_obspace
type(c_ptr), value, intent(in) :: c_conf

type(ufo_bgcheck), pointer :: self

call ufo_bgcheck_registry%setup(c_self, self)
call ufo_bgcheck_create(self, c_obspace, c_conf)

end subroutine ufo_bgcheck_create_c

! ------------------------------------------------------------------------------

subroutine ufo_bgcheck_delete_c(c_self) bind(c,name='ufo_bgcheck_delete_f90')
implicit none
integer(c_int), intent(inout) :: c_self

type(ufo_bgcheck), pointer :: self

call ufo_bgcheck_registry%get(c_self, self)
call ufo_bgcheck_delete(self)
call ufo_bgcheck_registry%delete(c_self, self)

end subroutine ufo_bgcheck_delete_c

! ------------------------------------------------------------------------------

subroutine ufo_bgcheck_prior_c(c_self, c_geovals) bind(c,name='ufo_bgcheck_prior_f90')
implicit none
integer(c_int), intent(in) :: c_self
integer(c_int), intent(in) :: c_geovals

type(ufo_bgcheck), pointer :: self
type(ufo_geovals), pointer :: geovals

call ufo_bgcheck_registry%get(c_self, self)
call ufo_geovals_registry%get(c_geovals, geovals)

call ufo_bgcheck_prior(self, geovals)

end subroutine ufo_bgcheck_prior_c

! ------------------------------------------------------------------------------

subroutine ufo_bgcheck_post_c(c_self, c_nobs, c_hofx) bind(c,name='ufo_bgcheck_post_f90')
implicit none
integer(c_int), intent(in) :: c_self
integer(c_int), intent(in) :: c_nobs
real(c_double), intent(in) :: c_hofx(c_nobs)

type(ufo_bgcheck), pointer :: self

call ufo_bgcheck_registry%get(c_self, self)

call ufo_bgcheck_post(self, c_hofx)

end subroutine ufo_bgcheck_post_c

! ------------------------------------------------------------------------------

end module ufo_bgcheck_mod_c
