!
! (C) Copyright 2017-2018 UCAR
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
!
module ufo_rogeorealitycheck_mod_c

use iso_c_binding
use ufo_rogeorealitycheck_mod
use ufo_geovals_mod
use ufo_geovals_mod_c,   only: ufo_geovals_registry

implicit none
private

#define LISTED_TYPE ufo_rogeorealitycheck

!> Linked list interface - defines registry_t type
#include "../../linkedList_i.f"

!> Global registry
type(registry_t) :: ufo_rogeorealitycheck_registry

! ------------------------------------------------------------------------------
contains
! ------------------------------------------------------------------------------
!> Linked list implementation
#include "../../linkedList_c.f"
! ------------------------------------------------------------------------------

subroutine ufo_rogeorealitycheck_create_c(c_self, c_obspace, c_conf) bind(c,name='ufo_rogeorealitycheck_create_f90')
implicit none
integer(c_int), intent(inout)  :: c_self
type(c_ptr), value, intent(in) :: c_obspace
type(c_ptr), value, intent(in) :: c_conf

type(ufo_rogeorealitycheck), pointer :: self
 
call ufo_rogeorealitycheck_registry%setup(c_self, self)
call ufo_rogeorealitycheck_create(self, c_obspace, c_conf)

end subroutine ufo_rogeorealitycheck_create_c

! ------------------------------------------------------------------------------

subroutine ufo_rogeorealitycheck_delete_c(c_self) bind(c,name='ufo_rogeorealitycheck_delete_f90')
implicit none
integer(c_int), intent(inout) :: c_self

type(ufo_rogeorealitycheck), pointer :: self

call ufo_rogeorealitycheck_registry%get(c_self, self)
call ufo_rogeorealitycheck_delete(self)
call ufo_rogeorealitycheck_registry%delete(c_self, self)

end subroutine ufo_rogeorealitycheck_delete_c

! ------------------------------------------------------------------------------

subroutine ufo_rogeorealitycheck_prior_c(c_self, c_geovals) bind(c,name='ufo_rogeorealitycheck_prior_f90')
implicit none
integer(c_int), intent(in) :: c_self
integer(c_int), intent(in) :: c_geovals

type(ufo_rogeorealitycheck), pointer :: self
type(ufo_geovals), pointer :: geovals

call ufo_rogeorealitycheck_registry%get(c_self, self)
call ufo_geovals_registry%get(c_geovals, geovals)
call ufo_rogeorealitycheck_prior(self, geovals)

end subroutine ufo_rogeorealitycheck_prior_c

! ------------------------------------------------------------------------------

subroutine ufo_rogeorealitycheck_post_c(c_self, c_nobs, c_hofx) bind(c,name='ufo_rogeorealitycheck_post_f90')
implicit none
integer(c_int), intent(in) :: c_self
integer(c_int), intent(in) :: c_nobs
real(c_double), intent(in) :: c_hofx(c_nobs)

type(ufo_rogeorealitycheck), pointer :: self

call ufo_rogeorealitycheck_registry%get(c_self, self)

call ufo_rogeorealitycheck_post(self, c_hofx)

end subroutine ufo_rogeorealitycheck_post_c

! ------------------------------------------------------------------------------

end module ufo_rogeorealitycheck_mod_c
