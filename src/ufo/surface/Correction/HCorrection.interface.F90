!
! (C) Copyright 2017-2018 UCAR
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
!
module ufo_hcorrection_mod_c

use fckit_configuration_module, only: fckit_configuration
use iso_c_binding
use ufo_hcorrection_mod
use ufo_geovals_mod
use ufo_geovals_mod_c,   only: ufo_geovals_registry
implicit none
private

#define LISTED_TYPE ufo_hcorrection

!> Linked list interface - defines registry_t type
#include "../../linkedList_i.f"

!> Global registry
type(registry_t) :: ufo_hcorrection_registry

! ------------------------------------------------------------------------------
contains
! ------------------------------------------------------------------------------
!> Linked list implementation
#include "../../linkedList_c.f"
! ------------------------------------------------------------------------------

subroutine ufo_hcorrection_create_c(c_self, c_conf, c_varlist) bind(c,name='ufo_hcorrection_create_f90')
use string_f_c_mod
implicit none
integer(c_int), intent(inout)  :: c_self
type(c_ptr), value, intent(in) :: c_conf
type(c_ptr), intent(in), value :: c_varlist ! list of geovals variables to be requested

type(ufo_hcorrection), pointer :: self
type(fckit_configuration)      :: f_conf

f_conf = fckit_configuration(c_conf)
call ufo_hcorrection_registry%setup(c_self, self)
call ufo_hcorrection_create(self, f_conf)

!> Update C++ ObsFilter with geovals variables list
if (allocated(self%geovars)) then
  call f_c_push_string_varlist(c_varlist, self%geovars)
endif

end subroutine ufo_hcorrection_create_c

! ------------------------------------------------------------------------------

subroutine ufo_hcorrection_delete_c(c_self) bind(c,name='ufo_hcorrection_delete_f90')
implicit none
integer(c_int), intent(inout) :: c_self

type(ufo_hcorrection), pointer :: self

call ufo_hcorrection_registry%get(c_self, self)
call ufo_hcorrection_delete(self)
call ufo_hcorrection_registry%delete(c_self, self)

end subroutine ufo_hcorrection_delete_c

! ------------------------------------------------------------------------------

subroutine ufo_hcorrection_prior_c(c_self, c_obspace, c_geovals) bind(c,name='ufo_hcorrection_prior_f90')
implicit none
integer(c_int), intent(in) :: c_self
type(c_ptr), value, intent(in) :: c_obspace
integer(c_int), intent(in) :: c_geovals

type(ufo_hcorrection), pointer :: self
type(ufo_geovals), pointer :: geovals

call ufo_hcorrection_registry%get(c_self, self)
call ufo_geovals_registry%get(c_geovals, geovals)

call ufo_hcorrection_prior(self, c_obspace, geovals)

end subroutine ufo_hcorrection_prior_c

! ------------------------------------------------------------------------------

subroutine ufo_hcorrection_post_c(c_self, c_obspace, c_nvars, c_nlocs, c_hofx) bind(c,name='ufo_hcorrection_post_f90')
implicit none
integer(c_int), intent(in) :: c_self
type(c_ptr), value, intent(in) :: c_obspace
integer(c_int), intent(in) :: c_nvars, c_nlocs
real(c_double), intent(in) :: c_hofx(c_nvars, c_nlocs)

type(ufo_hcorrection), pointer :: self

call ufo_hcorrection_registry%get(c_self, self)

call ufo_hcorrection_post(self, c_obspace, c_nvars, c_nlocs, c_hofx)

end subroutine ufo_hcorrection_post_c

! ------------------------------------------------------------------------------

end module ufo_hcorrection_mod_c
