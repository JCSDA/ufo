! (C) Copyright 2021.
!
! This software is developed by NOAA/NWS/EMC under the Apache 2.0 license
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

!> Fortran insitupm module for functions on the interface between C++ and Fortran
!  to handle tl/ad observation operators

module ufo_insitupm_tlad_mod_c

  use iso_c_binding
  use ufo_insitupm_tlad_mod
  implicit none
  private

#define LISTED_TYPE ufo_insitupm_tlad

  !> Linked list interface - defines registry_t type
#include "oops/util/linkedList_i.f"

  !> Global registry
  type(registry_t) :: ufo_insitupm_tlad_registry

contains

  ! ------------------------------------------------------------------------------
  !> Linked list implementation
#include "oops/util/linkedList_c.f"

! ------------------------------------------------------------------------------

subroutine ufo_insitupm_tlad_setup_c(c_key_self, c_conf, c_obsvars, c_geovars) bind(c,name='ufo_insitupm_tlad_setup_f90')
use fckit_configuration_module, only: fckit_configuration
use oops_variables_mod
use obs_variables_mod
implicit none
integer(c_int), intent(inout)  :: c_key_self
type(c_ptr), value, intent(in) :: c_conf
type(c_ptr), value, intent(in) :: c_obsvars ! variables to be simulated
type(c_ptr), value, intent(in) :: c_geovars ! variables requested from the model

type(ufo_insitupm_tlad), pointer :: self
type(fckit_configuration) :: f_conf

call ufo_insitupm_tlad_registry%setup(c_key_self, self)
f_conf = fckit_configuration(c_conf)

self%obsvars = obs_variables(c_obsvars)
self%geovars = oops_variables(c_geovars)

call self%setup(f_conf)

end subroutine ufo_insitupm_tlad_setup_c

! ------------------------------------------------------------------------------

subroutine ufo_insitupm_tlad_delete_c(c_key_self) bind(c,name='ufo_insitupm_tlad_delete_f90')
implicit none
integer(c_int), intent(inout) :: c_key_self

! type ufo_insitupm_tlad has allocatable data, but has a destructor marked final that
! should automatically deallocate (assuming compiler support...)
call ufo_insitupm_tlad_registry%remove(c_key_self)

end subroutine ufo_insitupm_tlad_delete_c

! ------------------------------------------------------------------------------

subroutine ufo_insitupm_tlad_settraj_c(c_key_self, c_key_geovals, c_obsspace) bind(c,name='ufo_insitupm_tlad_settraj_f90')
use ufo_geovals_mod_c, only: ufo_geovals_registry
use ufo_geovals_mod,   only: ufo_geovals
implicit none
integer(c_int),     intent(in) :: c_key_self
integer(c_int),     intent(in) :: c_key_geovals
type(c_ptr), value, intent(in) :: c_obsspace

type(ufo_insitupm_tlad), pointer :: self
type(ufo_geovals),      pointer :: geovals
type(ufo_geovals),      pointer :: hofxdiags

call ufo_insitupm_tlad_registry%get(c_key_self, self)
call ufo_geovals_registry%get(c_key_geovals, geovals)

call self%settraj(geovals, c_obsspace)

end subroutine ufo_insitupm_tlad_settraj_c

! ------------------------------------------------------------------------------

subroutine ufo_insitupm_simobs_tl_c(c_key_self, c_key_geovals, c_obsspace, c_nvars, c_nlocs, c_hofx) bind(c,name='ufo_insitupm_simobs_tl_f90')
use ufo_geovals_mod_c, only: ufo_geovals_registry
use ufo_geovals_mod,   only: ufo_geovals
implicit none
integer(c_int), intent(in) :: c_key_self
integer(c_int), intent(in) :: c_key_geovals
type(c_ptr), value, intent(in) :: c_obsspace
integer(c_int), intent(in) :: c_nvars, c_nlocs
real(c_double), intent(inout) :: c_hofx(c_nvars, c_nlocs)

type(ufo_insitupm_tlad), pointer :: self
type(ufo_geovals),      pointer :: geovals

call ufo_insitupm_tlad_registry%get(c_key_self, self)
call ufo_geovals_registry%get(c_key_geovals, geovals)
call self%simobs_tl(geovals, c_obsspace, c_nvars, c_nlocs, c_hofx)

end subroutine ufo_insitupm_simobs_tl_c

! ------------------------------------------------------------------------------

subroutine ufo_insitupm_simobs_ad_c(c_key_self, c_key_geovals, c_obsspace, c_nvars, c_nlocs, c_hofx) bind(c,name='ufo_insitupm_simobs_ad_f90')
use ufo_geovals_mod_c, only: ufo_geovals_registry
use ufo_geovals_mod,   only: ufo_geovals
implicit none
integer(c_int), intent(in) :: c_key_self
integer(c_int), intent(in) :: c_key_geovals
type(c_ptr), value, intent(in) :: c_obsspace
integer(c_int), intent(in) :: c_nvars, c_nlocs
real(c_double), intent(in) :: c_hofx(c_nvars, c_nlocs)

type(ufo_insitupm_tlad), pointer :: self
type(ufo_geovals),      pointer :: geovals

call ufo_insitupm_tlad_registry%get(c_key_self, self)
call ufo_geovals_registry%get(c_key_geovals, geovals)
call self%simobs_ad(geovals, c_obsspace, c_nvars, c_nlocs, c_hofx)

end subroutine ufo_insitupm_simobs_ad_c

! ------------------------------------------------------------------------------


end module ufo_insitupm_tlad_mod_c
