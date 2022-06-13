! (C) Copyright 2017-2018 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

!> Fortran module to handle atmvertinterplay observations

module ufo_atmvertinterplay_tlad_mod_c

  use fckit_configuration_module, only: fckit_configuration
  use ufo_atmvertinterplay_tlad_mod
  use ufo_geovals_mod_c, only: ufo_geovals_registry
  use ufo_geovals_mod,   only: ufo_geovals
  implicit none

  private

#define LISTED_TYPE ufo_atmvertinterplay_tlad

  !> Linked list interface - defines registry_t type
#include "oops/util/linkedList_i.f"

  !> Global registry
  type(registry_t) :: ufo_atmvertinterplay_tlad_registry

  ! ------------------------------------------------------------------------------
contains
  ! ------------------------------------------------------------------------------
  !> Linked list implementation
#include "oops/util/linkedList_c.f"

! ------------------------------------------------------------------------------

subroutine ufo_atmvertinterplay_tlad_setup_c(c_key_self, c_conf, c_obsvars, c_geovars) bind(c,name='ufo_atmvertinterplay_tlad_setup_f90')
use oops_variables_mod
implicit none
integer(c_int), intent(inout) :: c_key_self
type(c_ptr), value, intent(in) :: c_conf
type(c_ptr), value, intent(in) :: c_obsvars ! variables to be simulated
type(c_ptr), value, intent(in) :: c_geovars ! variables requested from the model

type(ufo_atmvertinterplay_tlad), pointer :: self
type(fckit_configuration) :: f_conf

call ufo_atmvertinterplay_tlad_registry%setup(c_key_self, self)
f_conf = fckit_configuration(c_conf)

self%obsvars = oops_variables(c_obsvars)
self%geovars = oops_variables(c_geovars)
call self%setup(f_conf)

end subroutine ufo_atmvertinterplay_tlad_setup_c

! ------------------------------------------------------------------------------

subroutine ufo_atmvertinterplay_tlad_delete_c(c_key_self) bind(c,name='ufo_atmvertinterplay_tlad_delete_f90')
implicit none
integer(c_int), intent(inout) :: c_key_self

type(ufo_atmvertinterplay_tlad), pointer :: self

call ufo_atmvertinterplay_tlad_registry%delete(c_key_self, self)

end subroutine ufo_atmvertinterplay_tlad_delete_c

! ------------------------------------------------------------------------------

subroutine ufo_atmvertinterplay_tlad_settraj_c(c_key_self, c_key_geovals, c_obsspace) bind(c,name='ufo_atmvertinterplay_tlad_settraj_f90')
implicit none
integer(c_int), intent(in) :: c_key_self
integer(c_int), intent(in) :: c_key_geovals
type(c_ptr), value, intent(in) :: c_obsspace

type(ufo_atmvertinterplay_tlad), pointer :: self
type(ufo_geovals),            pointer :: geovals

character(len=*), parameter :: myname_="ufo_atmvertinterplay_tlad_settraj_c"

call ufo_atmvertinterplay_tlad_registry%get(c_key_self, self)
call ufo_geovals_registry%get(c_key_geovals, geovals)

call self%settraj(geovals, c_obsspace)

end subroutine ufo_atmvertinterplay_tlad_settraj_c

! ------------------------------------------------------------------------------

subroutine ufo_atmvertinterplay_simobs_tl_c(c_key_self, c_key_geovals, c_obsspace, c_nvars, c_nlocs, c_hofx) &
           bind(c,name='ufo_atmvertinterplay_simobs_tl_f90')
implicit none
integer(c_int), intent(in) :: c_key_self
integer(c_int), intent(in) :: c_key_geovals
type(c_ptr), value, intent(in) :: c_obsspace
integer(c_int), intent(in) :: c_nvars, c_nlocs
real(c_double), intent(inout) :: c_hofx(c_nvars, c_nlocs)

type(ufo_atmvertinterplay_tlad), pointer :: self
type(ufo_geovals),            pointer :: geovals

character(len=*), parameter :: myname_="ufo_atmvertinterplay_simobs_tl_c"

call ufo_atmvertinterplay_tlad_registry%get(c_key_self, self)
call ufo_geovals_registry%get(c_key_geovals, geovals)

call self%simobs_tl(geovals, c_obsspace, c_nvars, c_nlocs, c_hofx)

end subroutine ufo_atmvertinterplay_simobs_tl_c

! ------------------------------------------------------------------------------

subroutine ufo_atmvertinterplay_simobs_ad_c(c_key_self, c_key_geovals, c_obsspace, c_nvars, c_nlocs, c_hofx) &
           bind(c,name='ufo_atmvertinterplay_simobs_ad_f90')
implicit none
integer(c_int), intent(in) :: c_key_self
integer(c_int), intent(in) :: c_key_geovals
type(c_ptr), value, intent(in) :: c_obsspace
integer(c_int), intent(in) :: c_nvars, c_nlocs
real(c_double), intent(in) :: c_hofx(c_nvars, c_nlocs)

type(ufo_atmvertinterplay_tlad), pointer :: self
type(ufo_geovals),            pointer :: geovals

character(len=*), parameter :: myname_="ufo_atmvertinterplay_simobs_ad_c"

call ufo_atmvertinterplay_tlad_registry%get(c_key_self, self)
call ufo_geovals_registry%get(c_key_geovals, geovals)

call self%simobs_ad(geovals, c_obsspace, c_nvars, c_nlocs, c_hofx)

end subroutine ufo_atmvertinterplay_simobs_ad_c

! ------------------------------------------------------------------------------

end module ufo_atmvertinterplay_tlad_mod_c
