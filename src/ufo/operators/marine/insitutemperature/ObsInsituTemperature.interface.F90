! (C) Copyright 2017-2022 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

!> Fortran insitutemperature module for functions on the interface between C++ and Fortran
!  to handle observation operators

module ufo_insitutemperature_mod_c

  use iso_c_binding
  use ufo_insitutemperature_mod
  use ufo_geovals_mod
  use ufo_geovals_mod_c,   only: ufo_geovals_registry
  use oops_variables_mod
  use obs_variables_mod
  implicit none
  private

  ! ------------------------------------------------------------------------------
#define LISTED_TYPE ufo_insitutemperature

  !> Linked list interface - defines registry_t type
#include "oops/util/linkedList_i.f"

  !> Global registry
  type(registry_t) :: ufo_insitutemperature_registry

  ! ------------------------------------------------------------------------------

contains

  ! ------------------------------------------------------------------------------
  !> Linked list implementation
#include "oops/util/linkedList_c.f"

! ------------------------------------------------------------------------------

subroutine ufo_insitutemperature_setup_c(c_key_self, c_conf, c_obsvars, c_obsvarindices, &
    c_nobsvars, c_geovars) bind(c,name='ufo_insitutemperature_setup_f90')
  integer(c_int),        intent(inout) :: c_key_self
  type(c_ptr),    value, intent(in)    :: c_conf
  type(c_ptr),    value, intent(in)    :: c_obsvars ! variables to be simulated...
  integer(c_int), value, intent(in)    :: c_nobsvars
  integer(c_int),        intent(in)    :: c_obsvarindices(c_nobsvars) ! ... and their global indices
  type(c_ptr),    value, intent(in)    :: c_geovars ! variables requested from the model

  type(ufo_insitutemperature), pointer :: self

  call ufo_insitutemperature_registry%setup(c_key_self, self)

  ! assuming that there is only 1 variable
  if (c_nobsvars /= 1) call abor1_ftn("InsituTemperature only works on a single variable")
  self%obsvars = obs_variables(c_obsvars)
  self%obsvaridx = c_obsvarindices(1) + 1
  self%geovars = oops_variables(c_geovars)

  call self%setup()

end subroutine ufo_insitutemperature_setup_c

! ------------------------------------------------------------------------------

subroutine ufo_insitutemperature_delete_c(c_key_self) bind(c,name='ufo_insitutemperature_delete_f90')
  integer(c_int), intent(inout) :: c_key_self

  type(ufo_insitutemperature), pointer :: self

  call ufo_insitutemperature_registry%delete(c_key_self, self)

end subroutine ufo_insitutemperature_delete_c

! ------------------------------------------------------------------------------

subroutine ufo_insitutemperature_simobs_c(c_key_self, c_key_geovals, c_obsspace, c_nvars, &
                                          c_nlocs, c_hofx) &
                                          bind(c,name='ufo_insitutemperature_simobs_f90')
  integer(c_int),        intent(in) :: c_key_self
  integer(c_int),        intent(in) :: c_key_geovals
  type(c_ptr),    value, intent(in) :: c_obsspace
  integer(c_int),        intent(in) :: c_nvars, c_nlocs
  real(c_double),     intent(inout) :: c_hofx(c_nvars, c_nlocs)

  type(ufo_insitutemperature), pointer :: self
  type(ufo_geovals), pointer :: geovals

  call ufo_insitutemperature_registry%get(c_key_self, self)
  call ufo_geovals_registry%get(c_key_geovals,geovals)
  call self%simobs(geovals, c_obsspace, c_nvars, c_nlocs, c_hofx)

end subroutine ufo_insitutemperature_simobs_c

! ------------------------------------------------------------------------------

end module ufo_insitutemperature_mod_c
