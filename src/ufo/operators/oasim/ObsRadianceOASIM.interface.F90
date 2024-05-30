! (C) Copyright 2022 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

!> Fortran oasim module for functions on the interface between C++ and Fortran
!  to handle observation operators

module ufo_radianceoasim_mod_c

  use fckit_configuration_module, only: fckit_configuration
  use iso_c_binding
  use ufo_oasim_mod
  use ufo_geovals_mod
  use ufo_geovals_mod_c,   only: ufo_geovals_registry

  implicit none
  private

  ! ------------------------------------------------------------------------------
#define LISTED_TYPE ufo_oasim

  !> Linked list interface - defines registry_t type
#include "oops/util/linkedList_i.f"

  !> Global registry
  type(registry_t) :: ufo_radianceoasim_registry

  ! ------------------------------------------------------------------------------

contains

  ! ------------------------------------------------------------------------------
  !> Linked list implementation
#include "oops/util/linkedList_c.f"

! ------------------------------------------------------------------------------

subroutine ufo_radianceoasim_setup_c(c_key_self, c_conf, c_nchan, c_channels) bind(c,name='ufo_radianceoasim_setup_f90')

integer(c_int), intent(inout)  :: c_key_self
type(c_ptr), value, intent(in) :: c_conf
integer(c_int), intent(in) :: c_nchan
integer(c_int), intent(in) :: c_channels(c_nchan)

type(ufo_oasim), pointer :: self
type(fckit_configuration) :: f_conf

call ufo_radianceoasim_registry%setup(c_key_self, self)
f_conf = fckit_configuration(c_conf)

call self%setup(f_conf, c_channels)

end subroutine ufo_radianceoasim_setup_c

! ------------------------------------------------------------------------------

subroutine ufo_radianceoasim_delete_c(c_key_self) bind(c,name='ufo_radianceoasim_delete_f90')
integer(c_int), intent(inout) :: c_key_self

type(ufo_oasim), pointer :: self

call ufo_radianceoasim_registry%get(c_key_self, self)
call self%delete()
call ufo_radianceoasim_registry%remove(c_key_self)

end subroutine ufo_radianceoasim_delete_c

! ------------------------------------------------------------------------------

subroutine ufo_radianceoasim_simobs_c(c_key_self, c_key_geovals, c_obsspace, c_nvars, c_nlocs, c_hofx) &
  bind(c,name='ufo_radianceoasim_simobs_f90')

integer(c_int), intent(in) :: c_key_self
integer(c_int), intent(in) :: c_key_geovals
type(c_ptr), value, intent(in) :: c_obsspace
integer(c_int), intent(in) :: c_nlocs, c_nvars
real(c_double), intent(inout) :: c_hofx(c_nvars,c_nlocs)

type(ufo_oasim), pointer :: self
type(ufo_geovals), pointer :: geovals

call ufo_radianceoasim_registry%get(c_key_self, self)
call ufo_geovals_registry%get(c_key_geovals,geovals)

call self%simobs(geovals, c_obsspace, c_nvars, c_nlocs, c_hofx)

end subroutine ufo_radianceoasim_simobs_c

end module ufo_radianceoasim_mod_c
