! (C) Copyright 2017-2018 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

!> Fortran atmvertinterplay module for functions on the interface between C++ and Fortran
!  to handle observation operators

module ufo_atmvertinterplay_mod_c

  use iso_c_binding

  use fckit_configuration_module, only: fckit_configuration
  use ufo_atmvertinterplay_mod
  use ufo_geovals_mod,   only: ufo_geovals
  use ufo_geovals_mod_c, only: ufo_geovals_registry
  implicit none
  private

  ! ------------------------------------------------------------------------------
#define LISTED_TYPE ufo_atmvertinterplay

  !> Linked list interface - defines registry_t type
#include "oops/util/linkedList_i.f"

  !> Global registry
  type(registry_t) :: ufo_atmvertinterplay_registry

  ! ------------------------------------------------------------------------------

contains

  ! ------------------------------------------------------------------------------
  !> Linked list implementation
#include "oops/util/linkedList_c.f"

! ------------------------------------------------------------------------------

subroutine ufo_atmvertinterplay_setup_c(c_key_self, c_conf, c_obsvars, c_geovars) bind(c,name='ufo_atmvertinterplay_setup_f90')
use oops_variables_mod
use obs_variables_mod
implicit none
integer(c_int), intent(inout)  :: c_key_self
type(c_ptr), value, intent(in) :: c_conf
type(c_ptr), value, intent(in) :: c_obsvars ! variables to be simulated
type(c_ptr), value, intent(in) :: c_geovars ! variables requested from the model

type(ufo_atmvertinterplay), pointer :: self
type(fckit_configuration) :: f_conf

call ufo_atmvertinterplay_registry%setup(c_key_self, self)

f_conf = fckit_configuration(c_conf)

self%obsvars = obs_variables(c_obsvars)
self%geovars = oops_variables(c_geovars)

call self%setup(f_conf)

end subroutine ufo_atmvertinterplay_setup_c

! ------------------------------------------------------------------------------

subroutine ufo_atmvertinterplay_delete_c(c_key_self) bind(c,name='ufo_atmvertinterplay_delete_f90')
implicit none
integer(c_int), intent(inout) :: c_key_self

type(ufo_atmvertinterplay), pointer :: self

call ufo_atmvertinterplay_registry%get(c_key_self, self)
call self%delete
call ufo_atmvertinterplay_registry%remove(c_key_self)

end subroutine ufo_atmvertinterplay_delete_c

! ------------------------------------------------------------------------------

subroutine ufo_atmvertinterplay_simobs_c(c_key_self, c_key_geovals, c_obsspace, c_nvars, c_nlocs, &
                                c_hofx) bind(c,name='ufo_atmvertinterplay_simobs_f90')

implicit none
integer(c_int), intent(in) :: c_key_self
integer(c_int), intent(in) :: c_key_geovals
type(c_ptr), value, intent(in) :: c_obsspace
integer(c_int), intent(in)     :: c_nvars, c_nlocs
real(c_double), intent(inout)  :: c_hofx(c_nvars, c_nlocs)

type(ufo_atmvertinterplay), pointer :: self
type(ufo_geovals),       pointer :: geovals
character(len=*), parameter :: myname_="ufo_atmvertinterplay_simobs_c"

call ufo_atmvertinterplay_registry%get(c_key_self, self)
call ufo_geovals_registry%get(c_key_geovals, geovals)
call self%simobs(geovals, c_obsspace, c_nvars, c_nlocs, c_hofx)

end subroutine ufo_atmvertinterplay_simobs_c

! ------------------------------------------------------------------------------

end module ufo_atmvertinterplay_mod_c
