! (C) Copyright 2017-2018 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

!> Fortran module to handle atmvertinterp observations

module ufo_atmvertinterp_mod_c

  use iso_c_binding
  use fckit_configuration_module, only: fckit_configuration
  use ufo_atmvertinterp_mod
  use ufo_geovals_mod,    only: ufo_geovals
  use ufo_geovals_mod_c,  only: ufo_geovals_registry

  implicit none

  private

#define LISTED_TYPE ufo_atmvertinterp

  !> Linked list interface - defines registry_t type
#include "oops/util/linkedList_i.f"

  !> Global registry
  type(registry_t) :: ufo_atmvertinterp_registry

  ! ------------------------------------------------------------------------------

contains
  ! ------------------------------------------------------------------------------
  !> Linked list implementation
#include "oops/util/linkedList_c.f"

! ------------------------------------------------------------------------------

subroutine ufo_atmvertinterp_setup_c(c_key_self, c_conf, c_obsvars, c_obsvarindices, c_nobsvars, &
                                     c_geovars) bind(c,name='ufo_atmvertinterp_setup_f90')
use oops_variables_mod
implicit none
integer(c_int), intent(inout)     :: c_key_self
type(c_ptr), intent(in), value    :: c_conf
type(c_ptr), intent(in), value    :: c_obsvars                    ! variables to be simulated...
integer(c_int), intent(in), value :: c_nobsvars
integer(c_int), intent(in)        :: c_obsvarindices(c_nobsvars)  ! ... and their global indices
type(c_ptr), intent(in), value    :: c_geovars  ! variables requested from the model

type(ufo_atmvertinterp), pointer :: self
type(fckit_configuration) :: f_conf

call ufo_atmvertinterp_registry%setup(c_key_self, self)
f_conf = fckit_configuration(c_conf)

self%obsvars = oops_variables(c_obsvars)
allocate(self%obsvarindices(self%obsvars%nvars()))
self%obsvarindices(:) = c_obsvarindices(:) + 1  ! Convert from C to Fortran indexing
self%geovars = oops_variables(c_geovars)

call self%setup(f_conf)

end subroutine ufo_atmvertinterp_setup_c

! ------------------------------------------------------------------------------

subroutine ufo_atmvertinterp_delete_c(c_key_self) bind(c,name='ufo_atmvertinterp_delete_f90')
implicit none
integer(c_int), intent(inout) :: c_key_self

type(ufo_atmvertinterp), pointer :: self

call ufo_atmvertinterp_registry%get(c_key_self, self)

! the obsvarindices array is allocated in the interface layer, so we deallocate here as well
if (allocated(self%obsvarindices)) deallocate(self%obsvarindices)

call ufo_atmvertinterp_registry%remove(c_key_self)

end subroutine ufo_atmvertinterp_delete_c

! ------------------------------------------------------------------------------

subroutine ufo_atmvertinterp_simobs_c(c_key_self, c_key_geovals, c_obsspace, c_nvars, c_nlocs, &
                                      c_hofx) bind(c,name='ufo_atmvertinterp_simobs_f90')
implicit none
integer(c_int), intent(in) :: c_key_self
integer(c_int), intent(in) :: c_key_geovals
type(c_ptr), value, intent(in) :: c_obsspace
integer(c_int), intent(in)     :: c_nvars, c_nlocs
real(c_double), intent(inout)  :: c_hofx(c_nvars, c_nlocs)

type(ufo_atmvertinterp), pointer :: self
type(ufo_geovals),       pointer :: geovals
character(len=*), parameter :: myname_="ufo_atmvertinterp_simobs_c"

call ufo_atmvertinterp_registry%get(c_key_self, self)
call ufo_geovals_registry%get(c_key_geovals, geovals)

call self%simobs(geovals, c_obsspace, c_nvars, c_nlocs, c_hofx)

end subroutine ufo_atmvertinterp_simobs_c

! ------------------------------------------------------------------------------

end module ufo_atmvertinterp_mod_c
