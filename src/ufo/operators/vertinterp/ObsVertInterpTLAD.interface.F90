! (C) Copyright 2017-2023 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

!> Fortran module to handle vertinterp observations

module ufo_vertinterp_tlad_mod_c

  use fckit_configuration_module, only: fckit_configuration
  use ufo_vertinterp_tlad_mod
  use ufo_geovals_mod_c, only: ufo_geovals_registry
  use ufo_geovals_mod,   only: ufo_geovals
  implicit none

  private

#define LISTED_TYPE ufo_vertinterp_tlad

  !> Linked list interface - defines registry_t type
#include "oops/util/linkedList_i.f"

  !> Global registry
  type(registry_t) :: ufo_vertinterp_tlad_registry

  ! ------------------------------------------------------------------------------
contains
  ! ------------------------------------------------------------------------------
  !> Linked list implementation
#include "oops/util/linkedList_c.f"

! ------------------------------------------------------------------------------

subroutine ufo_vertinterp_tlad_setup_c(c_key_self, c_conf, c_obsvars, &
                                          c_obsvarindices, c_nobsvars, c_geovars) bind(c,name='ufo_vertinterp_tlad_setup_f90')
use oops_variables_mod
use obs_variables_mod
implicit none
integer(c_int), intent(inout)     :: c_key_self
type(c_ptr), value, intent(in)    :: c_conf
type(c_ptr), intent(in), value    :: c_obsvars                    ! variables to be simulated...
integer(c_int), intent(in), value :: c_nobsvars
integer(c_int), intent(in)        :: c_obsvarindices(c_nobsvars)  ! ... and their global indices
type(c_ptr), value, intent(in)    :: c_geovars                    ! variables requested from the model

type(ufo_vertinterp_tlad), pointer :: self
type(fckit_configuration) :: f_conf

call ufo_vertinterp_tlad_registry%setup(c_key_self, self)
f_conf = fckit_configuration(c_conf)

self%obsvars = obs_variables(c_obsvars)
allocate(self%obsvarindices(self%obsvars%nvars()))
self%obsvarindices(:) = c_obsvarindices(:) + 1  ! Convert from C to Fortran indexing
self%geovars = oops_variables(c_geovars)
call self%setup(f_conf)

end subroutine ufo_vertinterp_tlad_setup_c

! ------------------------------------------------------------------------------

subroutine ufo_vertinterp_tlad_delete_c(c_key_self) bind(c,name='ufo_vertinterp_tlad_delete_f90')
implicit none
integer(c_int), intent(inout) :: c_key_self

type(ufo_vertinterp_tlad), pointer :: self

call ufo_vertinterp_tlad_registry%get(c_key_self, self)

! the obsvarindices array is allocated in the interface layer, so we deallocate here as well
if (allocated(self%obsvarindices)) deallocate(self%obsvarindices)

! type ufo_vertinterp_tlad has allocatable data, but has a destructor marked final that
! should automatically deallocate (assuming compiler support...)
call ufo_vertinterp_tlad_registry%remove(c_key_self)

end subroutine ufo_vertinterp_tlad_delete_c

! ------------------------------------------------------------------------------

subroutine ufo_vertinterp_tlad_settraj_c(c_key_self, c_key_geovals, c_obsspace) bind(c,name='ufo_vertinterp_tlad_settraj_f90')
implicit none
integer(c_int), intent(in) :: c_key_self
integer(c_int), intent(in) :: c_key_geovals
type(c_ptr), value, intent(in) :: c_obsspace

type(ufo_vertinterp_tlad), pointer :: self
type(ufo_geovals),            pointer :: geovals

character(len=*), parameter :: myname_="ufo_vertinterp_tlad_settraj_c"

call ufo_vertinterp_tlad_registry%get(c_key_self, self)
call ufo_geovals_registry%get(c_key_geovals, geovals)

call self%settraj(geovals, c_obsspace)

end subroutine ufo_vertinterp_tlad_settraj_c

! ------------------------------------------------------------------------------

subroutine ufo_vertinterp_simobs_tl_c(c_key_self, c_key_geovals, c_obsspace, c_nvars, c_nlocs, c_hofx) &
           bind(c,name='ufo_vertinterp_simobs_tl_f90')
implicit none
integer(c_int), intent(in) :: c_key_self
integer(c_int), intent(in) :: c_key_geovals
type(c_ptr), value, intent(in) :: c_obsspace
integer(c_int), intent(in) :: c_nvars, c_nlocs
real(c_double), intent(inout) :: c_hofx(c_nvars, c_nlocs)

type(ufo_vertinterp_tlad), pointer :: self
type(ufo_geovals),            pointer :: geovals

character(len=*), parameter :: myname_="ufo_vertinterp_simobs_tl_c"

call ufo_vertinterp_tlad_registry%get(c_key_self, self)
call ufo_geovals_registry%get(c_key_geovals, geovals)

call self%simobs_tl(geovals, c_obsspace, c_nvars, c_nlocs, c_hofx)

end subroutine ufo_vertinterp_simobs_tl_c

! ------------------------------------------------------------------------------

subroutine ufo_vertinterp_simobs_ad_c(c_key_self, c_key_geovals, c_obsspace, c_nvars, c_nlocs, c_hofx) &
           bind(c,name='ufo_vertinterp_simobs_ad_f90')
implicit none
integer(c_int), intent(in) :: c_key_self
integer(c_int), intent(in) :: c_key_geovals
type(c_ptr), value, intent(in) :: c_obsspace
integer(c_int), intent(in) :: c_nvars, c_nlocs
real(c_double), intent(in) :: c_hofx(c_nvars, c_nlocs)

type(ufo_vertinterp_tlad), pointer :: self
type(ufo_geovals),            pointer :: geovals

character(len=*), parameter :: myname_="ufo_vertinterp_simobs_ad_c"

call ufo_vertinterp_tlad_registry%get(c_key_self, self)
call ufo_geovals_registry%get(c_key_geovals, geovals)

call self%simobs_ad(geovals, c_obsspace, c_nvars, c_nlocs, c_hofx)

end subroutine ufo_vertinterp_simobs_ad_c

! ------------------------------------------------------------------------------

end module ufo_vertinterp_tlad_mod_c
