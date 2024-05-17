! (C) Copyright 2017-2018 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

!> Fortran lightning module for functions on the interface between C++ and Fortran
!  to handle tl/ad observation operators

module ufo_lightning_tlad_mod_c

  use iso_c_binding
  use ufo_lightning_tlad_mod
  use ufo_geovals_mod_c, only: ufo_geovals_registry
  use ufo_geovals_mod,   only: ufo_geovals
  implicit none
  private

#define LISTED_TYPE ufo_lightning_tlad

  !> Linked list interface - defines registry_t type
#include "oops/util/linkedList_i.f"

  !> Global registry
  type(registry_t) :: ufo_lightning_tlad_registry

contains

  ! ------------------------------------------------------------------------------
  !> Linked list implementation
#include "oops/util/linkedList_c.f"

! ------------------------------------------------------------------------------

subroutine ufo_lightning_tlad_setup_c(c_key_self, nhoriz_, c_obsvars, c_geovars) &
                                      bind(c,name='ufo_lightning_tlad_setup_f90')
use oops_variables_mod
use obs_variables_mod
use ufo_vars_mod, only: MAXVARLEN, var_qg, var_delp
implicit none
integer(c_int), intent(inout)  :: c_key_self
integer(c_int), intent(in)     :: nhoriz_
!logical(c_bool), intent(in)    :: l_fed_nonlinear_
type(c_ptr), value, intent(in) :: c_obsvars ! variables to be simulated
type(c_ptr), value, intent(in) :: c_geovars ! variables requested from the model
character(len=maxvarlen), dimension(2), parameter :: geovars_default = (/var_qg, var_delp/)

type(ufo_lightning_tlad), pointer :: self
call ufo_lightning_tlad_registry%setup(c_key_self, self)

self%obsvars = obs_variables(c_obsvars)
self%geovars = oops_variables(c_geovars)
call self%geovars%push_back(geovars_default)

self%n_horiz = nhoriz_
!self%l_fed_nonlinear = l_fed_nonlinear_

end subroutine ufo_lightning_tlad_setup_c

! ------------------------------------------------------------------------------

subroutine ufo_lightning_tlad_delete_c(c_key_self) bind(c,name='ufo_lightning_tlad_delete_f90')
implicit none
integer(c_int), intent(inout) :: c_key_self

call ufo_lightning_tlad_registry%remove(c_key_self)

end subroutine ufo_lightning_tlad_delete_c

! ------------------------------------------------------------------------------

subroutine ufo_lightning_tlad_settraj_c(c_key_self, c_key_geovals, c_obsspace) &
                                        bind(c,name='ufo_lightning_tlad_settraj_f90')

implicit none
integer(c_int),     intent(in) :: c_key_self
integer(c_int),     intent(in) :: c_key_geovals
type(c_ptr), value, intent(in) :: c_obsspace

type(ufo_lightning_tlad), pointer :: self
type(ufo_geovals),      pointer :: geovals

call ufo_lightning_tlad_registry%get(c_key_self, self)
call ufo_geovals_registry%get(c_key_geovals, geovals)

call self%settraj(geovals, c_obsspace)

end subroutine ufo_lightning_tlad_settraj_c

! ------------------------------------------------------------------------------

subroutine ufo_lightning_simobs_tl_c(c_key_self, c_key_geovals, c_obsspace, c_nvars, c_nlocs, c_hofx) &
                                     bind(c,name='ufo_lightning_simobs_tl_f90')

use ufo_geovals_mod_c, only: ufo_geovals_registry
use ufo_geovals_mod,   only: ufo_geovals
implicit none
integer(c_int), intent(in) :: c_key_self
integer(c_int), intent(in) :: c_key_geovals
type(c_ptr), value, intent(in) :: c_obsspace
integer(c_int), intent(in) :: c_nvars, c_nlocs
real(c_double), intent(inout) :: c_hofx(c_nvars, c_nlocs)


type(ufo_lightning_tlad), pointer :: self
type(ufo_geovals),      pointer :: geovals

call ufo_lightning_tlad_registry%get(c_key_self, self)
call ufo_geovals_registry%get(c_key_geovals, geovals)
call self%simobs_tl(geovals, c_obsspace, c_nvars, c_nlocs, c_hofx)

end subroutine ufo_lightning_simobs_tl_c

! ------------------------------------------------------------------------------

subroutine ufo_lightning_simobs_ad_c(c_key_self, c_key_geovals, c_obsspace, c_nvars, c_nlocs, c_hofx) &
                                     bind(c,name='ufo_lightning_simobs_ad_f90')

use ufo_geovals_mod_c, only: ufo_geovals_registry
use ufo_geovals_mod,   only: ufo_geovals
implicit none
integer(c_int), intent(in) :: c_key_self
integer(c_int), intent(in) :: c_key_geovals
type(c_ptr), value, intent(in) :: c_obsspace
integer(c_int), intent(in) :: c_nvars, c_nlocs
real(c_double), intent(in) :: c_hofx(c_nvars, c_nlocs)

type(ufo_lightning_tlad), pointer :: self
type(ufo_geovals),      pointer :: geovals

call ufo_lightning_tlad_registry%get(c_key_self, self)
call ufo_geovals_registry%get(c_key_geovals, geovals)
call self%simobs_ad(geovals, c_obsspace, c_nvars, c_nlocs, c_hofx)

end subroutine ufo_lightning_simobs_ad_c

! ------------------------------------------------------------------------------

end module ufo_lightning_tlad_mod_c
