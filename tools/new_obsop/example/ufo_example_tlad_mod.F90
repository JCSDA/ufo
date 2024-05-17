! (C) Copyright 2017-2018 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

!> Fortran module for example tl/ad observation operator

module ufo_example_tlad_mod

 use oops_variables_mod
 use obs_variables_mod
 use ufo_vars_mod

 implicit none
 private

 !> Fortran derived type for the tl/ad observation operator
 ! TODO: add to the below type what you need for your tl/ad observation operator
 !       this type can hold information on trajectory, for example
 type, public :: ufo_example_tlad
 private
  type(obs_variables), public :: obsvars
  type(oops_variables), public :: geovars
 contains
  procedure :: setup  => ufo_example_tlad_setup
  procedure :: settraj => ufo_example_tlad_settraj
  procedure :: simobs_tl  => ufo_example_simobs_tl
  procedure :: simobs_ad  => ufo_example_simobs_ad
  final :: destructor
 end type ufo_example_tlad

contains

! ------------------------------------------------------------------------------
! TODO: add setup of your TL/AD observation operator (optional)
subroutine ufo_example_tlad_setup(self, f_conf)
use fckit_configuration_module, only: fckit_configuration
implicit none
class(ufo_example_tlad), intent(inout) :: self

! TODO: consider whether passing the Configuration object to this function
! is necessary. If only a small number of parameters are used,
! you could pass them in directly instead. In that case you can modify the
! interface appropriately.
type(fckit_configuration), intent(in)  :: f_conf

! TODO: setup input variables varin (updated model variables)
!  self%geovars%push_back("variable name")

end subroutine ufo_example_tlad_setup

! ------------------------------------------------------------------------------
! TODO: add cleanup of your TL/AD observation operator (optional)
subroutine destructor(self)
implicit none
type(ufo_example_tlad), intent(inout) :: self

end subroutine destructor

! ------------------------------------------------------------------------------
! TODO: replace below function with your set trajectory for tl/ad code
subroutine ufo_example_tlad_settraj(self, geovals, obss, hofxdiags)
use iso_c_binding
use ufo_geovals_mod, only: ufo_geovals, ufo_geoval, ufo_geovals_get_var
use obsspace_mod
implicit none
class(ufo_example_tlad), intent(inout) :: self
type(ufo_geovals),       intent(in)    :: geovals
type(c_ptr), value,      intent(in)    :: obss
type(ufo_geovals),       intent(inout) :: hofxdiags    !non-h(x) diagnostics

end subroutine ufo_example_tlad_settraj

! ------------------------------------------------------------------------------
! TODO: replace below function with your tl observation operator.
! Note: this can use information saved from trajectory in your ufo_example_tlad type
! Input geovals parameter represents dx for tangent linear model
subroutine ufo_example_simobs_tl(self, geovals, obss, nvars, nlocs, hofx)
use iso_c_binding
use ufo_geovals_mod, only: ufo_geovals, ufo_geoval, ufo_geovals_get_var
use obsspace_mod
implicit none
class(ufo_example_tlad), intent(in)    :: self
type(ufo_geovals),       intent(in)    :: geovals
integer,                 intent(in)    :: nvars, nlocs
real(c_double),          intent(inout) :: hofx(nvars, nlocs)
type(c_ptr), value,      intent(in)    :: obss

end subroutine ufo_example_simobs_tl

! ------------------------------------------------------------------------------
! TODO: replace below function with your ad observation operator.
! Note: this can use information saved from trajectory in your ufo_example_tlad type
subroutine ufo_example_simobs_ad(self, geovals, obss, nvars, nlocs, hofx)
use iso_c_binding
use ufo_geovals_mod, only: ufo_geovals, ufo_geoval, ufo_geovals_get_var
use obsspace_mod
implicit none
class(ufo_example_tlad), intent(in)    :: self
type(ufo_geovals),       intent(inout) :: geovals
integer,                 intent(in)    :: nvars, nlocs
real(c_double),          intent(in)    :: hofx(nvars, nlocs)
type(c_ptr), value,      intent(in)    :: obss


end subroutine ufo_example_simobs_ad

! ------------------------------------------------------------------------------

end module ufo_example_tlad_mod
