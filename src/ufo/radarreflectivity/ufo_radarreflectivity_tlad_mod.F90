! (C) Copyright 2017-2018 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

!> Fortran module for radarreflectivity tl/ad observation operator

module ufo_radarreflectivity_tlad_mod

 use fckit_configuration_module, only: fckit_configuration
 use iso_c_binding
 use kinds

 use ufo_geovals_mod, only: ufo_geovals, ufo_geoval, ufo_geovals_get_var
 use ufo_vars_mod
 use obsspace_mod

 implicit none
 private
 integer, parameter :: max_string=800

 !> Fortran derived type for the tl/ad observation operator
 ! TODO: add to the below type what you need for your tl/ad observation operator
 !       this type can hold information on trajectory, for radarreflectivity
 type, public :: ufo_radarreflectivity_tlad
 private
  integer :: nvars_in
  character(len=max_string), public, allocatable :: varin(:)
 contains
  procedure :: setup  => ufo_radarreflectivity_tlad_setup
  procedure :: settraj => ufo_radarreflectivity_tlad_settraj
  procedure :: simobs_tl  => ufo_radarreflectivity_simobs_tl
  procedure :: simobs_ad  => ufo_radarreflectivity_simobs_ad
  final :: destructor
 end type ufo_radarreflectivity_tlad

contains

! ------------------------------------------------------------------------------
! TODO: add setup of your TL/AD observation operator (optional)
subroutine ufo_radarreflectivity_tlad_setup(self, f_conf, vars)
implicit none
class(ufo_radarreflectivity_tlad), intent(inout) :: self
type(fckit_configuration), intent(in)  :: f_conf
character(len=MAXVARLEN), dimension(:), intent(inout) :: vars ! variables to be simulated

! TODO: setup input variables varin (updated model variables)
  self%nvars_in = 0

end subroutine ufo_radarreflectivity_tlad_setup

! ------------------------------------------------------------------------------
! TODO: add cleanup of your TL/AD observation operator (optional)
subroutine destructor(self)
implicit none
type(ufo_radarreflectivity_tlad), intent(inout) :: self

  if (allocated(self%varin))   deallocate(self%varin)

end subroutine destructor

! ------------------------------------------------------------------------------
! TODO: replace below function with your set trajectory for tl/ad code
subroutine ufo_radarreflectivity_tlad_settraj(self, geovals, obss)
implicit none
class(ufo_radarreflectivity_tlad), intent(inout) :: self
type(ufo_geovals),       intent(in)    :: geovals
type(c_ptr), value,      intent(in)    :: obss

end subroutine ufo_radarreflectivity_tlad_settraj

! ------------------------------------------------------------------------------
! TODO: replace below function with your tl observation operator.
! Note: this can use information saved from trajectory in your ufo_radarreflectivity_tlad type
! Input geovals parameter represents dx for tangent linear model
subroutine ufo_radarreflectivity_simobs_tl(self, geovals, obss, nvars, nlocs, hofx)
implicit none
class(ufo_radarreflectivity_tlad), intent(in)    :: self
type(ufo_geovals),       intent(in)    :: geovals
integer,                 intent(in)    :: nvars, nlocs
real(c_double),          intent(inout) :: hofx(nvars, nlocs)
type(c_ptr), value,      intent(in)    :: obss

end subroutine ufo_radarreflectivity_simobs_tl

! ------------------------------------------------------------------------------
! TODO: replace below function with your ad observation operator.
! Note: this can use information saved from trajectory in your ufo_radarreflectivity_tlad type
subroutine ufo_radarreflectivity_simobs_ad(self, geovals, obss, nvars, nlocs, hofx)
implicit none
class(ufo_radarreflectivity_tlad), intent(in)    :: self
type(ufo_geovals),       intent(inout) :: geovals
integer,                 intent(in)    :: nvars, nlocs
real(c_double),          intent(in)    :: hofx(nvars, nlocs)
type(c_ptr), value,      intent(in)    :: obss


end subroutine ufo_radarreflectivity_simobs_ad

! ------------------------------------------------------------------------------

end module ufo_radarreflectivity_tlad_mod
