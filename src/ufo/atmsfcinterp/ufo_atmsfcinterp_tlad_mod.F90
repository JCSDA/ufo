! (C) Copyright 2017-2018 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

!> Fortran module for atmsfcinterp tl/ad observation operator

module ufo_atmsfcinterp_tlad_mod

 use iso_c_binding
 use config_mod
 use kinds

 use ufo_geovals_mod, only: ufo_geovals, ufo_geoval, ufo_geovals_get_var
 use ufo_basis_tlad_mod, only: ufo_basis_tlad
 use ufo_vars_mod
 use obsspace_mod

 implicit none
 private

 !> Fortran derived type for the tl/ad observation operator
 ! TODO: add to the below type what you need for your tl/ad observation operator
 !       this type can hold information on trajectory, for atmsfcinterp
 type, extends(ufo_basis_tlad), public :: ufo_atmsfcinterp_tlad
 private
   integer :: nvars_in
   character(len=MAXVARLEN), public, allocatable :: varin(:)
 contains
  procedure :: setup  => ufo_atmsfcinterp_tlad_setup
  procedure :: delete  => ufo_atmsfcinterp_tlad_delete
  procedure :: settraj => ufo_atmsfcinterp_tlad_settraj
  procedure :: simobs_tl  => ufo_atmsfcinterp_simobs_tl
  procedure :: simobs_ad  => ufo_atmsfcinterp_simobs_ad
 end type ufo_atmsfcinterp_tlad

contains

! ------------------------------------------------------------------------------
! TODO: add setup of your TL/AD observation operator (optional)
subroutine ufo_atmsfcinterp_tlad_setup(self, c_conf)
implicit none
class(ufo_atmsfcinterp_tlad), intent(inout) :: self
type(c_ptr),              intent(in)    :: c_conf

end subroutine ufo_atmsfcinterp_tlad_setup

! ------------------------------------------------------------------------------
! TODO: add cleanup of your TL/AD observation operator (optional)
subroutine ufo_atmsfcinterp_tlad_delete(self)
implicit none
class(ufo_atmsfcinterp_tlad), intent(inout) :: self

end subroutine ufo_atmsfcinterp_tlad_delete

! ------------------------------------------------------------------------------
! TODO: replace below function with your set trajectory for tl/ad code
subroutine ufo_atmsfcinterp_tlad_settraj(self, geovals, obss)
implicit none
class(ufo_atmsfcinterp_tlad), intent(inout) :: self
type(ufo_geovals),       intent(in)    :: geovals
type(c_ptr), value,      intent(in)    :: obss

end subroutine ufo_atmsfcinterp_tlad_settraj

! ------------------------------------------------------------------------------
! TODO: replace below function with your tl observation operator.
! Note: this can use information saved from trajectory in your ufo_atmsfcinterp_tlad type
! Input geovals parameter represents dx for tangent linear model
subroutine ufo_atmsfcinterp_simobs_tl(self, geovals, hofx, obss)
implicit none
class(ufo_atmsfcinterp_tlad), intent(in)    :: self
type(ufo_geovals),       intent(in)    :: geovals
real(c_double),          intent(inout) :: hofx(:)
type(c_ptr), value,      intent(in)    :: obss

end subroutine ufo_atmsfcinterp_simobs_tl

! ------------------------------------------------------------------------------
! TODO: replace below function with your ad observation operator.
! Note: this can use information saved from trajectory in your ufo_atmsfcinterp_tlad type
subroutine ufo_atmsfcinterp_simobs_ad(self, geovals, hofx, obss)
implicit none
class(ufo_atmsfcinterp_tlad), intent(in)    :: self
type(ufo_geovals),       intent(inout) :: geovals
real(c_double),          intent(in)    :: hofx(:)
type(c_ptr), value,      intent(in)    :: obss


end subroutine ufo_atmsfcinterp_simobs_ad

! ------------------------------------------------------------------------------

end module ufo_atmsfcinterp_tlad_mod
