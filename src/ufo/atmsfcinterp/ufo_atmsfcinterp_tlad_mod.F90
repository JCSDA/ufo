! (C) Copyright 2017-2018 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

!> Fortran module for atmsfcinterp tl/ad observation operator

module ufo_atmsfcinterp_tlad_mod

 use iso_c_binding
 use kinds
 use ufo_vars_mod
 use ufo_geovals_mod
 use ufo_geovals_mod_c, only: ufo_geovals_registry
 use obsspace_mod
 use missing_values_mod

 implicit none
 private

 !> Fortran derived type for the tl/ad observation operator
 ! TODO: add to the below type what you need for your tl/ad observation operator
 !       this type can hold information on trajectory, for atmsfcinterp
 type, public :: ufo_atmsfcinterp_tlad
 private
   integer :: nvars
   character(len=MAXVARLEN), public, allocatable :: varin(:)
   integer :: nval, nlocs
 contains
  procedure :: setup  => ufo_atmsfcinterp_tlad_setup_
  procedure :: cleanup  => ufo_atmsfcinterp_tlad_cleanup_
  procedure :: settraj => ufo_atmsfcinterp_tlad_settraj_
  procedure :: simobs_tl  => ufo_atmsfcinterp_simobs_tl_
  procedure :: simobs_ad  => ufo_atmsfcinterp_simobs_ad_
  final :: destructor
 end type ufo_atmsfcinterp_tlad

contains

! ------------------------------------------------------------------------------
! TODO: add setup of your TL/AD observation operator (optional)
subroutine ufo_atmsfcinterp_tlad_setup_(self, vars)
  implicit none
  class(ufo_atmsfcinterp_tlad), intent(inout) :: self
  character(len=MAXVARLEN), dimension(:), intent(inout) :: vars

end subroutine ufo_atmsfcinterp_tlad_setup_

! ------------------------------------------------------------------------------
! TODO: add cleanup of your TL/AD observation operator (optional)
subroutine ufo_atmsfcinterp_tlad_cleanup_(self)
  implicit none
  class(ufo_atmsfcinterp_tlad), intent(inout) :: self

end subroutine ufo_atmsfcinterp_tlad_cleanup_

! ------------------------------------------------------------------------------
! TODO: replace below function with your set trajectory for tl/ad code
subroutine ufo_atmsfcinterp_tlad_settraj_(self, geovals, obss)
  implicit none
  class(ufo_atmsfcinterp_tlad), intent(inout) :: self
  type(ufo_geovals),       intent(in)    :: geovals
  type(c_ptr), value,      intent(in)    :: obss

end subroutine ufo_atmsfcinterp_tlad_settraj_

! ------------------------------------------------------------------------------
! TODO: replace below function with your tl observation operator.
! Note: this can use information saved from trajectory in your ufo_atmsfcinterp_tlad type
! Input geovals parameter represents dx for tangent linear model
subroutine ufo_atmsfcinterp_simobs_tl_(self, geovals, obss, nvars, nlocs, hofx)
  implicit none
  class(ufo_atmsfcinterp_tlad), intent(in)    :: self
  type(ufo_geovals),       intent(in)    :: geovals
  integer,                 intent(in)    :: nvars, nlocs
  real(c_double),          intent(inout) :: hofx(nvars, nlocs)
  type(c_ptr), value,      intent(in)    :: obss

end subroutine ufo_atmsfcinterp_simobs_tl_

! ------------------------------------------------------------------------------
! TODO: replace below function with your ad observation operator.
! Note: this can use information saved from trajectory in your ufo_atmsfcinterp_tlad type
subroutine ufo_atmsfcinterp_simobs_ad_(self, geovals, obss, nvars, nlocs, hofx)
  implicit none
  class(ufo_atmsfcinterp_tlad), intent(in)    :: self
  type(ufo_geovals),       intent(inout) :: geovals
  integer,                 intent(in)    :: nvars, nlocs
  real(c_double),          intent(in)    :: hofx(nvars, nlocs)
  type(c_ptr), value,      intent(in)    :: obss


end subroutine ufo_atmsfcinterp_simobs_ad_

! ------------------------------------------------------------------------------

subroutine  destructor(self)
  type(ufo_atmsfcinterp_tlad), intent(inout)  :: self

  call self%cleanup()
  self%nvars = 0
  if (allocated(self%varin)) deallocate(self%varin)

end subroutine destructor


! ------------------------------------------------------------------------------

end module ufo_atmsfcinterp_tlad_mod
