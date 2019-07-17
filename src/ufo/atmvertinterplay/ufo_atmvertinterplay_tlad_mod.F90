! (C) Copyright 2017-2018 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

!> Fortran module for atmvertinterplay tl/ad observation operator

module ufo_atmvertinterplay_tlad_mod

 use iso_c_binding
 use kinds

 use ufo_geovals_mod
 use ufo_geovals_mod_c, only: ufo_geovals_registry
 use ufo_vars_mod
 use obsspace_mod

 implicit none
 private
 integer, parameter :: max_string=800

 !> Fortran derived type for the tl/ad observation operator
 ! TODO: add to the below type what you need for your tl/ad observation operator
 !       this type can hold information on trajectory, for atmvertinterplay
 type, public :: ufo_atmvertinterplay_tlad
 private
  integer :: nvars
  integer :: nval, nlocs
  character(len=max_string), public, allocatable :: varin(:)
 contains
  procedure :: setup  => atmvertinterplay_tlad_setup_
  procedure :: delete  => atmvertinterplay_tlad_delete_
  procedure :: settraj => atmvertinterplay_tlad_settraj_
  procedure :: simobs_tl  => atmvertinterplay_simobs_tl_
  procedure :: simobs_ad  => atmvertinterplay_simobs_ad_
 end type ufo_atmvertinterplay_tlad

contains

! ------------------------------------------------------------------------------
subroutine atmvertinterplay_tlad_setup_(self, vars)
implicit none
class(ufo_atmvertinterplay_tlad), intent(inout) :: self
character(len=MAXVARLEN), dimension(:), intent(inout) :: vars

  self%nvars = size(vars)
  !> Allocate varin
  allocate(self%varin(self%nvars))
  self%varin = vars

end subroutine atmvertinterplay_tlad_setup_

! ------------------------------------------------------------------------------
subroutine atmvertinterplay_tlad_delete_(self)
implicit none
class(ufo_atmvertinterplay_tlad), intent(inout) :: self

  if (allocated(self%varin))   deallocate(self%varin)

end subroutine atmvertinterplay_tlad_delete_

! ------------------------------------------------------------------------------
subroutine atmvertinterplay_tlad_settraj_(self, geovals, obss)
implicit none
class(ufo_atmvertinterplay_tlad), intent(inout) :: self
type(ufo_geovals),       intent(in)    :: geovals
type(c_ptr), value,      intent(in)    :: obss

end subroutine atmvertinterplay_tlad_settraj_

! ------------------------------------------------------------------------------
! Note: this can use information saved from trajectory in your ufo_atmvertinterplay_tlad type
! Input geovals parameter represents dx for tangent linear model
subroutine atmvertinterplay_simobs_tl_(self, geovals, obss, nvars, nlocs, hofx)
implicit none
class(ufo_atmvertinterplay_tlad), intent(in)    :: self
type(ufo_geovals),       intent(inout)    :: geovals
integer, intent(in) :: nvars, nlocs
real(c_double),          intent(in) :: hofx(nvars, nlocs)
type(c_ptr), value,      intent(in)    :: obss

end subroutine atmvertinterplay_simobs_tl_

! ------------------------------------------------------------------------------
! Note: this can use information saved from trajectory in your ufo_atmvertinterplay_tlad type
subroutine atmvertinterplay_simobs_ad_(self, geovals, obss, nvars, nlocs, hofx)
implicit none
class(ufo_atmvertinterplay_tlad), intent(in)    :: self
type(ufo_geovals),       intent(inout) :: geovals
integer, intent(in) :: nvars, nlocs
real(c_double),          intent(in)    :: hofx(nvars, nlocs)
type(c_ptr), value,      intent(in)    :: obss


end subroutine atmvertinterplay_simobs_ad_

! ------------------------------------------------------------------------------

end module ufo_atmvertinterplay_tlad_mod
