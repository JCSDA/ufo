! (C) Copyright 2019 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

!> Fortran module for geosaod tl/ad observation operator

module ufo_geosaod_tlad_mod

 use iso_c_binding
 use config_mod
 use kinds

 use ufo_geovals_mod, only: ufo_geovals, ufo_geoval, ufo_geovals_get_var
 use ufo_basis_tlad_mod, only: ufo_basis_tlad
 use ufo_vars_mod
 use obsspace_mod

 implicit none
 private
 integer, parameter :: max_string=800

 !> Fortran derived type for the tl/ad observation operator
 type, extends(ufo_basis_tlad), public :: ufo_geosaod_tlad
 private
  integer :: nvars_in
  character(len=max_string), public, allocatable :: varin(:)
 contains
  procedure :: setup  => ufo_geosaod_tlad_setup
  procedure :: delete  => ufo_geosaod_tlad_delete
  procedure :: settraj => ufo_geosaod_tlad_settraj
  procedure :: simobs_tl  => ufo_geosaod_simobs_tl
  procedure :: simobs_ad  => ufo_geosaod_simobs_ad
 end type ufo_geosaod_tlad

contains

! ------------------------------------------------------------------------------

subroutine ufo_geosaod_tlad_setup(self, c_conf)

implicit none
class(ufo_geosaod_tlad), intent(inout) :: self
type(c_ptr),             intent(in)    :: c_conf


end subroutine ufo_geosaod_tlad_setup

! ------------------------------------------------------------------------------

subroutine ufo_geosaod_tlad_delete(self)

implicit none
class(ufo_geosaod_tlad), intent(inout) :: self

  if (allocated(self%varin))   deallocate(self%varin)

end subroutine ufo_geosaod_tlad_delete

! ------------------------------------------------------------------------------
subroutine ufo_geosaod_tlad_settraj(self, geovals, obss)

implicit none
class(ufo_geosaod_tlad), intent(inout) :: self
type(ufo_geovals),       intent(in)    :: geovals
type(c_ptr), value,      intent(in)    :: obss

end subroutine ufo_geosaod_tlad_settraj

! ------------------------------------------------------------------------------

subroutine ufo_geosaod_simobs_tl(self, geovals, hofx, obss)

implicit none
class(ufo_geosaod_tlad), intent(in)    :: self
type(ufo_geovals),       intent(in)    :: geovals
real(c_double),          intent(inout) :: hofx(:)
type(c_ptr), value,      intent(in)    :: obss

end subroutine ufo_geosaod_simobs_tl

! ------------------------------------------------------------------------------

subroutine ufo_geosaod_simobs_ad(self, geovals, hofx, obss)

implicit none
class(ufo_geosaod_tlad), intent(in)    :: self
type(ufo_geovals),       intent(inout) :: geovals
real(c_double),          intent(in)    :: hofx(:)
type(c_ptr), value,      intent(in)    :: obss


end subroutine ufo_geosaod_simobs_ad

! ------------------------------------------------------------------------------

end module ufo_geosaod_tlad_mod
