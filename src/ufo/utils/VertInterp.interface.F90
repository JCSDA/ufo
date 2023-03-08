! (C) Copyright 2017-2019 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

!> Fortran sfcpcorrected module for functions on the interface between C++ and Fortran
!  to handle observation operators

module vert_interp_mod_c

  use iso_c_binding
  use vert_interp_mod
  implicit none
  private

contains

! ------------------------------------------------------------------------------

subroutine vert_interp_weights_c(c_nlev, c_obl, c_vec, c_wi, c_wf) &
  bind(c,name='vert_interp_weights_f90')

implicit none
integer(c_int), intent(in ) :: c_nlev         !Number of model levels
real(c_double), intent(in ) :: c_obl          !Observation location
real(c_double), intent(in ) :: c_vec(c_nlev)  !Structured vector of grid points
integer(c_int), intent(out) :: c_wi           !Index for interpolation
real(c_double), intent(out) :: c_wf           !Weight for interpolation

call vert_interp_weights(c_nlev, c_obl, c_vec, c_wi, c_wf)

end subroutine vert_interp_weights_c

! ------------------------------------------------------------------------------

subroutine vert_interp_apply_c(c_nlev, c_fvec, c_f, c_wi, c_wf) &
  bind(c,name='vert_interp_apply_f90')

implicit none
integer(c_int), intent(in ) :: c_nlev          !Number of model levels
real(c_double), intent(in ) :: c_fvec(c_nlev)  !Field at grid points
real(c_double), intent(out) :: c_f             !Output at obs location using linear interp
integer(c_int), intent(in ) :: c_wi            !Index for interpolation
real(c_double), intent(in ) :: c_wf            !Weight for interpolation

call vert_interp_apply(c_nlev, c_fvec, c_f, c_wi, c_wf)

end subroutine vert_interp_apply_c

! ------------------------------------------------------------------------------

subroutine nearestneighbor_interp_index_c(c_nlev, c_obl, c_vec, c_idx) &
  bind(c,name='nearestneighbor_interp_index_f90')

implicit none
integer(c_int), intent(in ) :: c_nlev         !Number of model levels
real(c_double), intent(in ) :: c_obl          !Observation location
real(c_double), intent(in ) :: c_vec(c_nlev)  !Structured vector of grid points
integer(c_int), intent(out) :: c_idx          !Index for interpolation

call nearestneighbor_interp_index(c_nlev, c_obl, c_vec, c_idx)

end subroutine nearestneighbor_interp_index_c

! ------------------------------------------------------------------------------

subroutine nearestneighbor_interp_apply_c(c_nlev, c_fvec, c_f, c_idx) &
  bind(c,name='nearestneighbor_interp_apply_f90')

implicit none
integer(c_int), intent(in ) :: c_nlev          !Number of model levels
real(c_double), intent(in ) :: c_fvec(c_nlev)  !Field at grid points
real(c_double), intent(out) :: c_f             !Output at obs location using linear interp
integer(c_int), intent(in ) :: c_idx           !Index for interpolation

call nearestneighbor_interp_apply(c_nlev, c_fvec, c_f, c_idx)

end subroutine nearestneighbor_interp_apply_c

! ------------------------------------------------------------------------------

end module vert_interp_mod_c
