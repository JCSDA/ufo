!
! (C) Copyright 2020 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

!>  Define interface for C++ ufo::SampledLocations code called from Fortran

!-------------------------------------------------------------------------------
interface
!-------------------------------------------------------------------------------
integer(kind=c_size_t) function c_sampled_locations_get_npaths(locs) &
              & bind(C,name='sampled_locations_get_npaths_f')
  use, intrinsic :: iso_c_binding
  implicit none

  type(c_ptr), value :: locs
end function c_sampled_locations_get_npaths

subroutine c_sampled_locations_get_lons(locs, npaths, lons) &
              & bind(C,name='sampled_locations_get_lons_f')
  use, intrinsic :: iso_c_binding, only : c_ptr,c_char,c_size_t,c_double
  implicit none
  type(c_ptr), value :: locs
  integer(c_size_t), intent(in) :: npaths
  real(c_double), intent(inout) :: lons(npaths)
end subroutine c_sampled_locations_get_lons

subroutine c_sampled_locations_get_lats(locs, npaths, lats) &
              & bind(C,name='sampled_locations_get_lats_f')
  use, intrinsic :: iso_c_binding, only : c_ptr,c_char,c_size_t,c_double
  implicit none
  type(c_ptr), value :: locs
  integer(c_size_t), intent(in) :: npaths
  real(c_double), intent(inout) :: lats(npaths)
end subroutine c_sampled_locations_get_lats

subroutine c_sampled_locations_get_timemask(locs, t1, t2, npaths, mask) &
              & bind(C,name='sampled_locations_get_timemask_f')
  use, intrinsic :: iso_c_binding, only : c_ptr,c_char,c_size_t,c_bool
  implicit none
  type(c_ptr), value :: locs, t1, t2
  integer(c_size_t), intent(in) :: npaths
  logical(c_bool), intent(inout) :: mask(npaths)
end subroutine c_sampled_locations_get_timemask

!-------------------------------------------------------------------------------
end interface
!-------------------------------------------------------------------------------
