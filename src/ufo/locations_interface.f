!
! (C) Copyright 2020 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

!>  Define interface for C++ ufo::Locations code called from Fortran

!-------------------------------------------------------------------------------
interface
!-------------------------------------------------------------------------------
integer(kind=c_size_t) function c_locations_get_nlocs(locs) bind(C,name='locations_get_nlocs_f')
  use, intrinsic :: iso_c_binding
  implicit none

  type(c_ptr), value :: locs
end function c_locations_get_nlocs

subroutine c_locations_get_lons(locs, nlocs, lons) &
              & bind(C,name='locations_get_lons_f')
  use, intrinsic :: iso_c_binding, only : c_ptr,c_char,c_size_t,c_double
  implicit none
  type(c_ptr), value :: locs
  integer(c_size_t), intent(in) :: nlocs
  real(c_double), intent(inout) :: lons(nlocs)
end subroutine c_locations_get_lons

subroutine c_locations_get_lats(locs, nlocs, lats) &
              & bind(C,name='locations_get_lats_f')
  use, intrinsic :: iso_c_binding, only : c_ptr,c_char,c_size_t,c_double
  implicit none
  type(c_ptr), value :: locs
  integer(c_size_t), intent(in) :: nlocs
  real(c_double), intent(inout) :: lats(nlocs)
end subroutine c_locations_get_lats

subroutine c_locations_get_timemask(locs, t1, t2, nlocs, mask) &
              & bind(C,name='locations_get_timemask_f')
  use, intrinsic :: iso_c_binding, only : c_ptr,c_char,c_size_t,c_bool
  implicit none
  type(c_ptr), value :: locs, t1, t2
  integer(c_size_t), intent(in) :: nlocs
  logical(c_bool), intent(inout) :: mask(nlocs)
end subroutine c_locations_get_timemask

!-------------------------------------------------------------------------------
end interface
!-------------------------------------------------------------------------------
