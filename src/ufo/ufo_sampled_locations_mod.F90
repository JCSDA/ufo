!
! (C) Copyright 2020 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

!>  Fortran interface to ufo::SampledLocations

module ufo_sampled_locations_mod

use iso_c_binding
implicit none

public :: ufo_sampled_locations

type ufo_sampled_locations
private
  type(c_ptr) :: ptr
contains
  procedure, public :: npaths
  procedure, public :: get_lons
  procedure, public :: get_lats
  procedure, public :: get_timemask
end type

interface ufo_sampled_locations
  module procedure ctor_from_ptr
end interface

private

#include "ufo/sampled_locations_interface.f"

contains

!-------------------------------------------------------------------------------
!>  Saves pointer to the C++ SampledLocations object
function ctor_from_ptr(ptr) result(this)
  type(ufo_sampled_locations)     :: this
  type(c_ptr), intent(in) :: ptr

  this%ptr = ptr
end function ctor_from_ptr

!-------------------------------------------------------------------------------
!>  Return the number of paths in this collection
integer function npaths(this)
  implicit none
  class(ufo_sampled_locations), intent(in) :: this

  ! Implicit conversion from c_size_t to integer which is safe in this case
  npaths = c_sampled_locations_get_npaths(this%ptr)
end function npaths

!-------------------------------------------------------------------------------

!> Get longitudes from the SampledLocations object
subroutine get_lons(this, lons)
  implicit none
  class(ufo_sampled_locations), intent(in) :: this
  real(c_double), intent(inout)    :: lons(:)

  integer(c_size_t) :: length

  length = size(lons)
  call c_sampled_locations_get_lons(this%ptr, length, lons)

end subroutine get_lons

!-------------------------------------------------------------------------------

!> Get latitudes from the SampledLocations object
subroutine get_lats(this, lats)
  implicit none
  class(ufo_sampled_locations), intent(in) :: this
  real(c_double), intent(inout)    :: lats(:)

  integer(c_size_t) :: length

  length = size(lats)
  call c_sampled_locations_get_lats(this%ptr, length, lats)

end subroutine get_lats

!-------------------------------------------------------------------------------

!> Get time mask (paths that are between t1 & t2) from the SampledLocations object
subroutine get_timemask(this, t1, t2, mask)
  use datetime_mod
  implicit none
  class(ufo_sampled_locations), intent(in) :: this
  type(datetime), intent(in)       :: t1, t2
  logical(c_bool), intent(inout)   :: mask(:)

  integer(c_size_t) :: length
  type(c_ptr) :: c_t1, c_t2

  length = size(mask)
  call f_c_datetime(t1, c_t1)
  call f_c_datetime(t2, c_t2)
  call c_sampled_locations_get_timemask(this%ptr, c_t1, c_t2, length, mask)

end subroutine get_timemask

!-------------------------------------------------------------------------------

end module ufo_sampled_locations_mod
