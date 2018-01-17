!
! (C) Copyright 2017 UCAR
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 

!> Fortran module handling observation locations

module ufo_locs_mod

use ufo_obs_vectors
use kinds

implicit none
private
public :: ufo_locs, ufo_locs_create, ufo_locs_setup, ufo_locs_delete

! ------------------------------------------------------------------------------

!> Fortran derived type to hold observation locations
type :: ufo_locs
  integer :: nlocs
  real(kind_real), allocatable, dimension(:) :: lat     !< latitude
  real(kind_real), allocatable, dimension(:) :: lon     !< longitude
  real(kind_real), allocatable, dimension(:) :: time    !< obs-time
end type ufo_locs

! ------------------------------------------------------------------------------
contains
! ------------------------------------------------------------------------------

subroutine ufo_locs_create(self, nlocs, lats, lons)
implicit none
type(ufo_locs), intent(inout) :: self
integer, intent(in)           :: nlocs
real(kind_real), intent(in) :: lats(nlocs)
real(kind_real), intent(in) :: lons(nlocs)

self%nlocs = nlocs
allocate(self%lat(nlocs), self%lon(nlocs), self%time(nlocs))
self%lat(:) = lats(:)
self%lon(:) = lons(:)
self%time(:) = 0.0

end subroutine ufo_locs_create

! ------------------------------------------------------------------------------

subroutine ufo_locs_setup(self, nlocs)
implicit none
type(ufo_locs), intent(inout) :: self
integer, intent(in)           :: nlocs

call ufo_locs_delete(self)

self%nlocs = nlocs
allocate(self%lat(nlocs), self%lon(nlocs), self%time(nlocs))
self%lat(:) = 0.0
self%lon(:) = 0.0
self%time(:) = 0.0

end subroutine ufo_locs_setup

! ------------------------------------------------------------------------------

subroutine ufo_locs_delete(self)
implicit none
type(ufo_locs), intent(inout) :: self

self%nlocs = 0
if (allocated(self%lat)) deallocate(self%lat)
if (allocated(self%lon)) deallocate(self%lon)
if (allocated(self%time)) deallocate(self%time)

end subroutine ufo_locs_delete

! ------------------------------------------------------------------------------

end module ufo_locs_mod
