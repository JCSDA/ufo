!
! (C) Copyright 2017 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

!> Fortran module handling observation locations

module ufo_locs_mod

use datetime_mod
use iso_c_binding
use kinds
use fckit_log_module, only : fckit_log
use obsspace_mod

implicit none
private
public :: ufo_locs, ufo_locs_create, ufo_locs_copy, ufo_locs_setup, ufo_locs_delete
public :: ufo_locs_init, ufo_locs_concatenate, ufo_locs_time_mask

! --------------------------------------------------------------------------------------------------

!> Fortran derived type to hold observation locations
type :: ufo_locs
  integer :: nlocs
  integer :: max_indx
  real(kind_real), allocatable, dimension(:) :: lat     !< latitude
  real(kind_real), allocatable, dimension(:) :: lon     !< longitude
  type(datetime),  allocatable, dimension(:) :: time    !< obs-time
  integer,         allocatable, dimension(:) :: indx    !< indices of locations in the full [geovals] array
end type ufo_locs

! --------------------------------------------------------------------------------------------------
contains
! --------------------------------------------------------------------------------------------------

subroutine ufo_locs_create(self, obss, nlocs, lats, lons)
implicit none
type(ufo_locs),     intent(inout) :: self
type(c_ptr), value, intent(in)    :: obss
integer,            intent(in)    :: nlocs
real(kind_real),    intent(in)    :: lats(nlocs)
real(kind_real),    intent(in)    :: lons(nlocs)

integer :: n
character(len=20), allocatable :: fstring(:)
type(datetime), dimension(:), allocatable :: date_time
self%nlocs = nlocs
self%max_indx = nlocs
allocate(self%lat(nlocs), self%lon(nlocs), self%time(nlocs))
allocate(self%indx(nlocs))
self%lat(:) = lats(:)
self%lon(:) = lons(:)

allocate(fstring(nlocs))

if (obsspace_has(obss,"MetaData", "datetime")) then
  allocate(date_time(nlocs))
  call obsspace_get_db(obss, "MetaData", "datetime", date_time)
  do n = 1, self%nlocs
    call datetime_to_string(date_time(n), fstring(n))
  enddo
  deallocate(date_time)
else
  fstring(:) = "9999-09-09T09:09:09Z"
endif

do n = 1, self%nlocs
  call datetime_create(fstring(n), self%time(n))
enddo
do n = 1, self%nlocs
  self%indx(n) = n
enddo

end subroutine ufo_locs_create

! --------------------------------------------------------------------------------------------------

subroutine ufo_locs_setup(self, nlocs)
implicit none
type(ufo_locs), intent(inout) :: self
integer, intent(in)           :: nlocs

character(len=20) :: fstring
integer :: n

call ufo_locs_delete(self)

self%max_indx = nlocs
self%nlocs = nlocs
allocate(self%lat(nlocs), self%lon(nlocs), self%time(nlocs), self%indx(nlocs))
self%lat(:) = 0.0
self%lon(:) = 0.0
fstring="9999-09-09T09:09:09Z"
do n = 1, self%nlocs
  call datetime_create(fstring, self%time(n))
enddo
self%indx(:) = 0

end subroutine ufo_locs_setup

! --------------------------------------------------------------------------------------------------

subroutine ufo_locs_copy(self, other)

implicit none
type(ufo_locs), intent(inout) :: self
type(ufo_locs), intent(in)    :: other

self%nlocs = other%nlocs
self%max_indx = other%max_indx

allocate(self%lat (self%nlocs))
allocate(self%lon (self%nlocs))
allocate(self%time(self%nlocs))
allocate(self%indx(self%nlocs))

self%lat  = other%lat
self%lon  = other%lon
self%time = other%time
self%indx = other%indx

end subroutine ufo_locs_copy

! --------------------------------------------------------------------------------------------------

subroutine ufo_locs_delete(self)
implicit none
type(ufo_locs), intent(inout) :: self

if (allocated(self%lat)) deallocate(self%lat)
if (allocated(self%lon)) deallocate(self%lon)
if (allocated(self%time)) deallocate(self%time)
if (allocated(self%indx)) deallocate(self%indx)
self%nlocs = 0
self%max_indx = -1 ! not set

end subroutine ufo_locs_delete

! --------------------------------------------------------------------------------------------------

subroutine ufo_locs_concatenate(self, other)
implicit none
type(ufo_locs), intent(inout) :: self
type(ufo_locs), intent(in) :: other

type(ufo_locs) :: temp_self, temp_other
character(255) :: message
integer :: n

if ((self%max_indx < 0) .OR. (other%max_indx < 0)) then
  write(message,'(A, A, I6, A, I6)') &
    'ufo_locs_concatenate: either self or other needs to be constructed valid indices ', &
    ' self%max_indx =', self%max_indx, ' other%max_indx = ', other%max_indx
  call fckit_log%info(message)
  stop
end if

! make a temporary copy of self
temp_self%nlocs = self%nlocs
temp_self%max_indx = self%max_indx
allocate(temp_self%lat(self%nlocs), temp_self%lon(self%nlocs), &
         temp_self%time(self%nlocs), temp_self%indx(self%nlocs))

temp_self%lat(:) = self%lat(:)
temp_self%lon(:) = self%lon(:)
temp_self%indx(:) = self%indx(:)
do n = 1, self%nlocs
  temp_self%time(n) = self%time(n)
end do

! make a temporary copy of other
temp_other%nlocs = other%nlocs
temp_other%max_indx = other%max_indx
allocate(temp_other%lat(other%nlocs), temp_other%lon(other%nlocs), &
         temp_other%time(other%nlocs), temp_other%indx(other%nlocs))

temp_other%lat(:) = other%lat(:)
temp_other%lon(:) = other%lon(:)
temp_other%indx(:) = other%indx(:)
do n = 1, other%nlocs
  temp_other%time(n) = other%time(n)
end do

! deallocate self
call ufo_locs_delete(self)

! reallocate self with combined concatenation
self%nlocs = temp_self%nlocs + temp_other%nlocs
allocate(self%lat(self%nlocs), self%lon(self%nlocs), &
         self%time(self%nlocs), self%indx(self%nlocs))

self%lat(1:temp_self%nlocs) = temp_self%lat(1:temp_self%nlocs)
self%lat(temp_self%nlocs+1:) = temp_other%lat(1:temp_other%nlocs)

self%lon(1:temp_self%nlocs) = temp_self%lon(1:temp_self%nlocs)
self%lon(temp_self%nlocs+1:) = temp_other%lon(1:temp_other%nlocs)

do n = 1, temp_self%nlocs
  self%time(n) = temp_self%time(n)
end do
do n = 1, temp_other%nlocs
  self%time(temp_self%nlocs + n) = temp_other%time(n)
end do

self%indx(1:temp_self%nlocs) = temp_self%indx(1:temp_self%nlocs)
do n = 1, temp_other%nlocs
  self%indx(temp_self%nlocs + n) = temp_other%indx(n) + &
                                  temp_self%max_indx
end do
self%max_indx = temp_self%max_indx + temp_other%max_indx

call ufo_locs_delete(temp_other)
call ufo_locs_delete(temp_self)

end subroutine ufo_locs_concatenate


! --------------------------------------------------------------------------------------------------

subroutine ufo_locs_init(self, obss)

  implicit none

  type(ufo_locs), intent(inout)   :: self
  type(c_ptr), value, intent(in)  :: obss

  integer :: nlocs, i

  nlocs = obsspace_get_nlocs(obss)
  call ufo_locs_setup(self, nlocs)

  call obsspace_get_db(obss, "MetaData", "datetime", self%time)
  call obsspace_get_db(obss, "MetaData", "longitude", self%lon)
  call obsspace_get_db(obss, "MetaData", "latitude", self%lat)

  !> index is left for back-compatibility; should be removed in the next
  !refactoring
  do i = 1, nlocs
    self%indx(i) = i
  enddo

  self%max_indx = obsspace_get_gnlocs(obss)

end subroutine ufo_locs_init

! --------------------------------------------------------------------------------------------------

subroutine ufo_locs_time_mask(self, t1, t2, time_mask)

type(ufo_locs),       intent(in)    :: self
type(datetime),       intent(in)    :: t1
type(datetime),       intent(in)    :: t2
logical, allocatable, intent(inout) :: time_mask(:)

! Locals
integer :: n

! Return a mask that is true where the location times are between t1 and t2

! Check for sensible inputs
if (t1>t2) call abor1_ftn("ufo_locs_mod.ufo_locs_time_mask t2 is not greater than or equal to t1")

! Allocate the array to output
if (.not.allocated(time_mask)) allocate(time_mask(self%nlocs))
time_mask = .false.

! Loop over times and check if between two times
do n = 1, self%nlocs

  ! Check if in the time range
  if (self%time(n) > t1 .and. self%time(n) <= t2 ) then
    time_mask(n) = .true.
  endif

enddo

end subroutine ufo_locs_time_mask

! --------------------------------------------------------------------------------------------------

end module ufo_locs_mod
