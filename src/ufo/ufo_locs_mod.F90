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
use type_distribution, only: random_distribution

implicit none
private
public :: ufo_locs, ufo_locs_create, ufo_locs_setup, ufo_locs_delete
public :: ufo_locs_init

! ------------------------------------------------------------------------------

!> Fortran derived type to hold observation locations
type :: ufo_locs
  integer :: nlocs
  real(kind_real), allocatable, dimension(:) :: lat     !< latitude
  real(kind_real), allocatable, dimension(:) :: lon     !< longitude
  type(datetime),  allocatable, dimension(:) :: time    !< obs-time
  integer,         allocatable, dimension(:) :: indx    !< indices of locations in the full [geovals] array
end type ufo_locs

! ------------------------------------------------------------------------------
contains
! ------------------------------------------------------------------------------

subroutine ufo_locs_create(self, nlocs, lats, lons, rdist)
implicit none
type(ufo_locs), intent(inout) :: self
integer, intent(in)           :: nlocs
real(kind_real), intent(in) :: lats(nlocs)
real(kind_real), intent(in) :: lons(nlocs)
integer, intent(in)         :: rdist

type(random_distribution) :: ran_dist
integer, allocatable :: dist_indx(:)
real(kind_real), allocatable, dimension(:) :: latsr, lonsr
type(datetime), allocatable, dimension(:) :: timer
integer :: n
character(len=20) :: fstring

self%nlocs = nlocs
allocate(self%lat(nlocs), self%lon(nlocs), self%time(nlocs))
allocate(self%indx(nlocs))
self%lat(:) = lats(:)
self%lon(:) = lons(:)
fstring="9999-09-09T09:09:09Z"
do n = 1, self%nlocs
  call datetime_create(fstring, self%time(n))
enddo
do n = 1, self%nlocs
  self%indx(n) = n
enddo

ran_dist = random_distribution(self%nlocs)
dist_indx = ran_dist%indx

if (rdist == 1) then

  !Redistribute randomly
  self%nlocs = ran_dist%nobs_pe()
  allocate(latsr(self%nlocs), lonsr(self%nlocs), timer(self%nlocs))
  do n = 1,self%nlocs
    latsr(n) =  self%lat(dist_indx(n))
    lonsr(n) =  self%lon(dist_indx(n))
    timer(n) = self%time(dist_indx(n))
  enddo
  do n = 1,self%nlocs
    call datetime_delete(self%time(n))
  enddo
  deallocate(self%lat, self%lon, self%time, self%indx)
  
  allocate(self%lat(self%nlocs), self%lon(self%nlocs), self%time(self%nlocs), self%indx(self%nlocs))
  self%lat(:) = latsr(:)
  self%lon(:) = lonsr(:)
  self%time(:) = timer(:)
  do n = 1, self%nlocs
    self%indx(n) = n
  enddo

  do n = 1,self%nlocs
    call datetime_delete(timer(n))
  enddo
  deallocate(latsr, lonsr, timer)

endif

end subroutine ufo_locs_create

! ------------------------------------------------------------------------------

subroutine ufo_locs_setup(self, nlocs)
implicit none
type(ufo_locs), intent(inout) :: self
integer, intent(in)           :: nlocs

character(len=20) :: fstring
integer :: n

call ufo_locs_delete(self)

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

! ------------------------------------------------------------------------------

subroutine ufo_locs_delete(self)
implicit none
type(ufo_locs), intent(inout) :: self

if (allocated(self%lat)) deallocate(self%lat)
if (allocated(self%lon)) deallocate(self%lon)
if (allocated(self%time)) deallocate(self%time)
if (allocated(self%indx)) deallocate(self%indx)
self%nlocs = 0

end subroutine ufo_locs_delete

! ------------------------------------------------------------------------------
    
subroutine ufo_locs_init(self, obss, t1, t2)
  use kinds
  use datetime_mod
  use twindow_utils_mod
  use obsspace_mod

  implicit none

  type(ufo_locs), intent(inout) :: self
  type(c_ptr), value, intent(in)              :: obss
  type(datetime), intent(in)                  :: t1, t2

  integer :: nlocs

  character(len=*),parameter:: &
     myname = "ufo_locs_init"
  integer :: i
  integer :: tw_nlocs
  integer, dimension(:), allocatable :: tw_indx
  real(kind_real), dimension(:), allocatable :: lon, lat
  type(datetime), dimension(:), allocatable :: date_time

  ! Local copies pre binning
  nlocs = obsspace_get_nlocs(obss)

  allocate(date_time(nlocs), lon(nlocs), lat(nlocs))

!TODO(JG): Add "MetaData" or similar group attribute to all ioda ObsSpace objects
  if (obsspace_has(obss,"MetaData", "datetime")) then
    call obsspace_get_db(obss, "MetaData", "datetime", date_time)
  else
    call obsspace_get_db(obss, "", "datetime", date_time)
  endif

  ! Generate the timing window indices
  allocate(tw_indx(nlocs))
  tw_nlocs = 0
  do i = 1, nlocs
    if (date_time(i) > t1 .and. date_time(i) <= t2) then
      tw_nlocs = tw_nlocs + 1
      tw_indx(tw_nlocs) = i
    endif
  enddo

!TODO(JG): Add "MetaData" or similar group attribute to all ioda ObsSpace objects
  if (obsspace_has(obss,"MetaData", "longitude")) then
    call obsspace_get_db(obss, "MetaData", "longitude", lon)
    call obsspace_get_db(obss, "MetaData", "latitude", lat)
  else
    call obsspace_get_db(obss, "", "longitude", lon)
    call obsspace_get_db(obss, "", "latitude", lat)
  endif

  !Setup ufo locations
  call ufo_locs_setup(self, tw_nlocs)
  do i = 1, tw_nlocs
    self%lon(i)  = lon(tw_indx(i))
    self%lat(i)  = lat(tw_indx(i))
    self%time(i) = date_time(tw_indx(i))
  enddo
  self%indx = tw_indx(1:tw_nlocs)

  do i = 1, nlocs
    call datetime_delete(date_time(i))
  enddo
  deallocate(date_time, lon, lat, tw_indx)

end subroutine ufo_locs_init

! ------------------------------------------------------------------------------

end module ufo_locs_mod
