!
! (C) Copyright 2017 UCAR
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 

!> Fortran module handling sea ice thickness observation space

module ufo_obs_seaicethick_mod

use kinds

implicit none
private
integer, parameter :: max_string=800

public ufo_obs_seaicethick
public ufo_obs_seaicethick_setup, ufo_obs_seaicethick_delete
public ufo_obs_seaicethick_read, ufo_obs_seaicethick_generate
public ufo_obs_seaicethick_getlocs

! ------------------------------------------------------------------------------

!> Fortran derived type to hold observation locations
type :: ufo_obs_seaicethick
  integer :: nobs
  real(kind_real), allocatable, dimension(:) :: lat      !< latitude
  real(kind_real), allocatable, dimension(:) :: lon      !< longitude
  real(kind_real), allocatable, dimension(:) :: icethick  !< ice thickness
  real(kind_real), allocatable, dimension(:) :: icethick_err  !< ice thickness
  real(kind_real), allocatable, dimension(:) :: freeboard  !< ice freeboard    
end type ufo_obs_seaicethick

! ------------------------------------------------------------------------------
contains

! ------------------------------------------------------------------------------

subroutine ufo_obs_seaicethick_setup(self, nobs)
implicit none

type(ufo_obs_seaicethick), intent(inout) :: self
integer, intent(in) :: nobs

call ufo_obs_seaicethick_delete(self)

self%nobs = nobs
allocate(self%lat(nobs), self%lon(nobs))
allocate(self%icethick(nobs),self%icethick_err(nobs))
self%lat = 0.
self%lon = 0.
self%icethick = 0.
self%icethick_err = 0.

end subroutine ufo_obs_seaicethick_setup

! ------------------------------------------------------------------------------

subroutine ufo_obs_seaicethick_delete(self)
implicit none
type(ufo_obs_seaicethick), intent(inout) :: self

self%nobs = 0
if (allocated(self%lat)) deallocate(self%lat)
if (allocated(self%lon)) deallocate(self%lon)
if (allocated(self%icethick)) deallocate(self%icethick)
if (allocated(self%icethick_err)) deallocate(self%icethick_err)

end subroutine ufo_obs_seaicethick_delete

! ------------------------------------------------------------------------------

subroutine ufo_obs_seaicethick_generate(self, nobs, lat, lon1, lon2)
implicit none
type(ufo_obs_seaicethick), intent(inout) :: self
integer, intent(in) :: nobs
real, intent(in) :: lat, lon1, lon2

integer :: i

call ufo_obs_seaicethick_setup(self, nobs)

self%icethick(:) = 3.0
self%icethick_err(:) = 0.1

self%lat(:)     = lat
do i = 1, nobs
  self%lon(i) = lon1 + (i-1)*(lon2-lon1)/(nobs-1)
enddo

end subroutine ufo_obs_seaicethick_generate

! ------------------------------------------------------------------------------

subroutine ufo_obs_seaicethick_read(filename, self)
use nc_diag_read_mod, only: nc_diag_read_get_var, nc_diag_read_get_dim
use nc_diag_read_mod, only: nc_diag_read_init, nc_diag_read_close
implicit none
character(max_string), intent(in)   :: filename
type(ufo_obs_seaicethick), intent(inout) :: self

integer :: iunit, nobs
real, dimension(:), allocatable :: field

call ufo_obs_seaicethick_delete(self)

call nc_diag_read_init(filename, iunit)
nobs = nc_diag_read_get_dim(iunit,'nobs')
self%nobs = nobs

allocate(field(nobs))
allocate(self%lon(nobs), self%lat(nobs), self%icethick(nobs), self%freeboard(nobs))
allocate(self%icethick_err(nobs))

call nc_diag_read_get_var(iunit, "lat", field)
self%lat = field
call nc_diag_read_get_var(iunit, "lon", field)
self%lon = field
call nc_diag_read_get_var(iunit, "thickness", field)
self%icethick = field
call nc_diag_read_get_var(iunit, "freeboard", field)
self%freeboard = field
call nc_diag_read_close(filename)
self%icethick_err = 0.1
end subroutine ufo_obs_seaicethick_read

! ------------------------------------------------------------------------------

subroutine ufo_obs_seaicethick_getlocs(self, locs)
use ufo_locs_mod
implicit none
type(ufo_obs_seaicethick), intent(in) :: self
type(ufo_locs), intent(inout) :: locs

call ufo_locs_setup(locs, self%nobs)
locs%lat = self%lat
locs%lon = self%lon
locs%time = 0.

end subroutine ufo_obs_seaicethick_getlocs

! ------------------------------------------------------------------------------

end module ufo_obs_seaicethick_mod
