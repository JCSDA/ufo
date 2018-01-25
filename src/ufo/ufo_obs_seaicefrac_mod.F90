!
! (C) Copyright 2017 UCAR
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 

!> Fortran module handling observation locations

module ufo_obs_seaicefrac_mod

use kinds

implicit none
private
integer, parameter :: max_string=800

public ufo_obs_seaicefrac
public ufo_obs_seaicefrac_setup, ufo_obs_seaicefrac_delete
public ufo_obs_seaicefrac_read, ufo_obs_seaicefrac_generate
public ufo_obs_seaicefrac_getlocs

! ------------------------------------------------------------------------------

!> Fortran derived type to hold observation locations
type :: ufo_obs_seaicefrac
  integer :: nobs
  real(kind_real), allocatable, dimension(:) :: lat      !< latitude
  real(kind_real), allocatable, dimension(:) :: lon      !< longitude
  real(kind_real), allocatable, dimension(:) :: icefrac  !< total ice concentration
  real(kind_real), allocatable, dimension(:) :: icetmp   !< ice temperature (?)
  integer,         allocatable, dimension(:) :: qc       !< QC flag (from file?)
end type ufo_obs_seaicefrac


! ------------------------------------------------------------------------------
contains

! ------------------------------------------------------------------------------

subroutine ufo_obs_seaicefrac_setup(self, nobs)
implicit none

type(ufo_obs_seaicefrac), intent(inout) :: self
integer, intent(in) :: nobs

call ufo_obs_seaicefrac_delete(self)

self%nobs = nobs
allocate(self%lat(nobs), self%lon(nobs))
allocate(self%icefrac(nobs), self%icetmp(nobs))
allocate(self%qc(nobs))
self%lat = 0.
self%lon = 0.
self%icefrac = 0.
self%icetmp  = 0.
self%qc = 0

end subroutine ufo_obs_seaicefrac_setup

! ------------------------------------------------------------------------------

subroutine ufo_obs_seaicefrac_delete(self)
implicit none
type(ufo_obs_seaicefrac), intent(inout) :: self

self%nobs = 0
if (allocated(self%lat)) deallocate(self%lat)
if (allocated(self%lon)) deallocate(self%lon)
if (allocated(self%icefrac)) deallocate(self%icefrac)
if (allocated(self%icetmp))  deallocate(self%icetmp)
if (allocated(self%qc))  deallocate(self%qc)

end subroutine ufo_obs_seaicefrac_delete

! ------------------------------------------------------------------------------

subroutine ufo_obs_seaicefrac_generate(self, nobs, lat, lon1, lon2)
implicit none
type(ufo_obs_seaicefrac), intent(inout) :: self
integer, intent(in) :: nobs
real, intent(in) :: lat, lon1, lon2

integer :: i

call ufo_obs_seaicefrac_setup(self, nobs)

self%icefrac(:) = 1.
self%icetmp(:)  = 0.
self%qc(:)      = 1.

self%lat(:)     = lat
do i = 1, nobs
  self%lon(i) = lon1 + (i-1)*(lon2-lon1)/(nobs-1)
enddo

print *, 'in random:', self%nobs, self%lon, self%lat

end subroutine ufo_obs_seaicefrac_generate

! ------------------------------------------------------------------------------

subroutine ufo_obs_seaicefrac_read(filename, self)
use nc_diag_read_mod, only: nc_diag_read_get_var, nc_diag_read_get_dim
use nc_diag_read_mod, only: nc_diag_read_init, nc_diag_read_close
implicit none
character(max_string), intent(in)   :: filename
type(ufo_obs_seaicefrac), intent(inout) :: self

integer :: iunit, nr, nc, nobs
real, allocatable, dimension(:,:)    :: field
integer, allocatable, dimension(:,:) :: ifield

call ufo_obs_seaicefrac_delete(self)

! open netcdf file and read dimensions
call nc_diag_read_init(filename, iunit)
nr = nc_diag_read_get_dim(iunit,'Rows')
nc = nc_diag_read_get_dim(iunit,'Columns')
nobs = nc * nr

! allocate geovals structure
call ufo_obs_seaicefrac_setup(self, nobs)

allocate(field(nc, nr), ifield(nc, nr))

call nc_diag_read_get_var(iunit, "Latitude", field)
self%lat = reshape(field, (/nobs/))

call nc_diag_read_get_var(iunit, "Longitude", field)
self%lon = reshape(field, (/nobs/))

call nc_diag_read_get_var(iunit, "IceConc", field)
self%icefrac = reshape(field, (/nobs/))

call nc_diag_read_get_var(iunit, "IceSrfTemp", field)
self%icetmp = reshape(field, (/nobs/))

call nc_diag_read_get_var(iunit, "QCFlags", ifield)
self%qc = reshape(ifield, (/nobs/))
deallocate(field, ifield)

call nc_diag_read_close(filename)

print *, 'in read: ', self%nobs

end subroutine ufo_obs_seaicefrac_read

! ------------------------------------------------------------------------------

subroutine ufo_obs_seaicefrac_getlocs(self, locs)
use ufo_locs_mod
implicit none
type(ufo_obs_seaicefrac), intent(in) :: self
type(ufo_locs), intent(inout) :: locs

call ufo_locs_setup(locs, self%nobs)
locs%lat = self%lat
locs%lon = self%lon
locs%time = 0.

end subroutine ufo_obs_seaicefrac_getlocs

! ------------------------------------------------------------------------------

end module ufo_obs_seaicefrac_mod
