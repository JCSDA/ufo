!
! (C) Copyright 2017 UCAR
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 

!> Fortran module handling observation locations

module ufo_obs_stericheight_mod

use kinds

implicit none
private
integer, parameter :: max_string=800

public ufo_obs_stericheight
public ufo_obs_stericheight_setup, ufo_obs_stericheight_delete
public ufo_obs_stericheight_read, ufo_obs_stericheight_generate
public ufo_obs_stericheight_getlocs

! ------------------------------------------------------------------------------

!> Fortran derived type to hold observation locations
type :: ufo_obs_stericheight
  integer :: nobs
  real(kind_real), allocatable, dimension(:) :: lat      !< latitude
  real(kind_real), allocatable, dimension(:) :: lon      !< longitude
  real(kind_real), allocatable, dimension(:) :: adt      !< ADT
  real(kind_real), allocatable, dimension(:) :: adt_err  !< ADT error  
  real(kind_real), allocatable, dimension(:) :: madt     !< Mean ADT
  integer,         allocatable, dimension(:) :: qc       !< QC flag (from file?)
end type ufo_obs_stericheight

! ------------------------------------------------------------------------------
contains

! ------------------------------------------------------------------------------

subroutine ufo_obs_stericheight_setup(self, nobs)
implicit none

type(ufo_obs_stericheight), intent(inout) :: self
integer, intent(in) :: nobs

call ufo_obs_stericheight_delete(self)

self%nobs = nobs
allocate(self%lat(nobs), self%lon(nobs))
allocate(self%adt(nobs), self%madt(nobs), self%adt_err(nobs))
allocate(self%qc(nobs))
self%lat = 0.
self%lon = 0.
self%adt = 0.
self%madt  = 0.
self%adt_err  = 0.
self%qc = 0

end subroutine ufo_obs_stericheight_setup

! ------------------------------------------------------------------------------

subroutine ufo_obs_stericheight_delete(self)
implicit none
type(ufo_obs_stericheight), intent(inout) :: self

self%nobs = 0
if (allocated(self%lat)) deallocate(self%lat)
if (allocated(self%lon)) deallocate(self%lon)
if (allocated(self%adt)) deallocate(self%adt)
if (allocated(self%adt_err)) deallocate(self%adt_err)
if (allocated(self%madt))  deallocate(self%madt)
if (allocated(self%qc))  deallocate(self%qc)

end subroutine ufo_obs_stericheight_delete

! ------------------------------------------------------------------------------

subroutine ufo_obs_stericheight_generate(self, nobs, lat, lon1, lon2)
implicit none
type(ufo_obs_stericheight), intent(inout) :: self
integer, intent(in) :: nobs
real, intent(in) :: lat, lon1, lon2

integer :: i

call ufo_obs_stericheight_setup(self, nobs)

self%adt(:) = 0.0
self%adt_err(:) = 0.1
self%madt(:)  = 0.
self%qc(:)      = 1.

self%lat(:)     = lat
do i = 1, nobs
  self%lon(i) = lon1 + (i-1)*(lon2-lon1)/(nobs-1)
enddo

end subroutine ufo_obs_stericheight_generate

! ------------------------------------------------------------------------------

subroutine ufo_obs_stericheight_read(filename, self)
use nc_diag_read_mod, only: nc_diag_read_get_var, nc_diag_read_get_dim
use nc_diag_read_mod, only: nc_diag_read_init, nc_diag_read_close
use ncd_kinds, only: i_short
implicit none
character(max_string), intent(in)   :: filename
type(ufo_obs_stericheight), intent(inout) :: self

integer :: iunit, nt, nr, nc, nobs, qcnobs
integer(i_short), allocatable, dimension(:)    :: ssh
real(kind_real), allocatable, dimension(:)    :: qc
integer, allocatable, dimension(:)    :: lon, lat, mdt

integer :: i, qci
real :: undef = -99999.9

call ufo_obs_stericheight_delete(self)

! open netcdf file and read dimensions
call nc_diag_read_init(filename, iunit)
nt = nc_diag_read_get_dim(iunit,'time')
nobs = nt
allocate(lon(nobs), lat(nobs), ssh(nobs), qc(nobs), mdt(nobs))
call nc_diag_read_get_var(iunit, "lat", lat)
call nc_diag_read_get_var(iunit, "lon", lon)
call nc_diag_read_get_var(iunit, "ssha", ssh)
call nc_diag_read_get_var(iunit, "mean_topography", mdt)
call nc_diag_read_close(filename)

qc = 1
where ( (lat.eq.0).or.(ssh.gt.9999).or.(mdt.gt.9999) )
   qc=0
end where

qcnobs = sum(qc)
self%nobs = qcnobs
print *,qcnobs, nobs

! allocate geovals structure
call ufo_obs_stericheight_setup(self, qcnobs)

qci = 0
do i = 1, nobs
   if ( qc(i).eq.1 ) then
      qci = qci + 1
      self%lat(qci) = lat(i)*1e-6      
      self%lon(qci) = lon(i)*1e-6
      if (self%lon(qci)>80.0) self%lon(qci) = self%lon(qci) - 360.0
      self%adt(qci) = ssh(i)*0.001 + mdt(i)*0.0001
      self%adt_err(qci) = 0.1
      write(101,*)self%lon(i),self%lat(i),self%adt(i)
   end if
end do

self%madt = 0.0
self%qc = 1
deallocate(lon, lat, ssh, qc, mdt)

print *, 'in read: ', self%nobs, nobs

end subroutine ufo_obs_stericheight_read

! ------------------------------------------------------------------------------

subroutine ufo_obs_stericheight_getlocs(self, locs)
use ufo_locs_mod
implicit none
type(ufo_obs_stericheight), intent(in) :: self
type(ufo_locs), intent(inout) :: locs

call ufo_locs_setup(locs, self%nobs)
locs%lat = self%lat
locs%lon = self%lon
locs%time = 0.

end subroutine ufo_obs_stericheight_getlocs

! ------------------------------------------------------------------------------

end module ufo_obs_stericheight_mod
