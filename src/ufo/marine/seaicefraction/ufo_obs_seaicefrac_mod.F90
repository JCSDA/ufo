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
public ufo_obs_seaicefrac_read_oisic
public ufo_obs_seaicefrac_getlocs

! ------------------------------------------------------------------------------

!> Fortran derived type to hold observation locations
type :: ufo_obs_seaicefrac
  integer :: nobs
  real(kind_real), allocatable, dimension(:) :: lat      !< latitude
  real(kind_real), allocatable, dimension(:) :: lon      !< longitude
  real(kind_real), allocatable, dimension(:) :: icefrac  !< total ice concentration
  real(kind_real), allocatable, dimension(:) :: icefrac_err  !< total ice concentration  
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
allocate(self%icefrac(nobs), self%icetmp(nobs), self%icefrac_err(nobs))
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
if (allocated(self%icefrac_err)) deallocate(self%icefrac_err)
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

self%icefrac(:) = 0.0
self%icefrac_err(:) = 0.1
!call random_number(self%icefrac(:))
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

integer :: iunit, nr, nc, nobs, qcnobs
real, allocatable, dimension(:,:)    :: field
real, allocatable, dimension(:)    :: lon, lat, icefrac, qc
integer, allocatable, dimension(:,:) :: ifield

integer :: i, qci
real :: undef = -99999.9

call ufo_obs_seaicefrac_delete(self)

! open netcdf file and read dimensions
call nc_diag_read_init(filename, iunit)
nr = nc_diag_read_get_dim(iunit,'Rows')
nc = nc_diag_read_get_dim(iunit,'Columns')
nobs = nc * nr

allocate(field(nc, nr), ifield(nc, nr), lon(nobs), lat(nobs), icefrac(nobs), qc(nobs))

call nc_diag_read_get_var(iunit, "Latitude", field)
lat = reshape(field, (/nobs/))

call nc_diag_read_get_var(iunit, "Longitude", field)
lon = reshape(field, (/nobs/))

call nc_diag_read_get_var(iunit, "IceConc", field)
icefrac = reshape(field, (/nobs/))

call nc_diag_read_close(filename)

qc = 1
where ( (icefrac.gt.100.0).or.(icefrac.lt.0.0))
   qc=0
end where
qcnobs = sum(qc)
self%nobs = qcnobs

! allocate geovals structure
call ufo_obs_seaicefrac_setup(self, qcnobs)

qci = 0
do i = 1, nobs
   if ( qc(i).eq.1 ) then
      qci = qci + 1
      self%lat(qci) = lat(i)
      self%lon(qci) = lon(i)      
      self%icefrac(qci) = icefrac(i)/100.0
      self%icefrac_err(qci) = 0.1  !< total ice concentration        
      write(101,*)lon(i),lat(i),icefrac(i)
   end if
end do

self%icetmp = 0.0
self%qc = 1
deallocate(field, ifield)

print *, 'in ufo_obs_seaicefrac_read : ', self%nobs, nobs
!call abor1_ftn("======================")
end subroutine ufo_obs_seaicefrac_read

! ------------------------------------------------------------------------------
! ------------------------------------------------------------------------------

subroutine ufo_obs_seaicefrac_read_oisic(filename, self)
use nc_diag_read_mod, only: nc_diag_read_get_var, nc_diag_read_get_dim
use nc_diag_read_mod, only: nc_diag_read_init, nc_diag_read_close
use netcdf
implicit none
character(max_string), intent(in)   :: filename
type(ufo_obs_seaicefrac), intent(inout) :: self

integer :: iunit, nlon, nlat, nobs, qcnobs
real, allocatable    :: tmp_lon(:), tmp_lat(:), tmp_sic(:,:,:), qc(:)
real, allocatable    :: lon(:,:), lat(:,:)

integer :: varid, ierr
integer :: dimids(3)
integer :: i, qci, ii, jj
real :: undef = -99999.9

call ufo_obs_seaicefrac_delete(self)

! open netcdf file and read dimensions
call nc_diag_read_init(filename, iunit)
nlon = nc_diag_read_get_dim(iunit,'lon')
nlat = nc_diag_read_get_dim(iunit,'lat')
print *,'nlon,nlat:',nlon,nlat
allocate(tmp_lon(nlon), tmp_lat(nlat), tmp_sic(nlon,nlat,1))
call nc_diag_read_get_var(iunit, "lon", tmp_lon)
call nc_diag_read_get_var(iunit, "lat", tmp_lat)
!call nc_diag_read_get_var(iunit, "icec", tmp_sic)
call nc_diag_read_close(filename)

ierr = nf90_open(filename, nf90_nowrite, iunit)
ierr = nf90_inq_varid(iunit, "icec", varid)
ierr = nf90_get_var(iunit, varid, tmp_sic)
ierr = nf90_close(iunit)

nobs = nlon * nlat
allocate(lon(nlon,nlat),lat(nlon,nlat))
allocate(qc(nobs))
do ii=1,nlon
   do jj=1,nlat
      lon(ii,jj)=tmp_lon(ii)
      lat(ii,jj)=tmp_lat(jj)
   end do
end do
tmp_sic=reshape(tmp_sic,(/1,1,nobs/))
lon=reshape(lon,(/1,nobs/))
lat=reshape(lat,(/1,nobs/))

qc = 1
!!$do i=1,nobs
!!$   !print *,tmp_sic(1,1,i)
!!$   if ((tmp_sic(1,1,i).gt.1.0).or.(tmp_sic(1,1,i).le.0.0)) then
!!$      qc(i)=0
!!$   end if
!!$end do
where ( (tmp_sic(1,1,:).gt.1.0).or.(tmp_sic(1,1,:).lt.0.0))
   qc=0
end where
qcnobs = sum(qc)
self%nobs = qcnobs

print *,'QC:',nobs,qcnobs

! allocate geovals structure
call ufo_obs_seaicefrac_setup(self, qcnobs)

qci = 0
do i = 1, nobs
   if ( qc(i).eq.1 ) then
      qci = qci + 1
      self%lat(qci) = lat(1,i)
      self%lon(qci) = lon(1,i)      
      self%icefrac(qci) = tmp_sic(1,1,i)
      self%icefrac_err(qci) = 0.01
      write(101,*)lon(1,i),lat(1,i),tmp_sic(1,1,i)
   end if
end do

self%icetmp = 0.0
self%qc = 1


end subroutine ufo_obs_seaicefrac_read_oisic

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
