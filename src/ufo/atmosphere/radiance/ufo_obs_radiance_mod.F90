!
! (C) Copyright 2017 UCAR
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 

!> Fortran module handling radiance observation space

module ufo_obs_radiance_mod

use kinds
use read_diag, only: set_radiag,&
                     diag_header_fix_list,&
                     diag_header_chan_list,&
                     diag_data_name_list,&
                     read_radiag_header,&
                     set_netcdf_read
use read_diag, only: read_radiag_data,&
                     diag_data_fix_list,&
                     diag_data_extra_list,&
                     diag_data_chan_list,&
                     open_radiag, &
                     close_radiag, &
                     read_all_radiag
use fckit_log_module, only : fckit_log

implicit none
private

character(len=*),parameter :: myname ="radNode_mod"
integer, parameter :: max_string=800

public ufo_obs_radiance
public ufo_obs_radiance_setup, ufo_obs_radiance_delete
public ufo_obs_radiance_read, ufo_obs_radiance_generate
public ufo_obs_radiance_getlocs

! ------------------------------------------------------------------------------

!> Fortran derived type to hold observation locations
type :: ufo_obs_radiance
  integer :: nobs
  integer :: nlocs
  type(diag_header_fix_list )              ::  header_fix
  type(diag_header_chan_list),allocatable  ::  header_chan(:)
  type(diag_data_name_list)                ::  header_name
  type(diag_data_fix_list)   ,allocatable  ::  datafix(:)
  type(diag_data_chan_list)  ,allocatable  ::  datachan(:,:)
  type(diag_data_extra_list) ,allocatable  ::  dataextra(:,:,:)
end type ufo_obs_radiance

! ------------------------------------------------------------------------------

contains

! ------------------------------------------------------------------------------

subroutine ufo_obs_radiance_setup(self, nobs)
implicit none
type(ufo_obs_radiance), intent(inout) :: self
integer, intent(in) :: nobs

call ufo_obs_radiance_delete(self)

self%nobs = nobs
!Allocatable arrays allocated in nc_diag

end subroutine ufo_obs_radiance_setup

! ------------------------------------------------------------------------------

subroutine ufo_obs_radiance_delete(self)
implicit none
type(ufo_obs_radiance), intent(inout) :: self

self%nobs = 0
if (allocated(self%header_chan)) deallocate(self%header_chan)
if (allocated(self%datafix)) deallocate(self%datafix)
if (allocated(self%datachan)) deallocate(self%datachan)
if (allocated(self%dataextra)) deallocate(self%dataextra)

end subroutine ufo_obs_radiance_delete

! ------------------------------------------------------------------------------

subroutine ufo_obs_radiance_generate(self, nobs, lat, lon1, lon2)
implicit none
type(ufo_obs_radiance), intent(inout) :: self
integer, intent(in) :: nobs
real, intent(in) :: lat, lon1, lon2

integer :: i

call ufo_obs_radiance_setup(self, nobs)

self%datafix(:)%Lat = lat
do i = 1, nobs
  self%datafix(:)%Lon = lon1 + (i-1)*(lon2-lon1)/(nobs-1)
enddo

print *, 'in random:', self%nobs, self%datafix(:)%Lat, self%datafix(:)%Lon

end subroutine ufo_obs_radiance_generate

! ------------------------------------------------------------------------------

subroutine ufo_obs_radiance_read(filename, self)
use nc_diag_read_mod, only: nc_diag_read_get_var, nc_diag_read_get_dim
use nc_diag_read_mod, only: nc_diag_read_init, nc_diag_read_close

use m_diag_raob, only: read_raob_diag_nc_header, read_raob_diag_nc_mass

implicit none
character(max_string), intent(in)   :: filename
type(ufo_obs_radiance), intent(inout) :: self

character(len=*),parameter :: myname_ =myname//"*rad_read"
integer :: ier
integer :: luin=0
integer :: npred = 7   
integer :: iversion=30303
logical :: lverbose  = .true.  ! control verbose
logical :: retrieval = .false. ! true when dealing with SST retrievals

call ufo_obs_radiance_delete(self)

call set_netcdf_read(.true.)
call open_radiag(filename, luin)
call set_radiag("version",iversion,ier)

call read_radiag_header(luin,npred,retrieval,self%header_fix,self%header_chan,self%header_name,ier,lverbose)

print*, myname_, ': Found this many channels: ', self%header_fix%nchan
print*, myname_, ': Observation type in file: ', self%header_fix%obstype
print*, myname_, ': Date of input file:       ', self%header_fix%idate


call read_all_radiag(luin, self%header_fix, retrieval, self%datafix, &
                     self%datachan, self%dataextra, self%nobs, ier)

self%nlocs = self%nobs
self%nobs  = self%nobs * self%header_fix%nchan
call close_radiag(filename,luin)
print *, myname_, ' Total number of observations in file: (nobs,nlocs) ', self%nobs, self%nlocs

end subroutine ufo_obs_radiance_read

! ------------------------------------------------------------------------------

subroutine ufo_obs_radiance_getlocs(self, locs)
use ufo_locs_mod
implicit none
type(ufo_obs_radiance), intent(in) :: self
type(ufo_locs), intent(inout) :: locs

character(len=*),parameter:: myname_=myname//"*rad_getlocs"
character(len=255) :: record
integer :: failed

call ufo_locs_setup(locs, self%nlocs)

failed=0
if(failed==0 .and. size(self%datafix(:)%Lat)==self%nlocs) then
  locs%lat(:) = self%datafix(:)%Lat
else
  failed=1
endif
if(failed==0 .and. size(self%datafix(:)%Lon)==self%nlocs) then
  locs%lon(:) = self%datafix(:)%Lon
else
  failed=2
endif
if(failed==0 .and. size(self%datafix(:)%obstime)==self%nlocs) then
  locs%time(:) = self%datafix(:)%obstime
else
  failed=3
endif
if(failed==0)then
  write(record,*)myname_,': allocated/assinged obs-data'
  call fckit_log%info(record)
else
  write(record,*)myname_,': failed allocation/assignment of obs-data, ier: ', failed
  call fckit_log%info(record)
  ! should exit in error here
endif

end subroutine ufo_obs_radiance_getlocs

! ------------------------------------------------------------------------------

end module ufo_obs_radiance_mod
