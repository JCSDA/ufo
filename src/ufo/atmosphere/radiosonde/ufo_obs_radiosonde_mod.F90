!
! (C) Copyright 2017 UCAR
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 

!> Fortran module handling radiosonde observation space

module ufo_obs_radiosonde_mod

use kinds
use m_diag_raob, only: diag_raob_header, diag_raob_mass
use fckit_log_module, only : fckit_log

implicit none
private
integer, parameter :: max_string=800

public ufo_obs_radiosonde
public ufo_obs_radiosonde_setup, ufo_obs_radiosonde_delete
public ufo_obs_radiosonde_read, ufo_obs_radiosonde_generate
public ufo_obs_radiosonde_getlocs

! ------------------------------------------------------------------------------

!> Fortran derived type to hold observation locations
type :: ufo_obs_radiosonde
  integer :: nobs
  integer :: nlocs
  type(diag_raob_header)                     :: header
  type(diag_raob_mass), pointer              :: mass(:)
end type ufo_obs_radiosonde

! ------------------------------------------------------------------------------

contains

! ------------------------------------------------------------------------------

subroutine ufo_obs_radiosonde_setup(self, nobs)
implicit none
type(ufo_obs_radiosonde), intent(inout) :: self
integer, intent(in) :: nobs

call ufo_obs_radiosonde_delete(self)

self%nobs = nobs
allocate(self%mass(nobs))

end subroutine ufo_obs_radiosonde_setup

! ------------------------------------------------------------------------------

subroutine ufo_obs_radiosonde_delete(self)
implicit none
type(ufo_obs_radiosonde), intent(inout) :: self

self%nobs = 0
self%nlocs = 0
if (associated(self%mass)) nullify(self%mass)

end subroutine ufo_obs_radiosonde_delete

! ------------------------------------------------------------------------------

subroutine ufo_obs_radiosonde_generate(self, nobs, lat, lon1, lon2)
implicit none
type(ufo_obs_radiosonde), intent(inout) :: self
integer, intent(in) :: nobs
real, intent(in) :: lat, lon1, lon2

character(len=*),parameter :: myname = "ufo_obs_radiosonde_generate"
integer :: i

call ufo_obs_radiosonde_setup(self, nobs)

self%mass(:)%Latitude = lat
do i = 1, nobs
  self%mass(i)%Longitude = lon1 + (i-1)*(lon2-lon1)/(nobs-1)
enddo

print *, myname, self%nobs, self%mass%longitude, self%mass%latitude

end subroutine ufo_obs_radiosonde_generate

! ------------------------------------------------------------------------------

subroutine ufo_obs_radiosonde_read(filename, self)

use nc_diag_read_mod, only: nc_diag_read_get_var, nc_diag_read_get_dim
use nc_diag_read_mod, only: nc_diag_read_init, nc_diag_read_close
use m_diag_raob, only: read_raob_diag_nc_header, read_raob_diag_nc_mass

implicit none
character(max_string), intent(in)   :: filename
type(ufo_obs_radiosonde), intent(inout), target :: self

character(len=*),parameter :: myname = "ufo_obs_radiosonde_read"
integer :: ier

call ufo_obs_radiosonde_delete(self)

call read_raob_diag_nc_header(filename,self%header)
self%nobs=self%header%n_Observations_Mass
allocate(self%mass(self%nobs))
call read_raob_diag_nc_mass(filename,self%header,self%mass,ier)
self%nobs=self%header%n_Observations_Mass
self%nlocs = self%nobs

print*, myname, ': Found this many observations: ', self%nobs
print*, myname, ': Found this many instances:    ', self%nlocs
print*, myname, ': Size of type holding RAOB:    ', size(self%mass)
print*, myname, ': Date of input file:           ', self%header%date
if(self%nobs>0)&
print*, myname, ': Mean observations:            ', sum(self%mass(:)%Observation)/self%nobs

end subroutine ufo_obs_radiosonde_read

! ------------------------------------------------------------------------------

subroutine ufo_obs_radiosonde_getlocs(self, locs)
use ufo_locs_mod
implicit none
type(ufo_obs_radiosonde), intent(in) :: self
type(ufo_locs), intent(inout) :: locs

character(len=*),parameter:: myname = "ufo_obs_radiosonde_getlocs"
character(len=255) :: record
integer :: failed

call ufo_locs_setup(locs, self%nlocs)

failed=0
if(failed==0 .and. size(self%mass(:)%Longitude)==self%nlocs) then
  locs%lon(:) = self%mass(:)%Longitude
else
  failed=1
endif
if(failed==0 .and. size(self%mass(:)%Latitude) ==self%nlocs) then
  locs%lat(:) = self%mass(:)%Latitude
else
  failed=2
endif
if(failed==0 .and. size(self%mass(:)%Time) == self%nlocs) then
  locs%time(:) = self%mass(:)%Time
else
  failed=3
endif
if(failed==0)then
  write(record,*)myname,': allocated/assinged obs-data'
  call fckit_log%info(record)
else
  write(record,*)myname,': failed allocation/assignment of obs-data, ier: ', failed
  call fckit_log%info(record)
  ! should exit in error here
endif

end subroutine ufo_obs_radiosonde_getlocs

! ------------------------------------------------------------------------------

end module ufo_obs_radiosonde_mod
