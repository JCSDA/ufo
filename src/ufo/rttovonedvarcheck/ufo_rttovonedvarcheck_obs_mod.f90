! (C) copyright 2020 Met Office
!
! this software is licensed under the terms of the apache licence version 2.0
! which can be obtained at http://www.apache.org/licenses/license-2.0.

!> Fortran module which contains the observation data for a all the obs space

module ufo_rttovonedvarcheck_obs_mod

use kinds
use iso_c_binding
use missing_values_mod
use obsspace_mod
use oops_variables_mod
use ufo_geovals_mod
use ufo_rttovonedvarcheck_constants_mod
use ufo_rttovonedvarcheck_utils_mod

implicit none
private

! Ob info type definition
type, public :: obs_type

integer                      :: iloc
integer, allocatable         :: QCflags(:,:)    ! current qc flags needed for channel selection
real(kind_real), allocatable :: yobs(:,:)       ! observation value from obs files
real(kind_real), allocatable :: yerr(:,:)       ! observation error from obs files
real(kind_real), allocatable :: ybias(:,:)      ! observation bias from obs files
real(kind_real), allocatable :: lat(:)          ! observation latitude
real(kind_real), allocatable :: lon(:)          ! observation longitude
real(kind_real), allocatable :: elevation(:)    ! observation elevation
real(kind_real), allocatable :: sat_zen(:)      ! observation satellite zenith angle
real(kind_real), allocatable :: sat_azi(:)      ! observation satellite azimuth angle
real(kind_real), allocatable :: sol_zen(:)      ! observation solar zenith angle
real(kind_real), allocatable :: sol_azi(:)      ! observation solar azimuth angle
real(kind_real), allocatable :: emissivity(:,:) ! initial surface emissivity

contains
  procedure :: setup  => ufo_rttovonedvarcheck_obs_setup
  procedure :: delete => ufo_rttovonedvarcheck_obs_delete

end type

contains

!-------------------------------------------------------------------------------
!> Initialize observation object
!!
!! \author Met Office
!!
!! \date 09/06/2020: Created
!!
subroutine ufo_rttovonedvarcheck_obs_setup(self,     & ! out
                                           obsspace, & ! in
                                           nchans,   & ! in
                                           vars,     & ! in
                                           mwemiss,  & ! in
                                           iremiss)    ! in

implicit none

! subroutine arguments:
class(obs_type), intent(out)     :: self     !< observation metadata type
type(c_ptr), value, intent(in)   :: obsspace !< observation database pointer
integer, intent(in)              :: nchans   !< total number of channels
type(oops_variables), intent(in) :: vars     !< channels for 1D-Var
logical                          :: mwemiss  !< flag for reading in microwave emissivity
logical                          :: iremiss  !< flag for reading in infrared emissivity

character(len=*), parameter :: routinename = "ufo_rttovonedvarcheck_obs_init"
real(kind_real)             :: missing
integer                     :: jvar    !< counters
character(len=max_string)   :: var
character(len=max_string)   :: varname
logical                     :: variable_present = .false.
logical                     :: model_surface_present = .false.

missing = missing_value(missing)
self % iloc = obsspace_get_nlocs(obsspace)

! allocate arrays
allocate(self % yobs(nchans, self % iloc))
allocate(self % yerr(nchans, self % iloc))
allocate(self % ybias(nchans, self % iloc))
allocate(self % QCflags(nchans, self % iloc))
allocate(self % lat(self % iloc))
allocate(self % lon(self % iloc))
allocate(self % elevation(self % iloc))
allocate(self % sat_zen(self % iloc))
allocate(self % sat_azi(self % iloc))
allocate(self % sol_zen(self % iloc))
allocate(self % sol_azi(self % iloc))

! initialize arrays
self % yobs(:,:) = 0.0
self % yerr(:,:) = 0.0
self % ybias(:,:) = 0.0
self % QCflags(:,:) = 0
self % lat(:) = 0.0
self % lon(:) = 0.0
self % elevation(:) = 0.0
self % sat_zen(:) = 0.0
self % sat_azi(:) = 0.0
self % sol_zen(:) = 0.0
self % sol_azi(:) = 0.0

! Setup optional arrays
if (mwemiss .or. iremiss) then
  allocate(self % emissivity(nchans, self % iloc))
  self % emissivity(:,:) = 0.0
end if

! read in observations and associated errors / biases for full ObsSpace
do jvar = 1, nchans

  var = vars % variable(jvar)
  call obsspace_get_db(obsspace, "FortranQC", trim(var), self % QCflags(jvar,:))
  call obsspace_get_db(obsspace, "ObsValue",  trim(var), self % yobs(jvar,:))
  call obsspace_get_db(obsspace, "ObsError",  trim(var), self % yerr(jvar,:))

  ! Optionally get the observation bias
  variable_present = obsspace_has(obsspace, "ObsBias", trim(var))
  if (variable_present) then
    call obsspace_get_db(obsspace, "ObsBias",   trim(var), self % ybias(jvar,:))
  end if

  ! If emissivity part of the state vector read it from the first guess from the db
  if (mwemiss) then
    call obsspace_get_db(obsspace, "MwEmiss",   trim(var), self % emissivity(jvar,:))
  else if (iremiss) then
    call obsspace_get_db(obsspace, "IREmiss",   trim(var), self % emissivity(jvar,:))
  end if

end do

if (.not. variable_present) write(*,*) "Using uncorrected brightness temperature"

! Add bias correction to the observations
self % yobs = self % yobs + self % ybias

! Read in prerequisites
call obsspace_get_db(obsspace, "MetaData", "latitude", self % lat(:))
call obsspace_get_db(obsspace, "MetaData", "longitude", self % lon(:))
call obsspace_get_db(obsspace, "MetaData", "sensor_zenith_angle", self % sat_zen(:))

! Read in optional angles
variable_present = obsspace_has(obsspace, "MetaData", "sensor_azimuth_angle")
if (variable_present) then
  call obsspace_get_db(obsspace, "MetaData", "sensor_azimuth_angle", self % sat_azi(:))
end if

variable_present = obsspace_has(obsspace, "MetaData", "solar_zenith_angle")
if (variable_present) then
  call obsspace_get_db(obsspace, "MetaData", "solar_zenith_angle", self % sol_zen(:))
end if

variable_present = obsspace_has(obsspace, "MetaData", "solar_azimuth_angle")
if (variable_present) then
  call obsspace_get_db(obsspace, "MetaData", "solar_azimuth_angle", self % sol_azi(:))
end if

! Read in elevation for all obs
variable_present = obsspace_has(obsspace, "MetaData", "elevation")
if (variable_present) then
  call obsspace_get_db(obsspace, "MetaData", "elevation", self % elevation(:))
else
  model_surface_present = obsspace_has(obsspace, "MetaData", "model_surface")
  if (model_surface_present) then
    call obsspace_get_db(obsspace, "MetaData", "model_surface", self % elevation(:))
  else
    self % elevation(:) = 0.0
  end if
end if

end subroutine ufo_rttovonedvarcheck_obs_setup

!------------------------------------------------------------------------------
!> Delete the observation object
!!
!! \author Met Office
!!
!! \date 09/06/2020: Created
!!
subroutine ufo_rttovonedvarcheck_obs_delete(self)    ! inout

implicit none

! subroutine arguments:
class(obs_type), intent(inout) :: self !< observation metadata type

character(len=*), parameter :: routinename = "ufo_rttovonedvarcheck_obs_delete"

! deallocate arrays
if (allocated(self % QCflags))    deallocate(self % QCflags)
if (allocated(self % yobs))       deallocate(self % yobs)
if (allocated(self % yerr))       deallocate(self % yerr)
if (allocated(self % ybias))      deallocate(self % ybias)
if (allocated(self % lat))        deallocate(self % lat)
if (allocated(self % lon))        deallocate(self % lon)
if (allocated(self % elevation))  deallocate(self % elevation)
if (allocated(self % sat_zen))    deallocate(self % sat_zen)
if (allocated(self % sat_azi))    deallocate(self % sat_azi)
if (allocated(self % sol_zen))    deallocate(self % sol_zen)
if (allocated(self % sol_azi))    deallocate(self % sol_azi)
if (allocated(self % emissivity)) deallocate(self % emissivity)

end subroutine ufo_rttovonedvarcheck_obs_delete

!-------------------------------------------------------------------------------

end module ufo_rttovonedvarcheck_obs_mod
