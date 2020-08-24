! (C) copyright 2020 Met Office
!
! this software is licensed under the terms of the apache licence version 2.0
! which can be obtained at http://www.apache.org/licenses/license-2.0.

!> Fortran module which contains the observation metadata for a single observation

module ufo_rttovonedvarcheck_ob_mod

use kinds
use missing_values_mod
use ufo_geovals_mod
use ufo_rttovonedvarcheck_constants_mod
use ufo_rttovonedvarcheck_utils_mod
use ufo_rttovonedvarcheck_pcemis_mod

implicit none
private

! Ob info type definition
type, public :: ufo_rttovonedvarcheck_ob

  character(len=max_string) :: forward_mod_name !< forward model name (RTTOV only one at the moment)
  integer              :: nlocs !< number of locations = 1
  integer              :: surface_type  !< surface type of observation
  integer, allocatable :: channels_used(:) !< channels used for this observation
  real(kind_real)      :: latitude !< latitude of observation
  real(kind_real)      :: longitude !< longitude of observation
  real(kind_real)      :: elevation  !< elevation above sea level of observation
  real(kind_real)      :: sensor_zenith_angle  !< sensor zenith of observation
  real(kind_real)      :: sensor_azimuth_angle  !< sensor azimuth of observation
  real(kind_real)      :: solar_zenith_angle !< solar zenith of observation
  real(kind_real)      :: solar_azimuth_angle  !< solar azimuth of observation
  real(kind_real)      :: cloudtopp !< cloud top pressure (used in if cloudy retrieval used)
  real(kind_real)      :: cloudfrac !< cloud fraction (used in if cloudy retrieval used)
  real(kind_real), allocatable :: yobs(:) !< satellite BTs
  real(kind_real), allocatable :: emiss(:) !< surface emissivity
  logical              :: retrievecloud  !< flag to turn on retrieve cloud
  logical              :: mwscatt !< flag to use rttov-scatt model through the interface
  logical              :: mwscatt_totalice !< flag to use total ice (rather then ciw) for rttov-scatt simulations
  logical, allocatable :: calc_emiss(:) !< flag to decide if RTTOV calculates emissivity
  type(ufo_rttovonedvarcheck_pcemis), pointer :: pcemis !< pointer to the IR pc emissivity object

contains
  procedure :: setup  => ufo_rttovonedvarcheck_InitOb
  procedure :: init_emiss  => ufo_rttovonedvarcheck_InitEmiss
  procedure :: delete => ufo_rttovonedvarcheck_DeleteOb
  procedure :: info => ufo_rttovonedvarcheck_PrintOb

end type ufo_rttovonedvarcheck_ob

contains

!-------------------------------------------------------------------------------
!> Initialize observation object
!!
!! \author Met Office
!!
!! \date 09/06/2020: Created
!!
subroutine ufo_rttovonedvarcheck_InitOb(self, & ! out
                                        nchans) ! in 

implicit none

! subroutine arguments:
class(ufo_rttovonedvarcheck_ob), intent(out) :: self !< observation metadata type
integer, intent(in) :: nchans !< number of channels used for this particular observation

character(len=*), parameter :: routinename = "ufo_rttovonedvarcheck_InitOb"
real(kind_real) :: missing

missing = missing_value(missing)

call self % delete()

allocate(self % yobs(nchans))
allocate(self % channels_used(nchans))
allocate(self % emiss(nchans))
allocate(self % calc_emiss(nchans))

self % yobs(nchans) = missing
self % emiss(:) = 0.0
self % calc_emiss(:) = .false.

end subroutine ufo_rttovonedvarcheck_InitOb

!------------------------------------------------------------------------------
!> Initialize the emissivity arrays within the single obs structure
!!
!! \author Met Office
!!
!! \date 09/06/2020: Created
!!
subroutine ufo_rttovonedvarcheck_InitEmiss(self, config, geovals)

implicit none

! subroutine arguments:
class(ufo_rttovonedvarcheck_ob), intent(inout) :: self !< observation metadata type
type(ufo_rttovonedvarcheck), intent(in) :: config !< main rttovonedvarcheck type
type(ufo_geovals), intent(in)           :: geovals !< model data at obs location

type(ufo_geoval), pointer    :: geoval

! ----------------------------------
! 1.0 Get surface type from geovals
!-----------------------------------
call ufo_geovals_get_var(geovals, "surface_type", geoval)
self % surface_type = geoval%vals(1, 1)

!-------------
! 2.1 Defaults
!-------------

! Only calculate in RTTOV over sea
if (self % surface_type == RTSea) then
  self % calc_emiss(:) = .true.
else
  self % calc_emiss(:) = .false.
end if

! The default emissivity for land is a very crude estimate - the same
! for all surface types and all frequencies. However, we do not use
! channels which see the surface over land where we rely on this default.
self % emiss(:) = 0.0
if (self % surface_type == RTLand) then
  self % emiss(:) = config % EmissLandDefault
else if (self % surface_type == RTIce) then
  self % emiss(:) = config % EmissSeaIceDefault
end if

end subroutine ufo_rttovonedvarcheck_InitEmiss

!------------------------------------------------------------------------------
!> Delete the single observation object
!!
!! \author Met Office
!!
!! \date 09/06/2020: Created
!!
subroutine ufo_rttovonedvarcheck_DeleteOb(self)    ! inout

implicit none

! subroutine arguments:
class(ufo_rttovonedvarcheck_ob), intent(inout) :: self !< observation metadata type

character(len=*), parameter :: routinename = "ufo_rttovonedvarcheck_DeleteOb"

self % nlocs = 1
self % latitude = 0.0
self % longitude = 0.0
self % elevation = 0.0
self % surface_type = 0
self % sensor_zenith_angle = 0.0
self % sensor_azimuth_angle = 0.0
self % solar_zenith_angle = 0.0
self % solar_azimuth_angle = 0.0
self % cloudtopp = 500.0
self % cloudfrac = 0.0
self % retrievecloud = .false.
self % mwscatt = .false.
self % mwscatt_totalice = .false.

if (allocated(self % yobs))          deallocate(self % yobs)
if (allocated(self % channels_used)) deallocate(self % channels_used)
if (allocated(self % emiss))         deallocate(self % emiss)
if (allocated(self % calc_emiss))    deallocate(self % calc_emiss)

self % pcemis => null()

end subroutine ufo_rttovonedvarcheck_DeleteOb

!------------------------------------------------------------------------------
!> Print information about the single observation object
!!
!! \author Met Office
!!
!! \date 09/06/2020: Created
!!
subroutine ufo_rttovonedvarcheck_PrintOb(self)

implicit none

! subroutine arguments:
class(ufo_rttovonedvarcheck_ob), intent(inout) :: self !< observation metadata type

character(len=20) :: surface_type

if (self % surface_type == RTLand) then
  surface_type = "land"
else if (self % surface_type == RTSea) then
  surface_type = "sea"
else if (self % surface_type == RTIce) then
  surface_type = "ice"
else
  surface_type = "unknown"
end if

write(*,"(A,2F8.2)") "Lat,Long:",self % latitude, self % longitude
write(*,*) "Surface type for RTTOV: ",surface_type
write(*,"(A,F8.2)") "Surface height:",self % elevation
write(*,"(A,F8.2)") "Satellite zenith angle: ",self % sensor_zenith_angle
write(*,"(A,F8.2)") "Solar zenith angle: ",self % solar_zenith_angle

end subroutine

!-------------------------------------------------------------------------------

end module ufo_rttovonedvarcheck_ob_mod
