! (C) copyright 2020 Met Office
!
! this software is licensed under the terms of the apache licence version 2.0
! which can be obtained at http://www.apache.org/licenses/license-2.0.

!> Fortran module which contains the observation metadata for a single observation

module ufo_rttovonedvarcheck_obinfo_mod

use kinds
use missing_values_mod
use ufo_geovals_mod
use ufo_rttovonedvarcheck_constants_mod
use ufo_rttovonedvarcheck_utils_mod

implicit none
private

! Ob info type definition
type, public :: obinfo_type

  character(len=max_string) :: forward_mod_name !< forward model name (RTTOV only one at the moment)
  integer         :: nlocs !< number of locations = 1
  real(kind_real) :: latitude !< latitude of observation
  real(kind_real) :: longitude !< longitude of observation
  real(kind_real) :: elevation  !< elevation above sea level of observation
  integer         :: surface_type  !< surface type of observation
  real(kind_real) :: sensor_zenith_angle  !< sensor zenith of observation
  real(kind_real) :: sensor_azimuth_angle  !< sensor azimuth of observation
  real(kind_real) :: solar_zenith_angle !< solar zenith of observation
  real(kind_real) :: solar_azimuth_angle  !< solar azimuth of observation
  real(kind_real) :: cloudtopp !< cloud top pressure (used in if cloudy retrieval used)
  real(kind_real) :: cloudfrac !< cloud fraction (used in if cloudy retrieval used)
  logical         :: retrievecloud  !< flag to turn on retrieve cloud
  logical         :: mwscatt !< flag to use rttov-scatt model through the interface
  logical         :: mwscatt_totalice !< flag to use total ice (rather then ciw) for rttov-scatt simulations
  real(kind_real), allocatable :: yobs(:) !< satellite BTs
  integer, allocatable :: channels_used(:) !< channels used for this observation
  real(kind_real), allocatable :: emiss(:) !< surface emissivity
  logical, allocatable :: calc_emiss(:) !< flag to decide if RTTOV calculates emissivity

contains
  procedure :: setup  => ufo_rttovonedvarcheck_InitObInfo
  procedure :: init_emiss  => ufo_rttovonedvarcheck_InitEmiss
  procedure :: delete => ufo_rttovonedvarcheck_DeleteObInfo

end type

contains

!-------------------------------------------------------------------------------
!> Initialize observation object
!!
!! \author Met Office
!!
!! \date 09/06/2020: Created
!!
subroutine ufo_rttovonedvarcheck_InitObInfo(ob_info, & ! out
                                            nchans)    ! in

implicit none

! subroutine arguments:
class(obinfo_type), intent(out) :: ob_info !< observation metadata type
integer, intent(in)  :: nchans !< number of channels used for this particular observation

character(len=*), parameter :: routinename = "ufo_rttovonedvarcheck_InitObInfo"
real(kind_real) :: missing

missing = missing_value(missing)

ob_info % nlocs = 1
ob_info % latitude = 0.0
ob_info % longitude = 0.0
ob_info % elevation = 0.0
ob_info % surface_type = 0
ob_info % sensor_zenith_angle = 0.0
ob_info % sensor_azimuth_angle = 0.0
ob_info % solar_zenith_angle = 0.0
ob_info % solar_azimuth_angle = 0.0
ob_info % cloudtopp = 500.0
ob_info % cloudfrac = 0.0
ob_info % retrievecloud = .false.
ob_info % mwscatt = .false.
ob_info % mwscatt_totalice = .false.

allocate(ob_info % yobs(nchans))
allocate(ob_info % channels_used(nchans))
allocate(ob_info % emiss(nchans))
allocate(ob_info % calc_emiss(nchans))

ob_info % yobs(nchans) = missing
ob_info % emiss(:) = 0.0
ob_info % calc_emiss(:) = .false.

end subroutine ufo_rttovonedvarcheck_InitObInfo

!------------------------------------------------------------------------------
!> Initialize the emissivity arrays within the single obs structure
!!
!! \author Met Office
!!
!! \date 09/06/2020: Created
!!
subroutine ufo_rttovonedvarcheck_InitEmiss(ob_info, self, geovals)

implicit none

! subroutine arguments:
class(obinfo_type), intent(inout)       :: ob_info !< observation metadata type
type(ufo_rttovonedvarcheck), intent(in) :: self !< main rttovonedvarcheck type
type(ufo_geovals), intent(in)           :: geovals !< model data at obs location

type(ufo_geoval), pointer    :: geoval

! ----------------------------------
! 1.0 Get surface type from geovals
!-----------------------------------
call ufo_geovals_get_var(geovals, "surface_type", geoval)
ob_info % surface_type = geoval%vals(1, 1)

!-------------
! 2.1 Defaults
!-------------

! Only calculate in RTTOV over sea
if (ob_info % surface_type == RTSea) then
  ob_info % calc_emiss(:) = .true.
else
  ob_info % calc_emiss(:) = .false.
end if

! The default emissivity for land is a very crude estimate - the same
! for all surface types and all frequencies. However, we do not use
! channels which see the surface over land where we rely on this default.
ob_info % emiss(:) = 0.0
if (ob_info % surface_type == RTLand) then
  ob_info % emiss(:) = self % EmissLandDefault
else if (ob_info % surface_type == RTIce) then
  ob_info % emiss(:) = self % EmissSeaIceDefault
end if

end subroutine ufo_rttovonedvarcheck_InitEmiss

!------------------------------------------------------------------------------
!> Delete the single observation object
!!
!! \author Met Office
!!
!! \date 09/06/2020: Created
!!
subroutine ufo_rttovonedvarcheck_DeleteObInfo(ob_info)    ! inout

implicit none

! subroutine arguments:
class(obinfo_type), intent(inout) :: ob_info !< observation metadata type

character(len=*), parameter :: routinename = "ufo_rttovonedvarcheck_InitObInfo"

ob_info % nlocs = 1
ob_info % latitude = 0.0
ob_info % longitude = 0.0
ob_info % elevation = 0.0
ob_info % surface_type = 0
ob_info % sensor_zenith_angle = 0.0
ob_info % sensor_azimuth_angle = 0.0
ob_info % solar_zenith_angle = 0.0
ob_info % solar_azimuth_angle = 0.0
ob_info % cloudtopp = 500.0
ob_info % cloudfrac = 0.0
ob_info % retrievecloud = .false.
ob_info % mwscatt = .false.
ob_info % mwscatt_totalice = .false.

if (allocated(ob_info % yobs))          deallocate(ob_info % yobs)
if (allocated(ob_info % channels_used)) deallocate(ob_info % channels_used)
if (allocated(ob_info % emiss))         deallocate(ob_info % emiss)
if (allocated(ob_info % calc_emiss))    deallocate(ob_info % calc_emiss)

end subroutine ufo_rttovonedvarcheck_DeleteObInfo

!-------------------------------------------------------------------------------

end module ufo_rttovonedvarcheck_obinfo_mod
