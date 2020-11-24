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
  integer, allocatable :: channels_all(:) !< all channels used for output
  real(kind_real)      :: latitude !< latitude of observation
  real(kind_real)      :: longitude !< longitude of observation
  real(kind_real)      :: elevation  !< elevation above sea level of observation
  real(kind_real)      :: sensor_zenith_angle  !< sensor zenith of observation
  real(kind_real)      :: sensor_azimuth_angle  !< sensor azimuth of observation
  real(kind_real)      :: solar_zenith_angle !< solar zenith of observation
  real(kind_real)      :: solar_azimuth_angle  !< solar azimuth of observation
  real(kind_real)      :: cloudtopp !< cloud top pressure (used in if cloudy retrieval used)
  real(kind_real)      :: cloudfrac !< cloud fraction (used in if cloudy retrieval used)
  real(kind_real)      :: final_cost !< final cost at solution
  real(kind_real), allocatable :: yobs(:) !< satellite BTs
  real(kind_real), allocatable :: emiss(:) !< surface emissivity
  real(kind_real), allocatable :: output_profile(:) !< retrieved state at converge as profile vector
  real(kind_real), allocatable :: output_BT(:) !< Brightness temperatures using retrieved state
  logical              :: retrievecloud  !< flag to turn on retrieve cloud
  logical              :: mwscatt !< flag to use rttov-scatt model through the interface
  logical              :: mwscatt_totalice !< flag to use total ice (rather then ciw) for rttov-scatt simulations
  logical, allocatable :: calc_emiss(:) !< flag to decide if RTTOV calculates emissivity
  type(ufo_rttovonedvarcheck_pcemis), pointer :: pcemis !< pointer to the IR pc emissivity object

contains
  procedure :: setup  => ufo_rttovonedvarcheck_InitOb
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
                                        nchans, &  ! in
                                        nprofelements, & ! in
                                        nchans_all ) ! in

implicit none

! subroutine arguments:
class(ufo_rttovonedvarcheck_ob), intent(out) :: self !< observation metadata type
integer, intent(in) :: nchans !< number of channels used for this particular observation
integer, intent(in) :: nprofelements !< number of profile elements used
integer :: nchans_all !< Size of all channels in ObsSpace

character(len=*), parameter :: routinename = "ufo_rttovonedvarcheck_InitOb"
real(kind_real) :: missing

missing = missing_value(missing)

call self % delete()

allocate(self % yobs(nchans))
allocate(self % channels_used(nchans))
allocate(self % channels_all(nchans_all))
allocate(self % emiss(nchans))
allocate(self % output_profile(nprofelements))
allocate(self % output_BT(nchans_all))
allocate(self % calc_emiss(nchans))

self % yobs(:) = missing
self % emiss(:) = 0.0
self % output_profile(:) = missing
self % output_BT(:) = missing
self % calc_emiss(:) = .true.

end subroutine ufo_rttovonedvarcheck_InitOb

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
real(kind_real) :: missing

missing = missing_value(missing)

self % nlocs = 1
self % latitude = missing
self % longitude = missing
self % elevation = missing
self % surface_type = 0
self % sensor_zenith_angle = missing
self % sensor_azimuth_angle = missing
self % solar_zenith_angle = missing
self % solar_azimuth_angle = missing
self % cloudtopp = 500.0
self % cloudfrac = 0.0
self % final_cost = missing
self % retrievecloud = .false.
self % mwscatt = .false.
self % mwscatt_totalice = .false.

if (allocated(self % yobs))           deallocate(self % yobs)
if (allocated(self % channels_used))  deallocate(self % channels_used)
if (allocated(self % channels_all))   deallocate(self % channels_all)
if (allocated(self % emiss))          deallocate(self % emiss)
if (allocated(self % output_profile)) deallocate(self % output_profile)
if (allocated(self % output_BT))      deallocate(self % output_BT)
if (allocated(self % calc_emiss))     deallocate(self % calc_emiss)

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
