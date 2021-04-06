! (C) copyright 2020 Met Office
!
! this software is licensed under the terms of the apache licence version 2.0
! which can be obtained at http://www.apache.org/licenses/license-2.0.

!> Fortran module constants used throughout the rttovonedvarcheck filter

module ufo_rttovonedvarcheck_constants_mod

use kinds

implicit none
private

!-----------------------------------------------------------------------------
! 1. miscellaneous definitions
!-----------------------------------------------------------------------------
integer, parameter, public :: max_string = 200 !< maximum string length

! RTTOV values for surface land, sea and seaice
integer, parameter, public :: &
  RTLand = 0, & !< integer for land surface type
  RTSea  = 1, & !< integer for sea surface type
  RTIce  = 2    !< integer for seaice surface type

!-----------------------------------------------------------------------------
! 2. Physical Constants
!-----------------------------------------------------------------------------

real(kind_real), parameter, public :: &
  MaxTemperature  =   340.0_kind_real,    & !< Maximum temperature ( K )
  MinTemperature  =    70.0_kind_real       !< Minimum temperature ( K )

!-----------------------------------------------------------------------------
! 3. Information for emissivity retrieval
! This is not currently available but has been left in for future development
!-----------------------------------------------------------------------------

!< Mapping for each of the 20 ATOVS channels (1-15 AMSU-A; 16-20 AMSU-B)
integer, parameter, public :: EmissElements(20) = &
  (/ 1,2,3,3,3,3,3,3,3,3,3,3,3,3,4, &  ! AMSU-A mapping
     4,5,5,5,5                      /) ! AMSU-B mapping

!< ATOVS Channel numbers for each of the 5 emissivity values
integer, parameter, public :: EmissMap(5) = (/ 1, 2, 3, 16, 17 /)

end module ufo_rttovonedvarcheck_constants_mod
