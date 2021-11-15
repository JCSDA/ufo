! (C) copyright 2020 Met Office
!
! this software is licensed under the terms of the apache licence version 2.0
! which can be obtained at http://www.apache.org/licenses/license-2.0.

!> Fortran module constants used throughout the rttovonedvarcheck filter

module ufo_rttovonedvarcheck_constants_mod

use kinds
use rttov_const, only: surftype_land, surftype_sea, surftype_seaice

implicit none
private

!-----------------------------------------------------------------------------
! 1. miscellaneous definitions
!-----------------------------------------------------------------------------
integer, parameter, public :: max_string = 200 !< maximum string length

! RTTOV values for surface land, sea and seaice
integer, parameter, public :: &
  RTLand = surftype_land, & !< integer for land surface type
  RTSea  = surftype_sea, &  !< integer for sea surface type
  RTIce  = surftype_seaice  !< integer for seaice surface type

!-----------------------------------------------------------------------------
! 2. Physical Constants
!-----------------------------------------------------------------------------

real(kind_real), parameter, public :: &
  MaxTemperature = 340.0_kind_real,   & !< Maximum temperature ( K )
  MinTemperature =  70.0_kind_real      !< Minimum temperature ( K )

end module ufo_rttovonedvarcheck_constants_mod
