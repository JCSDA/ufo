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
! 2. 1d-var profile elements
!-----------------------------------------------------------------------------

! define id codes for 1d-var retrieval fields.
! a list of these fieldtype codes is always present in the header of the bmatrix
! file and it's that list which decides the form of the retrieval vector.
!
! new definitions should be made in conjunction with the profileinfo_type
! structure found in ufo_rttovonedvarcheck_profindex_mod.F90.

integer, parameter, public :: nfieldtypes = 19 !< number of fieldtypes
integer, parameter, public :: &
  fieldtype_t          =  1, &   !< temperature
  fieldtype_q          =  2, &   !< specific humidity profile
  fieldtype_t2         =  3, &   !< surface air temperature
  fieldtype_q2         =  4, &   !< surface spec humidity
  fieldtype_tstar      =  5, &   !< surface skin temperature
  fieldtype_pstar      =  6, &   !< surface pressure
  fieldtype_o3total    =  7, &   !< total column ozone - not currently setup
  fieldtype_not_used   =  8, &   !< not currently in use - not currently setup
  fieldtype_ql         =  9, &   !< liquid water profile - not currently setup
  fieldtype_qt         = 10, &   !< total water profile
  fieldtype_windspeed  = 11, &   !< surface wind speed
  fieldtype_o3profile  = 12, &   !< ozone - not currently setup
  fieldtype_lwp        = 13, &   !< liquid water path - not currently setup
  fieldtype_mwemiss    = 14, &   !< microwave emissivity - not currently setup
  fieldtype_qi         = 15, &   !< ice profile - not currently setup
  fieldtype_cloudtopp  = 16, &   !< single-level cloud top pressure
  fieldtype_cloudfrac  = 17, &   !< effective cloud fraction
  fieldtype_emisspc    = 18, &   !< emissivity prinipal components - not currently setup
  fieldtype_cf         = 19      !< cloud fraction profile - not currently setup

character(len=*), parameter, public :: fieldtype_text(nfieldtypes) = &
  (/ 't                   ', &
     'q water vapour      ', &
     't2                  ', &
     'q2                  ', &
     'tstar               ', &
     'pstar               ', &
     'ozone (total column)', &
     '[unused field type] ', &
     'q liquid            ', &
     'q total             ', &
     'wind speed          ', &
     'ozone (profile)     ', &
     'liquid water path   ', &
     'microwave emissivity', &
     'q ice               ', &
     'cloud top pressure  ', &
     'cloud fraction      ', &
     'emissivity pcs      ', &
     'cloud fraction prof ' /)

!-----------------------------------------------------------------------------
! 3. Physical Constants
!-----------------------------------------------------------------------------

real(kind_real), parameter, public :: &
  MaxTotalOzone   =   650.0,    & !< Maximum total ozone ( Dobson units )
  MinTotalOzone   =    70.0,    & !< Minimum total ozone ( Dobson units )
  MaxSurfaceP     =  1200.0,    & !< Maximum surface pressure ( hPa )
  MinSurfaceP     =   300.0,    & !< Minimum surface pressure ( hPa )
  Min_q           =     3.0E-6, & !< Minimum humidity ( kg / kg )
  MaxTemperature  =   340.0,    & !< Maximum temperature ( K )
  MinTemperature  =    70.0,    & !< Minimum temperature ( K )
  IceShelfLimit   =   -72.0,    & !< Assumed limit of SH seaice
  WetLevelLid     =   115.0,    & !< Uppermost wet pressure level
  MWCloudLevelLid =   310.0

!-----------------------------------------------------------------------------
! 4. Information for emissivity retrieval
! This is not currently available but has been left in for future development
!-----------------------------------------------------------------------------

!< Mapping for each of the 20 ATOVS channels (1-15 AMSU-A; 16-20 AMSU-B)
integer, parameter, public :: EmissElements(20) = &
  (/ 1,2,3,3,3,3,3,3,3,3,3,3,3,3,4, &  ! AMSU-A mapping
     4,5,5,5,5                      /) ! AMSU-B mapping

!< ATOVS Channel numbers for each of the 5 emissivity values
integer, parameter, public :: EmissMap(5) = (/ 1, 2, 3, 16, 17 /)

end module ufo_rttovonedvarcheck_constants_mod
