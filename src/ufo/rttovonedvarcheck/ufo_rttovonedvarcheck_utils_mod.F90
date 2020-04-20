! (C) copyright 2020 Met Office
!
! this software is licensed under the terms of the apache licence version 2.0
! which can be obtained at http://www.apache.org/licenses/license-2.0.

!> fortran module to provide code shared between nonlinear and tlm/adm radiance calculations

module ufo_rttovonedvarcheck_utils_mod

use fckit_configuration_module, only: fckit_configuration
use iso_c_binding
use kinds

implicit none

! defaults are public
! defined types : ufo_rttovonedvarcheck, profileinfo_type, obinfo_type
!

integer, parameter :: max_string = 200

!===============================================================================
! type definitions
!===============================================================================

!---------------
! 1. 1DVar type
!---------------

type :: ufo_rttovonedvarcheck
  character(len=max_string)        :: qcname
  character(len=max_string)        :: rtype
  character(len=max_string)        :: b_matrix_path
  character(len=max_string)        :: forward_mod_name
  character(len=max_string), allocatable :: model_variables(:)
  type(c_ptr)                      :: obsdb
  type(fckit_configuration)        :: conf
  integer                          :: nlevels   ! number 1D-Var model levels
  integer                          :: nmvars
  integer                          :: nchans
  integer(c_int), allocatable      :: channels(:)
  logical                          :: qtotal
  logical                          :: RTTOV_mwscattSwitch
  logical                          :: use_totalice
  logical                          :: UseMLMinimization
  logical                          :: UseJforConvergence
  integer                          :: Max1DVarIterations
  integer                          :: JConvergenceOption
  integer                          :: IterNumForLWPCheck
  real(kind_real)                  :: ConvergenceFactor
  real(kind_real)                  :: Cost_ConvergenceFactor
  real(kind_real)                  :: Mqstart
  real(kind_real)                  :: Mqstep
  real(kind_real)                  :: MaxMLIterations
end type ufo_rttovonedvarcheck

!----------------------
! 2. Profile Information
!----------------------

!The profile variables locate a particular field, or element thereof, in the
!1d-var profile vector. absence of a field will be designated with zero. note
!that retrieval fields also require b matrix fieldtype definitions (see
!section 13 below).

!also, remember to update ufo_rttovonedvarcheck_initprofinfo and ufo_rttovonedvarcheck_mapprofiletob
!if adding fields to this structure.

type profileinfo_type

 !general

  integer :: nprofelements

 !atmosphere (locate start and end points of relevant fields)

  integer :: t(2)           ! temperature profile
  integer :: q(2)           ! water vapour profile (specific humidity)
  integer :: ql(2)          ! liquid water profile
  integer :: qt(2)          ! total water profile
  integer :: qi(2)          ! frozen ice  profile
  integer :: cf(2)          ! cloud fraction profile
  integer :: o3total        ! total column ozone
  integer :: o3profile(2)   ! ozone profile
  integer :: lwp            ! liquid water path

 !surface

  integer :: t2             ! screen temperature
  integer :: q2             ! screen specific humidity
  integer :: rh2            ! screen relative humidity
  integer :: tstar          ! skin temperature
  integer :: pstar          ! surface pressure
  integer :: mwemiss(2)     ! microwave emissivity
  integer :: emisspc(2)     ! emissivity principal components
  integer :: windspeed      ! surface windspeed

 !cloud (single-level grey cloud only)
  integer :: cloudtopp      ! single-level cloud top pressure
  integer :: cloudfrac      ! effective cloud fraction

 !other

  integer :: t70hpa         ! temperature at 70hpa
  integer :: t700hpa        ! temperature at 700hpa
  integer :: t950hpa        ! temperature at 950hpa
  integer :: t1000hpa       ! temperature at 1000hpa
  integer :: qsurf          ! surface specific humidity

end type

!---------------------------------------------------------
! 3. container for information about a single observation
!---------------------------------------------------------

type obinfo_type

  character(len=max_string) :: forward_mod_name
  integer         :: nlocs
  real(kind_real) :: latitude
  real(kind_real) :: longitude
  real(kind_real) :: elevation
  real(kind_real) :: sensor_zenith_angle
  real(kind_real) :: sensor_azimuth_angle
  real(kind_real) :: solar_zenith_angle
  real(kind_real) :: solar_azimuth_angle
  real(kind_real),allocatable :: yobs(:)

end type

!===============================================================================
! variables definitions
!===============================================================================

!---------------------------
! 1. 1d-var profile elements
!---------------------------

! define id codes for 1d-var retrieval fields.
! a list of these fieldtype codes is always present in the header of the bmatrix
! file and it's that list which decides the form of the retrieval vector.
!
! new definitions should be made in conjunction with the profileinfo_type
! structure found in section 3 above.

integer, parameter :: nfieldtypes = 19
integer, parameter :: &
  fieldtype_t          =  1, &   ! temperature
  fieldtype_q          =  2, &   ! specific humidity profile
  fieldtype_t2         =  3, &   ! surface air temperature
  fieldtype_q2         =  4, &   ! surface spec humidity
  fieldtype_tstar      =  5, &   ! surface skin temperature
  fieldtype_pstar      =  6, &   ! surface pressure
  fieldtype_o3total    =  7, &   ! total column ozone
  fieldtype_not_used   =  8, &   ! not currently in use
  fieldtype_ql         =  9, &   ! liquid water profile
  fieldtype_qt         = 10, &   ! total water profile
  fieldtype_windspeed  = 11, &   ! surface wind speed
  fieldtype_o3profile  = 12, &   ! ozone
  fieldtype_lwp        = 13, &   ! liquid water path
  fieldtype_mwemiss    = 14, &   ! microwave emissivity
  fieldtype_qi         = 15, &   ! ice profile
  fieldtype_cloudtopp  = 16, &   ! single-level cloud top pressure
  fieldtype_cloudfrac  = 17, &   ! effective cloud fraction
  fieldtype_emisspc    = 18, &   ! emissivity prinipal components
  fieldtype_cf         = 19      ! cloud fraction profile

character(len=*), parameter :: fieldtype_text(nfieldtypes) = &
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

!------------------------------
! 2. miscellaneous definitions
!------------------------------

!data source codes, used to assign background ozone information. the numbers
!are arbitrary.

integer, parameter :: &
  source_bg        = 1, &
  source_estimate  = 2, &
  source_product   = 3, &
  source_reference = 4

! for now default ozone profile to be provided by model
integer :: ozonesource = source_bg

!===============================================================================
!Physical Constants
!===============================================================================

real(kind_real), parameter :: &
  IR_Emiss       =      0.98,   & ! surf emiss to use for IR channels
  MaxTotalOzone  =    650.0,    & ! ( Dobson units )
  MinTotalOzone  =     70.0,    & ! ( Dobson units )
  MaxSurfaceP    =   1200.0,    & ! ( hPa )
  MinSurfaceP    =    300.0,    & ! ( hPa )
  Min_q          =      3.0E-6, & ! ( kg / kg )
  MaxTemperature =    340.0,    & ! ( K )
  MinTemperature =     70.0,    & ! ( K )
  IceShelfLimit  =    -72.0,    & ! assumed limit of SH seaice
  WetLevelLid    =    115.0,    & ! uppermost wet pressure level
  MWCloudLevelLid=    310.0,    & ! uppermost clw pressure lvl for mwcalcs(hpa)
  g              =      9.80665   ! gravity at the surface (m/s2)

end module ufo_rttovonedvarcheck_utils_mod
