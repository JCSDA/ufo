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
  logical                          :: FullDiagnostics
  integer                          :: Max1DVarIterations
  integer                          :: JConvergenceOption
  integer                          :: IterNumForLWPCheck
  real(kind_real)                  :: ConvergenceFactor
  real(kind_real)                  :: Cost_ConvergenceFactor
  real(kind_real)                  :: MaxMLIterations
end type ufo_rttovonedvarcheck

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

contains

!===============================================================================
! Subroutines
!===============================================================================

subroutine ufo_rttovonedvarcheck_setup(self, channels)

implicit none

! subroutine arguments
type(ufo_rttovonedvarcheck), intent(inout) :: self
integer(c_int), intent(in)                 :: channels(:)

! local variables
character(len=max_string) :: tmp
character(len=:), allocatable :: str
character(kind=c_char, len=max_string), allocatable :: char_array(:)
integer(c_size_t),parameter :: csize = max_string
integer :: ii

! Setup core paths and names
self % qcname = "rttovonedvarcheck"
call self % conf % get_or_die("BMatrix",str)
self % b_matrix_path = str
call self % conf % get_or_die("ModName",str)
self % forward_mod_name = str
call self % conf % get_or_die("nlevels",self % nlevels)

! Variables for profile (x,xb)
self % nmvars = self % conf % get_size("model_variables")
allocate(self % model_variables(self % nmvars))
call self % conf % get_or_die("model_variables", csize, char_array)
self % model_variables(1:self % nmvars) = char_array

! Satellite channels
self % nchans = size(channels)
allocate(self % channels(self % nchans))
self % channels(:) = channels(:)
write(*,*) "nchans setup = ",self%nchans
write(*,*) "channels setup = ",self%channels

! Set defaults for 1D-var
self % rtype = "diagonal"
self % qtotal = .false.
self % RTTOV_mwscattSwitch = .false.
self % use_totalice = .false.
self % UseMLMinimization = .false.
self % UseJforConvergence = .false.
self % FullDiagnostics = .false.
self % Max1DVarIterations = 7
self % JConvergenceOption = 1
self % IterNumForLWPCheck = 2
self % ConvergenceFactor = 0.40
self % Cost_ConvergenceFactor = 0.01
self % MaxMLIterations = 7

! R matrix type to use
if (self % conf % has("rtype")) then
  call self % conf % get_or_die("rtype",str)
  self % rtype = str
end if

! Flag for total humidity
if (self % conf % has("qtotal")) then
  call self % conf % get_or_die("qtotal", self % qtotal)
end if

! Flag for RTTOV MW scatt
if (self % conf % has("RTTOV_mwscattSwitch")) then
  call self % conf % get_or_die("RTTOV_mwscattSwitch", self % RTTOV_mwscattSwitch)
end if

! Flag for use of total ice in RTTOV MW scatt
if (self % conf % has("use_totalice")) then
  call self % conf % get_or_die("use_totalice", self % use_totalice)
end if

! Flag to turn on marquardt-levenberg minimiser
if (self % conf % has("UseMLMinimization")) then
  call self % conf % get_or_die("UseMLMinimization", self % UseMLMinimization)
end if

! Flag to Use J for convergence
if (self % conf % has("UseJforConvergence")) then
  call self % conf % get_or_die("UseJforConvergence", self % UseJforConvergence)
end if

! Flag to turn on full diagnostics
if (self % conf % has("FullDiagnostics")) then
  call self % conf % get_or_die("FullDiagnostics", self % FullDiagnostics)
end if

! Flag to specify if delta_x has to be negative for converg. to be true
if (self % conf % has("Max1DVarIterations")) then
  call self % conf % get_or_die("Max1DVarIterations", self % Max1DVarIterations)
end if

! Flag to specify if delta_x has to be negative for converg. to be true
if (self % conf % has("JConvergenceOption")) then
  call self % conf % get_or_die("JConvergenceOption", self % JConvergenceOption)
end if

! Choose which iteration to start checking LWP
if (self % conf % has("IterNumForLWPCheck")) then
  call self % conf % get_or_die("IterNumForLWPCheck", self % IterNumForLWPCheck)
end if

! Flag to specify if delta_x has to be negative for converg. to be true
if (self % conf % has("ConvergenceFactor")) then
  call self % conf % get_or_die("ConvergenceFactor", self % ConvergenceFactor)
end if

! Flag to specify if delta_x has to be negative for converg. to be true
if (self % conf % has("Cost_ConvergenceFactor")) then
  call self % conf % get_or_die("Cost_ConvergenceFactor", self % Cost_ConvergenceFactor)
end if

! Print self
if (self % FullDiagnostics) then
  write(*,*) "qcname = ",trim(self % qcname)
  write(*,*) "rtype = ",trim(self % rtype)
  write(*,*) "b_matrix_path = ",trim(self % b_matrix_path)
  write(*,*) "forward_mod_name = ",trim(self % forward_mod_name)
  write(*,*) "model_variables = "
  do ii = 1, self % nmvars
    write(*,*) trim(self % model_variables(ii))," "
  end do
  write(*,*) "nlevels = ",self %  nlevels
  write(*,*) "nmvars = ",self % nmvars
  write(*,*) "nchans = ",self % nchans
  write(*,*) "channels(:) = ",self % channels(:)
  write(*,*) "qtotal = ",self % qtotal
  write(*,*) "RTTOV_mwscattSwitch = ",self % RTTOV_mwscattSwitch
  write(*,*) "use_totalice = ",self % use_totalice
  write(*,*) "UseMLMinimization = ",self % UseMLMinimization
  write(*,*) "UseJforConvergence = ",self % UseJforConvergence
  write(*,*) "FullDiagnostics = ",self % FullDiagnostics
  write(*,*) "Max1DVarIterations = ",self % Max1DVarIterations
  write(*,*) "JConvergenceOption = ",self % JConvergenceOption
  write(*,*) "IterNumForLWPCheck = ",self % IterNumForLWPCheck
  write(*,*) "ConvergenceFactor = ",self % ConvergenceFactor
  write(*,*) "Cost_ConvergenceFactor = ",self % Cost_ConvergenceFactor
  write(*,*) "MaxMLIterations = ",self % MaxMLIterations
end if

end subroutine

!-------------------------------------------------------------------------------

subroutine ufo_rttovonedvarcheck_InitObInfo(ob_info, & ! out
                                            nchans)    ! in

implicit none

! subroutine arguments:
type(obinfo_type), intent(out) :: ob_info
integer, intent(in)  :: nchans

character(len=*), parameter :: routinename = "ufo_rttovonedvarcheck_InitObInfo"

ob_info % nlocs = 1
ob_info % latitude = 0.0
ob_info % longitude = 0.0
ob_info % elevation = 0.0
ob_info % sensor_zenith_angle = 0.0
ob_info % sensor_azimuth_angle = 0.0
ob_info % solar_zenith_angle = 0.0
ob_info % solar_azimuth_angle = 0.0

! Moved this to obs loop to allow channel selection
!allocate(ob_info%yobs(nchans))
!ob_info % yobs = 0.0

end subroutine ufo_rttovonedvarcheck_InitObInfo

!-------------------------------------------------------------------------------

end module ufo_rttovonedvarcheck_utils_mod
