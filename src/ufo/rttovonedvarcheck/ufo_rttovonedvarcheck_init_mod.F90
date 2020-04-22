! (C) Copyright 2020 Met Office
!
! This software is licensed under the terms of the apache licence version 2.0
! which can be obtained at http://www.apache.org/licenses/license-2.0.

!> Fortran module to provide code shared between nonlinear and tlm/adm radiance calculations

module ufo_rttovonedvarcheck_init_mod

use iso_c_binding
use kinds
use ufo_geovals_mod
use ufo_rttovonedvarcheck_utils_mod
use ufo_rttovonedvarcheck_bmatrix_mod, only: &
        bmatrix_type

implicit none

private

! subroutines - public
public ufo_rttovonedvarcheck_setup
public ufo_rttovonedvarcheck_InitObInfo
public ufo_rttovonedvarcheck_profile_setup

! Subroutines - private
private ufo_rttovonedvarcheck_InitProfInfo

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

subroutine ufo_rttovonedvarcheck_profile_setup( &
  bmatrix,  & ! in
  profinfo)   ! out

!-------------------------------------------------------------------------------
! use b matrix elements to define a 1d-var profile vector. the profinfo data
! structure contains variables specifying the locations of all possible fields.
!
! in order to allow as much flexibility as possible, the profile vector is
! undefined at the beginning of processing so that different combinations of
! fields can be used. the essential criterion for deciding what should be
! present is that we have suitable background error covariances and therefore
! the form of the bmatrix can be used to construct the 1d-var profile. when the
! matrix is read in, we keep a record of the fields that are present, and we use
! that record in here.
!-------------------------------------------------------------------------------

implicit none

! subroutine arguments:
type(bmatrix_type),     intent(in)   :: bmatrix   !background error covariances
type(profileinfo_type), intent(out)  :: profinfo  !profile information

! local constants:
character(len=*), parameter :: routinename = "ufo_rttovonedvarcheck_profile_setup"

! local variables:
integer :: i,j
integer :: firstelement
integer :: lastelement
integer :: nelements

! Initialise to defaults
call ufo_rttovonedvarcheck_InitProfInfo(profinfo)

!loop through fields in bmatrix. the number of elements that the field is composed of
!is held in bmatrix % fields(i,2).
nelements = 0

do j = 1, bmatrix % nfields
  firstelement = nelements + 1
  lastelement  = nelements + bmatrix % fields(j,2)

  !assign start and end points. if the field wasn't found then assign a value of
  !zero, which indicates absence.

  select case( bmatrix % fields(j,1) )

   !----------
   !atmosphere (set start and end points for multi-level fields)
   !----------

    case( fieldtype_t )
      profinfo % t(1)         = firstelement
      profinfo % t(2)         = lastelement

    case( fieldtype_q )
      profinfo % q(1)         = firstelement
      profinfo % q(2)         = lastelement

    case( fieldtype_ql )
      profinfo % ql(1)        = firstelement
      profinfo % ql(2)        = lastelement

    case( fieldtype_qi )
      profinfo % qi(1)        = firstelement
      profinfo % qi(2)        = lastelement

    case( fieldtype_cf )
      profinfo % cf(1)        = firstelement
      profinfo % cf(2)        = lastelement

    case( fieldtype_qt )
      profinfo % qt(1)        = firstelement
      profinfo % qt(2)        = lastelement

    case( fieldtype_o3profile )
      profinfo % o3profile(1) = firstelement
      profinfo % o3profile(2) = lastelement
      if ( ozonesource /= source_bg .and. profinfo % o3profile(1) > 0 ) then
        write(*,*) 'o3 profile retrieval requested but no cx % ozone available'
        write(*,*) 'o3 profile retrieval may result in high rejection rate'
      end if

    case( fieldtype_o3total )
      profinfo % o3total      = firstelement

    case( fieldtype_lwp )
      profinfo % lwp          = firstelement

   !-------
   !surface
   !-------

    case( fieldtype_t2 )
      profinfo % t2         = firstelement

    case( fieldtype_q2 )
      profinfo % q2         = firstelement

    case( fieldtype_tstar )
      profinfo % tstar      = firstelement

    case( fieldtype_pstar )
      profinfo % pstar      = firstelement

    case( fieldtype_windspeed )
      profinfo % windspeed  = firstelement

    case( fieldtype_mwemiss )
      profinfo % mwemiss(1) = firstelement
      profinfo % mwemiss(2) = lastelement

    case( fieldtype_emisspc )
      profinfo % emisspc(1) = firstelement
      profinfo % emisspc(2) = lastelement

   !------------------------------------
   !cloud (single-level grey cloud only)
   !------------------------------------

    case( fieldtype_cloudtopp )
      profinfo % cloudtopp         = firstelement

    case( fieldtype_cloudfrac )
      profinfo % cloudfrac         = firstelement

    case( fieldtype_not_used ) ! currently unused
      continue

    case default
      write(*,*) 'invalid field type in b matrix file: ',j
      cycle

  end select

  if ( firstelement /= 0 ) nelements = nelements + bmatrix % fields(j,2)

end do

profinfo % nprofelements = nelements

end subroutine ufo_rttovonedvarcheck_profile_setup

!-------------------------------------------------------------------------------

subroutine ufo_rttovonedvarcheck_InitProfInfo( profinfo ) ! out

!-------------------------------------------------------------------------------
! initialize profile indices to zero. the profinfo structure is used to store
! the locations of fields in the retrieval vector. zero will indicate absence of
! a field.
!-------------------------------------------------------------------------------

implicit none

! subroutine arguments:
type(profileinfo_type), intent(out) :: profinfo

character(len=*), parameter :: routinename = "ufo_rttovonedvarcheck_InitProfInfo"

!-------------------------------------------------------------------------------

profinfo % nprofelements = 0
profinfo % t = 0
profinfo % q = 0
profinfo % ql = 0
profinfo % qt = 0
profinfo % qi = 0
profinfo % cf = 0
profinfo % o3total = 0
profinfo % o3profile = 0
profinfo % t2 = 0
profinfo % q2 = 0
profinfo % rh2 = 0
profinfo % tstar = 0
profinfo % pstar = 0
profinfo % windspeed = 0
profinfo % t70hpa = 0
profinfo % t700hpa = 0
profinfo % t950hpa = 0
profinfo % t1000hpa = 0
profinfo % qsurf = 0
profinfo % lwp = 0
profinfo % mwemiss = 0
profinfo % cloudtopp = 0
profinfo % cloudfrac = 0
profinfo % emisspc = 0

end subroutine ufo_rttovonedvarcheck_InitProfInfo

!------------------------------------------------------------------------------

end module ufo_rttovonedvarcheck_init_mod
