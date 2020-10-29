! (C) copyright 2020 Met Office
!
! this software is licensed under the terms of the apache licence version 2.0
! which can be obtained at http://www.apache.org/licenses/license-2.0.

!> Fortran module containing main type, setup and utilities for the
!! main rttovonedvarcheck object

module ufo_rttovonedvarcheck_utils_mod

use fckit_configuration_module, only: fckit_configuration
use iso_c_binding
use kinds
use ufo_rttovonedvarcheck_constants_mod

implicit none
private

public ufo_rttovonedvarcheck_setup
public ufo_rttovonedvarcheck_iogetfreeunit

!===============================================================================
! type definitions
!===============================================================================

!---------------
! 1. 1DVar type
!---------------

type, public :: ufo_rttovonedvarcheck
  character(len=max_string)        :: qcname !< name of the filter
  character(len=max_string)        :: b_matrix_path !< path to the b-matrix file
  character(len=max_string)        :: r_matrix_path !< path to the r-matrix file
  character(len=max_string)        :: forward_mod_name !< forward model name (only RTTOV at the moment)
  character(len=max_string), allocatable :: retrieval_variables(:) !< list of variables which form the 1D-var retrieval vector
  type(c_ptr)                      :: obsdb !< pointer to the observation space
  type(fckit_configuration)        :: conf  !< contents of the yaml file
  integer(c_int)                   :: onedvarflag !< flag uased by the qc manager for a 1D-var check
  integer(c_int)                   :: passflag !< flag uased by the qc manager to flag good data
  integer                          :: nlevels ! number 1D-Var model levels
  integer                          :: nmvars !< number of variables being used in the retrieval
  integer                          :: nchans !< maximum number of channels (channels can be removed by previous qc checks)
  integer(c_int), allocatable      :: channels(:) !< integer list of channels
  logical                          :: qtotal !< flag to enable total humidity retrievals
  logical                          :: RTTOV_mwscattSwitch !< flag to switch on RTTOV-scatt
  logical                          :: RTTOV_usetotalice !< flag for use of total ice in RTTOV MW scatt
  logical                          :: UseMLMinimization !< flag to turn on marquardt-levenberg minimizer
  logical                          :: UseJforConvergence !< flag to Use J for convergence
  logical                          :: FullDiagnostics !< flag to turn on full diagnostics
  logical                          :: pcemiss !< flag gets turned off in emissivity eigen vector file is present
  integer                          :: Max1DVarIterations !< maximum number of iterations
  integer                          :: JConvergenceOption !< integer to select convergence option
  integer                          :: IterNumForLWPCheck !< choose which iteration to start checking LWP
  real(kind_real)                  :: ConvergenceFactor !< 1d-var convergence if using change in profile
  real(kind_real)                  :: Cost_ConvergenceFactor !< 1d-var convergence if using % change in cost
  real(kind_real)                  :: MaxMLIterations !< maximum number of iterations for internal Marquardt-Levenberg loop
  real(kind_real)                  :: EmissLandDefault !< default emissivity value to use over land
  real(kind_real)                  :: EmissSeaIceDefault !< default emissivity value to use over sea ice
  character(len=max_string)        :: EmisEigVecPath !< path to eigen vector file for IR PC emissivity
  character(len=max_string)        :: EmisAtlas !< path to the emissivity atlas for IR PC emissivity
end type ufo_rttovonedvarcheck

contains

!------------------------------------------------------------------------------
!> Setup the defaults for the main rttovonedvarcheck object and read in the 
!! contents of the yaml file.
!!
!! \author Met Office
!!
!! \date 09/06/2020: Created
!!
subroutine ufo_rttovonedvarcheck_setup(self, channels)

implicit none

! subroutine arguments
type(ufo_rttovonedvarcheck), intent(inout) :: self
integer(c_int), intent(in)                 :: channels(:)

! local variables
character(len=max_string) :: tmp
character(len=:), allocatable :: str
character(len=:), allocatable :: str_array(:)

! Setup core paths and names
self % qcname = "rttovonedvarcheck"
call self % conf % get_or_die("BMatrix",str)
self % b_matrix_path = str
call self % conf % get_or_die("RMatrix",str)
self % r_matrix_path = str
call self % conf % get_or_die("ModName",str)
self % forward_mod_name = str
call self % conf % get_or_die("nlevels",self % nlevels)

! Variables for profile (x,xb)
self % nmvars = self % conf % get_size("retrieval variables")
allocate(self % retrieval_variables(self % nmvars))
call self % conf % get_or_die("retrieval variables", str_array)
self % retrieval_variables(1:self % nmvars) = str_array

! Satellite channels
self % nchans = size(channels)
allocate(self % channels(self % nchans))
self % channels(:) = channels(:)

! Set defaults for 1D-var
self % qtotal = .false.
self % RTTOV_mwscattSwitch = .false.
self % RTTOV_usetotalice = .false.
self % UseMLMinimization = .false.
self % UseJforConvergence = .false.
self % FullDiagnostics = .false.
self % pcemiss = .false.
self % Max1DVarIterations = 7
self % JConvergenceOption = 1
self % IterNumForLWPCheck = 2
self % ConvergenceFactor = 0.40
self % Cost_ConvergenceFactor = 0.01
self % MaxMLIterations = 7
self % EmissLandDefault = 0.95    ! default land surface emissivity
self % EmissSeaIceDefault = 0.92  ! default seaice surface emissivity
self % EmisEigVecPath = ""
self % EmisAtlas = ""

! Flag for total humidity
if (self % conf % has("qtotal")) then
  call self % conf % get_or_die("qtotal", self % qtotal)
end if

! Flag for RTTOV MW scatt
if (self % conf % has("RTTOV_mwscattSwitch")) then
  call self % conf % get_or_die("RTTOV_mwscattSwitch", self % RTTOV_mwscattSwitch)
end if

! Flag for use of total ice in RTTOV MW scatt
if (self % conf % has("RTTOV_usetotalice")) then
  call self % conf % get_or_die("RTTOV_usetotalice", self % RTTOV_usetotalice)
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

! maximum number of iterations
if (self % conf % has("Max1DVarIterations")) then
  call self % conf % get_or_die("Max1DVarIterations", self % Max1DVarIterations)
end if

! integer to select convergence option
if (self % conf % has("JConvergenceOption")) then
  call self % conf % get_or_die("JConvergenceOption", self % JConvergenceOption)
end if

! Choose which iteration to start checking LWP
if (self % conf % has("IterNumForLWPCheck")) then
  call self % conf % get_or_die("IterNumForLWPCheck", self % IterNumForLWPCheck)
end if

! Cost threhsold for convergence check
if (self % conf % has("ConvergenceFactor")) then
  call self % conf % get_or_die("ConvergenceFactor", self % ConvergenceFactor)
end if

! Flag to specify if delta_x has to be negative for converg. to be true
if (self % conf % has("Cost_ConvergenceFactor")) then
  call self % conf % get_or_die("Cost_ConvergenceFactor", self % Cost_ConvergenceFactor)
end if

! Maximum number of iterations for internal Marquardt-Levenberg loop
if (self % conf % has("MaxMLIterations")) then
  call self % conf % get_or_die("MaxMLIterations", self % MaxMLIterations)
end if

! Default emissivity value to use over land
if (self % conf % has("EmissLandDefault")) then
  call self % conf % get_or_die("EmissLandDefault", self % EmissLandDefault)
end if

! Default emissivity value to use over seaice
if (self % conf % has("EmissSeaIceDefault")) then
  call self % conf % get_or_die("EmissSeaIceDefault", self % EmissSeaIceDefault)
end if

! Default eigen value path is blank but needs to be present if using PC emiss
if (self % conf % has("EmisEigVecPath")) then
  call self % conf % get_or_die("EmisEigVecPath",str)
  self % EmisEigVecPath = str
  self % pcemiss = .true.
end if

! Default emis atlas path is blank
if (self % conf % has("EmisAtlas")) then
  call self % conf % get_or_die("EmisAtlas",str)
  self % EmisAtlas = str
end if

! Print self
if (self % FullDiagnostics) then
  call ufo_rttovonedvarcheck_print(self)
end if

end subroutine ufo_rttovonedvarcheck_setup

!------------------------------------------------------------------------------
!> Print contents of rttovonedvarcheck object
!!
!! \author Met Office
!!
!! \date 03/08/2020: Created
!!
subroutine ufo_rttovonedvarcheck_print(self)

implicit none

type(ufo_rttovonedvarcheck), intent(in) :: self

integer :: ii

write(*,*) "qcname = ",trim(self % qcname)
write(*,*) "b_matrix_path = ",trim(self % b_matrix_path)
write(*,*) "r_matrix_path = ",trim(self % r_matrix_path)
write(*,*) "forward_mod_name = ",trim(self % forward_mod_name)
write(*,*) "retrieval_variables = "
do ii = 1, self % nmvars
  write(*,*) trim(self % retrieval_variables(ii))," "
end do
write(*,*) "nlevels = ",self %  nlevels
write(*,*) "nmvars = ",self % nmvars
write(*,*) "nchans = ",self % nchans
write(*,*) "channels(:) = ",self % channels(:)
write(*,*) "qtotal = ",self % qtotal
write(*,*) "RTTOV_mwscattSwitch = ",self % RTTOV_mwscattSwitch
write(*,*) "RTTOV_usetotalice = ",self % RTTOV_usetotalice
write(*,*) "UseMLMinimization = ",self % UseMLMinimization
write(*,*) "UseJforConvergence = ",self % UseJforConvergence
write(*,*) "FullDiagnostics = ",self % FullDiagnostics
write(*,*) "Max1DVarIterations = ",self % Max1DVarIterations
write(*,*) "JConvergenceOption = ",self % JConvergenceOption
write(*,*) "IterNumForLWPCheck = ",self % IterNumForLWPCheck
write(*,*) "ConvergenceFactor = ",self % ConvergenceFactor
write(*,*) "Cost_ConvergenceFactor = ",self % Cost_ConvergenceFactor
write(*,*) "MaxMLIterations = ",self % MaxMLIterations
write(*,*) "EmissLandDefault = ",self % EmissLandDefault
write(*,*) "EmissSeaIceDefault = ",self % EmissSeaIceDefault
write(*,*) "EmisEigVecPath = ",self % EmisEigVecPath
write(*,*) "EmisAtlas = ",self % EmisAtlas

end subroutine ufo_rttovonedvarcheck_print

!------------------------------------------------------------------------------
!> Find a free file unit.
!!
!! \author Met Office
!!
!! \date 09/06/2020: Created
!!
subroutine ufo_rttovonedvarcheck_iogetfreeunit(unit)

implicit none

integer, intent(out) :: unit

integer, parameter :: unit_min=10
integer, parameter :: unit_max=1000
logical            :: opened
integer            :: lun
integer            :: newunit

newunit=-1
do lun=unit_min,unit_max
  inquire(unit=lun,opened=opened)
  if (.not. opened) then
      newunit=lun
    exit
  end if
end do
unit=newunit

end subroutine ufo_rttovonedvarcheck_iogetfreeunit

! ------------------------------------------------------------------------------

end module ufo_rttovonedvarcheck_utils_mod
