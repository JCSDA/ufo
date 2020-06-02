! (C) copyright 2020 Met Office
!
! this software is licensed under the terms of the apache licence version 2.0
! which can be obtained at http://www.apache.org/licenses/license-2.0.

!> fortran module to provide code shared between nonlinear and tlm/adm radiance calculations

module ufo_rttovonedvarcheck_setup_mod

use iso_c_binding
use kinds
use ufo_geovals_mod
use ufo_rttovonedvarcheck_utils_mod
use ufo_rttovonedvarcheck_profindex_mod, only: profindex_type
use missing_values_mod

implicit none
private

public ufo_rttovonedvarcheck_setup
public ufo_rttovonedvarcheck_check_geovals
public ufo_rttovonedvarcheck_InitObInfo
public ufo_rttovonedvarcheck_InitEmiss
public ufo_rttovonedvarcheck_DeleteObInfo

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
self % EmissLandDefault = 0.95    ! default land surface emissivity
self % EmissSeaIceDefault = 0.92  ! default seaice surface emissivity
self % MwEmiss = .false.
self % IREmiss = .false.

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

! Default emissivity value to use over land
if (self % conf % has("EmissLandDefault")) then
  call self % conf % get_or_die("EmissLandDefault", self % EmissLandDefault)
end if

! Default emissivity value to use over seaice
if (self % conf % has("EmissSeaIceDefault")) then
  call self % conf % get_or_die("EmissSeaIceDefault", self % EmissSeaIceDefault)
end if

! Flag to specify if microwave emissivity is used
if (self % conf % has("MwEmiss")) then
  call self % conf % get_or_die("MwEmiss", self % MwEmiss)
end if

! Flag to specify if infrared emissivity is used
if (self % conf % has("IREmiss")) then
  call self % conf % get_or_die("IREmiss", self % IREmiss)
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
  write(*,*) "EmissLandDefault = ",self % EmissLandDefault
  write(*,*) "EmissSeaIceDefault = ",self % EmissSeaIceDefault
  write(*,*) "MwEmiss = ",self % MwEmiss
  write(*,*) "IREmiss = ",self % IREmiss
end if

end subroutine

!-------------------------------------------------------------------------------

subroutine ufo_rttovonedvarcheck_check_geovals(geovals, profindex)

! Heritage: Ops_SatRad_SetUpRTprofBg_RTTOV12.f90

use ufo_rttovonedvarcheck_minimize_utils_mod
implicit none

! subroutine arguments:
type(ufo_geovals), intent(inout) :: geovals
type(profindex_type), intent(in) :: profindex

character(len=*), parameter  :: routinename = "ufo_rttovonedvarcheck_check_geovals"
character(len=max_string)    :: varname
type(ufo_geoval), pointer    :: geoval
integer                      :: gv_index, i     ! counters
integer                      :: nlevels
real(kind_real), allocatable :: temperature(:)  ! temperature (K)
real(kind_real), allocatable :: pressure(:)     ! pressure (Pa)
real(kind_real), allocatable :: humidity_total(:)
real(kind_real), allocatable :: q(:)            ! specific humidity (kg/kg)
real(kind_real), allocatable :: ql(:)
real(kind_real), allocatable :: qi(:)

write(*,*) routinename, " : started"

!-------------------------
! Specific humidity total
!-------------------------

if (profindex % qt(1) > 0) then

  nlevels = profindex % qt(2) - profindex % qt(1) + 1
  allocate(temperature(nlevels))
  allocate(pressure(nlevels))
  allocate(humidity_total(nlevels))
  allocate(q(nlevels))
  allocate(ql(nlevels))
  allocate(qi(nlevels))

  ! Get temperature and pressure from geovals
  call ufo_geovals_get_var(geovals, "air_temperature", geoval)
  temperature(:) = geoval%vals(:, 1) ! K
  call ufo_geovals_get_var(geovals, "air_pressure", geoval)
  pressure(:) = geoval%vals(:, 1)    ! Pa

  ! Get humidity data from geovals
  humidity_total(:) = 0.0
  call ufo_geovals_get_var(geovals, "specific_humidity", geoval)
  humidity_total(:) = humidity_total(:) + geoval%vals(:, 1)
  call ufo_geovals_get_var(geovals, "mass_content_of_cloud_liquid_water_in_atmosphere_layer", geoval)
  humidity_total(:) = humidity_total(:) + geoval%vals(:, 1)

  ! Split qtotal to q(water_vapour), q(liquid), q(ice)
  call ufo_rttovonedvarcheck_Qsplit (1,      & ! in
                          temperature(:),    & ! in
                          pressure(:),       & ! in
                          nlevels,           & ! in
                          humidity_total(:), & ! in
                          q(:),              & ! out
                          ql(:),             & ! out
                          qi(:))               ! out

  ! Assign values to geovals
  varname = "specific_humidity"  ! kg/kg
  gv_index = 0
  do i=1,geovals%nvar
    if (varname == trim(geovals%variables(i))) gv_index = i
  end do
  geovals%geovals(gv_index) % vals(:,1) = q(:)

  varname = "mass_content_of_cloud_liquid_water_in_atmosphere_layer"  ! kg/kg
  gv_index = 0
  do i=1,geovals%nvar
    if (varname == trim(geovals%variables(i))) gv_index = i
  end do
  geovals%geovals(gv_index) % vals(:,1) = ql(:)

  varname = "mass_content_of_cloud_ice_in_atmosphere_layer"  ! kg/kg
  gv_index = 0
  do i=1,geovals%nvar
    if (varname == trim(geovals%variables(i))) gv_index = i
  end do
  geovals%geovals(gv_index)%vals(:,1) = qi(:)

  deallocate(temperature)
  deallocate(pressure)
  deallocate(humidity_total)
  deallocate(q)
  deallocate(ql)
  deallocate(qi)

end if

! Tidy up
if (allocated(temperature))    deallocate(temperature)
if (allocated(pressure))       deallocate(pressure)
if (allocated(humidity_total)) deallocate(humidity_total)
if (allocated(q))              deallocate(q)
if (allocated(ql))             deallocate(ql)
if (allocated(qi))             deallocate(qi)

write(*,*) routinename, " : ended"

end subroutine ufo_rttovonedvarcheck_check_geovals

!-------------------------------------------------------------------------------

subroutine ufo_rttovonedvarcheck_InitObInfo(ob_info, & ! out
                                            nchans)    ! in

implicit none

! subroutine arguments:
type(obinfo_type), intent(out) :: ob_info
integer, intent(in)  :: nchans

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

allocate(ob_info % yobs(nchans))
allocate(ob_info % emiss(nchans))
allocate(ob_info % calc_emiss(nchans))

ob_info % yobs(nchans) = missing
ob_info % emiss(:) = 0.0
ob_info % calc_emiss(:) = .false.

end subroutine ufo_rttovonedvarcheck_InitObInfo

!-------------------------------------------------------------------------------

subroutine ufo_rttovonedvarcheck_InitEmiss(self, geovals, ob_info)

implicit none

! subroutine arguments:
type(ufo_rttovonedvarcheck), intent(in) :: self
type(ufo_geovals), intent(inout)        :: geovals
type(obinfo_type), intent(inout)        :: ob_info

type(ufo_geoval), pointer    :: geoval

! ----------------------------------
! 1.0 Get surface type from geovals
!-----------------------------------
call ufo_geovals_get_var(geovals, "surface_type", geoval)
ob_info % surface_type = geoval%vals(1, 1)

!-------------
! 2.1 Defaults
!-------------

if (self % MwEmiss) then

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

end if

end subroutine ufo_rttovonedvarcheck_InitEmiss

!-------------------------------------------------------------------------------

subroutine ufo_rttovonedvarcheck_DeleteObInfo(ob_info)    ! inout

implicit none

! subroutine arguments:
type(obinfo_type), intent(inout) :: ob_info

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

if (allocated(ob_info % yobs))       deallocate(ob_info % yobs)
if (allocated(ob_info % emiss))      deallocate(ob_info % emiss)
if (allocated(ob_info % calc_emiss)) deallocate(ob_info % calc_emiss)

end subroutine ufo_rttovonedvarcheck_DeleteObInfo

!-------------------------------------------------------------------------------

end module ufo_rttovonedvarcheck_setup_mod
