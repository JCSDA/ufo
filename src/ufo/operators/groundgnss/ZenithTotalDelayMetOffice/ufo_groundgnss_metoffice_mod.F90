!-------------------------------------------------------------------------------
! (C) British Crown Copyright 2020 Met Office
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
!-------------------------------------------------------------------------------

!> Fortran module for ground based gnss Met Office forward operator

module ufo_groundgnss_metoffice_mod

use, intrinsic :: iso_c_binding
use kinds
use ufo_vars_mod
use ufo_geovals_mod
use ufo_geovals_mod_c, only: ufo_geovals_registry
use ufo_basis_mod,     only: ufo_basis
use obsspace_mod
use missing_values_mod
use ufo_groundgnss_ukmo_utils_mod
use fckit_log_module,  only : fckit_log
use ufo_utils_refractivity_calculator, only: ufo_calculate_refractivity

implicit none
public             :: ufo_groundgnss_metoffice
private

integer, parameter :: max_string = 800

  !> Fortran derived type for groundgnss trajectory
type, extends(ufo_basis) :: ufo_groundgnss_MetOffice
  logical :: vert_interp_ops
  logical :: pseudo_ops
  real(kind_real) :: min_temp_grad
  real(kind_real) :: dryRefractivityConstant
  real(kind_real) :: wetRefractivityConstant
  integer :: num_fov                    ! Number of field-of-view sampling points (1 = no averaging)
  real(kind_real) :: fov_radius         ! radius of the sampling points, in km
  contains
    procedure :: setup     => ufo_groundgnss_metoffice_setup
    procedure :: simobs    => ufo_groundgnss_metoffice_simobs
end type ufo_groundgnss_MetOffice

public :: ufo_groundgnss_metoffice_locs_init

contains

! ------------------------------------------------------------------------------
! Get the optional settings for the forward model, and save them in the object
! so that they can be used in the code.
! ------------------------------------------------------------------------------
subroutine ufo_groundgnss_metoffice_setup(self, &
                                         vertInterpOPS, &
                                         pseudoLevels, &
                                         minTempGrad, &
                                         dryRefractivityConstant, &
                                         wetRefractivityConstant, &
                                         numFOV, &
                                         fovRadius)

use fckit_configuration_module, only: fckit_configuration
implicit none
class(ufo_groundgnss_MetOffice), intent(inout) :: self
logical(c_bool), intent(in) :: vertInterpOPS
logical(c_bool), intent(in) :: pseudoLevels
real(c_float), intent(in) :: minTempGrad
real(c_float), intent(in) :: dryRefractivityConstant
real(c_float), intent(in) :: wetRefractivityConstant
integer(c_int), intent(in) :: numFOV
real(c_float), intent(in) :: fovRadius

character(max_string) :: err_msg

self % vert_interp_ops = vertInterpOPS
self % pseudo_ops = pseudoLevels
self % min_temp_grad = minTempGrad
self % dryRefractivityConstant = dryRefractivityConstant
self % wetRefractivityConstant = wetRefractivityConstant
self % num_fov = numFOV
self % fov_radius = fovRadius

! Reject unsupported values up front rather than silently sampling at undefined locations.
if (.NOT. ANY(numFOV == [1, 9])) then
  write(err_msg,"(A,I0)") "ufo_groundgnss_metoffice_setup: num_fov is set 1 or 9 currently, more options (e.g., num_fov=5) can be added "
  call abor1_ftn(err_msg)
end if

end subroutine ufo_groundgnss_metoffice_setup


! -----------------------------------------------------------------------------------------------
! Work out where the model state should be sampled for each observation.
! For num_fov = 9, 8 sample points are on the circle of radius (fov_radius) centered at obs locs.
! (that is, on the directions: 0=N, 45=NE, 90=E, 135=SE, 180=S, 225=SW, 270=W, 315=NW)
! The points of observation i occupy entries (i-1)*num_fov+1 : i*num_fov of lons/lats, which is
! the grouping that ObsGroundgnssMetOffice::locations() reports to OOPS.
!
! Another option of num_fov = 5 can be added in case of limited memory, using 4 sample points
! on the circle (that is, on the directions: 0=N, 90=E, 180=S, 270=W)
! -----------------------------------------------------------------------------------------------
subroutine ufo_groundgnss_metoffice_locs_init(self, obss, npaths, lons, lats)

implicit none

class(ufo_groundgnss_MetOffice), intent(in)    :: self
type(c_ptr), value,              intent(in)    :: obss     ! The observations and their meta-data
integer,                         intent(in)    :: npaths   ! nobs * num_fov
real(c_float),                   intent(inout) :: lons(npaths)
real(c_float),                   intent(inout) :: lats(npaths)

character(len=*), parameter  :: myname_ = "ufo_groundgnss_metoffice_locs_init"
integer                      :: iobs, nobs, num_fov
real(kind_real)              :: fov_radius               ! Sampling radius, km
real(kind_real), dimension(self%num_fov) :: lat_2d, lon_2d     ! Sample points
real(kind_real), allocatable :: lon(:), lat(:)           ! Observation positions
character(max_string)        :: err_msg

num_fov = self % num_fov
fov_radius = self % fov_radius
nobs = npaths / num_fov

if (num_fov /= 9) then
    write(err_msg,"(A,I0)") myname_ // ": unsupported num_fov ", num_fov
    call abor1_ftn(err_msg)
end if

allocate(lon(nobs), lat(nobs))
call obsspace_get_db(obss, "MetaData", "longitude", lon)   ! degrees
call obsspace_get_db(obss, "MetaData", "latitude", lat)    ! degrees

do iobs = 1, nobs
  call fov_sample_points(lat(iobs), lon(iobs), fov_radius, lat_2d, lon_2d, num_fov)
  lons( (iobs-1)*num_fov+1 : iobs*num_fov) = real(lon_2d(1:num_fov), c_float)
  lats( (iobs-1)*num_fov+1 : iobs*num_fov) = real(lat_2d(1:num_fov), c_float)
end do

deallocate(lon, lat)
end subroutine ufo_groundgnss_metoffice_locs_init

!-----------------------------------------------------------------------------
subroutine fov_sample_points(lat0, lon0, fov_radius, lat_2d, lon_2d, num_fov)
!-----------------------------------------------------------------------------
use ufo_constants_mod, only : mean_earth_rad_m, deg2rad, rad2deg
implicit none

integer, intent(in) :: num_fov
real(kind_real), intent(in) :: lat0, lon0     ! obs loc in degree
real(kind_real), intent(in) :: fov_radius     ! FOV average radius in km
real(kind_real), intent(out) :: lat_2d(num_fov), lon_2d(num_fov)    ! sample points

! For num_fov =9 (obs loc + 8 points along a circle of radius 25km centered at obs loc)
real(kind_real), parameter :: bearings(8) = [ &
                 0.0,   45.0,  90.0, 135.0, 180.0,  225.0, 270.0, 315.0 ]

real(kind_real) :: lat1, lon1, lat2, lon2
real(kind_real) :: theta, angular_dist, rearth
character(max_string)        :: err_msg
integer :: i

if( abs(lat0) > 89.0) then
  write(err_msg,*) " FOV sampling points not defined near the earth poles!"
  call abor1_ftn(err_msg)
end if

! Save the center point to the 1st of the sample points
   lat_2d(1) = lat0
   lon_2d(1) = lon0

   lat1 = lat0 * deg2rad
   lon1 = lon0 * deg2rad

! Distance on Earth's surface -> angular distance
rearth = mean_earth_rad_m *0.001       ! km
angular_dist = fov_radius / rearth

! Define the sample points along a circle centered at the obs location
   do i = 2, num_fov
      theta = bearings(i-1) * deg2rad
      lat2 = asin( sin(lat1)*cos(angular_dist) + &
                   cos(lat1)*sin(angular_dist)*cos(theta) )
      lon2 = lon1 + atan2(sin(theta)*sin(angular_dist)*cos(lat1), &
                   cos(angular_dist) - sin(lat1)*sin(lat2) )

      lat_2d(i) = lat2 * rad2deg
      lon_2d(i) = lon2 * rad2deg
   end do
end subroutine fov_sample_points

! ------------------------------------------------------------------------------
! Ground GNSS forward operator for the Met Office system
! ------------------------------------------------------------------------------
subroutine ufo_groundgnss_metoffice_simobs(self, geovals, hofx, obss)

  use fckit_exception_module, only: fckit_exception

  implicit none

  ! Arguments to this routine
  class(ufo_groundgnss_MetOffice), intent(in)    :: self     ! The object in which this operator is contained
  type(ufo_geovals),               intent(in)    :: geovals  ! The model values, interpolated to the obsevation locations
  real(kind_real),                 intent(inout) :: hofx(:)  ! The model forecast of the observations
  type(c_ptr), value,              intent(in)    :: obss     ! The observations, and meta-data for those observations

  character(len=*), parameter        :: myname_ = "ufo_groundgnss_metoffice_simobs"

  character(max_string)              :: err_msg         ! Error message for output
  character(max_string)              :: message         ! General message for output
  integer                            :: nobs            ! Number of observations
  integer                            :: iobs            ! Loop variable, observation number
  type(ufo_geoval), pointer          :: q               ! Model background values of specific humidity
  type(ufo_geoval), pointer          :: prs             ! Model background values of air pressure
  type(ufo_geoval), pointer          :: theta_heights   ! Model heights of levels containing specific humidity
  type(ufo_geoval), pointer          :: rho_heights     ! Model heights of levels containing air pressure

  real(kind_real), allocatable       :: zStation(:)

  ! Local variables
  INTEGER :: ilev, nlevp, nlevq, iflip
  REAL(kind_real), allocatable :: pressure(:)   ! Model background values of air pressure (monotonic order)
  REAL(kind_real), allocatable :: humidity(:)   ! Model background specific humidity  (in pressure monotonic order)
  REAL(kind_real), allocatable :: za(:)         ! Model heights of rho levs (in pressure monotonic order)
  REAL(kind_real), allocatable :: zb(:)         ! Model heights of theta levs (in pressure monotonic order)

  INTEGER :: num_fov, ifov, ipath             ! Field-of-view sampling counters

! get variables from geovals
  call ufo_geovals_get_var(geovals, var_q, q)               ! specific humidity
  call ufo_geovals_get_var(geovals, var_prsi, prs)          ! pressure
  call ufo_geovals_get_var(geovals, var_z, theta_heights)   ! Geopotential height of the normal model levels
  call ufo_geovals_get_var(geovals, var_zi, rho_heights)    ! Geopotential height of the pressure levels

  num_fov = self % num_fov

! check that geovals and hofx agree on how many obs there are.  With FOV averaging each
! observation is sampled along num_fov paths, so it is the profile count that must match.
  IF (num_fov > 1) THEN
    IF (q%nprofiles /= size(hofx)*num_fov .OR. prs%nprofiles /= size(hofx)*num_fov .OR. &
        theta_heights%nprofiles /= size(hofx)*num_fov .OR. &
        rho_heights%nprofiles /= size(hofx)*num_fov) THEN
      write(err_msg,*) myname_, " error: number of sampled profiles inconsistent!"
      call abor1_ftn(err_msg)
    END IF
  ELSE
    IF (geovals%nlocs /= size(hofx)) THEN
      write(err_msg,*) myname_, " error: nlocs inconsistent!"
      call abor1_ftn(err_msg)
    END IF
  END IF

  nlevp = prs % nval
  nlevq = q % nval

  nobs  = obsspace_get_nlocs(obss)
  allocate(zStation(nobs))

  call obsspace_get_db(obss, "MetaData", "stationElevation", zStation)

  hofx(:) = 0

  allocate(pressure(1:nlevp))
  allocate(humidity(1:nlevq))
  allocate(za(1:nlevp))
  allocate(zb(1:nlevq))

  if (nobs > 0) then
    obs_loop: do iobs = 1, nobs

      if (num_fov > 1) then
        ! Average the model state over the field-of-view sampling points before running the
        ! forward model, giving one representative profile per observation.  Equal weights.
        pressure(:) = 0.0
        humidity(:) = 0.0
        za(:) = 0.0
        zb(:) = 0.0
        do ifov = 1, num_fov
          ipath = (iobs-1)*num_fov + ifov
          pressure = pressure + prs % vals(:,ipath)
          humidity = humidity + q % vals(:,ipath)
          za = za + rho_heights % vals(:,ipath)
          zb = zb + theta_heights % vals(:,ipath)
        end do
        pressure = pressure / num_fov
        humidity = humidity / num_fov
        za = za / num_fov
        zb = zb / num_fov
      else
        pressure = prs % vals(:,iobs)
        humidity = q % vals(:,iobs)
        za = rho_heights % vals(:,iobs)
        zb = theta_heights % vals(:,iobs)
      end if

      call Ops_Groundgnss_ForwardModel(nlevp,                  &
                                       nlevq,                  &
                                       za(1:nlevp),            &
                                       zb(1:nlevq),            &
                                       pressure(1:nlevp),      &
                                       humidity(1:nlevq),      &
                                       self % vert_interp_ops, &
                                       self % pseudo_ops,      &
                                       self % min_temp_grad,   &
                                       self % dryRefractivityConstant, &
                                       self % wetRefractivityConstant, &
                                       1,                      &
                                       zStation(iobs),         &
                                       hofx(iobs))

      write(message,"(A,10I6)") "Size of hofx = ", shape(hofx)
      call fckit_log%debug(message)
      write(message,"(A,F12.4)") "hofx(iobs) = ", hofx(iobs)
      call fckit_log%debug(message)

    end do obs_loop
  end if

  deallocate(pressure)
  deallocate(humidity)
  deallocate(za)
  deallocate(zb)
  deallocate(zStation)

end subroutine ufo_groundgnss_metoffice_simobs
! ------------------------------------------------------------------------------

SUBROUTINE Ops_Groundgnss_ForwardModel(nlevp,                &
                                       nlevq,                &
                                       za,                   &
                                       zb,                   &
                                       pressure,             &
                                       humidity,             &
                                       vert_interp_ops,      &
                                       pseudo_ops,           &
                                       gbgnss_min_temp_grad, &
                                       dryRefractivityConstant, &
                                       wetRefractivityConstant, &
                                       nobs,                 &
                                       zStation,             &
                                       Model_ZTD)

INTEGER, INTENT(IN)            :: nlevp                    ! no. of p levels in state vec.
INTEGER, INTENT(IN)            :: nlevq                    ! no. of theta levels
REAL(kind_real), INTENT(IN)    :: za(1:nlevp)              ! heights of rho levs
REAL(kind_real), INTENT(IN)    :: zb(1:nlevq)              ! heights of theta levs
REAL(kind_real), INTENT(IN)    :: pressure(1:nlevp)        ! Model background pressure
REAL(kind_real), INTENT(IN)    :: humidity(1:nlevq)        ! Model background specific humidity
LOGICAL, INTENT(IN)            :: vert_interp_ops          ! Pressure varies exponentially with height?
LOGICAL, INTENT(IN)            :: pseudo_ops               ! Use pseudo-levels in calculation?
REAL(kind_real), INTENT(IN)    :: gbgnss_min_temp_grad     ! The minimum temperature gradient which is used
REAL(kind_real), INTENT(IN)    :: dryRefractivityConstant  ! Constant used to multiply the dry refractivity components
REAL(kind_real), INTENT(IN)    :: wetRefractivityConstant  ! Constant used to multiply the wet refractivity components
INTEGER, INTENT(IN)            :: nobs                     ! Number of observations

REAL(kind_real), INTENT(IN)    :: zStation               ! Station heights
REAL(kind_real), INTENT(INOUT) :: Model_ZTD              ! Model forecast of the observations

!
! Things that may need to be output, as they are used by the TL/AD calculation
!

REAL(kind_real)                  :: pN(nlevq)            ! Presure on theta levels
REAL(kind_real), ALLOCATABLE     :: refrac(:)            ! Model refractivity
LOGICAL                          :: refracerr            ! Refraction error
INTEGER                          :: nRefLevels           ! Number of levels in refractivity calculation
REAL(kind_real)                  :: TopCorrection        ! Zenith Total Delay Top of atmos correction
!
! Local parameters
!
integer, parameter           :: max_string = 800  ! Length of strings
character(len=*), parameter  :: myname_ = "Ops_Groundgnss_ForwardModel"
!
! Local variables
!
INTEGER                      :: nstate            ! no. of levels in state vec.
REAL(kind_real)              :: x(1:nlevp+nlevq)  ! state vector
character(max_string)        :: err_msg           ! Error message to be output
character(max_string)        :: message           ! General message for output

REAL(kind_real), ALLOCATABLE :: model_heights(:)  ! Geopotential heights of the refractivity levels (not needed for this oper)

! The model data must be on a staggered grid, with nlevp = nlevq+1
IF (nlevp /= nlevq + 1) THEN
    write(err_msg,*) myname_ // ":" // " Data must be on a staggered grid nlevp, nlevq = ", nlevp, nlevq
    call fckit_log % warning(err_msg)
    write(err_msg,*) myname_ // ":" // " error: number of levels inconsistent!"
    call abor1_ftn(err_msg)
END IF

! The ufo_utils_refractivity_calculator currently relies on the variables being provided
! bottom to top. The geovals are provided top to bottom and so the vertical variables need
! to be reveresed when they are passed into this routine and the outputs then need to be reversed
CALL ufo_calculate_refractivity(nlevp,                   &
                                nlevq,                   &
                                za(nlevp:1:-1),          &
                                zb(nlevq:1:-1),          &
                                pressure(nlevp:1:-1),    &
                                humidity(nlevq:1:-1),    &
                                vert_interp_ops,         &
                                pseudo_ops,              &
                                gbgnss_min_temp_grad,    &
                                dryRefractivityConstant, &
                                wetRefractivityConstant, &
                                refracerr,               &
                                nRefLevels,              &
                                refrac,                  &
                                model_heights)
! Flip the vertical direction of refrac and model_height back to top to bottom.
refrac = refrac(nRefLevels:1:-1)
model_heights = model_heights(nRefLevels:1:-1)

CALL Ops_groundgnss_TopCorrection(pressure,    &
                                  nlevq,       &
                                  za,          &
                                  zb,          &
                                  dryRefractivityConstant, &
                                  TopCorrection)

CALL Ops_Groundgnss_ZTD  (nlevq,     &
                          refrac,    &
                          zb,        &
                          zStation,  &
                          Model_ZTD)

Model_ZTD = Model_ZTD + TopCorrection

write(message,"(A,F16.14)") "Model_ZTD = ", Model_ZTD
call fckit_log%debug(message)

END SUBROUTINE ops_groundgnss_forwardmodel

END MODULE ufo_groundgnss_metoffice_mod
