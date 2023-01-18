!-------------------------------------------------------------------------------
! (C) British Crown Copyright 2020 Met Office
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
!-------------------------------------------------------------------------------

!> Fortran module for ground based gnss Met Office forward operator

module ufo_groundgnss_metoffice_mod

use iso_c_binding
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

  !> Fortran derived type for groundgnss trajectory
type, extends(ufo_basis) :: ufo_groundgnss_MetOffice
  logical :: vert_interp_ops
  logical :: pseudo_ops
  real(kind_real) :: min_temp_grad
  contains
    procedure :: setup     => ufo_groundgnss_metoffice_setup
    procedure :: simobs    => ufo_groundgnss_metoffice_simobs
end type ufo_groundgnss_MetOffice

contains

! ------------------------------------------------------------------------------
! Get the optional settings for the forward model, and save them in the object
! so that they can be used in the code.
! ------------------------------------------------------------------------------
subroutine ufo_groundgnss_metoffice_setup(self, f_conf)

use fckit_configuration_module, only: fckit_configuration
implicit none
class(ufo_groundgnss_MetOffice), intent(inout) :: self
type(fckit_configuration), intent(in)          :: f_conf

call f_conf%get_or_die("vert_interp_ops", self % vert_interp_ops)
call f_conf%get_or_die("pseudo_ops", self % pseudo_ops)
call f_conf%get_or_die("min_temp_grad", self % min_temp_grad)

end subroutine ufo_groundgnss_metoffice_setup


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
  integer, parameter                 :: max_string = 800

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


  write(err_msg,*) "TRACE: ufo_groundgnss_metoffice_simobs: begin"
  call fckit_log%info(err_msg)

! check if nlocs is consistent in geovals & hofx
  IF (geovals%nlocs /= size(hofx)) THEN
    write(err_msg,*) myname_, ' error: nlocs inconsistent!'
    call abor1_ftn(err_msg)
  END IF

! get variables from geovals
  call ufo_geovals_get_var(geovals, var_q, q)               ! specific humidity
  call ufo_geovals_get_var(geovals, var_prsi, prs)          ! pressure
  call ufo_geovals_get_var(geovals, var_z, theta_heights)   ! Geopotential height of the normal model levels
  call ufo_geovals_get_var(geovals, var_zi, rho_heights)    ! Geopotential height of the pressure levels

! Check that the geovals are ordered top to bottom
  if( prs%vals(1,1) > prs%vals(prs%nval,1) ) then
    write(err_msg,'(a)') 'Geovals should be ordered top to bottom'
    call fckit_exception%throw(err_msg)
  endif

  nlevp = prs % nval
  nlevq = q % nval

  nobs  = obsspace_get_nlocs(obss)
  allocate(zStation(nobs))

  call obsspace_get_db(obss, "MetaData", "stationElevation", zStation)

  write(err_msg,*) "TRACE: ufo_groundgnss_metoffice_simobs: begin observation loop, nobs =  ", nobs
  call fckit_log%info(err_msg)

  hofx(:) = 0

  allocate(pressure(1:nlevp)) 
  allocate(humidity(1:nlevq))
  allocate(za(1:nlevp))
  allocate(zb(1:nlevq))

  obs_loop: do iobs = 1, nobs

    pressure = prs % vals(:,iobs)
    humidity = q % vals(:,iobs)
    za = rho_heights % vals(:,iobs)
    zb = theta_heights % vals(:,iobs)

    call Ops_Groundgnss_ForwardModel(nlevp,                  &
                                     nlevq,                  &
                                     za(1:nlevp),            &
                                     zb(1:nlevq),            &
                                     pressure(1:nlevp),      &
                                     humidity(1:nlevq),      &
                                     self % vert_interp_ops, &
                                     self % pseudo_ops,      &
                                     self % min_temp_grad,   &
                                     1,                      &
                                     zStation(iobs),         &
                                     hofx(iobs))

    write(message,'(A,10I6)') "Size of hofx = ", shape(hofx)
    call fckit_log%debug(message)
    write(message,'(A,F12.4)') "hofx(iobs) = ", hofx(iobs)
    call fckit_log%debug(message)

  end do obs_loop

  write(err_msg,*) "TRACE: ufo_groundgnss_metoffice_simobs: completed"
  call fckit_log%info(err_msg)
  
  deallocate(pressure)
  deallocate(humidity)
  deallocate(za)
  deallocate(zb)

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
                                       nobs,                 &
                                       zStation,             &
                                       Model_ZTD)

INTEGER, INTENT(IN)            :: nlevp                  ! no. of p levels in state vec.
INTEGER, INTENT(IN)            :: nlevq                  ! no. of theta levels
REAL(kind_real), INTENT(IN)    :: za(1:nlevp)            ! heights of rho levs
REAL(kind_real), INTENT(IN)    :: zb(1:nlevq)            ! heights of theta levs
REAL(kind_real), INTENT(IN)    :: pressure(1:nlevp)      ! Model background pressure
REAL(kind_real), INTENT(IN)    :: humidity(1:nlevq)      ! Model background specific humidity
LOGICAL, INTENT(IN)            :: vert_interp_ops        ! Pressure varies exponentially with height?
LOGICAL, INTENT(IN)            :: pseudo_ops             ! Use pseudo-levels in calculation?
REAL(kind_real), INTENT(IN)    :: gbgnss_min_temp_grad   ! The minimum temperature gradient which is used
INTEGER, INTENT(IN)            :: nobs                   ! Number of observations

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
    write(err_msg,*) myname_ // ':' // ' Data must be on a staggered grid nlevp, nlevq = ', nlevp, nlevq
    call fckit_log % warning(err_msg)
    write(err_msg,*) myname_ // ':' // ' error: number of levels inconsistent!'
    call abor1_ftn(err_msg)
END IF

! The ufo_utils_refractivity_calculator currently relies on the variables being provided
! bottom to top. The geovals are provided top to bottom and so the vertical variables need
! to be reveresed when they are passed into this routine and the outputs then need to be reversed
CALL ufo_calculate_refractivity(nlevp,                  &
                                nlevq,                  &
                                za(nlevp:1:-1),         &
                                zb(nlevq:1:-1),         &
                                pressure(nlevp:1:-1),   &
                                humidity(nlevq:1:-1),   &
                                vert_interp_ops,        &
                                pseudo_ops,             &
                                gbgnss_min_temp_grad,   &
                                refracerr,              &
                                nRefLevels,             &
                                refrac,                 &
                                model_heights)
! Flip the vertical direction of refrac and model_height back to top to bottom.
refrac = refrac(nRefLevels:1:-1)
model_heights = model_heights(nRefLevels:1:-1)

CALL Ops_groundgnss_TopCorrection(pressure,    &
                                  nlevq,       &
                                  za,          &
                                  zb,          &
                                  TopCorrection)


CALL Ops_Groundgnss_ZTD  (nlevq,     &
                          refrac,    &
                          zb,        &
                          zStation,  &
                          Model_ZTD)

Model_ZTD = Model_ZTD + TopCorrection

write(message,'(A,F16.14)') "Model_ZTD = ", Model_ZTD
call fckit_log%debug(message)


END SUBROUTINE ops_groundgnss_forwardmodel

END MODULE ufo_groundgnss_metoffice_mod
