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
use ufo_refractivity_utils_mod
use fckit_log_module,  only : fckit_log

implicit none
public             :: ufo_groundgnss_metoffice
private

  !> Fortran derived type for groundgnss trajectory
type, extends(ufo_basis) :: ufo_groundgnss_MetOffice
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
type(fckit_configuration), intent(in)   :: f_conf

end subroutine ufo_groundgnss_metoffice_setup


! ------------------------------------------------------------------------------
! Ground GNSS forward operator for the Met Office system
! ------------------------------------------------------------------------------
subroutine ufo_groundgnss_metoffice_simobs(self, geovals, hofx, obss)

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

  write(err_msg,*) "TRACE: ufo_groundgnss_metoffice_simobs: begin"
  call fckit_log%info(err_msg)

! check if nlocs is consistent in geovals & hofx
  if (geovals%nlocs /= size(hofx)) then
      write(err_msg,*) myname_, ' error: nlocs inconsistent!'
      call abor1_ftn(err_msg)
  endif


! get variables from geovals
  call ufo_geovals_get_var(geovals, var_q, q)               ! specific humidity
  call ufo_geovals_get_var(geovals, var_prsi, prs)          ! pressure
  call ufo_geovals_get_var(geovals, var_z, theta_heights)   ! Geopotential height of the normal model levels
  call ufo_geovals_get_var(geovals, var_zi, rho_heights)    ! Geopotential height of the pressure levels

  nobs  = obsspace_get_nlocs(obss)
  allocate(zStation(nobs))

  call obsspace_get_db(obss, "MetaData", "station_height", zStation)

  write(err_msg,*) "TRACE: ufo_groundgnss_metoffice_simobs: begin observation loop, nobs =  ", nobs
  call fckit_log%info(err_msg)

  obs_loop: do iobs = 1, nobs 

    call Ops_Groundgnss_ForwardModel(prs % nval, &
                                     q % nval, &
                                     rho_heights % vals(:,iobs), &
                                     theta_heights % vals(:,iobs), &
                                     prs % vals(:,iobs), &
                                     q % vals(:,iobs), &
                                     1, &
                                     zStation(iobs), &
                                     hofx(iobs))

    write(message,'(A,10I6)') "Size of hofx = ", shape(hofx)
    call fckit_log%info(message)
    write(message,'(A,1F2.4)') "hofx(iobs) = ", hofx(iobs)
    call fckit_log%info(message)
  end do obs_loop

  write(err_msg,*) "TRACE: ufo_groundgnss_metoffice_simobs: completed"
  call fckit_log%info(err_msg)

end subroutine ufo_groundgnss_metoffice_simobs
! ------------------------------------------------------------------------------


SUBROUTINE Ops_Groundgnss_ForwardModel(nlevp, &
                                       nlevq, &
                                       za, &
                                       zb, &
                                       pressure, &
                                       humidity, &
                                       nobs, &
                                       zStation, &
                                       Model_ZTD)

INTEGER, INTENT(IN)            :: nlevp                  ! no. of p levels in state vec.
INTEGER, INTENT(IN)            :: nlevq                  ! no. of theta levels
REAL(kind_real), INTENT(IN)    :: za(1:nlevp)            ! heights of rho levs
REAL(kind_real), INTENT(IN)    :: zb(1:nlevq)            ! heights of theta levs
REAL(kind_real), INTENT(IN)    :: pressure(1:nlevp)      ! Model background pressure
REAL(kind_real), INTENT(IN)    :: humidity(1:nlevq)      ! Model background specific humidity
INTEGER, INTENT(IN)            :: nobs                   ! Number of observations

REAL(kind_real), INTENT(IN)    :: zStation
REAL(kind_real), INTENT(INOUT) :: Model_ZTD          ! Model forecast of the observations

! 
! Things that may need to be output, as they are used by the TL/AD calculation
! 

REAL(kind_real)                  :: pN(nlevq)              ! Presure on theta levels
REAL(kind_real)                  :: refrac(nlevq)
LOGICAL                          :: refracerr

REAL(kind_real)                  :: TopCorrection
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
character(max_string)        :: message         ! General message for output

! The model data must be on a staggered grid, with nlevp = nlevq+1
IF (nlevp /= nlevq + 1) THEN
    write(err_msg,*) myname_ // ':' // ' Data must be on a staggered grid nlevp, nlevq = ', nlevp, nlevq
    call fckit_log % warning(err_msg)
    write(err_msg,*) myname_ // ':' // ' error: number of levels inconsistent!'
    call abor1_ftn(err_msg)
END IF

nstate = nlevp + nlevq
x(1:nlevp) = pressure
x(nlevp+1:nstate) = humidity

CALL ufo_refractivity(nlevq,     &
                      nlevp,     &
                      za,        &
                      zb,        &
                      x,         &
                      pN,        &
                      refracerr, &
                      refrac)

CALL Ops_groundgnss_TopCorrection(pN,    &
                                  nlevq, &
                                  TopCorrection)

CALL Ops_Groundgnss_ZTD  (nlevq,  &
                          refrac, &
                          zb,     &
                          zStation,  &
                          Model_ZTD)

Model_ZTD = Model_ZTD + TopCorrection

write(message,'(A,F10.4)') "Model_ZTD = ", Model_ZTD
call fckit_log%info(message)

END SUBROUTINE ops_groundgnss_forwardmodel

END MODULE ufo_groundgnss_metoffice_mod
