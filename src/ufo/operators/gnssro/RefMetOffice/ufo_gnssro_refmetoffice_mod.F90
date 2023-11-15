!-------------------------------------------------------------------------------
! (C) British Crown Copyright 2021 Met Office
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
!-------------------------------------------------------------------------------

!> Fortran module for gnssro refractivity Met Office forward operator

module ufo_gnssro_refmetoffice_mod

use iso_c_binding
use kinds
use ufo_vars_mod
use ufo_geovals_mod
use ufo_geovals_mod_c, only: ufo_geovals_registry
use vert_interp_mod
use lag_interp_mod,    only: lag_interp_const, lag_interp_smthWeights
use obsspace_mod  
use missing_values_mod
use ufo_utils_refractivity_calculator, only: ufo_calculate_refractivity
use fckit_log_module,  only : fckit_log
use fckit_exception_module, only: fckit_exception
use ufo_constants_mod, only: &
    rd,                      &    ! Gas constant for dry air
    grav,                    &    ! Gravitational field strength
    mw_ratio,                &    ! Ratio of molecular weights of water and dry air
    c_virtual,               &    ! Related to mw_ratio
    n_alpha,                 &    ! Refractivity constant a
    n_beta                        ! Refractivity constant b


implicit none
public             :: ufo_gnssro_refmetoffice
private

  !> Fortran derived type for gnssro trajectory
type :: ufo_gnssro_RefMetOffice
  logical :: vert_interp_ops
  logical :: pseudo_ops
  real(kind_real) :: min_temp_grad
  contains
    procedure :: setup     => ufo_gnssro_refmetoffice_setup
    procedure :: simobs    => ufo_gnssro_refmetoffice_simobs
end type ufo_gnssro_RefMetOffice

contains

!-------------------------------------------------------------------------------
!> \brief Set up the Met Office GNSS-RO refractivity operator
!!
!! \details **ufo_gnssro_refmetoffice_setup**
!! * Get the optional settings for the forward model, and save them in the
!!   object so that they can be used in the code.
!!
!! \author Neill Bowler (Met Office)
!!
!! \date 20 March 2021
!!
!-------------------------------------------------------------------------------
subroutine ufo_gnssro_refmetoffice_setup(self, vert_interp_ops, pseudo_ops, min_temp_grad)

implicit none

class(ufo_gnssro_refmetoffice), intent(inout) :: self
logical(c_bool), intent(in)  :: vert_interp_ops
logical(c_bool), intent(in)  :: pseudo_ops
real(c_float), intent(in)  :: min_temp_grad

self % vert_interp_ops = vert_interp_ops
self % pseudo_ops = pseudo_ops
self % min_temp_grad = min_temp_grad

end subroutine ufo_gnssro_refmetoffice_setup

!-------------------------------------------------------------------------------
!> \brief Calculate the model forecast of the observations
!!
!! \details **ufo_gnssro_refmetoffice_simobs**
!! * 1-dimensional GNSS-RO forward operator for the Met Office system
!!
!! \author Neill Bowler (Met Office)
!!
!! \date 20 March 2021
!!
!-------------------------------------------------------------------------------
subroutine ufo_gnssro_refmetoffice_simobs(self, geovals, obss, hofx, obs_diags)

  implicit none

  ! Arguments to this routine
  class(ufo_gnssro_RefMetOffice), intent(in)    :: self       ! The object in which this operator is contained
  type(ufo_geovals),              intent(in)    :: geovals    ! The model values, interpolated to the observation locations
  real(kind_real),                intent(inout) :: hofx(:)    ! The model forecast of the observations
  type(c_ptr), value,             intent(in)    :: obss       ! The observations, and meta-data for those observations
  type(ufo_geovals),              intent(inout) :: obs_diags  ! Observation diagnostics

  character(len=*), parameter        :: myname_ = "ufo_gnssro_refmetoffice_simobs"
  integer, parameter                 :: max_string = 800

  character(max_string)              :: err_msg               ! Error message for output
  character(max_string)              :: message               ! General message for output
  integer                            :: nobs                  ! Number of observations
  integer                            :: ilev                  ! Loop variable, level number
  integer                            :: iobs                  ! Loop variable, observation number
  type(ufo_geoval), pointer          :: q                     ! Model background values of specific humidity
  type(ufo_geoval), pointer          :: prs                   ! Model background values of air pressure
  type(ufo_geoval), pointer          :: theta_heights         ! Model heights of levels containing specific humidity
  type(ufo_geoval), pointer          :: rho_heights           ! Model heights of levels containing air pressure
  real(kind_real), allocatable       :: obsLat(:)             ! Latitude of the observation
  real(kind_real), allocatable       :: obsLon(:)             ! Longitude of the observation
  real(kind_real), allocatable       :: obs_height(:)         ! Geopotential height of the observation
  logical                            :: BAErr                 ! Was there an error in the calculation?
  integer                            :: iVar                  ! Loop variable, obs diagnostics variable number
  real(kind_real), allocatable       :: refractivity(:)       ! Refractivity on various model levels
  real(kind_real), allocatable       :: model_heights(:)      ! Geopotential heights that refractivity is calculated on

  write(err_msg,*) "TRACE: ufo_gnssro_refmetoffice_simobs: begin"
  call fckit_log%info(err_msg)

  ! If output to refractivity (and heights of the refractivity levels) is needed,
  ! then use nval as a way to check whether the array has been initialised (since
  ! it is called in a loop).
  DO iVar = 1, obs_diags % nvar
    IF (obs_diags % variables(ivar) == "atmosphericRefractivity" .OR. &
        obs_diags % variables(ivar) == "geopotentialHeight") THEN
      write(err_msg,*) "TRACE: ufo_gnssro_refmetoffice_simobs: initialising obs_diags for " // &
        obs_diags % variables(ivar)
      call fckit_log%info(err_msg)
      obs_diags % geovals(iVar) % nval = 0
    END IF
  END DO

! check if nlocs is consistent in geovals & hofx
  if (geovals%nlocs /= size(hofx)) then
      write(err_msg,*) myname_, ' error: nlocs inconsistent!'
      call abor1_ftn(err_msg)
  endif

  write(message, *) myname_, ' Running Met Office GNSS-RO forward operator with:'
  call fckit_log%info(message)
  write(message, *) 'vert_interp_ops =', self % vert_interp_ops, &
    'pseudo_ops =', self % pseudo_ops, 'min_temp_grad =', self % min_temp_grad
  call fckit_log%info(message)

! get variables from geovals
  call ufo_geovals_get_var(geovals, var_q, q)               ! specific humidity
  call ufo_geovals_get_var(geovals, var_prsi, prs)          ! pressure
  call ufo_geovals_get_var(geovals, var_z, theta_heights)   ! Geopotential height of the normal model levels
  call ufo_geovals_get_var(geovals, var_zi, rho_heights)    ! Geopotential height of the pressure levels

! make sure that the geovals are in the correct vertical order (top-to-bottom)
  if (prs%vals(1,1) > prs%vals(prs%nval,1) ) then
    write(err_msg,'(a)') 'Geovals should be ordered top to bottom'
    call fckit_exception%throw(err_msg)
  end if

  write(message, '(A,10I6)') 'Q: ', q%nval, q%nprofiles, shape(q%vals)
  call fckit_log%info(message)
  write(message, '(A,10I6)') 'Pressure: ', prs%nval, prs%nprofiles, shape(prs%vals)
  call fckit_log%info(message)

  nobs  = obsspace_get_nlocs(obss)

! set obs space struture
  allocate(obsLon(nobs))
  allocate(obsLat(nobs))
  allocate(obs_height(nobs))

  call obsspace_get_db(obss, "MetaData", "longitude", obsLon)
  call obsspace_get_db(obss, "MetaData", "latitude", obsLat)
  call obsspace_get_db(obss, "MetaData", "height", obs_height)

  obs_loop: do iobs = 1, nobs 

    call RefMetOffice_ForwardModel(prs % nval, &
                                   q % nval, &
                                   rho_heights % vals(prs%nval:1:-1,iobs), &
                                   theta_heights % vals(q%nval:1:-1,iobs), &
                                   prs % vals(prs%nval:1:-1,iobs), &
                                   q % vals(q%nval:1:-1,iobs), &
                                   self % pseudo_ops, &
                                   self % vert_interp_ops, &
                                   self % min_temp_grad, &
                                   1, &
                                   obs_height(iobs:iobs), &
                                   obsLat(iobs), &
                                   hofx(iobs:iobs), &
                                   BAErr, &
                                   refractivity, &
                                   model_heights)

    if (BAErr) then
      write(err_msg,*) "Error with observation processing ", iobs
      call fckit_log % info(err_msg)
    end if

    ! If required, then save the refractivity and model heights to the obs diagnostics.
    ! Allocate the output arrays on the first iteration.
    DO iVar = 1, obs_diags % nvar
        IF (obs_diags % variables(ivar) == "atmosphericRefractivity") THEN
            IF (iobs == 1) THEN
                obs_diags % geovals(iVar) % nval = SIZE(refractivity)
                ALLOCATE(obs_diags % geovals(iVar) % vals(SIZE(refractivity), obs_diags % nlocs))
            END IF

            IF (BAerr) THEN
                obs_diags % geovals(iVar) % vals(:,iobs) = missing_value(obs_diags % geovals(iVar) % vals(1,1))
            ELSE
                ! Flip the order of the calculated refractivity, since the geovals are oriented
                ! top-to-bottom, but the code works bottom-to-top
                obs_diags % geovals(iVar) % vals(:,iobs) = refractivity(SIZE(refractivity):1:-1)
            END IF
        END IF

        IF (obs_diags % variables(ivar) == "geopotentialHeight") THEN
            IF (iobs == 1) THEN
                obs_diags % geovals(iVar) % nval = SIZE(model_heights)
                ALLOCATE(obs_diags % geovals(iVar) % vals(SIZE(model_heights), obs_diags % nlocs))
            END IF

            IF (BAerr) THEN
                obs_diags % geovals(iVar) % vals(:,iobs) = missing_value(obs_diags % geovals(iVar) % vals(1,1))
            ELSE
                ! Flip the order of the calculated refractivity, since the geovals are oriented
                ! top-to-bottom, but the code works bottom-to-top
                obs_diags % geovals(iVar) % vals(:,iobs) = model_heights(SIZE(model_heights):1:-1)
            END IF
        END IF
    END DO
  end do obs_loop

  deallocate(obsLat)
  deallocate(obsLon)
  deallocate(obs_height)

  write(err_msg,*) "TRACE: ufo_gnssro_refmetoffice_simobs: completed"
  call fckit_log%info(err_msg)

end subroutine ufo_gnssro_refmetoffice_simobs

!-------------------------------------------------------------------------------
!> \brief Interface routine for the GNSS-RO refractivity forward operator
!!
!! \details **RefMetOffice_ForwardModel**
!! * Calculate the refractivity on model or pseudo levels
!! * Vertically interpolate the model refractivity to the observation locations
!!
!! \author Neill Bowler (Met Office)
!!
!! \date 20 March 2021
!!
!-------------------------------------------------------------------------------
SUBROUTINE RefMetOffice_ForwardModel(nlevp, &
                                     nlevq, &
                                     za, &
                                     zb, &
                                     pressure, &
                                     humidity, &
                                     GPSRO_pseudo_ops, &
                                     GPSRO_vert_interp_ops, &
                                     GPSRO_min_temp_grad, &
                                     nobs, &
                                     zobs, &
                                     Latitude, &
                                     ycalc, &
                                     BAErr, &
                                     refractivity, &
                                     model_heights)

INTEGER, INTENT(IN)            :: nlevp                  ! no. of p levels in state vec.
INTEGER, INTENT(IN)            :: nlevq                  ! no. of theta levels
REAL(kind_real), INTENT(IN)    :: za(1:nlevp)            ! heights of rho levs
REAL(kind_real), INTENT(IN)    :: zb(1:nlevq)            ! heights of theta levs
REAL(kind_real), INTENT(IN)    :: pressure(1:nlevp)      ! Model background pressure
REAL(kind_real), INTENT(IN)    :: humidity(1:nlevq)      ! Model background specific humidity
LOGICAL, INTENT(IN)            :: GPSRO_pseudo_ops       ! Option: Use pseudo-levels in vertical interpolation?
LOGICAL, INTENT(IN)            :: GPSRO_vert_interp_ops  ! Option: Use ln(p) for vertical interpolation? (rather than exner)
REAL(kind_real), INTENT(IN)    :: GPSRO_min_temp_grad    ! The minimum temperature gradient which is used
INTEGER, INTENT(IN)            :: nobs                   ! Number of observations in the profile
REAL(kind_real), INTENT(IN)    :: zobs(1:nobs)           ! Impact parameter for the obs
REAL(kind_real), INTENT(IN)    :: Latitude               ! Latitude of this profile
REAL(kind_real), INTENT(INOUT) :: ycalc(1:nobs)          ! Model forecast of the observations
LOGICAL, INTENT(OUT)           :: BAErr                  ! Was an error encountered during the calculation?
REAL(kind_real), INTENT(INOUT), ALLOCATABLE :: refractivity(:)     ! Model refractivity on model/pseudo levels
REAL(kind_real), INTENT(INOUT), ALLOCATABLE :: model_heights(:)    ! Geopotential heights of the refractivity levels
!
! Things that may need to be output, as they are used by the TL/AD calculation
! 
INTEGER                      :: nRefLevels          ! Number of levels in refractivity calculation
REAL(kind_real), ALLOCATABLE :: nr(:)               ! Model calculation of impact parameters
! 
! Local parameters
! 
integer, parameter           :: max_string = 800  ! Length of strings
character(len=*), parameter  :: myname_ = "Ops_GPSRO_ForwardModel"
!
! Local variables
! 
INTEGER                      :: num_pseudo        ! Number of levels, including pseudo levels
REAL(kind_real)              :: x(1:nlevp+nlevq)  ! state vector
character(max_string)        :: err_msg           ! Error message to be output
INTEGER                      :: iObs              ! Loop counter, observation number
INTEGER                      :: iLevel            ! Number of model level just above observation
REAL(kind_real)              :: temperature(nlevq)! Calculated profile temperature
REAL(kind_real)              :: Pb(nlevq)         ! Air pressure, interpolated to temperature-levels
REAL(kind_real)              :: humidity_ob       ! Specific humidity interpolated to observation height
REAL(kind_real)              :: T_ob              ! Temperature interpolated to observation height
REAL(kind_real)              :: P_ob              ! Pressure interpolated to observation height
REAL(kind_real)              :: log_humid         ! Interpolation coefficient for humidity
REAL(kind_real)              :: temperature_grad  ! Interpolation coefficient for temperature
REAL(kind_real)              :: pressure_coeff    ! Interpolation coefficient for pressure
REAL(kind_real)              :: alpha_interp      ! Interpolation coefficient for refractivity
REAL(kind_real)              :: beta_interp       ! Interpolation coefficient for refractivity
REAL(kind_real)              :: model_height_diff ! Difference in height between two model levels
REAL(kind_real)              :: obs_height_diff   ! Height difference between observation and closest model level

! The model data must be on a staggered grid, with nlevp = nlevq+1
IF (nlevp /= nlevq + 1) THEN
    write(err_msg,*) myname_ // ':' // ' Data must be on a staggered grid nlevp, nlevq = ', nlevp, nlevq
    call fckit_log % warning(err_msg)
    write(err_msg,*) myname_ // ':' // ' error: number of levels inconsistent!'
    call abor1_ftn(err_msg)
END IF

BAErr = .FALSE.
ycalc(:) = missing_value(ycalc(1))

! Calculate the refractivity on model or pseudo levels
CALL ufo_calculate_refractivity (nlevp,                 &
                                 nlevq,                 &
                                 za,                    &
                                 zb,                    &
                                 pressure,              &
                                 humidity,              &
                                 GPSRO_pseudo_ops,      &
                                 GPSRO_vert_interp_ops, &
                                 GPSRO_min_temp_grad,   &
                                 BAerr,                 &
                                 nRefLevels,            &
                                 refractivity,          &
                                 model_heights,         &
                                 temperature=temperature, &
                                 interp_pressure=Pb)

! Vertically interpolate the model refractivity to the observation locations
DO iObs = 1, nobs

  ! Use separate interpolation of input quantities in refractivity calculation.
  ! This does not use the refractivity values calculated above, but uses the
  ! calculated temperature, and pressure on temperature-levels.
  IF (GPSRO_pseudo_ops) THEN
    IF (zb(nlevq) > zobs(iObs) .AND. zb(1) < zobs(iObs)) THEN
      ! Find the model level immediately above the observation
      iLevel = 1
      DO
        IF (zobs(iObs) < zb(iLevel + 1)) EXIT
        iLevel = iLevel + 1
      END DO
      
      model_height_diff = zb(iLevel + 1) - zb(iLevel)
      obs_height_diff = zobs(iObs) - zb(iLevel)

      ! Interpolate P,T,q separately
      IF (MIN (humidity(iLevel - 1), humidity(iLevel)) > 0.0) THEN
        ! humidity varies exponentially with height
        log_humid = LOG(humidity(iLevel) / humidity(iLevel + 1)) / model_height_diff
        humidity_ob = humidity(iLevel) * EXP(-log_humid * obs_height_diff)
        ! Assume linear variation if humidities are -ve
      ELSE
        humidity_ob = humidity(iLevel) + (humidity(iLevel + 1) - humidity(iLevel)) / &
            model_height_diff * obs_height_diff
      END IF

      ! T varies linearly with height
      temperature_grad = (temperature(iLevel + 1) - temperature(iLevel)) / &
                         model_height_diff
      T_ob = temperature(iLevel) + temperature_grad * obs_height_diff

      ! P varies to maintain hydrostatic balance
      IF (ABS (temperature(iLevel + 1) - temperature(iLevel)) > GPSRO_min_temp_grad) THEN
        pressure_coeff = ((Pb(iLevel + 1) / Pb(iLevel)) * (temperature(iLevel + 1) / temperature(iLevel)) ** &
                          (grav / (rd * temperature_grad)) - 1.0) / model_height_diff
        P_ob = (Pb(iLevel) * (T_ob / temperature(iLevel)) ** (-grav / (rd * temperature_grad))) * &
               (1.0 + pressure_coeff * obs_height_diff)
      ELSE
        ! If layer is isothermal, assume exponential variation to
        ! avoid singularity
        P_ob = Pb(iLevel) * EXP (LOG (Pb(iLevel + 1) / Pb(iLevel)) * &
               obs_height_diff / model_height_diff)
      END IF


      ! Calculate refractivity

      ycalc(iObs) = n_alpha * P_ob / T_ob + n_beta * P_ob * humidity_ob / (T_ob ** 2 * &
                (mw_ratio + (1.0 - mw_ratio) * humidity_ob))
    END IF

  ELSE

    IF (model_heights(nRefLevels) > zobs(iObs) .AND. model_heights(1) < zobs(iObs)) THEN
      ! Find the model level immediately above the observation
      iLevel = 1
      DO
        IF (zobs(iObs) < model_heights(iLevel + 1)) EXIT
        iLevel = iLevel + 1
      END DO

      ! Use simple assumption of exponentially varying refractivity
      alpha_interp = (model_heights(iLevel + 1) - zobs(iObs)) / &
                     (model_heights(iLevel + 1) - model_heights(iLevel))
      beta_interp = 1.0 - alpha_interp
      ycalc(iObs) = EXP(alpha_interp * LOG (refractivity(iLevel)) + &
                        beta_interp * LOG (refractivity(iLevel + 1)))
    ELSE IF (zobs(iObs) /= missing_value(zobs(iObs))) THEN
      PRINT*, 'Height out of range', iObs, nRefLevels, model_heights(nRefLevels), zobs(iObs)
    END IF
  END IF

END DO

END SUBROUTINE RefMetOffice_ForwardModel

end module ufo_gnssro_refmetoffice_mod
