!-------------------------------------------------------------------------------
! (C) British Crown Copyright 2025 Met Office
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
!-------------------------------------------------------------------------------

!> Fortran module for gnssro refractivity 1DVar check

module ufo_refractivityonedvarcheck_do1dvar_mod

use kinds, only: kind_real
use missing_values_mod, only: missing_value
use logger_mod, only: oops_log

private
public :: Refractivity_Do1DVar

contains

!-------------------------------------------------------------------------------
! Find a solution to the satellite sounding inverse problem
!-------------------------------------------------------------------------------
SUBROUTINE Refractivity_Do1DVar (nlevp,                  &
                              nlevq,                  &
                              BM1,                    &
                              Bsig,                   &
                              Back,                   &
                              Ob,                     &
                              RMatrix,                &
                              GPSRO_pseudo_ops,       &
                              GPSRO_vert_interp_ops,  &
                              GPSRO_min_temp_grad,    &
                              GPSRO_cost_funct_test,  &   ! Threshold value for the cost function convergence test
                              GPSRO_y_test,           &   ! Threshold value for the yobs-ysol tes
                              GPSRO_n_iteration_test, &   ! Maximum number of iterations
                              GPSRO_Delta_factor,     &   ! Delta
                              GPSRO_Delta_ct2,        &   ! Delta observations
                              GPSRO_OB_test,          &   ! Threshold value for the O-B test
                              minval_ytest,           &   ! Minimum value for RHS of y-test
                              maxval_ytest,           &   ! Maximum value for RHS of y-test
                              capsupersat,            &
                              RefracErr,              &
                              Tb,                     &
                              Ts,                     &
                              O_Bdiff,                &
                              DFS)

use ufo_gnssroonedvarcheck_setom1_mod, only: &
    Ops_GPSRO_setOM1

use ufo_refractivityonedvarcheck_utils_mod, only: &
    singlerefbg_type,             &
    singlerefob_type

use ufo_refractivityonedvarcheck_rootsolv_mod, only: &
    Ops_GPSRO_rootsolv

use ufo_roobserror_utils_mod, only: &
    Rmatrix_type

use ufo_utils_refractivity_calculator, only: &
    ufo_calculate_refractivity

IMPLICIT NONE

! Subroutine arguments:
INTEGER, INTENT(IN)                    :: nlevp                   !< Number of pressure levels
INTEGER, INTENT(IN)                    :: nlevq                   !< Number of specific humidity (or temperature) levels
REAL(kind_real), INTENT(IN)            :: BM1(:,:)                !< Inverse of the background-error covariance matrix
REAL(kind_real), INTENT(IN)            :: Bsig(:)                 !< Diagonal of the background-error covariance matrix
TYPE (SingleRefBg_type), INTENT(INOUT) :: Back                    !< Background state (i.e. geovals)
TYPE (SingleRefOb_type), INTENT(INOUT) :: Ob                      !< Observations
TYPE(RMatrix_type), INTENT(IN)         :: RMatrix                 !< Observation error covariance matrix information
LOGICAL, INTENT(IN)                    :: GPSRO_pseudo_ops        !< Whether to use pseudo-levels in the vertical interpolation
LOGICAL, INTENT(IN)                    :: GPSRO_vert_interp_ops   !< Use log(p) for vertical interpolation?
REAL(kind_real), INTENT(IN)            :: GPSRO_min_temp_grad     !< Threshold for vertical temperature gradient
REAL(kind_real), INTENT(IN)            :: GPSRO_cost_funct_test   !< Threshold value for the cost function convergence test
REAL(kind_real), INTENT(IN)            :: GPSRO_y_test            !< Threshold value for the yobs-ysol tes
INTEGER, INTENT(IN)                    :: GPSRO_n_iteration_test  !< Maximum number of iterations
REAL(kind_real), INTENT(IN)            :: GPSRO_Delta_ct2         !< Threshold for second convergence test (change in state divided by gradient)
REAL(kind_real), INTENT(IN)            :: GPSRO_Delta_factor      !< Threshold for change in the cost function for convergence
REAL(kind_real), INTENT(IN)            :: GPSRO_OB_test           !< Threshold value for the O-B test
REAL(kind_real), INTENT(IN)            :: minval_ytest            !< Minimum value for RHS of y-test
REAL(kind_real), INTENT(IN)            :: maxval_ytest            !< Maximum value for RHS of y-test
LOGICAL, INTENT(IN)                    :: capsupersat             !< Remove super-saturation?
LOGICAL, INTENT(OUT)                   :: RefracErr               !< Were there errors in the refractivity calculation?
REAL(kind_real), INTENT(INOUT)         :: Tb(nlevq)               !< Calculated background temperature
REAL(kind_real), INTENT(INOUT)         :: Ts(nlevq)               !< 1DVar solution temperature
REAL(kind_real), INTENT(INOUT)         :: O_Bdiff                 !< Measure of O-B for whole profile
REAL(kind_real), INTENT(INOUT)         :: DFS                     !< Measure of degrees of freedom of signal for whole profile

! Local parameters
INTEGER, PARAMETER                  :: max_string = 800

! Local variables
INTEGER                             :: nobs        ! size of the 1DVar observation vector
INTEGER                             :: i
INTEGER                             :: j
INTEGER                             :: it
LOGICAL                             :: OM1_error
LOGICAL                             :: converged
LOGICAL                             :: Do1DVar_error
LOGICAL                             :: ran_iteration
REAL(kind_real)                     :: J_pen
REAL(kind_real)                     :: xb(nlevp+nlevq)                 ! background profile used in the 1D-Var
REAL(kind_real)                     :: x(nlevp+nlevq)                  ! 1Dvar solution profile
REAL(kind_real)                     :: Amat(nlevp+nlevq,nlevp+nlevq)   ! solultion error cov matrix
REAL(kind_real), ALLOCATABLE        :: zobs(:)
REAL(kind_real), ALLOCATABLE        :: yobs(:)
REAL(kind_real), ALLOCATABLE        :: yb(:)
REAL(kind_real), ALLOCATABLE        :: ycalc(:)
REAL(kind_real), ALLOCATABLE        :: Om1(:,:)
REAL(kind_real), ALLOCATABLE        :: OSigma(:)
INTEGER, ALLOCATABLE                :: index_packed(:)
REAL(kind_real), ALLOCATABLE        :: model_heights(:)                ! Heights of model and pseudo-levels
REAL(kind_real), ALLOCATABLE        :: refractivity(:)                 ! Refractivity on model and pseudo_levels
INTEGER                             :: nRefLevels                      ! Number of levs to calculate ref on
CHARACTER(LEN=max_string)           :: message

!--------------
! 1. Initialise
!--------------

Do1DVar_error = .FALSE.
RefracErr = .FALSE.
ran_iteration = .FALSE.
OM1_error = .FALSE.

! Set all the PGE values to gross error

Ob % Refractivity(:) % PGEFinal = 1.0

! Set the background vector
xb(1:nlevp) = 1.0E-2 * Back % p(:)                 ! in hPa
xb(nlevp + 1:nlevp+nlevq) = 1.0E3 * Back % q(:)    ! in g/kg

nobs = COUNT (Ob % refractivity(:) % value /= missing_value(Ob % refractivity(1) % value) .AND. & ! not missing refractivity
              Ob % z(:) % value /= missing_value(Ob % z(1) % value) .AND. &                       ! not missing height
              Ob % qc_flags(:) == 0)                                                              ! not flagged as bad

WRITE (message, '(A,I0)') 'size of input obs vector ', SIZE (Ob % Refractivity(:) % value)
call oops_log % info(message)
WRITE (message, '(A,I0)') 'size of packed obs vector ', nobs
call oops_log % info(message)

! Only continue if we have some observations to process
IF (nobs > 0) THEN

  ! calculate an array of indices of the packed elements

  ALLOCATE (index_packed(nobs))                           ! allocate the packed index vector
  index_packed = missing_value(index_packed(1))           ! initialise

  ! Allocate arrays used in 1D-Var after test to stop allocating size nobs=0
  ALLOCATE (om1(nobs,nobs))
  ALLOCATE (yobs(nobs))
  ALLOCATE (zobs(nobs))
  ALLOCATE (yb(nobs))
  ALLOCATE (ycalc(nobs))
  ALLOCATE (OSigma(nobs))

  ! Pack observation arrays for valid values
  ! Note: This hard-codes the R-matrix to be diagonal, since that is all that
  ! is currently available in JEDI.  This will need to be revisted once the
  ! full capability is available.

  om1 = 0
  j = 1
  DO i = 1, SIZE (Ob % refractivity(:) % Value)
    IF (ABS(Ob % refractivity(i) % Value - missing_value(Ob % refractivity(i) % Value)) > 1 .AND. &
        Ob % z(i) % value /= missing_value(Ob % z(i) % value)           .AND. &
        Ob % qc_flags(i) == 0) THEN
      index_packed(j) = i
      zobs(j) = Ob % z(i) % value
      yobs(j) = Ob % refractivity(i) % Value
      om1(j,j) = (Ob % refractivity(i) % oberr)**(-2)
      j = j + 1
    END IF
  END DO
  
  CALL Ops_GPSRO_setOM1( &
    nobs,     &  ! Number of observations
    zobs,     &  ! Height of observations
    yobs,     &  ! Observed refractivity
    Rmatrix,  &  ! Observation error matrix information
    OSigma,   &  ! Standard deviation of observation errors
    OM1,      &  ! Inverse observation error matrix
    OM1_error &  ! Error flag
  )

  !-----------------------------------------------
  ! 2. If no errors so far, call the 1DVar routine
  !-----------------------------------------------

  IF (ALL(zobs(:) /= missing_value(zobs(1))) .AND.  &
      ALL(yobs(:) /= missing_value(yobs(1))) .AND.  &
      .NOT. OM1_error) THEN

    CALL Ops_GPSRO_rootsolv ( &
      nlevp+nlevq,               &    ! size of state vector
      nlevp,                     &    ! no. of press. levels
      nlevq,                     &    ! no. of theta levels
      nlevq,                     &    ! no. of q levels
      Nobs,                      &    ! no of obs
      Back % za,                 &    ! height of rho levels
      Back % zb,                 &    ! height of theta levels
      xb,                        &    ! background vector
      yobs,                      &    ! ob. vector
      zobs,                      &    ! ob. impact parameters
      Bsig,                      &    ! standard dev. of B errors
      Bm1,                       &    ! Inverse Back. cov matrix
      Om1,                       &    ! Inverse Ob cov matrix
      it,                        &    ! no of iterations
      x,                         &    ! solution vector
      yb,                        &    ! obs at first guess
      ycalc,                     &    ! obs at solution
      J_pen,                     &    ! penalty value
      Amat,                      &    ! solution cov. matrix
      converged,                 &    ! convergence flag
      RefracErr,                 &    ! error flag
      Do1DVar_error,             &    ! error flag
      GPSRO_n_iteration_test,    &    !
      GPSRO_Delta_factor,        &    !
      GPSRO_Delta_ct2,           &    !
      GPSRO_pseudo_ops,          &
      GPSRO_vert_interp_ops,     &
      GPSRO_min_temp_grad,       &
      capsupersat,               &
      O_Bdiff,                   &    ! observed -background BA value
      Tb,                        &
      Ts,                        &
      DFS                        &
    )
    ran_iteration = .TRUE.
  ELSE

    Do1DVar_error = .TRUE.

    Ob % refractivity(:) % PGEFinal = 0.99

  END IF

  IF (.NOT. Do1DVar_error) THEN

    ! store iteration and cost

    Ob % Niter = it

    Ob % Jcost = 2.0_kind_real * J_pen / REAL (nobs)

    ! map the retrieval information back into the ob structures

    Ob % p(:) % Value = 1.0E2 * x(1:nlevp)

    DO i = 1, nlevp

      Ob % p(i) % ObErr = 1.0E2 * SQRT(Amat(i,i))   ! error. est. from cov

    END DO

    Ob % q(:) % Value = 1.0E-3 * x(nlevp + 1:nlevp+nlevq)

    DO i = 1, nlevq

      j = i + nlevp

      Ob % q(i) % ObErr = 1.0E-3 *  SQRT(Amat(j,j))   ! error est. from cov

    END DO

    ! Refractivity calculated with solution

    Ob % SolutRefractivity(index_packed) = ycalc(1:nobs)

    ! PROBABILTY OF GROSS ERROR. Use the cost function value.

    IF (J_pen > GPSRO_cost_funct_test * REAL(NOBS, kind=kind_real)) THEN

      ! the cost function for profile is too high - GROSS ERROR -
      ! set all pge's in profile to 0.8

      Ob % refractivity(:) % PGEFinal = 0.8

    ELSE

      ! For each value in profile, estimate probability of gross error
      ! from difference between solution and observed value

      DO i = 1,nobs

        IF (ABS(Ob % refractivity(index_packed(i)) % Value - &
                Ob % SolutRefractivity(index_packed(i))) < &
                MIN(maxval_ytest, &
                    MAX(minval_ytest, &
                        GPSRO_y_test * Ob % refractivity(index_packed(i)) % ObErr))) THEN
          Ob % refractivity(index_packed(i)) % PGEFinal = 0.1
        ELSE
          Ob % refractivity(index_packed(i)) % PGEFinal = 0.7
        END IF

      END DO

    END IF

    ! if the initial 2J/m value exceeds read-in value then flag
    IF  (O_Bdiff > GPSRO_OB_test) THEN
      Ob % refractivity(:) % PGEFinal = 0.6
    END IF

    ! check for the RefracErr being set
    IF (ran_iteration .AND. RefracErr) THEN
      Ob % refractivity(:) % PGEFinal = 0.58     ! flag RefracErr
    END IF

  ELSE  ! Do 1D-Var_error

    IF (ran_iteration) THEN
      Ob % refractivity(:) % PGEFinal = 0.9     ! flag lack of convergence

      ! Refractivity calculated with solution
      Ob % SolutRefractivity(index_packed) = ycalc(1:nobs)

      ! map the retrieval information back into the ob structures
      Ob % p(:) % Value = 1.0E2 * x(1:nlevp)
      Ob % q(:) % Value = 1.0E-3 * x(nlevp + 1:nlevp+nlevq)

      ! store iterations and cost
      Ob % Niter = it
      Ob % Jcost = 2.0_kind_real * J_pen / nobs

    ELSE IF (OM1_error) THEN

      Ob % refractivity(:) % PGEFinal = 0.85    ! Flag error in getting
                                                ! observation error inverse
    END IF

  END IF

ELSE
  IF (nobs <= 10) THEN
    WRITE (message, '(A)') 'nobs is less than 10: exit Refractivity_Do1DVar'
    call oops_log % info(message)

    Ob % refractivity(:) % PGEFinal = 0.55     ! flag lack of observation data
  END IF

  IF (RefracErr) THEN
    WRITE (message, '(A)') 'Error in Ops_Refractivity: exit Refractivity_Do1DVar'
    call oops_log % info(message)
    Ob % refractivity(:) % PGEFinal = 0.58     ! flag RefracErr
  END IF
END IF

END SUBROUTINE Refractivity_Do1DVar

end module ufo_refractivityonedvarcheck_do1dvar_mod
