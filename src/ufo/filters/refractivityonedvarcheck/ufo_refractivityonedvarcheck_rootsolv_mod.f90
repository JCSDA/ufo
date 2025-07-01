!-------------------------------------------------------------------------------
! (C) British Crown Copyright 2025 Met Office
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
!-------------------------------------------------------------------------------

!> Fortran module for solver of the refractivity 1DVar with gnssro

module ufo_refractivityonedvarcheck_rootsolv_mod

use kinds, only: kind_real
use missing_values_mod, only: missing_value
use logger_mod, only: oops_log
use ufo_gnssro_refmetoffice_mod, only: RefMetOffice_ForwardModel
use ufo_gnssro_refmetoffice_tlad_mod, only: jacobian_interface

private
public :: Ops_GPSRO_rootsolv

contains

!-------------------------------------------------------------------------------
! Solve the 1dvar problem
!-------------------------------------------------------------------------------

SUBROUTINE Ops_GPSRO_rootsolv( &
  nstate,        &   ! size of state vector
  nlevp,         &   ! no. of press. levels
  nb,            &   ! no. of theta levels
  nlevq,         &   ! no. of q levels
  Nobs,          &   ! no of obs
  za,            &   ! height of rho levels
  zb,            &   ! height of theta levels
  xb,            &   ! background vector
  yobs,          &   ! ob. vector
  zobs,          &   ! ob. heights
  Bsig,          &   ! standard dev. of b errors
  Bm1,           &   ! Inverse Back. cov matrix
  Om1,           &   ! Inverse Ob cov matrix
  it,            &   ! no of iterations
  x,             &   ! solution vector
  yb,            &   ! obs first guess
  ycalc,         &   ! obs at solution
  J_pen,         &   ! penalty value
  Amat,          &   ! solution cov. matrix
  converged,     &   ! convergence flag
  RefracErr,     &   ! error flag
  Do1DVar_error, &   ! error flag
  Iter_max,      &   ! Max number of iterations
  Delta,         &   !
  Delta_ct2,     &   !
  GPSRO_pseudo_ops,      & ! Whether to use pseudo levels
  GPSRO_vert_interp_ops, & ! Whether to vertically interpolate using exner or ln(p)
  GPSRO_min_temp_grad, &   ! Minimum vertical temperature gradient allowed
  capsupersat,   &
  O_Bdiff,       &   ! observed-background refractivity value
  Tb,            &
  Ts,            &
  DFS            &
)


USE ufo_gnssroonedvarcheck_humidcheck_mod, only: &
    Ops_GPSRO_humidcheck

USE ufo_gnssroonedvarcheck_pen_mod, only: &
    Ops_GPSRO_pen

USE ufo_gnssroonedvarcheck_eval_derivs_mod, only: &
    Ops_GPSRO_eval_derivs

USE ufo_utils_mod, only: &
    InvertMatrix, &
    Ops_Cholesky

use ufo_utils_refractivity_calculator, only: &
    ufo_calculate_refractivity, &
    ufo_refractivity_kmat

IMPLICIT NONE

! Subroutine arguments:
INTEGER, INTENT(IN)            :: nstate                 !< Number of state variables (=nlevp+nlevq)
INTEGER, INTENT(IN)            :: nlevp                  !< Number of pressure levels
INTEGER, INTENT(IN)            :: nb                     !< Number of theta levels
INTEGER, INTENT(IN)            :: nlevq                  !< Number of specific humidity levels
INTEGER, INTENT(IN)            :: Nobs                   !< Number of observations
REAL(kind_real), INTENT(IN)    :: za(:)                  !< Geopotential height of pressure levels
REAL(kind_real), INTENT(IN)    :: zb(:)                  !< Geopotential height of theta levels
REAL(kind_real), INTENT(IN)    :: xb(:)                  !< Background state (combining pressure and specific humidity)
REAL(kind_real), INTENT(IN)    :: yobs(:)                !< Observed values
REAL(kind_real), INTENT(IN)    :: zobs(:)                !< Geopotential height of the observations
REAL(kind_real), INTENT(IN)    :: Bsig(:)                !< Diagonal of the background-error covariance matrix
REAL(kind_real), INTENT(IN)    :: Bm1(:,:)               !< Inverse of the background-error covariance matrix
REAL(kind_real), INTENT(IN)    :: Om1(:,:)               !< Inverse of the observation-error covariance matrix
INTEGER, INTENT(IN)            :: Iter_max               !< Maximum number of iterations
REAL(kind_real), INTENT(IN)    :: Delta                  !< Threshold for change in the cost function for convergence
REAL(kind_real), INTENT(IN)    :: Delta_ct2              !< Threshold for second convergence test (change in state divided by gradient)
LOGICAL, INTENT(IN)            :: GPSRO_pseudo_ops       !< Whether to use pseudo-levels in the vertical interpolation
LOGICAL, INTENT(IN)            :: GPSRO_vert_interp_ops  !< Use log(p) for vertical interpolation?
REAL(kind_real), INTENT(IN)    :: GPSRO_min_temp_grad    !< Threshold for vertical temperature gradient
LOGICAL, INTENT(IN)            :: capsupersat            !< Remove super-saturation?
INTEGER, INTENT(OUT)           :: it                     !< Number of iterations taken
REAL(kind_real), INTENT(OUT)   :: x(:)                   !< Model state at minimum
REAL(kind_real), INTENT(OUT)   :: yb(:)                  !< Calculated refractivity for the background state
REAL(kind_real), INTENT(OUT)   :: ycalc(:)               !< Calculated refractivity for the solution
REAL(kind_real), INTENT(OUT)   :: J_pen                  !< Cost function at solution
REAL(kind_real), INTENT(OUT)   :: Amat(:,:)              !< Analysis-error covariance matrix
LOGICAL, INTENT(OUT)           :: Converged              !< Did the 1DVar run to convergence?
LOGICAL, INTENT(OUT)           :: RefracErr              !< Was there an error in calculating the refractivity?
LOGICAL, INTENT(OUT)           :: Do1DVar_Error          !< Was there an error in the 1DVar step?
REAL(kind_real), INTENT(OUT)   :: O_Bdiff                !< Measure of O-B for whole profile
REAL(kind_real), INTENT(INOUT) :: Tb(nlevq)              !< Calculated background temperature
REAL(kind_real), INTENT(INOUT) :: Ts(nlevq)              !< 1DVar solution temperature
REAL(kind_real), INTENT(INOUT) :: DFS                    !< Measure of degrees of freesom of signal for whole profile

! Local declarations:
INTEGER, PARAMETER           :: max_string = 800
REAL(kind_real), PARAMETER   :: threshold = 50   ! Ignore differences of more then 50 N-units when calculating derivates

LOGICAL                      :: MARQ
INTEGER                      :: i
INTEGER                      :: it_marq
INTEGER                      :: ErrorCode
REAL(kind_real)              :: J_old
REAL(kind_real)              :: J_min
REAL(kind_real)              :: pen_back
REAL(kind_real)              :: pen_ob
REAL(kind_real)              :: xold(nstate)
REAL(kind_real)              :: x_temp(nstate)
REAL(kind_real)              :: xmin(nstate)
REAL(kind_real)              :: ymin(nobs)
REAL(kind_real)              :: lambda
REAL(kind_real)              :: lamp1
REAL(kind_real)              :: d2J_dx2(nstate,nstate)
REAL(kind_real)              :: dJ_dx(nstate)
REAL(kind_real)              :: diag_d2J(nstate)
REAL(kind_real)              :: Kmat(Nobs,nstate)
REAL(kind_real)              :: dx(nstate)
REAL(kind_real)              :: Conv_Test
REAL(kind_real)              :: ct2
REAL(kind_real)              :: ct3
REAL(kind_real)              :: d2                         ! measure of step taken
REAL(kind_real)              :: sdx(nstate)
REAL(kind_real)              :: OK(Nobs,nstate)            ! O^-1 * Kmat
REAL(kind_real)              :: KOK(nstate,nstate)         ! KT *O^-1 *K
REAL(kind_real)              :: AKOK(nstate,nstate)        ! Amat * above
REAL(kind_real)              :: T(nlevq)
REAL(kind_real)              :: pressure(1:nlevp)
REAL(kind_real)              :: humidity(1:nlevq)

REAL(kind_real), ALLOCATABLE :: model_heights(:)
REAL(kind_real), ALLOCATABLE :: refractivity(:)
INTEGER                      :: nRefLevels
CHARACTER(LEN=max_string)    :: message

!--------------
! 1. Initialise
!--------------

x(:) = xb(:)   ! first guess = background
xmin(:) = xb(:)

T(:) = missing_value(T(1))
Tb(:) = missing_value(Tb(1))
Ts(:) = missing_value(Ts(1))
yb(:) = missing_value(yb(1))

! Logicals

Converged = .FALSE.
Do1DVar_Error = .FALSE.
MARQ = .FALSE.

! Set convergence + penalty  values to a large no.

Conv_Test = 1.0E30
ct2 = 1.0E30
ct3 = 1.0E30
J_old = 1.0E30
J_pen = 1.0E30
J_min = 1.0E30

! Set lambda used in the Marqu. Lev. soln -initially a small value

lambda = 1.0E-4_kind_real

it = 0
it_marq = 0  ! no. of Marquardt iterations

O_Bdiff = missing_value(O_Bdiff)
DFS = missing_value(DFS)
ycalc(:) = missing_value(ycalc(1))
ymin(:) = missing_value(ymin(1))
ErrorCode = missing_value(ErrorCode)

!-----------------------
! 2. Main iteration loop
!-----------------------

! Data to stdout on convergence of iteration loop
call oops_log % info('J_pen|Conv_test|ct2|lambda|d2|(dJ/dx)^2|')

Iteration_loop: DO

  IF (it > Iter_max .OR. &
      (Conv_test < Delta  .AND.  &
      ct2 < Delta_ct2 * (Nobs / 200.0_kind_real))) EXIT   ! exit the iteration

  ! Count no. of iterations

  it = it + 1

  ! Calculate the refractivity profile for current guess x

  ! Call the 1D refractivity forward model

  !  1.  First calculate model refractivity on theta levels
  
  ! Unpack the solution to p and q, changing back to SI units
  pressure = 100 * x(1:nlevp)
  humidity = 0.001 * x(nlevp+1:nlevp+nlevq)

  CALL RefMetOffice_ForwardModel(nlevp, &
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
                                 ycalc, &
                                 RefracErr, &
                                 refractivity, &
                                 model_heights)

  ! Don't need model refractivity and heights, so get rid of them
  DEALLOCATE(refractivity, model_heights)

  ! no point proceeding further if ...
  IF (RefracErr) EXIT

  ! Store the refractivity values calculated with `x(:)=xb(:)'

  IF (it == 1) yb(:) = ycalc(:)

  ! For writing the background temperature to a file
  IF (it == 1) Tb(:) = T(:)

  ! Calculate the penalty (cost) function

  CALL Ops_GPSRO_pen (nstate,   &
                      nobs,     &
                      x,        &
                      xb,       &
                      yobs,     &
                      ycalc,    &
                      BM1,      &
                      OM1,      &
                      pen_ob,   &   !out obs penalty
                      pen_back, &   !out back penalty
                      J_pen)

  !store O-B value i.e. 2J/m on first iteration when dx=0
  IF (it == 1) O_Bdiff = 2.0_kind_real * J_pen / REAL (nobs)

  ! store the lowest cost value calculated

  IF (J_min > J_pen) THEN

     J_min = J_pen
     xmin(:) = x(:)
     ymin(:) = ycalc(:)

  END IF

  ! Convergence test

  Conv_test = ABS ((J_old - J_pen) / MIN ( J_old, REAL(NOBS, kind=kind_real)))

  IF (it == 1 .OR. J_pen - J_old < 0.1_kind_real .OR. lambda > 1.0E6_kind_real) THEN

    ! Normal Newtonian iteration

    J_old = J_pen

    IF (.NOT. MARQ) lambda = 0.1_kind_real * lambda

    MARQ = .FALSE.

    ! Evaluate the K matrix for current x
    ! 1.  Calculate the gradient of ref wrt p (on rho levels) and q (on theta levels)
    ! Note: pressure and humidity were unpacked from x earlier in the loop

    CALL jacobian_interface(nlevp, &                  ! Number of pressure levels
                            nlevq, &                  ! Number of specific humidity levels
                            za, &                     ! Heights of the pressure levels
                            zb, &                     ! Heights of the specific humidity levels
                            pressure, &               ! Values of the pressure
                            humidity, &               ! Values of the specific humidity
                            GPSRO_pseudo_ops, &       ! Whether to use pseudo-levels in the calculation
                            GPSRO_vert_interp_ops, &  ! Whether to interpolate using log(pressure)
                            GPSRO_min_temp_grad, &    ! Minimum allowed vertical temperature gradient
                            nobs, &                   ! Number of observations in the profile
                            zobs, &                   ! Height of the observations
                            Kmat)                     ! K-matrix (Jacobian of the observation with respect to the inputs)

    ! Store the state vector in xold

    xold(:) = x(:)

    ! Evaluate the -dJ_dx(vector NOTE SIGN!!) and d2J_dx2(matrix) at x
    CALL Ops_GPSRO_eval_derivs( &
      nstate,    &
      nobs,      &
      x,         &
      xb,        &
      yobs,      &
      ycalc,     &
      BM1,       &
      OM1,       &
      Kmat,      &
      threshold, &  ! Ignore differences larger than 50 N-units
      dJ_dx,     &
      d2J_dx2,   &
      diag_d2J   &
    )

    ! Store inverse of soln. cov matrix
    Amat(:,:) = d2J_dx2(:,:)

  ELSE

    ! M.Lev iteration.The previous increment increased the value
    ! of the penalty function  Use previous values of -dJ_dx, d2J_dx2 and xold

    MARQ = .TRUE.
    it_marq = it_marq + 1
    lambda = 10.0_kind_real * lambda

  END IF

  ! Marq.Lev adjustment to diagonal terms

  lamp1 = lambda + 1.0_kind_real

  DO i = 1,nstate

    d2J_dx2(i,i) = diag_d2J(i)

    ! Marq.Lev. modification to diagonals of Hessian

    d2J_dx2(i,i) = d2J_dx2(i,i) * lamp1

  END DO

  ! Solve the matrix equation
  !
  !  d2J_dx2(:,:) . dx (:) = - dJ_dx (:)
  !
  ! to find dx(:) using the Cholesky decomposition routine
  !

  CALL Ops_Cholesky (d2J_dx2,   &
                     dJ_dx,     &
                     nstate,    &
                     dx,        &   ! The answer, i.e. the "increment"
                     ErrorCode)

  IF (ErrorCode /= 0) EXIT

  ! Update estimate, but limit magnitude of increment with expected background
  ! error

  x(:) = xold(:) + SIGN (MIN (ABS (dx(:)), 2.0_kind_real * Bsig(:), xold(:) / 2.0_kind_real), dx(:))

  ! Check values of humidity -limit to supersat and set <0.0 to = 0.0

  CALL Ops_GPSRO_humidcheck (nstate, &
                             nlevp,  &
                             nlevq,  &
                             za,     &
                             zb,     &
                             capsupersat, &
                             x)

  ! Second convergence criteria in terms of gradient
  ct2 = DOT_PRODUCT (ABS (x(:) - xold(:)), ABS (dJ_dx(:)))
  ct3 = DOT_PRODUCT (dJ_dx(:), dJ_dx(:))

  !-----------------------------------------
  ! Save iteration info to standard output
  !-----------------------------------------
  sdx(:) = MATMUL (Amat(:,:) , (x(:) - xold(:)))  !S^-1.dx
  d2 = DOT_PRODUCT ((x(:) - xold(:)) , Sdx(:))    !d^2=dx(S^-1)dx, size of step normalized by error size

  WRITE (message,'(6E14.6)') J_pen, Conv_test, ct2, lambda, d2, ct3
  call oops_log % info(message)

END DO Iteration_loop

Ts(:) = T(:)                  !1DVAR solution temperature

WRITE (message, '(A,I0)') 'Number of iterations ', it   !write out number of iterations done
call oops_log % info(message)
WRITE (message, '(A,E16.8)') 'O-B size ', O_Bdiff
call oops_log % info(message)

! Output the x(:) that gave the lowest cost function

IF (J_min < J_pen) THEN

  J_pen = J_min
  x(:) = xmin(:)
  ycalc(:) = ymin(:)

END IF

IF (it <= Iter_max) Converged = .TRUE.

! Calculate the solution error cov. matrix

IF (ErrorCode == 0 .AND. it <= iter_max) THEN

  CALL InvertMatrix (nstate,    &
                     nstate,    &
                     Amat(:,:), &
                     ErrorCode)

  ! Check there isn't an error during the matrix inversion

  IF (ErrorCode /= 0) Do1DVar_Error = .TRUE.

  !calculate degrees of freedom for signal
  OK = MATMUL (Om1(:,:) , Kmat(:,:))                  ! O^-1K
  KOK = MATMUL (TRANSPOSE (Kmat(:,:)) , OK(:,:))       ! KTO^-1K
  AKOK = MATMUL (Amat(:,:) , KOK(:,:))                ! AKTO^-1K
  DFS = 0.0                                             ! initialise DFS
  !DFS is equal to trace of AKOK
  DO i = 1, nstate
    DFS = DFS + AKOK(i,i)
  END DO

  WRITE (message,'(A,F16.4)') 'DFS', DFS
  call oops_log % info(message)
ELSE

   Do1DVar_Error = .TRUE.

END IF

END SUBROUTINE Ops_GPSRO_rootsolv

end module ufo_refractivityonedvarcheck_rootsolv_mod
