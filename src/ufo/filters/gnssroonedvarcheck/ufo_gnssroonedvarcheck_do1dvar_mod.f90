!-------------------------------------------------------------------------------
! (C) British Crown Copyright 2020 Met Office
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
!-------------------------------------------------------------------------------

!> Fortran module for gnssro bending angle Met Office forward operator

module ufo_gnssroonedvarcheck_do1dvar_mod

use kinds, only: kind_real
use missing_values_mod, only: missing_value
use fckit_log_module, only: fckit_log

private
public :: Ops_GPSRO_Do1DVar_BA

contains

!-------------------------------------------------------------------------------
! Find a solution to the satellite sounding inverse problem
!-------------------------------------------------------------------------------
SUBROUTINE Ops_GPSRO_Do1DVar_BA (nlevp,                  &
                                 nlevq,                  &
                                 BM1,                    &
                                 Bsig,                   &
                                 Back,                   &
                                 Ob,                     &
                                 GPSRO_pseudo_ops,       &
                                 GPSRO_vert_interp_ops,  &
                                 GPSRO_Zmin,             &
                                 GPSRO_Zmax,             &
                                 GPSRO_cost_funct_test,  &   ! Threshold value for the cost function convergence test
                                 GPSRO_y_test,           &   ! Threshold value for the yobs-ysol tes
                                 GPSRO_n_iteration_test, &   ! Maximum number of iterations
                                 GPSRO_Delta_factor,     &   ! Delta
                                 GPSRO_Delta_ct2,        &   ! Delta observations
                                 GPSRO_OB_test,          &   ! Threshold value for the O-B test
                                 GPSRO_max_grad,         &   ! Max background vertical N gradient: N units/ metre
                                 GPSRO_hgt_shift,        &   ! Addition to imp_low if max N grad triggered: metre
                                 capsupersat,            &
                                 BAerr,                  &
                                 Tb,                     &
                                 Ts,                     &
                                 O_Bdiff,                &
                                 DFS)

use ufo_gnssroonedvarcheck_utils_mod, only: &
    singlebg_type,             &
    singleob_type

use ufo_gnssro_ukmo1d_utils_mod, only: &
    Ops_GPSRO_refrac

use ufo_gnssroonedvarcheck_deriv_mod, only: &
    Ops_GPSRO_deriv

use ufo_gnssroonedvarcheck_rootsolv_mod, only: &
    Ops_GPSRO_rootsolv_BA

IMPLICIT NONE

! Subroutine arguments:
INTEGER, INTENT(IN)                 :: nlevp
INTEGER, INTENT(IN)                 :: nlevq
REAL(kind_real), INTENT(IN)         :: BM1(:,:)
REAL(kind_real), INTENT(IN)         :: Bsig(:)
TYPE (SingleBg_type), INTENT(INOUT) :: Back
TYPE (SingleOb_type), INTENT(INOUT) :: Ob
LOGICAL, INTENT(IN)                 :: GPSRO_pseudo_ops
LOGICAL, INTENT(IN)                 :: GPSRO_vert_interp_ops
REAL(kind_real), INTENT(IN)         :: GPSRO_Zmin
REAL(kind_real), INTENT(IN)         :: GPSRO_Zmax
REAL(kind_real), INTENT(IN)         :: GPSRO_cost_funct_test
REAL(kind_real), INTENT(IN)         :: GPSRO_y_test
INTEGER, INTENT(IN)                 :: GPSRO_n_iteration_test
REAL(kind_real), INTENT(IN)         :: GPSRO_Delta_ct2
REAL(kind_real), INTENT(IN)         :: GPSRO_Delta_factor
REAL(kind_real), INTENT(IN)         :: GPSRO_OB_test
REAL(kind_real), INTENT(INOUT)      :: GPSRO_max_grad
REAL(kind_real), INTENT(INOUT)      :: GPSRO_hgt_shift
LOGICAL, INTENT(IN)                 :: capsupersat
LOGICAL, INTENT(OUT)                :: BAerr
REAL(kind_real), INTENT(INOUT)      :: Tb(nlevq)
REAL(kind_real), INTENT(INOUT)      :: Ts(nlevq)
REAL(kind_real), INTENT(INOUT)      :: O_Bdiff  ! measure of O-B for whole profile
REAL(kind_real), INTENT(INOUT)      :: DFS      ! measure of degrees of freedom of signal for whole profile

! Local parameters
CHARACTER(len=*), PARAMETER         :: RoutineName = "Ops_GPSRO_Do1DVar_BA"
INTEGER, PARAMETER                  :: max_string = 800

! Local variables
REAL(kind_real), PARAMETER          :: nr_gap = 500.0   ! padding for vertical range, metres
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
REAL(kind_real)                     :: Amat(nlevp+nlevq,nlevp+nlevq)        ! solultion error cov matrix
REAL(kind_real)                     :: Refmodel(nlevq)            ! model refractivity on theta levels
REAL(kind_real)                     :: lower_impact               ! lower impact height range
REAL(kind_real)                     :: upper_impact               ! upper impact height range
REAL(kind_real)                     :: Tdummy(nlevq)              ! to work create dummy
REAL(kind_real), ALLOCATABLE        :: dyb_dz(:)                  ! refractivity vertical gradient array
REAL(kind_real)                     :: imp_low                    ! lower impact height from gradient test
REAL(kind_real), ALLOCATABLE        :: zobs(:)
REAL(kind_real), ALLOCATABLE        :: yobs(:)
REAL(kind_real), ALLOCATABLE        :: yb(:)
REAL(kind_real), ALLOCATABLE        :: ycalc(:)
REAL(kind_real), ALLOCATABLE        :: Om1(:,:)
INTEGER, ALLOCATABLE                :: index_packed(:)
REAL(kind_real), ALLOCATABLE        :: z_pseudo(:)                ! Heights of model and pseudo-levels
REAL(kind_real), ALLOCATABLE        :: N_pseudo(:)                ! Refractivity on model and pseudo_levels
INTEGER                             :: nb_pseudo                  ! Number of levs to calculate ref on
CHARACTER(LEN=max_string)           :: message

REAL(kind_real) :: temp_rad_curv, temp_latitude, temp_undulation

!--------------
! 1. Initialise
!--------------

IF (GPSRO_pseudo_ops) THEN
  ALLOCATE (dyb_dz((2 * nlevq - 1) - 1))
ELSE
  ALLOCATE (dyb_dz(nlevq - 1))
END IF

Do1DVar_error = .FALSE.
BAerr = .FALSE.
ran_iteration = .FALSE.
OM1_error = .FALSE.

! Set all the PGE values to gross error

Ob % BendingAngle(:) % PGEFinal = 1.0

! Set the background vector xb

xb(1:nlevp) = Back % p(:)          ! in h/Pa
xb(nlevp + 1:nlevp+nlevq) = Back % q(:)    ! in g/kg

! Calculate refractivity on theta levels, to find appropriate
! impact height vertical range

CALL Ops_GPSRO_Refrac (nlevp,     &
                       nlevq,     &
                       Back % za, &
                       Back % zb, &
                       xb,        &
                       GPSRO_vert_interp_ops, &
                       GPSRO_pseudo_ops, &
                       BAerr,     &
                       Refmodel,  &
                       Tdummy,    &
                       z_pseudo,  &
                       N_pseudo,  &
                       nb_pseudo)

! Calculate refractivity vertical gradient
! Pseudo levels
IF (GPSRO_pseudo_ops) THEN
  CALL Ops_GPSRO_deriv (nb_pseudo, &    ! no of modelled levels
                        n_pseudo,  &    ! background refractivity on pseudo levels
                        z_pseudo,  &    ! refractivity height levels
                        dyb_dz)         ! output gradient
! Normal model levels
ELSE
  CALL Ops_GPSRO_deriv (nlevq,     &    ! no of modelled levels
                        Refmodel,  &    ! background refractivity on theta levels
                        Back % zb, &    ! refractivity height levels
                        dyb_dz)         ! output gradient
END IF

!PRINT*, 'Refmodel'
!WRITE(*,'(10F14.4)') Refmodel

xb(1:nlevp) = 1.0E-2 * Back % p(:)          ! in h/Pa
xb(nlevp + 1:nlevp+nlevq) = 1.0E3 * Back % q(:)    ! in g/kg

! initialise lower impact for gradient test
imp_low = 0.0

! Check that Obs % GPSRO_max_grad(i) and  Obs % GPSRO_hgt_shift(i) have been
! read in from the stationlist
IF (GPSRO_max_grad == missing_value(GPSRO_max_grad)) THEN

  ! set default value - a value that should never be triggered
  GPSRO_max_grad = -100000.0

END IF

IF (GPSRO_hgt_shift == missing_value(GPSRO_hgt_shift)) THEN

  ! set default value
  GPSRO_hgt_shift = 0.0

END IF

WRITE (message, '(A,F12.2)') 'GPSRO_max_grad', GPSRO_max_grad
CALL fckit_log % info(message)
WRITE (message, '(A,F12.2)') 'GPSRO_hgt_shift', GPSRO_hgt_shift
CALL fckit_log % info(message)

! Perform refractivity vertical gradient check, and assign lower impact height accordingly
! Pseudo levels
IF (GPSRO_pseudo_ops) THEN
  DO i = nb_pseudo - 1, 1, -1
    ! check gradient
    IF (dyb_dz(i) < GPSRO_max_grad .AND.  &
        dyb_dz(i) /= missing_value(dyb_dz(i))) THEN

      ! make sure imp. hght is > that of sharp gradient,(x(i) -RoC) +GPSRO_hgt_shift
      imp_low = (1.0E-6 * n_pseudo(i) + 1.0) * (z_pseudo(i) + Ob % RO_Rad_Curv % value) - &
                            Ob % RO_Rad_Curv % value + GPSRO_hgt_shift

      WRITE(message, '(A,F12.3)') 'large back refractivity gradient found:', dyb_dz(i)
      CALL fckit_log % info(message)

      ! exit the loop once GPSRO_max_grad exceeded
      EXIT
    END IF
  END DO
  ! lower vertical range for BA assimilation, maximum of x(1)-RoC, GPSRO_Zmin or imp_low
  lower_impact = MAX ((1.0E-6 * n_pseudo(1) + 1.0) * (z_pseudo(1) + Ob % RO_Rad_Curv % value) &
                      - Ob % RO_Rad_Curv % value + nr_gap, GPSRO_Zmin, imp_low)

  ! upper vertical range for BA assimilation, minimum of (x(nb)-RoC) or GPSRO_Zmax
  upper_impact = MIN ((1.0E-6 * n_pseudo(nb_pseudo) + 1.0) * (z_pseudo(nb_pseudo) + Ob % RO_Rad_Curv % value) &
                      - Ob % RO_Rad_Curv % value, GPSRO_Zmax)
! Normal model levels
ELSE
  DO i = nlevq - 1, 1, -1
    ! check gradient
    IF (dyb_dz(i) < GPSRO_max_grad .AND.  &
        dyb_dz(i) /= missing_value(dyb_dz(i))) THEN

      !make sure imp. hght is > that of sharp gradient,(x(i) -RoC) +GPSRO_hgt_shift
      imp_low = (1.0E-6 * Refmodel(i) + 1.0) * (Back % zb(i) + Ob % RO_Rad_Curv % value) &
                      - Ob % RO_Rad_Curv % value + GPSRO_hgt_shift

      WRITE(message,'(A,F12.2)') 'large back refractivity gradient found:', dyb_dz(i)
      CALL fckit_log % info(message)

      !exit the loop once GPSRO_max_grad exceeded
      EXIT
    END IF
  END DO

!  print*, Refmodel(1), Back % zb(1), Ob % RO_Rad_Curv % value, nr_gap, GPSRO_Zmin, imp_low

  ! lower vertical range for BA assimilation, maximum of x(1)-RoC, GPSRO_Zmin or imp_low
  lower_impact = MAX ((1.0E-6 * Refmodel(1) + 1.0) * (Back % zb(1) + Ob % RO_Rad_Curv % value) &
                      - Ob % RO_Rad_Curv % value + nr_gap, GPSRO_Zmin, imp_low)

  ! upper vertical range for BA assimilation, minimum of (x(nlevq)-RoC) or GPSRO_Zmax
  upper_impact = MIN ((1.0E-6 * Refmodel(nlevq) + 1.0) * (Back % zb(nlevq) + Ob % RO_Rad_Curv % value) &
                      - Ob % RO_Rad_Curv % value, GPSRO_Zmax)
END IF

WRITE (message, '(A,F12.2)') 'lower_impact', lower_impact
CALL fckit_log % info(message)
WRITE (message, '(A,F12.2)') 'upper_impact', upper_impact
CALL fckit_log % info(message)

! Set size of obs vector used in 1D- Var

nobs = COUNT (Ob % BendingAngle(:) % value /= missing_value(Ob % BendingAngle(1) % value) .AND. & ! not missing bending angle
              Ob % ImpactParam(:) % value /= missing_value(Ob % ImpactParam(1) % value)   .AND. & ! not missing impact parameter
              Ob % ImpactParam(:) % value - Ob % RO_Rad_Curv % value > lower_impact       .AND. & ! within lower impact height cut-off
              Ob % ImpactParam(:) % value - Ob % RO_Rad_Curv % value < upper_impact       .AND. & ! within upper impact height cut-off
              Ob % qc_flags(:) == 0)

WRITE (message, '(A,I0)') 'size of input obs vector ', SIZE (Ob % BendingAngle(:) % value)
CALL fckit_log % info(message)
WRITE (message, '(A,I0)') 'size of packed obs vector ', nobs
CALL fckit_log % info(message)

! Only continue if more than 10 obs in vector, note 10 is a arbitrary value
! and the BAerr is .FALSE.
IF (nobs > 10 .AND. .NOT. BAerr) THEN

  ! calculate an array of indices of the packed elements

  ALLOCATE (index_packed(nobs))                           ! allocate the packed index vector
  index_packed = missing_value(index_packed(1))           ! initialise

  ! Allocate arrays used in 1D- Var after test to stop allocating size nobs=0
  ALLOCATE (om1(nobs,nobs))
  ALLOCATE (yobs(nobs))
  ALLOCATE (zobs(nobs))
  ALLOCATE (yb(nobs))
  ALLOCATE (ycalc(nobs))


  ! Pack observation arrays for valid values
  ! Note: This hard-codes the R-matrix to be diagonal, since that is all that
  ! is currently available in JEDI.  This will need to be revisted once the
  ! full capability is available.

  om1 = 0
  j = 1
  DO i = 1, SIZE (Ob % BendingAngle(:) % Value)
    IF (Ob % BendingAngle(i) % Value /= missing_value(Ob % BendingAngle(i) % Value) .AND. &
        Ob % ImpactParam(i) % value /= missing_value(Ob % ImpactParam(i) % value)   .AND. &
        Ob % ImpactParam(i) % value - Ob % RO_Rad_Curv % value > lower_impact       .AND. &
        Ob % ImpactParam(i) % value - Ob % RO_Rad_Curv % value < upper_impact       .AND. &
        Ob % qc_flags(i) == 0) THEN
      index_packed(j) = i
      zobs(j) = Ob % ImpactParam(i) % value
      yobs(j) = Ob % BendingAngle(i) % Value
      om1(j,j) = (Ob % BendingAngle(i) % oberr)**(-2)
      j = j + 1
    END IF
  END DO

  !-----------------------------------------------
  ! 2. If no errors so far, call the 1DVar routine
  !-----------------------------------------------

  IF (ALL(zobs(:) /= missing_value(zobs(1))) .AND.  &
      ALL(yobs(:) /= missing_value(yobs(1))) .AND.  &
      .NOT. OM1_error) THEN

    temp_rad_curv = Ob % RO_Rad_Curv % Value
    temp_latitude = Ob % Latitude
    temp_undulation = Ob % RO_geoid_und % value
    CALL Ops_GPSRO_rootsolv_BA (nlevp+nlevq,               &    ! size of state vector
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
                                BAerr,                     &    ! error flag
                                Do1DVar_error,             &    ! error flag
                                GPSRO_n_iteration_test,    &    !
                                GPSRO_Delta_factor,        &    !
                                GPSRO_Delta_ct2,           &    !
                                GPSRO_pseudo_ops,          &
                                GPSRO_vert_interp_ops,     &
                                capsupersat,               &
                                O_Bdiff,                   &    ! observed -background BA value
                                temp_rad_curv,             &    ! Radius of curvature of ellipsoid
                                temp_latitude,             &    ! Latitude of occ
                                temp_undulation,           &    ! geoid undulation
                                Tb,                        &
                                Ts,                        &
                                DFS)
    ran_iteration = .TRUE.
  ELSE

    Do1DVar_error = .TRUE.

    Ob % BendingAngle(:) % PGEFinal = 0.99

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

    ! Bending angle calculated with solution

    Ob % SolutBendingAngle(index_packed) = ycalc(1:nobs)

!    ! store the bending angle calculated from the background

!    Back % BendingAngle(index_packed) = yb(1:nobs)

    ! PROBABILTY OF GROSS ERROR. Use the cost function value.

    IF (J_pen > GPSRO_cost_funct_test * REAL(NOBS, kind=kind_real)) THEN

      ! the cost function for profile is too high - GROSS ERROR -
      ! set all pge's in profile to 0.8

      Ob % BendingAngle(:) % PGEFinal = 0.8

    ELSE

      ! For each value in profile, estimate probability of gross error
      ! from difference between solution and observed value

      DO i = 1,nobs

        IF (ABS (Ob % BendingAngle(index_packed(i)) % Value - Ob % SolutBendingAngle(index_packed(i))) < &
              (GPSRO_y_test * Ob % BendingAngle(index_packed(i)) % ObErr)) THEN
          Ob % BendingAngle(index_packed(i)) % PGEFinal = 0.1
        ELSE
          Ob % BendingAngle(index_packed(i)) % PGEFinal = 0.7
        END IF

      END DO

    END IF

    ! if the initial 2J/m value exceeds read-in value then flag
    IF  (O_Bdiff > GPSRO_OB_test) THEN
      Ob % BendingAngle(:) % PGEFinal = 0.6
    END IF

    ! check for the BAerr being set
    IF (BAerr) THEN
      Ob % BendingAngle(:) % PGEFinal = 0.58     ! flag BAerr
    END IF

  ELSE  ! Do 1D-Var_error

    IF (ran_iteration) THEN
      Ob % BendingAngle(:) % PGEFinal = 0.9     ! flag lack of convergence

      ! Bending angle calculated with solution
      Ob % SolutBendingAngle(index_packed) = ycalc(1:nobs)

      ! map the retrieval information back into the ob structures
      Ob % p(:) % Value = 1.0E2 * x(1:nlevp)
      Ob % q(:) % Value = 1.0E-3 * x(nlevp + 1:nlevp+nlevq)

      ! store iterations and cost
      Ob % Niter = it
      Ob % Jcost = 2.0_kind_real * J_pen / nobs

    ELSE IF (OM1_error) THEN

      Ob % BendingAngle(:) % PGEFinal = 0.85    ! Flag error in getting
                                                ! observation error inverse
    END IF

  END IF

ELSE
  IF (nobs <= 10) THEN
    WRITE (message, '(A)') 'nobs is less than 10: exit Ops_GPSRO_Do1DVar_BA'
    CALL fckit_log % info(message)
    Ob % BendingAngle(:) % PGEFinal = 0.55     ! flag lack of observation data
  END IF

  IF (BAerr) THEN
    WRITE (message, '(A)') 'Error in Ops_Refractivity: exit Ops_GPSRO_Do1DVar_BA'
    CALL fckit_log % info(message)
    Ob % BendingAngle(:) % PGEFinal = 0.58     ! flag BAerr
  END IF
END IF

END SUBROUTINE Ops_GPSRO_Do1DVar_BA

end module ufo_gnssroonedvarcheck_do1dvar_mod
