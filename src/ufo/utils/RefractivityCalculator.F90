!-------------------------------------------------------------------------------
! (C) British Crown Copyright 2021 Met Office
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
!> \brief Module for containing a general refractivity forward operator and its K-matrix
!!
!! \author Neill Bowler (Met Office)
!!
!! \date 6 May 2021
!!
!-------------------------------------------------------------------------------

module ufo_utils_refractivity_calculator

use fckit_log_module, only: fckit_log
use missing_values_mod
use ufo_constants_mod, only: &
    rd,                      &    ! Gas constant for dry air
    cp,                      &    ! Heat capacity at constant pressure for air
    rd_over_cp,              &    ! Ratio of gas constant to heat capacity
    pref,                    &    ! Reference pressure for calculating exner
    pi,                      &    ! Something to do with circles...
    grav,                    &    ! Gravitational field strength
    ecc,                     &    ! eccentricity
    k_somig,                 &    ! Somigliana's constant
    g_equat,                 &    ! equatorial gravity (ms-2)
    a_earth,                 &    ! semi-major axis of earth (m)
    flatt,                   &    ! flattening
    m_ratio,                 &    ! gravity ratio
    mw_ratio,                &    ! Ratio of molecular weights of water and dry air
    c_virtual,               &    ! Related to mw_ratio
    n_alpha,                 &    ! Refractivity constant a
    n_beta                        ! Refractivity constant b
use kinds,            only: kind_real

implicit none

private
public :: ufo_calculate_refractivity
public :: ufo_refractivity_kmat
public :: ufo_refractivity_partial_derivatives

contains

!-------------------------------------------------------------------------------
!> \brief Calculation of the refractivity from the pressure and specific humidity
!!
!! \details **ufo_calculate_refractivity**
!! * Checks the inputs that the pressure is not missing, is monotonically
!!   decreasing and not negative.
!! * Perform the calculation on the model temperature levels.  It first calculates
!!   the exner function and virtual temperature from the inputs.  Then using
!!   these calculates the wet and dry components of the refractivity.
!! * If pseudo-level processing is being used, then calculate the temperature,
!!   pressure and specific humidity on the pseudo-levels by interpolation.  Then
!!   calculate the refractivity on all levels.
!!
!! \author Neill Bowler (Met Office)
!!
!! \date 6 May 2021
!!
!-------------------------------------------------------------------------------

SUBROUTINE ufo_calculate_refractivity (nlevP,           &
                                       nlevq,           &
                                       za,              &
                                       zb,              &
                                       P,               &
                                       q,               &
                                       vert_interp_ops, &
                                       pseudo_ops,      &
                                       min_temp_grad,   &
                                       refracerr,       &
                                       nRefLevels,      &
                                       refractivity,    &
                                       model_heights,   &
                                       temperature,     &
                                       interp_pressure, &
                                       tpseudo)

IMPLICIT NONE

! Subroutine arguments:
INTEGER, INTENT(IN)                                   :: nlevP                   !< no. of p levels in state vec.
INTEGER, INTENT(IN)                                   :: nlevq                   !< no. of theta levels
REAL(kind_real), INTENT(IN)                           :: za(nlevp)               !< heights of rho levs
REAL(kind_real), INTENT(IN)                           :: zb(nlevq)               !< heights of theta levs
REAL(kind_real), INTENT(IN)                           :: P(nlevp)                !< state vector
REAL(kind_real), INTENT(IN)                           :: q(nlevq)                !< state vector
LOGICAL, INTENT(IN)                                   :: vert_interp_ops         !< Use log(p) for vertical interpolation?
LOGICAL, INTENT(IN)                                   :: pseudo_ops              !< Use pseudo-levels to reduce errors?
REAL(kind_real), INTENT(IN)                           :: min_temp_grad           !< Minimum value for the vertical temperature gradient
LOGICAL, INTENT(OUT)                                  :: refracerr               !< refractivity error
INTEGER, INTENT(OUT)                                  :: nRefLevels              !< no. of pseudo levs
REAL(kind_real), ALLOCATABLE, INTENT(OUT)             :: refractivity(:)         !< Ref. on pseudo levs
REAL(kind_real), ALLOCATABLE, INTENT(OUT)             :: model_heights(:)        !< height of pseudo levs
REAL(kind_real), OPTIONAL, INTENT(OUT)                :: temperature(nlevq)      !< Calculated temperature on model
REAL(kind_real), OPTIONAL, INTENT(OUT)                :: interp_pressure(nlevq)  !< Model pressure, interpolated to temperature levels
REAL(kind_real), ALLOCATABLE, OPTIONAL, INTENT(INOUT) :: tpseudo(:)              !< Calculated temperature on pseudo levels (or model levels if pseudo_ops is false)

! Local declarations:
CHARACTER(len=*), PARAMETER               :: RoutineName = "ufo_calculate_refractivity"
INTEGER                                   :: i
INTEGER                                   :: counter
REAL(kind_real)                           :: Exner(nlevP)
REAL(kind_real)                           :: Pb(nlevq)
REAL(kind_real)                           :: T_local(nlevq)      ! Temp. on theta levs
REAL(kind_real)                           :: T_virtual
REAL(kind_real)                           :: Ex_theta
REAL(kind_real)                           :: pwt1
REAL(kind_real)                           :: pwt2
REAL(kind_real)                           :: Ndry
REAL(kind_real)                           :: Nwet
REAL(kind_real), ALLOCATABLE              :: P_pseudo(:)
REAL(kind_real), ALLOCATABLE              :: q_pseudo(:)
REAL(kind_real), ALLOCATABLE              :: T_pseudo(:)
REAL(kind_real)                           :: refracModel(1:nlevq)
REAL(kind_real)                           :: gamma
REAL(kind_real)                           :: beta
REAL(kind_real)                           :: c ! continuity constant for hydrostatic pressure
CHARACTER(LEN=200)                        :: message   ! Message to be output to user

! Allocate arrays for pseudo-level processing and output
IF (pseudo_ops) THEN
  nRefLevels = 2 * nlevq - 1
  ALLOCATE (P_pseudo(nRefLevels))
  ALLOCATE (q_pseudo(nRefLevels))
  ALLOCATE (T_pseudo(nRefLevels))
ELSE
  nRefLevels = nlevq
END IF
ALLOCATE (model_heights(nRefLevels))
ALLOCATE (refractivity(nRefLevels))

! Set up the P and q vectors from x

! Initialise refractivity arrays to missing Data
refracModel(:) = missing_value(refracModel(1))
refractivity(:) = missing_value(refractivity(1))
T_local(:) = missing_value(T_local(1))
refracerr = .FALSE.

DO i = 1, nlevq
  IF (P(i) == missing_value(P(i))) THEN  ! pressure missing
    refracerr = .TRUE.
    WRITE(message, *) RoutineName, " Input pressure missing", i
    CALL fckit_log % warning(message)
    EXIT
  END IF

  IF (P(i) - P(i + 1) < 0.0) THEN  ! or non-monotonic pressure
    refracerr = .TRUE.
    WRITE(message, *) RoutineName, " Input pressure non-monotonic", i, P(i), P(i+1)
    CALL fckit_log % warning(message)
    EXIT
  END IF
END DO

IF (ANY (P(:) <= 0.0)) THEN        ! pressure zero or negative
  refracerr = .TRUE.
  WRITE(message, *) RoutineName, " Input pressure not physical"
  CALL fckit_log % warning(message)
END IF

! only proceed if pressure is valid
IF (.NOT. refracerr) THEN

  ! Calculate exner on rho levels.

  Exner(:) = (P(:) / Pref) ** rd_over_cp

  ! Calculate the refractivity on the same model levels as specific humidity
  DO i = 1, nlevq
    ! Calc. pressure pb
    pwt1 = (za(i + 1) - zb(i)) / (za(i + 1) - za(i))

    pwt2 = 1.0 - pwt1

    ! Calculate the pressure on the theta level.
    IF (vert_interp_ops) THEN
      ! Assume ln(P) linear with height
      Pb(i) = EXP (pwt1 * LOG (P(i)) + pwt2 * LOG (P(i + 1)))
    ELSE
      ! Assume Exner varies linearly with height
      Pb(i) = Pref * (pwt1 * (P(i) / Pref) ** rd_over_cp + pwt2 * (P(i + 1) / Pref) ** rd_over_cp) ** (1.0 / rd_over_cp)
    END IF

    ! Calculate the Exner on the theta level.
    Ex_theta = (Pb(i) / Pref) ** rd_over_cp

    ! Calculate mean layer T_virtual on staggered vertical levels

    T_virtual = grav * (za(i + 1) - za(i)) * Ex_theta / &
        (Cp * (Exner(i) - Exner(i + 1)))

    T_local(i) = T_virtual / (1.0 + C_virtual * q(i))

    ! Wet component
    Nwet = n_beta * Pb(i) * q(i) / (T_local(i) ** 2 * (mw_ratio + (1.0 - mw_ratio) * q(i)))

    Ndry = n_alpha * Pb(i) / T_local(i)

    refracModel(i) = Ndry + Nwet

  END DO

  ! Do pseudo-level processing
  IF (pseudo_ops) THEN
    counter = 1
    DO i = 1, nRefLevels

      ! Odd 'i' (i.e. copies of actual model level values)
      IF (MOD (i, 2) > 0) THEN
        model_heights(i) = zb(counter)
        P_pseudo(i) = Pb(counter)
        q_pseudo(i) = q(counter)
        T_pseudo(i) = T_local(counter)
        counter = counter + 1

      ! Even 'i' (i.e. intermediate pseudo-levels)
      ELSE
        model_heights(i) = (zb(counter - 1) + zb(counter)) / 2.0

        ! Assume exponential variation when humidities are positive
        IF (MIN (q(counter - 1), q(counter)) > 0.0) THEN
          gamma = LOG (q(counter - 1) / q(counter)) / (zb(counter) - zb(counter - 1))
          q_pseudo(i) = q(counter - 1) * EXP (-gamma * (model_heights(i) - model_heights(i - 1)))

        ! Assume linear variation if humidities are -ve
        ELSE
          q_pseudo(i) = q(counter - 1) + (q(counter) - q(counter - 1)) / (zb(counter) - &
                        zb(counter - 1)) * (model_heights(i) - zb(counter - 1))
        END IF

        ! T varies linearly with height
        beta = (T_local(counter) - T_local(counter - 1)) / (zb(counter) - zb(counter - 1))
        T_pseudo(i) = T_local(counter - 1) + beta * (model_heights(i) - zb(counter - 1))

        ! Pressure varies to maintain hydrostatic balance
        IF (ABS (T_local(counter) - T_local(counter - 1)) > min_temp_grad) THEN
          c = ((Pb(counter) / Pb(counter - 1)) * (T_local(counter) / T_local(counter - 1)) ** (grav / (rd * beta)) - &
              1.0) / (zb(counter) - zb(counter - 1))
          P_pseudo(i) = (Pb(counter - 1) * (T_pseudo(i) / T_local(counter - 1)) ** &
                      (-grav / (rd * beta))) * (1.0 + c * (model_heights(i) - zb(counter - 1)))
        ELSE
          ! If layer is isothermal, explicitly force P to vary exponentially to avoid singularity
          P_pseudo(i) = Pb(counter - 1) * EXP (LOG (Pb(counter) / Pb(counter - 1)) * &
                      ((model_heights(i) - zb(counter - 1)) / (zb(counter) - zb(counter - 1))))
        END IF
      END IF
    END DO

    refractivity = n_alpha * P_pseudo / T_pseudo + n_beta * P_pseudo * q_pseudo / &
               (T_pseudo ** 2 * (mw_ratio + (1.0 - mw_ratio) * q_pseudo))
  ELSE
    refractivity = refracModel
    model_heights = zb
  END IF

END IF

IF (PRESENT(temperature)) THEN
  temperature(:) = T_local(:)
END IF

IF (PRESENT(tpseudo)) THEN
  IF (pseudo_ops) THEN
    ALLOCATE (tpseudo(1:nRefLevels))
    tpseudo(:) = T_pseudo(:)
  ELSE
    ALLOCATE (tpseudo(1:nlevq))
    tpseudo(:) = T_local(:)
  END IF
END IF

IF (PRESENT(interp_pressure)) THEN
  interp_pressure(:) = Pb(:)
END IF

IF (pseudo_ops) THEN
  IF (ALLOCATED (P_pseudo)) DEALLOCATE (P_pseudo)
  IF (ALLOCATED (q_pseudo)) DEALLOCATE (q_pseudo)
  IF (ALLOCATED (T_pseudo)) DEALLOCATE (T_pseudo)
END IF

END SUBROUTINE ufo_calculate_refractivity

!-------------------------------------------------------------------------------
!> \brief Calculate general refractivity K matrix.
!!
!! \details **ufo_refractivity_kmat**
!! * Allocate intermediate matrices used in calculation.
!! * Calculate the refractivity and basic gradients on the temperature/theta
!!   levels.
!! * If using pseudo-levels, then calculate the partial matrices on these
!!   levels, and use these to calculate the final matrices (dref_dp and dref_dq).
!! * If not using pseudo-levels, then calculate the final matrices.
!! * Deallocate the temporary matrices.
!!
!! \author Neill Bowler (Met Office)
!!
!! \date 6 May 2021
!!
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
!> Calculate general refractivity K matrix.
!-------------------------------------------------------------------------------

SUBROUTINE ufo_refractivity_kmat(nlevP,           &
                                 nlevq,           &
                                 nRefLevels,      &
                                 za,              &
                                 zb,              &
                                 P,               &
                                 q,               &
                                 pseudo_ops,      &
                                 vert_interp_ops, &
                                 min_temp_grad,   &
                                 dref_dP,         &
                                 dref_dq,         &
                                 dPb_dP,          &
                                 refractivity)

IMPLICIT NONE

! Subroutine arguments:
INTEGER, INTENT(IN)                       :: nlevP                     !< no. of pressure levels
INTEGER, INTENT(IN)                       :: nlevq                     !< no. of specific humidity levels
INTEGER, INTENT(IN)                       :: nRefLevels                !< no. of refractivity levels
REAL(kind_real), INTENT(IN)               :: za(nlevp)                 !< Height of the pressure levels
REAL(kind_real), INTENT(IN)               :: zb(nlevq)                 !< Height of the specific humidity levels
REAL(kind_real), INTENT(IN)               :: P(nlevp)                  !< Pressure
REAL(kind_real), INTENT(IN)               :: q(nlevq)                  !< Specific humidity
LOGICAL, INTENT(IN)                       :: pseudo_ops                !< Whether to use pseudo-levels in calculation
LOGICAL, INTENT(IN)                       :: vert_interp_ops           !< Whether to interpolate vertically using exner or ln(p)
REAL(kind_real), INTENT(IN)               :: min_temp_grad             !< Minimum value for the vertical temperature gradient
REAL(kind_real), ALLOCATABLE, INTENT(OUT) :: dref_dP(:,:)              !< kmatrix for p
REAL(kind_real), ALLOCATABLE, INTENT(OUT) :: dref_dq(:,:)              !< kmatrix for q
REAL(kind_real), OPTIONAL, INTENT(OUT)    :: dPb_dP(nlevq,nlevP)       !< Gradient of pressure on theta levels wrt pressure on pressure levels
REAL(kind_real), OPTIONAL, INTENT(OUT), ALLOCATABLE :: refractivity(:) !< Calculated refractivity

! Local declarations:
INTEGER                                   :: i
INTEGER                                   :: counter
REAL(kind_real)                           :: Pb(nlevq)
REAL(kind_real)                           :: T_virtual(nlevq)
REAL(kind_real)                           :: T(nlevq)
REAL(kind_real)                           :: refrac(nlevq)                   ! Refractivity on model levels
REAL(kind_real)                           :: Extheta
REAL(kind_real)                           :: pwt1
REAL(kind_real)                           :: pwt2
REAL(kind_real)                           :: Ndry
REAL(kind_real)                           :: Nwet
REAL(kind_real)                           :: dEx_dP(nlevP,nlevP)
REAL(kind_real)                           :: dPb_dP_local(nlevq,nlevP)
REAL(kind_real)                           :: dExtheta_dPb(nlevq,nlevq)
REAL(kind_real)                           :: dTv_dExtheta(nlevq,nlevq)
REAL(kind_real)                           :: dTv_dEx(nlevq,nlevP)
REAL(kind_real)                           :: dT_dTv(nlevq,nlevq)
REAL(kind_real)                           :: dT_dq(nlevq,nlevq)
REAL(kind_real)                           :: dref_dPb(nlevq,nlevq)
REAL(kind_real)                           :: dref_dT(nlevq,nlevq)
REAL(kind_real)                           :: dref_dq_model(nlevq,nlevq)      ! Gradient of refractivity wrt specific humidity on model levels
REAL(kind_real)                           :: m1(nRefLevels,nlevq)
REAL(kind_real)                           :: m2(nRefLevels,nlevP)
REAL(kind_real)                           :: m3(nRefLevels,nlevq)
REAL(kind_real)                           :: m4(nRefLevels,nlevq)
REAL(kind_real), ALLOCATABLE              :: P_pseudo(:)
REAL(kind_real), ALLOCATABLE              :: q_pseudo(:)
REAL(kind_real), ALLOCATABLE              :: T_pseudo(:)
REAL(kind_real), ALLOCATABLE              :: model_heights(:)
REAL(kind_real)                           :: gamma
REAL(kind_real)                           :: beta
REAL(kind_real)                           :: c           ! continuity constant for hydrostatic pressure
REAL(kind_real)                           :: g_RB        ! Frequently used term
REAL(kind_real)                           :: c_ZZ        ! Frequently used term
REAL(kind_real)                           :: dPp_dT1     ! dP_pseudo / dT_below
REAL(kind_real)                           :: dPp_dTp     ! dP_pseudo / dT_pseudo
REAL(kind_real)                           :: dTp_dT1     ! dT_pseudo / dT_below
REAL(kind_real)                           :: dPp_dbeta   ! dP_pseudo / dbeta
REAL(kind_real)                           :: dbeta_dT1   ! dbeta / dT_below
REAL(kind_real)                           :: dTp_dbeta   ! dT_pseudo / dbeta
REAL(kind_real)                           :: dbeta_dT2   ! dbeta / dT_above
REAL(kind_real)                           :: dPp_dc      ! dP_pseudo / dc
REAL(kind_real)                           :: dc_dT1      ! dc / dT_below
REAL(kind_real)                           :: dc_dbeta    ! dc / dbeta
REAL(kind_real)                           :: dc_dT2      ! dc / dT_above
REAL(kind_real)                           :: dPp_dP1     ! dP_pseudo / dP_below
REAL(kind_real)                           :: dc_dP1      ! dc / dP_below
REAL(kind_real)                           :: dc_dP2      ! dc / dP_above
REAL(kind_real), ALLOCATABLE              :: dref_dPpseudo(:,:)
REAL(kind_real), ALLOCATABLE              :: dref_dTpseudo(:,:)
REAL(kind_real), ALLOCATABLE              :: dref_dqpseudo(:,:)
REAL(kind_real), ALLOCATABLE              :: dPpseudo_dPb(:,:)
REAL(kind_real), ALLOCATABLE              :: dTpseudo_dTb(:,:)
REAL(kind_real), ALLOCATABLE              :: dqpseudo_dqb(:,:)
REAL(kind_real), ALLOCATABLE              :: dPpseudo_dT(:,:)
REAL(kind_real), ALLOCATABLE              :: dTb_dP(:,:)
REAL(kind_real), ALLOCATABLE              :: dPpseudo_dP(:,:)
REAL(kind_real), ALLOCATABLE              :: dTpseudo_dP(:,:)

ALLOCATE (dref_dP(nRefLevels,nlevP))
ALLOCATE (dref_dq(nRefLevels,nlevq))

IF (pseudo_ops) THEN
  ALLOCATE (P_pseudo(nRefLevels))
  ALLOCATE (q_pseudo(nRefLevels))
  ALLOCATE (T_pseudo(nRefLevels))
  ALLOCATE (model_heights(nRefLevels))
  IF (PRESENT(refractivity)) ALLOCATE (refractivity(nRefLevels))

  ALLOCATE (dref_dPpseudo(nRefLevels,nRefLevels))
  ALLOCATE (dref_dTpseudo(nRefLevels,nRefLevels))
  ALLOCATE (dref_dqpseudo(nRefLevels,nRefLevels))
  ALLOCATE (dPpseudo_dPb(nRefLevels,nlevq))
  ALLOCATE (dPpseudo_dT(nRefLevels,nlevq))
  ALLOCATE (dTpseudo_dTb(nRefLevels,nlevq))
  ALLOCATE (dqpseudo_dqb(nRefLevels,nlevq))

  ALLOCATE (dTb_dP(nlevq,nlevP))
  ALLOCATE (dPpseudo_dP(nRefLevels,nlevP))
  ALLOCATE (dTpseudo_dP(nRefLevels,nlevP))

  dref_dPpseudo(:,:) = 0.0
  dref_dTpseudo(:,:) = 0.0
  dref_dqpseudo(:,:) = 0.0
  dPpseudo_dPb(:,:) = 0.0
  dPpseudo_dT(:,:) = 0.0
  dTpseudo_dTb(:,:) = 0.0
  dqpseudo_dqb(:,:) = 0.0
  dPpseudo_dP(:,:) = 0.0
  dTpseudo_dP(:,:) = 0.0
END IF

!-----------------------
! 1. Initialise matrices
!-----------------------

dref_dp(:,:) = 0.0

call ufo_refractivity_partial_derivatives(nlevP,           &
                                          nlevq,           &
                                          za,              &
                                          zb,              &
                                          P,               &
                                          q,               &
                                          vert_interp_ops, &
                                          dT_dTv,          &
                                          dT_dq,           &
                                          dref_dPb,        &
                                          dref_dT,         &
                                          dref_dq_model,   &
                                          refrac,          &
                                          T,               &
                                          Pb,              &
                                          dEx_dP,          &
                                          dExtheta_dPb,    &
                                          dTv_dExtheta,    &
                                          dPb_dP_local,    &
                                          dTv_dEx)

IF (pseudo_ops) THEN
  !----------------------------------!
  !- Add intermediate pseudo-levels -!
  !----------------------------------!
  dTb_dP = MATMUL (dT_dTv, MATMUL (MATMUL (dTv_dExtheta, dExtheta_dPb), dPb_dP_local) + MATMUL (dTv_dEx, dEx_dP))
  counter = 1
  DO i = 1, nRefLevels
    ! Odd 'i' (i.e. copies of actual model level values)
    IF (MOD (i, 2) > 0) THEN
      model_heights(i) = zb(counter)
      P_pseudo(i) = Pb(counter)
      q_pseudo(i) = q(counter)
      T_pseudo(i) = T(counter)
      IF (PRESENT(refractivity)) refractivity(i) = n_alpha * P_pseudo(i) / T_pseudo(i) + &
          n_beta * P_pseudo(i) * q_pseudo(i) / &
          (T_pseudo(i) ** 2 * (mw_ratio + (1.0 - mw_ratio) * q_pseudo(i)))

      dref_dPpseudo(i,i) = dref_dPb(counter,counter)
      dref_dTpseudo(i,i) = dref_dT(counter,counter)
      dref_dqpseudo(i,i) = dref_dq_model(counter,counter)

      dPpseudo_dPb(i,counter) = 1.0
      dTpseudo_dTb(i,counter) = 1.0
      dqpseudo_dqb(i,counter) = 1.0

      counter = counter + 1

    ! Even 'i' (i.e. intermediate pseudo-levels)
    ELSE
      model_heights(i) = (zb(counter - 1) + zb(counter)) / 2.0
      IF (MIN (q(counter - 1), q(counter)) > 0.0) THEN
        gamma = LOG (q(counter - 1) / q(counter)) / (zb(counter) - zb(counter - 1))
        q_pseudo(i) = q(counter - 1) * EXP (-gamma * (model_heights(i) - model_heights(i - 1)))
      ELSE
        q_pseudo(i) = q(counter - 1) + (q(counter) - q(counter - 1)) / &
                     (zb(counter) - zb(counter - 1)) * (model_heights(i) - zb(counter - 1))
      END IF

      beta = (T(counter) - T(counter - 1)) / (zb(counter) - zb(counter - 1))
      T_pseudo(i) = T(counter - 1) + beta * (model_heights(i) - zb(counter - 1))
      IF (ABS (T(counter) - T(counter - 1)) > min_temp_grad) THEN
        c = ((Pb(counter) / Pb(counter - 1)) * (T(counter) / T(counter - 1)) ** (grav / (rd * beta)) - 1.0) / &
             (zb(counter) - zb(counter - 1))
        P_pseudo(i) = (Pb(counter - 1) * (T_pseudo(i) / T(counter - 1)) ** &
                    (-grav / (rd * beta))) * (1.0 + c * (model_heights(i) - zb(counter - 1)))
      ELSE
        P_pseudo(i) = Pb(counter - 1) * EXP (LOG (Pb(counter) / Pb(counter - 1)) * &
                    ((model_heights(i) - zb(counter - 1)) / (zb(counter) - zb(counter - 1))))
      END IF

      Ndry = n_alpha * P_pseudo(i) / T_pseudo(i)
      Nwet = n_beta * P_pseudo(i) * q_pseudo(i) / (T_pseudo(i) ** 2 &
                    *(mw_ratio + (1.0 - mw_ratio) * q_pseudo(i)))

      IF (PRESENT(refractivity)) refractivity(i) = Ndry + Nwet

      dref_dPpseudo(i,i) = (Ndry + Nwet) / P_pseudo(i)
      dref_dTpseudo(i,i) = -(Ndry + 2.0 * Nwet) / T_pseudo(i)
      dref_dqpseudo(i,i) = n_beta * P_pseudo(i) * mw_ratio / (T_pseudo(i) * (mw_ratio + &
                           (1.0 - mw_ratio) * q_pseudo(i))) ** 2

      ! Transform P, q and T to pseudo-levels
      IF (MIN (q(counter - 1), q(counter)) > 0.0) THEN
        dqpseudo_dqb(i,counter - 1) = (q_pseudo(i) / q(counter - 1)) * (1.0 - &
             (model_heights(i) - zb(counter - 1)) / (zb(counter) - zb(counter - 1)))

        dqpseudo_dqb(i,counter) = (q_pseudo(i) / q(counter)) * ((model_heights(i) - zb(counter - 1)) / &
                                  (zb(counter) - zb(counter - 1)))
      ELSE
        dqpseudo_dqb(i,counter - 1) = (zb(counter) - model_heights(i)) / (zb(counter) - zb(counter - 1))

        dqpseudo_dqb(i,counter) = (model_heights(i) - zb(counter - 1)) / (zb(counter) - zb(counter - 1))
      END IF

      dTpseudo_dTb(i,counter - 1) = (zb(counter) - model_heights(i)) / (zb(counter) - zb(counter - 1))

      dTpseudo_dTb(i,counter) = (model_heights(i) - zb(counter - 1)) / (zb(counter) - zb(counter - 1))

      ! Items that need changing for dPpseudo/dTb:

      g_RB = grav / (rd * beta)
      c_ZZ = c + 1.0 / (zb(counter) - zb(counter - 1))

      dPp_dT1 = (g_RB / T(counter - 1)) * P_pseudo(i)
      dPp_dTp = -(g_RB / T_pseudo(i)) * P_pseudo(i)
      dTp_dT1 = dTpseudo_dTb(i,counter - 1)
      dPp_dbeta = (g_RB / beta) * LOG (T_pseudo(i) / T(counter - 1)) * P_pseudo(i)
      dbeta_dT1 = -1.0 / (zb(counter) - zb(counter - 1))
      dTp_dbeta = model_heights(i) - zb(counter - 1)
      dbeta_dT2 = 1.0 / (zb(counter) - zb(counter - 1))

      IF (ABS (T(counter) - T(counter - 1)) > min_temp_grad) THEN
        ! Incomplete computation of dPpseudo_dPb. Temperature derivatives done below
        dPp_dc = (model_heights(i) - zb(counter - 1)) * Pb(counter - 1) * (T_pseudo(i) / T(counter - 1)) ** (-g_RB)
        dc_dT1 = -(g_RB / (T(counter - 1))) * c_ZZ
        dc_dbeta = -(g_RB / beta) * LOG (T(counter) / T(counter - 1)) * c_ZZ
        dc_dT2 = (g_RB / T(counter)) * c_ZZ
        dPp_dP1 = P_pseudo(i) / Pb(counter - 1)
        dPpseudo_dT(i,counter - 1) = dPp_dT1 + dPp_dTp * dTp_dT1 + dPp_dbeta * dbeta_dT1 + dPp_dc * &
                                   (dc_dT1 + dc_dbeta * dbeta_dT1)
        dPpseudo_dT(i,counter) = (dPp_dbeta + dPp_dTp * dTp_dbeta + dPp_dc * dc_dbeta) * &
                                  dbeta_dT2 + dPp_dc * dc_dT2

        dc_dP1 = -(1.0 / Pb(counter - 1)) * c_ZZ
        dc_dP2 =  (1.0 / Pb(counter)) * c_ZZ

        dPpseudo_dPb(i,counter - 1) = dPp_dP1 + dPp_dc * dc_dP1
        dPpseudo_dPb(i,counter) = dPp_dc * dc_dP2
      ELSE
        dPpseudo_dPb(i,counter - 1) = EXP (LOG (Pb(counter) / Pb(counter - 1)) * ((model_heights(i) - zb(counter - 1)) / &
          (zb(counter) - zb(counter - 1)))) * (1.0 - ((model_heights(i) - zb(counter - 1)) / (zb(counter) - zb(counter - 1))))
        dPpseudo_dPb(i,counter) = (Pb(counter - 1) / Pb(counter)) * ((model_heights(i) - zb(counter - 1)) / &
          (zb(counter) - zb(counter - 1))) * EXP (LOG (Pb(counter) / Pb(counter - 1)) * ((model_heights(i) - zb(counter - 1)) / &
          (zb(counter) - zb(counter - 1))))
      END IF

    END IF

  END DO

  ! temperature derivatives:
  dPpseudo_dP = MATMUL (dPpseudo_dPb, dPb_dP_local) + MATMUL (dPpseudo_dT, dTb_dP)
  dTpseudo_dP = MATMUL (dTpseudo_dTb, dTb_dP)

  !-------------------------------------------------
  ! 3. Evaluate the Kmatrix by matrix multiplication
  !-------------------------------------------------

  ! calc K matrix for P on rho levels
  dref_dP(:,:) = MATMUL (dref_dPpseudo, dPpseudo_dP) + MATMUL (dref_dTpseudo, dTpseudo_dP)

  ! calc Kmatrix for q on theta levels
  dref_dq(:,:) = MATMUL (dref_dqpseudo, dqpseudo_dqb) + MATMUL (MATMUL (dref_dTpseudo, dTpseudo_dTb), dT_dq) + &
                 MATMUL (MATMUL (dref_dPpseudo, dPpseudo_dT), dT_dq)

! Normal model levels
ELSE

  if (PRESENT(refractivity)) refractivity = refrac

  !-------------------------------------------------
  ! 3. Evaluate the Kmatrix by matrix multiplication
  !-------------------------------------------------

  ! calc K matrix for P on rho levels
  !  dNmod/dP = (dNmod/dPb * dPb/dP)   + ....
  dref_dP(:,:) = MATMUL (dref_dPb, dPb_dP_local)

  !  .... (dNmod/dT * dT/dTv * dTv/dEx *dEx/dP) + .....
  m1(:,:) = MATMUL (dref_dT, dT_dTv)
  m2(:,:) = MATMUL (m1, dTv_dEx)
  dref_dP(:,:) = dref_dP(:,:) + MATMUL (m2, dEx_dP)

  !  .... (dNmod/dT * dT/dTv * dTv/dExtheta *dExtheta/dPb*dPb/dP)
  m3(:,:) = MATMUL (m1, dTv_dExtheta)
  m4(:,:) = MATMUL (m3, dExtheta_dPb)
  dref_dP(:,:) = dref_dP(:,:) + MATMUL (m4, dPb_dP_local)

  ! calc Kmatrix for q on theta levels
  !  dNmod/dq = (dNmod/dq) + (dNmod/dT*dT/dq)
  dref_dq(:,:) = dref_dq_model(:,:) + MATMUL (dref_dT, dT_dq)

END IF

if (present(dPb_dP)) dPb_dP(:,:) = dPb_dP_local(:,:)

IF (ALLOCATED (P_pseudo)) DEALLOCATE (P_pseudo)
IF (ALLOCATED (q_pseudo)) DEALLOCATE (q_pseudo)
IF (ALLOCATED (T_pseudo)) DEALLOCATE (T_pseudo)
IF (ALLOCATED (model_heights)) DEALLOCATE (model_heights)
IF (ALLOCATED (dref_dPpseudo)) DEALLOCATE (dref_dPpseudo)
IF (ALLOCATED (dref_dTpseudo)) DEALLOCATE (dref_dTpseudo)
IF (ALLOCATED (dref_dqpseudo)) DEALLOCATE (dref_dqpseudo)
IF (ALLOCATED (dPpseudo_dPb)) DEALLOCATE (dPpseudo_dPb)
IF (ALLOCATED (dTpseudo_dTb)) DEALLOCATE (dTpseudo_dTb)
IF (ALLOCATED (dqpseudo_dqb)) DEALLOCATE (dqpseudo_dqb)
IF (ALLOCATED (dPpseudo_dT)) DEALLOCATE (dPpseudo_dT)
IF (ALLOCATED (dTb_dP)) DEALLOCATE (dTb_dP)
IF (ALLOCATED (dPpseudo_dP)) DEALLOCATE (dPpseudo_dP)
IF (ALLOCATED (dTpseudo_dP)) DEALLOCATE (dTpseudo_dP)

END SUBROUTINE ufo_refractivity_kmat

!-------------------------------------------------------------------------------
!> \brief Calculate some partial derivatives of refractivity on model levels
!!
!! \details **ufo_refractivity_partial_derivatives**
!! * Calculate the pressure on model theta levels
!! * Calculate exner on model theta level
!! * Calculate mean layer T_virtual, and then various partial derivatives
!!
!! \author Neill Bowler (Met Office)
!!
!! \date 26 May 2021
!!
!-------------------------------------------------------------------------------

subroutine ufo_refractivity_partial_derivatives(nlevP,           &
                                                nlevq,           &
                                                za,              &
                                                zb,              &
                                                P,               &
                                                q,               &
                                                vert_interp_ops, &
                                                dT_dTv,          &
                                                dT_dq,           &
                                                dref_dPb,        &
                                                dref_dT,         &
                                                dref_dq,         &
                                                refractivity,    &
                                                T,               &
                                                Pb,              &
                                                dEx_dP,          &
                                                dExtheta_dPb,    &
                                                dTv_dExtheta,    &
                                                dPb_dP_local,    &
                                                dTv_dEx)

IMPLICIT NONE

! Subroutine arguments:
INTEGER, INTENT(IN)                       :: nlevP                     !< no. of pressure levels
INTEGER, INTENT(IN)                       :: nlevq                     !< no. of specific humidity levels
REAL(kind_real), INTENT(IN)               :: za(nlevp)                 !< Height of the pressure levels
REAL(kind_real), INTENT(IN)               :: zb(nlevq)                 !< Height of the specific humidity levels
REAL(kind_real), INTENT(IN)               :: P(nlevp)                  !< Pressure
REAL(kind_real), INTENT(IN)               :: q(nlevq)                  !< Specific humidity
LOGICAL, INTENT(IN)                       :: vert_interp_ops           !< Whether to interpolate vertically using exner or ln(p)
REAL(kind_real), INTENT(OUT)              :: dT_dTv(nlevq,nlevq)       !< Partial derivative of temperature wrt virtual temperature
REAL(kind_real), INTENT(OUT)              :: dT_dq(nlevq,nlevq)        !< Partial derivative of temperature wrt specific humidity
REAL(kind_real), INTENT(OUT)              :: dref_dPb(nlevq,nlevq)     !< Partial derivative of refractivity wrt pressure on model theta levels
REAL(kind_real), INTENT(OUT)              :: dref_dT(nlevq,nlevq)      !< Partial derivative of refractivity wrt temperature
REAL(kind_real), INTENT(OUT)              :: dref_dq(nlevq,nlevq)      !< Partial derivative of refractivity wrt specific humidity
REAL(kind_real), INTENT(OUT)              :: refractivity(nlevq)       !< Calculated refractivity
REAL(kind_real), INTENT(OUT)              :: T(nlevq)                  !< Calculated temperature
REAL(kind_real), INTENT(OUT)              :: Pb(nlevq)                 !< Pressure on model theta levels
REAL(kind_real), INTENT(OUT)              :: dEx_dP(nlevP,nlevP)       !< Partial derivative of exner wrt pressure
REAL(kind_real), INTENT(OUT)              :: dExtheta_dPb(nlevq,nlevq) !< Partial derivative of refractivity wrt pressure (at ob location)
REAL(kind_real), INTENT(OUT)              :: dTv_dExtheta(nlevq,nlevq) !< Virtual temperature divided by exner on theta levels
REAL(kind_real), INTENT(OUT)              :: dPb_dP_local(nlevq,nlevP) !< Partial derivative of pressure on theta levels wrt pressure on pressure levels
REAL(kind_real), INTENT(OUT)              :: dTv_dEx(nlevq,nlevP)      !< Partial derivative of virtual temperature wrt exner

! Local declarations:
INTEGER                                   :: i                         ! Loop variable
REAL(kind_real)                           :: Exner(nlevP)              ! Exner on model pressure levels
REAL(kind_real)                           :: T_virtual(nlevq)          ! Virtual temperature
REAL(kind_real)                           :: Extheta                   ! Exner on model theta levels
REAL(kind_real)                           :: pwt1                      ! Weight given to the lower model level in interpolation
REAL(kind_real)                           :: pwt2                      ! Weight given to the upper model level in interpolation
REAL(kind_real)                           :: Ndry                      ! Contribution to refractivity from dry terms
REAL(kind_real)                           :: Nwet                      ! Contribution to refractivity from wet terms

!-----------------------
! 1. Initialise matrices
!-----------------------

dPb_dP_local(:,:) = 0.0
dExtheta_dPb(:,:) = 0.0
dTv_dExtheta(:,:) = 0.0
dTv_dEx(:,:) = 0.0
dT_dTv(:,:) = 0.0
dT_dq(:,:) = 0.0
dref_dpb(:,:) = 0.0
dref_dT(:,:) = 0.0
dref_dq(:,:) = 0.0
dEx_dP(:,:) = 0.0

! Calculate exner on rho levels.
Exner(:) = (P(:) / Pref) ** rd_over_cp
DO i = 1,nlevp
  dEx_dP(i,i) = rd_over_cp / Pref * (P(i) / Pref) ** (rd_over_cp - 1.0)
END DO

!----------------------------------------------
! 2. Calculate the refractivity on the temperature/theta levels
!----------------------------------------------

DO i = 1, nlevq

  pwt1 = (za(i + 1) - zb(i)) / (za(i + 1) - za(i))

  pwt2 = 1.0 - pwt1

  ! calculate the pressure on the theta level.
  IF (vert_interp_ops) THEN
    Pb(i) = EXP (pwt1 * LOG (P(i)) + pwt2 * LOG (P(i + 1)))

    dPb_dP_local(i,i) = Pb(i) * pwt1 / P(i)
    dPb_dP_local(i,i + 1) = Pb(i) * pwt2 / P(i + 1)
  ELSE
    ! Assume Exner varies linearly with height
    Pb(i) = Pref * (pwt1 * (P(i) / Pref) ** rd_over_cp + pwt2 * (P(i + 1) / Pref) ** rd_over_cp) ** (1.0 / rd_over_cp)

    dPb_dP_local(i,i) = pwt1 * (pwt1 * (P(i) / Pref) ** rd_over_cp + pwt2 * &
     (P(i + 1) / Pref) ** rd_over_cp) ** (1.0 / rd_over_cp - 1.0) * (P(i) / Pref) ** (rd_over_cp - 1.0)
    dPb_dP_local(i,i + 1) = pwt2 * (pwt1 * (P(i) / Pref) ** rd_over_cp + pwt2 * &
     (P(i + 1) / Pref) ** rd_over_cp) ** (1.0 / rd_over_cp - 1.0) * (P(i + 1) / Pref) ** (rd_over_cp-1.0)
  END IF

  ! calculate Exner on the theta level.

  Extheta = (Pb(i) / Pref) ** rd_over_cp

  dExtheta_dPb(i,i) = rd_over_cp * (Pb(i) ** (rd_over_cp - 1.0)) / (Pref ** rd_over_cp)

  ! Calculate mean layer T_virtual on staggered vertical levels

  T_virtual(i) = grav * (za(i + 1) - za(i)) * Extheta / (Cp * (Exner(i) - Exner(i + 1)))

  dTv_dExtheta(i,i) = T_virtual(i) / Extheta

  dTv_dEx(i,i) = -T_virtual(i) / (Exner(i) - Exner(i + 1))

  dTv_dEx(i,i + 1) = T_virtual(i) / (Exner(i) - Exner(i + 1))

  T(i) = T_virtual(i) / (1.0 + C_virtual * q(i))

  dT_dTv(i,i) = 1.0 / (1.0 + C_virtual * q(i))

  dT_dq(i,i) = -C_virtual * T(i) / (1.0 + C_virtual * q(i))

  ! wet component

  Nwet = n_beta * Pb(i) * q(i) / (T(i) ** 2 * (mw_ratio + (1.0 - mw_ratio) * q(i)))

  dref_dq(i,i) = n_beta * Pb(i) * mw_ratio / (T(i) * (mw_ratio + (1.0 - mw_ratio) * q(i))) ** 2

  Ndry = n_alpha * Pb(i) / T(i)

  refractivity(i) = Ndry + Nwet

  dref_dPb(i,i) = (Ndry + Nwet) / Pb(i)

  dref_dT(i,i) = -(Ndry + 2.0 * Nwet) / T(i)

END DO

end subroutine ufo_refractivity_partial_derivatives


end module ufo_utils_refractivity_calculator
