!-------------------------------------------------------------------------------
! (C) British Crown Copyright 2020 Met Office
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
!-------------------------------------------------------------------------------
module ufo_gnssro_bendmetoffice_tlad_utils_mod

!use iso_c_binding
use fckit_log_module, only: fckit_log
use kinds,            only: kind_real

! Generic routines from elsewhere in jedi
use missing_values_mod
use ufo_constants_mod, only: &
    rd,                     &    ! Gas constant for dry air
    cp,                     &    ! Heat capacity at constant pressure for air
    rd_over_cp,             &    ! Ratio of gas constant to heat capacity
    pref,                   &    ! Reference pressure for calculating exner
    pi,                     &    ! Something to do with circles...
    grav,                   &    ! Gravitational field strength
    ecc,                    &    ! eccentricity
    k_somig,                &    ! Somigliana's constant
    g_equat,                &    ! equatorial gravity (ms-2)
    a_earth,                &    ! semi-major axis of earth (m)
    flatt,                  &    ! flattening
    m_ratio,                &    ! gravity ratio
    mw_ratio,               &    ! Ratio of molecular weights of water and dry air
    c_virtual,              &    ! Related to mw_ratio
    n_alpha,                &    ! Refractivity constant a
    n_beta                       ! Refractivity constant b

implicit none
public             :: Ops_GPSROcalc_alphaK
public             :: Ops_GPSROcalc_nrK
public             :: Ops_GPSRO_refracK
private

contains

!-------------------------------------------------------------------------------
! (C) Crown copyright Met Office. All rights reserved.
!     Refer to COPYRIGHT.txt of this distribution for details.
!-------------------------------------------------------------------------------
! Calculate GPSRO refractivity K matrix.
!-------------------------------------------------------------------------------

SUBROUTINE Ops_GPSRO_refracK (nstate,  &
                              nlevP,   &
                              nb,      &
                              nlevq,   &
                              za,      &
                              zb,      &
                              x,       &
                              GPSRO_pseudo_ops, &
                              GPSRO_vert_interp_ops, &
                              dref_dP, &
                              dref_dq)

IMPLICIT NONE

! Subroutine arguments:
INTEGER, INTENT(IN)                       :: nstate         ! Number of variables in full model state
INTEGER, INTENT(IN)                       :: nlevP          ! Number of p levels in state vector
INTEGER, INTENT(IN)                       :: nb             ! Number of theta levels
INTEGER, INTENT(IN)                       :: nlevq          ! Number of specific humidity levels (=nb)
REAL(kind_real), INTENT(IN)               :: za(:)          ! Heights of the pressure levels
REAL(kind_real), INTENT(IN)               :: zb(:)          ! Heights of the theta (specific humidity) levels
REAL(kind_real), INTENT(IN)               :: x(:)           ! Input model state
LOGICAL, INTENT(IN)                       :: GPSRO_pseudo_ops       ! Whether to use pseudo-levels in the calculation
LOGICAL, INTENT(IN)                       :: GPSRO_vert_interp_ops  ! Whether to use log(p) for the vertical interpolation
REAL(kind_real), ALLOCATABLE, INTENT(OUT) :: dref_dP(:,:)   ! K-matrix for p
REAL(kind_real), ALLOCATABLE, INTENT(OUT) :: dref_dq(:,:)   ! K-matrix for q

! Local declarations:
CHARACTER(len=*), PARAMETER  :: RoutineName = "Ops_GPSRO_refracK"
INTEGER                      :: i
INTEGER                      :: counter
REAL(kind_real)              :: P(nlevP)
REAL(kind_real)              :: Exner(nlevP)
REAL(kind_real)              :: q(nlevq)
REAL(kind_real)              :: Pb(nlevq)
REAL(kind_real)              :: Tv(nlevq)
REAL(kind_real)              :: T(nlevq)
REAL(kind_real)              :: Extheta
REAL(kind_real)              :: pwt1
REAL(kind_real)              :: pwt2
REAL(kind_real)              :: Ndry
REAL(kind_real)              :: Nwet
REAL(kind_real)              :: refrac(nlevq)
REAL(kind_real)              :: dEx_dP(nlevP,nlevP)
REAL(kind_real)              :: dPb_dP(nlevq,nlevP)
REAL(kind_real)              :: dExtheta_dPb(nlevq,nlevq)
REAL(kind_real)              :: dTv_dExtheta(nlevq,nlevq)
REAL(kind_real)              :: dTv_dEx(nlevq,nlevP)
REAL(kind_real)              :: dT_dTv(nlevq,nlevq)
REAL(kind_real)              :: dT_dq(nlevq,nlevq)
REAL(kind_real)              :: dref_dPb(nlevq,nlevq)
REAL(kind_real)              :: dref_dT(nlevq,nlevq)
REAL(kind_real)              :: m1(nb,nlevq)
REAL(kind_real)              :: m2(nb,nlevP)
REAL(kind_real)              :: m3(nb,nlevq)
REAL(kind_real)              :: m4(nb,nlevq)
REAL(kind_real), ALLOCATABLE :: P_pseudo(:)
REAL(kind_real), ALLOCATABLE :: q_pseudo(:)
REAL(kind_real), ALLOCATABLE :: T_pseudo(:)
REAL(kind_real), ALLOCATABLE :: z_pseudo(:)
REAL(kind_real), ALLOCATABLE :: N_pseudo(:)
REAL(kind_real)              :: gamma
REAL(kind_real)              :: beta
REAL(kind_real)              :: c !continuity constant for hydrostatic pressure
REAL(kind_real)              :: g_RB     ! Frequently used term
REAL(kind_real)              :: c_ZZ     ! Frequently used term
REAL(kind_real)              :: dPp_dT1     ! dP_pseudo / dT_below
REAL(kind_real)              :: dPp_dTp     ! dP_pseudo / dT_pseudo
REAL(kind_real)              :: dTp_dT1     ! dT_pseudo / dT_below
REAL(kind_real)              :: dPp_dbeta   ! dP_pseudo / dbeta
REAL(kind_real)              :: dbeta_dT1   ! dbeta / dT_below
REAL(kind_real)              :: dTp_dbeta   ! dT_pseudo / dbeta
REAL(kind_real)              :: dbeta_dT2   ! dbeta / dT_above
REAL(kind_real)              :: dPp_dc      ! dP_pseudo / dc
REAL(kind_real)              :: dc_dT1      ! dc / dT_below
REAL(kind_real)              :: dc_dbeta    ! dc / dbeta
REAL(kind_real)              :: dc_dT2      ! dc / dT_above
REAL(kind_real)              :: dPp_dP1     ! dP_pseudo / dP_below
REAL(kind_real)              :: dc_dP1      ! dc / dP_below
REAL(kind_real)              :: dc_dP2      ! dc / dP_above
REAL(kind_real), ALLOCATABLE :: dref_dPpseudo(:,:)
REAL(kind_real), ALLOCATABLE :: dref_dTpseudo(:,:)
REAL(kind_real), ALLOCATABLE :: dref_dqpseudo(:,:)
REAL(kind_real), ALLOCATABLE :: dPpseudo_dPb(:,:)
REAL(kind_real), ALLOCATABLE :: dTpseudo_dTb(:,:)
REAL(kind_real), ALLOCATABLE :: dqpseudo_dqb(:,:)
REAL(kind_real), ALLOCATABLE :: dPpseudo_dT(:,:)
REAL(kind_real), ALLOCATABLE :: dTb_dP(:,:)
REAL(kind_real), ALLOCATABLE :: dPpseudo_dP(:,:)
REAL(kind_real), ALLOCATABLE :: dTpseudo_dP(:,:)

IF (GPSRO_pseudo_ops) THEN
  ALLOCATE (dref_dP(nb,nlevP))
  ALLOCATE (dref_dq(nb,nlevq))

  ALLOCATE (P_pseudo(nb))
  ALLOCATE (q_pseudo(nb))
  ALLOCATE (T_pseudo(nb))
  ALLOCATE (z_pseudo(nb))
  ALLOCATE (N_pseudo(nb))

  ALLOCATE (dref_dPpseudo(nb,nb))
  ALLOCATE (dref_dTpseudo(nb,nb))
  ALLOCATE (dref_dqpseudo(nb,nb))
  ALLOCATE (dPpseudo_dPb(nb,nlevq))
  ALLOCATE (dPpseudo_dT(nb,nlevq))
  ALLOCATE (dTpseudo_dTb(nb,nlevq))
  ALLOCATE (dqpseudo_dqb(nb,nlevq))

  ALLOCATE (dTb_dP(nlevq,nlevP))
  ALLOCATE (dPpseudo_dP(nb,nlevP))
  ALLOCATE (dTpseudo_dP(nb,nlevP))

  dref_dPpseudo(:,:) = 0.0
  dref_dTpseudo(:,:) = 0.0
  dref_dqpseudo(:,:) = 0.0
  dPpseudo_dPb(:,:) = 0.0
  dPpseudo_dT(:,:) = 0.0
  dTpseudo_dTb(:,:) = 0.0
  dqpseudo_dqb(:,:) = 0.0
  dPpseudo_dP(:,:) = 0.0
  dTpseudo_dP(:,:) = 0.0
ELSE
  ALLOCATE (dref_dP(nlevq,nlevP))
  ALLOCATE (dref_dq(nlevq,nlevq))
END IF

!-----------------------
! 1. Initialise matrices
!-----------------------

dPb_dP(:,:) = 0.0
dExtheta_dPb(:,:) = 0.0
dEx_dP(:,:) = 0.0
dTv_dExtheta(:,:) = 0.0
dTv_dEx(:,:) = 0.0
dT_dTv(:,:) = 0.0
dT_dq(:,:) = 0.0
dref_dpb(:,:) = 0.0
dref_dT(:,:) = 0.0
dref_dq(:,:) = 0.0
dref_dp(:,:) = 0.0

! Set up the P and q vectors from x

P(:) = 1.0E2 * x(1:nlevP)
q(:) = 1.0E-3 * x(nlevP + 1:nstate)

! Calculate exner on rho levels.

Exner(:) = (P(:) / Pref) ** rd_over_cp

DO i = 1,nlevp

  dEx_dP(i,i) = rd_over_cp / Pref * (P(i) / Pref) ** (rd_over_cp - 1.0)

END DO

!----------------------------------------------
! 2. Calculate the refractivity on the b levels
!----------------------------------------------

DO i = 1, nlevq

  ! Calc. pressure on b levels

  pwt1 = (za(i + 1) - zb(i)) / (za(i + 1) - za(i))

  pwt2 = 1.0 - pwt1

  ! calculate the pressure on the theta level.
  IF (GPSRO_vert_interp_ops) THEN
    Pb(i) = EXP (pwt1 * LOG (P(i)) + pwt2 * LOG (P(i + 1)))

    dPb_dP(i,i) = Pb(i) * pwt1 / P(i)
    dPb_dP(i,i + 1) = Pb(i) * pwt2 / P(i + 1)
  ELSE
    ! Assume Exner varies linearly with height
    Pb(i) = Pref * (pwt1 * (P(i) / Pref) ** rd_over_cp + pwt2 * (P(i + 1) / Pref) ** rd_over_cp) ** (1.0 / rd_over_cp)

    dPb_dP(i,i) = pwt1 * (pwt1 * (P(i) / Pref) ** rd_over_cp + pwt2 * &
     (P(i + 1) / Pref) ** rd_over_cp) ** (1.0 / rd_over_cp - 1.0) * (P(i) / Pref) ** (rd_over_cp - 1.0)
    dPb_dP(i,i + 1) = pwt2 * (pwt1 * (P(i) / Pref) ** rd_over_cp + pwt2 * &
     (P(i + 1) / Pref) ** rd_over_cp) ** (1.0 / rd_over_cp - 1.0) * (P(i + 1) / Pref) ** (rd_over_cp-1.0)
  END IF

  ! calculate Exner on the theta level.

  Extheta = (Pb(i) / Pref) ** rd_over_cp

  dExtheta_dPb(i,i) = rd_over_cp * (Pb(i) ** (rd_over_cp - 1.0)) / (Pref ** rd_over_cp)

  ! Calculate mean layer Tv using ND definition

  Tv(i) = grav * (za(i + 1) - za(i)) * Extheta / (Cp * (Exner(i) - Exner(i + 1)))

  dTv_dExtheta(i,i) = Tv(i) / Extheta

  dTv_dEx(i,i) = -Tv(i) / (Exner(i) - Exner(i + 1))

  dTv_dEx(i,i + 1) = Tv(i) / (Exner(i) - Exner(i + 1))

  IF (i > nlevq) THEN

    T(i) = Tv(i)

    dT_dTv(i,i) = 1.0

    ! no wet component

    Nwet = 0.0

  ELSE

    T(i) = Tv(i) / (1.0 + C_virtual * q(i))

    dT_dTv(i,i) = 1.0 / (1.0 + C_virtual * q(i))

    dT_dq(i,i) = -C_virtual * T(i) / (1.0 + C_virtual * q(i))

    ! wet compontent

    Nwet = n_beta * Pb(i) * q(i) / (T(i) ** 2 * (mw_ratio + (1.0 - mw_ratio) * q(i)))

    dref_dq(i,i) = n_beta * Pb(i) * mw_ratio / (T(i) * (mw_ratio + (1.0 - mw_ratio) * q(i))) ** 2

  END IF

  Ndry = n_alpha * Pb(i) / T(i)

  refrac(i) = Ndry + Nwet

  dref_dPb(i,i) = refrac(i) / Pb(i)

  dref_dT(i,i) = -(Ndry + 2.0 * Nwet) / T(i)

END DO

IF (GPSRO_pseudo_ops) THEN
  !----------------------------------!
  !- Add intermediate pseudo-levels -!
  !----------------------------------!
  dTb_dP = MATMUL (dT_dTv, MATMUL (MATMUL (dTv_dExtheta, dExtheta_dPb), dPb_dP) + MATMUL (dTv_dEx, dEx_dP))
  counter = 1
  DO i = 1, nb
    ! Odd 'i' (i.e. copies of actual model level values)
    IF (MOD (i, 2) > 0) THEN
      z_pseudo(i) = zb(counter)
      P_pseudo(i) = Pb(counter)
      q_pseudo(i) = q(counter)
      T_pseudo(i) = T(counter)
      N_pseudo(i) = n_alpha * P_pseudo(i) / T_pseudo(i) + n_beta * P_pseudo(i) * q_pseudo(i) / &
                    (T_pseudo(i) ** 2 * (mw_ratio + (1.0 - mw_ratio) * q_pseudo(i)))

      dref_dPpseudo(i,i) = dref_dPb(counter,counter)
      dref_dTpseudo(i,i) = dref_dT(counter,counter)
      dref_dqpseudo(i,i) = dref_dq(counter,counter)

      dPpseudo_dPb(i,counter) = 1.0
      dTpseudo_dTb(i,counter) = 1.0
      dqpseudo_dqb(i,counter) = 1.0

      counter = counter + 1

    ! Even 'i' (i.e. intermediate pseudo-levels)
    ELSE
      z_pseudo(i) = (zb(counter - 1) + zb(counter)) / 2.0
      IF (MIN (q(counter - 1), q(counter)) > 0.0) THEN
        gamma = LOG (q(counter - 1) / q(counter)) / (zb(counter) - zb(counter - 1))
        q_pseudo(i) = q(counter - 1) * EXP (-gamma * (z_pseudo(i) - z_pseudo(i - 1)))
      ELSE
        q_pseudo(i) = q(counter - 1) + (q(counter) - q(counter - 1)) / &
                     (zb(counter) - zb(counter - 1)) * (z_pseudo(i) - zb(counter - 1))
      END IF

      beta = (T(counter) - T(counter - 1)) / (zb(counter) - zb(counter - 1))
      T_pseudo(i) = T(counter - 1) + beta * (z_pseudo(i) - zb(counter - 1))
      IF (ABS (T(counter) - T(counter - 1)) > 1.0e-10) THEN
        c = ((Pb(counter) / Pb(counter - 1)) * (T(counter) / T(counter - 1)) ** (grav / (rd * beta)) - 1.0) / &
             (zb(counter) - zb(counter - 1))
        P_pseudo(i) = (Pb(counter - 1) * (T_pseudo(i) / T(counter - 1)) ** &
                    (-grav / (rd * beta))) * (1.0 + c * (z_pseudo(i) - zb(counter - 1)))
      ELSE
        P_pseudo(i) = Pb(counter - 1) * EXP (LOG (Pb(counter) / Pb(counter - 1)) * &
                    ((z_pseudo(i) - zb(counter - 1)) / (zb(counter) - zb(counter - 1))))
      END IF

      Ndry = n_alpha * P_pseudo(i) / T_pseudo(i)
      Nwet = n_beta * P_pseudo(i) * q_pseudo(i) / (T_pseudo(i) ** 2 &
                    *(mw_ratio + (1.0 - mw_ratio) * q_pseudo(i)))

      N_pseudo(i) = Ndry + Nwet

      dref_dPpseudo(i,i) = N_pseudo(i) / P_pseudo(i)
      dref_dTpseudo(i,i) = -(Ndry + 2.0 * Nwet) / T_pseudo(i)
      dref_dqpseudo(i,i) = n_beta * P_pseudo(i) * mw_ratio / (T_pseudo(i) * (mw_ratio + &
                           (1.0 - mw_ratio) * q_pseudo(i))) ** 2

      ! Transform P, q and T to pseudo-levels
      IF (MIN (q(counter - 1), q(counter)) > 0.0) THEN
        dqpseudo_dqb(i,counter - 1) = (q_pseudo(i) / q(counter - 1)) * (1.0 - &
             (z_pseudo(i) - zb(counter - 1)) / (zb(counter) - zb(counter - 1)))

        dqpseudo_dqb(i,counter) = (q_pseudo(i) / q(counter)) * ((z_pseudo(i) - zb(counter - 1)) / &
                                  (zb(counter) - zb(counter - 1)))
      ELSE
        dqpseudo_dqb(i,counter - 1) = (zb(counter) - z_pseudo(i)) / (zb(counter) - zb(counter - 1))

        dqpseudo_dqb(i,counter) = (z_pseudo(i) - zb(counter - 1)) / (zb(counter) - zb(counter - 1))
      END IF

      dTpseudo_dTb(i,counter - 1) = (zb(counter) - z_pseudo(i)) / (zb(counter) - zb(counter - 1))

      dTpseudo_dTb(i,counter) = (z_pseudo(i) - zb(counter - 1)) / (zb(counter) - zb(counter - 1))

      ! Items that need changing for dPpseudo/dTb:

      g_RB = grav / (rd * beta)
      c_ZZ = c + 1.0 / (zb(counter) - zb(counter - 1))

      dPp_dT1 = (g_RB / T(counter - 1)) * P_pseudo(i)
      dPp_dTp = -(g_RB / T_pseudo(i)) * P_pseudo(i)
      dTp_dT1 = dTpseudo_dTb(i,counter - 1)
      dPp_dbeta = (g_RB / beta) * LOG (T_pseudo(i) / T(counter - 1)) * P_pseudo(i)
      dbeta_dT1 = -1.0 / (zb(counter) - zb(counter - 1))
      dTp_dbeta = z_pseudo(i) - zb(counter - 1)
      dbeta_dT2 = 1.0 / (zb(counter) - zb(counter - 1))

      IF (ABS (T(counter) - T(counter - 1)) > 1.0e-10) THEN
        ! Incomplete computation of dPpseudo_dPb. Temperature derivatives done below
        dPp_dc = (z_pseudo(i) - zb(counter - 1)) * Pb(counter - 1) * (T_pseudo(i) / T(counter - 1)) ** (-g_RB)
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
        dPpseudo_dPb(i,counter - 1) = EXP (LOG (Pb(counter) / Pb(counter - 1)) * ((z_pseudo(i) - zb(counter - 1)) / &
          (zb(counter) - zb(counter - 1)))) * (1.0 - ((z_pseudo(i) - zb(counter - 1)) / (zb(counter) - zb(counter - 1))))
        dPpseudo_dPb(i,counter) = (Pb(counter - 1) / Pb(counter)) * ((z_pseudo(i) - zb(counter - 1)) / &
          (zb(counter) - zb(counter - 1))) * EXP (LOG (Pb(counter) / Pb(counter - 1)) * ((z_pseudo(i) - zb(counter - 1)) / &
          (zb(counter) - zb(counter - 1))))
      END IF

    END IF

  END DO

  ! temperature derivatives:
  dPpseudo_dP = MATMUL (dPpseudo_dPb, dPb_dp) + MATMUL (dPpseudo_dT, dTb_dP)
  dTpseudo_dP = MATMUL (dTpseudo_dTb, dTb_dP)

  !-------------------------------------------------
  ! 3. Evaluate the Kmatrix by matrix multiplication
  !-------------------------------------------------

  ! calc K matrix for P on rho levels
  dref_dP(:,:) = MATMUL (dref_dPpseudo, dPpseudo_dP) + MATMUL (dref_dTpseudo, dTpseudo_dP)

  ! calc Kmatrix for q on theta levels
  dref_dq(:,:) = MATMUL (dref_dqpseudo, dqpseudo_dqb) + MATMUL (MATMUL (dref_dTpseudo, dTpseudo_dTb), dT_dq) + &
                 MATMUL (MATMUL (dref_dPpseudo, dPpseudo_dT), dT_dq)

  ! Put the K matrices in correct units

  dref_dP(:,:) = 1.0E2 * dref_dP(:,:)   ! hPa

  dref_dq(:,:) = 1.0E-3 * dref_dq(:,:)  ! g/kg

! Normal model levels
ELSE

  !-------------------------------------------------
  ! 3. Evaluate the Kmatrix by matrix multiplication
  !-------------------------------------------------

  ! calc K matrix for P on rho levels
  !  dNmod/dP = (dNmod/dPb * dPb/dP)   + ....
  dref_dP(:,:) = MATMUL (dref_dPb, dPb_dP)

  !  .... (dNmod/dT * dT/dTv * dTv/dEx *dEx/dP) + .....
  m1(:,:) = MATMUL (dref_dT, dT_dTv)
  m2(:,:) = MATMUL (m1, dTv_dEx)
  dref_dP(:,:) = dref_dP(:,:) + MATMUL (m2, dEx_dP)

  !  .... (dNmod/dT * dT/dTv * dTv/dExtheta *dExtheta/dPb*dPb/dP)
  m3(:,:) = MATMUL (m1, dTv_dExtheta)
  m4(:,:) = MATMUL (m3, dExtheta_dPb)
  dref_dP(:,:) = dref_dP(:,:) + MATMUL (m4, dPb_dp)

  ! calc Kmatrix for q on theta levels
  !  dNmod/dq = (dNmod/dq) + (dNmod/dT*dT/dq)
  dref_dq(:,:) = dref_dq(:,:) + MATMUL (dref_dT, dT_dq)

  ! Put the K matrices in correct units

  dref_dP(:,:) = 1.0E2 * dref_dP(:,:)   ! hPa

  dref_dq(:,:) = 1.0E-3 * dref_dq(:,:)  ! g/kg

END IF

IF (ALLOCATED (P_pseudo)) DEALLOCATE (P_pseudo)
IF (ALLOCATED (q_pseudo)) DEALLOCATE (q_pseudo)
IF (ALLOCATED (T_pseudo)) DEALLOCATE (T_pseudo)
IF (ALLOCATED (z_pseudo)) DEALLOCATE (z_pseudo)
IF (ALLOCATED (N_pseudo)) DEALLOCATE (N_pseudo)
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

END SUBROUTINE Ops_GPSRO_refracK


SUBROUTINE Ops_GPSROcalc_nrK (zb,       & ! geopotential heights of model levels
                              nb,       & ! number of levels in zb
                              Rad,      & ! radius of curvature of earth at observation
                              lat,      & ! latitude at observation
                              und,      & ! geoid undulation above WGS-84
                              refrac,   & ! refractivity of model on model levels
                              dnr_dref)   ! Calculated gradient of nr


USE ufo_gnssro_ukmo1d_utils_mod, only: Ops_GPSRO_geop_geom

IMPLICIT NONE

! Subroutine arguments:
INTEGER, INTENT(IN)          :: nb              ! number of levels in zb
REAL(kind_real), INTENT(IN)  :: zb(nb)          ! geopotential heights on zb levels /m
REAL(kind_real), INTENT(IN)  :: Rad             ! local radius of curvature of earth /m
REAL(kind_real), INTENT(IN)  :: lat             ! latitude at observation/ degrees
REAL(kind_real), INTENT(IN)  :: und             ! geoid undulation
REAL(kind_real), INTENT(IN)  :: refrac(nb)      ! refractivity on model levels / N
REAL(kind_real), INTENT(OUT) :: dnr_dref(nb,nb) ! Calculated gradient of nr wrt ref

! Local declarations:
CHARACTER(len=*), PARAMETER :: RoutineName = "Ops_GPSROcalc_nrK"
REAL(kind_real)             :: r(nb)           ! radius of model levels /m
REAL(kind_real)             :: z(nb)           ! geopotential heights on zb levels /m, local copy of zb
INTEGER                     :: i

!Initialise the matrix
dnr_dref = 0.0

!----------------------------------------------
! 1. Convert zb values to geometric altitudes
!---------------------------------------------
z= zb + und                 ! approx. convert to geopotential above WGS-84 ellipsoid
CALL Ops_GPSRO_geop_geom (lat, &
                          z)

!--------------------------------------------------
! 2. Calculate dnr/dref
!--------------------------------------------------

r = Rad + z

DO i = 1,nb
  IF (zb(i) > 0.0 .AND. refrac(i) > 0.0) THEN
    dnr_dref(i,i) = 1.0E-6 * r(i)
  END IF
END DO

END SUBROUTINE Ops_GPSROcalc_nrK

!-------------------------------------------------------------------------------
! GPSRO 1D bending angle operator K code.
!-------------------------------------------------------------------------------
SUBROUTINE Ops_GPSROcalc_alphaK (nobs,     &
                                 nlev,     &
                                 a,        &
                                 refrac,   &
                                 nr,       &
                                 Kmat_ref, &
                                 Kmat_nr)

IMPLICIT NONE

! Subroutine arguments:
INTEGER, INTENT(IN)          :: nobs                ! size of ob. vector
INTEGER, INTENT(IN)          :: nlev                ! no. of refractivity levels
REAL(kind_real), INTENT(IN)  :: a(nobs)             ! observation impact parameters
REAL(kind_real), INTENT(IN)  :: refrac(nlev)        ! refractivity values on model levels
REAL(kind_real), INTENT(IN)  :: nr(nlev)            ! refractive index * radius product
REAL(kind_real), INTENT(OUT) :: Kmat_ref(nobs,nlev) ! BA gradient wrt refractivity
REAL(kind_real), INTENT(OUT) :: kmat_nr(nobs,nlev)  ! BA gradient wrt index * radius product

! Local declarations:
CHARACTER(len=*), PARAMETER :: RoutineName = "Ops_GPSROcalc_alphaK"
INTEGER                     :: i
INTEGER                     :: n
INTEGER                     :: ibot
INTEGER                     :: jbot
INTEGER                     :: kbot
REAL(kind_real)             :: kval(nlev - 1)      ! exponential decay rate between levels
REAL(kind_real)             :: root_2pia
REAL(kind_real)             :: ref_low             ! refractivity at lower level
REAL(kind_real)             :: nr_low              ! impact parameter lower level
REAL(kind_real)             :: tup
REAL(kind_real)             :: tlow                ! upper/lower bounds to error function
REAL(kind_real)             :: dalpha              ! delta bending angle
REAL(kind_real)             :: diff_erf            ! error function result
REAL(kind_real)             :: erf_tup
REAL(kind_real)             :: erf_tlow            ! error function of tup and tlow
REAL(kind_real)             :: t                   ! int. step in approx. to error function
REAL(kind_real), PARAMETER  :: a1 = 0.3480242      ! consts used in error function approx.
REAL(kind_real), PARAMETER  :: a2 = -0.0958798
REAL(kind_real), PARAMETER  :: a3 = 0.7478556
REAL(kind_real), PARAMETER  :: p = 0.47047
REAL(kind_real)             :: dkval_dref(nlev - 1,2)
REAL(kind_real)             :: dkval_dnr(nlev - 1,2)
REAL(kind_real)             :: dalpha_dref(2)
REAL(kind_real)             :: dalpha_dnr(2)
REAL(kind_real)             :: dalpha_dk
REAL(kind_real)             :: dalpha_drlow
REAL(kind_real)             :: dalpha_derf
REAL(kind_real)             :: dalpha_dnrlow
REAL(kind_real)             :: dnrlow_dref(2)
REAL(kind_real)             :: dnrlow_dnr(2)
REAL(kind_real)             :: drlow_dref(2)
REAL(kind_real)             :: drlow_dnr(2)
REAL(kind_real)             :: drlow_dk
REAL(kind_real)             :: derf_dref(2)
REAL(kind_real)             :: derf_dnr(2)
REAL(kind_real)             :: derf_dtup
REAL(kind_real)             :: derf_dtlow
REAL(kind_real)             :: dtup_dnr(2)
REAL(kind_real)             :: dtlow_dnr(2)
REAL(kind_real)             :: dtup_dref(2)
REAL(kind_real)             :: dtlow_dref(2)
REAL(kind_real)             :: dtup_dk
REAL(kind_real)             :: dtlow_dk

!-------------------------------------------------------------------------------
! Initialise the K matrices
!-------------------------------------------------------------------------------

Kmat_ref(:,:) = 0.0
Kmat_nr(:,:) = 0.0
dkval_dref(:,:) = 0.0
dkval_dnr(:,:)  = 0.0

jbot = 1

DO

  IF (refrac(jbot) > 0.0 .AND. nr(jbot) > 0.0) EXIT

  jbot = jbot + 1

END DO

!-------------------------------------------------------------------------------
! Calculate lowest usable level (because of superrefraction)
!-------------------------------------------------------------------------------

kbot = nlev

DO i = nlev,jbot + 1,-1

  ! to avoid large gradients
  IF ((nr(kbot) - nr(kbot-1)) < 10.0) EXIT

  kbot = kbot - 1

END DO

jbot = MAX (jbot,kbot)

!-------------------------------------------------------------------------------
! Calculate the exponential decay rate between levels
!-------------------------------------------------------------------------------

DO i = jbot,nlev - 1

  kval(i) = LOG (refrac(i) / refrac(i + 1)) / &
              MAX (1.0,(nr(i + 1) - nr(i)))

  kval(i) = MAX (1.0E-6,kval(i))

  IF (kval(i) > 1.0E-6) THEN

    dkval_dref(i,1) = 1.0 / (refrac(i) * MAX (1.0,(nr(i + 1) - nr(i))))
    dkval_dref(i,2) = -1.0 / (refrac(i + 1) * MAX (1.0,(nr(i + 1) - nr(i))))

    dkval_dnr(i,1) = kval(i) / MAX (1.0,(nr(i + 1) - nr(i)))
    dkval_dnr(i,2) = -kval(i) / MAX (1.0,(nr(i + 1) - nr(i)))

  END IF

END DO

!-------------------------------------------------------------------------------
! Calculate the bending angle gradients
!-------------------------------------------------------------------------------

DO n = 1,nobs

  IF (a(n) < nr(jbot) .OR. a(n) > nr(nlev)) CYCLE

  Root_2PIa = SQRT (2.0 * pi * a(n))

  ibot = jbot

  ! Find bottom state vector level
  !----------------------------------
  DO

    ! check more than 1 metre apart to stop large gradients in K code
    ! ---------------------------------------------------------------
    IF (((nr(ibot + 1) - a(n)) > 1.0) .OR. ibot == nlev - 1) EXIT

    ibot = ibot + 1

  END DO

  tlow = 0.0

  DO i = ibot, nlev - 1

    ! initialise matrices
    !---------------------

    dalpha_dref(:) = 0.0
    dalpha_dnr(:) = 0.0
    drlow_dref(:) = 0.0
    drlow_dnr(:) = 0.0
    dnrlow_dref(:) = 0.0
    dnrlow_dnr(:) = 0.0
    dtup_dnr(:) = 0.0
    dtup_dref(:) = 0.0
    dtlow_dnr(:) = 0.0
    dtlow_dref(:) = 0.0
    derf_dtup = 0.0
    derf_dtlow = 0.0
    dalpha_drlow = 0.0
    dalpha_dnrlow = 0.0
    dalpha_dk = 0.0
    dalpha_derf = 0.0
    drlow_dk = 0.0
    dtup_dk = 0.0
    dtlow_dk = 0.0

    ! Values of refractivity and impact parameter at lower level
    !-----------------------------------------------------------
    IF (i == ibot) THEN

      ref_low = refrac(i) * EXP (-kval(i) * (a(n) - nr(i)))

      drlow_dref(1)= ref_low / refrac(i)
      drlow_dk  = -ref_low * (a(n) - nr(i))
      drlow_dnr(1) = ref_low * kval(i)

      nr_low = a(n)

      dnrlow_dnr(1) = 0.0

    ELSE

      ref_low = refrac(i)

      drlow_dref(1) = 1.0

      nr_low = nr(i)

      dnrlow_dnr(1) = 1.0

    END IF

    drlow_dref(:) = drlow_dref(:) + drlow_dk * dkval_dref(i,:)
    drlow_dnr(:) = drlow_dnr(:) + drlow_dk * dkval_dnr(i,:)


    ! Limits used in the error function
    !----------------------------------
    IF (i == nlev - 1) THEN

      ! simple extrapolation 100km above the uppermost level.
      !-----------------------------------------------------
      tup = SQRT (kval(i) * (nr(i + 1) + 1.0E5 - a(n)))

      dtup_dk = 0.5 * (nr(i + 1) + 1.0E5 - a(n)) / tup
      dtup_dnr(2) = 0.5 * kval(i) / tup

    ELSE

      tup = SQRT (kval(i) * (nr(i + 1) - a(n)))

      dtup_dk = 0.5 * (nr(i + 1) - a(n)) / tup
      dtup_dnr(2) = 0.5 * kval(i) / tup

    END IF

    dtup_dref(:) = dtup_dref(:) + dtup_dk * dkval_dref(i,:)
    dtup_dnr(:) = dtup_dnr(:) + dtup_dk * dkval_dnr(i,:)

    tlow = 0.0

    IF (i > ibot) THEN
      tlow = SQRT (kval(i) * (nr(i) - a(n)))
      dtlow_dk = 0.5 * (nr(i) - a(n)) / tlow
      dtlow_dnr(1) = 0.5 * kval(i) / tlow
    END IF

    dtlow_dref(:) = dtlow_dref(:) + dtlow_dk * dkval_dref(i,:)
    dtlow_dnr(:) = dtlow_dnr(:) + dtlow_dk * dkval_dnr(i,:)

    ! Abramowitz and Stegun approx. to error function
    !------------------------------------------------
    t = 1.0 / (1.0 + p * tup)
    erf_tup = 1.0 - EXP (-(tup ** 2))  * (a1 + (a2 + a3 * t) * t) * t
    derf_dtup = EXP (-(tup ** 2)) * ((2.0 * tup * (a1 + (a2 + a3 * t) * t) * t) + &
                                 ((a1 + (2.0 * a2 + 3.0 * a3 * t) * t) * p * t ** 2))

    t = 1.0 / (1.0 + p * tlow)
    erf_tlow = 1.0 - EXP (-(tlow ** 2)) * (a1 + (a2 + a3 * t) * t) * t

    ! Multiplied by -1.0 to account for the usage in derf_dref and derf_dnr
    !---------------------------------------------------------------------
    derf_dtlow = -1.0 * EXP (-(tlow ** 2)) * ((2.0 * tlow * (a1 + (a2 + a3 * t) * t) * t) + &
                              ((a1 + (2.0 * a2 + 3.0 * a3 * t) * t) * p * t ** 2))

    diff_erf = erf_tup - erf_tlow

    derf_dref(:) = derf_dtup * dtup_dref(:) + &
                    derf_dtlow * dtlow_dref(:)

    derf_dnr(:) = derf_dtup * dtup_dnr(:) + &
                   derf_dtlow * dtlow_dnr(:)

    dalpha = 1.0E-6 * Root_2PIa * SQRT (kval(i)) * &
                   ref_low * EXP (kval(i) * (nr_low - a(n))) * diff_erf

    dalpha_drlow = dalpha / MAX (1.0E-10,ref_low)

    dalpha_derf = dalpha / MAX (1.0E-10,diff_erf)

    dalpha_dnrlow = dalpha * kval(i)

    dalpha_dk = dalpha * (nr_low - a(n) + 0.5 / kval(i))

    ! Now apply chain rule
    !----------------------

    dalpha_dref(:) = dalpha_dref(:) + &
                     dalpha_drlow * drlow_dref(:) + &
                     dalpha_derf * derf_dref(:) + &
                     dalpha_dnrlow * dnrlow_dref(:) + &
                     dalpha_dk * dkval_dref(i,:)

    dalpha_dnr(:) = dalpha_dnr(:) + &
                    dalpha_drlow * drlow_dnr(:) + &
                    dalpha_derf * derf_dnr(:) + &
                    dalpha_dnrlow * dnrlow_dnr(:) + &
                    dalpha_dk * dkval_dnr(i,:)

    ! Now update matrices
    !---------------------

    Kmat_ref(n,i) = Kmat_ref(n,i) + dalpha_dref(1)
    Kmat_nr(n,i) = Kmat_nr(n,i) + dalpha_dnr(1)

    Kmat_ref(n,i + 1) = Kmat_ref(n,i + 1) + dalpha_dref(2)
    Kmat_nr(n,i + 1) = Kmat_nr(n,i + 1) + dalpha_dnr(2)

  END DO

END DO

END SUBROUTINE Ops_GPSROcalc_alphaK

end module ufo_gnssro_bendmetoffice_tlad_utils_mod
