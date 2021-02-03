!-------------------------------------------------------------------------------
! (C) British Crown Copyright 2020 Met Office
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
!-------------------------------------------------------------------------------

module ufo_refractivity_utils_mod


!use iso_c_binding
use fckit_log_module, only: fckit_log
use kinds,            only: kind_real

! Generic routines from elsewhere in jedi
use missing_values_mod
use ufo_constants_mod, only: &
    rd,                      &    ! Gas constant for dry air
    cp,                      &    ! Heat capacity at constant pressure for air
    rd_over_cp,              &    ! Ratio of gas constant to heat capacity
    pref,                    &    ! Reference pressure for calculating exner
    grav,                    &    ! Gravitational field strength
    n_alpha,                 &    ! Refractivity constant a
    n_beta,                  &    ! Refractivity constant b
    mw_ratio,                &    ! Ratio of molecular weights of water and dry air
    C_virtual                     ! Related to mw ratio

implicit none

public :: ufo_refractivity
public :: ufo_refractivityDeriv
private

contains


SUBROUTINE ufo_refractivity(nlevq,     &
                            nlevp,     &
                            za,        &
                            zb,        &
                            x,         &
                            pN,        &
                            refracerr, &
                            refrac)
			    
IMPLICIT NONE

! Subroutine arguments:

INTEGER, INTENT(IN)                    :: nlevq        ! no. of levels of wet refractivity required
INTEGER, INTENT(IN)                    :: nlevP        ! no. of levels of dry refractivity required
REAL(kind_real), INTENT(IN)            :: za(:)        ! heights of pressure levels


REAL(kind_real), INTENT(IN)            :: zb(:)        ! Heights of theta levels
REAL(kind_real), INTENT(IN)            :: x(:)         ! state vector
REAL(kind_real), INTENT(INOUT)         :: pN(:)        ! pressure on refractivity (theta) levels
LOGICAL, INTENT(OUT)                   :: refracerr    ! errors in refractivity calculation
REAL(kind_real), INTENT(INOUT)         :: refrac(:)    ! refrac on refractivity (theta) levels


! Local declarations:
CHARACTER(len=*), PARAMETER :: RoutineName = "ufo_refractivity"
INTEGER                     :: i           ! Loop counter
INTEGER                     :: Level       ! Loop counter
REAL, ALLOCATABLE           :: ExnerN(:)   ! Exner on refractivity (theta) levels
REAL, ALLOCATABLE           :: Exner(:)    ! Exner on pressure (rho) levels
REAL(kind_real)             :: T(nlevq)    ! Temp. on theta levs
REAL                        :: Tv          ! Virtual temperature on refractivity (theta) level
REAL                        :: Ndry        ! Dry refractivity
REAL                        :: Nwet        ! Wet refractivity


LOGICAL                     :: nonmon      ! non-monotonic pressure warning
LOGICAL                     :: unphys      ! zero or negative pressure warning
LOGICAL                     :: levelerr    ! nlevq greater than nlevp error
integer, parameter          :: max_string = 800
CHARACTER(len=max_string)   :: message     ! General message for output

REAL(kind_real)             :: P(nlevP)    ! Pressure
REAL(kind_real)             :: q(nlevq)    ! Humidity
REAL(kind_real)             :: pwt1        ! Weighting variable
REAL(kind_real)             :: pwt2        ! Weighting variable
INTEGER                     :: nstate      ! Number of states (P and q)


! Get constants into right units

nstate = nlevp + nlevq
P(:) = x(1:nlevp)
q(:) = x(nlevp + 1:nstate)

ALLOCATE (ExnerN(nlevq))
ALLOCATE (Exner(nlevp))

refrac(:) = missing_value(refrac(1))
T(:) = missing_value(T(1))
nonmon = .FALSE.
unphys = .FALSE.
levelerr = .FALSE.
refracerr = .FALSE.

DO i = 1, nlevp
  IF (P(i) == missing_value(P(i))) THEN  !pressure missing
    refracerr = .TRUE.
    WRITE(message, *) RoutineName, "Missing value P", i
    CALL fckit_log % warning(message)
    EXIT
  END IF
END DO


DO i = 1, nlevp-1
  IF (P(i) - P(i + 1) < 0.0) THEN       !or non-monotonic pressure
    refracerr = .TRUE.
    nonmon = .TRUE.
    WRITE(message,*) "Non monotonic", i, P(i), P(i+1)
    CALL fckit_log % warning(message)
    EXIT
  END IF
END DO

IF (ANY (P(:) <= 0.0)) THEN           !pressure zero or negative
  refracerr = .TRUE.
  unphys = .TRUE.
END IF

IF (nlevq >= nlevp) THEN  ! nlevq must be < than nlevp
  refracerr = .TRUE.
  levelerr = .TRUE.
END IF

! only proceed if pressure is valid
IF (refracerr) THEN
  IF (nonmon) THEN
    CALL fckit_log%warning (RoutineName // ": Pressure non-monotonic")
  ELSE IF (unphys) THEN
    CALL fckit_log%warning (RoutineName // ": Pressure <= zero")
  ELSE
    CALL fckit_log%warning (RoutineName // ": Pressure missing")
  END IF
  IF (levelerr) THEN
    CALL fckit_log%warning (RoutineName // ": Too many wet levels")
  END IF

ELSE

  ! Calculate exner on rho(pressure) levels.

  Exner(:) = (P(:) / pref) ** rd_over_cp

  ! Calculate the refractivity on the b levels


  DO Level = 1, nlevp - 1
	
    pwt1 = (za(Level + 1) - zb(Level)) / (za(Level + 1) - za(Level))
    pwt2 = 1.0 - pwt1

    pN(Level) = EXP (pwt1 * LOG (P(Level)) + pwt2 * LOG (P(Level + 1)))

  END DO

  DO i = 1, nlevq

    ! Calculate Exner on the refractivity level.
    ExnerN(i) = (pN(i) / pref) ** rd_over_cp

    ! Calculate mean layer virtual temp using ND definition

    Tv = grav * (za(i + 1) - za(i)) * ExnerN(i) / &
        (cp * (Exner(i) - Exner(i + 1)))

    IF (i > nlevq) THEN

      T(i) = Tv

      ! No wet component

      Nwet = 0.0

    ELSE

      T(i) = Tv / (1.0 + C_virtual * q(i))

      ! Wet component

      Nwet = n_beta * pN(i) * q(i) / (T(i) ** 2 * (mw_ratio + (1.0 - mw_ratio) * q(i)))


    END IF

    IF (i > nlevP ) THEN
      ! No dry component

      Ndry = 0.0

    ELSE
      ! Dry component

      Ndry = n_alpha * pN(i) / T(i)

    END IF

    refrac(i) = Ndry + Nwet

  END DO

END IF

END SUBROUTINE ufo_refractivity



SUBROUTINE ufo_refractivityDeriv(nlevq,     &
                                 nlevP,     &
                                 za,        &
                                 zb,        &
                                 pN,        &
                                 x,         &
                                 refracerr, &
                                 refrac,    &
                                 dpN_dP,    &
                                 dref_dP,   &
                                 dref_dq)

IMPLICIT NONE

! Subroutine arguments:

INTEGER, INTENT(IN)                    :: nlevq         ! no. of levels of wet refractivity required
INTEGER, INTENT(IN)                    :: nlevP         ! no. of levels of dry refractivity required
REAL(kind_real), INTENT(IN)            :: za(:)         ! heights of pressure (rho) levels
REAL(kind_real), INTENT(IN)            :: zb(:)         ! Heights of theta levels
REAL(kind_real), INTENT(INOUT)         :: pN(:)         ! pressure on refractivity (theta) levels
REAL(kind_real), INTENT(IN)            :: x(:)          ! state vector
LOGICAL, INTENT(OUT)                   :: refracerr     ! errors in refractivity calculation
REAL(kind_real), INTENT(INOUT)         :: refrac(:)     ! refractivity
REAL(kind_real), INTENT(INOUT)         :: dpN_dP(:,:)   ! p_on_theta derivatives wrt P_on_rho levels
REAL(kind_real), INTENT(INOUT)         :: dref_dp(:,:)  ! refrac derivatives wrt P on rho levels
REAL(kind_real), INTENT(INOUT)         :: dref_dq(:,:)  ! refrac derivatites wrt q on theta levels

! Local declarations:
CHARACTER(len=*), PARAMETER            :: RoutineName = "Ops_RefractivityDeriv"
INTEGER                                :: i             ! loop iterator
INTEGER                                :: maxlevs       ! number of levels to loop over
INTEGER                                :: Level         ! loop iterator
REAL(kind_real), ALLOCATABLE           :: ExnerN(:)     ! Exner on refractivity (theta) levels
REAL(kind_real), ALLOCATABLE           :: Exner(:)      ! Exner on pressure (rho) levels
REAL(kind_real)                        :: T(nlevq)      ! Temperature on refractivity (theta) level
REAL(kind_real)                        :: Tv            ! virtual temperature on refractivity (theta) level
REAL(kind_real)                        :: Ndry          ! Dry refractivity
REAL(kind_real)                        :: Nwet          ! Wet refractivity
REAL(kind_real),ALLOCATABLE            :: dEx_dP(:,:)   ! derivatives of exner wrt p on rho levels
REAL(kind_real),ALLOCATABLE            :: dExN_dpN(:,:) ! derivatives of exner wrt p on theta levels
REAL(kind_real),ALLOCATABLE            :: dTv_dEx(:,:)  ! derivatives of Tv wrt exner on rho levels
REAL(kind_real),ALLOCATABLE            :: dTv_dExN(:,:) ! derivatives of Tv wrt exner on theta levels
REAL(kind_real),ALLOCATABLE            :: dT_dTv(:,:)   ! derivatives of T wrt Tv on theta levels
REAL(kind_real),ALLOCATABLE            :: dT_dq(:,:)    ! derivatives of T wrt q on theta levels
REAL(kind_real),ALLOCATABLE            :: dref_dpN(:,:) ! derivatives of refrac wrt p on theta levels
REAL(kind_real),ALLOCATABLE            :: dref_dT(:,:)  ! derivatives of refrac wrt T on theta levels
REAL(kind_real),ALLOCATABLE            :: m1(:,:)       ! Matrix placeholder
REAL(kind_real),ALLOCATABLE            :: m2(:,:)       ! Matrix placeholder
REAL(kind_real),ALLOCATABLE            :: m3(:,:)       ! Matrix placeholder
REAL(kind_real),ALLOCATABLE            :: m4(:,:)       ! Matrix placeholder
LOGICAL                                :: nonmon        ! non-monotonic pressure warning
LOGICAL                                :: unphys        ! zero or negative pressure warning
LOGICAL                                :: levelerr      ! nlevq greater than nlevP error
integer, parameter                     :: max_string = 800
CHARACTER(len=max_string)              :: message       ! General message for output

REAL(kind_real)                        :: P(nlevP)      ! Pressure
REAL(kind_real)                        :: q(nlevq)      ! Humidity
REAL(kind_real)                        :: pwt1          ! Weighting variable
REAL(kind_real)                        :: pwt2          ! Weighting variable
INTEGER                                :: nstate        ! Number of states (P and q)


! Set arrays

nstate = nlevP + nlevq
P(:) = x(1:nlevP)
q(:) = x(nlevP + 1:nstate)

! initialise matrices

maxlevs = MAX (nlevP, nlevq)

ALLOCATE (Exner(nlevP))
ALLOCATE (ExnerN(nlevq))
ALLOCATE (dEx_dP(nlevP,nlevP))
ALLOCATE (dExN_dpN(nlevq,nlevq))
ALLOCATE (dTv_dEx(nlevq,nlevP))
ALLOCATE (dTv_dExN(nlevq,nlevq))
ALLOCATE (dT_dTv(nlevq,nlevq))
ALLOCATE (dT_dq(nlevq,nlevq))
ALLOCATE (dref_dpN(nlevq,nlevq))
ALLOCATE (dref_dT(nlevq,nlevq))
ALLOCATE (m1(nlevq,nlevq))
ALLOCATE (m2(nlevq,nlevP))
ALLOCATE (m3(nlevq,nlevq))
ALLOCATE (m4(nlevq,nlevq))

Exner(:) = 0.0
ExnerN(:) = 0.0
dEx_dP(:,:) = 0.0
dExN_dpN(:,:) = 0.0
dTv_dEx(:,:) = 0.0
dTv_dExN(:,:) = 0.0
dT_dTv(:,:) = 0.0
dT_dq(:,:) = 0.0
dref_dpN(:,:) = 0.0
dref_dT(:,:) = 0.0
dref_dp(:,:) = 0.0
dref_dq(:,:) = 0.0
m1(:,:) = 0.0
m2(:,:) = 0.0
m3(:,:) = 0.0
m4(:,:) = 0.0

T(:) = missing_value(T(1))
nonmon = .FALSE.
unphys = .FALSE.
levelerr = .FALSE.

DO i = 1, nlevP
  IF (P(i) == missing_value(P(i))) THEN  ! pressure missing
    refracerr = .TRUE.
    WRITE(message, *) RoutineName, "Missing value P", i
    CALL fckit_log % warning(message)
    EXIT
  END IF
END DO

DO i=1, nlevP-1
  IF (P(i) - P(i + 1) < 0.0) THEN  ! or non-monotonic pressure
    refracerr = .TRUE.
    nonmon = .TRUE.
    WRITE(message,*) "Non monotonic", i, P(i), P(i+1)
    CALL fckit_log % warning(message)
    EXIT
  END IF
END DO

IF (ANY (P(:) <= 0.0)) THEN       ! pressure zero or negative
  refracerr = .TRUE.
  unphys = .TRUE.
END IF

IF (nlevq >= nlevP) THEN  ! nlevq must be < than nlevP
  refracerr = .TRUE.
  levelerr = .TRUE.
END IF

! only proceed if pressure is valid
IF (refracerr) THEN
  IF (nonmon) THEN
   CALL fckit_log%warning (RoutineName // ": Pressure non-monotonic")
  ELSE IF (unphys) THEN
    CALL fckit_log%warning (RoutineName // ": Pressure <= zero")
  ELSE
    CALL fckit_log%warning (RoutineName // ": Pressure missing")
  END IF
  IF (levelerr) THEN
    CALL fckit_log%warning (RoutineName // ": Too many wet levels")
  END IF
ELSE

  ! Calculate exner on intermediate levels.

  Exner(:) = (P(:) / pref) ** rd_over_cp

  ! Calculate derivative w.r.t P
  DO i = 1, maxlevs

    dEx_dP(i,i) = rd_over_cp * (P(i)) ** (rd_over_cp - 1.0) / (pref ** rd_over_cp)

  END DO

  ! Calculate the refractivity on the b levels

  DO Level = 1, nlevP - 1

    pwt1 = (za(Level + 1) - zb(Level)) / (za(Level + 1) - za(Level))
    pwt2 = 1.0 - pwt1

    pN(Level) = EXP (pwt1 * LOG (P(Level)) + pwt2 * LOG (P(Level + 1)))

    dpN_dP(Level,Level) = pN(Level) * pwt1 / P(Level)
    dpN_dP(Level,Level+1) = pN(Level) * pwt2 / P(Level+1)


  END DO


  DO i = 1, maxlevs - 1

    ! Calculate Exner on the refractivity level.
    ExnerN(i) = (pN(i) / pref) ** rd_over_cp

    ! Calculate derivative w.r.t pN
    dExN_dpN(i,i) = rd_over_cp * (pN(i) ** (rd_over_cp - 1.0)) / (pref ** rd_over_cp)

    ! Calculate mean layer virtual temp using ND definition

    Tv = grav * (za(i + 1) - za(i)) * ExnerN(i) / &
        (cp * (Exner(i) - Exner(i + 1)))

    ! Calculate derivative w.r.t. Exners

    dTv_dExN(i,i) = Tv / ExnerN(i)
    dTv_dEx(i,i) = -Tv / (Exner(i) - Exner(i + 1))
    dTv_dEx(i,i + 1) = Tv / (Exner(i) - Exner(i + 1))

    IF (i > nlevq) THEN

      T(i) = Tv

      ! Calculate derivative w.r.t. Tv

      dT_dTv(i,i) = 1.0

      ! No wet component

      Nwet = 0.0

    ELSE

      T(i) = Tv / (1.0 + C_virtual * q(i))

      ! Calculate derivative w.r.t. Tv and q

      dT_dTv(i,i) = 1.0 / (1.0 + C_virtual * q(i))
      dT_dq (i,i) = -C_virtual * T(i) / (1.0 + C_virtual * q(i))

      ! Wet component

      Nwet = n_beta * pN(i) * q(i) / (T(i) ** 2 * (mw_ratio + (1.0 - mw_ratio) * q(i)))

      ! Calculate derivative with regards to q on theta

      dref_dq(i,i) = n_beta * pN(i) * mw_ratio / (T(i) * (mw_ratio + (1.0 - mw_ratio) * q(i))) ** 2

    END IF

    IF (i > nlevP) THEN

      ! No dry component

      Ndry = 0.0

    ELSE

      ! Dry component

      Ndry = n_alpha * pN(i) / T(i)

   END IF

     refrac(i) = Ndry + Nwet

     ! Calculate derivative w.r.t p on theta and T

     dref_dpN(i,i) = refrac(i) / pN(i)
     dref_dT (i,i) = -(Ndry + 2.0 * Nwet) / T(i)

  END DO

  dref_dp(:,:) = MATMUL (dref_dpN, dpN_dP)

  m1(:,:) = MATMUL (dref_dT, dT_dTv)
  m2(:,:) = MATMUL (m1, dTv_dEx)
  dref_dp(:,:) = dref_dp(:,:) + MATMUL (m2, dEx_dP)

  m3(:,:) = MATMUL (m1, dTv_dExN)
  m4(:,:) = MATMUL (m3, dExN_dpN)
  dref_dp(:,:) = dref_dp(:,:) + MATMUL (m4, dpN_dP)

  dref_dq(:,:) = dref_dq(:,:) + MATMUL (dref_dT, dT_dq)

END IF


DEALLOCATE (Exner)
DEALLOCATE (ExnerN)
DEALLOCATE (dEx_dP)
DEALLOCATE (dExN_dpN)
DEALLOCATE (dTv_dEx)
DEALLOCATE (dTv_dExN)
DEALLOCATE (dT_dTv)
DEALLOCATE (dT_dq)
DEALLOCATE (dref_dpN)
DEALLOCATE (dref_dT)
DEALLOCATE (m1)
DEALLOCATE (m2)
DEALLOCATE (m3)
DEALLOCATE (m4)

END SUBROUTINE ufo_refractivityDeriv

END MODULE ufo_refractivity_utils_mod
