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
    rd,                     &    ! Gas constant for dry air
    cp,                     &    ! Heat capacity at constant pressure for air
    rd_over_cp,             &    ! Ratio of gas constant to heat capacity
    pref,                   &    ! Reference pressure for calculating exner
    grav,                   &    ! Gravitational field strength
    n_alpha,                &    ! Refractivity constant a
    n_beta,                  &    ! Refractivity constant b
    mw_ratio,               &    ! Ratio of molecular weights of water and dry air
    C_virtual                    ! Related to mw ratio

implicit none

public :: ufo_refractivity
private

contains


SUBROUTINE ufo_refractivity(nlevq,     &
                            nlevP,     &
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
REAL(kind_real), INTENT(IN)            :: za(:)    ! heights of pressure levels

REAL(kind_real), INTENT(IN)            :: zb(:)             ! Heights of theta levels
REAL(kind_real), INTENT(IN)            :: x(:) ! state vector
REAL(kind_real), INTENT(INOUT)         :: pN(:)        ! pressure on refractivity (theta) levels
LOGICAL, INTENT(OUT)                   :: refracerr    ! errors in refractivity calculation
REAL(kind_real), INTENT(INOUT)         :: refrac(:)    ! refrac on refractivity (theta) levels


! Local declarations:
CHARACTER(len=*), PARAMETER :: RoutineName = "ufo_refractivity"
INTEGER                     :: i
INTEGER                     :: Level
REAL, ALLOCATABLE           :: ExnerN(:)      ! Exner on refractivity (theta) levels
REAL, ALLOCATABLE           :: Exner(:)  ! Exner on pressure (rho) levels
REAL(kind_real)             :: T(nlevq)       ! Temp. on theta levs
REAL                        :: Tv          ! Tv on refractivity (theta) level
REAL                        :: Ndry        ! Dry refractivity
REAL                        :: Nwet        ! Wet refractivity

INTEGER                     :: numPlevs    ! Number of levels of pInter
INTEGER                     :: numqlevs    ! Numbers of levels of qN
LOGICAL                     :: nonmon      ! non-monotonic pressure warning
LOGICAL                     :: unphys      ! zero or negative pressure warning
LOGICAL                     :: levelerr    ! nlevq greater than nlevp error
integer, parameter          :: max_string = 800
CHARACTER(len=max_string)   :: message

REAL(kind_real)             :: P(nlevP)
REAL(kind_real)             :: q(nlevq)
REAL(kind_real)             :: pwt1
REAL(kind_real)             :: pwt2
INTEGER                     :: nstate


! Get constants into right units

nstate = nlevP + nlevq
P(:) = x(1:nlevP)
q(:) = x(nlevP + 1:nstate)

ALLOCATE (ExnerN(nlevq))
ALLOCATE (Exner(nlevP))

refrac(:) = missing_value(refrac(1))
T(:) = missing_value(T(1))
nonmon = .FALSE.
unphys = .FALSE.
levelerr = .FALSE.
refracerr = .FALSE.

DO i = 1, nlevP
  IF (P(i) == missing_value(P(i))) THEN  !pressure missing
    refracerr = .TRUE.
    WRITE(message, *) RoutineName, "Missing value P", i
    CALL fckit_log % warning(message)
    EXIT
  END IF
END DO


DO i = 1, nlevP-1
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

  ! Calculate exner on rho(pressure) levels.

  Exner(:) = (P(:) / pref) ** rd_over_cp

  ! Calculate the refractivity on the b levels


  DO Level = 1, nlevP - 1
	
    pwt1 = (za(Level + 1) - zb(Level)) / (za(Level + 1) - za(Level))
    pwt2 = 1.0 - pwt1

    pN(Level) = EXP (pwt1 * LOG (P(Level)) + pwt2 * LOG (P(Level + 1)))

  END DO

  DO i = 1, nlevq

    ! Calculate Exner on the refractivity level.
    ExnerN(i) = (pN(i) / pref) ** rd_over_cp

    ! Calculate mean layer Tv using ND definition

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

    IF (i > nlevp ) THEN
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

END MODULE ufo_refractivity_utils_mod
