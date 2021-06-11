!-------------------------------------------------------------------------------
! (C) Crown copyright Met Office. All rights reserved.
!     Refer to COPYRIGHT.txt of this distribution for details.
!-------------------------------------------------------------------------------
! Check values of humidity -limit to supersat and set <0.0 to = 0.0
!-------------------------------------------------------------------------------

module ufo_gnssroonedvarcheck_humidcheck_mod

use kinds

private
public :: Ops_GPSRO_humidcheck

contains

SUBROUTINE Ops_GPSRO_humidcheck (nstate, &
                                 nlevp,  &
                                 nlevq,  &
                                 za,     &
                                 zb,     &
                                 capsupersat, &
                                 x)

use ufo_constants_mod, only: &
    rd,                      &    ! Gas constant for dry air
    grav,                    &    ! Gravitational field strength
    c_virtual                     ! Related to mw_ratio

USE ufo_utils_mod, ONLY: &
    Ops_Qsat,            &
    Ops_QsatWat

IMPLICIT NONE

! Subroutine arguments:
INTEGER, INTENT(IN)            :: nstate
INTEGER, INTENT(IN)            :: nlevp
INTEGER, INTENT(IN)            :: nlevq
REAL(kind_real), INTENT(IN)    :: za(:)
REAL(kind_real), INTENT(IN)    :: zb(:)
LOGICAL, INTENT(IN)            :: capsupersat   ! Whether to remove super-saturation (wrt ice?)
REAL(kind_real), INTENT(INOUT) :: x(:)

! Local declarations:
CHARACTER(len=*), PARAMETER :: RoutineName = "Ops_GPSRO_humidcheck"
INTEGER                     :: i
REAL(kind_real)                        :: P(nlevp)
REAL(kind_real)                        :: Q(nlevq)
REAL(kind_real)                        :: T(nlevq)
REAL(kind_real)                        :: pb(nlevq)
REAL(kind_real)                        :: qsaturated(nlevq)
REAL(kind_real)                        :: Tv
REAL(kind_real)                        :: pwt1
REAL(kind_real)                        :: pwt2

!---------------------------------------------------------------------
! 1. Check values of humidity -limit to supersat and set <0.0 to = 0.0
!---------------------------------------------------------------------

! Set up the P and Q vectors from x

P(:) = 1.0E2 * x(1:nlevp)           ! in Pa
Q(:) = 1.0E-3 * x(nlevp + 1:nstate)    ! in kg/kg

DO i = 1, nlevq

  ! Calculate `mean P' for layer

  pwt1 = (za(i + 1) - zb(i)) / (za(i + 1) - za(i))

  pwt2 = 1.0 - pwt1

  pb(i) = EXP (pwt1 * LOG (P(i)) + pwt2 * LOG (P(i + 1)))

  ! Derive the layer mean virtual temp. using the hydrostatic relationship

  Tv = grav * (za(i + 1) - za(i)) / (rd * LOG (P(i) / P(i + 1)))

  ! Calculate the temperature

  T(i) = Tv / (1.0 + C_virtual * Q(i))

END DO

! Calculate the super.sat for T and Pmean

! For T < 0 Ops_Qsat returns the saturated specific humidity over ice.
! Supersaturation with respect to ice is possible, although is largely
! suppressed when CapSupersat is true. When CapSupersat is false the
! humidity check should be done with respect to water.

IF (CapSupersat) THEN
  CALL Ops_Qsat (qsaturated, &  ! out
                 T,          &
                 Pb,         &
                 nlevq)
ELSE
  CALL Ops_QsatWat (qsaturated, &  ! out
                    T,          &
                    Pb,         &
                    nlevq)
END IF

! Check no values have gone -ve
WHERE (x(nlevp + 1:nstate) < 0.0)
  x(nlevp + 1:nstate) = 1.0E-4
END WHERE

! Limit saturated value
WHERE (x(nlevp + 1:nstate) > 1.0E3 * qsaturated(1:nlevq))
  x(nlevp + 1:nstate) = 1.0E3 * qsaturated(1:nlevq)  ! in g/kg
END WHERE

END SUBROUTINE Ops_GPSRO_humidcheck

end module ufo_gnssroonedvarcheck_humidcheck_mod
