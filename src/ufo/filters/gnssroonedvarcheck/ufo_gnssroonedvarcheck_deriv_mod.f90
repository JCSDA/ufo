!-------------------------------------------------------------------------------
! (C) British Crown Copyright 2020 Met Office
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
!-------------------------------------------------------------------------------

!> Fortran module for gnssro bending angle Met Office forward operator

module ufo_gnssroonedvarcheck_deriv_mod

use kinds, only: kind_real
use missing_values_mod, only: missing_value

private
public :: Ops_GPSRO_deriv

contains

!-------------------------------------------------------------------------------
! Calculate the simple gradient of the input
!-------------------------------------------------------------------------------
SUBROUTINE Ops_GPSRO_deriv (nobs,   & ! number of levels
                            yb,     & ! background refractivity on obs levels
                            zobs,   & ! obs height levels
                            dyb_dz)   ! gradient

IMPLICIT NONE

! Subroutine arguments:
INTEGER, INTENT(IN)          :: nobs              ! number of levels in yb
REAL(kind_real), INTENT(IN)  :: yb(nobs)          ! input profile to be differentiated, y
REAL(kind_real), INTENT(IN)  :: zobs(nobs)        ! independent co-ordinate for differentiation, x
REAL(kind_real), INTENT(OUT) :: dyb_dz(nobs-1)    ! Calculated dy/dx

! Local declarations:
CHARACTER(len=*), PARAMETER :: RoutineName = "Ops_GPSRO_deriv"
INTEGER                     :: i

! initialise
dyb_dz(:) = missing_value(dyb_dz(1))

! calculate simple gradient delta (y)/delta(x)
DO i = 1,nobs - 1
  IF (yb(i+1) /= missing_value(yb(i+1)) .AND. &
      yb(i) /= missing_value(yb(i)) .AND. &
      zobs(i+1) /= missing_value(zobs(i+1)) .AND. &
      zobs(i) /= missing_value(zobs(i))) THEN
    dyb_dz(i) = (yb(i + 1) - yb(i)) / (zobs(i + 1) - zobs(i))
  END IF
END DO

END SUBROUTINE Ops_GPSRO_deriv

end module ufo_gnssroonedvarcheck_deriv_mod
