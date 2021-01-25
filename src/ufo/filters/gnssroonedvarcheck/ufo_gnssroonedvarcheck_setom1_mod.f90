!-------------------------------------------------------------------------------
! (C) Crown copyright Met Office. All rights reserved.
!     Refer to COPYRIGHT.txt of this distribution for details.
!-------------------------------------------------------------------------------
! Find a solution to the GPSRO inverse problem.
!-------------------------------------------------------------------------------

module ufo_gnssroonedvarcheck_setom1_mod

use kinds
use missing_values_mod
use ufo_roobserror_utils_mod, only: Rmatrix_type

private
public :: Ops_GPSRO_setOM1

contains

SUBROUTINE Ops_GPSRO_setOM1 (nobs,      &
                             zobs,      &
                             yobs,      &
                             Rmatrix,   &
                             OSigma,    &
                             OM1,       &
                             OM1_error)

USE ufo_utils_mod, ONLY: InvertMatrix

IMPLICIT NONE

! Subroutine arguments:
INTEGER, INTENT(IN)             :: nobs
REAL(kind_real), INTENT(IN)     :: zobs(:)
REAL(kind_real), INTENT(IN)     :: yobs(:)
TYPE (rmatrix_type), INTENT(IN) :: Rmatrix
REAL(kind_real), INTENT(OUT)    :: OSigma(:)
REAL(kind_real), INTENT(OUT)    :: OM1(:,:)
LOGICAL, INTENT(OUT)            :: OM1_error

! Local declarations:
CHARACTER(len=*), PARAMETER     :: RoutineName = "Ops_GPSRO_setOM1"
INTEGER                         :: n
INTEGER                         :: i
INTEGER                         :: ReturnCode
REAL(kind_real)                 :: frac_err
REAL(kind_real)                 :: omat(nobs,nobs)
REAL(kind_real), ALLOCATABLE    :: gradient(:)

IF (Rmatrix % satid <= 0) THEN
  ! Rmatrix has not been set up
  OM1_error = .TRUE.
  OSigma(:) = 1
  OM1(:,:) = 1
ELSE
  ALLOCATE (gradient(Rmatrix % num_heights - 1))

  DO i = 1, Rmatrix % num_heights - 1

    ! Calculate the gradient of fractional error with height

    gradient(i) = (Rmatrix % frac_err(i + 1) - Rmatrix % frac_err(i)) / &
                  (Rmatrix % height(i + 1) - Rmatrix % height(i))

  END DO

  OM1_error = .FALSE.

  ! Initialise covariance matrix

  omat(:,:) = 0.0

  ! Calculate the variance values

  DO n = 1, nobs

    i = 1

    DO

      IF (zobs(n) < Rmatrix % height(i + 1) .OR. &
          i + 1 >= Rmatrix % num_heights) THEN
        EXIT
      END IF

      i = i + 1

    END DO

    ! Fractional error

    frac_err = Rmatrix % frac_err(i) + &
               gradient(i) * (zobs(n) - Rmatrix % height(i))

    ! Standard deviation

    OSigma(n) = MAX (frac_err * yobs(n), Rmatrix % min_error)

    ! Variance

    omat(n,n) = OSigma(n) ** 2

  END DO

  ! Calculate the covariances

  DO n = 1,nobs

    DO i = n + 1,nobs

      omat(n,i) = SQRT (omat(n,n) * omat(i,i)) * &
                  EXP (-Rmatrix % clen * ABS (zobs(i) - zobs(n)))

      omat(i,n) = omat(n,i)

    END DO

  END DO

  ! Invert the matrix

  CALL InvertMatrix (nobs,   &
                     nobs,   &
                     omat,   &
                     ReturnCode)

  ! Set the inverse

  OM1(:,:) = omat(:,:)

  IF (ReturnCode /= 0) OM1_error = .TRUE.

  DEALLOCATE (gradient)
END IF

END SUBROUTINE Ops_GPSRO_setOM1

end module ufo_gnssroonedvarcheck_setom1_mod
