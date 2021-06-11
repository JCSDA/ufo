!-------------------------------------------------------------------------------
! (C) British Crown Copyright 2020 Met Office
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
!-------------------------------------------------------------------------------

!> Evaluate the 1st and 2nd deriv. of the cost function.

module ufo_gnssroonedvarcheck_eval_derivs_mod

use kinds, only: kind_real
use missing_values_mod, only: missing_value

private
public :: Ops_GPSRO_eval_derivs_BA

contains

SUBROUTINE Ops_GPSRO_eval_derivs_BA (Nstate,  &  ! size ot state vector
                                     Nobs,    &  ! no of obs
                                     x,       &  ! current est. of soln
                                     xb,      &  ! background vector
                                     yobs,    &  ! obs. vector
                                     ycalc,   &  ! y(x)
                                     BM1,     &  ! inverse .of bsck cov matrix
                                     OM1,     &  ! inv. of obs+forw cov matrix
                                     Kmat,    &  ! gradient matrix
                                     dJ_dx,   &  ! -ve of first deriv. of cost function
                                     d2J_dx2,  & ! second deriv. of cost function.
                                     diag_d2J)   ! vector containing the diagonal values of the matrix above


IMPLICIT NONE

! Subroutine arguments:
INTEGER, INTENT(IN)         :: Nstate
INTEGER, INTENT(IN)         :: Nobs
REAL(kind_real), INTENT(IN)            :: x(:)
REAL(kind_real), INTENT(IN)            :: xb(:)
REAL(kind_real), INTENT(IN)            :: yobs(:)
REAL(kind_real), INTENT(IN)            :: ycalc(:)
REAL(kind_real), INTENT(IN)            :: BM1(:,:)
REAL(kind_real), INTENT(IN)            :: OM1(:,:)
REAL(kind_real), INTENT(IN)            :: Kmat(:,:)
REAL(kind_real), INTENT(OUT)           :: dJ_dx(:)
REAL(kind_real), INTENT(OUT)           :: d2J_dx2(:,:)
REAL(kind_real), INTENT(OUT)           :: diag_d2J(:)

! Local declarations:
CHARACTER(len=*), PARAMETER :: RoutineName = "Ops_GPSRO_eval_derivs_BA"
INTEGER                     :: i
REAL(kind_real)                        :: dx(Nstate)
REAL(kind_real)                        :: dy(Nobs)
REAL(kind_real)                        :: Bdx(Nstate)
REAL(kind_real)                        :: KO(Nstate,Nobs)

!--------------------------------------------------------
! 1. Evaluate the 1st and 2nd deriv. of the cost function
!--------------------------------------------------------

! Deviation from background

dx(:) = x(:) - xb(:)

! Obs. meas-calc

dy(:) = yobs(:) - ycalc(:)

! If absolute difference is greater than 1 radian set difference to 0

WHERE (ABS (dy(:)) > 1.0)

  dy(:) = 0.0

END WHERE

! calc. Bdx matrix   ie B^-1(x-xb)

Bdx(:) = MATMUL (BM1(:,:), dx(:))

! K^T O^-1 matrix

KO(:,:) = MATMUL (TRANSPOSE (Kmat(:,:)), OM1(:,:))

! Calculate -dJ_dx vector    -note the NEGATIVE sign

dJ_dx(:) = MATMUL (KO(:,:), dy(:)) - Bdx(:)
!print*, 'dJ_dx components'
!write(*,'(10E15.8)') dJ_dx(11), Bdx(11)
!write(*,'(10E15.8)') KO(11,:)
!write(*,'(10E15.8)') dy(:)

! d2J_dx2 MATRIX

d2J_dx2(:,:) = MATMUL (KO(:,:), Kmat(:,:)) + BM1(:,:)

! Save the diagonal terms as a vector.

DO i = 1,Nstate

  diag_d2J(i) = d2J_dx2(i,i)

END DO

END SUBROUTINE Ops_GPSRO_eval_derivs_BA

end module ufo_gnssroonedvarcheck_eval_derivs_mod
