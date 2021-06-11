!-------------------------------------------------------------------------------
! (C) British Crown Copyright 2020 Met Office
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
!-------------------------------------------------------------------------------

!> Calculates GPSRO penalty.

module ufo_gnssroonedvarcheck_pen_mod

use kinds, only: kind_real
use missing_values_mod, only: missing_value

private
public :: Ops_GPSRO_pen

contains

SUBROUTINE Ops_GPSRO_pen (Nstate,   &   ! size of state vec.
                          Nobs,     &   ! size of obs vec.
                          x,        &   ! current estimate of soltution
                          xb,       &   ! background
                          yobs,     &   ! observed values
                          ycalc,    &   ! y(x)
                          BM1,      &   ! inverse .of bsck cov matrix
                          OM1,      &   ! inv. of obs+forw cov matrix
                          pen_ob,   &   ! obs penalty
                          pen_back, &   ! back penalty
                          pen_func)     ! total penalty


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
REAL(kind_real),  INTENT(OUT)          :: pen_ob
REAL(kind_real),  INTENT(OUT)          :: pen_back
REAL(kind_real),  INTENT(OUT)          :: pen_func

! Local declarations:
REAL(kind_real)                        :: dx(Nstate)
REAL(kind_real)                        :: dy(Nobs)
REAL(kind_real)                        :: Bdx(Nstate)
REAL(kind_real)                        :: Ody(Nobs)
REAL(kind_real)                        :: J_back
REAL(kind_real)                        :: J_obs
CHARACTER(len=*), PARAMETER :: RoutineName = "Ops_GPSRO_pen"

! deviation from background

dx(:) = x(:) - xb(:)

! calc. Bdx matrix   ie B^-1(x-xb)

Bdx(:) = MATMUL (BM1(:,:), dx(:))

! background term (scalar)

J_back = DOT_PRODUCT (dx(:), Bdx(:))

! obs. meas-calc

dy(:) = yobs(:) - ycalc(:)

!make sure missing data is not included
WHERE (ycalc(:) == missing_value(ycalc(1)) .OR. &
       yobs(:) == missing_value(yobs(1)))

  dy(:) = 0.0
END WHERE

! Ody   O^-1 (ymeas-ycalc)

Ody(:) = MATMUL (OM1(:,:), dy (:))

! observation term. (scalar)

J_obs = DOT_PRODUCT (dy(:), Ody(:))

! SCALAR value required

pen_func = 0.5 * (J_back + J_obs)

pen_ob = 0.5 * J_obs
pen_back = 0.5 * J_back

END SUBROUTINE Ops_GPSRO_pen

end module ufo_gnssroonedvarcheck_pen_mod
