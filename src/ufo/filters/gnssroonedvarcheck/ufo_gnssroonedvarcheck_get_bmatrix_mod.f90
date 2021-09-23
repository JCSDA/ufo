!-------------------------------------------------------------------------------
! (C) British Crown Copyright 2020 Met Office
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
!-------------------------------------------------------------------------------

!> Set up the background error covariance matrix (B matrix).

module ufo_gnssroonedvarcheck_get_bmatrix_mod

use kinds, only: kind_real
use fckit_log_module, only : fckit_log
use ufo_utils_mod, only: ufo_utils_iogetfreeunit

private
public :: Bmatrix_type

type Bmatrix_type
  INTEGER       :: nlevp
  INTEGER       :: nlevq
  INTEGER       :: nstate
  INTEGER       :: nband
  INTEGER       :: nseason
  REAL(kind_real), POINTER :: band_up_lim(:)      ! band_up_lim(nband)
  REAL(kind_real), POINTER :: sigma(:,:,:)        ! sigma(nseason,nband,nstate)
  REAL(kind_real), POINTER :: inverse(:,:,:,:)    ! inverse(nseason,nband,nstate,nstate)
  CONTAINS
    procedure :: get => Ops_GPSRO_GetBmatrix
end type

contains

SUBROUTINE Ops_GPSRO_GetBmatrix (Bmatrix, &
                                 filename, &
                                 cx_nlevp, &
                                 cx_nlevq)

IMPLICIT NONE

! Subroutine arguments:
CLASS(Bmatrix_type), INTENT(OUT) :: Bmatrix    !< The background errors read in
CHARACTER(LEN=*)                 :: filename   !< The name of the file to be read in
INTEGER, INTENT(IN)              :: cx_nlevp   !< The number of pressure levels in the model
INTEGER, INTENT(IN)              :: cx_nlevq   !< The number of temperature levels in the model

! Local declarations:
CHARACTER(len=*), PARAMETER      :: RoutineName = "Ops_GPSRO_GetBmatrix"
INTEGER                          :: i
INTEGER                          :: j
INTEGER                          :: m
INTEGER                          :: n
INTEGER                          :: nlevp
INTEGER                          :: nlevq
INTEGER                          :: nstate
INTEGER                          :: nband
INTEGER                          :: nseason
INTEGER                          :: fileunit
CHARACTER(len=*), PARAMETER      :: filetype_name = "Bmatrix"
CHARACTER(len=20)                :: prefix
CHARACTER(len=256)               :: ErrorMessage
INTEGER                          :: return_code

!-----------------------------------------------
! 0. Determine BMatrix environment variable name
!-----------------------------------------------

prefix = 'GPSRO_'
fileunit = ufo_utils_iogetfreeunit()

OPEN(UNIT=fileunit, FILE=filename, ACTION='READ', STATUS='OLD', IOSTAT=return_code)
if (return_code /= 0) then
  WRITE(ErrorMessage, '(3A,I0)') "Error opening ", TRIM(filename), &
    ", return code = ", return_code
  call abor1_ftn(ErrorMessage)
end if

!---------------------
! 2. Read in the file
!---------------------

READ (fileunit, '(5I5)') nlevp, nlevq, nstate, nband, nseason

IF (cx_nlevp /= nlevp) THEN

  WRITE (ErrorMessage, '(A,I0,A,I0)')'nlevp = ', nlevp, ' cx_nlevp = ', cx_nlevp
  call fckit_log % error(ErrorMessage)
  ErrorMessage = 'no. of pressure levels in vector and bmatrix not the same'
  call abor1_ftn(ErrorMessage)

END IF

IF (cx_nlevq /= nlevq) THEN

  WRITE (ErrorMessage, '(A,I0,A,I0)') 'nlevq = ', nlevq, ' cx_nlevq = ', cx_nlevq
  call fckit_log % error(ErrorMessage)
  ErrorMessage = 'no. of humidity levels in vector and bmatrix not the same'
  call abor1_ftn(ErrorMessage)

END IF

! Allocate storage variables

Bmatrix % nlevp = nlevp
Bmatrix % nlevq = nlevq
Bmatrix % nstate = nstate
Bmatrix % nband = nband
Bmatrix % nseason = nseason

! Allocate the arrays in Bmatrix type

ALLOCATE (Bmatrix % band_up_lim(nband))
ALLOCATE (Bmatrix % sigma(nseason,nband,nstate))
ALLOCATE (Bmatrix % inverse(nseason,nband,nstate,nstate))

! Read the band upper limit

READ (fileunit, '(3F5.1)') (Bmatrix % band_up_lim(i), i = 1, nband)

DO n = 1,nseason

  DO m = 1,nband

    READ (fileunit, *)  ! space

    ! Read in the sigma values

    READ (fileunit, '(10E15.6)') (Bmatrix % sigma (n,m,i), i = 1, nstate)

    ! Read in the inverse B matrix

    DO i = 1,nstate

      READ (fileunit, *)  ! space
      READ (fileunit, '(10E15.6)') (Bmatrix % inverse (n,m,i,j), j = 1, nstate)

    END DO ! each B matrix

  END DO ! each band

END DO  ! each season

! Close the file

CLOSE(fileunit)

END SUBROUTINE Ops_GPSRO_GetBmatrix

end module ufo_gnssroonedvarcheck_get_bmatrix_mod
