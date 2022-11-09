!-------------------------------------------------------------------------------
! (C) British Crown Copyright 2020 Met Office
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
!-------------------------------------------------------------------------------

!> Fortran module for utility routines used in calculating the observation
!> uncertainties, using the Met Office calculations

module ufo_roobserror_utils_mod

use kinds
use fckit_log_module, only: fckit_log
use ufo_utils_mod, only: ufo_utils_iogetfreeunit
use missing_values_mod

implicit none

private
public  :: Rmatrix_type
public  :: ufo_roobserror_getrmatrix
public  :: ufo_roobserror_interpolate_rmatrix
public  :: ufo_roobserror_findnearest_rmatrix

type Rmatrix_type
  INTEGER       :: satid = 0             ! Sat ID for each R matrix
  INTEGER       :: origctr = 0           ! Originating centre for each R matrix
  REAL          :: av_temp               ! Average temperature of the troposphere (if used)
  INTEGER       :: max_height            ! The maximum height used in calculating the average temperature
  REAL          :: latitude              ! The latitude of this R matrix (if used)
  INTEGER       :: num_heights           ! Number of heights in the R matrix specification (i.e. the size of the matrix)
  REAL          :: clen                  ! Correlation length-scale for the R matrix
  REAL          :: min_error             ! Minimum value of the observation error (in radians)
  REAL, POINTER :: height(:)             ! The height for each specified observation error
  REAL, POINTER :: frac_err(:)           ! The fractional error at this height (obs error / background measurement)
end type

contains

!-------------------------------------------------------------------------------
! Read data needed to set up the Rmatrix.
!-------------------------------------------------------------------------------

SUBROUTINE ufo_roobserror_getrmatrix(max_mats, &     ! The maximum number of matrices
                                     filename, &     ! The name of the file to be read in
                                     Rmatrix, &      ! The R matrices read from the files
                                     R_num_mats)     ! Number of R matrices read in

IMPLICIT NONE

! Subroutine arguments:
INTEGER, INTENT(IN)                           :: max_mats    ! The maximum number of matrices possible
CHARACTER(LEN=*), INTENT(IN)                  :: filename    ! The name of the file to be read in
TYPE (Rmatrix_type), ALLOCATABLE, INTENT(OUT) :: Rmatrix(:)  ! The R matrices read from the files
INTEGER, INTENT(OUT)                          :: R_num_mats  ! Number of sats/processing R matrices in these files

! Local declarations:
CHARACTER(len=*), PARAMETER                   :: RoutineName = "ufo_roobserror_getrmatrix"
TYPE (Rmatrix_type)                           :: temp_rmat(max_mats)   ! Temporary store of the matrices read in
INTEGER                                       :: i                     ! Loop variable when copying matrices
INTEGER                                       :: fileunit              ! The unit number of the input file
CHARACTER(len=200)                            :: ErrorMessage          ! The error message to be output
INTEGER                                       :: return_code           ! Return code from routines called

!-----------------------------------------------
! 0. Open the R-matrix file
!-----------------------------------------------

R_num_mats = 0

fileunit = ufo_utils_iogetfreeunit()
OPEN(UNIT=fileunit, FILE=filename, ACTION='READ', STATUS='OLD', IOSTAT=return_code)
if (return_code /= 0) then
  WRITE(ErrorMessage, '(3A,I0)') "Error opening ", TRIM(filename), &
    ", return code = ", return_code
  call abor1_ftn(ErrorMessage)
end if

!-----------------------------------------------
! 1. Read all the namelists that are in the file - read these into a temporary
! array to store the data
!-----------------------------------------------

DO
  IF (R_num_mats >= max_mats) THEN
    WRITE (ErrorMessage, '(A,I0,A,I0)') &
      "Trying to read too many R-matrices ", &
      R_num_mats, " out of ", max_mats
    CALL abor1_ftn(ErrorMessage)
  END IF
  
  CALL ufo_roobserror_read_rmatrix (fileunit,                  &
                                    temp_rmat(R_num_mats + 1), &
                                    return_code)
  IF (return_code < 0) THEN
    ! End of file reached
    EXIT
  ELSE IF (return_code > 0) THEN
    WRITE(ErrorMessage, '(A,I0)') &
      "Trouble reading R-matrix namelist, error number ", &
      return_code
    CALL abor1_ftn(ErrorMessage)
  ELSE
    R_num_mats = R_num_mats + 1
  END IF
END DO

!-----------------------------------------------
! 2. Copy the R-matrices from the temporary array to the output variable
!-----------------------------------------------

ALLOCATE (RMatrix(R_num_mats))
DO i = 1, R_num_mats
  CALL ufo_roobserror_copy_rmatrix (temp_rmat(i), &
                                    RMatrix(i))
END DO

CLOSE(fileunit)

!-----------------------------------------------
! 3. Check that all the matrices have the same max_height value
!-----------------------------------------------

DO i = 2, R_num_mats
  IF (RMatrix(i) % max_height /= RMatrix(1) % max_height) THEN
    WRITE (ErrorMessage, '(A,I0,A,I0,A,I0)') &
      "All the R-matrices should have the same value for max_height Rmatrix(", &
      i, ") % max_height = ", RMatrix(i) % max_height, &
      ", RMatrix(1) % max_height = ", RMatrix(1) % max_height
    CALL abor1_ftn(ErrorMessage)
  END IF
END DO

END SUBROUTINE ufo_roobserror_getrmatrix


!-------------------------------------------------------------------------------
! (C) Crown copyright Met Office. All rights reserved.
!     Refer to COPYRIGHT.txt of this distribution for details.
!-------------------------------------------------------------------------------
! Read data needed to set up the Rmatrix.
!-------------------------------------------------------------------------------

SUBROUTINE ufo_roobserror_read_rmatrix (fileunit,    &  ! The maximum number of matrices
                                 Rmatrix,     &  ! The R matrices read from the files
                                 return_code)    ! Number of R matrices read in

IMPLICIT NONE

! Subroutine arguments:
INTEGER, INTENT(IN)              :: fileunit     ! The unit number of the open file
TYPE (Rmatrix_type), INTENT(OUT) :: Rmatrix      ! The R matrix read from the file
INTEGER, INTENT(OUT)             :: return_code  ! Return code from this routine

! Local declarations:
CHARACTER(len=*), PARAMETER      :: RoutineName = "ufo_roobserror_read_rmatrix"
CHARACTER(len=200)               :: ErrorMessage
REAL                             :: av_temp          ! Average troposphere temperature
INTEGER                          :: max_height       ! The maximum height for calculating the average temperature
REAL                             :: latitude         ! Latitude
REAL                             :: heights(600)     ! Heights of the observation errors
REAL                             :: obs_errors(600)  ! Observation errors (on vertical levels)
INTEGER                          :: satid            ! Satellite identifier
INTEGER                          :: origc            ! Originating centre for the RO data
REAL                             :: clen             ! Correlation length-scale for the R-matrix
REAL                             :: min_error        ! Minimum observation error

NAMELIST / GPSRO_ob_error / &
  av_temp,                  &
  heights,                  &
  obs_errors,               &
  satid,                    &
  origc,                    &
  clen,                     &
  min_error,                &
  latitude,                 &
  max_height

! Initialise values to blank / sensible values

satid = missing_value(satid)
origc = missing_value(origc)
min_error = 0
heights = missing_value(heights(1))
obs_errors = missing_value(obs_errors(1))
clen = 1.0E10
av_temp = missing_value(av_temp)
latitude = missing_value(latitude)
max_height = 0

! Read the namelist, and decode this into the R-matrix structure

READ (fileunit, NML = GPSRO_ob_error, IOSTAT = return_code)
IF (return_code == 0) THEN

  RMatrix % num_heights = COUNT (heights /= missing_value(heights(1)))
  IF (COUNT(obs_errors /= missing_value(obs_errors(1))) /= RMatrix % num_heights) THEN
    WRITE (ErrorMessage, '(A,I0,1X,I0)') "Counts do not match ", &
                                         COUNT(obs_errors /= missing_value(obs_errors(1))), RMatrix % num_heights
    CALL abor1_ftn(ErrorMessage)
  ELSE
    ALLOCATE (RMatrix % height(1:RMatrix % num_heights))
    ALLOCATE (RMatrix % frac_err(1:RMatrix % num_heights))
    RMatrix % height(:) = PACK (heights, heights /= missing_value(heights(1)))
    RMatrix % frac_err(:) = PACK (obs_errors, obs_errors /= missing_value(obs_errors(1)))
    RMatrix % frac_err(:) = RMatrix % frac_err(:)
  END IF

  RMatrix % av_temp = av_temp
  RMatrix % max_height = max_height
  RMatrix % latitude = latitude
  RMatrix % satid = satid
  RMatrix % origctr = origc
  RMatrix % clen = clen
  RMatrix % min_error = min_error
END IF

END SUBROUTINE ufo_roobserror_read_rmatrix


!-------------------------------------------------------------------------------
! (C) Crown copyright Met Office. All rights reserved.
!     Refer to COPYRIGHT.txt of this distribution for details.
!-------------------------------------------------------------------------------
! Choose the R matrix to apply for this observation.  The matrix will be
! be selected to have the same satellite identifier and originating centre as
! requested.  The matrices will be interpolated from those with the nearest
! average temperature value.
!-------------------------------------------------------------------------------

SUBROUTINE ufo_roobserror_interpolate_rmatrix (satid,        &
                                               origc,        &
                                               av_temp,      &
                                               R_num_sats,   &
                                               Rmatrix_list, &
                                               out_matrix)


IMPLICIT NONE

! Subroutine arguments:
INTEGER, INTENT(IN)              :: satid                     ! The satellite identifier | These three entries will
INTEGER, INTENT(IN)              :: origc                     ! The originating centre   | be used to find the
REAL(kind_real), INTENT(IN)      :: av_temp                   ! The average temperature  | R matrix from the list
INTEGER, INTENT(IN)              :: R_num_sats                ! The number of matrices in the list
TYPE (Rmatrix_type), INTENT(IN)  :: RMatrix_list(R_num_sats)  ! The list of R matrices to select from
TYPE (Rmatrix_type), INTENT(OUT) :: out_matrix                ! Output R matrix

! Local parameters
logical, parameter               :: verboseOutput = .FALSE.   ! Whether to output extra debugging information

! Local declarations:
CHARACTER(len=*), PARAMETER      :: RoutineName = "ufo_roobserror_interpolate_rmatrix"
INTEGER                          :: i                         ! Loop variable
INTEGER                          :: lower_match               ! The largest av_temp which is smaller than the one being searched for
INTEGER                          :: upper_match               ! The smallest av_temp which is larger than the one being searched for
REAL                             :: weight                    ! Weight to give to the upper matched av_temp
CHARACTER(len=200)               :: ErrorMessage              ! Error message to be output
REAL, ALLOCATABLE                :: interp_errors(:)          ! Observation errors, interpolated to the required av_temp

lower_match = 0
upper_match = 0
DO i = 1, R_num_sats
  IF (satid == RMatrix_list(i) % satid .AND. &
      origc == RMatrix_list(i) % origctr) THEN
    IF (RMatrix_list(i) % av_temp /= missing_value(RMatrix_list(i) % av_temp)) THEN
      IF (RMatrix_list(i) % av_temp > av_temp) THEN
        IF (upper_match == 0) THEN
          upper_match = i
        ELSE IF (RMatrix_list(i) % av_temp < RMatrix_list(upper_match) % av_temp) THEN
          upper_match = i
        END IF
      ELSE
        IF (lower_match == 0) THEN
          lower_match = i
        ELSE IF (RMatrix_list(i) % av_temp > RMatrix_list(lower_match) % av_temp) THEN
          lower_match = i
        END IF
      END IF
    END IF
  END IF
END DO

IF (upper_match > 0) THEN
  IF (lower_match > 0) THEN
    if (verboseOutput) then
      WRITE (ErrorMessage, '(A,3F10.2)') "Interpolating between locations.. ", &
                                         RMatrix_list(lower_match) % av_temp, &
                                         av_temp, RMatrix_list(upper_match) % av_temp
    end if
    weight = (av_temp - RMatrix_list(lower_match) % av_temp) / &
      (RMatrix_list(upper_match) % av_temp - RMatrix_list(lower_match) % av_temp)
    
    out_matrix % num_heights = RMatrix_list(lower_match) % num_heights
    ALLOCATE (out_matrix % height(1:out_matrix % num_heights))
    ALLOCATE (out_matrix % frac_err(1:out_matrix % num_heights))
    ALLOCATE (interp_errors(1:out_matrix % num_heights))
    CALL ufo_roobserror_interpolate_heights(RMatrix_list(lower_match) % num_heights, &
                                            RMatrix_list(lower_match) % height,      &
                                            RMatrix_list(upper_match) % num_heights, &
                                            RMatrix_list(upper_match) % height,      &
                                            RMatrix_list(upper_match) % frac_err,    &
                                            interp_errors)
    out_matrix % height(:) = RMatrix_list(lower_match) % height(:)
    out_matrix % frac_err(:) = weight * interp_errors + &
                                 (1 - weight) * RMatrix_list(lower_match) % frac_err

    out_matrix % av_temp = av_temp

    ! max_height should be the same for every matrix, so do not interpolate
    out_matrix % max_height = RMatrix_list(lower_match) % max_height
    out_matrix % satid = satid
    out_matrix % origctr = origc
    out_matrix % clen = weight * RMatrix_list(upper_match) % clen + &
                             (1 - weight) * RMatrix_list(lower_match) % clen
    out_matrix % min_error = weight * RMatrix_list(upper_match) % min_error + &
                             (1 - weight) * RMatrix_list(lower_match) % min_error
    out_matrix % latitude = missing_value(out_matrix % latitude)
  ELSE
    CALL ufo_roobserror_copy_rmatrix(RMatrix_list(upper_match), &
                                     out_matrix)
  END IF
ELSE
  IF (lower_match > 0) THEN
    CALL ufo_roobserror_copy_rmatrix(RMatrix_list(lower_match), &
                                     out_matrix)
  ELSE
    WRITE (ErrorMessage, '(A,I0,1X,I0,1X)') "Did not match any matrices ", satid, origc
    CALL fckit_log % warning(ErrorMessage)
    CALL ufo_roobserror_copy_rmatrix(RMatrix_list(1), &
                                     out_matrix)
  END IF
END IF

END SUBROUTINE ufo_roobserror_interpolate_rmatrix


!-------------------------------------------------------------------------------
! (C) Crown copyright Met Office. All rights reserved.
!     Refer to COPYRIGHT.txt of this distribution for details.
!-------------------------------------------------------------------------------
! Interpolate the R-matrix values between the source and target heights
!-------------------------------------------------------------------------------

SUBROUTINE ufo_roobserror_interpolate_heights(target_nheights, &
                                              target_heights,  &
                                              source_nheights, &
                                              source_heights,  &
                                              source_values,   &
                                              output_values)

IMPLICIT NONE

! Subroutine arguments:
INTEGER, INTENT(IN)         :: target_nheights
REAL, INTENT(IN)            :: target_heights(1:target_nheights)
INTEGER, INTENT(IN)         :: source_nheights
REAL, INTENT(IN)            :: source_heights(1:source_nheights)
REAL, INTENT(IN)            :: source_values(1:source_nheights)
REAL, INTENT(OUT)           :: output_values(1:target_nheights)

! Local declarations:
CHARACTER(len=*), PARAMETER :: RoutineName = "ufo_roobserror_interpolate_heights"
INTEGER                     :: i
INTEGER                     :: j
LOGICAL                     :: MatchFound
REAL                        :: weight     ! Weight to give to the upper level value

DO i = 1, target_nheights
  IF (ALL (target_heights(i) < source_heights)) THEN
    ! The target height is less than all in the source, so take the bottom value
    output_values(i) = source_values(1)
  ELSE
    ! Search for a match, where the target height is between two adjacent values
    MatchFound = .FALSE.
    DO j = 1, source_nheights - 1
      IF (target_heights(i) >= source_heights(j) .AND. &
          target_heights(i) < source_heights(j + 1)) THEN
        MatchFound = .TRUE.
        EXIT
      END IF
    END DO
    IF (MatchFound) THEN
      ! If we have found a match, then interpolate between the heights in
      ! question
      weight = (target_heights(i) - source_heights(j)) / &
               (source_heights(j + 1) - source_heights(j))
      output_values(i) = weight * source_values(j + 1) + &
                           (1 - weight) * source_values(j)
    ELSE
      ! If no match has been found, then the target is greater than all the
      ! source heights
      output_values(i) = source_values(source_nheights)
    END IF
  END IF
END DO

END SUBROUTINE ufo_roobserror_interpolate_heights


!-------------------------------------------------------------------------------
! (C) Crown copyright Met Office. All rights reserved.
!     Refer to COPYRIGHT.txt of this distribution for details.
!-------------------------------------------------------------------------------
! Read data needed to set up the Rmatrix.
!-------------------------------------------------------------------------------

SUBROUTINE ufo_roobserror_copy_rmatrix(in_mat,  &    ! The RMatrix to copy
                                       out_mat)      ! The target for the copy

IMPLICIT NONE

! Subroutine arguments:
TYPE (Rmatrix_type), INTENT(IN)  :: in_mat     ! The R matrix to copy
TYPE (Rmatrix_type), INTENT(OUT) :: out_mat    ! The target for the copy

! Local declarations:
CHARACTER(len=*), PARAMETER      :: RoutineName = "ufo_roobserror_copy_rmatrix"

out_mat % num_heights = in_mat % num_heights
ALLOCATE (out_mat % height(out_mat % num_heights))
ALLOCATE (out_mat % frac_err(out_mat % num_heights))
out_mat % height(:) = in_mat % height(:)
out_mat % frac_err(:) = in_mat % frac_err(:)

out_mat % av_temp = in_mat % av_temp
out_mat % max_height = in_mat % max_height
out_mat % latitude = in_mat % latitude
out_mat % satid = in_mat % satid
out_mat % origctr = in_mat % origctr
out_mat % clen = in_mat % clen
out_mat % min_error = in_mat % min_error

END SUBROUTINE ufo_roobserror_copy_rmatrix

!-------------------------------------------------------------------------------
! Choose the R matrix to apply for this observation.  The matrix will be
! be selected to have the same satellite identifier and originating centre as
! requested.  The matrix with the nearest latitude to that specified will be
! selected (without interpolation).
!-------------------------------------------------------------------------------

SUBROUTINE ufo_roobserror_findnearest_rmatrix(satid,        &
                                              origc,        &
                                              latitude,     &
                                              R_num_sats,   &
                                              Rmatrix_list, &
                                              out_matrix)

IMPLICIT NONE

! Subroutine arguments:
INTEGER, INTENT(IN)              :: satid                       ! The satellite identifier | These three entries will
INTEGER, INTENT(IN)              :: origc                       ! The originating centre   | be used to find the
REAL(kind_real), INTENT(IN)      :: latitude                    ! The latitude             | R matrix from the list
INTEGER, INTENT(IN)              :: R_num_sats                  ! The number of matrices in the list
TYPE (Rmatrix_type), INTENT(IN)  :: RMatrix_list(1:R_num_sats)  ! The list of R matrices to select from
TYPE (Rmatrix_type), INTENT(OUT) :: out_matrix                  ! Output R matrix

! Local declarations:
CHARACTER(len=*), PARAMETER      :: RoutineName = "ufo_roobserror_findnearest_rmatrix"
INTEGER                          :: i                           ! Loop variable
INTEGER                          :: closest_match               ! The R matrix whose latitude is closest to the specification
CHARACTER(len=200)               :: ErrorMessage                ! Error message to be output

closest_match = -1
DO i = 1, R_num_sats
  IF (satid == RMatrix_list(i) % satid .AND. &
      origc == RMatrix_list(i) % origctr) THEN
    IF (RMatrix_list(i) % latitude /= missing_value(RMatrix_list(i) % latitude)) THEN
      IF (closest_match < 1) THEN
        closest_match = i
      ELSE IF (ABS (latitude - RMatrix_list(i) % latitude) < &
               ABS (latitude - RMatrix_list(closest_match) % latitude)) THEN
        closest_match = i
      END IF
    END IF
  END IF
END DO

IF (closest_match > 0) THEN
  CALL ufo_roobserror_copy_rmatrix(RMatrix_list(closest_match), &
                                   out_matrix)
ELSE
  WRITE (ErrorMessage, '(A,I0,1X,I0)') "Did not match any matrices ", satid, origc
  CALL fckit_log % warning(ErrorMessage)
  CALL ufo_roobserror_copy_rmatrix(RMatrix_list(1), &
                                   out_matrix)
END IF

END SUBROUTINE ufo_roobserror_findnearest_rmatrix

end module ufo_roobserror_utils_mod
