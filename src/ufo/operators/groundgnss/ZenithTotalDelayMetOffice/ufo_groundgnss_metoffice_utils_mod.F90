!-------------------------------------------------------------------------------
! (C) British Crown Copyright 2020 Met Office
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
!-------------------------------------------------------------------------------

module ufo_groundgnss_ukmo_utils_mod

!use iso_c_binding
use fckit_log_module, only: fckit_log
use kinds,            only: kind_real

! Generic routines from elsewhere in jedi
use missing_values_mod
use ufo_constants_mod, only: &
    rd,                      &   ! Gas constant for dry air
    grav,                    &   ! Gravitational field strength
    n_alpha                      ! Refractivity constant a

implicit none
public                      :: Ops_Groundgnss_ZTD
public                      :: Ops_groundgnss_TopCorrection
REAL(kind_real), PARAMETER  :: refrac_scale = 1.0E-6
private

contains

!-------------------------------------------------------------------------------
! Ground based GNSS Observation operator
!-------------------------------------------------------------------------------
SUBROUTINE Ops_Groundgnss_ZTD  (nlevq,    &
                                refrac,   &
                                zb,       &
                                zStation, &
                                Model_ZTD)

IMPLICIT NONE

    INTEGER, INTENT(IN)               :: nlevq           ! no. of temperature/theta levels
    REAL(kind_real), INTENT(IN)       :: refrac(:)       ! refractivty
    REAL(kind_real), INTENT(IN)       :: zb(:)           ! The geometric height of the model theta levels
    REAL(kind_real), INTENT(IN)       :: zStation        ! Station heights
    REAL(kind_real), INTENT(INOUT)    :: Model_ZTD       ! Model background ZTD


    REAL(kind_real)             :: LocalZenithDelay      ! Zenith total delay
    INTEGER                     :: Level                 ! level iterator
    REAL(kind_real)             :: StationRefrac         ! refrac at station
    REAL(kind_real)             :: c                     ! scale height
    REAL(kind_real)             :: const                 ! refrac constant
    REAL(kind_real)             :: term1                 ! refrac term1
    REAL(kind_real)             :: term2                 ! refrac term2
    INTEGER                     :: Lowest_Level          ! lowest height level

    !------------------------------------------------------------
    ! Calculate the zenith delay for each layer and add to total
    !------------------------------------------------------------

    StationRefrac = 0.0

    DO Level = nlevq, 1, -1
      IF (zb(Level) > zStation) THEN
        Lowest_Level = Level
        EXIT
      END IF
    END DO

    ! Start at bottom level
    ! The routine works from the bottom to the top
    ! Geovals are order top to bottom.

    DO Level = Lowest_Level, 1, -1

      LocalZenithDelay = 0.0

      IF (Level == Lowest_Level .AND. Level /= nlevq) THEN

        ! If station lies above the lowest model level, interpolate refractivity
        ! to station height

        c = (LOG (refrac(Level) / refrac(Level + 1))) / (zb(Level + 1) - zb(Level))
        StationRefrac = refrac(Level + 1) * EXP (-c * (zStation - zb(Level + 1)))
        const = -StationRefrac / c * EXP (c * zStation)
        term1 = EXP (-c * (zb(Level)))
        term2 = EXP (-c * zStation)
        LocalZenithDelay = refrac_scale * const * (term1 - term2)

      ELSE IF (Level == nlevq) THEN

        ! If station lies below model level 1 (ie. the lowest level for which refractivity is
        ! calculated, then use c from the first full layer, but integrate down to height of
        ! station

        c = (LOG (refrac(Level - 1) / refrac(Level))) / (zb(Level) - zb(Level - 1))
        const = -refrac(Level) / c * EXP (c * (zb(Level)))
        term1 = EXP (-c * (zb(Level - 1)))
        term2 = EXP (-c * ZStation)
        LocalZenithDelay = refrac_scale * const * (term1 - term2)

      ELSE IF (Level >= 1 .AND. Level < nlevq-1) THEN

        ! If not at top level

        c = (LOG (refrac(Level) / refrac(Level + 1))) / (zb(Level + 1) - zb(Level))
        const = -refrac(Level + 1) / c * EXP (c * (zb(Level + 1)))
        term1 = EXP (-c * (zb(Level)))
        term2 = EXP (-c * (zb(Level + 1)))
        LocalZenithDelay = refrac_scale * const * (term1 - term2)

      END IF

      Model_ZTD = Model_ZTD + LocalZenithDelay
  END DO
END SUBROUTINE Ops_Groundgnss_ZTD


SUBROUTINE Ops_groundgnss_TopCorrection(P,    &
                                        nlevq, &
                                        za,    &
                                        zb,    &
                                        TopCorrection)

    IMPLICIT NONE

    REAL(kind_real), INTENT(IN)      :: P(:)               ! Pressure on pressure (rho) levels
    INTEGER, INTENT(IN)              :: nlevq              ! no. of temperature/theta levels
    REAL(kind_real), INTENT(IN)      :: za(:)              ! heights of pressure (rho) levels
    REAL(kind_real), INTENT(IN)      :: zb(:)              ! Heights of temperature/theta levels
    REAL(kind_real), INTENT(INOUT)   :: TopCorrection      ! ZTD Top of atmos correction


    INTEGER                          :: Level              ! Loop counter
    REAL(kind_real)                  :: pN(nlevq)          ! Pressure on theta levels
    REAL(kind_real)                  :: pwt1               ! Weighting variable
    REAL(kind_real)                  :: pwt2               ! Weighting variable
    REAL(kind_real)                  :: TCconstant         ! Top correction constant

    DO Level = 1, nlevq

        pwt1 = (za(Level) - zb(Level)) / (za(Level) - za(Level+1))
        pwt2 = 1.0 - pwt1

        pN(Level) = EXP (pwt1 * LOG (P(Level + 1)) + pwt2 * LOG (P(Level)))

    END DO

    TCconstant = (refrac_scale * n_alpha * rd)/ grav

    TopCorrection = TCconstant * pN(1) !pN at the model top

END SUBROUTINE Ops_groundgnss_TopCorrection

END MODULE ufo_groundgnss_ukmo_utils_mod
