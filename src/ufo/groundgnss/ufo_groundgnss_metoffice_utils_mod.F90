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
    rd,                     &    ! Gas constant for dry air
    grav,                   &    ! Gravitational field strength
    n_alpha,                &    ! Refractivity constant a
    n_beta                       ! Refractivity constant b

implicit none
public                      :: Ops_Groundgnss_ZTD
public                      :: Ops_groundgnss_TopCorrection
REAL(kind_real), PARAMETER  :: refrac_scale = 1.0E-6
private

contains

!-------------------------------------------------------------------------------
! Ground based GNSS Observation operator
!-------------------------------------------------------------------------------
SUBROUTINE Ops_Groundgnss_ZTD  (nlevq,  &
                                refrac, &
                                zb,     &
                                zStation, &
                                Model_ZTD)

IMPLICIT NONE
    
    INTEGER, INTENT(IN)               :: nlevq           ! no. of theta levels
    REAL(kind_real), INTENT(IN)       :: refrac(:)
    REAL(kind_real), INTENT(IN)       :: zb(:)
    REAL(kind_real), INTENT(IN)       :: zStation
    REAL(kind_real), INTENT(INOUT)    :: Model_ZTD
    
    
    REAL(kind_real)             :: LocalZenithDelay
    INTEGER                     :: Level
    REAL(kind_real)             :: StationRefrac
    REAL(kind_real)             :: c
    REAL(kind_real)             :: const
    REAL(kind_real)             :: term1
    REAL(kind_real)             :: term2
    INTEGER                     :: Lowest_Level

    !------------------------------------------------------------
    ! Calculate the zenith delay for each layer and add to total
    !------------------------------------------------------------

    StationRefrac = 0.0

    DO Level = 1, nlevq
      IF (zb(Level) > zStation) THEN
        Lowest_Level = Level
        EXIT
      END IF
    END DO

    ! Start at bottom level

    DO Level = Lowest_Level, nlevq

      LocalZenithDelay = 0.0

      IF (Level == Lowest_Level .AND. Level /= 1) THEN

        ! If station lies above the lowest model level, interpolate refractivity
        ! to station height

        c = (LOG (refrac(Level) / refrac(Level - 1))) / (zb(Level - 1) - zb(Level))
        StationRefrac = refrac(Level - 1) * EXP (-c * (zStation - zb(Level - 1)))
        const = -StationRefrac / c * EXP (c * zStation)
        term1 = EXP (-c * (zb(Level)))
        term2 = EXP (-c * zStation)
        LocalZenithDelay = refrac_scale * const * (term1 - term2)

      ELSE IF (Level == 1) THEN

        ! If station lies below model level 1 (ie. the lowest level for which refractivity is
        ! calculated, then use c from the first full layer, but integrate down to height of
        ! station

        c = (LOG (Refrac(Level + 1) / Refrac(Level))) / (zb(Level) - zb(Level + 1))
        const = -refrac(Level) / c * EXP (c * (zb(Level)))
        term1 = EXP (-c * (zb(Level + 1)))
        term2 = EXP (-c * ZStation)
        LocalZenithDelay = refrac_scale * const * (term1 - term2)

      ELSE IF (Level <= nlevq .AND. Level > 2) THEN

        ! If not at top level

        c = (LOG (refrac(Level) / refrac(Level - 1))) / (zb(Level - 1) - zb(Level))
        const = -refrac(Level - 1) / c * EXP (c * (zb(Level - 1)))
        term1 = EXP (-c * (zb(Level)))
        term2 = EXP (-c * (zb(Level - 1)))
        LocalZenithDelay = refrac_scale * const * (term1 - term2)

      END IF

      Model_ZTD = Model_ZTD + LocalZenithDelay
  END DO
END SUBROUTINE Ops_Groundgnss_ZTD


SUBROUTINE Ops_groundgnss_TopCorrection(pN,    &
                                        nlevq, &
                                        TopCorrection)
					
    IMPLICIT NONE
    
    REAL(kind_real), INTENT(IN)      :: pN(:)
    INTEGER, INTENT(IN)              :: nlevq
    REAL(kind_real), INTENT(INOUT)   :: TopCorrection

    REAL(kind_real)                  :: TCconstant
    REAL(kind_real), PARAMETER       :: hpa_to_pa = 100.0

    TCconstant = (refrac_scale * n_alpha * rd)/ (hpa_to_pa * grav)
    
    TopCorrection = TCconstant * pN(nlevq)
    
END SUBROUTINE Ops_groundgnss_TopCorrection

END MODULE ufo_groundgnss_ukmo_utils_mod
