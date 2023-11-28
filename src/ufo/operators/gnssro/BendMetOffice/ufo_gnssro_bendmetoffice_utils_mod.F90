!-------------------------------------------------------------------------------
! (C) British Crown Copyright 2020 Met Office
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
!-------------------------------------------------------------------------------
module ufo_gnssro_ukmo1d_utils_mod

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
    pi,                     &    ! Something to do with circles...
    grav,                   &    ! Gravitational field strength
    ecc,                    &    ! eccentricity
    k_somig,                &    ! Somigliana's constant
    g_equat,                &    ! equatorial gravity (ms-2)
    a_earth,                &    ! semi-major axis of earth (m)
    flatt,                  &    ! flattening
    m_ratio,                &    ! gravity ratio
    mw_ratio,               &    ! Ratio of molecular weights of water and dry air
    c_virtual,              &    ! Related to mw_ratio
    n_alpha,                &    ! Refractivity constant a
    n_beta                       ! Refractivity constant b

implicit none
public             :: Ops_GPSROcalc_alpha
public             :: Ops_GPSROcalc_nr
public             :: Ops_GPSRO_geop_geom
private

contains

!-------------------------------------------------------------------------------
! GPSRO 1D bending angle operator.
!-------------------------------------------------------------------------------
SUBROUTINE Ops_GPSROcalc_alpha (nobs,   &
                                nlev,   &
                                a,      &
                                refrac, &
                                nr,     &
                                alpha,  &
                                noCheck)

IMPLICIT NONE

! Subroutine arguments:
INTEGER, INTENT(IN)          :: nobs           ! size of ob. vector
INTEGER, INTENT(IN)          :: nlev           ! no. of refractivity levels
REAL(kind_real), INTENT(IN)  :: a(nobs)        ! observation impact parameters
REAL(kind_real), INTENT(IN)  :: refrac(nlev)   ! refractivity values on model levels
REAL(kind_real), INTENT(IN)  :: nr(nlev)       ! refractive index * radius product
REAL(kind_real), INTENT(OUT) :: alpha(nobs)    ! bending angle
LOGICAL, INTENT(IN)          :: noCheck        ! If true, do not apply super-refraction check

! Local declarations:
CHARACTER(len=*), PARAMETER :: RoutineName = "Ops_GPSROcalc_alpha"
INTEGER                     :: i
INTEGER                     :: n
INTEGER                     :: jbot
INTEGER                     :: kbot
REAL(kind_real)             :: kval(nlev - 1) ! exponential decay rate between levels
REAL(kind_real)             :: root_2pia
REAL(kind_real)             :: ref_low        ! refractivity at lower level
REAL(kind_real)             :: nr_low         ! impact parameter lower level
REAL(kind_real)             :: tup
REAL(kind_real)             :: tlow           ! upper/lower bounds to error function
REAL(kind_real)             :: diff_erf       ! error function result
REAL(kind_real)             :: dalpha         ! delta bending angle
REAL(kind_real)             :: erf_tup
REAL(kind_real)             :: erf_tlow       ! error function of tup and tlow
REAL(kind_real)             :: t              ! int. step in approx. to error function
REAL(kind_real), PARAMETER  :: a1 = 0.3480242 ! consts used in error function approx.
REAL(kind_real), PARAMETER  :: a2 = -0.0958798
REAL(kind_real), PARAMETER  :: a3 = 0.7478556
REAL(kind_real), PARAMETER  :: p = 0.47047
INTEGER                     :: indices(nlev)  ! The indices to use in the refractivity profile
INTEGER                     :: ngood          ! The number of refractivity values to use
REAL(kind_real)             :: minNR          ! Minimum value of nr found so far
INTEGER                     :: iMin           ! Minimum value of indices to use in profile

jbot = 1

DO

  IF (refrac(jbot) > 0.0 .AND. nr(jbot) > 0.0) EXIT

  jbot = jbot + 1

END DO

!-------------------------------------------------------------------------------
! Calculate the indices of the model levels which can be used in calculating the
! bending angles.
! If no super-refraction check, then search downwards and create a profile where
! the impact parameter is monotonically decreasing.
! If using the super-refraction check, then search downwards; if the impact
! parameter decreases by less than 10 metres, then reject all model levels below
! this point.
!-------------------------------------------------------------------------------

if (noCheck) then
  ! Remove regions where the IP reduces in the model
  nGood = 0
  indices = 0
  minNR = nr(nlev)
  do i = nlev, jbot, -1
    if (nr(i) <= minNR) then
      minNR = nr(i)
      nGood = nGood + 1
      indices(nGood) = i
    end if
  end do
  indices(1:nGood) = indices(nGood:1:-1)
else
  ! Remove regions below points where the impact parameter decreases by less
  ! than 10 metres.
  kbot = nlev
  DO i = nlev, jbot + 1, -1
    ! to avoid large gradients
    IF ((nr(kbot) - nr(kbot - 1)) < 10.0) EXIT
    kbot = kbot - 1
  END DO

  jbot = MAX (jbot,kbot)
  nGood = 0
  do i = jbot, nlev
    nGood = nGood + 1
    indices(nGood) = i
  end do
end if

!-------------------------------------------------------------------------------
! Calculate the exponential decay rate between levels
!-------------------------------------------------------------------------------

DO i = 1, nGood - 1

  kval(i) = LOG(refrac(indices(i)) / refrac(indices(i+1))) / &
               MAX(1.0,(nr(indices(i+1)) - nr(indices(i))))
  kval(i) = MAX(1.0E-6, kval(i))

END DO

!-------------------------------------------------------------------------------
! Calculate the bending angle values
!-------------------------------------------------------------------------------

alpha(:) = missing_value(alpha(1))

DO n = 1, nobs

  IF (a(n) < nr(indices(1)) .OR. a(n) > nr(indices(nGood))) CYCLE

  Root_2PIa = SQRT (2.0 * pi * a(n))

  ! Find bottom state vector level
  !----------------------------------

  iMin = 1
  DO
    ! check more than 1 metre apart to stop large gradients in K code
    IF (((nr(indices(iMin + 1)) - a(n)) > 1.0) .OR. nGood < iMin + 2) EXIT
    iMin = iMin + 1
  END DO

  ! Initialise bending angle value
  !-------------------------------

  alpha(n) = 0.0

  ! Values of refractivity and impact parameter at lower level
  !-----------------------------------------------------------
!The following line is a compiler directive to avoid bugs with intel compilers
!DIR$ NOVECTOR
  DO i = iMin, nGood - 1

    IF (i == iMin) THEN
      ref_low = refrac(indices(i)) * EXP (-kval(i) * (a(n) - nr(indices(i))))
      nr_low = a(n)
    ELSE
      ref_low = refrac(indices(i))
      nr_low = nr(indices(i))
    END IF

    ! Limits used in the error function
    !----------------------------------

    IF (i == nGood - 1) THEN
      ! Simple extrapolation 100km above the uppermost level
      !-----------------------------------------------------
      tup = SQRT(kval(i) * (nr(indices(i+1)) + 1.0E5 - a(n)))
    ELSE
      tup = SQRT(kval(i) * (nr(indices(i+1)) - a(n)))
    END IF

    tlow = 0.0
    IF (i > iMin) tlow = SQRT(kval(i) * (nr(indices(i)) - a(n)))

    ! Abramowitz and Stegun approx. to error function
    !------------------------------------------------
    t = 1.0 / (1.0 + p * tup)
    erf_tup = 1.0 - EXP (-(tup ** 2)) * (a1 + (a2 + a3 * t) * t) * t

    t = 1.0 / (1.0 + p * tlow)
    erf_tlow = 1.0 - EXP (-(tlow ** 2)) * (a1 + (a2 + a3 * t) * t) * t

    diff_erf = erf_tup - erf_tlow

    dalpha =  1.0E-6 * Root_2PIa * SQRT (kval(i)) * &
                  ref_low * EXP (kval(i) * (nr_low - a(n))) * diff_erf

    ! Bending angle value
    !--------------------
    alpha(n) = alpha(n) + dalpha

  END DO

END DO

END SUBROUTINE Ops_GPSROcalc_alpha


!-------------------------------------------------------------------------------
! (C) Crown copyright Met Office. All rights reserved.
!     Refer to COPYRIGHT.txt of this distribution for details.
!-------------------------------------------------------------------------------
SUBROUTINE Ops_GPSROcalc_nr (nb,     & ! number of levels in zb
                             zb,     & ! geopotential heights of model levels
                             refrac, & ! refractivity of model on model levels
                             Rad,    & ! radius of curvature of earth at observation
                             lat,    & ! latitude at observation
                             und,    & ! geoid undulation above WGS-84
                             nr)       ! Calculated model impact parameters

IMPLICIT NONE

! Subroutine arguments:
INTEGER, INTENT(IN)          :: nb          ! number of levels in zb
REAL(kind_real), INTENT(IN)  :: zb(nb)      ! geopotential heights on zb levels /m
REAL(kind_real), INTENT(IN)  :: Rad         ! local radius of curvature of earth /m
REAL(kind_real), INTENT(IN)  :: lat         ! latitude at observation/ degrees
REAL(kind_real), INTENT(IN)  :: und         ! geoid undulation
REAL(kind_real), INTENT(IN)  :: refrac(nb)  ! refractivity on model levels / N
REAL(kind_real), INTENT(OUT) :: nr(nb)      ! Calculated model impact parameters

! Local declarations:
CHARACTER(len=*), PARAMETER  :: RoutineName = "Ops_GPSROcalc_nr"
REAL(kind_real)              :: r(nb)       ! radius of model levels /m
REAL(kind_real)              :: z(nb)       ! geopotential heights on zb levels /m, local copy of zb
INTEGER                      :: i

nr(:) = missing_value(nr(1))

!----------------------------------------------
! 1. Convert zb values to geometric altitudes
!---------------------------------------------
z = zb + und               ! approx. convert to geopotential above WGS-84 ellipsoid
CALL Ops_GPSRO_geop_geom (lat, &
                          z)

!--------------------------------------------------
! 2. Calculate x =nr, i.e. model impact parameters
!--------------------------------------------------

r = Rad + z

DO i = 1,nb
  IF (zb(i) > 0.0 .AND. refrac(i) > 0.0) THEN
    nr(i) = (1.0E0 + 1.0E-6 * refrac(i)) * r(i)
  END IF
END DO

END SUBROUTINE Ops_GPSROcalc_nr


!-------------------------------------------------------------------------------
! Convert geopotential height above ellispoid to geometric height on WGS-84
! (ellipsoid).  The method for this calculation is available from
! http://mtp.jpl.nasa.gov/notes/altitude/altitude.html.
!-------------------------------------------------------------------------------
SUBROUTINE Ops_GPSRO_geop_geom (lat, &
                                z)

IMPLICIT NONE

! Subroutine arguments:
REAL(kind_real), INTENT(IN)    :: lat         ! latitude of observation
REAL(kind_real), INTENT(INOUT) :: z(:)        ! geopotential height in, geometric height out

! Local declarations:
CHARACTER(len=*), PARAMETER :: RoutineName = "Ops_GPSRO_geop_geom"
REAL(kind_real)             :: r_eff   ! effective radius of Earth, formulas from MJ Mahoney (2005)
REAL(kind_real)             :: g_somig ! Somigliana's equation for normal gravity on the surface of an ellipsoid of revolution
REAL(kind_real)             :: latrad  ! latitude in radians

!------------------------
! 1. correct units
!-----------------------

latrad = lat * (Pi / 180.0)             !convert latitude from degrees to radians

!-------------------------
! 2. Calculate r and g
!-------------------------

r_eff = a_earth / (1.0 + flatt + m_ratio - 2.0 * flatt * (SIN (latrad)) ** 2)

g_somig = g_equat * (1.0 + k_somig * (SIN (latrad)) ** 2) / (SQRT (1.0 - (ecc ** 2) * (SIN (latrad)) ** 2))

!-------------------------------------------------------------------
! 3. convert z (in geopotential height) to geometric wrt ellipsoid
!-------------------------------------------------------------------

z(:) = (r_eff * z(:)) / ((g_somig / grav) * r_eff - z(:))

END SUBROUTINE Ops_GPSRO_geop_geom

!-------------------------------------------------------------------------
end module ufo_gnssro_ukmo1d_utils_mod
