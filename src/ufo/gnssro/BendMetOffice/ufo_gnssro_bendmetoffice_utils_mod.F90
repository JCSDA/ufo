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
public             :: Ops_GPSRO_refrac
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
                                alpha)

IMPLICIT NONE

! Subroutine arguments:
INTEGER, INTENT(IN)          :: nobs           ! size of ob. vector
INTEGER, INTENT(IN)          :: nlev           ! no. of refractivity levels
REAL(kind_real), INTENT(IN)  :: a(nobs)        ! observation impact parameters
REAL(kind_real), INTENT(IN)  :: refrac(nlev)   ! refractivity values on model levels
REAL(kind_real), INTENT(IN)  :: nr(nlev)       ! refractive index * radius product
REAL(kind_real), INTENT(OUT) :: alpha(nobs)    ! bending angle

! Local declarations:
CHARACTER(len=*), PARAMETER :: RoutineName = "Ops_GPSROcalc_alpha"
INTEGER                     :: i
INTEGER                     :: n
INTEGER                     :: ibot
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

jbot = 1

DO

  IF (refrac(jbot) > 0.0 .AND. nr(jbot) > 0.0) EXIT

  jbot = jbot + 1

END DO

!-------------------------------------------------------------------------------
! Calculate lowest usable level (because of superrefraction)
!-------------------------------------------------------------------------------

kbot = nlev

DO i = nlev,jbot + 1,-1

  ! to avoid large gradients
  IF ((nr(kbot) - nr(kbot - 1)) < 10.0) EXIT

  kbot = kbot - 1

END DO

jbot = MAX (jbot,kbot)

!-------------------------------------------------------------------------------
! Calculate the exponential decay rate between levels
!-------------------------------------------------------------------------------

DO i = jbot,nlev - 1

  kval(i) = LOG (refrac(i) / refrac(i + 1)) / &
               MAX (1.0,(nr(i + 1) - nr(i)))
  kval(i) = MAX (1.0E-6,kval(i))

END DO

!-------------------------------------------------------------------------------
! Calculate the bending angle values
!-------------------------------------------------------------------------------

alpha(:) = missing_value(alpha(1))

DO n = 1,nobs

  IF (a(n) < nr(jbot) .OR. a(n) > nr(nlev)) CYCLE

  Root_2PIa = SQRT (2.0 * pi * a(n))

  ibot = jbot

  ! Find bottom state vector level
  !----------------------------------

  DO
    ! check more than 1 metre apart to stop large gradients in K code
    IF (((nr(ibot + 1) - a(n)) > 1.0) .OR. ibot == nlev - 1) EXIT

    ibot = ibot + 1

  END DO

  ! Initialise bending angle value
  !-------------------------------

  alpha(n) = 0.0

  ! Values of refractivity and impact parameter at lower level
  !-----------------------------------------------------------
!The following line is a compiler directive to avoid bugs with intel compilers
!DIR$ NOVECTOR
  DO i = ibot, nlev - 1

    IF (i == ibot) THEN

      ref_low = refrac(ibot) * EXP (-kval(ibot) * (a(n) - nr(ibot)))
      nr_low = a(n)

    ELSE

      ref_low = refrac(i)
      nr_low = nr(i)

    END IF

    ! Limits used in the error function
    !----------------------------------

    IF (i == nlev - 1) THEN

      ! Simple extrapolation 100km above the uppermost level
      !-----------------------------------------------------
      tup = SQRT (kval(i) * (nr(i + 1) + 1.0E5 - a(n)))

    ELSE

      tup = SQRT (kval(i) * (nr(i + 1) - a(n)))

    END IF

    tlow = 0.0

    IF (i > ibot) tlow = SQRT (kval(i) * (nr(i) - a(n)))

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
SUBROUTINE Ops_GPSROcalc_nr (zb,     & ! geopotential heights of model levels
                             nb,     & ! number of levels in zb
                             Rad,    & ! radius of curvature of earth at observation
                             lat,    & ! latitude at observation
                             und,    & ! geoid undulation above WGS-84
                             refrac, & ! refractivity of model on model levels
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


!-------------------------------------------------------------------------------
! (C) Crown copyright Met Office. All rights reserved.
!     Refer to COPYRIGHT.txt of this distribution for details.
!-------------------------------------------------------------------------------
! GPSRO refractivity forward operator
!-------------------------------------------------------------------------------

SUBROUTINE Ops_GPSRO_refrac (nlevp,     &
                             nlevq,     &
                             za,        &
                             zb,        &
                             x,         &
                             GPSRO_pseudo_ops, &
                             GPSRO_vert_interp_ops, &
                             refracerr, &
                             refrac,    &
                             T,         &
                             z_pseudo,  &
                             N_pseudo,  &
                             nb_pseudo)

IMPLICIT NONE

! Subroutine arguments:
INTEGER, INTENT(IN)                      :: nlevp                  ! no. of p levels in state vec.
INTEGER, INTENT(IN)                      :: nlevq                  ! no. of theta levels
REAL(kind_real), INTENT(IN)              :: za(:)                  ! heights of rho levs
REAL(kind_real), INTENT(IN)              :: zb(:)                  ! heights of theta levs
REAL(kind_real), INTENT(IN)              :: x(:)                   ! state vector
LOGICAL, INTENT(IN)                      :: GPSRO_pseudo_ops       ! Use pseudo levels?
LOGICAL, INTENT(IN)                      :: GPSRO_vert_interp_ops  ! Interpolate using ln(p)
LOGICAL, INTENT(OUT)                     :: refracerr              ! refractivity error
REAL(kind_real), INTENT(OUT)             :: refrac(nlevq)          ! refrac on theta levs
REAL(kind_real), INTENT(OUT)             :: T(nlevq)               ! Temp. on theta levs
REAL(kind_real), ALLOCATABLE, INTENT(OUT), OPTIONAL :: z_pseudo(:) ! height of pseudo levs
REAL(kind_real), ALLOCATABLE, INTENT(OUT), OPTIONAL :: N_pseudo(:) ! Ref. on pseudo levs
INTEGER, INTENT(OUT), OPTIONAL           :: nb_pseudo              ! no. of pseudo levs

! Local declarations:
integer, parameter           :: max_string = 800
CHARACTER(len=*), PARAMETER  :: RoutineName = "Ops_GPSRO_refrac"
CHARACTER(len=max_string)    :: message
INTEGER                      :: i
INTEGER                      :: counter
INTEGER                      :: search_lev    ! The vertical level to start searching from to
                                              ! find matching temperature-level heights
INTEGER                      :: this_lev      ! Matching level to temperature-level height
REAL(kind_real)              :: P(nlevp)
REAL(kind_real)              :: Exner(nlevp)
REAL(kind_real)              :: q(nlevq)
REAL(kind_real)              :: Pb(nlevq)
REAL(kind_real)              :: Tv
REAL(kind_real)              :: Ex_theta
REAL(kind_real)              :: pwt1
REAL(kind_real)              :: pwt2
REAL(kind_real)              :: Ndry
REAL(kind_real)              :: Nwet
REAL(kind_real), ALLOCATABLE :: P_pseudo(:)
REAL(kind_real), ALLOCATABLE :: q_pseudo(:)
REAL(kind_real), ALLOCATABLE :: T_pseudo(:)
REAL(kind_real)              :: gamma
REAL(kind_real)              :: beta
REAL(kind_real)              :: c ! continuity constant for hydrostatic pressure
LOGICAL                      :: nonmon
LOGICAL                      :: unphys

! Allocate arrays for pseudo-level processing
IF (GPSRO_pseudo_ops) THEN
  nb_pseudo = 2 * nlevq - 1
  ALLOCATE (P_pseudo(nb_pseudo))
  ALLOCATE (q_pseudo(nb_pseudo))
  ALLOCATE (T_pseudo(nb_pseudo))
  ALLOCATE (z_pseudo(nb_pseudo))
  ALLOCATE (N_pseudo(nb_pseudo))
END IF

! Set up the P and q vectors from x
P(:) = x(1:nlevp)
q(:) = x(nlevp+1:nlevp+nlevq)

! Initialise refractivity arrays to missing Data
refrac(:) = missing_value(refrac(1))
T(:) = missing_value(T(1))
nonmon = .FALSE.
unphys = .FALSE.
refracerr = .FALSE.

DO i = 1, nlevp
  IF (P(i) == missing_value(P(i))) THEN  !pressure missing
    refracerr = .TRUE.
    WRITE(message, *) RoutineName, "Missing value P", i
    CALL fckit_log % warning(message)
    EXIT
  END IF
END DO

DO i = 1, nlevp-1
  IF (P(i) - P(i + 1) < 0.0) THEN       !or non-monotonic pressure
    refracerr = .TRUE.
    nonmon = .TRUE.
    WRITE(message,*) "Non monotonic", i, P(i), P(i+1)
    CALL fckit_log % warning(message)
    EXIT
  END IF
END DO

IF (ANY (P(:) <= 0.0)) THEN           !pressure zero or negative
  refracerr = .TRUE.
  unphys = .TRUE.
END IF

IF (ANY (q(:) < 0.0)) THEN           !humidity negative
  refracerr = .TRUE.
  unphys = .TRUE.
END IF

! only proceed if pressure is valid
IF (refracerr) THEN
  IF (nonmon) THEN
    CALL fckit_log%warning (RoutineName // ": GPSRO Pressure non-monotonic")
  ELSE IF (unphys) THEN
    CALL fckit_log%warning (RoutineName // ": GPSRO Pressure <= zero")
  ELSE
    CALL fckit_log%warning (RoutineName // ": GPSRO Pressure missing")
  END IF
ELSE

  ! Calculate exner on rho levels.
  Exner(:) = (P(:) / Pref) ** rd_over_cp

!  PRINT*, 'Exner'
!  WRITE(*, '(7E18.10)') exner

  ! Calculate the refractivity on the b levels
  search_lev = 1
  DO i = 1, nlevq
    ! Search for the first pressure (rho) level which has a height greater than
    ! the temperature (theta) level being considered.  Interpolate the pressure
    ! to the temperature levels using the level which has been found.
    DO this_lev = search_lev, nlevp
      IF (za(this_lev) > zb(i)) THEN
        EXIT
      END IF
    END DO
    IF (this_lev == 1) THEN
      ! Calc. pressure pb
      Pb(i) = P(i)
      CALL fckit_log%warning (RoutineName // ": Bottom temperature level is below bottom pressure level!  Extrapolating!")
      CALL fckit_log%warning (RoutineName // ": Results could be very inaccurate")
    ELSE
      ! Calc. pressure pb
      pwt1 = (za(this_lev) - zb(i)) / (za(this_lev) - za(this_lev-1))
      search_lev = this_lev
      pwt2 = 1.0 - pwt1

      ! Calculate the pressure on the theta level.
      IF (GPSRO_vert_interp_ops) THEN
        ! Assume ln(P) linear with height
        Pb(i) = EXP (pwt1 * LOG (P(this_lev-1)) + pwt2 * LOG (P(this_lev)))
      ELSE
        ! Assume Exner varies linearly with height
        Pb(i) = Pref * (pwt1 * (P(this_lev-1) / Pref) ** rd_over_cp + pwt2 * &
          (P(this_lev) / Pref) ** rd_over_cp) ** (1.0 / rd_over_cp)
      END IF
    END IF

    ! Calculate the Exner on the theta level.
    Ex_theta = (Pb(i) / Pref) ** rd_over_cp

    ! Calculate mean layer Tv (virtual temperature) using ND definition
    Tv = grav * (za(this_lev) - za(this_lev-1)) * Ex_theta / &
        (Cp * (Exner(this_lev-1) - Exner(this_lev)))

    IF (i > nlevq) THEN

      T(i) = Tv

      ! No wet component

      Nwet = 0.0

    ELSE

      T(i) = Tv / (1.0 + C_virtual * q(i))

      ! Wet component

      Nwet = n_beta * Pb(i) * q(i) / (T(i) ** 2 * (mw_ratio + (1.0 - mw_ratio) * q(i)))

    END IF

    Ndry = n_alpha * Pb(i) / T(i)

    refrac(i) = Ndry + Nwet

  END DO

  ! Do pseudo-level processing
  IF (GPSRO_pseudo_ops) THEN
    counter = 1
    DO i = 1, nb_pseudo

      ! Odd 'i' (i.e. copies of actual model level values)
      IF (MOD (i, 2) > 0) THEN
        z_pseudo(i) = zb(counter)
        P_pseudo(i) = Pb(counter)
        q_pseudo(i) = q(counter)
        T_pseudo(i) = T(counter)
        counter = counter + 1

      ! Even 'i' (i.e. intermediate pseudo-levels)
      ELSE
        z_pseudo(i) = (zb(counter - 1) + zb(counter)) / 2.0

        ! Assume exponential variation when humidities are positive
        IF (MIN (q(counter - 1), q(counter)) > 0.0) THEN
          gamma = LOG (q(counter - 1) / q(counter)) / (zb(counter) - zb(counter - 1))
          q_pseudo(i) = q(counter - 1) * EXP (-gamma * (z_pseudo(i) - z_pseudo(i - 1)))

        ! Assume linear variation if humidities are -ve
        ELSE
          q_pseudo(i) = q(counter - 1) + (q(counter) - q(counter - 1)) / (zb(counter) - &
                        zb(counter - 1)) * (z_pseudo(i) - zb(counter - 1))
        END IF

        ! T varies linearly with height
        beta = (T(counter) - T(counter - 1)) / (zb(counter) - zb(counter - 1))
        T_pseudo(i) = T(counter - 1) + beta * (z_pseudo(i) - zb(counter - 1))

        ! Pressure varies to maintain hydrostatic balance
        IF (ABS (T(counter) - T(counter - 1)) > 1.0E-10) THEN
          c = ((Pb(counter) / Pb(counter - 1)) * (T(counter) / T(counter - 1)) ** (grav / (rd * beta)) - &
              1.0) / (zb(counter) - zb(counter - 1))
          P_pseudo(i) = (Pb(counter - 1) * (T_pseudo(i) / T(counter - 1)) ** &
                      (-grav / (rd * beta))) * (1.0 + c * (z_pseudo(i) - zb(counter - 1)))
        ELSE
          ! If layer is isothermal, explicitly force P to vary exponentially to avoid singularity
          P_pseudo(i) = Pb(counter - 1) * EXP (LOG (Pb(counter) / Pb(counter - 1)) * &
                      ((z_pseudo(i) - zb(counter - 1)) / (zb(counter) - zb(counter - 1))))
        END IF
      END IF
    END DO

    N_pseudo = n_alpha * P_pseudo / T_pseudo + n_beta * P_pseudo * q_pseudo / &
               (T_pseudo ** 2 * (mw_ratio + (1.0 - mw_ratio) * q_pseudo))
  END IF

END IF

IF (GPSRO_pseudo_ops) THEN
  IF (ALLOCATED (P_pseudo)) DEALLOCATE (P_pseudo)
  IF (ALLOCATED (q_pseudo)) DEALLOCATE (q_pseudo)
  IF (ALLOCATED (T_pseudo)) DEALLOCATE (T_pseudo)
END IF

END SUBROUTINE Ops_GPSRO_refrac

!-------------------------------------------------------------------------
end module ufo_gnssro_ukmo1d_utils_mod
