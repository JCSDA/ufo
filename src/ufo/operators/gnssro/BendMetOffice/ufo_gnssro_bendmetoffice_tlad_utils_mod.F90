!-------------------------------------------------------------------------------
! (C) British Crown Copyright 2020 Met Office
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
!-------------------------------------------------------------------------------
module ufo_gnssro_bendmetoffice_tlad_utils_mod

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
public             :: Ops_GPSROcalc_alphaK
public             :: Ops_GPSROcalc_nrK
private

contains

SUBROUTINE Ops_GPSROcalc_nrK (zb,       & ! geopotential heights of model levels
                              nb,       & ! number of levels in zb
                              Rad,      & ! radius of curvature of earth at observation
                              lat,      & ! latitude at observation
                              und,      & ! geoid undulation above WGS-84
                              refrac,   & ! refractivity of model on model levels
                              dnr_dref)   ! Calculated gradient of nr


USE ufo_gnssro_ukmo1d_utils_mod, only: Ops_GPSRO_geop_geom

IMPLICIT NONE

! Subroutine arguments:
INTEGER, INTENT(IN)          :: nb              ! number of levels in zb
REAL(kind_real), INTENT(IN)  :: zb(nb)          ! geopotential heights on zb levels /m
REAL(kind_real), INTENT(IN)  :: Rad             ! local radius of curvature of earth /m
REAL(kind_real), INTENT(IN)  :: lat             ! latitude at observation/ degrees
REAL(kind_real), INTENT(IN)  :: und             ! geoid undulation
REAL(kind_real), INTENT(IN)  :: refrac(nb)      ! refractivity on model levels / N
REAL(kind_real), INTENT(OUT) :: dnr_dref(nb,nb) ! Calculated gradient of nr wrt ref

! Local declarations:
CHARACTER(len=*), PARAMETER :: RoutineName = "Ops_GPSROcalc_nrK"
REAL(kind_real)             :: r(nb)           ! radius of model levels /m
REAL(kind_real)             :: z(nb)           ! geopotential heights on zb levels /m, local copy of zb
INTEGER                     :: i

!Initialise the matrix
dnr_dref = 0.0

!----------------------------------------------
! 1. Convert zb values to geometric altitudes
!---------------------------------------------
z= zb + und                 ! approx. convert to geopotential above WGS-84 ellipsoid
CALL Ops_GPSRO_geop_geom (lat, &
                          z)

!--------------------------------------------------
! 2. Calculate dnr/dref
!--------------------------------------------------

r = Rad + z

DO i = 1,nb
  IF (zb(i) > 0.0 .AND. refrac(i) > 0.0) THEN
    dnr_dref(i,i) = 1.0E-6 * r(i)
  END IF
END DO

END SUBROUTINE Ops_GPSROcalc_nrK

!-------------------------------------------------------------------------------
! GPSRO 1D bending angle operator K code.
!-------------------------------------------------------------------------------
SUBROUTINE Ops_GPSROcalc_alphaK (nobs,     &
                                 nlev,     &
                                 a,        &
                                 refrac,   &
                                 nr,       &
                                 Kmat_ref, &
                                 Kmat_nr,  &
                                 noCheck)

IMPLICIT NONE

! Subroutine arguments:
INTEGER, INTENT(IN)          :: nobs                ! size of ob. vector
INTEGER, INTENT(IN)          :: nlev                ! no. of refractivity levels
REAL(kind_real), INTENT(IN)  :: a(nobs)             ! observation impact parameters
REAL(kind_real), INTENT(IN)  :: refrac(nlev)        ! refractivity values on model levels
REAL(kind_real), INTENT(IN)  :: nr(nlev)            ! refractive index * radius product
REAL(kind_real), INTENT(OUT) :: Kmat_ref(nobs,nlev) ! BA gradient wrt refractivity
REAL(kind_real), INTENT(OUT) :: kmat_nr(nobs,nlev)  ! BA gradient wrt index * radius product
LOGICAL, INTENT(IN)          :: noCheck             ! If true, do not apply super-refraction check

! Local declarations:
CHARACTER(len=*), PARAMETER :: RoutineName = "Ops_GPSROcalc_alphaK"
INTEGER                     :: i
INTEGER                     :: n
INTEGER                     :: jbot
INTEGER                     :: kbot
REAL(kind_real)             :: kval(nlev - 1)      ! exponential decay rate between levels
REAL(kind_real)             :: root_2pia
REAL(kind_real)             :: ref_low             ! refractivity at lower level
REAL(kind_real)             :: nr_low              ! impact parameter lower level
REAL(kind_real)             :: tup
REAL(kind_real)             :: tlow                ! upper/lower bounds to error function
REAL(kind_real)             :: dalpha              ! delta bending angle
REAL(kind_real)             :: diff_erf            ! error function result
REAL(kind_real)             :: erf_tup
REAL(kind_real)             :: erf_tlow            ! error function of tup and tlow
REAL(kind_real)             :: t                   ! int. step in approx. to error function
REAL(kind_real), PARAMETER  :: a1 = 0.3480242      ! consts used in error function approx.
REAL(kind_real), PARAMETER  :: a2 = -0.0958798
REAL(kind_real), PARAMETER  :: a3 = 0.7478556
REAL(kind_real), PARAMETER  :: p = 0.47047
REAL(kind_real)             :: dkval_dref(nlev - 1,2)
REAL(kind_real)             :: dkval_dnr(nlev - 1,2)
REAL(kind_real)             :: dalpha_dref(2)
REAL(kind_real)             :: dalpha_dnr(2)
REAL(kind_real)             :: dalpha_dk
REAL(kind_real)             :: dalpha_drlow
REAL(kind_real)             :: dalpha_derf
REAL(kind_real)             :: dalpha_dnrlow
REAL(kind_real)             :: dnrlow_dref(2)
REAL(kind_real)             :: dnrlow_dnr(2)
REAL(kind_real)             :: drlow_dref(2)
REAL(kind_real)             :: drlow_dnr(2)
REAL(kind_real)             :: drlow_dk
REAL(kind_real)             :: derf_dref(2)
REAL(kind_real)             :: derf_dnr(2)
REAL(kind_real)             :: derf_dtup
REAL(kind_real)             :: derf_dtlow
REAL(kind_real)             :: dtup_dnr(2)
REAL(kind_real)             :: dtlow_dnr(2)
REAL(kind_real)             :: dtup_dref(2)
REAL(kind_real)             :: dtlow_dref(2)
REAL(kind_real)             :: dtup_dk
REAL(kind_real)             :: dtlow_dk
INTEGER                     :: indices(nlev)  ! The indices to use in the refractivity profile
INTEGER                     :: ngood          ! The number of refractivity values to use
REAL(kind_real)             :: minNR          ! Minimum value of nr found so far
INTEGER                     :: iMin           ! Minimum value of indices to use in profile

!-------------------------------------------------------------------------------
! Initialise the K matrices
!-------------------------------------------------------------------------------

Kmat_ref(:,:) = 0.0
Kmat_nr(:,:) = 0.0
dkval_dref(:,:) = 0.0
dkval_dnr(:,:)  = 0.0

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

  IF (kval(i) > 1.0E-6) THEN
    dkval_dref(i,1) =  1.0 / (refrac(indices(i)) * MAX(1.0, (nr(indices(i+1)) - nr(indices(i)))))
    dkval_dref(i,2) = -1.0 / (refrac(indices(i+1)) * MAX(1.0, (nr(indices(i+1)) - nr(indices(i)))))

    dkval_dnr(i,1) =  kval(i) / MAX(1.0, (nr(indices(i+1)) - nr(indices(i))))
    dkval_dnr(i,2) = -kval(i) / MAX(1.0, (nr(indices(i+1)) - nr(indices(i))))
  ELSE
    kval(i) = 1.0E-6
  END IF

END DO

!-------------------------------------------------------------------------------
! Calculate the bending angle gradients
!-------------------------------------------------------------------------------

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

  tlow = 0.0

  DO i = iMin, nGood - 1

    ! initialise matrices
    !---------------------

    dalpha_dref(:) = 0.0
    dalpha_dnr(:) = 0.0
    drlow_dref(:) = 0.0
    drlow_dnr(:) = 0.0
    dnrlow_dref(:) = 0.0
    dnrlow_dnr(:) = 0.0
    dtup_dnr(:) = 0.0
    dtup_dref(:) = 0.0
    dtlow_dnr(:) = 0.0
    dtlow_dref(:) = 0.0
    derf_dtup = 0.0
    derf_dtlow = 0.0
    dalpha_drlow = 0.0
    dalpha_dnrlow = 0.0
    dalpha_dk = 0.0
    dalpha_derf = 0.0
    drlow_dk = 0.0
    dtup_dk = 0.0
    dtlow_dk = 0.0

    ! Values of refractivity and impact parameter at lower level
    !-----------------------------------------------------------
    IF (i == iMin) THEN

      ref_low = refrac(indices(i)) * EXP (-kval(i) * (a(n) - nr(indices(i))))

      drlow_dref(1) =  ref_low / refrac(indices(i))
      drlow_dk      = -ref_low * (a(n) - nr(indices(i)))
      drlow_dnr(1)  =  ref_low * kval(i)

      nr_low = a(n)

      dnrlow_dnr(1) = 0.0

    ELSE

      ref_low = refrac(indices(i))

      drlow_dref(1) = 1.0

      nr_low = nr(indices(i))

      dnrlow_dnr(1) = 1.0

    END IF

    drlow_dref(:) = drlow_dref(:) + drlow_dk * dkval_dref(i,:)
    drlow_dnr(:) = drlow_dnr(:) + drlow_dk * dkval_dnr(i,:)


    ! Limits used in the error function
    !----------------------------------
    IF (i == nGood - 1) THEN

      ! simple extrapolation 100km above the uppermost level.
      !-----------------------------------------------------
      tup = SQRT(kval(i) * (nr(indices(i+1)) + 1.0E5 - a(n)))

      dtup_dk = 0.5 * (nr(indices(i+1)) + 1.0E5 - a(n)) / tup
      dtup_dnr(2) = 0.5 * kval(i) / tup

    ELSE

      tup = SQRT (kval(i) * (nr(indices(i+1)) - a(n)))

      dtup_dk = 0.5 * (nr(indices(i+1)) - a(n)) / tup
      dtup_dnr(2) = 0.5 * kval(i) / tup

    END IF

    dtup_dref(:) = dtup_dref(:) + dtup_dk * dkval_dref(i,:)
    dtup_dnr(:) = dtup_dnr(:) + dtup_dk * dkval_dnr(i,:)

    tlow = 0.0

    IF (i > iMin) THEN
      tlow = SQRT(kval(i) * (nr(indices(i)) - a(n)))
      dtlow_dk = 0.5 * (nr(indices(i)) - a(n)) / tlow
      dtlow_dnr(1) = 0.5 * kval(i) / tlow
    END IF

    dtlow_dref(:) = dtlow_dref(:) + dtlow_dk * dkval_dref(i,:)
    dtlow_dnr(:) = dtlow_dnr(:) + dtlow_dk * dkval_dnr(i,:)

    ! Abramowitz and Stegun approx. to error function
    !------------------------------------------------
    t = 1.0 / (1.0 + p * tup)
    erf_tup = 1.0 - EXP (-(tup ** 2))  * (a1 + (a2 + a3 * t) * t) * t
    derf_dtup = EXP (-(tup ** 2)) * ((2.0 * tup * (a1 + (a2 + a3 * t) * t) * t) + &
                                 ((a1 + (2.0 * a2 + 3.0 * a3 * t) * t) * p * t ** 2))

    t = 1.0 / (1.0 + p * tlow)
    erf_tlow = 1.0 - EXP (-(tlow ** 2)) * (a1 + (a2 + a3 * t) * t) * t

    ! Multiplied by -1.0 to account for the usage in derf_dref and derf_dnr
    !---------------------------------------------------------------------
    derf_dtlow = -1.0 * EXP (-(tlow ** 2)) * ((2.0 * tlow * (a1 + (a2 + a3 * t) * t) * t) + &
                              ((a1 + (2.0 * a2 + 3.0 * a3 * t) * t) * p * t ** 2))

    diff_erf = erf_tup - erf_tlow

    derf_dref(:) = derf_dtup * dtup_dref(:) + &
                    derf_dtlow * dtlow_dref(:)

    derf_dnr(:) = derf_dtup * dtup_dnr(:) + &
                   derf_dtlow * dtlow_dnr(:)

    dalpha = 1.0E-6 * Root_2PIa * SQRT (kval(i)) * &
                   ref_low * EXP (kval(i) * (nr_low - a(n))) * diff_erf

    dalpha_drlow = dalpha / MAX (1.0E-10,ref_low)

    dalpha_derf = dalpha / MAX (1.0E-10,diff_erf)

    dalpha_dnrlow = dalpha * kval(i)

    dalpha_dk = dalpha * (nr_low - a(n) + 0.5 / kval(i))

    ! Now apply chain rule
    !----------------------

    dalpha_dref(:) = dalpha_dref(:) + &
                     dalpha_drlow * drlow_dref(:) + &
                     dalpha_derf * derf_dref(:) + &
                     dalpha_dnrlow * dnrlow_dref(:) + &
                     dalpha_dk * dkval_dref(i,:)

    dalpha_dnr(:) = dalpha_dnr(:) + &
                    dalpha_drlow * drlow_dnr(:) + &
                    dalpha_derf * derf_dnr(:) + &
                    dalpha_dnrlow * dnrlow_dnr(:) + &
                    dalpha_dk * dkval_dnr(i,:)

    ! Now update matrices
    !---------------------

    Kmat_ref(n,indices(i)) = Kmat_ref(n,indices(i)) + dalpha_dref(1)
    Kmat_nr(n,indices(i)) = Kmat_nr(n,indices(i)) + dalpha_dnr(1)

    Kmat_ref(n,indices(i+1)) = Kmat_ref(n,indices(i+1)) + dalpha_dref(2)
    Kmat_nr(n,indices(i+1)) = Kmat_nr(n,indices(i+1)) + dalpha_dnr(2)

  END DO

END DO

END SUBROUTINE Ops_GPSROcalc_alphaK

end module ufo_gnssro_bendmetoffice_tlad_utils_mod
