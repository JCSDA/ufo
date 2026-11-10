! (c) copyright 2018 ucar
!
! this software is licensed under the terms of the apache licence version 2.0
! which can be obtained at http://www.apache.org/licenses/license-2.0.

!> fortran module to provide code shared between nonlinear and tlm/adm radiance calculations

MODULE ufo_luts_utils_mod

  USE fckit_configuration_module, ONLY: fckit_configuration
  USE iso_c_binding
  USE kinds

  USE crtm_module

  USE ufo_vars_mod
  USE ufo_geovals_mod, ONLY: ufo_geovals, ufo_geoval, ufo_geovals_get_var
  USE obsspace_mod

  USE ufo_crtm_utils_mod, ONLY: assign_aerosol_names => ufo_crtm_assign_aerosol_names, &
       &ufo_crtm_sensor_wavelengths
  USE cf_mieobs_mod, ONLY: get_rc_wavelengths
  USE ufo_utils_mod, ONLY: upper2lower
  USE ufo_constants_mod, ONLY: grav, rd_over_rv, c_virtual, t0c,&
       &aerosol_concentration_minvalue,&
       &aerosol_concentration_minvalue_layer

  IMPLICIT NONE
  PRIVATE

  PUBLIC luts_conf
  PUBLIC luts_conf_setup
  PUBLIC luts_conf_delete
  PUBLIC luts_set_wavelengths
  PUBLIC :: calculate_aero_layers

  INTEGER, PARAMETER, PUBLIC :: max_string=800

! Saturation vapor pressure table for qsmith: es over water every 1/es_table_per_k K,
! starting es_table_below_triple K below the triple point (113.16 K to 375.16 K)
  INTEGER, PARAMETER :: es_table_per_k = 10
  INTEGER, PARAMETER :: es_table_n = 2621
  REAL,    PARAMETER :: es_table_below_triple = 160.
! es table and its forward differences, built on first use
  REAL, ALLOCATABLE, SAVE :: es_table(:), es_table_diff(:)


!type for general config
  TYPE luts_conf
     CHARACTER(len=maxvarlen) :: aerosol_option
     LOGICAL :: aaod
     CHARACTER(len=max_string) :: rcfile
     CHARACTER(len=255), ALLOCATABLE :: sensor_id(:)
     INTEGER :: n_sensors
     REAL(kind_real), ALLOCATABLE :: wavelengths(:)
     LOGICAL :: use_crtm
     CHARACTER(len=255) :: endian_type
     CHARACTER(len=255) :: coefficient_path
     REAL(kind_real) :: convert_factor_model
     LOGICAL :: dry_mixr_model
  END TYPE luts_conf

CONTAINS

! ------------------------------------------------------------------------------

  SUBROUTINE luts_conf_setup(conf, f_confopts, f_confoper)

    IMPLICIT NONE
    TYPE(luts_conf),            INTENT(inout) :: conf
    TYPE(fckit_configuration),  INTENT(in)    :: f_confopts
    TYPE(fckit_configuration),  INTENT(in)    :: f_confoper

    CHARACTER(len=:), ALLOCATABLE :: str

    CALL f_confopts%get_or_die("AerosolOption",str)
    conf%aerosol_option = upper2lower(str)
    CALL f_confopts%get_or_die("RCFile",str)
    conf%rcfile = str
    conf%n_sensors = 1
    ALLOCATE(conf%sensor_id(conf%n_sensors))
    CALL f_confopts%get_or_die("Sensor_ID",str)
    conf%sensor_id(conf%n_sensors) = str

    IF (f_confOpts%has("CoefficientPath")) THEN
       CALL f_confopts%get_or_die("CoefficientPath",str)
       conf%coefficient_path = str
       CALL f_confopts%get_or_die("EndianType",str)
       conf%endian_type = str
       conf%use_crtm=.TRUE.
    ELSE
       conf%use_crtm=.FALSE.
    END IF

    IF (f_confOpts%has("model units coeff")) THEN
       CALL f_confopts%get_or_die("model units coeff",conf%convert_factor_model)
    END IF

    IF (f_confOpts%has("dry mix ratio")) THEN
       CALL f_confopts%get_or_die("dry mix ratio",conf%dry_mixr_model)
    END IF

    IF (f_confOpts%has("AbsorptionAod")) THEN
       CALL f_confopts%get_or_die("AbsorptionAod", conf%aaod)
    END IF

    DEALLOCATE(str)

  END SUBROUTINE luts_conf_setup

! -----------------------------------------------------------------------------

  SUBROUTINE luts_conf_delete(conf)

    IMPLICIT NONE
    TYPE(luts_conf), INTENT(inout) :: conf

    IF (ALLOCATED(conf%sensor_id)) DEALLOCATE(conf%sensor_id)

  END SUBROUTINE luts_conf_delete

! ------------------------------------------------------------------------------

!> Set conf%wavelengths for the channels used, from the sensor's CRTM SpcCoeff file
!> (use_crtm) or from the rc file. Reads the file directly, with no CRTM_Init.
  SUBROUTINE luts_set_wavelengths(conf, channels)

    IMPLICIT NONE
    TYPE(luts_conf), INTENT(inout) :: conf
    INTEGER,         INTENT(in)    :: channels(:)

    REAL(kind_real), ALLOCATABLE :: wavelengths_all(:)
    INTEGER :: rc

    IF (conf%use_crtm) THEN
       CALL ufo_crtm_sensor_wavelengths(conf%coefficient_path, conf%sensor_id(1), &
            &wavelengths_all)
    ELSE
       CALL get_rc_wavelengths(conf%rcfile, wavelengths_all, rc)
       IF (rc /= 0) CALL abor1_ftn("luts_set_wavelengths: error getting wavelengths from rcfile")
    END IF

    conf%wavelengths = wavelengths_all(channels)

  END SUBROUTINE luts_set_wavelengths

! ------------------------------------------------------------------------------

  SUBROUTINE calculate_aero_layers(conf,&
       &n_aerosols, n_profiles, n_layers, &
       &geovals, aero_layers, rh, layer_factors)

    IMPLICIT NONE

    TYPE(luts_conf), INTENT(in) :: conf
    INTEGER, INTENT(in) :: n_aerosols, n_profiles, n_layers
    TYPE(ufo_geovals), INTENT(in) :: geovals
    REAL(kind_real), OPTIONAL, INTENT(out) :: &
         &aero_layers(n_aerosols,n_layers,n_profiles)
    REAL(kind_real), OPTIONAL, INTENT(out) :: rh(n_layers,n_profiles)
    REAL(kind_real), OPTIONAL, INTENT(out) :: layer_factors(n_layers,n_profiles)

! local variables
    INTEGER :: m, ivar,k
    TYPE(ufo_geoval), POINTER :: geoval
    CHARACTER(max_string) :: err_msg
    REAL(kind_real), DIMENSION(n_layers,n_profiles) :: t,sphum,pmid
    REAL(kind_real), DIMENSION(n_layers+1,n_profiles) :: pint
    CHARACTER(len=maxvarlen), ALLOCATABLE :: var_aerosols(:)
    REAL(kind_real) :: factors(n_layers,n_profiles)

    CALL ufo_geovals_get_var(geovals, var_ts, geoval)
    IF (geoval%nval /= n_layers) THEN
       WRITE(err_msg,*) "get_atm_aero_data error: layers inconsistent!"
       CALL abor1_ftn(err_msg)
    END IF
    t=geoval%vals

    CALL ufo_geovals_get_var(geovals, var_prs, geoval)
    pmid=geoval%vals

    CALL ufo_geovals_get_var(geovals, var_prsi, geoval)
    pint=geoval%vals

    CALL ufo_geovals_get_var(geovals, var_q, geoval)
    sphum=geoval%vals

!three versions of GOCART:

!1) the original GOCART parameterization that exists
!   in CRTM (bc1-2,oc1-2,sulf,dust1-5,seas1-4)

!2) GOCART parameterization
!(bc1-2,oc1-2,sulf,dust1-5,seas1-5) that was implemented
!in NOAA's GEFS-Aerosols model

!3) GOCART parameterization
!(bc1-2,oc1-2,sulf,dust1-5,seas1-5,nitrate1-3) that was implemented
!in NOAA's UFS-Aerosols model

    IF (TRIM(conf%aerosol_option) /= "aerosols_gocart_gefs" .AND. &
         &TRIM(conf%aerosol_option) /= "aerosols_gocart_ufs") THEN
       WRITE(err_msg,*) "this aerosol not implemented - check later"
       CALL abor1_ftn(err_msg)
    END IF

    CALL assign_aerosol_names(conf%aerosol_option,var_aerosols)

    IF (conf%dry_mixr_model) THEN

       DO k=1,n_layers
          DO m = 1, n_profiles
             ! Convert to a dry-air mass basis. delta_p/g is the moist-air mass, so use
             ! the mass fraction 1/(1+mixing_ratio), not the dry partial-pressure ratio
             ! of Cotton eq. (2.4), p_dry=p_total/(1+r_v/r_d*mixing_ratio) (ufo#4450).
             factors(k,m)=(pint(k+1,m)-pint(k,m))/grav/&
                  &(1.0_kind_real+sphum(k,m)/(1.0_kind_real-sphum(k,m)))
          END DO
       END DO

    ELSE

       DO k=1,n_layers
          DO m = 1, n_profiles
             factors(k,m)=(pint(k+1,m)-pint(k,m))/grav
          END DO
       END DO

    END IF

    IF ( PRESENT(aero_layers) ) THEN
       DO ivar=1,n_aerosols
          CALL ufo_geovals_get_var(geovals, var_aerosols(ivar), geoval)
          aero_layers(ivar,:,:)=conf%convert_factor_model*geoval%vals*factors
       END DO
    END IF

    IF ( PRESENT(rh) ) THEN
       CALL qsmith(t,sphum,pmid,rh)
       WHERE (rh > 1.0_kind_real) rh=1.0_kind_real
    END IF

    IF ( PRESENT(layer_factors) ) layer_factors=conf%convert_factor_model*factors

  END SUBROUTINE calculate_aero_layers

! ------------------------------------------------------------------------------

!> Relative humidity rh = q/qs (not capped at 1) from a saturation vapor pressure
!> lookup table (ported from GFDL FV3 qsmith, which returns qs). es(T) is linearly
!> interpolated in the table and qs = (rd/rv) * es * (1 + c_virtual*q) / p.
!>   t [K], sphum = specific humidity [kg/kg], p [Pa], all (n_layers, n_profiles).
  SUBROUTINE qsmith(t,sphum,p,rh)

    REAL(kind_real), DIMENSION(:,:), INTENT(in) :: t,sphum,p
    REAL(kind_real), DIMENSION(:,:), INTENT(out) :: rh

    REAL :: es, qs, q
    REAL :: ap1
    REAL :: Tmin
    INTEGER :: i, k, it, n_layers, n_profiles

    n_layers=SIZE(t,1)
    n_profiles=SIZE(t,2)

! The lookup starts at t0c, the table at the triple point (0.01 K apart), as in
! the original code
    Tmin = t0c-es_table_below_triple

    IF ( .NOT. ALLOCATED(es_table) ) CALL qsmith_init()

    DO k=1,n_layers
       DO i=1,n_profiles
          ap1 = REAL(es_table_per_k)*DIM(t(k,i), Tmin) + 1.
          ap1 = MIN(REAL(es_table_n), ap1)
          it = ap1
          es = es_table(it) + (ap1-it)*es_table_diff(it)
          q=sphum(k,i)
          qs = rd_over_rv*es*(1.+c_virtual*q)/p(k,i)
          rh(k,i) = q/qs
       END DO
    END DO

  END SUBROUTINE qsmith

! ------------------------------------------------------------------------------

!> Build the es table and its forward differences
!> es_table_diff(i) = es_table(i+1) - es_table(i), used for linear interpolation.
  SUBROUTINE qsmith_init()

    INTEGER :: i

    ALLOCATE(es_table(es_table_n), es_table_diff(es_table_n))

    CALL qs_table(es_table_n, es_table)

    DO i=1,es_table_n-1
       es_table_diff(i) = es_table(i+1) - es_table(i)
    END DO
    es_table_diff(es_table_n) = es_table_diff(es_table_n-1)

  END SUBROUTINE qsmith_init

! ------------------------------------------------------------------------------

!> Saturation vapor pressure over water [Pa] on n points from es_table_below_triple K
!> below the triple point, every 1/es_table_per_k K (Goff-Gratch formula, Smithsonian
!> Meteorological Tables p. 350).
  SUBROUTINE qs_table(n,table)

    INTEGER, INTENT(in):: n
    REAL, INTENT(out) :: table (n)

    REAL, PARAMETER :: esbasw = 1013246.0  ! es at the steam point (dyn cm**-2)
    REAL, PARAMETER :: tbasw  = 373.16     ! steam point of water (K)
    REAL, PARAMETER :: tbasi  = 273.16     ! triple point of water (K)
    REAL, PARAMETER :: dt = 1./es_table_per_k
    REAL :: Tmin, tem, aa, b, c, d, e
    INTEGER :: i

    Tmin = tbasi - es_table_below_triple
    DO  i=1,n
       tem = Tmin+dt*REAL(i-1)
       aa  = -7.90298*(tbasw/tem-1)
       b   =  5.02808*log10(tbasw/tem)
       c   = -1.3816e-07*(10**((1-tem/tbasw)*11.344)-1)
       d   =  8.1328e-03*(10**((tbasw/tem-1)*(-3.49149))-1)
       e   =  log10(esbasw)
       table(i)  = 0.1*10**(aa+b+c+d+e)  ! dyn cm**-2 to Pa
    END DO

  END SUBROUTINE qs_table

! ------------------------------------------------------------------------------

END MODULE ufo_luts_utils_mod

