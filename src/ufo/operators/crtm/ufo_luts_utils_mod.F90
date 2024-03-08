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

  USE ufo_crtm_utils_mod, ONLY: qsmith, assign_aerosol_names, upper2lower
  USE ufo_crtm_utils_mod, ONLY: grav,rv_rd,&
       &aerosol_concentration_minvalue,&
       &aerosol_concentration_minvalue_layer

  IMPLICIT NONE
  PRIVATE

  PUBLIC luts_conf
  PUBLIC luts_conf_setup
  PUBLIC luts_conf_delete
  PUBLIC :: calculate_aero_layers

  INTEGER, PARAMETER, PUBLIC :: max_string=800


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

    CHARACTER(*), PARAMETER :: routine_name = 'luts_conf_setup'
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
    ENDIF

    IF (f_confOpts%has("model units coeff")) THEN
       CALL f_confopts%get_or_die("model units coeff",conf%convert_factor_model)
    ENDIF

    IF (f_confOpts%has("dry mix ratio")) THEN
       CALL f_confopts%get_or_die("dry mix ratio",conf%dry_mixr_model)
    ENDIF

    IF (f_confOpts%has("AbsorptionAod")) THEN
       CALL f_confopts%get_or_die("AbsorptionAod", conf%aaod)
    ENDIF

    DEALLOCATE(str)

  END SUBROUTINE luts_conf_setup

! -----------------------------------------------------------------------------

  SUBROUTINE luts_conf_delete(conf)

    IMPLICIT NONE
    TYPE(luts_conf), INTENT(inout) :: conf

    IF (ALLOCATED(conf%sensor_id)) DEALLOCATE(conf%sensor_id)

  END SUBROUTINE luts_conf_delete

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
    CHARACTER(len=maxvarlen) :: varname
    REAL(kind_real), DIMENSION(n_layers,n_profiles) :: t,sphum,pmid
    REAL(kind_real), DIMENSION(n_layers+1,n_profiles) :: pint
    CHARACTER(len=maxvarlen), ALLOCATABLE :: var_aerosols(:)
    REAL(kind_real) :: factors(n_layers,n_profiles)

    CALL ufo_geovals_get_var(geovals, var_ts, geoval)
    IF (geoval%nval /= n_layers) THEN
       WRITE(err_msg,*) 'get_atm_aero_data error: layers inconsistent!'
       CALL abor1_ftn(err_msg)
    ENDIF
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
       WRITE(err_msg,*) 'this aerosol not implemented - check later'
       CALL abor1_ftn(err_msg)
    ENDIF

    CALL assign_aerosol_names(conf%aerosol_option,var_aerosols)

    IF (conf%dry_mixr_model) THEN

       DO k=1,n_layers
          DO m = 1, n_profiles
!correct for mixing ratio factor 
!being calculated from dry pressure, cotton eq. (2.4)
!p_dry=p_total/(1+r_v/r_d*mixing_ratio)
             factors(k,m)=(pint(k+1,m)-pint(k,m))/grav/&
                  &(1.0_kind_real+rv_rd*sphum(k,m)/(1.0_kind_real-sphum(k,m)))
          ENDDO
       ENDDO

    ELSE
   
       DO k=1,n_layers
          DO m = 1, n_profiles
             factors(k,m)=(pint(k+1,m)-pint(k,m))/grav
          ENDDO
       ENDDO

    ENDIF

    IF ( PRESENT(aero_layers) ) THEN
       DO ivar=1,n_aerosols
          CALL ufo_geovals_get_var(geovals, var_aerosols(ivar), geoval)
          aero_layers(ivar,:,:)=conf%convert_factor_model*geoval%vals*factors
       ENDDO
    ENDIF

    IF ( PRESENT(rh) ) THEN
       CALL qsmith(t,sphum,pmid,rh)
       WHERE (rh > 1.0_kind_real) rh=1.0_kind_real
    ENDIF

    IF ( PRESENT(layer_factors) ) layer_factors=conf%convert_factor_model*factors

  END SUBROUTINE calculate_aero_layers

! ------------------------------------------------------------------------------

END MODULE ufo_luts_utils_mod

