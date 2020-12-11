! (c) copyright 2017-2018 ucar
!
! this software is licensed under the terms of the apache licence version 2.0
! which can be obtained at http://www.apache.org/licenses/license-2.0.

!> fortran module to handle tl/ad for aod observations

MODULE ufo_aodluts_tlad_mod

  USE fckit_configuration_module, ONLY: fckit_configuration
  USE iso_c_binding
  USE kinds
  USE missing_values_mod

  USE ufo_geovals_mod, ONLY: ufo_geovals, ufo_geoval, ufo_geovals_get_var
  USE ufo_vars_mod
  USE ufo_crtm_utils_mod, ONLY: assign_aerosol_names, max_string
  USE ufo_luts_utils_mod, ONLY: luts_conf, luts_conf_setup, &
       &luts_conf_delete, calculate_aero_layers
  USE crtm_module
  USE crtm_spccoeff, ONLY: sc
  USE obsspace_mod

  USE cf_mieobs_mod, ONLY: get_cf_aod
  USE geos_mieobs_mod, ONLY: get_geos_aod_tl, get_geos_aod_ad

  IMPLICIT NONE
  PRIVATE

!> fortran derived type for aod trajectory
  TYPE, PUBLIC :: ufo_aodluts_tlad
     PRIVATE
     CHARACTER(len=maxvarlen), PUBLIC, ALLOCATABLE :: varin(:)  ! variablesrequested from the model
     INTEGER, ALLOCATABLE                          :: channels(:)
     REAL(kind_real), ALLOCATABLE                  :: wavelengths(:)
     TYPE(luts_conf) :: conf
     INTEGER :: n_profiles
     INTEGER :: n_layers
     INTEGER :: n_channels
     INTEGER :: n_aerosols
     REAL(kind_real), ALLOCATABLE  :: bext(:,:,:,:)  
     REAL(kind_real), ALLOCATABLE  :: layer_factors(:,:)
     LOGICAL :: ltraj
   CONTAINS
     PROCEDURE :: setup  => ufo_aodluts_tlad_setup
     PROCEDURE :: delete  => ufo_aodluts_tlad_delete
     PROCEDURE :: settraj => ufo_aodluts_tlad_settraj
     PROCEDURE :: simobs_tl  => ufo_aodluts_simobs_tl
     PROCEDURE :: simobs_ad  => ufo_aodluts_simobs_ad
  END TYPE ufo_aodluts_tlad

CONTAINS

! ------------------------------------------------------------------------------

  SUBROUTINE ufo_aodluts_tlad_setup(self, f_confoper, channels)

    IMPLICIT NONE
    CLASS(ufo_aodluts_tlad),   INTENT(inout) :: self
    TYPE(fckit_configuration), INTENT(in)    :: f_confoper
    INTEGER(c_int),               INTENT(in)    :: channels(:)  !list of channels to use

    TYPE(fckit_configuration) :: f_confopts

    CHARACTER(len=maxvarlen), ALLOCATABLE :: var_aerosols(:)
    CHARACTER(len=:), ALLOCATABLE :: str

    CALL f_confoper%get_or_die("obs options",f_confopts)

    CALL luts_conf_setup(self%conf, f_confopts, f_confoper)

    CALL assign_aerosol_names(self%conf%aerosol_option,var_aerosols)

    self%n_aerosols=SIZE(var_aerosols)
    ALLOCATE(self%varin(self%n_aerosols))
    self%varin(1:self%n_aerosols) = var_aerosols

    ALLOCATE(self%channels(SIZE(channels)))
    ALLOCATE(self%wavelengths(SIZE(channels)))

    self%channels(:) = channels(:)

    DEALLOCATE(var_aerosols)

  END SUBROUTINE ufo_aodluts_tlad_setup

! ------------------------------------------------------------------------------

  SUBROUTINE ufo_aodluts_tlad_delete(self)

    IMPLICIT NONE
    CLASS(ufo_aodluts_tlad), INTENT(inout) :: self

    self%ltraj = .FALSE.
    CALL luts_conf_delete(self%conf)

!deallocate all arrays

    IF (ALLOCATED(self%bext)) DEALLOCATE(self%bext)
    IF (ALLOCATED(self%layer_factors)) DEALLOCATE(self%layer_factors)
    IF (ALLOCATED(self%varin)) DEALLOCATE(self%varin)
    IF (ALLOCATED(self%channels)) DEALLOCATE(self%channels)
    IF (ALLOCATED(self%wavelengths)) DEALLOCATE(self%wavelengths)

  END SUBROUTINE ufo_aodluts_tlad_delete

! ------------------------------------------------------------------------------

  SUBROUTINE ufo_aodluts_tlad_settraj(self, geovals, obss)

    IMPLICIT NONE

    CLASS(ufo_aodluts_tlad), INTENT(inout) :: self
    TYPE(ufo_geovals),        INTENT(in)    :: geovals
    TYPE(c_ptr), VALUE,       INTENT(in)    :: obss

! local variables
    CHARACTER(*), PARAMETER :: program_name = 'ufo_aodluts_tlad_mod.f90'
    CHARACTER(255) :: message, version
    INTEGER        :: err_stat, alloc_stat
    INTEGER        :: n,l,m
    TYPE(ufo_geoval), POINTER :: temp

! define the "non-demoninational" arguments
    TYPE(crtm_channelinfo_type)             :: chinfo(self%conf%n_sensors)

    CHARACTER(len=maxvarlen), ALLOCATABLE :: var_aerosols(:)
    REAL(kind_real), ALLOCATABLE :: aero_layers(:,:,:),rh(:,:)
    REAL(kind_real), ALLOCATABLE :: wavelengths_all(:)

    INTEGER :: rc,nvars

! get number of profile and layers from geovals

    CALL assign_aerosol_names(self%conf%aerosol_option,var_aerosols)

    self%n_profiles = geovals%nlocs
    CALL ufo_geovals_get_var(geovals, var_aerosols(1), temp)
    self%n_layers = temp%nval
    NULLIFY(temp)

    err_stat = crtm_init( self%conf%sensor_id, &
         chinfo, &
         file_path=TRIM(self%conf%coefficient_path), &
         quiet=.TRUE.)
    IF ( err_stat /= success ) THEN
       message = 'error initializing crtm (settraj)'
       CALL display_message( program_name, message, failure )
       STOP
    END IF

    ALLOCATE(aero_layers(self%n_aerosols,self%n_layers,self%n_profiles),&
         &rh(self%n_layers,self%n_profiles))

    ALLOCATE(self%layer_factors(self%n_layers,self%n_profiles))

    nvars=SIZE(self%channels)

    sensor_loop:DO n = 1, self%conf%n_sensors

       self%n_channels = crtm_channelinfo_n_channels(chinfo(n))

       IF (ALLOCATED(wavelengths_all)) DEALLOCATE(wavelengths_all)

       ALLOCATE(wavelengths_all(self%n_channels), stat = alloc_stat)

       wavelengths_all=1.e7/sc(chinfo(n)%sensor_index)%wavenumber(:)

       self%wavelengths=wavelengths_all(self%channels)

       CALL calculate_aero_layers(self%conf%aerosol_option,&
            &self%n_aerosols, self%n_profiles, self%n_layers,&
            &geovals, aero_layers=aero_layers, rh=rh, &
            &layer_factors=self%layer_factors)

       ALLOCATE(self%bext(self%n_layers, nvars, self%n_aerosols, self%n_profiles))

       CALL get_cf_aod(self%n_layers, self%n_profiles, nvars, &
            &self%n_aerosols, self%conf%rcfile,  &
            &self%wavelengths, var_aerosols, aero_layers, rh, &
            &ext=self%bext, rc = rc)

       IF (rc /= 0) THEN
          message = 'error on exit from get_cf_aod'
          CALL display_message( program_name, message, failure )
          STOP
       END IF

       DEALLOCATE(rh)
       DEALLOCATE(aero_layers)
       DEALLOCATE(wavelengths_all)

    END DO sensor_loop

    err_stat = crtm_destroy( chinfo )
    IF ( err_stat /= success ) THEN
       message = 'error destroying crtm (settraj)'
       CALL display_message( program_name, message, failure )
       STOP
    END IF


! set flag that the tracectory was set

    self%ltraj = .TRUE.

  END SUBROUTINE ufo_aodluts_tlad_settraj

! ------------------------------------------------------------------------------

  SUBROUTINE ufo_aodluts_simobs_tl(self, geovals, obss, nvars, nlocs, hofx)

    IMPLICIT NONE
    CLASS(ufo_aodluts_tlad), INTENT(in) :: self
    TYPE(ufo_geovals),        INTENT(in) :: geovals
    TYPE(c_ptr), VALUE,    INTENT(in)    :: obss
    INTEGER,                  INTENT(in)    :: nvars, nlocs
    REAL(c_double),        INTENT(inout) :: hofx(nvars, nlocs)

    CHARACTER(len=*), PARAMETER :: myname_="ufo_aodluts_simobs_tl"
    CHARACTER(max_string) :: err_msg
    INTEGER :: jaero
    TYPE(ufo_geoval), POINTER :: var_p

    CHARACTER(len=maxvarlen), ALLOCATABLE :: var_aerosols(:)

    REAL(kind_real), ALLOCATABLE :: aero_layers(:,:,:)

! initial checks

! check if trajectory was set
    IF (.NOT. self%ltraj) THEN
       WRITE(err_msg,*) myname_, ' trajectory wasnt set!'
       CALL abor1_ftn(err_msg)
    ENDIF

! check if nlocs is consistent in geovals & hofx
    IF (geovals%nlocs /= self%n_profiles) THEN
       WRITE(err_msg,*) myname_, ' error: nlocs inconsistent!'
       CALL abor1_ftn(err_msg)
    ENDIF

    CALL assign_aerosol_names(self%conf%aerosol_option,var_aerosols)

    CALL ufo_geovals_get_var(geovals, var_aerosols(1), var_p)

! check model levels is consistent in geovals & crtm
    IF (var_p%nval /= self%n_layers) THEN
       WRITE(err_msg,*) myname_, ' error: layers inconsistent!'
       CALL abor1_ftn(err_msg)
    ENDIF

    ALLOCATE(aero_layers(self%n_aerosols,self%n_layers,self%n_profiles))

    DO jaero=1,self%n_aerosols
       CALL ufo_geovals_get_var(geovals, var_aerosols(jaero), var_p)
       aero_layers(jaero,:,:)=var_p%vals*self%layer_factors
    ENDDO
    
    CALL get_geos_aod_tl(self%n_layers,self%n_profiles, nvars,&
         &self%n_aerosols, self%bext, aero_layers, aod_tot_tl=hofx)

    NULLIFY(var_p)

    DEALLOCATE(aero_layers)
    DEALLOCATE(var_aerosols)

  END SUBROUTINE ufo_aodluts_simobs_tl

! ------------------------------------------------------------------------------

  SUBROUTINE ufo_aodluts_simobs_ad(self, geovals, obss, nvars, nlocs, hofx)

    IMPLICIT NONE
    CLASS(ufo_aodluts_tlad), INTENT(in) :: self
    TYPE(ufo_geovals),     INTENT(inout) :: geovals
    TYPE(c_ptr), VALUE,    INTENT(in)    :: obss
    INTEGER,                  INTENT(in)    :: nvars, nlocs
    REAL(c_double),           INTENT(in) :: hofx(nvars, nlocs)

    CHARACTER(len=*), PARAMETER :: myname_="ufo_aodluts_simobs_ad"
    CHARACTER(max_string) :: err_msg
    TYPE(ufo_geoval), POINTER :: var_p
    INTEGER :: jaero

    REAL(kind_real), DIMENSION(:,:,:), ALLOCATABLE :: qm_ad
    CHARACTER(len=maxvarlen), ALLOCATABLE :: var_aerosols(:)
    REAL(kind_real), DIMENSION(:,:), ALLOCATABLE :: layer_factors


! initial checks

! check if trajectory was set
    IF (.NOT. self%ltraj) THEN
       WRITE(err_msg,*) myname_, ' trajectory wasnt set!'
       CALL abor1_ftn(err_msg)
    ENDIF

! check if nlocs is consistent in geovals & hofx
    IF (geovals%nlocs /= self%n_profiles) THEN
       WRITE(err_msg,*) myname_, ' error: nlocs inconsistent!'
       CALL abor1_ftn(err_msg)
    ENDIF

    CALL assign_aerosol_names(self%conf%aerosol_option,var_aerosols)

    ALLOCATE(qm_ad(self%n_aerosols, self%n_layers, nlocs))

    CALL get_geos_aod_ad(self%n_layers, nlocs, nvars, self%n_aerosols, &
         &self%bext, hofx, qm_ad)

    DO jaero=self%n_aerosols,1,-1

       CALL ufo_geovals_get_var(geovals, var_aerosols(jaero), var_p)
       IF (.NOT. ALLOCATED(var_p%vals)) THEN
          var_p%nlocs = self%n_profiles
          var_p%nval = self%n_layers
          ALLOCATE(var_p%vals(var_p%nval,var_p%nlocs))
          var_p%vals = 0.0_kind_real
       ENDIF

       qm_ad(jaero,:,:) = qm_ad(jaero,:,:) * self%layer_factors
       var_p%vals=qm_ad(jaero,:,:)

    ENDDO

    NULLIFY(var_p)

    DEALLOCATE(qm_ad)
    DEALLOCATE(var_aerosols)

    IF (.NOT. geovals%linit ) geovals%linit=.TRUE.

  END SUBROUTINE ufo_aodluts_simobs_ad

! ------------------------------------------------------------------------------

END MODULE ufo_aodluts_tlad_mod
