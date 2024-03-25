! (c) copyright 2017-2018 ucar
!
! this software is licensed under the terms of the apache licence version 2.0
! which can be obtained at http://www.apache.org/licenses/license-2.0.

!> fortran module to handle aod observations

MODULE ufo_aodluts_mod

  USE fckit_configuration_module, ONLY: fckit_configuration
  USE iso_c_binding
  USE kinds
  USE missing_values_mod

  USE ufo_geovals_mod, ONLY: ufo_geovals, ufo_geoval, ufo_geovals_get_var
  USE ufo_vars_mod
  USE ufo_crtm_utils_mod, ONLY: assign_aerosol_names, max_string, upper2lower
  USE ufo_luts_utils_mod, ONLY: luts_conf, luts_conf_setup, &
       &luts_conf_delete, calculate_aero_layers
  USE crtm_module
  USE crtm_spccoeff, ONLY: sc
  USE obsspace_mod

  USE cf_mieobs_mod, ONLY: get_cf_aod,get_rc_wavelengths

  IMPLICIT NONE
  PRIVATE

!> fortran derived type for aod trajectory
  TYPE, PUBLIC :: ufo_aodluts
     PRIVATE
     CHARACTER(len=maxvarlen), PUBLIC, ALLOCATABLE :: varin(:)  ! variablesrequested from the model
     INTEGER, ALLOCATABLE                          :: channels(:)
     INTEGER :: n_aerosols
     TYPE(luts_conf) :: conf
   CONTAINS
     PROCEDURE :: setup  => ufo_aodluts_setup
     PROCEDURE :: delete => ufo_aodluts_delete
     PROCEDURE :: simobs => ufo_aodluts_simobs
  END TYPE ufo_aodluts
  CHARACTER(len=maxvarlen), DIMENSION(4), PARAMETER :: varin_default = &
       (/var_ts, var_q, var_prs, var_prsi/)

  CHARACTER(maxvarlen), PARAMETER :: varname_tmplate="aerosol_optical_depth"

CONTAINS

! ------------------------------------------------------------------------------

  SUBROUTINE ufo_aodluts_setup(self, f_confoper, channels)

    IMPLICIT NONE
    CLASS(ufo_aodluts),        INTENT(inout) :: self
    TYPE(fckit_configuration), INTENT(in)    :: f_confoper
    INTEGER(c_int),            INTENT(in)    :: channels(:)  !list of channels to use

    INTEGER :: nvars_in, rc
    CHARACTER(len=max_string) :: err_msg
    TYPE(fckit_configuration) :: f_confopts

    CHARACTER(len=maxvarlen), ALLOCATABLE :: var_aerosols(:)

    CALL f_confoper%get_or_die("obs options",f_confopts)

    CALL luts_conf_setup(self%conf, f_confopts, f_confoper)

    CALL assign_aerosol_names(self%conf%aerosol_option,var_aerosols)

    self%n_aerosols=SIZE(var_aerosols)
    nvars_in = SIZE(varin_default)+self%n_aerosols

    ALLOCATE(self%varin(nvars_in))
    self%varin(1:SIZE(varin_default)) = varin_default
    self%varin(SIZE(varin_default)+1:) = var_aerosols

    ALLOCATE(self%channels(SIZE(channels)))
    ALLOCATE(self%conf%wavelengths(SIZE(channels)))

    self%channels(:) = channels(:)

    DEALLOCATE(var_aerosols)

  END SUBROUTINE ufo_aodluts_setup

! ------------------------------------------------------------------------------

  SUBROUTINE ufo_aodluts_delete(self)

    IMPLICIT NONE
    CLASS(ufo_aodluts), INTENT(inout) :: self

    CALL luts_conf_delete(self%conf)

    IF (ALLOCATED(self%varin)) DEALLOCATE(self%varin)
    IF (ALLOCATED(self%channels)) DEALLOCATE(self%channels)
    IF (ALLOCATED(self%conf%wavelengths)) DEALLOCATE(self%conf%wavelengths)

  END SUBROUTINE ufo_aodluts_delete

! ------------------------------------------------------------------------------

  SUBROUTINE ufo_aodluts_simobs(self, geovals, obss, nvars, nlocs, hofx)

    IMPLICIT NONE
    CLASS(ufo_aodluts),       INTENT(inout) :: self
    TYPE(ufo_geovals),        INTENT(in) :: geovals
    INTEGER,                  INTENT(in) :: nvars, nlocs
    REAL(c_double),           INTENT(inout) :: hofx(nvars, nlocs)
    TYPE(c_ptr), VALUE,       INTENT(in) :: obss

! local variables
    CHARACTER(*), PARAMETER :: program_name = 'ufo_aodluts_mod.f90'
    CHARACTER(255) :: message, version
    INTEGER        :: err_stat, alloc_stat
    INTEGER        :: l, m, n, i
    TYPE(ufo_geoval), POINTER :: temp
    REAL(c_double) :: missing

    INTEGER :: n_profiles
    INTEGER :: n_layers
    INTEGER :: n_channels
    INTEGER :: n_aerosols

! define the "non-demoninational" arguments
    TYPE(crtm_channelinfo_type)             :: chinfo(self%conf%n_sensors)

    REAL(kind_real), ALLOCATABLE :: wavelengths_all(:)
    REAL(kind_real), ALLOCATABLE :: aero_layers(:,:,:),rh(:,:)
    
    CHARACTER(len=maxvarlen), ALLOCATABLE :: var_aerosols(:)

    INTEGER :: rc

    CALL assign_aerosol_names(self%conf%aerosol_option,var_aerosols)

    n_profiles = geovals%nlocs
    CALL ufo_geovals_get_var(geovals, var_aerosols(1), temp)
    n_layers = temp%nval
    NULLIFY(temp)

    n_aerosols=self%n_aerosols

    ALLOCATE(aero_layers(n_aerosols,n_layers,n_profiles),&
         &rh(n_layers,n_profiles))
    
    IF (self%conf%use_crtm) THEN
       err_stat = crtm_init( self%conf%sensor_id, &
            chinfo, &
            file_path=TRIM(self%conf%coefficient_path), &
            quiet=.TRUE.)
       
       IF ( err_stat /= success ) THEN
          message = 'error initializing crtm'
          CALL display_message( program_name, message, failure )
          STOP
       END IF
    ENDIF

    sensor_loop:DO n = 1, self%conf%n_sensors

       IF (ALLOCATED(wavelengths_all)) DEALLOCATE(wavelengths_all)
       
       IF (self%conf%use_crtm) THEN

          n_channels = crtm_channelinfo_n_channels(chinfo(n))
         
          ALLOCATE(wavelengths_all(n_channels), stat = alloc_stat)
          
          IF ( alloc_stat /= 0 ) THEN
             message = 'error allocating wavelengths_all'
             CALL display_message( program_name, message, failure )
             STOP
          END IF
          
          wavelengths_all=1.e7/sc(chinfo(n)%sensor_index)%wavenumber(:)

       ELSE

          CALL get_rc_wavelengths(self%conf%rcfile,wavelengths_all,rc)

          IF ( rc /= 0 ) THEN
             message = 'error getting wavelengths from rcfile'
             CALL display_message( program_name, message, failure )
             STOP
          END IF
          
       ENDIF

       self%conf%wavelengths=wavelengths_all(self%channels)

       CALL calculate_aero_layers(self%conf,&
            &n_aerosols, n_profiles, n_layers,&
            &geovals, aero_layers=aero_layers, rh=rh)
 
!note the difference in the "if" construct:
!a) aaod_tot denotes Absorption AOD 
!b) aod_tot denotes AOD
!that are alternatively output from 
!subroutine get_cf_aod in geos-aero/src/CF_MieObs.F90 

!aaod_tot option is chosen with the following configuration   

!- obs space:
!    name: Aod
!..........
!     channels: 3,5,6,7
!  obs operator:
!    name: AodLUTs
!    obs options:
!      Sensor_ID: aeronet
!      AerosolOption: aerosols_gocart_gefs
!      RCFile: [geosaod_aeronet.rc]
!      AbsorptionAod: true

!It can only exist for AERONET observations
!see 
!https://amt.copernicus.org/articles/13/3375/2020/#section4
!for details (note that only available for "channels": 3,5,6,7)
!i.e. 440,675,870,1020nm
!An illustrative example is given in
!fv3-jedi/test/testinput/hofx_gfs_aero.yaml

       IF (self%conf%aaod) THEN
          CALL get_cf_aod(n_layers, n_profiles, nvars, n_aerosols, &
               &self%conf%rcfile,  &
               &self%conf%wavelengths, var_aerosols, aero_layers, rh,&
               &aaod_tot = hofx, rc = rc)  
       ELSE
          CALL get_cf_aod(n_layers, n_profiles, nvars, n_aerosols, &
               &self%conf%rcfile,  &
               &self%conf%wavelengths, var_aerosols, aero_layers, rh,&
               &aod_tot = hofx, rc = rc)  
       ENDIF

       DEALLOCATE(aero_layers,rh,wavelengths_all)

       IF (rc /= 0) THEN
          message = 'error on exit from get_cf_aod'
          CALL display_message( program_name, message, failure )
          STOP
       END IF
          
    END DO sensor_loop
    
    IF (self%conf%use_crtm) THEN
       err_stat = crtm_destroy( chinfo )
       IF ( err_stat /= success ) THEN
          message = 'error destroying crtm (settraj)'
          CALL display_message( program_name, message, failure )
          STOP
       END IF
    ENDIF

  END SUBROUTINE ufo_aodluts_simobs

! ------------------------------------------------------------------------------

END MODULE ufo_aodluts_mod
