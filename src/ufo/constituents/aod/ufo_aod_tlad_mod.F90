! (C) Copyright 2017-2018 UCAR
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 

!> Fortran module to handle aod observations

MODULE ufo_aod_tlad_mod
  USE iso_c_binding
  USE ufo_vars_mod
  USE ufo_locs_mod
  USE ufo_geovals_mod
  USE kinds  
  USE CRTM_module 
  USE ufo_aod_utils_mod
  USE ufo_basis_tlad_mod, ONLY: ufo_basis_tlad
  USE obsspace_mod
  
  IMPLICIT NONE
  PRIVATE
  
!> Fortran derived type for aod trajectory
  TYPE, extends(ufo_basis_tlad), PUBLIC :: ufo_aod_tlad
    PRIVATE
    TYPE(aod_conf) :: rc
    INTEGER :: n_Profiles
    INTEGER :: n_Layers
    INTEGER :: n_Channels
    TYPE(CRTM_atmosphere_type), ALLOCATABLE :: atm_k(:,:)
  CONTAINS
    PROCEDURE :: delete  => ufo_aod_tlad_delete
    PROCEDURE :: settraj => ufo_aod_tlad_settraj 
    PROCEDURE :: simobs_tl  => ufo_aod_simobs_tl
    PROCEDURE :: simobs_ad  => ufo_aod_simobs_ad
  END TYPE ufo_aod_tlad

  CHARACTER(len=MAXVARLEN) :: varname
  CHARACTER(max_string) :: message, version

CONTAINS
  
  SUBROUTINE ufo_aod_tlad_setup(self)
!optimally
!  SUBROUTINE ufo_aod_tlad_setup(self, c_conf)
    
    IMPLICIT NONE
    class(ufo_aod_tlad), INTENT(inout) :: self
!    TYPE(c_ptr),              INTENT(in)    :: c_conf
    
!    CALL aod_conf_setup(self%rc,c_conf)
 
!    self%rc
   
  END SUBROUTINE ufo_aod_tlad_setup
  

  SUBROUTINE ufo_aod_tlad_delete(self)
    IMPLICIT NONE
    class(ufo_aod_tlad), INTENT(inout)  :: self
    
! Nothing here yet
    
  END SUBROUTINE ufo_aod_tlad_delete


  SUBROUTINE ufo_aod_tlad_settraj(self, geovals, obss)

    IMPLICIT NONE
    class(ufo_aod_tlad), INTENT(inout) :: self
    TYPE(ufo_geovals),   INTENT(in)    :: geovals
    TYPE(c_ptr), value,  INTENT(in)    :: obss
    CHARACTER(len=*), PARAMETER :: program_name = "ufo_aod_tlad_settraj"

    TYPE(CRTM_channelinfo_type)             :: chinfo(self%rc%n_sensors)
    TYPE(CRTM_atmosphere_type), ALLOCATABLE :: atm(:)
    TYPE(CRTM_rtsolution_type), ALLOCATABLE :: rts(:,:)
    TYPE(CRTM_RTSolution_type), ALLOCATABLE :: rts_k(:,:)

    TYPE(ufo_geoval), POINTER :: geoval
    INTEGER        :: err_stat, alloc_stat    
    INTEGER :: l,m,n

    CONTINUE

    IF (self%rc%n_aerosols /= n_aerosols_gocart_default) THEN
       message = 'Only default GOCART with 14 species allowed for now'
       CALL display_message( program_name, message, failure )
       STOP
    ENDIF
       
    self%n_profiles = geovals%nobs
    CALL ufo_geovals_get_var(geovals, var_t, geoval)
    self%n_layers = geoval%nval
    NULLIFY(geoval)

    CALL CRTM_version( version )
    CALL program_message( program_name, &
         'Check/example program for the CRTM forward and &
         &k-matrix (settraj) functions using '//&
         &TRIM(self%rc%endian_type)//' coefficient datafiles',&
         &'CRTM version: '//TRIM(version) )
    WRITE( *,'(/5x,"initializing the CRTM (settraj) ...")' )
    err_stat = CRTM_init( self%rc%sensor_id, &
         chinfo, &
         file_path=TRIM(self%rc%coefficient_path), &
         quiet=.TRUE.)

    IF ( err_stat /= success ) THEN
       message = 'error initializing CRTM (settraj)'
       CALL display_message( program_name, message, failure )
       STOP
    END IF
    
    sensor_loop:DO n = 1, self%rc%n_sensors

       self%n_channels = CRTM_channelinfo_n_channels(chinfo(n))

! 5b. allocate the arrays
! -----------------------
       ALLOCATE( atm( self%n_Profiles ),&
            rts( self%n_channels, self%n_profiles ), &
            self%atm_k( self%n_channels, self%n_profiles ), &
            rts_k( self%n_channels, self%n_profiles ), &
            stat = alloc_stat )
       IF ( alloc_stat /= 0 ) THEN
          message = 'error allocating structure arrays'
          CALL display_message( program_name, message, failure )
          STOP
       END IF

! Create the input FORWARD structure (atm)
       CALL CRTM_atmosphere_create( atm, self%n_layers, self%rc%n_absorbers, &
            &self%rc%n_clouds,self%rc%n_aerosols )
       IF ( ANY(.NOT. CRTM_atmosphere_associated(atm)) ) THEN
          message = 'error allocating CRTM forward atmosphere structure'
          CALL display_message( program_name, message, failure )
          STOP
       END IF

! Create the output k-matrix structure
       CALL CRTM_atmosphere_create( self%atm_k, self%n_layers, self%rc%n_absorbers, &
            &self%rc%n_clouds, self%rc%n_aerosols )
       IF ( ANY(.NOT. CRTM_atmosphere_associated(self%atm_k)) ) THEN
          message = 'error allocating CRTM k-matrix atmosphere structure'
          CALL display_message( program_name, message, failure )
          STOP
       END IF
       
!       CALL CRTM_rtsolution_create(rts, self%n_layers )
       
       CALL load_atm_data(self%n_profiles,self%n_layers,geovals,atm)
       IF (self%rc%n_aerosols > 0)  &
            &CALL load_aerosol_data(self%n_profiles,self%n_layers,geovals,atm)
       
       CALL CRTM_atmosphere_zero(self%atm_k )
       
       rts_k%radiance               = zero
       rts_k%brightness_temperature = zero

       DO m = 1, self%n_profiles
          DO l = 1, self%n_channels
             rts_k(l,m)%layer_optical_depth = one
          ENDDO
       ENDDO
       
! 8b.1 the k-matrix model for aod
! ----------------------
       err_stat = CRTM_aod_k( atm,    &  ! forward  input
            rts_k                   , &  ! k-matrix input
            chinfo(n:n)             , &  ! input
            rts                     , &  ! forward  output
            self%atm_k        )          ! k-matrix output
       IF ( err_stat /= success ) THEN
          message = 'error calling CRTM k-matrix model for '&
               &//TRIM(self%rc%sensor_id(n))
          CALL display_message( program_name, message, failure )
          STOP
       END IF
       
       call CRTM_atmosphere_destroy(atm)
       call CRTM_rtsolution_destroy(rts_k)
       call CRTM_rtsolution_destroy(rts)

   ! deallocate all arrays
   ! ---------------------
       deallocate(atm, rts, rts_k, stat = alloc_stat)
       IF ( alloc_stat /= 0 ) THEN
          message = 'error deallocating structure arrays (settraj)'
          CALL display_message( program_name, message, failure )
          STOP
       END IF
       
    END DO sensor_loop
       
  END SUBROUTINE ufo_aod_tlad_settraj

! ------------------------------------------------------------------------------

  SUBROUTINE ufo_aod_simobs_tl(self, geovals, hofx, obss)
    IMPLICIT NONE
    class(ufo_aod_tlad), INTENT(in)     :: self
    TYPE(ufo_geovals),   INTENT(in)     :: geovals
    REAL(kind_real),     INTENT(inout)  :: hofx(:)
    TYPE(c_ptr), value,  INTENT(in)     :: obss

    CHARACTER(len=*), PARAMETER :: myname_="ufo_aod_simobs_tl"

    ! nothing here yet

  END SUBROUTINE ufo_aod_simobs_tl

! ------------------------------------------------------------------------------

  SUBROUTINE ufo_aod_simobs_ad(self, geovals, hofx, obss)

    IMPLICIT NONE
    class(ufo_aod_tlad), INTENT(in)    :: self
    TYPE(ufo_geovals),   INTENT(inout) :: geovals
    REAL(kind_real),     INTENT(in)    :: hofx(:)
    TYPE(c_ptr), value,  INTENT(in)    :: obss

    CHARACTER(len=*), PARAMETER :: myname_="ufo_aod_simobs_tl"

    ! nothing here yet

  END SUBROUTINE ufo_aod_simobs_ad

! ------------------------------------------------------------------------------

END MODULE ufo_aod_tlad_mod
