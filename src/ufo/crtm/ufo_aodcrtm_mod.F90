! (C) Copyright 2017-2018 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

!> Fortran module to handle aod observations

module ufo_aodcrtm_mod

 use fckit_configuration_module, only: fckit_configuration
 use iso_c_binding
 use kinds
 use missing_values_mod

 use ufo_geovals_mod, only: ufo_geovals, ufo_geoval, ufo_geovals_get_var
 use ufo_vars_mod
 use ufo_crtm_utils_mod
 use crtm_module
 use obsspace_mod

 implicit none
 private

 !> Fortran derived type for aod trajectory
 type, public :: ufo_aodcrtm
 private
  character(len=MAXVARLEN), public, allocatable :: varin(:)  ! variablesrequested from the model
  integer, allocatable                          :: channels(:)
  type(crtm_conf) :: conf
 contains
   procedure :: setup  => ufo_aodcrtm_setup
   procedure :: delete => ufo_aodcrtm_delete
   procedure :: simobs => ufo_aodcrtm_simobs
 end type ufo_aodcrtm
 CHARACTER(len=maxvarlen), DIMENSION(5), PARAMETER :: varin_default = (/var_ts, var_mixr, var_rh, var_prs, var_prsi/) 

 CHARACTER(MAXVARLEN), PARAMETER :: varname_tmplate="aerosol_optical_depth"

contains

! ------------------------------------------------------------------------------

subroutine ufo_aodcrtm_setup(self, f_confOper, channels)

implicit none
class(ufo_aodcrtm),        intent(inout) :: self
type(fckit_configuration), intent(in)    :: f_confOper
integer(c_int),            intent(in)    :: channels(:)  !List of channels to use

integer :: nvars_in
character(len=max_string) :: err_msg
type(fckit_configuration) :: f_confOpts

CHARACTER(len=MAXVARLEN), ALLOCATABLE :: var_aerosols(:)

 call f_confOper%get_or_die("obs options",f_confOpts)

 call crtm_conf_setup(self%conf, f_confOpts, f_confOper)
 if ( ufo_vars_getindex(self%conf%Absorbers, var_mixr) /= 1 ) then
   write(err_msg,*) 'ufo_aodcrtm_setup error: H2O must be first in CRTM Absorbers for AOD'
   call abor1_ftn(err_msg)
 end if
 if ( ufo_vars_getindex(self%conf%Absorbers, var_oz) < 2 ) then
   write(err_msg,*) 'ufo_aodcrtm_setup error: O3 must be included in CRTM Absorbers'
   call abor1_ftn(err_msg)
 end if

 CALL assign_aerosol_names(self%conf%aerosol_option,var_aerosols)

 nvars_in = SIZE(varin_default)+SIZE(var_aerosols)
 allocate(self%varin(nvars_in))
 self%varin(1:size(varin_default)) = varin_default
 self%varin(SIZE(varin_default)+1:) = var_aerosols

 ! save channels
 allocate(self%channels(size(channels)))
 self%channels(:) = channels(:)

end subroutine ufo_aodcrtm_setup

! ------------------------------------------------------------------------------

subroutine ufo_aodcrtm_delete(self)

implicit none
class(ufo_aodcrtm), intent(inout) :: self

 call crtm_conf_delete(self%conf)

end subroutine ufo_aodcrtm_delete

! ------------------------------------------------------------------------------

SUBROUTINE ufo_aodcrtm_simobs(self, geovals, obss, nvars, nlocs, hofx)

implicit none
class(ufo_aodcrtm),      intent(in) :: self
type(ufo_geovals),        intent(in) :: geovals
integer,                  intent(in) :: nvars, nlocs
real(c_double),        intent(inout) :: hofx(nvars, nlocs)
type(c_ptr), value,       intent(in) :: obss

! Local Variables
character(*), parameter :: PROGRAM_NAME = 'ufo_aodcrtm_mod.F90'
character(255) :: message, version
integer        :: err_stat, alloc_stat
integer        :: l, m, n, i
type(ufo_geoval), pointer :: temp
real(c_double) :: missing

integer :: n_Profiles
integer :: n_Layers
integer :: n_Channels

! Define the "non-demoninational" arguments
type(CRTM_ChannelInfo_type)             :: chinfo(self%conf%n_Sensors)
type(CRTM_Geometry_type),   allocatable :: geo(:)

! Define the FORWARD variables
type(CRTM_Atmosphere_type), allocatable :: atm(:)
type(CRTM_Surface_type),    allocatable :: sfc(:)
type(CRTM_RTSolution_type), allocatable :: rts(:,:)

! Define the K-MATRIX variables - necessary for AOD call
! ---------------------------------
TYPE(CRTM_Atmosphere_type), ALLOCATABLE :: atm_K(:,:)
TYPE(CRTM_RTSolution_type), ALLOCATABLE :: rts_K(:,:)

 ! Get number of profile and layers from geovals
 ! ---------------------------------------------
 n_Profiles = geovals%nlocs
 call ufo_geovals_get_var(geovals, var_ts, temp)
 n_Layers = temp%nval
 nullify(temp)

 ! Program header
 ! --------------
 ! call CRTM_Version( Version )
 ! call Program_Message( PROGRAM_NAME, &
 !                       'Check/example program for the CRTM Forward and K-Matrix functions using '//&
 !                       trim(self%conf%ENDIAN_type)//' coefficient datafiles', &
 !                       'CRTM Version: '//TRIM(Version) )


 ! Initialise all the sensors at once
 ! ----------------------------------
 !** NOTE: CRTM_Init points to the various binary files needed for CRTM.  See the
 !**       CRTM_Lifecycle.f90 for more details.

 ! write( *,'(/5x,"Initializing the CRTM...")' )
 err_stat = CRTM_Init( self%conf%SENSOR_ID, &
            chinfo, &
            File_Path=trim(self%conf%COEFFICIENT_PATH), &
            Quiet=.TRUE.)
 if ( err_stat /= SUCCESS ) THEN
   message = 'Error initializing CRTM'
   call Display_Message( PROGRAM_NAME, message, FAILURE )
   stop
 end if


 ! Loop over all sensors. Not necessary if we're calling CRTM for each sensor
 ! ----------------------------------------------------------------------------
 Sensor_Loop:DO n = 1, self%conf%n_Sensors


   ! Determine the number of channels for the current sensor
   ! -------------------------------------------------------
   N_Channels = CRTM_ChannelInfo_n_Channels(chinfo(n))


   ! Allocate the ARRAYS
   ! -------------------
   allocate( geo( n_Profiles ),               &
             atm( n_Profiles ),               &
             sfc( n_Profiles ),               &
             rts( N_Channels, n_Profiles ),   &
             STAT = alloc_stat )
   if ( alloc_stat /= 0 ) THEN
      message = 'Error allocating structure arrays'
      call Display_Message( PROGRAM_NAME, message, FAILURE )
      stop
   end if


   ! Create the input FORWARD structure (atm)
   ! ----------------------------------------
   call CRTM_Atmosphere_Create( atm, n_Layers, self%conf%n_Absorbers, self%conf%n_Clouds, self%conf%n_Aerosols )
   if ( ANY(.NOT. CRTM_Atmosphere_Associated(atm)) ) THEN
      message = 'Error allocating CRTM Forward Atmosphere structure'
      CALL Display_Message( PROGRAM_NAME, message, FAILURE )
      STOP
   END IF

   ALLOCATE( atm_K( n_channels, N_PROFILES ), &
        rts_K( n_channels, N_PROFILES ), &
        STAT = alloc_stat )
   IF ( alloc_stat /= 0 ) THEN
      message = 'Error allocating structure arrays'
      CALL Display_Message( PROGRAM_NAME, message, FAILURE )
      STOP
   END IF

! The output K-MATRIX structure
   CALL CRTM_Atmosphere_Create( atm_K, n_layers, self%conf%n_Absorbers, self%conf%n_Clouds, self%conf%n_Aerosols)
   IF ( ANY(.NOT. CRTM_Atmosphere_Associated(atm_K)) ) THEN
      message = 'Error allocating CRTM K-matrix Atmosphere structure'
      CALL Display_Message( PROGRAM_NAME, message, FAILURE )
      STOP
   END IF
   
   CALL CRTM_RTSolution_Create(rts, n_Layers )
   CALL CRTM_RTSolution_Create(rts_k, n_Layers )

   !Assign the data from the GeoVaLs
   !--------------------------------
   CALL Load_Atm_Data(n_Profiles,n_Layers,geovals,atm,self%conf)

   IF (TRIM(self%conf%aerosol_option) /= "") &
        &CALL load_aerosol_data(n_profiles,n_layers,geovals,&
        &self%conf%aerosol_option,atm)

   ! Call THE CRTM inspection
   ! ------------------------
   IF (self%conf%inspect > 0) THEN
      CALL CRTM_Atmosphere_Inspect(atm(self%conf%inspect))
      CALL CRTM_ChannelInfo_Inspect(chinfo(n))
   ENDIF

! 8b.1 The K-matrix model for AOD
! ----------------------
   err_stat = CRTM_AOD_K( atm,    &  ! FORWARD  Input
        rts_K                   , &  ! K-MATRIX Input
        chinfo(n:n)             , &  ! Input
        rts                     , &  ! FORWARD  Output
        atm_k        )               ! K-MATRIX Output

   if ( err_stat /= SUCCESS ) THEN
      message = 'Error calling CRTM Forward Model for '//TRIM(self%conf%SENSOR_ID(n))
      call Display_Message( PROGRAM_NAME, message, FAILURE )
      stop
   end if


   ! Put simulated brightness temperature into hofx
   ! ----------------------------------------------

   ! Set missing value
   missing = missing_value(missing)
   hofx = missing

   DO m = 1, n_profiles
       DO l = 1, size(self%channels)
         hofx(l,m) = SUM(rts(self%channels(l),m)%layer_optical_depth)
      END DO
   END DO

   ! Deallocate the structures
   ! -------------------------
   call CRTM_Atmosphere_Destroy(atm)
   call CRTM_RTSolution_Destroy(rts)

   call CRTM_Atmosphere_Destroy(atm_k)
   call CRTM_RTSolution_Destroy(rts_k)

   ! Deallocate all arrays
   ! ---------------------
   DEALLOCATE(geo, atm, sfc, rts, atm_k, rts_k, STAT = alloc_stat)
   if ( alloc_stat /= 0 ) THEN
      message = 'Error deallocating structure arrays'
      call Display_Message( PROGRAM_NAME, message, FAILURE )
      stop
   end if

end do Sensor_Loop

 ! Destroy CRTM instance
 ! ---------------------
 ! write( *, '( /5x, "Destroying the CRTM..." )' )
 err_stat = CRTM_Destroy( chinfo )
 if ( err_stat /= SUCCESS ) THEN
    message = 'Error destroying CRTM'
    call Display_Message( PROGRAM_NAME, message, FAILURE )
    stop
 end if

end subroutine ufo_aodcrtm_simobs

! ------------------------------------------------------------------------------

end module ufo_aodcrtm_mod
