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
  character(len=MAXVARLEN), public, allocatable :: varin(:), varin_aero(:) ! variablesrequested from the model
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

 CALL assign_aerosol_names(self%conf%aerosol_option, var_aerosols)

 nvars_in = SIZE(varin_default)+SIZE(var_aerosols)
 allocate(self%varin(nvars_in))
 self%varin(1:size(varin_default)) = varin_default
 self%varin(SIZE(varin_default)+1:) = var_aerosols

 
 allocate(self%varin_aero(SIZE(var_aerosols)))
 self%varin_aero(:) = var_aerosols(:)

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
class(ufo_aodcrtm),       intent(in)    :: self
type(ufo_geovals),        intent(in)    :: geovals
integer,                  intent(in)    :: nvars, nlocs
real(c_double),           intent(inout) :: hofx(nvars, nlocs)
type(c_ptr), value,       intent(in)    :: obss

! Local Variables
character(*), parameter :: PROGRAM_NAME = 'ufo_aodcrtm_mod.F90'
character(len=MAXVARLEN) :: def_aero_mod
character(255) :: message, version
integer        :: err_stat, alloc_stat
integer        :: l, m, n, i
type(ufo_geoval), pointer :: temp
real(c_double) :: missing

integer :: n_Profiles
integer :: n_Layers
integer :: n_Channels

logical :: jacobian_needed

! Define the Channel Info and Geometry  arguments
type(CRTM_ChannelInfo_type)             :: chinfo(self%conf%n_Sensors)
type(CRTM_Geometry_type),   allocatable :: geo(:)

! Define the FORWARD variables
type(CRTM_Atmosphere_type), allocatable :: atm(:)
type(CRTM_Surface_type),    allocatable :: sfc(:)
type(CRTM_RTSolution_type), allocatable :: rts(:,:)

! Define the K-MATRIX variables - necessary for AOD call
! ---------------------------------
type(CRTM_Atmosphere_type), allocatable :: atm_K(:,:)
type(CRTM_RTSolution_type), allocatable :: rts_K(:,:)

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
 call define_aerosol_model(self%conf%AerosolCoeff_File, def_aero_mod)
 err_stat = CRTM_Init( self%conf%SENSOR_ID, &
            chinfo, &
            File_Path           = trim(self%conf%COEFFICIENT_PATH), &
            NC_File_Path        = trim(self%conf%NC_COEFFICIENT_PATH), &
            Aerosol_Model       = trim(def_aero_mod), &
            AerosolCoeff_Format = trim(self%conf%AerosolCoeff_Format), &
            AerosolCoeff_File   = trim(self%conf%AerosolCoeff_File), &
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
   if ( any(.not. CRTM_Atmosphere_Associated(atm)) ) then
      message = 'Error allocating CRTM Forward Atmosphere structure'
      call Display_Message( PROGRAM_NAME, message, FAILURE )
      stop
   end if

   call CRTM_RTSolution_Create(rts, n_Layers )


   !Assign the data from the GeoVaLs
   !--------------------------------
   call Load_Atm_Data(n_Profiles,n_Layers,geovals,atm,self%conf)

   if (trim(self%conf%aerosol_option) /= "") &
       & call load_aerosol_data(n_profiles, n_layers, geovals, &
       & self%conf, self%varin_aero, trim(def_aero_mod), atm)

   ! Call THE CRTM inspection
   ! ------------------------
   if (self%conf%inspect > 0) then
     call CRTM_Atmosphere_Inspect(atm(self%conf%inspect))
     call CRTM_ChannelInfo_Inspect(chinfo(n))
   end if

   ! Start processing the CRTM AOD 
   ! |-> Jacobian K-Matrix .OR.
   ! |-> Forward Operator
   ! ------------------------------------------------
   jacobian_needed = .false.
   kmatrix : select case ( jacobian_needed )
     
     case(.true.)
       !
       ! Description:
       ! ===========       
       ! For the time being this branch is NEVER executed because
       ! jacobian_needed = .false. . The Jacobian may be needed in the future
       ! for QC and bias correction but running the Forward operator alone for
       ! ufo simobs requires less resources. The CRTM AOD Jacobian code will
       ! remain unexecuted until the necessary configuration code is added.
       !

       ! The output K-MATRIX structure:
 
       ! Allocate the k-matrix arrays     
       allocate( atm_K( n_channels, N_PROFILES ), &
                 rts_K( n_channels, N_PROFILES ), &
                 STAT = alloc_stat )
       if ( alloc_stat /= 0 ) then
         message = 'Error allocating structure arrays'
         call Display_Message( PROGRAM_NAME, message, FAILURE )
         stop
       end if

       ! Call the constructor for the k-matrix output structure
       call CRTM_Atmosphere_Create( atm_K, n_layers, self%conf%n_Absorbers, self%conf%n_Clouds, self%conf%n_Aerosols)
       if ( any(.not. CRTM_Atmosphere_Associated(atm_K)) ) then
         message = 'Error allocating CRTM K-matrix Atmosphere structure'
         call Display_Message( PROGRAM_NAME, message, FAILURE )
         stop
       end if
   
       ! Call the constructor for the k-matrix input structure
       call CRTM_RTSolution_Create(rts_k, n_Layers )

       ! ============================================================================
       !    **** INITIALIZE THE K-MATRIX ARGUMENTS ****
       !
       !     Zero the K-matrix OUTPUT structures
       ! ---------------------------------------
       call CRTM_Atmosphere_Zero( atm_k )

       !     Initialize the K-matrix INPUT
       ! ----------------------------------
       do m = 1, n_profiles
         do l = 1, SIZE(self%channels)
           rts_K(self%channels(l),m)%Layer_Optical_Depth = ONE
         end do
       end do
       ! ============================================================================


       !     The K-matrix model y = K.x for AOD 
       ! -----------------------------------------
       err_stat = CRTM_AOD_K( atm,    &  ! FORWARD  Input
            rts_K                   , &  ! K-MATRIX Input
            chinfo(n:n)             , &  ! Input
            rts                     , &  ! FORWARD  Output
            atm_k        )               ! K-MATRIX Output

       if ( err_stat /= SUCCESS ) then
         message = 'Error calling CRTM K-Matrix Model for '//TRIM(self%conf%SENSOR_ID(n))
         call Display_Message( PROGRAM_NAME, message, FAILURE )
         stop
       end if

        
       call CRTM_Atmosphere_Destroy(atm_k)
       call CRTM_RTSolution_Destroy(rts_k)
        
       deallocate( atm_K, rts_K )
              
       if ( alloc_stat /= 0 ) THEN
         message = 'Error deallocating Jacobian structure arrays'
         call Display_Message( PROGRAM_NAME, message, FAILURE )
         stop
       end if

     case default
 
       !     The Forward Operator y = H(x) for AOD
       ! ------------------------------------------
       err_stat = CRTM_AOD( atm          , &  ! FORWARD  Input
                            chinfo(n:n)  , &  ! Input
                            rts          )    ! FORWARD  Output

       if ( err_stat /= SUCCESS ) then
         message = 'Error calling CRTM Forward Model for '//TRIM(self%conf%SENSOR_ID(n))
         call Display_Message( PROGRAM_NAME, message, FAILURE )
         stop
       end if

   end select kmatrix


   ! Put simulated brightness temperature into hofx
   ! ----------------------------------------------

   ! Set missing value
   missing = missing_value(missing)
   hofx = missing

   do m = 1, n_profiles
       do l = 1, size(self%channels)
         hofx(l,m) = sum(rts(self%channels(l),m)%Layer_Optical_Depth)
      end do
   end do

   ! Deallocate the structures
   ! -------------------------
   call CRTM_Atmosphere_Destroy(atm)
   call CRTM_RTSolution_Destroy(rts)

   ! Deallocate all arrays
   ! ---------------------
   deallocate(geo, atm, sfc, rts, STAT = alloc_stat)
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
 if ( err_stat /= SUCCESS ) then
    message = 'Error destroying CRTM'
    call Display_Message( PROGRAM_NAME, message, FAILURE )
    stop
 end if

end subroutine ufo_aodcrtm_simobs

! ------------------------------------------------------------------------------

end module ufo_aodcrtm_mod
