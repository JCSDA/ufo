! (C) Copyright 2017-2018 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

!> Fortran module to handle radiancecrtm observations

module ufo_radiancecrtm_mod

 use iso_c_binding
 use config_mod
 use kinds

 use ufo_geovals_mod, only: ufo_geovals, ufo_geoval, ufo_geovals_get_var
 use ufo_vars_mod
 use ufo_crtm_utils_mod
 use crtm_module
 use obsspace_mod

 implicit none
 private

 !> Fortran derived type for radiancecrtm trajectory
 type, public :: ufo_radiancecrtm
 private
   character(len=max_string), public, allocatable :: varin(:)  ! variables requested from the model
   character(len=max_string), public, allocatable :: varout(:) ! variables simulated by CRTM
   integer, allocatable                           :: channels(:)
   type(crtm_conf) :: conf
 contains
   procedure :: setup  => ufo_radiancecrtm_setup
   procedure :: delete => ufo_radiancecrtm_delete
   procedure :: simobs => ufo_radiancecrtm_simobs
 end type ufo_radiancecrtm

 character(len=maxvarlen), dimension(24), parameter :: varin_default = &
                            (/var_ts, var_mixr, var_prs, var_prsi, var_oz, var_co2,       &
                              var_sfc_wfrac, var_sfc_lfrac, var_sfc_ifrac, var_sfc_sfrac, &
                              var_sfc_wtmp,  var_sfc_ltmp,  var_sfc_itmp,  var_sfc_stmp,  &
                              var_sfc_vegfrac, var_sfc_wspeed, var_sfc_wdir, var_sfc_lai, &
                              var_sfc_soilm, var_sfc_soilt, var_sfc_landtyp,              &
                              var_sfc_vegtyp, var_sfc_soiltyp, var_sfc_sdepth/)
contains

! ------------------------------------------------------------------------------

subroutine ufo_radiancecrtm_setup(self, c_conf, channels)

implicit none
class(ufo_radiancecrtm), intent(inout) :: self
type(c_ptr),             intent(in)    :: c_conf
integer(c_int),          intent(in)    :: channels(:)  !List of channels to use

integer :: nvars_in, nvars_out
integer :: ind, ich

 call crtm_conf_setup(self%conf,c_conf)

 ! request from the model all the hardcoded atmospheric & surface variables + 
 ! 2 * n_clouds (for mass content and effective radius)
 nvars_in = size(varin_default) + self%conf%n_clouds * 2
 allocate(self%varin(nvars_in))
 self%varin(1:size(varin_default)) = varin_default
 ind = size(varin_default) + 1
 if (self%conf%n_clouds > 0) then
   self%varin(ind)   = var_clw
   self%varin(ind+1) = var_clwefr
   ind = ind + 2
 endif
 if (self%conf%n_clouds > 1) then
   self%varin(ind)   = var_cli
   self%varin(ind+1) = var_cliefr
   ind = ind + 2
 endif

 ! output variables: all requested channels
 nvars_out = size(channels)
 allocate(self%varout(nvars_out))
 allocate(self%channels(nvars_out))
 self%channels(:) = channels(:)
 do ich = 1, size(channels)
   call get_var_name(ich, self%varout(ich))
 enddo

end subroutine ufo_radiancecrtm_setup

! ------------------------------------------------------------------------------

subroutine ufo_radiancecrtm_delete(self)

implicit none
class(ufo_radiancecrtm), intent(inout) :: self

 call crtm_conf_delete(self%conf)

end subroutine ufo_radiancecrtm_delete

! ------------------------------------------------------------------------------

subroutine ufo_radiancecrtm_simobs(self, geovals, obss, nvars, nlocs, hofx)

implicit none
class(ufo_radiancecrtm),  intent(in) :: self         !Radiance object
type(ufo_geovals),        intent(in) :: geovals      !Inputs from the model
integer,                  intent(in) :: nvars, nlocs
real(c_double),        intent(inout) :: hofx(nvars, nlocs) !h(x) to return
type(c_ptr), value,       intent(in) :: obss         !ObsSpace

! Local Variables
character(*), parameter :: PROGRAM_NAME = 'ufo_radiancecrtm_mod.F90'
character(255) :: message, version
integer        :: err_stat, alloc_stat
integer        :: l, m, n, s
type(ufo_geoval), pointer :: temp

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


 ! Get number of profile and layers from geovals
 ! ---------------------------------------------
 n_Profiles = geovals%nobs
 call ufo_geovals_get_var(geovals, var_ts, temp)
 n_Layers = temp%nval
 nullify(temp)


 ! Program header
 ! --------------
 call CRTM_Version( Version )
 call Program_Message( PROGRAM_NAME, &
                       'Check/example program for the CRTM Forward and K-Matrix functions using '//&
                       trim(self%conf%ENDIAN_type)//' coefficient datafiles', &
                       'CRTM Version: '//TRIM(Version) )


 ! Initialise all the sensors at once
 ! ----------------------------------
 !** NOTE: CRTM_Init points to the various binary files needed for CRTM.  See the
 !**       CRTM_Lifecycle.f90 for more details.

 write( *,'(/5x,"Initializing the CRTM...")' )
 err_stat = CRTM_Init( self%conf%SENSOR_ID, chinfo, &
                       File_Path=trim(self%conf%COEFFICIENT_PATH), &
                       Quiet=.TRUE.)
 if ( err_stat /= SUCCESS ) THEN
   message = 'Error initializing CRTM'
   call Display_Message( PROGRAM_NAME, message, FAILURE )
   stop
 end if


 ! Loop over all sensors. Not necessary if we're calling CRTM for each sensor
 ! ----------------------------------------------------------------------------
 Sensor_Loop:do n = 1, self%conf%n_Sensors


   ! Pass channel list to CRTM
   ! -------------------------
   !err_stat = CRTM_ChannelInfo_Subset(chinfo(n), channels, reset=.false.)
   !if ( err_stat /= SUCCESS ) THEN
   !   message = 'Error subsetting channels'
   !   call Display_Message( PROGRAM_NAME, message, FAILURE )
   !   stop
   !end if


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


   ! Create the input FORWARD structure (sfc)
   ! ----------------------------------------
   call CRTM_Surface_Create(sfc, N_Channels)
   IF ( ANY(.NOT. CRTM_Surface_Associated(sfc)) ) THEN
      message = 'Error allocating CRTM Surface structure'
      CALL Display_Message( PROGRAM_NAME, message, FAILURE )
      STOP
   END IF


   !Assign the data from the GeoVaLs
   !--------------------------------
   call Load_Atm_Data(n_Profiles,n_Layers,geovals,atm,self%conf)
   call Load_Sfc_Data(n_Profiles,n_Layers,n_Channels,self%channels,geovals,sfc,chinfo,obss)
   call Load_Geom_Data(obss,geo)


   ! Call THE CRTM inspection
   ! ------------------------
   if (self%conf%inspect > 0) then
     call CRTM_Atmosphere_Inspect(atm(self%conf%inspect))
     call CRTM_Surface_Inspect(sfc(self%conf%inspect))
     call CRTM_Geometry_Inspect(geo(self%conf%inspect))
     call CRTM_ChannelInfo_Inspect(chinfo(n))
   endif

   ! Call the forward model call for each sensor
   ! -------------------------------------------
   err_stat = CRTM_Forward( atm        , &  ! Input
                            sfc        , &  ! Input
                            geo        , &  ! Input
                            chinfo(n:n), &  ! Input
                            rts          )  ! Output
   if ( err_stat /= SUCCESS ) THEN
      message = 'Error calling CRTM Forward Model for '//TRIM(self%conf%SENSOR_ID(n))
      call Display_Message( PROGRAM_NAME, message, FAILURE )
      stop
   end if

   !call CRTM_RTSolution_Inspect(rts)

   ! Put simulated brightness temperature into hofx
   ! ----------------------------------------------

   !Set to zero and initialize counter
   hofx(:,:) = 0.0_kind_real

   do m = 1, n_Profiles
     do l = 1, size(self%channels)

       hofx(l, m) = rts(self%channels(l),m)%Brightness_Temperature

     end do
   end do


   ! Deallocate the structures
   ! -------------------------
   call CRTM_Geometry_Destroy(geo)
   call CRTM_Atmosphere_Destroy(atm)
   call CRTM_RTSolution_Destroy(rts)
   call CRTM_Surface_Destroy(sfc)


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
 write( *, '( /5x, "Destroying the CRTM..." )' )
 err_stat = CRTM_Destroy( chinfo )
 if ( err_stat /= SUCCESS ) THEN
    message = 'Error destroying CRTM'
    call Display_Message( PROGRAM_NAME, message, FAILURE )
    stop
 end if

end subroutine ufo_radiancecrtm_simobs

! ------------------------------------------------------------------------------

end module ufo_radiancecrtm_mod
