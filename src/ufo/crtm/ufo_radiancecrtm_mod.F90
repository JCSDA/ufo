! (C) Copyright 2017-2018 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

!> Fortran module to handle radiancecrtm observations

module ufo_radiancecrtm_mod

 use fckit_configuration_module, only: fckit_configuration
 use iso_c_binding
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
   integer, allocatable                           :: channels(:)
   type(crtm_conf) :: conf
 contains
   procedure :: setup  => ufo_radiancecrtm_setup
   procedure :: delete => ufo_radiancecrtm_delete
   procedure :: simobs => ufo_radiancecrtm_simobs
 end type ufo_radiancecrtm

 character(len=maxvarlen), dimension(21), parameter :: varin_default = &
                            (/var_ts, var_prs, var_prsi,                                  &
                              var_sfc_wfrac, var_sfc_lfrac, var_sfc_ifrac, var_sfc_sfrac, &
                              var_sfc_wtmp,  var_sfc_ltmp,  var_sfc_itmp,  var_sfc_stmp,  &
                              var_sfc_vegfrac, var_sfc_wspeed, var_sfc_wdir, var_sfc_lai, &
                              var_sfc_soilm, var_sfc_soilt, var_sfc_landtyp,              &
                              var_sfc_vegtyp, var_sfc_soiltyp, var_sfc_sdepth/)

contains

! ------------------------------------------------------------------------------

subroutine ufo_radiancecrtm_setup(self, f_confOpts, f_confOper, channels)

implicit none
class(ufo_radiancecrtm), intent(inout) :: self
type(fckit_configuration),  intent(in)    :: f_confOpts
type(fckit_configuration),  intent(in)    :: f_confOper
integer(c_int),          intent(in)    :: channels(:)  !List of channels to use

integer :: nvars_in, nvars_out
integer :: ind, jspec, ich
character(len=max_string) :: err_msg

 call crtm_conf_setup(self%conf,f_confOpts,f_confOper)
 if ( ufo_vars_getindex(self%conf%Absorbers, var_mixr) < 1 ) then
   write(err_msg,*) 'ufo_radiancecrtm_setup error: H2O must be included in CRTM Absorbers'
   call abor1_ftn(err_msg)
 end if
 if ( ufo_vars_getindex(self%conf%Absorbers, var_oz) < 1 ) then
   write(err_msg,*) 'ufo_radiancecrtm_setup error: O3 must be included in CRTM Absorbers'
   call abor1_ftn(err_msg)
 end if

 ! request from the model all the hardcoded atmospheric & surface variables + 
 ! 1 * n_Absorbers
 ! 2 * n_Clouds (mass content and effective radius)
 nvars_in = size(varin_default) + self%conf%n_Absorbers + 2 * self%conf%n_Clouds
 allocate(self%varin(nvars_in))
 self%varin(1:size(varin_default)) = varin_default
 ind = size(varin_default) + 1
 !Use list of Absorbers and Clouds from conf
 do jspec = 1, self%conf%n_Absorbers
   self%varin(ind) = self%conf%Absorbers(jspec)
   ind = ind + 1
 end do
 do jspec = 1, self%conf%n_Clouds
   self%varin(ind) = self%conf%Clouds(jspec,1)
   ind = ind + 1
   self%varin(ind) = self%conf%Clouds(jspec,2)
   ind = ind + 1
 end do

 ! save channels
 allocate(self%channels(size(channels)))
 self%channels(:) = channels(:)

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
 n_Profiles = geovals%nlocs
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
                       IRwaterCoeff_File=trim(self%conf%IRwaterCoeff_File), &
                       IRlandCoeff_File=trim(self%conf%IRlandCoeff_File), &
                       IRsnowCoeff_File=trim(self%conf%IRsnowCoeff_File), &
                       IRiceCoeff_File=trim(self%conf%IRiceCoeff_File), &
                       VISwaterCoeff_File=trim(self%conf%VISwaterCoeff_File), &
                       VISlandCoeff_File=trim(self%conf%VISlandCoeff_File), &
                       VISsnowCoeff_File=trim(self%conf%VISsnowCoeff_File), &
                       VISiceCoeff_File=trim(self%conf%VISiceCoeff_File), &
                       MWwaterCoeff_File=trim(self%conf%MWwaterCoeff_File), &
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
   err_stat = CRTM_ChannelInfo_Subset(chinfo(n), self%channels, reset=.false.)
   if ( err_stat /= SUCCESS ) THEN
      message = 'Error subsetting channels!'
      call Display_Message( PROGRAM_NAME, message, FAILURE )
      stop
   end if


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
   call Load_Sfc_Data(n_Profiles,n_Layers,n_Channels,self%channels,geovals,sfc,chinfo,obss,self%conf)
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
       hofx(l,m) = rts(l,m)%Brightness_Temperature
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
