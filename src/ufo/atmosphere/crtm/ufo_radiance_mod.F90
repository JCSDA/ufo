! (C) Copyright 2017-2018 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

!> Fortran module to handle radiance observations

module ufo_radiance_mod

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

 !> Fortran derived type for radiance trajectory
 type, public :: ufo_radiance
 private
  type(crtm_conf) :: rc
 contains
   procedure :: setup  => ufo_radiance_setup
   procedure :: delete => ufo_radiance_delete
   procedure :: simobs => ufo_radiance_simobs
 end type ufo_radiance

 CHARACTER(MAXVARLEN), PARAMETER :: varname_tmplate="brightness_temperature"

contains

! ------------------------------------------------------------------------------

subroutine ufo_radiance_setup(self, c_conf)

implicit none
class(ufo_radiance), intent(inout) :: self
type(c_ptr),         intent(in)    :: c_conf

 call crtm_conf_setup(self%rc,c_conf)

end subroutine ufo_radiance_setup

! ------------------------------------------------------------------------------

subroutine ufo_radiance_delete(self)

implicit none
class(ufo_radiance), intent(inout) :: self

 call crtm_conf_delete(self%rc)

end subroutine ufo_radiance_delete

! ------------------------------------------------------------------------------

subroutine ufo_radiance_simobs(self, geovals, hofx, obss, channels)

implicit none
class(ufo_radiance),      intent(in) :: self         !Radiance object
type(ufo_geovals),        intent(in) :: geovals      !Inputs from the model
real(c_double),        intent(inout) :: hofx(:)      !h(x) to return
type(c_ptr), value,       intent(in) :: obss         !ObsSpace
integer(c_int),           intent(in) :: channels(:)  !List of channels to use

! Local Variables
character(*), parameter :: PROGRAM_NAME = 'ufo_radiance_mod.F90'
character(255) :: message, version
integer        :: err_stat, alloc_stat
integer        :: l, m, n, i, s
type(ufo_geoval), pointer :: temp

integer :: n_Profiles
integer :: n_Layers
integer :: n_Channels

! Define the "non-demoninational" arguments
type(CRTM_ChannelInfo_type)             :: chinfo(self%rc%n_Sensors)
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
                       trim(self%rc%ENDIAN_type)//' coefficient datafiles', &
                       'CRTM Version: '//TRIM(Version) )


 ! Initialise all the sensors at once
 ! ----------------------------------
 !** NOTE: CRTM_Init points to the various binary files needed for CRTM.  See the
 !**       CRTM_Lifecycle.f90 for more details.

 write( *,'(/5x,"Initializing the CRTM...")' )
 err_stat = CRTM_Init( self%rc%SENSOR_ID, chinfo, &
                       File_Path=trim(self%rc%COEFFICIENT_PATH), &
                       Quiet=.TRUE.)
 if ( err_stat /= SUCCESS ) THEN
   message = 'Error initializing CRTM'
   call Display_Message( PROGRAM_NAME, message, FAILURE )
   stop
 end if


 ! Loop over all sensors. Not necessary if we're calling CRTM for each sensor
 ! ----------------------------------------------------------------------------
 Sensor_Loop:do n = 1, self%rc%n_Sensors


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
   call CRTM_Atmosphere_Create( atm, n_Layers, self%rc%n_Absorbers, self%rc%n_Clouds, self%rc%n_Aerosols )
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
   call Load_Atm_Data(n_Profiles,n_Layers,geovals,atm,self%rc)
   call Load_Sfc_Data(n_Profiles,n_Layers,n_Channels,channels,geovals,sfc,chinfo,obss)
   call Load_Geom_Data(obss,geo)


   ! Call THE CRTM inspection
   ! ------------------------
   if (self%rc%inspect > 0) then
     call CRTM_Atmosphere_Inspect(atm(self%rc%inspect))
     call CRTM_Surface_Inspect(sfc(self%rc%inspect))
     call CRTM_Geometry_Inspect(geo(self%rc%inspect))
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
      message = 'Error calling CRTM Forward Model for '//TRIM(self%rc%SENSOR_ID(n))
      call Display_Message( PROGRAM_NAME, message, FAILURE )
      stop
   end if

   !call CRTM_RTSolution_Inspect(rts)

   ! Put simulated brightness temperature into hofx
   ! ----------------------------------------------

   !Set to zero and initialize counter
   hofx(:) = 0.0_kind_real
   i = 1

   do m = 1, n_Profiles
     do l = 1, size(channels)

       hofx(i) = rts(channels(l),m)%Brightness_Temperature
       i = i + 1

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

end subroutine ufo_radiance_simobs

! ------------------------------------------------------------------------------

end module ufo_radiance_mod
