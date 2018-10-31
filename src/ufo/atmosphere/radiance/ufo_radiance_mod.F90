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
 use ufo_basis_mod, only: ufo_basis
 use ufo_vars_mod
 use ufo_radiance_utils_mod
 use crtm_module
 use obsspace_mod

 implicit none
 private

 !> Fortran derived type for radiance trajectory
 type, extends(ufo_basis), public :: ufo_radiance
 private
  type(rad_conf) :: rc
 contains
   procedure :: setup  => ufo_radiance_setup
   procedure :: delete => ufo_radiance_delete
   procedure :: simobs => ufo_radiance_simobs
 end type ufo_radiance

contains

! ------------------------------------------------------------------------------

subroutine ufo_radiance_setup(self, c_conf)

implicit none
class(ufo_radiance), intent(inout) :: self
type(c_ptr),         intent(in)    :: c_conf

 call rad_conf_setup(self%rc,c_conf)

end subroutine ufo_radiance_setup

! ------------------------------------------------------------------------------

subroutine ufo_radiance_delete(self)

implicit none
class(ufo_radiance), intent(inout) :: self

 call rad_conf_delete(self%rc)

end subroutine ufo_radiance_delete

! ------------------------------------------------------------------------------

subroutine ufo_radiance_simobs(self, geovals, hofx, obss)

implicit none
class(ufo_radiance),      intent(in) :: self
type(ufo_geovals),        intent(in) :: geovals
real(c_double),        intent(inout) :: hofx(:)
type(c_ptr), value,       intent(in) :: obss

! Local Variables
character(*), parameter :: PROGRAM_NAME = 'ufo_radiance_mod.F90'
character(255) :: message, version
integer        :: err_stat, alloc_stat
integer        :: l, m, n, i, s, ierr
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
 call ufo_geovals_get_var(geovals, var_tv, temp, status=ierr)
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
 err_stat = CRTM_Init( self%rc%SENSOR_ID, &
            chinfo, &
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
   call Load_Atm_Data(n_Profiles,n_Layers,geovals,atm)
   call Load_Sfc_Data(n_Profiles,n_Layers,n_Channels,geovals,sfc,chinfo,obss)
   call Load_Geom_Data(obss,geo)


   ! Call THE CRTM inspection
   ! ------------------------
   call CRTM_Atmosphere_Inspect(atm(12))
   call CRTM_Surface_Inspect(sfc(12))
   call CRTM_Geometry_Inspect(geo(12))
   call CRTM_ChannelInfo_Inspect(chinfo(1))


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


   ! Put simulated brightness temperature into hofx
   ! ----------------------------------------------

   !Set to zero and initializ counter
   hofx(:) = 0.0_kind_real
   i = 1

   do m = 1, n_Profiles
     do l = 1, N_Channels

       hofx(i) = rts(l,m)%Brightness_Temperature
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
