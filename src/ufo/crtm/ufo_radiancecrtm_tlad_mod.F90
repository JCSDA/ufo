! (C) Copyright 2017-2018 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

!> Fortran module to handle tl/ad for radiancecrtm observations

module ufo_radiancecrtm_tlad_mod

 use iso_c_binding
 use config_mod
 use kinds

 use ufo_geovals_mod, only: ufo_geovals, ufo_geoval, ufo_geovals_get_var
 use ufo_vars_mod
 use ufo_crtm_utils_mod
 use crtm_module
 use obsspace_mod
 use missing_values_mod

 implicit none
 private

 !> Fortran derived type for radiancecrtm trajectory
 type, public :: ufo_radiancecrtm_tlad
 private
  character(len=max_string), public, allocatable :: varin(:)  ! variables requested from the model
  integer, allocatable                           :: channels(:)
  type(crtm_conf) :: conf
  integer :: n_Profiles
  integer :: n_Layers
  integer :: n_Channels
  type(CRTM_Atmosphere_type), allocatable :: atm_K(:,:)
  type(CRTM_Surface_type), allocatable :: sfc_K(:,:)
  logical :: ltraj
 contains
  procedure :: setup  => ufo_radiancecrtm_tlad_setup
  procedure :: delete  => ufo_radiancecrtm_tlad_delete
  procedure :: settraj => ufo_radiancecrtm_tlad_settraj
  procedure :: simobs_tl  => ufo_radiancecrtm_simobs_tl
  procedure :: simobs_ad  => ufo_radiancecrtm_simobs_ad
 end type ufo_radiancecrtm_tlad

contains

! ------------------------------------------------------------------------------

subroutine ufo_radiancecrtm_tlad_setup(self, c_conf, channels)

implicit none
class(ufo_radiancecrtm_tlad), intent(inout) :: self
type(c_ptr),                  intent(in)    :: c_conf
integer(c_int),               intent(in)    :: channels(:)  !List of channels to use

integer :: nvars_in

 call crtm_conf_setup(self%conf,c_conf)

 ! only air_temperature adjusted for now
 nvars_in = 1
 allocate(self%varin(nvars_in))
 self%varin(1) = var_ts

 ! save channels
 allocate(self%channels(size(channels)))
 self%channels(:) = channels(:)

end subroutine ufo_radiancecrtm_tlad_setup

! ------------------------------------------------------------------------------

subroutine ufo_radiancecrtm_tlad_delete(self)

implicit none
class(ufo_radiancecrtm_tlad), intent(inout) :: self

 self%ltraj = .false.

 call crtm_conf_delete(self%conf)

 if (allocated(self%atm_k)) then
   call CRTM_Atmosphere_Destroy(self%atm_K)
   deallocate(self%atm_k)
 endif

 if (allocated(self%sfc_k)) then
   call CRTM_Surface_Destroy(self%sfc_K)
   deallocate(self%sfc_k)
 endif

end subroutine ufo_radiancecrtm_tlad_delete

! ------------------------------------------------------------------------------

subroutine ufo_radiancecrtm_tlad_settraj(self, geovals, obss)

implicit none

class(ufo_radiancecrtm_tlad), intent(inout) :: self
type(ufo_geovals),        intent(in)    :: geovals
type(c_ptr), value,       intent(in)    :: obss

! Local Variables
character(*), parameter :: PROGRAM_NAME = 'ufo_radiancecrtm_mod.F90'
character(255) :: message, version
integer        :: err_stat, alloc_stat
integer        :: n, k1
type(ufo_geoval), pointer :: temp

! Define the "non-demoninational" arguments
type(CRTM_ChannelInfo_type)             :: chinfo(self%conf%n_Sensors)
type(CRTM_Geometry_type),   allocatable :: geo(:)

! Define the FORWARD variables
type(CRTM_Atmosphere_type), allocatable :: atm(:)
type(CRTM_Surface_type),    allocatable :: sfc(:)
type(CRTM_RTSolution_type), allocatable :: rts(:,:)

! Define the K-MATRIX variables
type(CRTM_RTSolution_type), allocatable :: rts_K(:,:)


 ! Get number of profile and layers from geovals
 ! ---------------------------------------------
 self%n_Profiles = geovals%nlocs
 call ufo_geovals_get_var(geovals, var_ts, temp)
 self%n_Layers = temp%nval
 nullify(temp)

 ! Program header
 ! --------------
 call CRTM_Version( Version )
 call Program_Message( PROGRAM_NAME, &
                       'Check/example program for the CRTM Forward and K-Matrix (setTraj) functions using '//&
                       trim(self%conf%ENDIAN_type)//' coefficient datafiles', &
                       'CRTM Version: '//TRIM(Version) )


 ! Initialise all the sensors at once
 ! ----------------------------------
 !** NOTE: CRTM_Init points to the various binary files needed for CRTM.  See the
 !**       CRTM_Lifecycle.f90 for more details.

 write( *,'(/5x,"Initializing the CRTM (setTraj) ...")' )
 err_stat = CRTM_Init( self%conf%SENSOR_ID, &
            chinfo, &
            File_Path=trim(self%conf%COEFFICIENT_PATH), &
            Quiet=.TRUE.)
 if ( err_stat /= SUCCESS ) THEN
   message = 'Error initializing CRTM (setTraj)'
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
   self%N_Channels = CRTM_ChannelInfo_n_Channels(chinfo(n))


   ! Allocate the ARRAYS
   ! -------------------
   allocate( geo( self%n_Profiles )                         , &
             atm( self%n_Profiles )                         , &
             sfc( self%n_Profiles )                         , &
             rts( self%N_Channels, self%n_Profiles )        , &
             self%atm_K( self%N_Channels, self%n_Profiles ) , &
             self%sfc_K( self%N_Channels, self%n_Profiles ) , &
             rts_K( self%N_Channels, self%n_Profiles )      , &
             STAT = alloc_stat                                )
   if ( alloc_stat /= 0 ) THEN
      message = 'Error allocating structure arrays (setTraj)'
      call Display_Message( PROGRAM_NAME, message, FAILURE )
      stop
   end if


   ! Create the input FORWARD structure (atm)
   ! ----------------------------------------
   call CRTM_Atmosphere_Create( atm, self%n_Layers, self%conf%n_Absorbers, self%conf%n_Clouds, self%conf%n_Aerosols )
   if ( ANY(.NOT. CRTM_Atmosphere_Associated(atm)) ) THEN
      message = 'Error allocating CRTM Forward Atmosphere structure (setTraj)'
      CALL Display_Message( PROGRAM_NAME, message, FAILURE )
      STOP
   END IF


   ! Create the input FORWARD structure (sfc)
   ! ----------------------------------------
   call CRTM_Surface_Create(sfc, self%N_Channels)
   IF ( ANY(.NOT. CRTM_Surface_Associated(sfc)) ) THEN
      message = 'Error allocating CRTM Surface structure (setTraj)'
      CALL Display_Message( PROGRAM_NAME, message, FAILURE )
      STOP
   END IF


   ! Create output K-MATRIX structure (atm)
   ! --------------------------------------
   call CRTM_Atmosphere_Create( self%atm_K, self%n_Layers, self%conf%n_Absorbers, self%conf%n_Clouds, self%conf%n_Aerosols )
   if ( ANY(.NOT. CRTM_Atmosphere_Associated(self%atm_K)) ) THEN
      message = 'Error allocating CRTM K-matrix Atmosphere structure (setTraj)'
      CALL Display_Message( PROGRAM_NAME, message, FAILURE )
      STOP
   END IF


   ! Create output K-MATRIX structure (sfc)
   ! --------------------------------------
   call CRTM_Surface_Create(self%sfc_K, self%N_Channels)
   IF ( ANY(.NOT. CRTM_Surface_Associated(self%sfc_K)) ) THEN
      message = 'Error allocating CRTM K-matrix Surface structure (setTraj)'
      CALL Display_Message( PROGRAM_NAME, message, FAILURE )
      STOP
   END IF


   !Assign the data from the GeoVaLs
   !--------------------------------
   call Load_Atm_Data(self%N_PROFILES,self%N_LAYERS,geovals,atm,self%conf)
   call Load_Sfc_Data(self%N_PROFILES,self%N_LAYERS,self%N_Channels,self%channels,geovals,sfc,chinfo,obss)
   call Load_Geom_Data(obss,geo)


   ! Zero the K-matrix OUTPUT structures
   ! -----------------------------------
   call CRTM_Atmosphere_Zero( self%atm_K )
   call CRTM_Surface_Zero( self%sfc_K )


   ! Inintialize the K-matrix INPUT so that the results are dTb/dx
   ! -------------------------------------------------------------
   rts_K%Radiance               = ZERO
   rts_K%Brightness_Temperature = ONE


   ! Call the K-matrix model
   ! -----------------------
   err_stat = CRTM_K_Matrix( atm         , &  ! FORWARD  Input
                             sfc         , &  ! FORWARD  Input
                             rts_K       , &  ! K-MATRIX Input
                             geo         , &  ! Input
                             chinfo(n:n) , &  ! Input
                             self%atm_K  , &  ! K-MATRIX Output
                             self%sfc_K  , &  ! K-MATRIX Output
                             rts           )  ! FORWARD  Output
   if ( err_stat /= SUCCESS ) THEN
      message = 'Error calling CRTM (setTraj) K-Matrix Model for '//TRIM(self%conf%SENSOR_ID(n))
      call Display_Message( PROGRAM_NAME, message, FAILURE )
      stop
   end if

   !call CRTM_RTSolution_Inspect(rts)


   ! Deallocate the structures
   ! -------------------------
   call CRTM_Geometry_Destroy(geo)
   call CRTM_Atmosphere_Destroy(atm)
   call CRTM_RTSolution_Destroy(rts_K)
   call CRTM_RTSolution_Destroy(rts)
   call CRTM_Surface_Destroy(sfc)


   ! Deallocate all arrays
   ! ---------------------
   deallocate(geo, atm, sfc, rts, rts_K, STAT = alloc_stat)
   if ( alloc_stat /= 0 ) THEN
      message = 'Error deallocating structure arrays (setTraj)'
      call Display_Message( PROGRAM_NAME, message, FAILURE )
      stop
   end if

 end do Sensor_Loop


 ! Destroy CRTM instance
 ! ---------------------
 write( *, '( /5x, "Destroying the CRTM (setTraj)..." )' )
 err_stat = CRTM_Destroy( chinfo )
 if ( err_stat /= SUCCESS ) THEN
    message = 'Error destroying CRTM (setTraj)'
    call Display_Message( PROGRAM_NAME, message, FAILURE )
    stop
 end if


 ! Set flag that the tracectory was set
 ! ------------------------------------
 self%ltraj = .true.

end subroutine ufo_radiancecrtm_tlad_settraj

! ------------------------------------------------------------------------------

subroutine ufo_radiancecrtm_simobs_tl(self, geovals, obss, nvars, nlocs, hofx)

implicit none
class(ufo_radiancecrtm_tlad), intent(in)    :: self
type(ufo_geovals),        intent(in)    :: geovals
type(c_ptr), value,       intent(in)    :: obss
integer,                  intent(in)    :: nvars, nlocs
real(c_double),           intent(inout) :: hofx(nvars, nlocs)

character(len=*), parameter :: myname_="ufo_radiancecrtm_simobs_tl"
character(max_string) :: err_msg
integer :: jprofile, jchannel, jlevel
type(ufo_geoval), pointer :: tv_d


 ! Initial checks
 ! --------------

 ! Check if trajectory was set
 if (.not. self%ltraj) then
   write(err_msg,*) myname_, ' trajectory wasnt set!'
   call abor1_ftn(err_msg)
 endif

 ! Check if nlocs is consistent in geovals & hofx
 if (geovals%nlocs /= self%n_Profiles) then
   write(err_msg,*) myname_, ' error: nlocs inconsistent!'
   call abor1_ftn(err_msg)
 endif

 ! Initialize hofx
 ! ---------------
 hofx(:,:) = 0.0_kind_real

 ! Temperature
 ! -----------

 ! Get t from geovals
 call ufo_geovals_get_var(geovals, var_ts, tv_d)

 ! Check model levels is consistent in geovals & crtm
 if (tv_d%nval /= self%n_Layers) then
   write(err_msg,*) myname_, ' error: layers inconsistent!'
   call abor1_ftn(err_msg)
 endif

 ! Multiply by Jacobian and add to hofx
 do jprofile = 1, self%n_Profiles
   do jchannel = 1, size(self%channels)
     do jlevel = 1, tv_d%nval
       hofx(jchannel, jprofile) = hofx(jchannel, jprofile) + &
                    self%atm_K(self%channels(jchannel),jprofile)%Temperature(jlevel) * &
                    tv_d%vals(jlevel,jprofile)
     enddo
   enddo
 enddo


end subroutine ufo_radiancecrtm_simobs_tl


! ------------------------------------------------------------------------------

subroutine ufo_radiancecrtm_simobs_ad(self, geovals, obss, nvars, nlocs, hofx)

implicit none
class(ufo_radiancecrtm_tlad), intent(in)    :: self
type(ufo_geovals),        intent(inout) :: geovals
type(c_ptr), value,       intent(in)    :: obss
integer,                  intent(in)    :: nvars, nlocs
real(c_double),           intent(in)    :: hofx(nvars, nlocs)

character(len=*), parameter :: myname_="ufo_radiancecrtm_simobs_ad"
character(max_string) :: err_msg
integer :: jprofile, jchannel, jlevel
type(ufo_geoval), pointer :: tv_d
real(c_double) :: missing


 ! Initial checks
 ! --------------

 ! Check if trajectory was set
 if (.not. self%ltraj) then
   write(err_msg,*) myname_, ' trajectory wasnt set!'
   call abor1_ftn(err_msg)
 endif

 ! Check if nlocs is consistent in geovals & hofx
 if (geovals%nlocs /= self%n_Profiles) then
   write(err_msg,*) myname_, ' error: nlocs inconsistent!'
   call abor1_ftn(err_msg)
 endif

 ! Set missing value
 missing = missing_value(missing)

 ! Temperature
 ! -----------

 ! Get t from geovals
 call ufo_geovals_get_var(geovals, var_ts, tv_d)

 ! allocate if not yet allocated
 if (.not. allocated(tv_d%vals)) then
    tv_d%nlocs = self%n_Profiles
    tv_d%nval = self%n_Layers
    allocate(tv_d%vals(tv_d%nval,tv_d%nlocs))
    tv_d%vals = 0.0_kind_real
 endif


 ! Multiply by Jacobian and add to hofx (adjoint)
 do jprofile = 1, self%n_Profiles
   do jchannel = 1, size(self%channels)
     do jlevel = 1, tv_d%nval
       if (hofx(jchannel, jprofile) /= missing) then
         tv_d%vals(jlevel,jprofile) = tv_d%vals(jlevel,jprofile) + &
                                      self%atm_K(self%channels(jchannel),jprofile)%Temperature(jlevel) * &
                                      hofx(jchannel, jprofile)
       endif
     enddo
   enddo
 enddo


 ! Once all geovals set replace flag
 ! ---------------------------------
 if (.not. geovals%linit ) geovals%linit=.true.


end subroutine ufo_radiancecrtm_simobs_ad

! ------------------------------------------------------------------------------

end module ufo_radiancecrtm_tlad_mod
