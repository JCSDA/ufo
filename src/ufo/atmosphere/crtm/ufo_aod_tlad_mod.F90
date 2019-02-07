! (C) Copyright 2017-2018 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

!> Fortran module to handle tl/ad for aod observations

module ufo_aod_tlad_mod

 use iso_c_binding
 use config_mod
 use kinds

 use ufo_geovals_mod, only: ufo_geovals, ufo_geoval, ufo_geovals_get_var
 use ufo_basis_tlad_mod, only: ufo_basis_tlad
 use ufo_vars_mod
 use ufo_crtm_utils_mod
 use crtm_module
 use obsspace_mod

 implicit none
 private

 !> Fortran derived type for aod trajectory
 type, extends(ufo_basis_tlad), public :: ufo_aod_tlad
 private
  type(crtm_conf) :: rc
  integer :: n_Profiles
  integer :: n_Layers
  integer :: n_Channels
  type(CRTM_Atmosphere_type), allocatable :: atm_K(:,:)
  type(CRTM_Surface_type), allocatable :: sfc_K(:,:)
 contains
  procedure :: setup  => ufo_aod_tlad_setup
  procedure :: delete  => ufo_aod_tlad_delete
  procedure :: settraj => ufo_aod_tlad_settraj
  procedure :: simobs_tl  => ufo_aod_simobs_tl
  procedure :: simobs_ad  => ufo_aod_simobs_ad
 end type ufo_aod_tlad

contains

! ------------------------------------------------------------------------------

subroutine ufo_aod_tlad_setup(self, c_conf)

implicit none
class(ufo_aod_tlad), intent(inout) :: self
type(c_ptr),              intent(in)    :: c_conf

 call crtm_conf_setup(self%rc,c_conf)

end subroutine ufo_aod_tlad_setup

! ------------------------------------------------------------------------------

subroutine ufo_aod_tlad_delete(self)

implicit none
class(ufo_aod_tlad), intent(inout) :: self

 self%ltraj = .false.
 call crtm_conf_delete(self%rc)

 if (allocated(self%atm_k)) then
   call CRTM_Atmosphere_Destroy(self%atm_K)
   deallocate(self%atm_k)
 endif

 if (allocated(self%sfc_k)) then
   call CRTM_Surface_Destroy(self%sfc_K)
   deallocate(self%sfc_k)
 endif

end subroutine ufo_aod_tlad_delete

! ------------------------------------------------------------------------------

subroutine ufo_aod_tlad_settraj(self, geovals, obss)

implicit none

class(ufo_aod_tlad), intent(inout) :: self
type(ufo_geovals),        intent(in)    :: geovals
type(c_ptr), value,       intent(in)    :: obss

! Local Variables
character(*), parameter :: PROGRAM_NAME = 'ufo_aod_tlad_mod.F90'
character(255) :: message, version
integer        :: err_stat, alloc_stat
INTEGER        :: n, k1,l,m
type(ufo_geoval), pointer :: temp

! Define the "non-demoninational" arguments
type(CRTM_ChannelInfo_type)             :: chinfo(self%rc%n_Sensors)
type(CRTM_Geometry_type),   allocatable :: geo(:)

! Define the FORWARD variables
type(CRTM_Atmosphere_type), allocatable :: atm(:)
type(CRTM_Surface_type),    allocatable :: sfc(:)
type(CRTM_RTSolution_type), allocatable :: rts(:,:)

! Define the K-MATRIX variables
type(CRTM_RTSolution_type), allocatable :: rts_K(:,:)


 ! Get number of profile and layers from geovals
 ! ---------------------------------------------
 self%n_Profiles = geovals%nobs
 call ufo_geovals_get_var(geovals, var_ts, temp)
 self%n_Layers = temp%nval
 nullify(temp)

 ! Program header
 ! --------------
 call CRTM_Version( Version )
 call Program_Message( PROGRAM_NAME, &
                       'Check/example program for the CRTM Forward and K-Matrix (setTraj) functions using '//&
                       trim(self%rc%ENDIAN_type)//' coefficient datafiles', &
                       'CRTM Version: '//TRIM(Version) )


 ! Initialise all the sensors at once
 ! ----------------------------------
 !** NOTE: CRTM_Init points to the various binary files needed for CRTM.  See the
 !**       CRTM_Lifecycle.f90 for more details.

 write( *,'(/5x,"Initializing the CRTM (setTraj) ...")' )
 err_stat = CRTM_Init( self%rc%SENSOR_ID, &
            chinfo, &
            File_Path=trim(self%rc%COEFFICIENT_PATH), &
            Quiet=.TRUE.)
 if ( err_stat /= SUCCESS ) THEN
   message = 'Error initializing CRTM (setTraj)'
   call Display_Message( PROGRAM_NAME, message, FAILURE )
   stop
 end if

 ! Loop over all sensors. Not necessary if we're calling CRTM for each sensor
 ! ----------------------------------------------------------------------------
 Sensor_Loop:do n = 1, self%rc%n_Sensors


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
   call CRTM_Atmosphere_Create( atm, self%n_Layers, self%rc%n_Absorbers, self%rc%n_Clouds, self%rc%n_Aerosols )
   if ( ANY(.NOT. CRTM_Atmosphere_Associated(atm)) ) THEN
      message = 'Error allocating CRTM Forward Atmosphere structure (setTraj)'
      CALL Display_Message( PROGRAM_NAME, message, FAILURE )
      STOP
   END IF

   ! Create the input FORWARD structure (sfc)
   ! ----------------------------------------
!   call CRTM_Surface_Create(sfc, self%N_Channels)
!   IF ( ANY(.NOT. CRTM_Surface_Associated(sfc)) ) THEN
!      message = 'Error allocating CRTM Surface structure (setTraj)'
!      CALL Display_Message( PROGRAM_NAME, message, FAILURE )
!      STOP
!   END IF


   ! Create output K-MATRIX structure (atm)
   ! --------------------------------------
   call CRTM_Atmosphere_Create( self%atm_K, self%n_Layers, self%rc%n_Absorbers, self%rc%n_Clouds, self%rc%n_Aerosols )
   if ( ANY(.NOT. CRTM_Atmosphere_Associated(self%atm_K)) ) THEN
      message = 'Error allocating CRTM K-matrix Atmosphere structure (setTraj)'
      CALL Display_Message( PROGRAM_NAME, message, FAILURE )
      STOP
   END IF

   ! Create output K-MATRIX structure (sfc)
   ! --------------------------------------
!   call CRTM_Surface_Create(self%sfc_K, self%N_Channels)
!   IF ( ANY(.NOT. CRTM_Surface_Associated(self%sfc_K)) ) THEN
!      message = 'Error allocating CRTM K-matrix Surface structure (setTraj)'
!      CALL Display_Message( PROGRAM_NAME, message, FAILURE )
!      STOP
!   END IF


   !Assign the data from the GeoVaLs
   !--------------------------------
   CALL Load_Atm_Data(self%N_PROFILES,self%N_LAYERS,geovals,atm,self%rc)
!   call Load_Sfc_Data(self%N_PROFILES,self%N_LAYERS,self%N_Channels,geovals,sfc,chinfo,obss)
!   call Load_Geom_Data(obss,geo)

   IF (TRIM(self%rc%aerosol_option) /= "") &
        &CALL load_aerosol_data(self%n_profiles,self%n_layers,geovals,&
        &self%rc%aerosol_option,atm)

   call CRTM_Atmosphere_Inspect(atm(12))

   CALL CRTM_RTSolution_Create(rts, self%n_Layers )
   CALL CRTM_RTSolution_Create(rts_k, self%n_Layers )

   DO m = 1, self%n_Profiles
      DO l = 1, self%N_Channels
         rts_k(l,m)%layer_optical_depth = one
      ENDDO
   ENDDO

   ! Zero the K-matrix OUTPUT structures
   ! -----------------------------------
   call CRTM_Atmosphere_Zero( self%atm_K )
!   call CRTM_Surface_Zero( self%sfc_K )


   ! Inintialize the K-matrix INPUT so that the results are daero/dx
   ! -------------------------------------------------------------


   FORALL (m=1:self%N_PROFILES,l=1:self%N_Channels) rts_k(l,m)%layer_optical_depth = one

   ! Call the K-matrix model
   ! -----------------------
   err_stat = CRTM_AOD_K( atm         , &  ! FORWARD  Input
                             rts_K       , &  ! K-MATRIX Input
                             chinfo(n:n) , &  ! Input
                             rts         , &  ! FORWARD  Output
                             self%atm_K    )  ! K-MATRIX Output
   if ( err_stat /= SUCCESS ) THEN
      message = 'Error calling CRTM (setTraj) K-Matrix Model for '//TRIM(self%rc%SENSOR_ID(n))
      call Display_Message( PROGRAM_NAME, message, FAILURE )
      stop
   end if

   ! Deallocate the structures
   ! -------------------------
!   call CRTM_Geometry_Destroy(geo)
   call CRTM_Atmosphere_Destroy(atm)
   call CRTM_RTSolution_Destroy(rts_K)
   call CRTM_RTSolution_Destroy(rts)
!   call CRTM_Surface_Destroy(sfc)



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

end subroutine ufo_aod_tlad_settraj

! ------------------------------------------------------------------------------

subroutine ufo_aod_simobs_tl(self, geovals, hofx, obss)

implicit none
class(ufo_aod_tlad), intent(in) :: self
type(ufo_geovals),        intent(in) :: geovals
real(c_double),        intent(inout) :: hofx(:)
type(c_ptr), value,    intent(in)    :: obss

character(len=*), parameter :: myname_="ufo_aod_simobs_tl"
character(max_string) :: err_msg
integer :: job, jprofile, jchannel, jlevel
type(ufo_geoval), pointer :: tv_d


 ! Initial checks
 ! --------------

 ! Check if trajectory was set
 if (.not. self%ltraj) then
   write(err_msg,*) myname_, ' trajectory wasnt set!'
   call abor1_ftn(err_msg)
 endif

 ! Check if nobs is consistent in geovals & hofx
 if (geovals%nobs /= self%n_Profiles) then
   write(err_msg,*) myname_, ' error: nobs inconsistent!'
   call abor1_ftn(err_msg)
 endif

 ! Initialize hofx
 ! ---------------
 hofx(:) = 0.0_kind_real

 PRINT *,'@@@2',SIZE(hofx),self%n_Profiles,self%n_Channels
 WRITE(err_msg,*) myname_, ' hofx size'
 CALL abor1_ftn(err_msg)


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
 job = 0
 do jprofile = 1, self%n_Profiles
   do jchannel = 1, self%n_Channels
     job = job + 1
     do jlevel = 1, tv_d%nval
       hofx(job) = hofx(job) + &
                    self%atm_K(jchannel,jprofile)%Temperature(jlevel) * tv_d%vals(jlevel,jprofile)
     enddo
   enddo
 enddo


end subroutine ufo_aod_simobs_tl

! ------------------------------------------------------------------------------

subroutine ufo_aod_simobs_ad(self, geovals, hofx, obss)

implicit none
class(ufo_aod_tlad), intent(in) :: self
type(ufo_geovals),     intent(inout) :: geovals
real(c_double),           intent(in) :: hofx(:)
type(c_ptr), value,    intent(in)    :: obss

character(len=*), parameter :: myname_="ufo_aod_simobs_ad"
character(max_string) :: err_msg
integer :: job, jprofile, jchannel, jlevel
type(ufo_geoval), pointer :: tv_d


 ! Initial checks
 ! --------------

 ! Check if trajectory was set
 if (.not. self%ltraj) then
   write(err_msg,*) myname_, ' trajectory wasnt set!'
   call abor1_ftn(err_msg)
 endif

 ! Check if nobs is consistent in geovals & hofx
 if (geovals%nobs /= self%n_Profiles) then
   write(err_msg,*) myname_, ' error: nobs inconsistent!'
   call abor1_ftn(err_msg)
 endif


 ! Temperature
 ! -----------

 ! Get t from geovals
 call ufo_geovals_get_var(geovals, var_ts, tv_d)

 ! allocate if not yet allocated
 if (.not. allocated(tv_d%vals)) then
    tv_d%nobs = self%n_Profiles
    tv_d%nval = self%n_Layers
    allocate(tv_d%vals(tv_d%nval,tv_d%nobs))
    tv_d%vals = 0.0_kind_real
 endif


 ! Multiply by Jacobian and add to hofx (adjoint)
 job = 0
 do jprofile = 1, self%n_Profiles
   do jchannel = 1, self%n_Channels
     job = job + 1
     do jlevel = 1, tv_d%nval
       tv_d%vals(jlevel,jprofile) = tv_d%vals(jlevel,jprofile) + &
                                    self%atm_K(jchannel,jprofile)%Temperature(jlevel) * hofx(job)
     enddo
   enddo
 enddo


 ! Once all geovals set replace flag
 ! ---------------------------------
 if (.not. geovals%linit ) geovals%linit=.true.


end subroutine ufo_aod_simobs_ad

! ------------------------------------------------------------------------------

end module ufo_aod_tlad_mod
