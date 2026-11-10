! (C) Copyright 2017-2018 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

!> Fortran module to handle tl/ad for aod observations

module ufo_aodcrtm_tlad_mod

 use fckit_configuration_module, only: fckit_configuration
 use fckit_mpi_module,   only: fckit_mpi_comm
 use, intrinsic :: iso_c_binding
 use kinds
 use missing_values_mod

 use ufo_geovals_mod, only: ufo_geovals, ufo_geoval, ufo_geovals_get_var, ufo_geovals_nlevs
 use ufo_vars_mod
 use ufo_crtm_utils_mod
 use CRTM_Module
 use obsspace_mod

 implicit none
 private

 !> Fortran derived type for aod trajectory
 type, public :: ufo_aodcrtm_tlad
 private
  character(len=MAXVARLEN), public, allocatable :: varin(:)  ! variables requested from the model
  integer, allocatable                          :: channels(:)
  type(crtm_conf) :: conf
  integer :: n_Profiles
  integer :: n_Layers
  integer :: n_Channels
  ! Aerosol concentration Jacobians for the requested channels, packed from the
  ! CRTM K-matrix output in settraj: (n_Layers, size(channels), n_Profiles, n_Aerosols)
  real(kind_real), allocatable :: K_aero(:,:,:,:)
  real(kind_real), allocatable  :: scaling_factor(:,:)
  type(CRTM_Options_type), allocatable :: Options(:)
  logical :: ltraj
 contains
  procedure :: setup  => ufo_aodcrtm_tlad_setup
  procedure :: delete  => ufo_aodcrtm_tlad_delete
  procedure :: settraj => ufo_aodcrtm_tlad_settraj
  procedure :: simobs_tl  => ufo_aodcrtm_simobs_tl
  procedure :: simobs_ad  => ufo_aodcrtm_simobs_ad
 end type ufo_aodcrtm_tlad

contains

! ------------------------------------------------------------------------------

subroutine ufo_aodcrtm_tlad_setup(self, f_confOper, channels, midPointJulday, comm)

implicit none
class(ufo_aodcrtm_tlad),   intent(inout) :: self
type(fckit_configuration), intent(in)    :: f_confOper
integer(c_int),            intent(in)    :: channels(:)  !List of channels to use
integer(c_int64_t),        intent(in)    :: midPointJulday
type(fckit_mpi_comm),      intent(in)    :: comm

type(fckit_configuration) :: f_confOpts

 call f_confOper%get_or_die("obs options",f_confOpts)

 call ufo_crtm_conf_setup(self%conf, f_confOpts, f_confOper, midPointJulday, comm)
 call ufo_crtm_check_absorbers(self%conf, .true., "ufo_aodcrtm_tlad_setup")
 if (self%conf%n_Aerosols == 0) then
   call abor1_ftn("ufo_aodcrtm_tlad_setup error: AerosolOption must be set")
 end if

 ! increment variables requested from the model: the aerosol species
 self%varin = self%conf%Aerosols
 self%channels = channels

 call f_confOpts%final()

 ! use the shared CRTM coefficient load in settraj
 call ufo_crtm_register_aerosol(self%conf)

end subroutine ufo_aodcrtm_tlad_setup

! ------------------------------------------------------------------------------

subroutine ufo_aodcrtm_tlad_delete(self)

implicit none
class(ufo_aodcrtm_tlad), intent(inout) :: self

 self%ltraj = .false.
 call ufo_crtm_deregister_aerosol(self%conf)
 call ufo_crtm_conf_delete(self%conf)

 if (allocated(self%K_aero)) deallocate(self%K_aero)

 if (allocated(self%scaling_factor)) then
   deallocate(self%scaling_factor)
 end if

 if (allocated(self%Options)) deallocate(self%Options)

end subroutine ufo_aodcrtm_tlad_delete

! ------------------------------------------------------------------------------

subroutine ufo_aodcrtm_tlad_settraj(self, geovals, obss)

implicit none

class(ufo_aodcrtm_tlad), intent(inout) :: self
type(ufo_geovals),        intent(in)    :: geovals
type(c_ptr), value,       intent(in)    :: obss

! Local Variables
character(len=*), parameter :: PROGRAM_NAME = "ufo_aodcrtm_tlad_settraj"
character(len=MAXVARLEN) :: def_aero_mod
type(fckit_mpi_comm)  :: f_comm
character(255) :: message
integer        :: err_stat
integer        :: n,l,m,jchannel,jaero

! Define the Channel Info argument
type(CRTM_ChannelInfo_type)             :: chinfo(self%conf%n_Sensors)

! Define the FORWARD variables
type(CRTM_Atmosphere_type), allocatable :: atm(:)
type(CRTM_RTSolution_type), allocatable :: rts(:,:)

! Define the K-MATRIX variables
type(CRTM_Atmosphere_type), allocatable :: atm_K(:,:)
type(CRTM_RTSolution_type), allocatable :: rts_K(:,:)


 ! Get number of profile and layers from geovals
 ! ---------------------------------------------
 self%n_Profiles = geovals%nlocs
 self%n_Layers = ufo_geovals_nlevs(geovals, var_ts)

 ! ChannelInfo for this operator's sensors, from the shared CRTM coefficient
 ! load (CRTM_Init runs only on first use; see ufo_crtm_utils_mod)
 ! ---------------------------------------------------------------------------

 call obsspace_get_comm(obss, f_comm)
 call ufo_crtm_init_aerosol(self%conf, chinfo, PROGRAM_NAME, f_comm, def_aero_mod)

 ! Loop over sensors (one per operator in practice)
 ! ----------------------------------------------------------------------------
 Sensor_Loop:do n = 1, self%conf%n_Sensors


   ! Determine the number of channels for the current sensor
   ! -------------------------------------------------------
   self%n_Channels = CRTM_ChannelInfo_n_Channels(chinfo(n))


   ! Allocate the ARRAYS
   ! -------------------
   allocate(atm(self%n_Profiles))
   allocate(rts(self%n_Channels, self%n_Profiles))
   allocate(atm_K(self%n_Channels, self%n_Profiles))
   allocate(self%scaling_factor(self%n_Layers, self%n_Profiles))
   allocate(self%Options(self%n_Profiles))
   allocate(rts_K(self%n_Channels, self%n_Profiles))


   ! Create the input FORWARD structure (atm)
   ! ----------------------------------------
   call CRTM_Atmosphere_Create( atm, self%n_Layers, self%conf%n_Absorbers, self%conf%n_Clouds, &
                               self%conf%n_Aerosols )
   if ( any(.not. CRTM_Atmosphere_Associated(atm)) ) then
      message = "Error allocating CRTM Forward Atmosphere structure (setTraj)"
      call ufo_crtm_comm_stat_check(FAILURE, PROGRAM_NAME, message, f_comm)
   end if

   ! Create output K-MATRIX structure (atm)
   ! --------------------------------------
   call CRTM_Atmosphere_Create( atm_K, self%n_Layers, self%conf%n_Absorbers, &
                               self%conf%n_Clouds, self%conf%n_Aerosols )
   if ( any(.not. CRTM_Atmosphere_Associated(atm_K)) ) then
      message = "Error allocating CRTM K-matrix Atmosphere structure (setTraj)"
      call ufo_crtm_comm_stat_check(FAILURE, PROGRAM_NAME, message, f_comm)
   end if

   !Assign the data from the GeoVaLs
   !--------------------------------
   call ufo_crtm_load_atm_data(self%n_Profiles,self%n_Layers,geovals,atm,self%conf)

   do m = 1, self%n_Profiles
     self%Options(m)%Skip_Profile = .not. ufo_crtm_pressure_is_monotonic(atm(m))
   end do

   if (self%conf%aerosol_option /= "") then
     call ufo_crtm_load_aerosol_data(self%n_profiles, self%n_layers, geovals, &
                                     self%conf, self%conf%Aerosols, trim(def_aero_mod), atm)
   end if

   call CRTM_RTSolution_Create(rts, self%n_layers )
   call CRTM_RTSolution_Create(rts_k, self%n_layers )

   ! (atm_K, the K-matrix OUTPUT, is already zeroed by CRTM_Atmosphere_Create)
   call ufo_crtm_calculate_aero_layer_factor(atm,self%scaling_factor)

   ! Initialize the K-matrix INPUT so that the results are daero/dx
   ! -------------------------------------------------------------

   forall (l=1:self%n_channels,m=1:self%n_profiles) rts_k(l,m)%layer_optical_depth = one

   ! Call the K-matrix model
   ! -----------------------
   err_stat = CRTM_AOD_K( atm         , &  ! FORWARD  Input
                             rts_K       , &  ! K-MATRIX Input
                             chinfo(n:n) , &  ! Input
                             rts         , &  ! FORWARD  Output
                             atm_K       , &  ! K-MATRIX Output
                             self%Options  )  ! Input
   message = "Error calling CRTM (setTraj) K-Matrix Model for "//self%conf%SENSOR_ID(n)
   call ufo_crtm_comm_stat_check(err_stat, PROGRAM_NAME, message, f_comm)

   ! Keep only the aerosol Jacobians of the requested channels, in a contiguous array
   ! ---------------------------------------------------------------------------------
   if (allocated(self%K_aero)) deallocate(self%K_aero)
   allocate(self%K_aero(self%n_Layers, size(self%channels), self%n_Profiles, self%conf%n_Aerosols))
   do jaero = 1, self%conf%n_Aerosols
     do m = 1, self%n_Profiles
       do jchannel = 1, size(self%channels)
         self%K_aero(:,jchannel,m,jaero) = &
             atm_K(self%channels(jchannel),m)%Aerosol(jaero)%Concentration(1:self%n_Layers)
       end do
     end do
   end do

   ! Deallocate the structures
   ! -------------------------
   call CRTM_Atmosphere_Destroy(atm)
   call CRTM_Atmosphere_Destroy(atm_K)
   call CRTM_RTSolution_Destroy(rts_K)
   call CRTM_RTSolution_Destroy(rts)



   ! Deallocate all arrays
   ! ---------------------
   deallocate(atm)
   deallocate(atm_K)
   deallocate(rts)
   deallocate(rts_K)

 end do Sensor_Loop


 ! Destroy CRTM instance
 ! ---------------------
 call ufo_crtm_destroy(chinfo, PROGRAM_NAME, f_comm)
 call f_comm%final()


 ! Set flag that the trajectory was set
 ! ------------------------------------
 self%ltraj = .true.

end subroutine ufo_aodcrtm_tlad_settraj

! ------------------------------------------------------------------------------

subroutine ufo_aodcrtm_simobs_tl(self, geovals, obss, nvars, nlocs, hofx)

implicit none
class(ufo_aodcrtm_tlad), intent(in) :: self
type(ufo_geovals),        intent(in) :: geovals
type(c_ptr), value,    intent(in)    :: obss
integer,                  intent(in)    :: nvars, nlocs
real(c_double),        intent(inout) :: hofx(nvars, nlocs)

character(len=*), parameter :: myname_="ufo_aodcrtm_simobs_tl"
character(max_string) :: err_msg
integer :: jprofile, jchannel, jlevel, jaero
type(ufo_geoval), pointer :: var_p


 ! Initial checks
 ! --------------

 ! Check if trajectory was set
 if (.not. self%ltraj) then
   err_msg = myname_ // " trajectory wasn't set!"
   call abor1_ftn(err_msg)
 end if

 ! Check if nlocs is consistent in geovals & hofx
 if (geovals%nlocs /= self%n_Profiles) then
   err_msg = myname_ // " error: nlocs inconsistent!"
   call abor1_ftn(err_msg)
 end if

 call ufo_geovals_get_var(geovals, self%conf%Aerosols(1), var_p)

 ! Check model levels is consistent in geovals & crtm
 if (var_p%nval /= self%n_Layers) then
   err_msg = myname_ // " error: layers inconsistent!"
   call abor1_ftn(err_msg)
 end if


 ! Initialize hofx
 ! ---------------
 hofx(:,:) = 0.0_kind_real

 ! Multiply by Jacobian and add to hofx
 ! (aerosols outermost so each GeoVaL is looked up once; every hofx(jchannel,jprofile)
 !  still accumulates aerosol by aerosol, level by level)
 do jaero = 1, self%conf%n_aerosols
   call ufo_geovals_get_var(geovals, self%conf%Aerosols(jaero), var_p)
   do jprofile = 1, self%n_profiles
     if (self%Options(jprofile)%Skip_Profile) cycle
     do jchannel = 1, size(self%channels)
       do jlevel = 1, var_p%nval
         hofx(jchannel, jprofile) = hofx(jchannel, jprofile) + &
             self%K_aero(jlevel,jchannel,jprofile,jaero) * &
             var_p%vals(jlevel,jprofile) * self%scaling_factor(jlevel,jprofile) * &
             self%conf%unit_coef
       end do
     end do
   end do
 end do

end subroutine ufo_aodcrtm_simobs_tl

! ------------------------------------------------------------------------------

subroutine ufo_aodcrtm_simobs_ad(self, geovals, obss, nvars, nlocs, hofx)

implicit none
class(ufo_aodcrtm_tlad), intent(in) :: self
type(ufo_geovals),     intent(inout) :: geovals
type(c_ptr), value,    intent(in)    :: obss
integer,                  intent(in)    :: nvars, nlocs
real(c_double),           intent(in) :: hofx(nvars, nlocs)

character(len=*), parameter :: myname_="ufo_aodcrtm_simobs_ad"
character(max_string) :: err_msg
integer :: jprofile, jchannel, jlevel
type(ufo_geoval), pointer :: var_p
real(c_double) :: missing

integer :: jaero

 ! Initial checks
 ! --------------

 ! Check if trajectory was set
 if (.not. self%ltraj) then
   err_msg = myname_ // " trajectory wasn't set!"
   call abor1_ftn(err_msg)
 end if

 ! Check if nlocs is consistent in geovals & hofx
 if (geovals%nlocs /= self%n_Profiles) then
   err_msg = myname_ // " error: nlocs inconsistent!"
   call abor1_ftn(err_msg)
 end if

 ! Set missing value
 missing = missing_value(missing)

 do jaero=1,self%conf%n_aerosols

    call ufo_geovals_get_var(geovals, self%conf%Aerosols(jaero), var_p)

! Multiply by Jacobian and add to the GeoVaLs (adjoint): dx = dx + K^T hofx,
! with the layer scaling applied to the added term only
    do jprofile = 1, self%n_Profiles
      if (self%Options(jprofile)%Skip_Profile) cycle
       do jchannel = 1, size(self%channels)
          if (hofx(jchannel, jprofile) /= missing) then
            do jlevel = 1, var_p%nval
               var_p%vals(jlevel,jprofile) = var_p%vals(jlevel,jprofile) + &
                   self%K_aero(jlevel,jchannel,jprofile,jaero) * &
                   self%scaling_factor(jlevel,jprofile) * self%conf%unit_coef * &
                   hofx(jchannel, jprofile)
            end do
          end if
       end do
    end do

 end do

 ! Once all geovals set replace flag
 ! ---------------------------------
 if (.not. geovals%linit ) geovals%linit=.true.


end subroutine ufo_aodcrtm_simobs_ad

! ------------------------------------------------------------------------------

end module ufo_aodcrtm_tlad_mod
