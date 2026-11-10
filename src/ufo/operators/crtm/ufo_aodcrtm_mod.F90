! (C) Copyright 2017-2018 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

!> Fortran module to handle aod observations

module ufo_aodcrtm_mod

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
 type, public :: ufo_aodcrtm
 private
  ! variables requested from the model
  character(len=MAXVARLEN), public, allocatable :: varin(:)
  integer, allocatable                          :: channels(:)
  type(crtm_conf) :: conf
 contains
   procedure :: setup  => ufo_aodcrtm_setup
   procedure :: delete => ufo_aodcrtm_delete
   procedure :: simobs => ufo_aodcrtm_simobs
 end type ufo_aodcrtm
 character(len=maxvarlen), dimension(5), &
           parameter :: varin_default = [var_ts, var_mixr, var_rh, var_prs, var_prsi]

 character(MAXVARLEN), parameter :: varname_tmplate="aerosol_optical_depth"

contains

! ------------------------------------------------------------------------------

subroutine ufo_aodcrtm_setup(self, f_confOper, channels, midPointJulday, comm)

implicit none
class(ufo_aodcrtm),        intent(inout) :: self
type(fckit_configuration), intent(in)    :: f_confOper
integer(c_int),            intent(in)    :: channels(:)  !List of channels to use
integer(c_int64_t),        intent(in)    :: midPointJulday
type(fckit_mpi_comm),      intent(in)    :: comm

type(fckit_configuration) :: f_confOpts

 call f_confOper%get_or_die("obs options",f_confOpts)

 call ufo_crtm_conf_setup(self%conf, f_confOpts, f_confOper, midPointJulday, comm)
 call ufo_crtm_check_absorbers(self%conf, .true., "ufo_aodcrtm_setup")
 if (self%conf%n_Aerosols == 0) then
   call abor1_ftn("ufo_aodcrtm_setup error: AerosolOption must be set")
 end if
 self%varin = [character(len=MAXVARLEN) :: varin_default, self%conf%Aerosols]
 self%channels = channels

 call f_confOpts%final()

 ! use the shared CRTM coefficient load in simobs
 call ufo_crtm_register_aerosol(self%conf)

end subroutine ufo_aodcrtm_setup

! ------------------------------------------------------------------------------

subroutine ufo_aodcrtm_delete(self)

implicit none
class(ufo_aodcrtm), intent(inout) :: self

 call ufo_crtm_deregister_aerosol(self%conf)
 call ufo_crtm_conf_delete(self%conf)

end subroutine ufo_aodcrtm_delete

! ------------------------------------------------------------------------------

subroutine ufo_aodcrtm_simobs(self, geovals, obss, nvars, nlocs, hofx)

implicit none
class(ufo_aodcrtm),       intent(in)    :: self
type(ufo_geovals),        intent(in)    :: geovals
integer,                  intent(in)    :: nvars, nlocs
real(c_double),           intent(inout) :: hofx(nvars, nlocs)
type(c_ptr), value,       intent(in)    :: obss

! Local Variables
character(len=*), parameter :: PROGRAM_NAME = "ufo_aodcrtm_simobs"
character(len=MAXVARLEN) :: def_aero_mod
type(fckit_mpi_comm)  :: f_comm
character(255) :: message
integer        :: err_stat
integer        :: l, m, n
real(c_double) :: missing

integer :: n_Profiles
integer :: n_Layers
integer :: n_Channels

type(CRTM_Options_type), allocatable :: Options(:)

! Define the Channel Info argument
type(CRTM_ChannelInfo_type)             :: chinfo(self%conf%n_Sensors)

! Define the FORWARD variables
type(CRTM_Atmosphere_type), allocatable :: atm(:)
type(CRTM_RTSolution_type), allocatable :: rts(:,:)

 ! Get number of profile and layers from geovals
 ! ---------------------------------------------
 n_Profiles = geovals%nlocs
 n_Layers = ufo_geovals_nlevs(geovals, var_ts)

 ! ChannelInfo for this operator's sensors, from the shared CRTM coefficient
 ! load (CRTM_Init runs only on first use; see ufo_crtm_utils_mod)
 ! ---------------------------------------------------------------------------

 call obsspace_get_comm(obss, f_comm)
 call ufo_crtm_init_aerosol(self%conf, chinfo, PROGRAM_NAME, f_comm, def_aero_mod)


 ! Set missing value; profiles that are skipped keep it
 missing = missing_value(missing)
 hofx = missing

 ! Loop over sensors (one per operator in practice)
 ! ----------------------------------------------------------------------------
 Sensor_Loop:do n = 1, self%conf%n_Sensors


   ! Determine the number of channels for the current sensor
   ! -------------------------------------------------------
   n_Channels = CRTM_ChannelInfo_n_Channels(chinfo(n))


   ! Allocate the ARRAYS
   ! -------------------
   allocate(atm(n_Profiles))
   allocate(rts(n_Channels, n_Profiles))


   ! Create the input FORWARD structure (atm)
   ! ----------------------------------------
   call CRTM_Atmosphere_Create( atm, n_Layers, self%conf%n_Absorbers, self%conf%n_Clouds, &
                               self%conf%n_Aerosols )
   if ( any(.not. CRTM_Atmosphere_Associated(atm)) ) then
      message = "Error allocating CRTM Forward Atmosphere structure"
      call ufo_crtm_comm_stat_check(FAILURE, PROGRAM_NAME, message, f_comm)
   end if

   call CRTM_RTSolution_Create(rts, n_Layers )


   !Assign the data from the GeoVaLs
   !--------------------------------
   call ufo_crtm_load_atm_data(n_Profiles,n_Layers,geovals,atm,self%conf)

   allocate(Options(n_Profiles))
   do m = 1, n_Profiles
     Options(m)%Skip_Profile = .not. ufo_crtm_pressure_is_monotonic(atm(m))
   end do

   if (self%conf%aerosol_option /= "") then
     call ufo_crtm_load_aerosol_data(n_profiles, n_layers, geovals, &
                                     self%conf, self%conf%Aerosols, trim(def_aero_mod), atm)
   end if

   ! Call THE CRTM inspection
   ! ------------------------
   if (self%conf%inspect > 0) then
     call CRTM_Atmosphere_Inspect(atm(self%conf%inspect))
     call CRTM_ChannelInfo_Inspect(chinfo(n))
   end if

   ! The Forward Operator y = H(x) for AOD
   ! (the AOD K-matrix, CRTM_AOD_K, is only needed by the TL/AD operator)
   ! ------------------------------------------
   err_stat = CRTM_AOD( atm          , &  ! FORWARD  Input
                        chinfo(n:n)  , &  ! Input
                        rts          , &  ! FORWARD  Output
                        Options      )    ! Input
   message = "Error calling CRTM Forward Model for "//self%conf%SENSOR_ID(n)
   call ufo_crtm_comm_stat_check(err_stat, PROGRAM_NAME, message, f_comm)


   ! Put simulated aerosol optical depth into hofx
   ! ----------------------------------------------
   do m = 1, n_profiles
       if (Options(m)%Skip_Profile) cycle
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
   deallocate(atm)
   deallocate(rts)
   deallocate(Options)

end do Sensor_Loop

 ! Destroy CRTM instance
 ! ---------------------
 call ufo_crtm_destroy(chinfo, PROGRAM_NAME, f_comm)
 call f_comm%final()

end subroutine ufo_aodcrtm_simobs

! ------------------------------------------------------------------------------

end module ufo_aodcrtm_mod
