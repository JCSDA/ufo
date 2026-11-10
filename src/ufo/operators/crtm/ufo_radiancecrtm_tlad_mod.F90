! (C) Copyright 2017-2026 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

!> Fortran module to handle tl/ad for radiancecrtm observations

module ufo_radiancecrtm_tlad_mod
 use,intrinsic :: iso_c_binding
 use,intrinsic :: iso_fortran_env
 use CRTM_Module

 use fckit_configuration_module, only: fckit_configuration
 use fckit_mpi_module,   only: fckit_mpi_comm
 use, intrinsic :: iso_c_binding
 use kinds
 use missing_values_mod

 use obsspace_mod
 use obsdatavector_mod
 use ufo_geovals_mod, only: ufo_geovals, ufo_geoval, ufo_geovals_get_var, ufo_geovals_nlevs
 use ufo_vars_mod
 use ufo_crtm_utils_mod

 use ufo_constants_mod, only: deg2rad, kg_to_g
 use ufo_crtm_passive_mod
 use ufo_crtm_active_mod
 use ufo_reconradop_crtm_mod
 implicit none
 private

 !> Fortran derived type for radiancecrtm trajectory
 type, public :: ufo_radiancecrtm_tlad
 private
  type(ufo_reconradop_crtm) :: reconradop_crtm
  character(len=MAXVARLEN), public, allocatable :: varin(:)  ! variables requested from the model
  integer, allocatable                          :: channels(:)
  type(crtm_conf) :: conf
  type(crtm_conf) :: conf_traj
  integer :: n_Profiles
  integer :: n_Layers
  integer :: n_Channels
  ! Jacobians used by simobs_tl/ad, packed from the CRTM K-matrix output in settraj.
  ! Species are in the order of conf (the linear operator's lists).
  real(kind_real), allocatable :: K_temp(:,:,:)   ! (n_Layers, n_Channels, n_Profiles)
  real(kind_real), allocatable :: K_abs(:,:,:,:)  ! (n_Layers, n_Channels, n_Profiles, n_Absorbers)
  real(kind_real), allocatable :: K_cld(:,:,:,:)  ! (n_Layers, n_Channels, n_Profiles, n_Clouds)
  real(kind_real), allocatable :: K_sfc(:,:,:)    ! (n_Channels, n_Profiles, n_Surfaces)
  logical :: ltraj
  type(CRTM_Options_type), allocatable :: Options(:)
  logical :: use_qc_flags
 contains
  procedure :: setup  => ufo_radiancecrtm_tlad_setup
  procedure :: delete  => ufo_radiancecrtm_tlad_delete
  procedure :: settraj => ufo_radiancecrtm_tlad_settraj
  procedure :: simobs_tl  => ufo_radiancecrtm_simobs_tl
  procedure :: simobs_ad  => ufo_radiancecrtm_simobs_ad
 end type ufo_radiancecrtm_tlad

 character(len=maxvarlen), dimension(1), parameter :: varin_default = &
                            [var_ts]

contains

! ------------------------------------------------------------------------------

subroutine ufo_radiancecrtm_tlad_setup(self, f_confOper, channels, midPointJulday, comm)

implicit none
class(ufo_radiancecrtm_tlad), intent(inout) :: self
type(fckit_configuration),    intent(in)    :: f_confOper
integer(c_int),               intent(in)    :: channels(:)  !List of channels to use
integer(c_int64_t),           intent(in)    :: midPointJulday
type(fckit_mpi_comm),         intent(in)    :: comm

type(fckit_configuration) :: f_confOpts,f_confLinOper
character(max_string) :: err_msg

 call f_confOper%get_or_die("obs options",f_confOpts)
 call ufo_crtm_conf_setup(self%conf_traj, f_confOpts, f_confOper, midPointJulday, comm)
 call f_confOper%get_or_die("UseQCFlagsToSkipHofX",self%use_qc_flags)

 call ufo_crtm_check_absorbers(self%conf_traj, .false., "ufo_radiancecrtm_tlad_setup")

 if ( f_confOper%has("linear obs operator") ) then
    call f_confOper%get_or_die("linear obs operator",f_confLinOper)
    call ufo_crtm_conf_setup(self%conf, f_confOpts, f_confLinOper, midPointJulday, comm)
    call f_confLinOper%final()
 else
    call ufo_crtm_conf_setup(self%conf, f_confOpts, f_confOper, midPointJulday, comm)
 end if

 ! Increment variables requested from the model: varin_default, the absorbers, the cloud
 ! mass contents and the surface variables of the linear operator's conf
 self%varin = [character(len=MAXVARLEN) :: varin_default, self%conf%Absorbers, &
               self%conf%Clouds(:,1), self%conf%Surfaces]
 if ( (ufo_vars_getindex(self%varin, var_sfc_wspeed) > 0 .or. &
       ufo_vars_getindex(self%varin, var_sfc_wdir) > 0) .and. &
      self%conf_traj%sfc_wind_geovars /= "vector") then
   err_msg = "ufo_radiancecrtm_tlad_setup error: " // &
             "sfc_wind_geovars not supported in tlad --> " // self%conf_traj%sfc_wind_geovars
   call abor1_ftn(err_msg)
 end if


 self%channels = channels
 call f_confOpts%final()

 ! use the shared CRTM coefficient load in settraj (CRTM_Init uses conf_traj)
 call ufo_crtm_register_radiance(self%conf_traj)

end subroutine ufo_radiancecrtm_tlad_setup

! ------------------------------------------------------------------------------

subroutine ufo_radiancecrtm_tlad_delete(self)

implicit none
class(ufo_radiancecrtm_tlad), intent(inout) :: self

 self%ltraj = .false.

 call ufo_crtm_deregister_radiance(self%conf_traj)
 call ufo_crtm_conf_delete(self%conf)
 call ufo_crtm_conf_delete(self%conf_traj)

 call free_packed_jacobians(self)

 if (allocated(self%Options)) then
    call CRTM_Options_Destroy(self%Options)
    deallocate(self%Options)
 end if


end subroutine ufo_radiancecrtm_tlad_delete

! ------------------------------------------------------------------------------

subroutine ufo_radiancecrtm_tlad_settraj(self, geovals, obss, hofxdiags, qcf_p)
use fckit_mpi_module,   only: fckit_mpi_comm
use fckit_log_module,   only: fckit_log
use obsdatavector_mod,  only: obsdatavector_int
use, intrinsic :: ieee_arithmetic,    only: ieee_is_nan
use CRTM_SpcCoeff, only: SC, &
                         SpcCoeff_IsMicrowaveSensor , &
                         SpcCoeff_IsInfraredSensor  , &
                         SpcCoeff_IsVisibleSensor   , &
                         SpcCoeff_IsUltravioletSensor

implicit none

class(ufo_radiancecrtm_tlad), intent(inout) :: self
type(ufo_geovals),        intent(in)    :: geovals
type(c_ptr), value,       intent(in)    :: obss
type(ufo_geovals),        intent(inout) :: hofxdiags    !non-h(x) diagnostics
type(c_ptr), value,       intent(in) :: qcf_p
type(obsdatavector_int) :: qc_flags

! Local Variables
character(len=*), parameter :: PROGRAM_NAME = "ufo_radiancecrtm_tlad_settraj"
character(255) :: message
character(max_string) :: err_msg
integer        :: err_stat
integer        :: isc  ! this sensor's entry in CRTM's SpcCoeff table (SC)
integer        :: n
integer :: jvar, jchannel, jprofile, jlevel, jspec
type(fckit_mpi_comm)  :: f_comm

! Define the Channel Info and Geometry arguments
type(CRTM_ChannelInfo_type)             :: chinfo(self%conf_traj%n_Sensors)
type(CRTM_Geometry_type),   allocatable :: geo(:)

! Define the FORWARD variables
type(CRTM_Atmosphere_type), allocatable :: atm(:)
type(CRTM_Surface_type),    allocatable :: sfc(:)
type(CRTM_RTSolution_type), allocatable :: rts(:,:)

! Define the K-MATRIX variables
type(CRTM_Atmosphere_type), allocatable :: atm_K(:,:)
type(CRTM_Surface_type),    allocatable :: sfc_K(:,:)
type(CRTM_RTSolution_type), allocatable :: rts_K(:,:)

!for gmi
type(CRTM_Geometry_type),   allocatable :: geo_hf(:)
type(CRTM_ChannelInfo_type) :: chinfo_gmi(2)  ! channels 1-9 and 10-13
integer :: n_lo                               ! number of gmi channels 1-9

character(len=MAXVARLEN),  dimension(hofxdiags%nvar) :: &
                          ystr_diags, xstr_diags
integer :: str_pos(4), ch_diags(hofxdiags%nvar),numNaN
integer :: lch
! Used to parse hofxdiags
character(len=MAXVARLEN) :: varstr
character(10), parameter :: jacobianstr = "_jacobian_"
integer(c_size_t) :: nvars, nlocs
integer :: ispec
integer :: ivar

! set a local boolean variable for whether we are in vis or ultraviolet channels
logical :: Is_Vis_or_UV

integer, allocatable :: zeroCloudInCRTM0(:)

 call obsspace_get_comm(obss, f_comm)

 ! Get number of profile and layers from geovals
 ! ---------------------------------------------
 self%n_Profiles = geovals%nlocs
 self%n_Layers = ufo_geovals_nlevs(geovals, var_ts)

 qc_flags%data_ptr = qcf_p

 ! ChannelInfo for this operator's sensors, from the shared CRTM coefficient
 ! load (CRTM_Init runs only on first use; see ufo_crtm_utils_mod)
 ! ---------------------------------------------------------------------------

 call ufo_crtm_init_radiance(self%conf_traj, chinfo, self%channels, PROGRAM_NAME, f_comm)

 ! Loop over sensors (one per operator in practice)
 ! ----------------------------------------------------------------------------
 Sensor_Loop:do n = 1, self%conf_traj%n_Sensors


   ! Determine the number of channels for the current sensor
   ! -------------------------------------------------------
   self%n_Channels = CRTM_ChannelInfo_n_Channels(chinfo(n))
   isc = chinfo(n)%Sensor_Index


   ! Allocate the ARRAYS
   ! -------------------
   allocate(geo(self%n_Profiles))
   allocate(atm(self%n_Profiles))
   allocate(sfc(self%n_Profiles))
   allocate(rts(self%n_Channels, self%n_Profiles))
   allocate(atm_K(self%n_Channels, self%n_Profiles))
   allocate(sfc_K(self%n_Channels, self%n_Profiles))
   allocate(rts_K(self%n_Channels, self%n_Profiles))
   allocate(self%Options(self%n_Profiles))

   ! Create the input FORWARD structure (atm)
   ! ----------------------------------------
   call CRTM_Atmosphere_Create( atm, self%n_Layers, self%conf_traj%n_Absorbers, &
                               self%conf_traj%n_Clouds, self%conf_traj%n_Aerosols )
   if ( any(.not. CRTM_Atmosphere_Associated(atm)) ) then
      message = "Error allocating CRTM Forward Atmosphere structure (setTraj)"
      call ufo_crtm_comm_stat_check(FAILURE, PROGRAM_NAME, message, f_comm)
   end if

   if (self%n_Layers > 0) call CRTM_RTSolution_Create(rts, self%n_Layers )

   ! Create the input FORWARD structure (sfc)
   ! ----------------------------------------
   call CRTM_Surface_Create(sfc, self%n_Channels)
   if ( any(.not. CRTM_Surface_Associated(sfc)) ) then
      message = "Error allocating CRTM Surface structure (setTraj)"
      call ufo_crtm_comm_stat_check(FAILURE, PROGRAM_NAME, message, f_comm)
   end if


   ! Create output K-MATRIX structure (atm)
   ! --------------------------------------
   call CRTM_Atmosphere_Create( atm_K, self%n_Layers, self%conf_traj%n_Absorbers, &
                                self%conf_traj%n_Clouds, self%conf_traj%n_Aerosols )
   if ( any(.not. CRTM_Atmosphere_Associated(atm_K)) ) then
      message = "Error allocating CRTM K-matrix Atmosphere structure (setTraj)"
      call ufo_crtm_comm_stat_check(FAILURE, PROGRAM_NAME, message, f_comm)
   end if


   ! Create output K-MATRIX structure (sfc)
   ! --------------------------------------
   ! No SensorData: CRTM reads it only from the forward sfc, never from sfc_K
   call CRTM_Surface_Create(sfc_K, 0)

   ! Some special treatment for Vis or UV as compared to IR or MW.
   ! ----------------------------------------
   if (SpcCoeff_IsVisibleSensor(SC(isc)) .or. SpcCoeff_IsUltravioletSensor(SC(isc))) then
      Is_Vis_or_UV = .true.
   else
      Is_Vis_or_UV = .false.
   end if

   !Assign the data from the GeoVaLs
   !--------------------------------
   allocate(zeroCloudInCRTM0(self%n_Profiles))
   call ufo_crtm_load_sfc_data(self%n_Profiles, self%n_Channels, self%channels, geovals, sfc, &
                               chinfo, obss, self%conf_traj, SC(isc)%Is_Active_Sensor, &
                               Is_Vis_or_UV, zeroCloudInCRTM0)
   call ufo_crtm_load_atm_data(self%n_Profiles, self%n_Layers, geovals, atm, self%conf_traj, &
                               SC(isc)%Is_Active_Sensor, zeroCloudInCRTM0)
   deallocate(zeroCloudInCRTM0)
   if (self%conf%SENSOR_ID(n) == "gmi_gpm") then
      allocate( geo_hf( self%n_Profiles ))
      call ufo_crtm_load_geom_data(obss,geo,geo_hf,self%conf%SENSOR_ID(n))
   else
      call ufo_crtm_load_geom_data(obss,geo)
   end if

   ! Zero the K-matrix OUTPUT structures
   ! (atm_K is already zeroed by CRTM_Atmosphere_Create)
   ! -----------------------------------
   call CRTM_Surface_Zero( sfc_K )

   ! Initialize the K-matrix INPUT so that the results are dTb/dx or dR/dx
   ! -------------------------------------------------------------
   if (SC(isc)%Is_Active_Sensor) then
       do jchannel = 1, self%n_Channels
          do jprofile = 1, self%n_Profiles
             do jlevel = 1, self%n_Layers
                rts_K(jchannel,jprofile)%Reflectivity(jlevel)            = ZERO
                rts_K(jchannel,jprofile)%Reflectivity_Attenuated(jlevel) = ONE
             end do
          end do
       end do

       rts_K%Radiance                = ZERO
       rts_K%Brightness_Temperature  = ZERO
   else if (Is_Vis_or_UV .or. self % conf % read_Cmatrix ) then
       rts_K%Radiance                = ONE
       rts_K%Brightness_Temperature  = ZERO
   else
       rts_K%Radiance                = ZERO
       rts_K%Brightness_Temperature  = ONE
   end if


   call ufo_crtm_skip_profiles(self%n_Profiles,self%n_Channels,self%channels,obss,atm,sfc, &
                               SC(isc)%Is_Active_Sensor, Is_Vis_or_UV, self%Options)

   ! eliminate remaining profiles that are "QCed" out
   if ( self%use_qc_flags ) then
     call ufo_crtm_skip_qc_profiles(qc_flags, size(self%channels), &
                                    self%Options(1:self%n_Profiles), "DEBUG (TLAD)")
   end if

   ! Call the K-matrix model
   ! -----------------------
   if (self%conf%SENSOR_ID(n) == "gmi_gpm") then
      ! channels 1-9 with geo, channels 10-13 with geo_hf: one call per block
      call ufo_crtm_gmi_split_chinfo(chinfo(n), self%channels, chinfo_gmi, n_lo, &
                                     PROGRAM_NAME, f_comm)
      if (n_lo > 0) then
         err_stat = CRTM_K_Matrix( atm, sfc, rts_K(1:n_lo,:), geo, chinfo_gmi(1:1), &
                                   atm_K(1:n_lo,:), sfc_K(1:n_lo,:), rts(1:n_lo,:), &
                                   self%Options )
         message = "Error calling CRTM (setTraj) K-Matrix Model for gmi_gpm channels 1-9"
         call ufo_crtm_comm_stat_check(err_stat, PROGRAM_NAME, message, f_comm)
      end if
      if (n_lo < self%n_Channels) then
         ! K-matrix input so that the results are dTb/dx
         rts_K(n_lo+1:,:)%Radiance               = ZERO
         rts_K(n_lo+1:,:)%Brightness_Temperature = ONE
         err_stat = CRTM_K_Matrix( atm, sfc, rts_K(n_lo+1:,:), geo_hf, chinfo_gmi(2:2), &
                                   atm_K(n_lo+1:,:), sfc_K(n_lo+1:,:), rts(n_lo+1:,:), &
                                   self%Options )
         message = "Error calling CRTM (setTraj, geo_hf) K-Matrix Model for gmi_gpm channels 10-13"
         call ufo_crtm_comm_stat_check(err_stat, PROGRAM_NAME, message, f_comm)
      end if
      call CRTM_ChannelInfo_Destroy(chinfo_gmi)
   else
      err_stat = CRTM_K_Matrix( atm         , &  ! FORWARD  Input
                                sfc         , &  ! FORWARD  Input
                                rts_K       , &  ! K-MATRIX Input
                                geo         , &  ! Input
                                chinfo(n:n) , &  ! Input
                                atm_K       , &  ! K-MATRIX Output
                                sfc_K       , &  ! K-MATRIX Output
                                rts         , &  ! FORWARD  Output
                                self%Options  )  ! Input
      message = "Error calling CRTM (setTraj) K-Matrix Model for "//self%conf_traj%SENSOR_ID(n)
      call ufo_crtm_comm_stat_check(err_stat, PROGRAM_NAME, message, f_comm)
   end if

   ! check for NaN values in atm_k: skip the profile at its first NaN
   numNaN = 0
   NaN_Profile_Loop: do jprofile = 1, self%n_Profiles
      if (self%Options(jprofile)%Skip_Profile) cycle NaN_Profile_Loop
      do jchannel = 1, size(self%channels)
         if (any(ieee_is_nan(atm_K(jchannel,jprofile)% &
                             Temperature(1:atm_K(jchannel,jprofile)%n_layers)))) then
            self%Options(jprofile)%Skip_Profile = .true.
            numNaN = numNaN + 1
            write(message,"(i0,a)") numNaN, "th NaN in Jacobian Profiles"
            call fckit_log%info(message)
            cycle NaN_Profile_Loop
         end if
      end do
   end do NaN_Profile_Loop

   !! Parse hofxdiags%variables into independent/dependent variables and channel
   !! assumed formats:
   !!   jacobian var -->     <ystr>_jacobian_<xstr>_<chstr>
   !!   non-jacobian var --> <ystr>_<chstr>

   ch_diags = -9999
   do jvar = 1, hofxdiags%nvar
      varstr = hofxdiags%variables(jvar)
      str_pos(4) = len_trim(varstr)
      if (str_pos(4) < 1) cycle
      str_pos(3) = index(varstr,"_",back=.true.)        !final "_" before channel
      read(varstr(str_pos(3)+1:str_pos(4)),*, err=999) ch_diags(jvar)
 999  str_pos(1) = index(varstr,jacobianstr) - 1        !position before jacobianstr
      if (str_pos(1) == 0) then
         err_msg = "ufo_radiancecrtm_tlad_settraj: _jacobian_ must be preceded by dependent " // &
                   "variable in config: " // hofxdiags%variables(jvar)
         call abor1_ftn(err_msg)
      else if (str_pos(1) > 0) then
         !Diagnostic is a Jacobian member (dy/dx)
         ystr_diags(jvar) = varstr(1:str_pos(1))
         str_pos(2) = str_pos(1) + len(jacobianstr) + 1 !begin xstr_diags
         str_pos(4) = str_pos(3) - str_pos(2)
         xstr_diags(jvar)(1:str_pos(4)) = varstr(str_pos(2):str_pos(3)-1)
         xstr_diags(jvar)(str_pos(4)+1:) = ""
      else !null
         !Diagnostic is a dependent variable (y)
         xstr_diags(jvar) = ""
         ystr_diags(jvar)(1:str_pos(3)-1) = varstr(1:str_pos(3)-1)
         ystr_diags(jvar)(str_pos(3):) = ""
         if (ch_diags(jvar) < 0) ystr_diags(jvar) = varstr
      end if
   end do

   nvars = hofxdiags%nvar
   nlocs = geovals%nlocs
   ! put CRTM output into hofxdiags
   if (SC(isc)%Is_Active_Sensor) then
      call ufo_crtm_active_diag(rts, &
                                rts_K, &
                                atm, &
                                atm_K, &
                                sfc_K, &
                                self%conf, &
                                n, &
                                self%Options, &
                                self%channels, &
                                geovals, &
                                obss, &
                                nvars, &
                                nlocs, &
                                self%n_Profiles, &
                                self%n_Layers, &
                                xstr_diags, &
                                ystr_diags, &
                                ch_diags, &
                                hofxdiags,&
                                err_stat)
   else
      if (self % conf % read_Cmatrix) then
         call self%reconradop_crtm%apply(self%conf % Cmatrix_path,&
                                         rts, &
                                         isc, &
                                         self%n_Profiles, &
                                         self%n_Channels, &
                                         self%conf%n_Absorbers,&
                                         self%n_Layers,&
                                         self%channels, &
                                         rts_K, &
                                         atm_K, &
                                         sfc_K, .true. )
      end if

      call ufo_crtm_passive_diag(rts, &
                                 rts_K, &
                                 atm, &
                                 atm_K, &
                                 sfc_K, &
                                 self%conf, &
                                 n, &
                                 self%Options, &
                                 self%channels, &
                                 geovals, &
                                 obss, &
                                 nvars, &
                                 nlocs, &
                                 self%n_Profiles, &
                                 self%n_Layers, &
                                 xstr_diags, &
                                 ystr_diags, &
                                 ch_diags, &
                                 hofxdiags,&
                                 err_stat)
   end if

   ! check for error from either passive or active
   if (err_stat > 0) then
       err_msg = "ufo_radiancecrtm_tlad_settraj error: " // &
                 "failed to put simulated diagnostics into hofxdiags"
       call abor1_ftn(err_msg)
    end if

   ! If desired, nullify out Jacobian for a user-defined per-variable channel list
   ! ---------------------------------------------------------------------------------------
   if (size(self%conf%nullifyJacVars) > 0) then
     do jspec = 1, self%conf%n_Absorbers
       ivar = ufo_vars_getindex(self%conf%nullifyJacVars, self%conf%Absorbers(jspec))
       if (ivar > 0) then
          ispec = ufo_vars_getindex(self%conf_traj%Absorbers, self%conf%Absorbers(jspec))
          do jprofile = 1, self%n_Profiles
            if (.not.self%Options(jprofile)%Skip_Profile) then
               do lch = 1, size(self%channels)
                  if (any(self%conf%nullifyJacChannels(ivar)%channels == self%channels(lch))) then
                     atm_K(lch,jprofile)%Absorber(:,ispec) = ZERO
                  end if
               end do
            end if
          end do
       end if
     end do
   end if

   ! Keep only the Jacobians simobs_tl/ad need, in contiguous arrays
   ! ---------------------------------------------------------------
   call pack_jacobians(self, atm_K, sfc_K)

   ! Deallocate the structures
   ! -------------------------
   call CRTM_Geometry_Destroy(geo)
   call CRTM_Atmosphere_Destroy(atm)
   call CRTM_Atmosphere_Destroy(atm_K)
   call CRTM_RTSolution_Destroy(rts_K)
   call CRTM_RTSolution_Destroy(rts)
   call CRTM_Surface_Destroy(sfc)
   call CRTM_Surface_Destroy(sfc_K)


   ! Deallocate all arrays
   ! ---------------------
   deallocate(geo)
   deallocate(atm)
   deallocate(sfc)
   deallocate(rts)
   deallocate(rts_K)
   deallocate(atm_K)
   deallocate(sfc_K)
   if(allocated(geo_hf)) deallocate(geo_hf)

 end do Sensor_Loop


 ! Destroy CRTM instance
 ! ---------------------
 call ufo_crtm_destroy(chinfo, PROGRAM_NAME, f_comm)

 ! Set flag that the trajectory was set
 ! ------------------------------------
 self%ltraj = .true.
 call f_comm%final()

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
integer :: jprofile, jchannel, jspec, nval
type(ufo_geoval), pointer :: geoval_d
real(kind_real) :: geoval_unit_rescale

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

 ! Initialize hofx
 ! ---------------
 hofx(:,:) = 0.0_kind_real

 ! In every block below: hofx = hofx + K . dx (Jacobian times increment), per channel and profile

 ! Temperature
 ! -----------
 call ufo_geovals_get_var(geovals, var_ts, geoval_d)

 ! Check model levels is consistent in geovals & crtm
 if (geoval_d%nval /= self%n_Layers) then
   err_msg = myname_ // " error: layers inconsistent!"
   call abor1_ftn(err_msg)
 end if

 nval = geoval_d%nval
 do jprofile = 1, self%n_Profiles
   if (self%Options(jprofile)%Skip_Profile) cycle
   do jchannel = 1, size(self%channels)
     hofx(jchannel, jprofile) = hofx(jchannel, jprofile) + &
         dot_product(self%K_temp(1:nval,jchannel,jprofile), geoval_d%vals(:,jprofile))
   end do
 end do

 ! Absorbers
 ! ---------
 do jspec = 1, self%conf%n_Absorbers
   call ufo_geovals_get_var(geovals, self%conf%Absorbers(jspec), geoval_d)
   geoval_unit_rescale = absorber_unit_rescale(self%conf%Absorbers(jspec))

   nval = geoval_d%nval
   do jprofile = 1, self%n_Profiles
     if (self%Options(jprofile)%Skip_Profile) cycle
     do jchannel = 1, size(self%channels)
       hofx(jchannel, jprofile) = hofx(jchannel, jprofile) + geoval_unit_rescale * &
           dot_product(self%K_abs(1:nval,jchannel,jprofile,jspec), &
                       geoval_d%vals(:,jprofile))
     end do
   end do
 end do

 ! Clouds (mass content only)
 ! --------------------------
 do jspec = 1, self%conf%n_Clouds
   call ufo_geovals_get_var(geovals, self%conf%Clouds(jspec,1), geoval_d)

   nval = geoval_d%nval
   do jprofile = 1, self%n_Profiles
     if (self%Options(jprofile)%Skip_Profile) cycle
     do jchannel = 1, size(self%channels)
       hofx(jchannel, jprofile) = hofx(jchannel, jprofile) + &
           dot_product(self%K_cld(1:nval,jchannel,jprofile,jspec), &
                       geoval_d%vals(:,jprofile))
     end do
   end do
 end do

 ! Surface Variables (single level): all channels at once
 ! --------------------------------
 do jspec = 1, self%conf%n_Surfaces
   call ufo_geovals_get_var(geovals, self%conf%Surfaces(jspec), geoval_d)
   do jprofile = 1, self%n_Profiles
     if (self%Options(jprofile)%Skip_Profile) cycle
     hofx(1:size(self%channels), jprofile) = hofx(1:size(self%channels), jprofile) + &
         self%K_sfc(1:size(self%channels), jprofile, jspec) * &
         geoval_d%vals(1, jprofile)
   end do
 end do

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
integer :: jprofile, jchannel, jspec, nval, nchan
type(ufo_geoval), pointer :: geoval_d
real(c_double) :: missing
real(kind_real) :: geoval_unit_rescale

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
 nchan = size(self%channels)

 ! In every block below: dx = dx + K^T hofx (adjoint), skipping missing hofx values

 ! Temperature
 ! -----------
 call ufo_geovals_get_var(geovals, var_ts, geoval_d)

 nval = geoval_d%nval
 do jprofile = 1, self%n_Profiles
   if (self%Options(jprofile)%Skip_Profile) cycle
   do jchannel = 1, nchan
     if (hofx(jchannel, jprofile) == missing) cycle
     geoval_d%vals(:,jprofile) = geoval_d%vals(:,jprofile) + &
         self%K_temp(1:nval,jchannel,jprofile) * hofx(jchannel, jprofile)
   end do
 end do

 ! Absorbers
 ! ---------
 do jspec = 1, self%conf%n_Absorbers
   call ufo_geovals_get_var(geovals, self%conf%Absorbers(jspec), geoval_d)
   geoval_unit_rescale = absorber_unit_rescale(self%conf%Absorbers(jspec))

   nval = geoval_d%nval
   do jprofile = 1, self%n_Profiles
     if (self%Options(jprofile)%Skip_Profile) cycle
     do jchannel = 1, nchan
       if (hofx(jchannel, jprofile) == missing) cycle
       geoval_d%vals(:,jprofile) = geoval_d%vals(:,jprofile) + &
           self%K_abs(1:nval,jchannel,jprofile,jspec) * &
           geoval_unit_rescale * hofx(jchannel, jprofile)
     end do
   end do
 end do

 ! Clouds (mass content only)
 ! --------------------------
 do jspec = 1, self%conf%n_Clouds
   call ufo_geovals_get_var(geovals, self%conf%Clouds(jspec,1), geoval_d)

   nval = geoval_d%nval
   do jprofile = 1, self%n_Profiles
     if (self%Options(jprofile)%Skip_Profile) cycle
     do jchannel = 1, nchan
       if (hofx(jchannel, jprofile) == missing) cycle
       geoval_d%vals(:,jprofile) = geoval_d%vals(:,jprofile) + &
           self%K_cld(1:nval,jchannel,jprofile,jspec) * hofx(jchannel, jprofile)
     end do
   end do
 end do

 ! Surface Variables (single level): sum over all non-missing channels at once
 ! --------------------------------
 do jspec = 1, self%conf%n_Surfaces
   call ufo_geovals_get_var(geovals, self%conf%Surfaces(jspec), geoval_d)
   do jprofile = 1, self%n_Profiles
     if (self%Options(jprofile)%Skip_Profile) cycle
     geoval_d%vals(1, jprofile) = geoval_d%vals(1, jprofile) + &
         sum(self%K_sfc(1:nchan, jprofile, jspec) * hofx(1:nchan, jprofile), &
             mask = hofx(1:nchan, jprofile) /= missing)
   end do
 end do

end subroutine ufo_radiancecrtm_simobs_ad

! ------------------------------------------------------------------------------

!> Jacobian of every channel of one profile w.r.t. the surface GeoVaL `varname`.
!> Surfaces are validated against UFO_Surfaces in ufo_crtm_conf_setup.
!> Copy the Jacobians simobs_tl/ad use out of the CRTM K-matrix structures into
!> contiguous arrays (levels first), in the order of the linear operator's conf.
subroutine pack_jacobians(self, atm_K, sfc_K)
  type(ufo_radiancecrtm_tlad), intent(inout) :: self
  type(CRTM_Atmosphere_type),  intent(in)    :: atm_K(:,:)   !< (n_Channels, n_Profiles)
  type(CRTM_Surface_type),     intent(in)    :: sfc_K(:,:)   !< (n_Channels, n_Profiles)

  integer :: jprofile, jchannel, jspec, ispec, nlay, nch, nprof

  nlay  = self%n_Layers
  nch   = self%n_Channels
  nprof = self%n_Profiles

  call free_packed_jacobians(self)
  allocate(self%K_temp(nlay, nch, nprof),                        &
           self%K_abs (nlay, nch, nprof, self%conf%n_Absorbers), &
           self%K_cld (nlay, nch, nprof, self%conf%n_Clouds),    &
           self%K_sfc (nch, nprof, self%conf%n_Surfaces))

  do jprofile = 1, nprof
    do jchannel = 1, nch
      self%K_temp(:,jchannel,jprofile) = atm_K(jchannel,jprofile)%Temperature(1:nlay)
    end do
  end do

  do jspec = 1, self%conf%n_Absorbers
    ispec = ufo_vars_getindex(self%conf_traj%Absorbers, self%conf%Absorbers(jspec))
    do jprofile = 1, nprof
      do jchannel = 1, nch
        self%K_abs(:,jchannel,jprofile,jspec) = atm_K(jchannel,jprofile)%Absorber(1:nlay,ispec)
      end do
    end do
  end do

  do jspec = 1, self%conf%n_Clouds
    ispec = ufo_vars_getindex(self%conf_traj%Clouds(:,1), self%conf%Clouds(jspec,1))
    do jprofile = 1, nprof
      do jchannel = 1, nch
        self%K_cld(:,jchannel,jprofile,jspec) = &
            atm_K(jchannel,jprofile)%Cloud(ispec)%Water_Content(1:nlay)
      end do
    end do
  end do

  do jspec = 1, self%conf%n_Surfaces
    do jprofile = 1, nprof
      self%K_sfc(:,jprofile,jspec) = sfc_jacobian(sfc_K(1:nch,jprofile), self%conf%Surfaces(jspec))
    end do
  end do

end subroutine pack_jacobians

! ------------------------------------------------------------------------------

subroutine free_packed_jacobians(self)
  type(ufo_radiancecrtm_tlad), intent(inout) :: self

  if (allocated(self%K_temp)) deallocate(self%K_temp)
  if (allocated(self%K_abs))  deallocate(self%K_abs)
  if (allocated(self%K_cld))  deallocate(self%K_cld)
  if (allocated(self%K_sfc))  deallocate(self%K_sfc)

end subroutine free_packed_jacobians

! ------------------------------------------------------------------------------

function sfc_jacobian(sfc_K, varname) result(jac)
  type(CRTM_Surface_type), intent(in) :: sfc_K(:)   !< one per channel
  character(len=*),        intent(in) :: varname
  real(kind_real) :: jac(size(sfc_K))

  select case (varname)
    case (var_sfc_wtmp);   jac = sfc_K(:)%Water_Temperature
    case (var_sfc_ltmp);   jac = sfc_K(:)%Land_Temperature
    case (var_sfc_itmp);   jac = sfc_K(:)%Ice_Temperature
    case (var_sfc_stmp);   jac = sfc_K(:)%Snow_Temperature
    case (var_sfc_wspeed); jac = sfc_K(:)%Wind_Speed
    case (var_sfc_wdir);   jac = sfc_K(:)%Wind_Direction
    case (var_sfc_sss);    jac = sfc_K(:)%Salinity
    case default
      call abor1_ftn("ufo_radiancecrtm_tlad: unsupported surface variable "//trim(varname))
  end select
end function sfc_jacobian

! ------------------------------------------------------------------------------

!> Factor converting the absorber GeoVaL to CRTM units:
!> "water_vapor_mixing_ratio_wrt_dry_air" is kg/kg in JEDI but g/kg in CRTM.
real(kind_real) function absorber_unit_rescale(varname)
  character(len=*), intent(in) :: varname
  absorber_unit_rescale = one
  if (varname == var_mixr) absorber_unit_rescale = kg_to_g
end function absorber_unit_rescale

! ------------------------------------------------------------------------------

end module ufo_radiancecrtm_tlad_mod
