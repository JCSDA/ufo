! (C) Copyright 2017-2018 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

!> Fortran module to handle radiancecrtm observations

module ufo_radiancecrtm_mod
 use,intrinsic :: iso_c_binding
 use,intrinsic :: iso_fortran_env
 use CRTM_Module

 use fckit_configuration_module, only: fckit_configuration
 use fckit_mpi_module,   only: fckit_mpi_comm
 use, intrinsic :: iso_c_binding
 use kinds
 use missing_values_mod
 use ufo_reconradop_crtm_mod
 use obsspace_mod
 use obsdatavector_mod
 use ufo_geovals_mod, only: ufo_geovals, ufo_geoval, ufo_geovals_get_var, ufo_geovals_nlevs
 use ufo_vars_mod
 use ufo_crtm_utils_mod
 use ufo_crtm_passive_mod
 use ufo_crtm_active_mod

 use ufo_constants_mod, only: deg2rad

 implicit none
 private

 !> Fortran derived type for radiancecrtm trajectory
 type, public :: ufo_radiancecrtm
 private
   character(len=MAXVARLEN), public, allocatable :: varin(:)  ! variables requested from the model
   integer, allocatable                          :: channels(:)
   type(crtm_conf) :: conf
   logical :: use_qc_flags
   type(ufo_reconradop_crtm) :: reconradop_crtm
 contains
   procedure :: setup  => ufo_radiancecrtm_setup
   procedure :: delete => ufo_radiancecrtm_delete
   procedure :: simobs => ufo_radiancecrtm_simobs
 end type ufo_radiancecrtm

 character(len=maxvarlen), dimension(16), parameter :: varin_default = &
                            [var_ts, var_prs, var_prsi,                                  &
                              var_sfc_wfrac, var_sfc_lfrac, var_sfc_ifrac, var_sfc_sfrac, &
                              var_sfc_wtmp,  var_sfc_ltmp,  var_sfc_itmp,  var_sfc_stmp,  &
                              var_sfc_vegfrac, var_sfc_lai,                               &
                              var_sfc_soilm, var_sfc_soilt, var_sfc_sdepth]

contains

! ------------------------------------------------------------------------------

subroutine ufo_radiancecrtm_setup(self, f_confOper, channels, midPointJulday, comm)

implicit none
class(ufo_radiancecrtm),   intent(inout) :: self
type(fckit_configuration), intent(in)    :: f_confOper
integer(c_int),            intent(in)    :: channels(:)  !List of channels to use
integer(c_int64_t),        intent(in)    :: midPointJulday
type(fckit_mpi_comm),      intent(in)    :: comm
integer :: js, jspec
integer :: sensor_type
character(len=max_string) :: err_msg
type(fckit_configuration) :: f_confOpts
character(len=maxvarlen), dimension(10), parameter :: more_hydrometeors = &
                            [var_qsat, var_airdens, var_q,  var_ni, var_nr, &
                              var_qc,   var_qi,      var_qr, var_qs, var_qg]

 ! 1. Read the obs operator configuration
 call f_confOper%get_or_die("obs options",f_confOpts)
 call f_confOper%get_or_die("UseQCFlagsToSkipHofX",self%use_qc_flags)
 call ufo_crtm_conf_setup(self%conf,f_confOpts,f_confOper,midPointJulday,comm)
 call ufo_crtm_check_absorbers(self%conf, .false., "ufo_radiancecrtm_setup")

 ! 2. Variables requested from the model (GeoVaLs)

 ! 2a. Atmospheric and surface state
 self%varin = varin_default

 ! 2b. Surface classification, by sensor type: vegetation & soil type for MW, land type for
 !     VIS/IR. The sensor type is read from the sensor's SpcCoeff file into a local
 !     structure, so CRTM's loaded coefficients (shared between operators) are untouched.
 do js = 1, self%conf%n_Sensors
   sensor_type = ufo_crtm_sensor_type(self%conf, js)
   select case (sensor_type)
     case (MICROWAVE_SENSOR)
       call add_varin([var_sfc_vegtyp, var_sfc_soiltyp])
     case (VISIBLE_SENSOR, INFRARED_SENSOR)
       call add_varin([ufo_crtm_landtyp_var(self%conf, "ufo_radiancecrtm_setup")])
     case default
       write(err_msg,"(a,i0,2a)") "ufo_radiancecrtm_setup error: unsupported Sensor_Type = ", &
                                  sensor_type, " from Sensor_ID = ", trim(self%conf%SENSOR_ID(js))
       call abor1_ftn(err_msg)
   end select
 end do

 ! 2c. Absorbers, and the mass content & effective radius of each cloud species
 call add_varin(self%conf%Absorbers)
 do jspec = 1, self%conf%n_Clouds
   call add_varin(self%conf%Clouds(jspec,1:2))
 end do

 ! 2d. Cloud fraction, unless it is a constant in the conf or computed in the FOV
 if (self%conf%n_Clouds > 0 .and. self%conf%Cloud_Fraction < 0.0 .and. &
     .not. self%conf%cal_cloud_frac_in_fov) then
   call add_varin([var_cldfrac])
   if (comm%rank() == 0) call Display_Message("ufo_radiancecrtm_setup", &
     "Cloud_Fraction is not provided in conf. Will request as a geoval.", WARNING)
 end if

 ! 2e. Sea surface salinity
 if (self%conf%salinity_option == "on") then
   call add_varin([var_sfc_sss])
 end if

 ! 2f. Surface wind: speed & direction, or u & v
 if (self%conf%sfc_wind_geovars == "vector") then
   call add_varin([var_sfc_wspeed, var_sfc_wdir])
 else if (self%conf%sfc_wind_geovars == "uv") then
   call add_varin([var_sfc_u, var_sfc_v])
 end if

 ! 2g. Hydrometeors for the in-FOV (Thompson) cloud fraction or effective radius
 if (self%conf%cal_cloud_frac_in_fov .or. self%conf%cal_cloud_reff_in_fov) then
   call add_varin(more_hydrometeors)
 end if

 ! 3. Channels, and registration for the shared CRTM coefficient load used in simobs
 self%channels = channels
 call f_confOpts%final()
 call ufo_crtm_register_radiance(self%conf)

contains

  !> Append the names not yet in the list of variables requested from the model
  subroutine add_varin(names)
    character(len=*), intent(in) :: names(:)
    integer :: i
    do i = 1, size(names)
      if (ufo_vars_getindex(self%varin, names(i)) < 1) then
        self%varin = [character(len=maxvarlen) :: self%varin, names(i)]
      end if
    end do
  end subroutine add_varin

end subroutine ufo_radiancecrtm_setup

! ------------------------------------------------------------------------------

subroutine ufo_radiancecrtm_delete(self)

implicit none
class(ufo_radiancecrtm), intent(inout) :: self

 call ufo_crtm_deregister_radiance(self%conf)
 call ufo_crtm_conf_delete(self%conf)

end subroutine ufo_radiancecrtm_delete

! ------------------------------------------------------------------------------

subroutine ufo_radiancecrtm_simobs(self, geovals, obss, nvars, nlocs, hofx, hofxdiags, qcf_p)
use fckit_mpi_module,   only: fckit_mpi_comm
use obsdatavector_mod,  only: obsdatavector_int
use CRTM_SpcCoeff, only: SC, &
                         SpcCoeff_IsMicrowaveSensor , &
                         SpcCoeff_IsInfraredSensor  , &
                         SpcCoeff_IsVisibleSensor   , &
                         SpcCoeff_IsUltravioletSensor

implicit none

class(ufo_radiancecrtm),  intent(in) :: self      !Radiance object
type(ufo_geovals),        intent(in) :: geovals      !Inputs from the model
integer(c_size_t),        intent(in) :: nvars, nlocs
real(c_double),        intent(inout) :: hofx(nvars, nlocs) !h(x) to return
type(ufo_geovals),     intent(inout) :: hofxdiags    !non-h(x) diagnostics
type(c_ptr), value,       intent(in) :: obss         !ObsSpace
type(c_ptr), value,       intent(in) :: qcf_p
type(obsdatavector_int) :: qc_flags
! Local Variables
character(len=*), parameter :: PROGRAM_NAME = "ufo_radiancecrtm_simobs"
character(255) :: message
character(max_string) :: err_msg
integer        :: err_stat
integer        :: isc  ! this sensor's entry in CRTM's SpcCoeff table (SC)
integer        :: n
integer :: jvar, jprofile, jlevel, jchannel
type(fckit_mpi_comm)  :: f_comm

integer :: n_Profiles, n_Layers, n_Channels
! Define the Channel Info and Geometry arguments
type(CRTM_ChannelInfo_type)             :: chinfo(self%conf%n_Sensors)
type(CRTM_Geometry_type),   allocatable :: geo(:)

! Define the FORWARD variables
type(CRTM_Atmosphere_type), allocatable :: atm(:)
type(CRTM_Surface_type),    allocatable :: sfc(:)
type(CRTM_RTSolution_type), allocatable :: rts(:,:)

! Define the K-MATRIX variables for hofxdiags
type(CRTM_Atmosphere_type), allocatable :: atm_K(:,:)
type(CRTM_Surface_type),    allocatable :: sfc_K(:,:)
type(CRTM_RTSolution_type), allocatable :: rts_K(:,:)
type(CRTM_Options_type),    allocatable :: Options(:)

!for gmi
type(CRTM_Geometry_type),   allocatable :: geo_hf(:)
type(CRTM_ChannelInfo_type) :: chinfo_gmi(2)  ! gmi channels 1-9 and 10-13
integer :: n_lo                               ! number of gmi channels 1-9

! Used to parse hofxdiags
character(len=MAXVARLEN) :: varstr
character(len=MAXVARLEN), dimension(hofxdiags%nvar) :: &
                          ystr_diags, xstr_diags
character(10), parameter :: jacobianstr = "_jacobian_"
integer :: str_pos(4), ch_diags(hofxdiags%nvar)
logical :: jacobian_needed

! set a local boolean variable for whether we are in vis or ultraviolet channels
logical        :: Is_Vis_or_UV

integer, allocatable :: zeroCloudInCRTM0(:)

 call obsspace_get_comm(obss, f_comm)

 ! Get number of profile and layers from geovals
 ! ---------------------------------------------
 n_Profiles = geovals%nlocs
 n_Layers = ufo_geovals_nlevs(geovals, var_ts)

 qc_flags%data_ptr = qcf_p
 ! ChannelInfo for this operator's sensors, from the shared CRTM coefficient
 ! load (CRTM_Init runs only on first use; see ufo_crtm_utils_mod)
 ! ---------------------------------------------------------------------------

 call ufo_crtm_init_radiance(self%conf, chinfo, self%channels, PROGRAM_NAME, f_comm)

 !! Parse hofxdiags%variables into independent/dependent variables and channel
 !! assumed formats:
 !!   jacobian var -->     <ystr>_jacobian_<xstr>_<chstr>
 !!   non-jacobian var --> <ystr>_<chstr>

 jacobian_needed = .false.
 ch_diags = -9999
 do jvar = 1, hofxdiags%nvar
   varstr = hofxdiags%variables(jvar)
   str_pos(4) = len_trim(varstr)
   if (str_pos(4) < 1) cycle
   str_pos(3) = index(varstr,"_",back=.true.)        !final "_" before channel
   read(varstr(str_pos(3)+1:str_pos(4)),*, err=999) ch_diags(jvar)
999 str_pos(1) = index(varstr,jacobianstr) - 1        !position before jacobianstr
   if (str_pos(1) == 0) then
     err_msg = "ufo_radiancecrtm_simobs: _jacobian_ must be preceded by dependent " // &
               "variable in config: " // hofxdiags%variables(jvar)
     call abor1_ftn(err_msg)
   else if (str_pos(1) > 0) then
     !Diagnostic is a Jacobian member (dy/dx)
     ystr_diags(jvar) = varstr(1:str_pos(1))
     str_pos(2) = str_pos(1) + len(jacobianstr) + 1 !begin xstr_diags
     jacobian_needed = .true.
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

 ! Loop over sensors (one per operator in practice)
 ! ----------------------------------------------------------------------------
 Sensor_Loop:do n = 1, self%conf%n_Sensors


   ! Determine the number of channels for the current sensor
   ! -------------------------------------------------------
   n_Channels = CRTM_ChannelInfo_n_Channels(chinfo(n))
   isc = chinfo(n)%Sensor_Index

   ! Allocate the ARRAYS (for CRTM_Forward)
   ! --------------------------------------
   allocate(geo(n_Profiles))
   allocate(atm(n_Profiles))
   allocate(sfc(n_Profiles))
   allocate(rts(n_Channels, n_Profiles))
   allocate(Options(n_Profiles))

   ! Create the input FORWARD structure (atm)
   ! ----------------------------------------
   call CRTM_Atmosphere_Create( atm, n_Layers, self%conf%n_Absorbers, self%conf%n_Clouds, &
                               self%conf%n_Aerosols )
   if ( any(.not. CRTM_Atmosphere_Associated(atm)) ) then
      message = "Error allocating CRTM Forward Atmosphere structure"
      call ufo_crtm_comm_stat_check(FAILURE, PROGRAM_NAME, message, f_comm)
   end if


   ! Create the input FORWARD structure (sfc)
   ! ----------------------------------------
   call CRTM_Surface_Create(sfc, n_Channels)
   if ( any(.not. CRTM_Surface_Associated(sfc)) ) then
      message = "Error allocating CRTM Surface structure"
      call ufo_crtm_comm_stat_check(FAILURE, PROGRAM_NAME, message, f_comm)
   end if

   call CRTM_RTSolution_Create(rts, n_Layers )

   ! Some special treatment for Vis or UV as compared to IR or MW.
   ! ----------------------------------------
   if (SpcCoeff_IsVisibleSensor(SC(isc)) .or. SpcCoeff_IsUltravioletSensor(SC(isc))) then
      Is_Vis_or_UV = .true.
   else
      Is_Vis_or_UV = .false.
   end if

   !Assign the data from the GeoVaLs
   !--------------------------------
   allocate(zeroCloudInCRTM0(n_Profiles))
   call ufo_crtm_load_sfc_data(n_Profiles, n_Channels, self%channels, geovals, sfc, chinfo, obss, &
                               self%conf, SC(isc)%Is_Active_Sensor, Is_Vis_or_UV, zeroCloudInCRTM0)

   call ufo_crtm_load_atm_data(n_Profiles, n_Layers, geovals, atm, self%conf, &
                               SC(isc)%Is_Active_Sensor, zeroCloudInCRTM0)
   deallocate(zeroCloudInCRTM0)
   if (self%conf%SENSOR_ID(n) == "gmi_gpm") then
     allocate( geo_hf( n_Profiles ))
     call ufo_crtm_load_geom_data(obss,geo,geo_hf,self%conf%SENSOR_ID(n))
   else
     call ufo_crtm_load_geom_data(obss,geo)
   end if

   ! Call THE CRTM inspection
   ! ------------------------
   if (self%conf%inspect > 0) then
     call CRTM_Atmosphere_Inspect(atm(self%conf%inspect))
     call CRTM_Surface_Inspect(sfc(self%conf%inspect))
     call CRTM_Geometry_Inspect(geo(self%conf%inspect))
     call CRTM_ChannelInfo_Inspect(chinfo(n))
   end if

   ! set profiles that should be skipped
   call ufo_crtm_skip_profiles(n_Profiles,n_Channels,self%channels,obss,atm,sfc,  &
                               SC(isc)%Is_Active_Sensor,Is_Vis_or_UV,Options)
   ! eliminate remaining profiles that are "QCed" out
   if ( self%use_qc_flags ) then
     call ufo_crtm_skip_qc_profiles(qc_flags, size(self%channels), Options(1:n_Profiles), "DEBUG")
   end if

   if (jacobian_needed) then
      ! Allocate the ARRAYS (for CRTM_K_Matrix)
      ! --------------------------------------
      allocate(atm_K(n_Channels, n_Profiles))
      allocate(sfc_K(n_Channels, n_Profiles))
      allocate(rts_K(n_Channels, n_Profiles))

      ! Create output K-MATRIX structure (atm)
      ! --------------------------------------
      call CRTM_Atmosphere_Create( atm_K, n_Layers, self%conf%n_Absorbers, self%conf%n_Clouds, &
                                  self%conf%n_Aerosols )
      if ( any(.not. CRTM_Atmosphere_Associated(atm_K)) ) then
         message = "Error allocating CRTM K-matrix Atmosphere structure (setTraj)"
         call ufo_crtm_comm_stat_check(FAILURE, PROGRAM_NAME, message, f_comm)
      end if

      ! Create output K-MATRIX structure (sfc)
      ! --------------------------------------
      ! No SensorData: CRTM reads it only from the forward sfc, never from sfc_K
      call CRTM_Surface_Create( sfc_K, 0)

      ! Zero the K-matrix OUTPUT structures
      ! (atm_K is already zeroed by CRTM_Atmosphere_Create)
      ! -----------------------------------
      call CRTM_Surface_Zero( sfc_K )

      ! Initialize the K-matrix INPUT so that the results are dTb/dx or dR/dx
      ! -------------------------------------------------------------
      if (SC(isc)%Is_Active_Sensor) then
         do jchannel = 1, n_Channels
            do jprofile = 1, n_Profiles
               do jlevel = 1, n_Layers
                  rts_K(jchannel,jprofile)%Reflectivity(jlevel)            = ZERO
                  rts_K(jchannel,jprofile)%Reflectivity_Attenuated(jlevel) = ONE
               end do
            end do
         end do
         rts_K%Radiance                = ZERO
         rts_K%Brightness_Temperature  = ZERO
      else if (Is_Vis_or_UV .or. self % conf % read_Cmatrix) then
         rts_K%Radiance                = ONE
         rts_K%Brightness_Temperature  = ZERO
      else
         rts_K%Radiance                = ZERO
         rts_K%Brightness_Temperature  = ONE
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
                                      Options )
            message = "Error calling CRTM (setTraj) K-Matrix Model for gmi_gpm channels 1-9"
            call ufo_crtm_comm_stat_check(err_stat, PROGRAM_NAME, message, f_comm)
         end if
         if (n_lo < n_Channels) then
            ! K-matrix input so that the results are dTb/dx
            rts_K(n_lo+1:,:)%Radiance               = ZERO
            rts_K(n_lo+1:,:)%Brightness_Temperature = ONE
            err_stat = CRTM_K_Matrix( atm, sfc, rts_K(n_lo+1:,:), geo_hf, chinfo_gmi(2:2), &
                                      atm_K(n_lo+1:,:), sfc_K(n_lo+1:,:), rts(n_lo+1:,:), &
                                      Options )
            message = "Error calling CRTM (setTraj, geo_hf) K-Matrix Model for " // &
                      "gmi_gpm channels 10-13"
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
                                   Options       )  ! Input
         message = "Error calling CRTM (setTraj) K-Matrix Model for "//self%conf%SENSOR_ID(n)
         call ufo_crtm_comm_stat_check(err_stat, PROGRAM_NAME, message, f_comm)
      end if
   else
      ! Call the forward model call for each sensor
      ! -------------------------------------------
      if (self%conf%SENSOR_ID(n) == "gmi_gpm") then
         ! channels 1-9 with geo, channels 10-13 with geo_hf: one call per block
         call ufo_crtm_gmi_split_chinfo(chinfo(n), self%channels, chinfo_gmi, n_lo, &
                                        PROGRAM_NAME, f_comm)
         if (n_lo > 0) then
            err_stat = CRTM_Forward( atm, sfc, geo, chinfo_gmi(1:1), rts(1:n_lo,:), Options )
            message = "Error calling CRTM Forward Model for gmi_gpm channels 1-9"
            call ufo_crtm_comm_stat_check(err_stat, PROGRAM_NAME, message, f_comm)
         end if
         if (n_lo < n_Channels) then
            err_stat = CRTM_Forward( atm, sfc, geo_hf, chinfo_gmi(2:2), rts(n_lo+1:,:), Options )
            message = "Error calling CRTM Forward Model for gmi_gpm channels 10-13"
            call ufo_crtm_comm_stat_check(err_stat, PROGRAM_NAME, message, f_comm)
         end if
         call CRTM_ChannelInfo_Destroy(chinfo_gmi)
      else
         err_stat = CRTM_Forward( atm         , &  ! Input
                                  sfc         , &  ! Input
                                  geo         , &  ! Input
                                  chinfo(n:n) , &  ! Input
                                  rts         , &  ! Output
                                  Options       )  ! Input
         message = "Error calling CRTM Forward Model for "//self%conf%SENSOR_ID(n)
         call ufo_crtm_comm_stat_check(err_stat, PROGRAM_NAME, message, f_comm)
      end if
   end if ! jacobian_needed

   ! put CRTM output into hofxdiags
   if (SC(isc)%Is_Active_Sensor) then
      call ufo_crtm_active_sim(rts, &
                               Options, &
                               nvars, &
                               nlocs, &
                               n_Profiles, &
                               n_Channels, &
                               hofx, &
                               obss)

      call ufo_crtm_active_diag(rts, &
                                rts_K, &
                                atm, &
                                atm_K, &
                                sfc_K, &
                                self%conf, &
                                n, &
                                Options, &
                                self%channels, &
                                geovals, &
                                obss, &
                                nvars, &
                                nlocs, &
                                n_Profiles, &
                                n_Layers, &
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
                                         n_Profiles, &
                                         n_Channels, &
                                         self%conf%n_Absorbers,&
                                         n_Layers,&
                                         self%channels, &
                                         rts_K, &
                                         atm_K, &
                                         sfc_K, jacobian_needed )
      end if

     call ufo_crtm_passive_sim(rts, &
                               Options, &
                               Is_Vis_or_UV, &
                               nvars, &
                               nlocs, &
                               n_Profiles, &
                               n_Channels, &
                               hofx)

      call ufo_crtm_passive_diag(rts, &
                                 rts_K, &
                                 atm, &
                                 atm_K, &
                                 sfc_K, &
                                 self%conf, &
                                 n, &
                                 Options, &
                                 self%channels, &
                                 geovals, &
                                 obss, &
                                 nvars, &
                                 nlocs, &
                                 n_Profiles, &
                                 n_Layers, &
                                 xstr_diags, &
                                 ystr_diags, &
                                 ch_diags, &
                                 hofxdiags,&
                                 err_stat)
   end if

   ! check for error from either passive or active
   if (err_stat > 0) then
       err_msg = "ufo_radiancecrtm_simobs error: " // &
                 "failed to put simulated diagnostics into hofxdiags"
       call abor1_ftn(err_msg)
   end if

   ! Deallocate the structures
   ! -------------------------
   call CRTM_Geometry_Destroy(geo)
   call CRTM_Atmosphere_Destroy(atm)
   call CRTM_Surface_Destroy(sfc)
   call CRTM_RTSolution_Destroy(rts)

   ! Deallocate all arrays
   ! ---------------------
   deallocate(geo)
   deallocate(atm)
   deallocate(sfc)
   deallocate(rts)
   deallocate(Options)
   if(allocated(geo_hf)) deallocate(geo_hf)

   if (jacobian_needed) then
      ! Deallocate the K structures
      ! ---------------------------
      call CRTM_Atmosphere_Destroy(atm_K)
      call CRTM_Surface_Destroy(sfc_K)
      call CRTM_RTSolution_Destroy(rts_K)

      ! Deallocate all K arrays
      ! -----------------------
      deallocate(atm_K)
      deallocate(sfc_K)
      deallocate(rts_K)
   end if

 end do Sensor_Loop


 ! Destroy CRTM instance
 ! ---------------------
 call ufo_crtm_destroy(chinfo, PROGRAM_NAME, f_comm)
 call f_comm%final()

end subroutine ufo_radiancecrtm_simobs

! ------------------------------------------------------------------------------

end module ufo_radiancecrtm_mod
