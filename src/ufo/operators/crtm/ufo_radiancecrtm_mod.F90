! (C) Copyright 2017-2018 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

!> Fortran module to handle radiancecrtm observations

module ufo_radiancecrtm_mod
 use,intrinsic :: iso_c_binding
 use,intrinsic :: iso_fortran_env
 use crtm_module

 use fckit_configuration_module, only: fckit_configuration
 use fckit_mpi_module,   only: fckit_mpi_comm
 use iso_c_binding
 use kinds
 use missing_values_mod

 use obsspace_mod
 use obsdatavector_mod
 use ufo_geovals_mod, only: ufo_geovals, ufo_geoval, ufo_geovals_get_var
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
 contains
   procedure :: setup  => ufo_radiancecrtm_setup
   procedure :: delete => ufo_radiancecrtm_delete
   procedure :: simobs => ufo_radiancecrtm_simobs
 end type ufo_radiancecrtm

 character(len=maxvarlen), dimension(16), parameter :: varin_default = &
                            (/var_ts, var_prs, var_prsi,                                  &
                              var_sfc_wfrac, var_sfc_lfrac, var_sfc_ifrac, var_sfc_sfrac, &
                              var_sfc_wtmp,  var_sfc_ltmp,  var_sfc_itmp,  var_sfc_stmp,  &
                              var_sfc_vegfrac, var_sfc_lai,                               &
                              var_sfc_soilm, var_sfc_soilt, var_sfc_sdepth/)

contains

! ------------------------------------------------------------------------------

subroutine ufo_radiancecrtm_setup(self, f_confOper, channels, comm)
use ufo_utils_mod, only: cmp_strings
use CRTM_SpcCoeff, only: CRTM_SpcCoeff_Load, SC

implicit none
class(ufo_radiancecrtm),   intent(inout) :: self
type(fckit_configuration), intent(in)    :: f_confOper
integer(c_int),            intent(in)    :: channels(:)  !List of channels to use
type(fckit_mpi_comm),      intent(in)    :: comm
integer :: nvars_in
integer :: ind, js, jspec
integer :: err_stat
character(len=max_string) :: err_msg
type(fckit_configuration) :: f_confOpts
logical :: request_mw_vegtyp_soiltyp_data, request_visir_landtyp_data
logical :: request_cldfrac, request_salinity


 call f_confOper%get_or_die("obs options",f_confOpts)
 call f_confOper%get_or_die("UseQCFlagsToSkipHofX",self%use_qc_flags)  
 call crtm_conf_setup(self%conf,f_confOpts,f_confOper,comm)
 if ( ufo_vars_getindex(self%conf%Absorbers, var_mixr) < 1 ) then
   write(err_msg,*) 'ufo_radiancecrtm_setup error: H2O must be included in CRTM Absorbers'
   call abor1_ftn(err_msg)
 end if
 if ( ufo_vars_getindex(self%conf%Absorbers, var_oz) < 1 ) then
   write(err_msg,*) 'ufo_radiancecrtm_setup error: O3 must be included in CRTM Absorbers'
   call abor1_ftn(err_msg)
 end if

 ! Check which surface type data CRTM will need to request from the model
 request_mw_vegtyp_soiltyp_data = .false.
 request_visir_landtyp_data = .false.
 ! Use CRTM_SpcCoeff_Load to fill the CRTM "SC" structure with data
 ! Perhaps in the future a simpler lookup table could be used, but this would require changes
 ! to the CRTM interface to expose the Sensor_ID -> Sensor_Type mapping in an easier way.
 err_stat = CRTM_SpcCoeff_Load(Sensor_ID = self%conf%Sensor_ID, &
                               File_Path = trim(self%conf%COEFFICIENT_PATH), &
                               Quiet     = .true.)
 if (err_stat /= success) then
   write(err_msg,*) 'ufo_radiancecrtm_setup error: failed CRTM_Load_SpcCoeff'
   call abor1_ftn(err_msg)
 end if
 do js=1, self%conf%n_Sensors
   if (SC(js)%Sensor_Type == MICROWAVE_SENSOR) then
     request_mw_vegtyp_soiltyp_data = .true.
   else if (SC(js)%Sensor_Type == VISIBLE_SENSOR .or. SC(js)%Sensor_Type == INFRARED_SENSOR) then
     request_visir_landtyp_data = .true.
   else
     write(err_msg,*) 'ufo_radiancecrtm_setup error: unsupported Sensor_Type =', &
                      SC(js)%Sensor_Type, ', from Sensor_ID =', SC(js)%Sensor_ID
     call abor1_ftn(err_msg)
   end if
 end do
 deallocate(SC)
 ! check that at least one surface type is requested
 if (.not. request_mw_vegtyp_soiltyp_data .and. .not. request_visir_landtyp_data) then
   write(err_msg,*) 'ufo_radiancecrtm_setup error: did not request surface data for Sensor_ID =', &
                    SC(js)%Sensor_ID
   call abor1_ftn(err_msg)
 end if

 request_cldfrac = self%conf%n_Clouds > 0 .and. self%conf%Cloud_Fraction < 0.0

 request_salinity = cmp_strings(self%conf%salinity_option, "on")

 ! request from the model all the hardcoded atmospheric & surface variables +
 ! 1-3 surface classification variables (depending on the instrument wavelength)
 ! 1 * n_Absorbers
 ! 2 * n_Clouds (mass content and effective radius)
 ! 2 for sfc_wind_geovars parsing
 ! 0-1 cloud fraction
 ! 0-1 sea surface salinity
 nvars_in = size(varin_default) + self%conf%n_Absorbers + 2 * self%conf%n_Clouds + 2
 if (request_mw_vegtyp_soiltyp_data) then
   nvars_in = nvars_in + 2
 else if (request_visir_landtyp_data) then
   nvars_in = nvars_in + 1
 end if
 if (request_cldfrac) then
   nvars_in = nvars_in + 1
 end if
 if (request_salinity) then
   nvars_in = nvars_in + 1
 end if

 allocate(self%varin(nvars_in))
 self%varin(1:size(varin_default)) = varin_default

 ! varin_default specifies the NPOESS land_type classification to match the CRTM default.
 ! Here we check if an alternate land_type classification has been requested, and, if so, we
 ! swap out the variable in varin.
 ind = size(varin_default) + 1
 if (request_mw_vegtyp_soiltyp_data) then
   self%varin(ind) = var_sfc_vegtyp
   ind = ind + 1
   self%varin(ind) = var_sfc_soiltyp
   ind = ind + 1
 end if
 if (request_visir_landtyp_data) then
   if (index(self%conf%IRlandCoeff_File, "IGBP") == 1) then
     self%varin(ind) = var_sfc_landtyp_igbp
   else if (index(self%conf%IRlandCoeff_File, "USGS") == 1) then
     self%varin(ind) = var_sfc_landtyp_usgs
   else if (index(self%conf%IRlandCoeff_File, "NPOESS") == 1) then
     self%varin(ind) = var_sfc_landtyp_npoess
   else
     write(err_msg,*) "ufo_radiancecrtm_setup error: cannot infer land type classification from " &
                      // "IRlandCoeff_File: " // trim(self%conf%IRlandCoeff_File)
     call abor1_ftn(err_msg)
   end if
   ind = ind + 1
 endif

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
 if (request_cldfrac) then
   self%varin(ind) = var_cldfrac
   ind = ind + 1
 end if
 if (request_salinity) then
   self%varin(ind) = var_sfc_sss
   ind = ind + 1
 end if

 if (trim(self%conf%sfc_wind_geovars) == 'vector') then
   self%varin(ind) = var_sfc_wspeed
   ind = ind + 1
   self%varin(ind) = var_sfc_wdir
   ind = ind + 1
 else if (trim(self%conf%sfc_wind_geovars) == 'uv') then
   self%varin(ind) = var_sfc_u
   ind = ind + 1
   self%varin(ind) = var_sfc_v
   ind = ind + 1
 end if

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

subroutine ufo_radiancecrtm_simobs(self, geovals, obss, nvars, nlocs, hofx, hofxdiags, qcf_p)
use fckit_mpi_module,   only: fckit_mpi_comm
use ufo_utils_mod,      only: cmp_strings
use CRTM_SpcCoeff, only: SC, &
                         SpcCoeff_IsMicrowaveSensor , &
                         SpcCoeff_IsInfraredSensor  , &
                         SpcCoeff_IsVisibleSensor   , &
                         SpcCoeff_IsUltravioletSensor

implicit none

class(ufo_radiancecrtm),  intent(in) :: self         !Radiance object
type(ufo_geovals),        intent(in) :: geovals      !Inputs from the model
integer(c_size_t),        intent(in) :: nvars, nlocs
real(c_double),        intent(inout) :: hofx(nvars, nlocs) !h(x) to return
type(ufo_geovals),     intent(inout) :: hofxdiags    !non-h(x) diagnostics
type(c_ptr), value,       intent(in) :: obss         !ObsSpace
type(c_ptr), value,       intent(in) :: qcf_p
type(obsdatavector_int) :: qc_flags
! Local Variables
character(*), parameter :: PROGRAM_NAME = 'ufo_radiancecrtm_simobs'
character(255) :: message, version
character(max_string) :: err_msg
integer        :: err_stat, alloc_stat
integer        :: l, m, n
type(ufo_geoval), pointer :: temp
integer :: jvar, jprofile, jlevel, jchannel, ichannel, jspec
real(c_double) :: missing
type(fckit_mpi_comm)  :: f_comm
integer :: n_skipped
integer :: qc_ff
real(kind_real) :: total_od, secant_term, wfunc_max
real(kind_real), allocatable :: TmpVar(:)
real(kind_real), allocatable :: Tao(:)
real(kind_real), allocatable :: Wfunc(:)

integer :: n_Profiles, n_Layers, n_Channels
integer(int64) :: ll,ml
! Define the "non-demoninational" arguments
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
type(CRTM_Atmosphere_type), allocatable :: atm_Ka(:,:)
type(CRTM_Surface_type),    allocatable :: sfc_Ka(:,:)
type(CRTM_RTSolution_type), allocatable :: rts_Ka(:,:)
type(CRTM_RTSolution_type), allocatable :: rtsa(:,:)

! Used to parse hofxdiags
character(len=MAXVARLEN) :: varstr
character(len=MAXVARLEN), dimension(hofxdiags%nvar) :: &
                          ystr_diags, xstr_diags
character(10), parameter :: jacobianstr = "_jacobian_"
integer :: str_pos(4), ch_diags(hofxdiags%nvar)
logical :: jacobian_needed, skip_prof
! For gmi_gpm geophysical angles at channels 10-13.
character(len=1) :: angle_hf

! set a local boolean variable for whether we are in vis or ultraviolet channels
logical        :: Is_Vis_or_UV = .false.

 call obsspace_get_comm(obss, f_comm)

 ! Get number of profile and layers from geovals
 ! ---------------------------------------------
 n_Profiles = geovals%nlocs
 call ufo_geovals_get_var(geovals, var_ts, temp)
 n_Layers = temp%nval
 nullify(temp)

 qc_flags%data_ptr = qcf_p
 ! Program header
 ! --------------
 ! call CRTM_Version( Version )
 ! call Program_Message( PROGRAM_NAME, &
 !                       'UFO interface for the CRTM Forward and K-Matrix functions using '//&
 !                       trim(self%conf%ENDIAN_type)//' coefficient datafiles', &
 !                       'CRTM Version: '//TRIM(Version) )


 ! Initialise all the sensors at once
 ! ----------------------------------
 !** NOTE: CRTM_Init points to the various binary files needed for CRTM.  See the
 !**       CRTM_Lifecycle.f90 for more details.

 ! write( *,'(/5x,"Initializing the CRTM...")' )
 err_stat = CRTM_Init( self%conf%SENSOR_ID                                      , &
                       chinfo                                                   , &
                       File_Path           = trim(self%conf%COEFFICIENT_PATH)   , &
                       NC_File_Path        = trim(self%conf%NC_COEFFICIENT_PATH), &
                       Aerosol_Model       = trim(self%conf%Aerosol_Model)      , &
                       AerosolCoeff_Format = trim(self%conf%AerosolCoeff_Format), &
                       AerosolCoeff_File   = trim(self%conf%AerosolCoeff_File)  , &
                       Cloud_Model         = trim(self%conf%Cloud_Model)        , &
                       CloudCoeff_Format   = trim(self%conf%CloudCoeff_Format)  , &
                       CloudCoeff_File     = trim(self%conf%CloudCoeff_File)    , &
                       IRwaterCoeff_File   = trim(self%conf%IRwaterCoeff_File)  , &
                       IRlandCoeff_File    = trim(self%conf%IRlandCoeff_File)   , &
                       IRsnowCoeff_File    = trim(self%conf%IRsnowCoeff_File)   , &
                       IRiceCoeff_File     = trim(self%conf%IRiceCoeff_File)    , &
                       VISwaterCoeff_File  = trim(self%conf%VISwaterCoeff_File) , &
                       VISlandCoeff_File   = trim(self%conf%VISlandCoeff_File)  , &
                       VISsnowCoeff_File   = trim(self%conf%VISsnowCoeff_File)  , &
                       VISiceCoeff_File    = trim(self%conf%VISiceCoeff_File)   , &
                       MWwaterCoeff_File   = trim(self%conf%MWwaterCoeff_File)  , &
                       Quiet               = .TRUE.)

 message = 'Error initializing CRTM'
 call crtm_comm_stat_check(err_stat, PROGRAM_NAME, message, f_comm)
 
 ! Loop over all sensors. Not necessary if we're calling CRTM for each sensor
 ! ----------------------------------------------------------------------------
 Sensor_Loop:do n = 1, self%conf%n_Sensors


   ! Pass channel list to CRTM
   ! -------------------------
   err_stat = CRTM_ChannelInfo_Subset(chinfo(n), self%channels, reset=.false.)
   message = 'Error subsetting channels!'
   call crtm_comm_stat_check(err_stat, PROGRAM_NAME, message, f_comm)

   ! Determine the number of channels for the current sensor
   ! -------------------------------------------------------
   n_Channels = CRTM_ChannelInfo_n_Channels(chinfo(n))

   ! Allocate the ARRAYS (for CRTM_Forward)
   ! --------------------------------------
   allocate( geo( n_Profiles ),               &
             atm( n_Profiles ),               &
             sfc( n_Profiles ),               &
             rts( n_Channels, n_Profiles ),   &
             Options( n_Profiles ),           &
             STAT = alloc_stat )
   message = 'Error allocating structure arrays'
   call crtm_comm_stat_check(alloc_stat, PROGRAM_NAME, message, f_comm)

   if (n_Layers > 0) call CRTM_RTSolution_Create (rts, n_Layers)

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
   call CRTM_Surface_Create(sfc, n_Channels)
   IF ( ANY(.NOT. CRTM_Surface_Associated(sfc)) ) THEN
      message = 'Error allocating CRTM Surface structure'
      CALL Display_Message( PROGRAM_NAME, message, FAILURE )
      STOP
   END IF

   CALL CRTM_RTSolution_Create(rts, n_Layers )

   ! Some special treatment for Vis or UV as compared to IR or MW.
   ! ----------------------------------------
   if (SpcCoeff_IsVisibleSensor(SC(n)) .or. SpcCoeff_IsUltravioletSensor(SC(n))) then
      Is_Vis_or_UV = .true.
   else
      Is_Vis_or_UV = .false.
   endif

   !Assign the data from the GeoVaLs
   !--------------------------------
   call Load_Atm_Data(n_Profiles,n_Layers,geovals,atm,self%conf, SC(n)%Is_Active_Sensor)

   call Load_Sfc_Data(n_Profiles,n_Channels,self%channels,geovals,sfc,chinfo,obss,self%conf, &
                      SC(n)%Is_Active_Sensor, Is_Vis_or_UV)
   if (cmp_strings(self%conf%SENSOR_ID(n),'gmi_gpm')) then
     allocate( geo_hf( n_Profiles ))
     call Load_Geom_Data(obss,geo,geo_hf,self%conf%SENSOR_ID(n))
   else
     call Load_Geom_Data(obss,geo)
   endif

   ! Call THE CRTM inspection
   ! ------------------------
   if (self%conf%inspect > 0) then
     call CRTM_Atmosphere_Inspect(atm(self%conf%inspect))
     call CRTM_Surface_Inspect(sfc(self%conf%inspect))
     call CRTM_Geometry_Inspect(geo(self%conf%inspect))
     call CRTM_ChannelInfo_Inspect(chinfo(n))
   endif

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
 999  str_pos(1) = index(varstr,jacobianstr) - 1        !position before jacobianstr
      if (str_pos(1) == 0) then
         write(err_msg,*) 'ufo_radiancecrtm_simobs: _jacobian_ must be // &
                           & preceded by dependent variable in config: ', &
                           & hofxdiags%variables(jvar)
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

   ! set profiles that should be skipped
   call ufo_crtm_skip_profiles(n_Profiles,n_Channels,self%channels,obss,atm,sfc,  &
                               SC(n)%Is_Active_Sensor,Is_Vis_or_UV,Options)
   if ( self%use_qc_flags ) then
   ! eliminate remaining profiles that are "QCed" out
     n_skipped = 0
     do m = 1, n_Profiles
       if (.not.Options(m)%Skip_Profile) then
         ml = m
         skip_prof = .true.
! can't use all here since the qc flags are out of order and not an array
         do l = 1, size(self%channels)
           ll = l
           qc_ff = qc_flags%get(ll,ml)
           if ( qc_ff < 2) then 
             skip_prof = .false.
             exit
           end if
         end do
         if ( skip_prof ) then
           n_skipped = n_skipped + 1
         end if
         Options(m)%Skip_Profile = skip_prof
       end if
     end do
   end if
   
   if (jacobian_needed) then
      ! Allocate the ARRAYS (for CRTM_K_Matrix)
      ! --------------------------------------
      allocate( atm_K( n_Channels, n_Profiles ),               &
                sfc_K( n_Channels, n_Profiles ),   &
                rts_K( n_Channels, n_Profiles ),   &
                STAT = alloc_stat )
      message = 'Error allocating K structure arrays'
      call crtm_comm_stat_check(alloc_stat, PROGRAM_NAME, message, f_comm)

      ! Create output K-MATRIX structure (atm)
      ! --------------------------------------
      call CRTM_Atmosphere_Create( atm_K, n_Layers, self%conf%n_Absorbers, self%conf%n_Clouds, self%conf%n_Aerosols )
      if ( ANY(.NOT. CRTM_Atmosphere_Associated(atm_K)) ) THEN
         message = 'Error allocating CRTM K-matrix Atmosphere structure (setTraj)'
         CALL Display_Message( PROGRAM_NAME, message, FAILURE )
         STOP
      END IF

      ! Create output K-MATRIX structure (sfc)
      ! --------------------------------------
      call CRTM_Surface_Create( sfc_K, n_Channels)
      IF ( ANY(.NOT. CRTM_Surface_Associated(sfc_K)) ) THEN
         message = 'Error allocating CRTM K-matrix Surface structure (setTraj)'
         CALL Display_Message( PROGRAM_NAME, message, FAILURE )
         STOP
      END IF

      ! Zero the K-matrix OUTPUT structures
      ! -----------------------------------
      call CRTM_Atmosphere_Zero( atm_K )
      call CRTM_Surface_Zero( sfc_K )

      ! Inintialize the K-matrix INPUT so that the results are dTb/dx or dR/dx
      ! -------------------------------------------------------------
      if (SC(n)%Is_Active_Sensor) then
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
      else if (Is_Vis_or_UV) then
         rts_K%Radiance                = ONE
         rts_K%Brightness_Temperature  = ZERO
      else
         rts_K%Radiance                = ZERO
         rts_K%Brightness_Temperature  = ONE
      end if

      ! Call the K-matrix model
      ! -----------------------
      err_stat = CRTM_K_Matrix( atm         , &  ! FORWARD  Input
                                sfc         , &  ! FORWARD  Input
                                rts_K       , &  ! K-MATRIX Input
                                geo         , &  ! Input
                                chinfo(n:n) , &  ! Input
                                atm_K       , &  ! K-MATRIX Output
                                sfc_K       , &  ! K-MATRIX Output
                                rts         , &  ! FORWARD  Output
                                Options       )  ! Input
      message = 'Error calling CRTM (setTraj) K-Matrix Model for '//TRIM(self%conf%SENSOR_ID(n))
      call crtm_comm_stat_check(err_stat, PROGRAM_NAME, message, f_comm)
      if (cmp_strings(self%conf%SENSOR_ID(n),'gmi_gpm')) then
         allocate( atm_Ka( n_Channels, n_Profiles ),               &
                   sfc_Ka( n_Channels, n_Profiles ),   &
                   rts_Ka( n_Channels, n_Profiles ),   &
                   rtsa( n_Channels, n_Profiles ),     &
                   STAT = alloc_stat )
         message = 'Error allocating K structure arrays rtsa, atm_Ka ......'
         call crtm_comm_stat_check(alloc_stat, PROGRAM_NAME, message, f_comm)
         !! save resutls for gmi channels 1-9.
         atm_Ka = atm_K
         sfc_Ka = sfc_K
         rts_Ka = rts_K
         rtsa   = rts
         !! call CRTM_K_Matrix again for geo_hf which has view angle for gmi channels 10-13.
         call CRTM_Atmosphere_Zero( atm_K )
         call CRTM_Surface_Zero( sfc_K )
         rts_K%Radiance               = ZERO
         rts_K%Brightness_Temperature = ONE
         ! Call the K-matrix model
         ! -----------------------
         err_stat = CRTM_K_Matrix( atm         , &  ! FORWARD  Input
                                   sfc         , &  ! FORWARD  Input
                                   rts_K       , &  ! K-MATRIX Input
                                   geo_hf        , &  ! Input
                                   chinfo(n:n) , &  ! Input
                                   atm_K       , &  ! K-MATRIX Output
                                   sfc_K       , &  ! K-MATRIX Output
                                   rts         , &  ! FORWARD  Output
                                   Options       )  ! Input
         message = 'Error calling CRTM (setTraj, geo_hf) K-Matrix Model for ' &
                   //TRIM(self%conf%SENSOR_ID(n))
         call crtm_comm_stat_check(err_stat, PROGRAM_NAME, message, f_comm)
         !! replace data for gmi channels 1-9 by early results calculated with geo.
         do l = 1, size(self%channels)
            if ( self%channels(l) <= 9 ) then
               atm_K(l,:) = atm_Ka(l,:)
               sfc_K(l,:) = sfc_Ka(l,:)
               rts_K(l,:) = rts_Ka(l,:)
               rts(l,:)   = rtsa(l,:)
            endif
         enddo
         deallocate(atm_Ka,sfc_Ka,rts_Ka,rtsa)
      endif ! cmp_strings(self%conf%SENSOR_ID(n),'gmi_gpm')
   else
      ! Call the forward model call for each sensor
      ! -------------------------------------------
      err_stat = CRTM_Forward( atm         , &  ! Input
                               sfc         , &  ! Input
                               geo         , &  ! Input
                               chinfo(n:n) , &  ! Input
                               rts         , &  ! Output
                               Options       )  ! Input
      message = 'Error calling CRTM Forward Model for '//TRIM(self%conf%SENSOR_ID(n))
      call crtm_comm_stat_check(err_stat, PROGRAM_NAME, message, f_comm)
      if (cmp_strings(self%conf%SENSOR_ID(n),'gmi_gpm')) then
         allocate( rtsa( n_Channels, n_Profiles ),     &
                   STAT = alloc_stat )
         message = 'Error allocating K structure arrays rtsa.'
         call crtm_comm_stat_check(alloc_stat, PROGRAM_NAME, message, f_comm)
         !! save resutls for gmi channels 1-9.
         rtsa = rts
         !! call crtm again for gmi channels 10-13 with geo_hf.
         ! -----------------------
         err_stat = CRTM_Forward( atm         , &  ! Input
                                  sfc         , &  ! Input
                                  geo_hf        , &  ! Input
                                  chinfo(n:n) , &  ! Input
                                  rts         , &  ! Output
                                  Options       )  ! Input
         message = 'Error calling CRTM Forward Model for gmi_gpm channels 10-13'
         call crtm_comm_stat_check(err_stat, PROGRAM_NAME, message, f_comm)
         !! replace data for gmi channels 1-9 by results calculated with geo.
         do l = 1, size(self%channels)
            if ( self%channels(l) <= 9 ) then
               rts(l,:)   = rtsa(l,:)
            endif
         enddo
         deallocate(rtsa)
      endif ! cmp_strings(self%conf%SENSOR_ID(n),'gmi_gpm')
   end if ! jacobian_needed

   !call CRTM_RTSolution_Inspect(rts)

   ! put CRTM output into hofxdiags
   if (SC(n)%Is_Active_Sensor) then
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
      call ufo_crtm_passive_sim(rts, &
                                Options, &
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
       write(err_msg,*) 'ufo_radiancecrtm_sim error: failed to put simulated diagnostics into hofxdiags'
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
   deallocate(geo, atm, sfc, rts, Options, STAT = alloc_stat)
   if(allocated(geo_hf)) deallocate(geo_hf)
   message = 'Error deallocating structure arrays'
   call crtm_comm_stat_check(alloc_stat, PROGRAM_NAME, message, f_comm)

   if (jacobian_needed) then
      ! Deallocate the K structures
      ! ---------------------------
      call CRTM_Atmosphere_Destroy(atm_K)
      call CRTM_Surface_Destroy(sfc_K)
      call CRTM_RTSolution_Destroy(rts_K)

      ! Deallocate all K arrays
      ! -----------------------
      deallocate(atm_K, sfc_K, rts_K, STAT = alloc_stat)
      message = 'Error deallocating K structure arrays'
      call crtm_comm_stat_check(alloc_stat, PROGRAM_NAME, message, f_comm)
   end if

 end do Sensor_Loop


 ! Destroy CRTM instance
 ! ---------------------
 ! write( *, '( /5x, "Destroying the CRTM..." )' )
 err_stat = CRTM_Destroy( chinfo )
 message = 'Error destroying CRTM'
 call crtm_comm_stat_check(err_stat, PROGRAM_NAME, message, f_comm)

end subroutine ufo_radiancecrtm_simobs

! ------------------------------------------------------------------------------

end module ufo_radiancecrtm_mod
