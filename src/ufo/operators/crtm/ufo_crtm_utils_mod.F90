! (C) Copyright 2018-2026 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

!> Fortran module to provide code shared between nonlinear and tlm/adm radiance calculations

module ufo_crtm_utils_mod

use fckit_configuration_module, only: fckit_configuration
use fckit_exception_module, only: fckit_exception
use fckit_mpi_module,   only: fckit_mpi_comm
use, intrinsic :: iso_c_binding
use kinds

use CRTM_Module

use ufo_geovals_mod, only: ufo_geovals, ufo_geoval, ufo_geovals_get_var
use ufo_vars_mod
use ufo_utils_mod, only: upper2lower
use ufo_ghg_utils_mod, only: co2_mauna_loa_rise
use obsspace_mod
use ufo_constants_mod, only: kg_to_g, co2_rescale_to_ppmv, co2_ppmv_value, midpoint_julday, &
                             t0c, grav, qsmall, &
                             cloudfractionsmall, aerosol_concentration_minvalue_layer, &
                             ozone_default_value, threshold_reflectivity

implicit none
private

public ufo_crtm_comm_stat_check
public crtm_conf
public nullify_jac_channels_type
public ufo_crtm_conf_setup
public ufo_crtm_conf_delete
public ufo_crtm_init_aerosol
public ufo_crtm_init_radiance
public ufo_crtm_register_radiance
public ufo_crtm_register_aerosol
public ufo_crtm_deregister_radiance
public ufo_crtm_deregister_aerosol
public ufo_crtm_sensor_type
public ufo_crtm_sensor_wavelengths
public ufo_crtm_check_absorbers
public ufo_crtm_landtyp_var
public ufo_crtm_destroy
public ufo_crtm_gmi_split_chinfo
public ufo_crtm_get_var_name
public ufo_crtm_load_atm_data
public ufo_crtm_load_sfc_data
public ufo_crtm_load_geom_data
public ufo_crtm_pressure_is_monotonic
public ufo_crtm_skip_profiles
public ufo_crtm_skip_qc_profiles

public ufo_crtm_load_aerosol_data
public ufo_crtm_assign_aerosol_names
public ufo_crtm_define_aerosol_model
public ufo_crtm_calculate_aero_layer_factor

integer, parameter, public :: max_string=800

! Set once the CRTM library version has been logged (once per process)
logical :: crtm_version_logged = .false.
! Set once missing gmi_gpm channel 10-13 geometry has been warned about (once per process)
logical :: gmi_hf_geometry_warned = .false.

! Shared CRTM coefficient loads.
! Instead of a CRTM_Init/CRTM_Destroy on every simobs/settraj call, each operator
! registers its sensors at setup (ufo_crtm_register_radiance/_aerosol). Operators
! with the same coefficient settings form a group, loaded with one CRTM_Init for
! the union of their sensors. CRTM holds one load at a time, so only one group is
! in memory; using another group replaces it. Each operator deregisters its sensors
! when deleted (ufo_crtm_deregister_radiance/_aerosol); a sensor leaves its group
! when its last operator goes, and the load is destroyed when its group is empty.
! Operators with ShareCoefficients: false do not register and load their own
! coefficients for each call, as before shared loads.
! The variables below hold this state for the whole process.
! CRTM_Init and CRTM_Destroy must not be called anywhere else in the executable:
! another load while this one is in memory fails.
type shared_load_group
  character(len=:),   allocatable :: key            ! coefficient settings
  logical                         :: aerosol        ! aerosol-only (vs radiance) CRTM_Init
  character(len=255), allocatable :: sensor_ids(:)  ! registered sensors
  integer,            allocatable :: n_users(:)     ! operators using each sensor
end type shared_load_group
type(shared_load_group),     allocatable :: shared_groups(:)
type(CRTM_ChannelInfo_type), allocatable :: shared_chinfo(:)  ! from CRTM_Init
! group whose load is in CRTM's memory (0 = none)
integer :: loaded_group = 0
! .true. while the CRTM coefficients in memory are a shared load
logical :: shared_loaded = .false.
! Set once more than one group has been warned about (once per process)
logical :: shared_groups_warned = .false.

!Type to hold the list of channels to zero for one nullify jacobian params entry
type nullify_jac_channels_type
 integer, allocatable :: channels(:)
end type nullify_jac_channels_type

!Type for general config
type crtm_conf
 integer :: n_Sensors
 integer :: n_Absorbers
 integer :: n_Clouds
 integer :: n_Aerosols
 integer :: n_Surfaces
 character(len=MAXVARLEN), allocatable :: Absorbers(:)
 integer, allocatable :: Absorber_Id(:)
 integer, allocatable :: Absorber_Units(:)
 character(len=MAXVARLEN), allocatable :: Clouds(:,:)
 integer, allocatable :: Cloud_Id(:)
 character(len=MAXVARLEN), allocatable :: Surfaces(:)
 character(len=MAXVARLEN), allocatable :: Aerosols(:)  ! aerosol species of aerosol_option

 character(len=255), allocatable :: SENSOR_ID(:)
 character(len=255) :: ENDIAN_TYPE
 character(len=255) :: COEFFICIENT_PATH, NC_COEFFICIENT_PATH
 character(len=255) :: TauCoeff_Format, SpcCoeff_Format
 character(len=255) :: CloudCoeff_Format, AerosolCoeff_Format
 character(len=255) :: Aerosol_Model, Cloud_Model
 character(len=255) :: &
    IRwaterCoeff_File, IRlandCoeff_File, IRsnowCoeff_File, IRiceCoeff_File, &
    VISwaterCoeff_File, VISlandCoeff_File, VISsnowCoeff_File, VISiceCoeff_File, &
    MWwaterCoeff_File, CloudCoeff_File, AerosolCoeff_File
 integer, allocatable :: Land_WSI(:)
 real(kind_real) :: Cloud_Fraction = -1.0_kind_real
 integer :: inspect
 character(len=MAXVARLEN) :: aerosol_option
 character(len=255) :: salinity_option
 character(len=MAXVARLEN) :: sfc_wind_geovars
 character(len=MAXVARLEN) :: co2_method_geovars
 real(kind_real) :: unit_coef = 1.0_kind_real
 real(kind_real) :: co2_rescale_to_ppmv = co2_rescale_to_ppmv
 real(kind_real) :: co2_ppmv_value = co2_ppmv_value
 real(kind_real) :: midpoint_julday = midpoint_julday
 ! nullifyJacVars(i) is zeroed for channels in nullifyJacChannels(i)%channels
 character(len=MAXVARLEN), allocatable :: nullifyJacVars(:)
 type(nullify_jac_channels_type), allocatable :: nullifyJacChannels(:)
 logical :: Cloud_Seeding = .false.
 logical :: cal_cloud_frac_in_fov = .false.
 logical :: cal_cloud_reff_in_fov  = .false.
 logical :: precip_hydro = .false.
 logical :: flag_deep_conv_mass_flux = .true.
 logical :: read_Cmatrix = .false.
 character(len=max_string) :: Cmatrix_path
 logical :: share_coefficients = .true.  ! register for the shared CRTM load
end type crtm_conf


interface ufo_crtm_calculate_aero_layer_factor

   module procedure calculate_aero_layer_factor_atm_profile, &
                    calculate_aero_layer_factor_atm

end interface


 ! Add more UFO_Absorbers as needed
 ! Note: must have same ordering as CRTM_Absorbers,
 !       CRTM_Absorber_Id, and CRTM_Absorber_Units
 character(len=MAXVARLEN), parameter :: &
      UFO_Absorbers(3) = &
         [ var_mixr, var_co2, var_oz ]

 ! copy of ABSORBER_ID_NAME defined in CRTM_Atmosphere_Define
 character(len=MAXVARLEN), parameter :: &
      CRTM_Absorbers(N_VALID_ABSORBER_IDS) = &
         ABSORBER_ID_NAME(1:N_VALID_ABSORBER_IDS)
 integer, parameter :: &
      CRTM_Absorber_Id(N_VALID_ABSORBER_IDS) = &
         [ H2O_ID,  CO2_ID,   O3_ID,  N2O_ID, &
            CO_ID,  CH4_ID,   O2_ID,   NO_ID, &
           SO2_ID,  NO2_ID,  NH3_ID, HNO3_ID, &
            OH_ID,   HF_ID,  HCl_ID,  HBr_ID, &
            HI_ID,  ClO_ID,  OCS_ID, H2CO_ID, &
          HOCl_ID,   N2_ID,  HCN_ID, CH3l_ID, &
          H2O2_ID, C2H2_ID, C2H6_ID,  PH3_ID, &
          COF2_ID,  SF6_ID,  H2S_ID,HCOOH_ID ]
 integer, parameter :: &
      CRTM_Absorber_Units(3) = [ &
         MASS_MIXING_RATIO_UNITS   & !H2O
       , VOLUME_MIXING_RATIO_UNITS & !CO2
       , VOLUME_MIXING_RATIO_UNITS & !O3
        ]

 ! I. Moradi: MAX_VALID_CLOUD_CATEGORIES is defined to declare the following variables that can
 ! work CRTM V3.0 and greater of the CRTM interface
 ! generally requires some considerations in terms of how clouds are defined
 integer, parameter :: MAX_VALID_CLOUD_CATEGORIES = 25  ! 6 default + 19 from DDA table

 character(len=MAXVARLEN), parameter :: &
      UFO_CLOUDS(MAX_VALID_CLOUD_CATEGORIES, 2) = &
       reshape( &
         [ var_clw_wp,  var_cli_wp,  var_clr_wp,  var_cls_wp,  var_clg_wp, &    ! 1- 5
           var_clh_wp,  var_cls_wp,  var_cls_wp,  var_cls_wp,  var_cls_wp, &    ! 5-10
           var_cls_wp,  var_cli_wp,  var_cls_wp,  var_cls_wp,  var_cls_wp, &    ! 11-15
           var_cls_wp,  var_cls_wp,  var_cls_wp,  var_cls_wp,  var_clh_wp, &    ! 16-20
           var_clg_wp,  var_cls_wp,  var_clh_wp,  var_cli_wp,  var_clw_wp, &    ! 21-25
           var_clwefr, var_cliefr, var_clrefr, var_clsefr, var_clgefr, &   ! 1- 5
           var_clhefr, var_clsefr, var_clsefr, var_clsefr, var_clsefr, &   ! 5-10
           var_clsefr, var_cliefr, var_clsefr, var_clsefr, var_clsefr, &   ! 11-15
           var_clsefr, var_clsefr, var_clsefr, var_clsefr, var_clhefr, &   ! 16-20
           var_clgefr, var_clsefr, var_clhefr, var_cliefr, var_clwefr] &   ! 21-25
           , [MAX_VALID_CLOUD_CATEGORIES,2] )

 ! copy of CLOUD_CATEGORY_NAME defined in CRTM_Cloud_Define
 character(len=MAXVARLEN), parameter :: &
      CRTM_Clouds(N_VALID_CLOUD_CATEGORIES) = &
         CLOUD_CATEGORY_NAME(1:N_VALID_CLOUD_CATEGORIES)

 ! The hydrometeor types from the DDA tables are defined using corresponding integer
 ! values so  the ufo interface does not break for the current version of CRTM
 ! implemented in the UFO
 integer, parameter :: &
      CRTM_Cloud_Id(MAX_VALID_CLOUD_CATEGORIES) = &
         [ WATER_CLOUD                   , &
           ICE_CLOUD                     , &
           RAIN_CLOUD                    , &
           SNOW_CLOUD                    , &
           GRAUPEL_CLOUD                 , &
           HAIL_CLOUD                    , &
           PlateType1                    , &
           ColumnType1                   , &
           SixBulletRosette              , &
           Perpendicular4_BulletRosette  , &
           Flat3_BulletRosette           , &
           IconCloudIce                  , &
           SectorSnowflake               , &
           EvansSnowAggregate            , &
           EightColumnAggregate          , &
           LargePlateAggregate           , &
           LargeColumnAggregate          , &
           LargeBlockAggregate           , &
           IconSnow                      , &
           IconHail                      , &
           GemGraupel                    , &
           GemSnow                       , &
           GemHail                       , &
           IceSphere                     , &
           LiquidSphere                  ]

! Clouds used in the Thompson microphysics scheme to calculate cloud fraction or particle radius
 integer, parameter :: &
      thompson_Cloud_Id(5) = &
         [ WATER_CLOUD                   , &
           ICE_CLOUD                     , &
           RAIN_CLOUD                    , &
           SNOW_CLOUD                    , &
           GRAUPEL_CLOUD                 ]
 character(len=MAXVARLEN), parameter :: &
      thompson_Cloud_var(5) = &
         [ var_qc, var_qi, var_qr, var_qs, var_qg ]
! A parameter for mass flux deep convection in the Thompson microphysics scheme
! to calculate cloud fraction
 real(kind_real) :: xrc3 = 200.0_kind_real

! Surface Variables

 character(len=MAXVARLEN), parameter :: &
      UFO_Surfaces(7) = &
         [ var_sfc_wtmp, var_sfc_ltmp, var_sfc_itmp, var_sfc_stmp, var_sfc_wspeed, var_sfc_wdir, &
           var_sfc_sss ]

 character(len=MAXVARLEN), parameter :: &
      CRTM_Surfaces(7) = &
         [  character(len=MAXVARLEN):: "Water_Temperature", "Land_Temperature", "Ice_Temperature", &
                                       "Snow_Temperature", "Wind_Speed", &
                                       "Wind_Direction", "Salinity" ]

 character(len=MAXVARLEN), parameter :: &
      ValidSurfaceWindGeoVars(2) = [character(len=MAXVARLEN) :: "vector", "uv"]

 character(len=MAXVARLEN), parameter :: &
      ValidCO2AbsorberMethod(3) = [character(len=MAXVARLEN) :: "Background", &
                                                               "EternalConstant", &
                                                               "MaunaLoa"]



contains

! ------------------------------------------------------------------------------

logical function ufo_crtm_pressure_is_monotonic(atm)

implicit none

type(CRTM_Atmosphere_type), intent(in) :: atm

integer :: jlevel

 ufo_crtm_pressure_is_monotonic = .true.
 do jlevel = atm%n_layers, 1, -1
   if (atm%level_pressure(jlevel) <= atm%level_pressure(jlevel-1)) then
     ufo_crtm_pressure_is_monotonic = .false.
     return
   end if
 end do

end function ufo_crtm_pressure_is_monotonic

! ------------------------------------------------------------------------------

subroutine ufo_crtm_conf_setup(conf, f_confOpts, f_confOper, midPointJulday, comm)

implicit none

type(crtm_conf),                intent(inout) :: conf
type(fckit_configuration),      intent(in)    :: f_confOpts
type(fckit_configuration),      intent(in)    :: f_confOper
integer(c_int64_t),             intent(in)    :: midPointJulday
type(fckit_mpi_comm), optional, intent(in)    :: comm

character(len=*), parameter :: routine_name = "ufo_crtm_conf_setup"
character(len=255) :: IRwaterCoeff, VISwaterCoeff, &
                      IRVISlandCoeff, IRVISsnowCoeff, IRVISiceCoeff, &
                      MWwaterCoeff
integer :: jspec, ivar
character(len=max_string) :: message
character(len=:), allocatable :: str
character(len=:), allocatable :: str_array(:)
integer, allocatable :: int_array(:)
type(fckit_configuration), allocatable :: nullifyJacParams_confs(:)

logical, save :: message_flag = .true.

 !Some config needs to come from user
 !-----------------------------------

 !Number of sensors, each call to CRTM will be for a single sensor
 !type (zenith/scan angle will be different)
 conf%n_Sensors = 1

 ! absorbers, clouds and aerosols (should match what model will provide)

 ! Set print rank
 ! --------------
 if (present(comm)) then
   if (comm%rank() /= 0) message_flag = .false.
 end if

 ! Absorbers
 !----------
 conf%n_Absorbers = 0
 if (f_confOper%has("Absorbers")) conf%n_Absorbers = f_confOper%get_size("Absorbers")

 allocate( conf%Absorbers     ( conf%n_Absorbers ), &
           conf%Absorber_Id   ( conf%n_Absorbers ), &
           conf%Absorber_Units( conf%n_Absorbers ) )

 if (conf%n_Absorbers > 0) then
   call f_confOper%get_or_die("Absorbers",str_array)
   conf%Absorbers(1:conf%n_Absorbers) = str_array
 end if

 call conf_check_duplicates(conf%Absorbers, "Absorbers", routine_name)

 ! convert from CRTM names to UFO CF names and define Id and Units
 do jspec = 1, conf%n_Absorbers
   ivar = ufo_vars_getindex(CRTM_Absorbers, conf%Absorbers(jspec))
   if (ivar < 1 .or. ivar > size(UFO_Absorbers)) then
     message = trim(routine_name) // " error: " // trim(conf%Absorbers(jspec)) // &
               " not supported by UFO_Absorbers"
     call fckit_exception % throw(message)
   end if
   conf%Absorbers(jspec) = UFO_Absorbers(ivar)
   conf%Absorber_Id(jspec) = CRTM_Absorber_Id(ivar)
   ! NOTE JEDI's units for the "water_vapor_mixing_ratio_wrt_dry_air" is kg/kg,
   ! but here we're passing CRTM's g/kg
   conf%Absorber_Units(jspec) = CRTM_Absorber_Units(ivar)
 end do


 ! Clouds
 !-------
 conf%n_Clouds = 0
 if (f_confOper%has("Clouds")) conf%n_Clouds = f_confOper%get_size("Clouds")
 allocate( conf%Clouds  ( conf%n_Clouds,2), &
           conf%Cloud_Id( conf%n_Clouds ) )
 if (conf%n_Clouds > 0) then
   conf%cal_cloud_frac_in_fov = conf_thompson_option(f_confOper, &
                                  "method for cloud fraction within fov", routine_name)
   !  Scale-aware mass-flux deep conv scheme flag used by the Thompson cloud fraction
   if (conf%cal_cloud_frac_in_fov .and. f_confOper%has("convection_mass_flux_flag")) then
     call f_confOper%get_or_die("convection_mass_flux_flag",conf%flag_deep_conv_mass_flux)
   end if
   conf%cal_cloud_reff_in_fov = conf_thompson_option(f_confOper, &
                                  "method for hydrometeor effective radii within fov", routine_name)

   call f_confOper%get_or_die("Clouds",str_array)
   conf%Clouds(1:conf%n_Clouds,1) = str_array

   call conf_cloud_fraction(conf, f_confOper, routine_name, message_flag)

   if (f_confOper%has("Cloud_Seeding")) then
     call f_confOper%get_or_die("Cloud_Seeding",conf%Cloud_Seeding)
     if ( conf%Cloud_Seeding ) then
       message = trim(routine_name) // ": Cloud Seeding is activated"
     else
       message = trim(routine_name) // ": Cloud Seeding is not activated"
     end if
     if (message_flag) call Display_Message(routine_name, trim(message), INFORMATION )
   end if
 end if

 call conf_check_duplicates(conf%Clouds(:,1), "Clouds", routine_name)

 ! convert from CRTM names to UFO CF names and define Id
 do jspec = 1, conf%n_Clouds
   ivar = ufo_vars_getindex(CRTM_Clouds, conf%Clouds(jspec,1))
   if (ivar < 1 .or. ivar > size(UFO_Clouds)) then
     message = trim(routine_name) // " error: " // trim(conf%Clouds(jspec,1)) // &
               " not supported by UFO_Clouds"
     call fckit_exception % throw(message)
   end if

   conf%Clouds(jspec,1:2) = UFO_Clouds(ivar,1:2)
   conf%Cloud_Id(jspec)   = CRTM_Cloud_Id(ivar)
   if (conf%Cloud_Id(jspec) == RAIN_CLOUD .or. conf%Cloud_Id(jspec) == SNOW_CLOUD .or. &
        conf%Cloud_Id(jspec) == GRAUPEL_CLOUD .or. conf%Cloud_Id(jspec) == HAIL_CLOUD) then
     conf%precip_hydro = .true.
   end if
 end do

 ! Aerosols
 !---------
 if (f_confOpts%get("AerosolOption",str)) then
    conf%aerosol_option = upper2lower(str)
    call ufo_crtm_assign_aerosol_names(conf%aerosol_option, conf%Aerosols)
 else
    conf%aerosol_option = ""
    allocate(conf%Aerosols(0))
 end if
 conf%n_Aerosols = size(conf%Aerosols)
 call f_confOpts%get_or_die("model units coeff", conf%unit_coef)

 ! Surface variables
 !----------
 conf%n_Surfaces = 0
 if (f_confOper%has("Surfaces")) conf%n_Surfaces = f_confOper%get_size("Surfaces")

 allocate( conf%Surfaces    ( conf%n_Surfaces ))

 if (conf%n_Surfaces > 0) then
   call f_confOper%get_or_die("Surfaces",str_array)
   conf%Surfaces(1:conf%n_Surfaces) = str_array
 end if

 call conf_check_duplicates(conf%Surfaces, "Surfaces", routine_name)

 ! convert from CRTM names to UFO CF names
 do jspec = 1, conf%n_Surfaces
   ivar = ufo_vars_getindex(CRTM_Surfaces, conf%Surfaces(jspec))
   if (ivar < 1 .or. ivar > size(UFO_Surfaces)) then
     message = trim(routine_name) // " error: " // trim(conf%Surfaces(jspec)) // &
               " not supported by UFO_Surfaces"
     call fckit_exception % throw(message)
   end if
   conf%Surfaces(jspec) = UFO_Surfaces(ivar)
 end do

 ! select between two surface wind geovals options
 ! valid options: vector [default], uv
 conf%sfc_wind_geovars = "vector"
 if (f_confOper%get("SurfaceWindGeoVars", str)) conf%sfc_wind_geovars = str
 if (ufo_vars_getindex(ValidSurfaceWindGeoVars, conf%sfc_wind_geovars) < 1) then
   message = trim(routine_name) // " error: invalid SurfaceWindGeoVars " // &
             conf%sfc_wind_geovars
   call fckit_exception % throw(message)
 end if

 ! Select among one of different methods to handling CO2 Absorber
 ! --------------------------------------------------------------
 conf%co2_method_geovars = "Background"
 if (f_confOper%get("CO2AbsorberMethod", str)) conf%co2_method_geovars = str
 if (ufo_vars_getindex(ValidCO2AbsorberMethod, conf%co2_method_geovars) < 1) then
   message = trim(routine_name) // " error: invalid CO2 method " // conf%co2_method_geovars
   call fckit_exception % throw(message)
 end if

 ! Per-variable list of channels whose Jacobian is zeroed
 ! ---------------------------------------------------------------------------------------
 if (f_confOper%has("nullify jacobian params")) then
   call f_confOper%get_or_die("nullify jacobian params", nullifyJacParams_confs)

   allocate(conf%nullifyJacVars(size(nullifyJacParams_confs)))
   allocate(conf%nullifyJacChannels(size(nullifyJacParams_confs)))

   do jspec = 1, size(nullifyJacParams_confs)
     call nullifyJacParams_confs(jspec)%get_or_die("variable", str)
     conf%nullifyJacVars(jspec) = str

     call nullifyJacParams_confs(jspec)%get_or_die("channels", int_array)
     conf%nullifyJacChannels(jspec)%channels = int_array
   end do
 else
   allocate(conf%nullifyJacVars(0))
   allocate(conf%nullifyJacChannels(0))
 end if

 ! observation midpoint Julian day in proleptic Julian calendar (from 01Jan 4713BC)
 ! --------------------------------------------------------------------------------
 conf%midpoint_julday = midPointJulday

 ! import a GeoVaLs scaling factor
 ! -------------------------------
 if (f_confOper%has("CO2gvConvertUnit")) then
   call f_confOper%get_or_die("CO2gvConvertUnit", conf%co2_rescale_to_ppmv)
 end if

 ! import a ppmv value of CO2, only used if EternalConstant is selected
 ! --------------------------------------------------------------------
 if (f_confOper%has("CO2ppmvValue")) then
   call f_confOper%get_or_die("CO2ppmvValue", conf%co2_ppmv_value)
 end if

 ! Sea_Surface_Salinity
 !---------
 conf%salinity_option = "off"
 if (f_confOpts%get("Salinity",str)) conf%salinity_option = str

 !Allocate SENSOR_ID
 allocate(conf%SENSOR_ID(conf%n_Sensors))

 !Get sensor ID from config
 call f_confOpts%get_or_die("Sensor_ID",str)
 conf%SENSOR_ID(conf%n_Sensors) = str

 !ENDIAN type
 call f_confOpts%get_or_die("EndianType",str)
 conf%ENDIAN_TYPE = str

 !Path to coefficient files
 call f_confOpts%get_or_die("CoefficientPath",str)
 conf%COEFFICIENT_PATH = str

 !Path to NetCDF coefficient files
 conf%NC_COEFFICIENT_PATH = conf%COEFFICIENT_PATH
 if (f_confOpts%get("NC_CoefficientPath",str)) conf%NC_COEFFICIENT_PATH = str

 ! Spc and Tau coefficient file format
 conf%SpcCoeff_Format = "Binary"
 if (f_confOpts%get("SpcCoeff_Format",str)) conf%SpcCoeff_Format = str
 conf%TauCoeff_Format = "Binary"
 if (f_confOpts%get("TauCoeff_Format",str)) conf%TauCoeff_Format = str

 ! Cloud coefficient model, file, and format
 conf%Cloud_Model = "CRTM"
 if (f_confOpts%get("Cloud_Model",str)) conf%Cloud_Model = str
 conf%CloudCoeff_File = "CloudCoeff.bin"
 if (f_confOpts%get("CloudCoeff_File",str)) conf%CloudCoeff_File = str
 conf%CloudCoeff_Format = "Binary"
 if (f_confOpts%get("CloudCoeff_Format",str)) conf%CloudCoeff_Format = str

 ! Aerosol coefficient model, file, and format
 conf%Aerosol_Model = "CRTM"
 if (f_confOpts%get("Aerosol_Model",str)) conf%Aerosol_Model = str
 conf%AerosolCoeff_File = "AerosolCoeff.bin"
 if (f_confOpts%get("AerosolCoeff_File",str)) conf%AerosolCoeff_File = str
 conf%AerosolCoeff_Format = "Binary"
 if (f_confOpts%get("AerosolCoeff_Format",str)) conf%AerosolCoeff_Format = str

 ! Coefficient file prefixes
 IRwaterCoeff = "Nalli"
 if (f_confOpts%get("IRwaterCoeff",str)) IRwaterCoeff = str
 VISwaterCoeff = "NPOESS"
 if (f_confOpts%get("VISwaterCoeff",str)) VISwaterCoeff = str
 IRVISlandCoeff = "NPOESS"
 if (f_confOpts%get("IRVISlandCoeff",str)) IRVISlandCoeff = str
 IRVISsnowCoeff = "NPOESS"
 if (f_confOpts%get("IRVISsnowCoeff",str)) IRVISsnowCoeff = str
 IRVISiceCoeff = "NPOESS"
 if (f_confOpts%get("IRVISiceCoeff",str)) IRVISiceCoeff = str
 MWwaterCoeff = "FASTEM6"
 if (f_confOpts%get("MWwaterCoeff",str)) MWwaterCoeff = str

 ! Define water, snow, ice (WSI) categories
 select case (IRVISlandCoeff)
    case ("USGS")
       allocate(conf%Land_WSI(2))
       conf%Land_WSI(1:2) = [16,24]
    case ("IGBP")
       allocate(conf%Land_WSI(2))
       conf%Land_WSI(1:2) = [15,17]
    case ("NPOESS")
       allocate(conf%Land_WSI(1))
       conf%Land_WSI(1) = -1
    case default
       message = trim(routine_name) // " error: unknown IR/vis land coeff classification " // &
                 IRVISlandCoeff
       call fckit_exception % throw(message)
 end select

 ! IR emissivity coeff files
 conf%IRwaterCoeff_File = trim(IRwaterCoeff)//".IRwater.EmisCoeff.bin"
 conf%IRlandCoeff_File  = trim(IRVISlandCoeff)//".IRland.EmisCoeff.bin"
 conf%IRsnowCoeff_File  = trim(IRVISsnowCoeff)//".IRsnow.EmisCoeff.bin"
 conf%IRiceCoeff_File   = trim(IRVISiceCoeff)//".IRice.EmisCoeff.bin"

 !VIS emissivity coeff files
 conf%VISwaterCoeff_File = trim(VISwaterCoeff)//".VISwater.EmisCoeff.bin"
 conf%VISlandCoeff_File  = trim(IRVISlandCoeff)//".VISland.EmisCoeff.bin"
 conf%VISsnowCoeff_File  = trim(IRVISsnowCoeff)//".VISsnow.EmisCoeff.bin"
 conf%VISiceCoeff_File   = trim(IRVISiceCoeff)//".VISice.EmisCoeff.bin"

 ! MW water emissivity coeff file
 conf%MWwaterCoeff_File = trim(MWwaterCoeff)//".MWwater.EmisCoeff.bin"

 conf%inspect = 0
 if (f_confOpts%has("InspectProfileNumber")) then
   call f_confOpts%get_or_die("InspectProfileNumber",conf%inspect)
 end if

 if (f_confOpts%has("ShareCoefficients")) then
   call f_confOpts%get_or_die("ShareCoefficients",conf%share_coefficients)
 end if

 ! reconstruction operator cmatrix path
 conf % read_Cmatrix = .false.
 if (f_confOpts % has("ReconstructedRadianceCorrection")) then
    call f_confOpts % get_or_die("ReconstructedRadianceCorrection", conf % read_Cmatrix)
 end if

 if (conf % read_Cmatrix) then
    call f_confOpts % get_or_die("CMatrixPath", str)
    conf % Cmatrix_path = trim(str)
 end if

 call log_version()

end subroutine ufo_crtm_conf_setup

! ------------------------------------------------------------------------------

!> Log the CRTM library version, once per process (it is the same for every
!> CRTM operator).
subroutine log_version()

use fckit_log_module, only: fckit_log

implicit none

character(len=256) :: version

 if (crtm_version_logged) return
 call CRTM_Version(version)
 call fckit_log%info("UFO CRTM operators using CRTM version: "//trim(version))
 crtm_version_logged = .true.

end subroutine log_version

! ------------------------------------------------------------------------------

!> Throw if any entry of names appears more than once.
subroutine conf_check_duplicates(names, label, routine_name)

implicit none
character(len=*), intent(in) :: names(:)
character(len=*), intent(in) :: label
character(len=*), intent(in) :: routine_name

character(len=max_string) :: message
integer :: j

 do j = 1, size(names) - 1
   if (any(names(j) == names(j+1:))) then
     message = trim(routine_name) // " error: " // trim(names(j)) // " is duplicated in " // label
     call fckit_exception % throw(message)
   end if
 end do

end subroutine conf_check_duplicates

! ------------------------------------------------------------------------------

!> Read an optional cloud method option that may be "thompson" or "none"
!> (case-insensitive). Returns .true. for "thompson"; absent means "none".
logical function conf_thompson_option(f_conf, key, routine_name) result(use_thompson)

implicit none
type(fckit_configuration), intent(in) :: f_conf
character(len=*),          intent(in) :: key
character(len=*),          intent(in) :: routine_name

character(len=:), allocatable :: str
character(len=max_string) :: message

 use_thompson = .false.
 if (.not. f_conf%get(key, str)) return
 select case (upper2lower(str))
   case ("thompson")
     use_thompson = .true.
   case ("none")
   case default
     message = trim(routine_name) // ' error: "' // key // &
               '" can only be "thompson" or "none" (case-insensitive).'
     call fckit_exception % throw(message)
 end select

end function conf_thompson_option

! ------------------------------------------------------------------------------

!> Read the optional fixed Cloud_Fraction, validate it, and print the related
!> warnings. Requires conf%cal_cloud_frac_in_fov to be set already.
subroutine conf_cloud_fraction(conf, f_conf, routine_name, message_flag)

implicit none
type(crtm_conf),           intent(inout) :: conf
type(fckit_configuration), intent(in)    :: f_conf
character(len=*),          intent(in)    :: routine_name
logical,                   intent(in)    :: message_flag

! To save the value read from yaml, since Cloud_Fraction is derived from Obs Opearator and used for
! Linear Obs Operator too.
real(kind_real), save :: cloudFraction = -1.0_kind_real
character(len=max_string) :: message

 if (f_conf%has("Cloud_Fraction")) then
   if (conf%cal_cloud_frac_in_fov) then
     message = trim(routine_name) // &
         ": Cloud_Fraction input is ignored since it will be calculated within UFO."
     if (message_flag) call Display_Message(routine_name, trim(message), WARNING )
   else
     call f_conf%get_or_die("Cloud_Fraction",conf%Cloud_Fraction)
     cloudFraction = conf%Cloud_Fraction ! saving from Yaml file to be used for linear-obs-operator
     if ( conf%Cloud_Fraction < 0.0 .or. &
          conf%Cloud_Fraction > 1.0 ) then
       message = trim(routine_name) // " error: must specify" // &
                 " 0.0 <= Cloud_Fraction <= 1.0" // &
                 " or remove Cloud_Fraction from conf" // &
                 " and provide as a geoval"
       call fckit_exception % throw(message)
     end if
   end if
 else if (.not. conf%cal_cloud_frac_in_fov .and. cloudFraction <= 0.0) then
   ! If Cloud_Fraction is not provided as parameter in Obs Operator then provide warning
   message = trim(routine_name) // &
           ": Cloud_Fraction is not provided in conf." // &
           " Will request as a geoval."
   if (message_flag) call Display_Message(routine_name, trim(message), WARNING )
 end if

 if (conf%cal_cloud_frac_in_fov) then
   message = trim(routine_name) // &
         ': Using the "method for cloud fraction within fov: thompson"' // &
         " in areas where clouds are not zeroed-out."
   if (message_flag) call Display_Message(routine_name, trim(message), WARNING )
   if (conf%flag_deep_conv_mass_flux) then
     message = trim(routine_name) // &
         ": convection_mass_flux_flag is TRUE for " // &
         '"method for cloud fraction within fov: thompson".'
     if (message_flag) call Display_Message(routine_name, trim(message), WARNING )
   end if
 end if

end subroutine conf_cloud_fraction

! -----------------------------------------------------------------------------

subroutine ufo_crtm_conf_delete(conf)

implicit none
type(crtm_conf), intent(inout) :: conf

 deallocate(conf%SENSOR_ID)
 deallocate(conf%Land_WSI)
 deallocate(conf%Absorbers)
 deallocate(conf%Absorber_Id)
 deallocate(conf%Absorber_Units)
 deallocate(conf%Clouds)
 deallocate(conf%Cloud_Id)
 deallocate(conf%Aerosols)
 deallocate(conf%nullifyJacVars)
 deallocate(conf%nullifyJacChannels)

end subroutine ufo_crtm_conf_delete

! ------------------------------------------------------------------------------

!> Initialize CRTM for the aerosol operators (ufo_aodcrtm_mod, ufo_aodcrtm_tlad_mod
!> and ufo_extcoeffprofcrtm_mod): Spc/Tau and aerosol coefficients only, no
!> Cloud/IR/VIS/MW. The aerosol model is derived from conf%AerosolCoeff_File and
!> optionally returned in aerosol_model for the caller's later use.
!> Sensors registered with ufo_crtm_register_aerosol use the shared load.
subroutine ufo_crtm_init_aerosol(conf, chinfo, program_name, f_comm, aerosol_model)

implicit none
type(crtm_conf),             intent(in)    :: conf
type(CRTM_ChannelInfo_type), intent(inout) :: chinfo(:)
character(len=*),            intent(in)    :: program_name
type(fckit_mpi_comm),        intent(in)    :: f_comm
character(len=*), optional,  intent(out)   :: aerosol_model

character(len=MAXVARLEN)  :: def_aero_mod

 call ufo_crtm_define_aerosol_model(conf%AerosolCoeff_File, def_aero_mod)
 if (present(aerosol_model)) aerosol_model = def_aero_mod

 call init_coeffs(conf, aerosol_coeff_key(conf), .true., chinfo, program_name, f_comm)

end subroutine ufo_crtm_init_aerosol

! ------------------------------------------------------------------------------

!> Initialize CRTM for the full-physics radiance operators
!> the full coefficient set (Cloud/IR/VIS/MW) plus the per-sensor channel subset
!> Sensors registered with ufo_crtm_register_radiance use the shared load.
subroutine ufo_crtm_init_radiance(conf, chinfo, channels, program_name, f_comm)

implicit none
type(crtm_conf),             intent(in)    :: conf
type(CRTM_ChannelInfo_type), intent(inout) :: chinfo(:)
integer,                     intent(in)    :: channels(:)
character(len=*),            intent(in)    :: program_name
type(fckit_mpi_comm),        intent(in)    :: f_comm

character(len=max_string) :: message
integer :: err_stat, n

 call init_coeffs(conf, radiance_coeff_key(conf), .false., chinfo, program_name, f_comm)

 do n = 1, conf%n_Sensors
   err_stat = CRTM_ChannelInfo_Subset(chinfo(n), channels, reset=.false.)
   message = "Error subsetting channels!"
   call ufo_crtm_comm_stat_check(err_stat, program_name, message, f_comm)
 end do

end subroutine ufo_crtm_init_radiance

! ------------------------------------------------------------------------------

!> chinfo for conf's sensors. If they are registered under key, from the shared
!> load of their group (loaded here on first use, replacing any other group's
!> load); otherwise from a CRTM_Init of their own (freed by ufo_crtm_destroy).
subroutine init_coeffs(conf, key, aerosol, chinfo, program_name, f_comm)

implicit none
type(crtm_conf),             intent(in)    :: conf
character(len=*),            intent(in)    :: key
logical,                     intent(in)    :: aerosol
type(CRTM_ChannelInfo_type), intent(inout) :: chinfo(:)
character(len=*),            intent(in)    :: program_name
type(fckit_mpi_comm),        intent(in)    :: f_comm

integer :: g, n, k

 g = find_group(key, aerosol)
 if (g > 0) then
   do n = 1, conf%n_Sensors
     if (.not. any(shared_groups(g)%sensor_ids == conf%SENSOR_ID(n))) g = 0
   end do
 end if

 if (g == 0) then
   call drop_shared_load(program_name, f_comm)
   call crtm_init_coeffs(conf, aerosol, conf%SENSOR_ID, chinfo, program_name, f_comm)
   return
 end if

 if (loaded_group /= g) then
   call drop_shared_load(program_name, f_comm)
   allocate(shared_chinfo(size(shared_groups(g)%sensor_ids)))
   call crtm_init_coeffs(conf, aerosol, shared_groups(g)%sensor_ids, shared_chinfo, &
                         program_name, f_comm)
   loaded_group = g
   shared_loaded = .true.
 end if

 ! This operator's ChannelInfo is a copy of its sensor's entry in the shared load
 do n = 1, conf%n_Sensors
   k = findloc(shared_chinfo%Sensor_ID, conf%SENSOR_ID(n), 1)
   chinfo(n) = shared_chinfo(k)
 end do

end subroutine init_coeffs

! ------------------------------------------------------------------------------

!> CRTM_Init for sensor_ids with the radiance or aerosol coefficient settings of conf.
subroutine crtm_init_coeffs(conf, aerosol, sensor_ids, chinfo, program_name, f_comm)

implicit none
type(crtm_conf),             intent(in)    :: conf
logical,                     intent(in)    :: aerosol
character(len=*),            intent(in)    :: sensor_ids(:)
type(CRTM_ChannelInfo_type), intent(inout) :: chinfo(:)
character(len=*),            intent(in)    :: program_name
type(fckit_mpi_comm),        intent(in)    :: f_comm

character(len=MAXVARLEN)  :: def_aero_mod
character(len=max_string) :: message
integer :: err_stat

 if (aerosol) then
   call ufo_crtm_define_aerosol_model(conf%AerosolCoeff_File, def_aero_mod)
   err_stat = CRTM_Init( sensor_ids                                          , &
                         chinfo                                              , &
                         File_Path           = trim(conf%COEFFICIENT_PATH)   , &
                         NC_File_Path        = trim(conf%NC_COEFFICIENT_PATH), &
                         SpcCoeff_Format     = trim(conf%SpcCoeff_Format)    , &
                         TauCoeff_Format     = trim(conf%TauCoeff_Format)    , &
                         Aerosol_Model       = trim(def_aero_mod)            , &
                         AerosolCoeff_Format = trim(conf%AerosolCoeff_Format), &
                         AerosolCoeff_File   = trim(conf%AerosolCoeff_File)  , &
                         Quiet               = .true.)
 else
   err_stat = CRTM_Init( sensor_ids                                          , &
                         chinfo                                              , &
                         File_Path           = trim(conf%COEFFICIENT_PATH)   , &
                         NC_File_Path        = trim(conf%NC_COEFFICIENT_PATH), &
                         SpcCoeff_Format     = trim(conf%SpcCoeff_Format)    , &
                         TauCoeff_Format     = trim(conf%TauCoeff_Format)    , &
                         Aerosol_Model       = trim(conf%Aerosol_Model)      , &
                         AerosolCoeff_Format = trim(conf%AerosolCoeff_Format), &
                         AerosolCoeff_File   = trim(conf%AerosolCoeff_File)  , &
                         Cloud_Model         = trim(conf%Cloud_Model)        , &
                         CloudCoeff_Format   = trim(conf%CloudCoeff_Format)  , &
                         CloudCoeff_File     = trim(conf%CloudCoeff_File)    , &
                         IRwaterCoeff_File   = trim(conf%IRwaterCoeff_File)  , &
                         IRlandCoeff_File    = trim(conf%IRlandCoeff_File)   , &
                         IRsnowCoeff_File    = trim(conf%IRsnowCoeff_File)   , &
                         IRiceCoeff_File     = trim(conf%IRiceCoeff_File)    , &
                         VISwaterCoeff_File  = trim(conf%VISwaterCoeff_File) , &
                         VISlandCoeff_File   = trim(conf%VISlandCoeff_File)  , &
                         VISsnowCoeff_File   = trim(conf%VISsnowCoeff_File)  , &
                         VISiceCoeff_File    = trim(conf%VISiceCoeff_File)   , &
                         MWwaterCoeff_File   = trim(conf%MWwaterCoeff_File)  , &
                         Quiet               = .true.)
 end if
 message = "Error initializing CRTM"
 call ufo_crtm_comm_stat_check(err_stat, program_name, message, f_comm)

end subroutine crtm_init_coeffs

! ------------------------------------------------------------------------------

!> Register the sensors of a radiance operator for the shared CRTM load. Called at
!> operator setup.
subroutine ufo_crtm_register_radiance(conf)

implicit none
type(crtm_conf), intent(in) :: conf

 call register_sensors(conf, radiance_coeff_key(conf), .false.)

end subroutine ufo_crtm_register_radiance

! ------------------------------------------------------------------------------

!> Register the sensors of an aerosol operator for the shared CRTM load. Called at
!> operator setup.
subroutine ufo_crtm_register_aerosol(conf)

implicit none
type(crtm_conf), intent(in) :: conf

 call register_sensors(conf, aerosol_coeff_key(conf), .true.)

end subroutine ufo_crtm_register_aerosol

! ------------------------------------------------------------------------------

!> Undo ufo_crtm_register_radiance. Called at operator delete, before conf is deleted.
subroutine ufo_crtm_deregister_radiance(conf)

implicit none
type(crtm_conf), intent(in) :: conf

 call deregister_sensors(conf, radiance_coeff_key(conf), .false.)

end subroutine ufo_crtm_deregister_radiance

! ------------------------------------------------------------------------------

!> Undo ufo_crtm_register_aerosol. Called at operator delete, before conf is deleted.
subroutine ufo_crtm_deregister_aerosol(conf)

implicit none
type(crtm_conf), intent(in) :: conf

 call deregister_sensors(conf, aerosol_coeff_key(conf), .true.)

end subroutine ufo_crtm_deregister_aerosol

! ------------------------------------------------------------------------------

!> Abort unless conf's Absorbers include H2O and O3; with h2o_first, H2O must also be
!> the first absorber (as CRTM_AOD requires).
subroutine ufo_crtm_check_absorbers(conf, h2o_first, routine_name)

implicit none
type(crtm_conf),  intent(in) :: conf
logical,          intent(in) :: h2o_first
character(len=*), intent(in) :: routine_name

character(len=max_string) :: err_msg
integer :: ih2o, ioz

 ih2o = ufo_vars_getindex(conf%Absorbers, var_mixr)
 ioz  = ufo_vars_getindex(conf%Absorbers, var_oz)
 if (h2o_first .and. ih2o /= 1) then
   err_msg = routine_name // " error: H2O must be first in CRTM Absorbers for CRTM_AOD"
   call abor1_ftn(err_msg)
 else if (ih2o < 1) then
   err_msg = routine_name // " error: H2O must be included in CRTM Absorbers"
   call abor1_ftn(err_msg)
 end if
 if (ioz < 1) then
   err_msg = routine_name // " error: O3 must be included in CRTM Absorbers"
   call abor1_ftn(err_msg)
 end if

end subroutine ufo_crtm_check_absorbers

! ------------------------------------------------------------------------------

!> Land type GeoVaL of the classification (NPOESS, IGBP or USGS) that conf's
!> IRlandCoeff_File is named after.
function ufo_crtm_landtyp_var(conf, routine_name) result(varname)

implicit none
type(crtm_conf),  intent(in) :: conf
character(len=*), intent(in) :: routine_name
character(len=MAXVARLEN) :: varname

character(len=max_string) :: message

 if (index(conf%IRlandCoeff_File, "NPOESS") == 1) then
   ! this is the default if CRTM not configured otherwise
   varname = var_sfc_landtyp_npoess
 else if (index(conf%IRlandCoeff_File, "IGBP") == 1) then
   varname = var_sfc_landtyp_igbp
 else if (index(conf%IRlandCoeff_File, "USGS") == 1) then
   varname = var_sfc_landtyp_usgs
 else
   message = routine_name // " error: cannot infer land type classification from " // &
             "IRlandCoeff_File " // conf%IRlandCoeff_File
   call fckit_exception % throw(message)
 end if

end function ufo_crtm_landtyp_var

! ------------------------------------------------------------------------------

!> Add conf's sensors to the group for key (created if new), unless conf opts out
!> of shared loads. Adding a sensor to the group that is in memory makes its next
!> use reload the group.
subroutine register_sensors(conf, key, aerosol)

use fckit_log_module, only: fckit_log

implicit none
type(crtm_conf),  intent(in) :: conf
character(len=*), intent(in) :: key
logical,          intent(in) :: aerosol

type(shared_load_group) :: new_group
character(len=max_string) :: message
integer :: g, n, k

 if (.not. conf%share_coefficients) return

 if (.not. allocated(shared_groups)) allocate(shared_groups(0))
 g = find_group(key, aerosol)
 if (g == 0) then
   new_group%key = key
   new_group%aerosol = aerosol
   allocate(new_group%sensor_ids(0), new_group%n_users(0))
   shared_groups = [shared_groups, new_group]
   g = size(shared_groups)
 end if

 do n = 1, conf%n_Sensors
   k = findloc(shared_groups(g)%sensor_ids, conf%SENSOR_ID(n), 1)
   if (k > 0) then
     shared_groups(g)%n_users(k) = shared_groups(g)%n_users(k) + 1
     cycle
   end if
   if (size(shared_groups(g)%sensor_ids) == MAX_N_SENSORS) then
     write(message, "(3a,i0,a)") "ufo_crtm: cannot add ", trim(conf%SENSOR_ID(n)), &
       " to a shared CRTM load, which holds at most ", MAX_N_SENSORS, &
       " sensors; set ShareCoefficients: false in the obs options of some operators"
     call fckit_exception%throw(message)
   end if
   shared_groups(g)%sensor_ids = [shared_groups(g)%sensor_ids, conf%SENSOR_ID(n)]
   shared_groups(g)%n_users = [shared_groups(g)%n_users, 1]
   if (loaded_group == g) loaded_group = -g   ! in memory, but now missing a sensor
 end do

 if (.not. shared_groups_warned .and. count(group_sizes() > 0) > 1) then
   call fckit_log%warning("ufo_crtm: CRTM operators have different coefficient " // &
     "settings, so they use more than one shared CRTM load. CRTM holds one at a " // &
     "time, and switching between them reloads all of a load's coefficients. If " // &
     "this is slow, set ShareCoefficients: false in the obs options of the operators " // &
     "with uncommon settings.")
   shared_groups_warned = .true.
 end if

end subroutine register_sensors

! ------------------------------------------------------------------------------

!> Remove one user of conf's sensors from the group for key. A sensor with no users
!> left leaves the group; when the group is empty its load, if in memory, is destroyed.
subroutine deregister_sensors(conf, key, aerosol)

implicit none
type(crtm_conf),  intent(in) :: conf
character(len=*), intent(in) :: key
logical,          intent(in) :: aerosol

type(fckit_mpi_comm) :: f_comm
integer :: g, n, k

 if (.not. conf%share_coefficients) return
 g = find_group(key, aerosol)
 if (g == 0) return

 do n = 1, conf%n_Sensors
   k = findloc(shared_groups(g)%sensor_ids, conf%SENSOR_ID(n), 1)
   if (k == 0) cycle
   shared_groups(g)%n_users(k) = shared_groups(g)%n_users(k) - 1
   if (shared_groups(g)%n_users(k) > 0) cycle
   ! a load in memory keeps the sensor until the group is next loaded
   shared_groups(g)%sensor_ids = [shared_groups(g)%sensor_ids(:k-1), shared_groups(g)%sensor_ids(k+1:)]
   shared_groups(g)%n_users = [shared_groups(g)%n_users(:k-1), shared_groups(g)%n_users(k+1:)]
 end do

 if (size(shared_groups(g)%sensor_ids) == 0 .and. abs(loaded_group) == g) then
   f_comm = fckit_mpi_comm()
   call drop_shared_load("ufo_crtm_deregister", f_comm)
 end if

end subroutine deregister_sensors

! ------------------------------------------------------------------------------

!> Number of sensors registered in each group.
function group_sizes() result(sizes)

implicit none
integer, allocatable :: sizes(:)

integer :: g

 allocate(sizes(size(shared_groups)))
 do g = 1, size(shared_groups)
   sizes(g) = size(shared_groups(g)%sensor_ids)
 end do

end function group_sizes

! ------------------------------------------------------------------------------

!> Index of the group with this key and kind in shared_groups (0 if none).
integer function find_group(key, aerosol)

implicit none
character(len=*), intent(in) :: key
logical,          intent(in) :: aerosol

integer :: g

 find_group = 0
 if (.not. allocated(shared_groups)) return
 do g = 1, size(shared_groups)
   if ((shared_groups(g)%aerosol .eqv. aerosol) .and. shared_groups(g)%key == key) then
     find_group = g
     return
   end if
 end do

end function find_group

! ------------------------------------------------------------------------------

!> Coefficient settings that must match for radiance operators to share one load.
function radiance_coeff_key(conf) result(key)

implicit none
type(crtm_conf), intent(in) :: conf
character(len=:), allocatable :: key

 key = trim(conf%COEFFICIENT_PATH)    // "|" // trim(conf%NC_COEFFICIENT_PATH) // "|" // &
       trim(conf%SpcCoeff_Format)     // "|" // trim(conf%TauCoeff_Format)     // "|" // &
       trim(conf%Aerosol_Model)       // "|" // trim(conf%AerosolCoeff_Format) // "|" // &
       trim(conf%AerosolCoeff_File)   // "|" // trim(conf%Cloud_Model)         // "|" // &
       trim(conf%CloudCoeff_Format)   // "|" // trim(conf%CloudCoeff_File)     // "|" // &
       trim(conf%IRwaterCoeff_File)   // "|" // trim(conf%IRlandCoeff_File)    // "|" // &
       trim(conf%IRsnowCoeff_File)    // "|" // trim(conf%IRiceCoeff_File)     // "|" // &
       trim(conf%VISwaterCoeff_File)  // "|" // trim(conf%VISlandCoeff_File)   // "|" // &
       trim(conf%VISsnowCoeff_File)   // "|" // trim(conf%VISiceCoeff_File)    // "|" // &
       trim(conf%MWwaterCoeff_File)

end function radiance_coeff_key

! ------------------------------------------------------------------------------

!> Coefficient settings that must match for aerosol operators to share one load
!> (the aerosol model follows from AerosolCoeff_File).
function aerosol_coeff_key(conf) result(key)

implicit none
type(crtm_conf), intent(in) :: conf
character(len=:), allocatable :: key

 key = trim(conf%COEFFICIENT_PATH)    // "|" // trim(conf%NC_COEFFICIENT_PATH) // "|" // &
       trim(conf%SpcCoeff_Format)     // "|" // trim(conf%TauCoeff_Format)     // "|" // &
       trim(conf%AerosolCoeff_Format) // "|" // trim(conf%AerosolCoeff_File)

end function aerosol_coeff_key

! ------------------------------------------------------------------------------

!> Release the shared load, if it is in memory, so another CRTM_Init can be made.
subroutine drop_shared_load(program_name, f_comm)

implicit none
character(len=*),     intent(in) :: program_name
type(fckit_mpi_comm), intent(in) :: f_comm

character(len=max_string) :: message
integer :: err_stat

 ! shared_chinfo can be allocated without a finished load if CRTM_Init failed
 if (shared_loaded) then
   err_stat = CRTM_Destroy( shared_chinfo )
   message = "Error destroying shared CRTM load"
   call ufo_crtm_comm_stat_check(err_stat, program_name, message, f_comm)
 end if
 if (allocated(shared_chinfo)) deallocate(shared_chinfo)
 loaded_group = 0
 shared_loaded = .false.

end subroutine drop_shared_load

! ------------------------------------------------------------------------------

!> Destroy CRTM. Shared teardown for every CRTM operator. While the shared load is
!> in memory, chinfo is a copy of its entries, so only the copy is freed.
subroutine ufo_crtm_destroy(chinfo, program_name, f_comm)

implicit none
type(CRTM_ChannelInfo_type), intent(inout) :: chinfo(:)
character(len=*),            intent(in)    :: program_name
type(fckit_mpi_comm),        intent(in)    :: f_comm

character(len=max_string) :: message
integer :: err_stat

 if (shared_loaded) then
   call CRTM_ChannelInfo_Destroy( chinfo )
 else
   err_stat = CRTM_Destroy( chinfo )
   message = "Error destroying CRTM"
   call ufo_crtm_comm_stat_check(err_stat, program_name, message, f_comm)
 end if

end subroutine ufo_crtm_destroy

! ------------------------------------------------------------------------------

!> GMI channels 1-9 use the geometry geo and channels 10-13 the geometry geo_hf
!> (see ufo_crtm_load_geom_data). Split chinfo (already subset to channels, with
!> channels <= 9 listed first) into chinfo_gmi(1) = channels <= 9 and chinfo_gmi(2) =
!> channels >= 10, so that each block runs CRTM once with its own geometry. n_lo is the size
!> of the first block: CRTM output for it is (1:n_lo,:), for the second (n_lo+1:,:).
!> A block with no channels is left empty and must not be passed to CRTM.
subroutine ufo_crtm_gmi_split_chinfo(chinfo, channels, chinfo_gmi, n_lo, program_name, f_comm)

implicit none
type(CRTM_ChannelInfo_type), intent(in)    :: chinfo
integer,                     intent(in)    :: channels(:)
type(CRTM_ChannelInfo_type), intent(inout) :: chinfo_gmi(2)
integer,                     intent(out)   :: n_lo
character(len=*),            intent(in)    :: program_name
type(fckit_mpi_comm),        intent(in)    :: f_comm

character(len=max_string) :: message
integer :: err_stat

 n_lo = count(channels <= 9)
 if (any(channels(1:n_lo) > 9)) then
   message = "gmi_gpm channels 1-9 must be listed before channels 10-13"
   call ufo_crtm_comm_stat_check(FAILURE, program_name, message, f_comm)
 end if

 message = "Error subsetting gmi_gpm channels"
 if (n_lo > 0) then
   chinfo_gmi(1) = chinfo
   err_stat = CRTM_ChannelInfo_Subset(chinfo_gmi(1), channels(1:n_lo))
   call ufo_crtm_comm_stat_check(err_stat, program_name, message, f_comm)
 end if
 if (n_lo < size(channels)) then
   chinfo_gmi(2) = chinfo
   err_stat = CRTM_ChannelInfo_Subset(chinfo_gmi(2), channels(n_lo+1:))
   call ufo_crtm_comm_stat_check(err_stat, program_name, message, f_comm)
 end if

end subroutine ufo_crtm_gmi_split_chinfo

! ------------------------------------------------------------------------------

!> CRTM sensor type (MICROWAVE_SENSOR, INFRARED_SENSOR, ...) of sensor n of conf,
!> from its SpcCoeff file.
integer function ufo_crtm_sensor_type(conf, n)

use SpcCoeff_Define, only: SpcCoeff_type, SpcCoeff_Destroy

implicit none
type(crtm_conf), intent(in) :: conf
integer,         intent(in) :: n

type(SpcCoeff_type) :: spccoeff

 call read_spccoeff(conf%COEFFICIENT_PATH, conf%SENSOR_ID(n), conf%SpcCoeff_Format, spccoeff)
 ufo_crtm_sensor_type = spccoeff%Sensor_Type
 call SpcCoeff_Destroy(spccoeff)

end function ufo_crtm_sensor_type

! ------------------------------------------------------------------------------

!> Channel central wavelengths (nm) of sensor_id, from its binary SpcCoeff file.
subroutine ufo_crtm_sensor_wavelengths(coefficient_path, sensor_id, wavelengths)

use SpcCoeff_Define, only: SpcCoeff_type, SpcCoeff_Destroy

implicit none
character(len=*),             intent(in)  :: coefficient_path
character(len=*),             intent(in)  :: sensor_id
real(kind_real), allocatable, intent(out) :: wavelengths(:)

type(SpcCoeff_type) :: spccoeff

 call read_spccoeff(coefficient_path, sensor_id, "Binary", spccoeff)
 wavelengths = 1.0e7_kind_real / spccoeff%Wavenumber
 call SpcCoeff_Destroy(spccoeff)

end subroutine ufo_crtm_sensor_wavelengths

! ------------------------------------------------------------------------------

!> Read the SpcCoeff file of sensor_id (spccoeff_format "Binary" or "netCDF") into a
!> local structure. CRTM's loaded coefficients (SC) are left untouched, so this
!> needs no CRTM_Init.
subroutine read_spccoeff(coefficient_path, sensor_id, spccoeff_format, spccoeff)

use SpcCoeff_Define, only: SpcCoeff_type
use SpcCoeff_IO,     only: SpcCoeff_ReadFile

implicit none
character(len=*),    intent(in)    :: coefficient_path
character(len=*),    intent(in)    :: sensor_id
character(len=*),    intent(in)    :: spccoeff_format
type(SpcCoeff_type), intent(inout) :: spccoeff

character(len=:), allocatable :: path
character(len=max_string) :: err_msg
logical :: netcdf
integer :: err_stat

 ! same file name as CRTM_SpcCoeff_Load: <COEFFICIENT_PATH>/<Sensor_ID>.SpcCoeff.bin,
 ! or .SpcCoeff.nc for netCDF (also from COEFFICIENT_PATH, not NC_COEFFICIENT_PATH)
 netcdf = (spccoeff_format == "netCDF")
 path = trim(adjustl(coefficient_path))
 if (len(path) > 0) then
   if (path(len(path):len(path)) /= "/") path = path // "/"
 end if
 if (netcdf) then
   path = path // trim(adjustl(sensor_id)) // ".SpcCoeff.nc"
 else
   path = path // trim(adjustl(sensor_id)) // ".SpcCoeff.bin"
 end if
 err_stat = SpcCoeff_ReadFile(path, spccoeff, netCDF = netcdf, Quiet = .true.)
 if (err_stat /= SUCCESS) then
   err_msg = "ufo_crtm: failed to read SpcCoeff file " // path
   call fckit_exception%throw(err_msg)
 end if

end subroutine read_spccoeff


! ------------------------------------------------------------------------------

subroutine ufo_crtm_comm_stat_check(stat, PROGRAM_NAME, message, f_comm)

use fckit_mpi_module,   only: fckit_mpi_comm

implicit none

integer,              intent(in) :: stat
character(len=*),     intent(in) :: PROGRAM_NAME
character(len=*),     intent(in) :: message
type(fckit_mpi_comm), intent(in) :: f_comm

character(len=max_string) :: rank_message, fmessage

 if ( stat /= SUCCESS ) then
    write(rank_message,"(2a,i0)") trim(message), " on rank ", f_comm%rank()
    call Display_Message( PROGRAM_NAME, rank_message, FAILURE )
    fmessage = "Abort from "//PROGRAM_NAME
    call fckit_exception % throw(fmessage)
 end if

end subroutine ufo_crtm_comm_stat_check

! ------------------------------------------------------------------------------

subroutine ufo_crtm_skip_profiles(n_Profiles,n_Channels,channels,obss,atm,sfc,  &
                                  Is_Active_Sensor, Is_Vis_or_UV, Options)
! Profiles are skipped when the ObsValue of all channels is missing,
! if the pressure doesn't increase with levels, and
! if SST is missing for profiles with non-zero sea coverage.
! TODO: Use complete QC information
! It would be more comprehensive to use EffectiveQC or EffectiveError. That
! would require those ObsSpace values to be initialized before calls to
! this subroutine within ufo_radiancecrtm_simobs+ufo_radiancecrtm_tlad_settraj.
use missing_values_mod
use ufo_constants_mod, only: one

implicit none

integer,              intent(in)    :: n_Profiles, n_Channels
type(c_ptr), value,   intent(in)    :: obss
integer(c_int),       intent(in)    :: channels(:)
type(CRTM_Atmosphere_type), intent(in) :: atm(:)
type(CRTM_Surface_type),    intent(in) :: sfc(:)
logical, intent(in):: Is_Active_Sensor
logical, intent(in):: Is_Vis_or_UV
type(CRTM_Options_type),    intent(inout) :: Options(:)

integer :: jprofile, jchannel
character(len=MAXVARLEN) :: varname
character(len=64) :: obsGroupName
character(len=max_string) :: message
real(kind_real)  :: ObsVal(n_Profiles,n_Channels)
integer :: extendedObs(n_Profiles), actObsAvgQC(n_Profiles)
integer :: seqNum(n_Profiles)



real(c_double)  :: missing_d
real(kind_real) :: missing_r
real(kind_real), parameter  :: lowest_albedo = 0.001
real(kind_real), allocatable :: alon(:),alat(:)

 ! Set missing values
 missing_d = missing_value(missing_d)
 missing_r = missing_value(missing_r)

 ObsVal = missing_d

 ! Do a quick test to set the group name for ObsValue
 call ufo_crtm_get_var_name(channels(1),varname, Is_Active_Sensor, Is_Vis_or_UV)
 if (obsspace_has(obss, "DerivedObsValue", varname)) then
   obsGroupName = "DerivedObsValue"
 else if (obsspace_has(obss, "ObsValue", varname)) then
   obsGroupName = "ObsValue"
 else
   message = "Group name for observed values is neither ObsValue nor DerivedObsValue"
   call fckit_exception % throw(message)
 end if

 do jchannel = 1, n_Channels
   call ufo_crtm_get_var_name(channels(jchannel),varname, Is_Active_Sensor, Is_Vis_or_UV)
   call obsspace_get_db(obss, trim(obsGroupName), varname, ObsVal(:,jchannel))
end do

 if (Is_Active_Sensor) then
    call obsspace_get_db(obss, "MetaData", "sequenceNumber", seqNum)
    if (obsspace_has(obss, "MetaData", "extendedObsSpace")) then
       call obsspace_get_db(obss, "MetaData", "extendedObsSpace", extendedObs)
    else
       extendedObs = 1
    end if
    if (obsspace_has(obss, "MetaData", "actObsAvgQC")) then
       call obsspace_get_db(obss, "MetaData", "actObsAvgQC", actObsAvgQC)
    else
       actObsAvgQC = 1
    end if
 end if

 allocate(alon(n_Profiles),alat(n_Profiles))
 if (obsspace_has(obss, "MetaData", "latitude")) call obsspace_get_db( &
                                    obss, "MetaData", "latitude", alat)
 if (obsspace_has(obss, "MetaData", "longitude")) call obsspace_get_db( &
                                    obss, "MetaData", "longitude", alon)

 !Loop over all n_Profiles, i.e. number of locations.
 !Each check that skips a profile sets Skip_Profile and moves to the next profile.
 profile_loop: do jprofile = 1, n_Profiles
   Options(jprofile)%Skip_Profile = .false.

   ! check whether observations for all channels are missing in the input file
   if (all(ObsVal(jprofile,:) == missing_d)) then
     Options(jprofile)%Skip_Profile = .true.
     cycle profile_loop
   end if

   ! check for pressure monotonicity
   if (.not. ufo_crtm_pressure_is_monotonic(atm(jprofile))) then
     Options(jprofile)%Skip_Profile = .true.
     cycle profile_loop
   end if

   ! Sum of the coverage types can be found less than 1 near the boundary of regional
   ! configurations. CRTM needs them to sum to 1, so abort.
   if ( abs(sfc(jprofile)%Water_Coverage + sfc(jprofile)%Land_Coverage + &
        sfc(jprofile)%Snow_Coverage + sfc(jprofile)%Ice_Coverage - one ) > 1.0e-6_kind_real) then
    write(message,"(a,6(a,g0))") &
       "Abort! Surface coverage fractions do not sum to 1 for CRTM. Please check ob at ", &
       "Lon: ", alon(jprofile), ", Lat: ", alat(jprofile), &
       ", Water cover: ", sfc(jprofile)%Water_Coverage, &
       ", Land cover: ", sfc(jprofile)%Land_Coverage, &
       ", Snow cover: ", sfc(jprofile)%Snow_Coverage, &
       ", Ice cover: ", sfc(jprofile)%Ice_Coverage
    call fckit_exception % throw(message)
   end if

   ! check for missing values in water surface temperature when the mask
   ! indicates there is water.
   ! (this can happen with generically coupled atm and ocean components if
   ! the land/sea masks don't match)
   if ((sfc(jprofile)%Water_Temperature == missing_r) .and.   &
       (sfc(jprofile)%Water_Coverage > 0.0) .and. (.not. Is_Vis_or_UV)) then
     Options(jprofile)%Skip_Profile = .true.
     cycle profile_loop
   end if

   ! check for all channels in Vis/UV profiles that have ObsValue/albedo
   ! that are below minimum threshold. Skip those.
   if (Is_Vis_or_UV) then
      if (all(ObsVal(jprofile,:) < lowest_albedo)) then
         Options(jprofile)%Skip_Profile = .true.
         cycle profile_loop
      end if
   end if

   ! check for all channels in active profiles that have ObsValue/Reflectivity
   ! that are beyond threshold. Skip those.
   if (Is_Active_Sensor) then
      if (extendedObs(jprofile) == 0 .or. actObsAvgQC(jprofile) == 0) then
         Options(jprofile)%Skip_Profile = .true.
         cycle profile_loop
      end if
      ! the second dimension is for channels so if any channel is missing then skip it
      if (any(abs(ObsVal(jprofile,:)) >= threshold_reflectivity)) then
         Options(jprofile)%Skip_Profile = .true.
      end if
   end if
 end do profile_loop

end subroutine ufo_crtm_skip_profiles

! ------------------------------------------------------------------------------

!> Skip profiles for which no channel needs H(x): keep a profile only if at
!> least one of its channels has a QC flag of pass (0) or passive (1).
!> Profiles already skipped stay skipped.
subroutine ufo_crtm_skip_qc_profiles(qc_flags, n_channels, options, label)

use obsdatavector_mod,  only: obsdatavector_int
use fckit_log_module,   only: fckit_log
use, intrinsic :: iso_fortran_env, only: int64
implicit none
type(obsdatavector_int), intent(in)    :: qc_flags
integer,                 intent(in)    :: n_channels   !< number of simulated channels
type(CRTM_Options_type), intent(inout) :: options(:)   !< one per profile
character(len=*),        intent(in)    :: label        !< prefix for the debug message

character(len=max_string) :: dbg_msg
integer(int64) :: jchan, jprof
integer :: n_skipped

 do jprof = 1, size(options)
   if (options(jprof)%Skip_Profile) cycle
   ! qc_flags is not a Fortran array (and channels are not contiguous), so look up one by one
   options(jprof)%Skip_Profile = .true.
   do jchan = 1, n_channels
     if (qc_flags%get(jchan, jprof) < 2) then
       options(jprof)%Skip_Profile = .false.
       exit
     end if
   end do
 end do

 n_skipped = count(options(:)%Skip_Profile)
 write(dbg_msg,"(a,a,i9,a,i9,a)") trim(label), ", total of ", n_skipped, &
       " profiles being skipped. Using ", size(options) - n_skipped, " good profiles"
 call fckit_log%debug(dbg_msg)

end subroutine ufo_crtm_skip_qc_profiles

! ------------------------------------------------------------------------------

subroutine ufo_crtm_load_atm_data(n_Profiles, n_Layers, geovals, atm, conf, Is_Active_Sensor, &
                                  zeroCloudInCRTM, row_indices)

use ufo_constants_mod, only: zero
implicit none

integer, intent(in) :: n_Profiles, n_Layers
type(ufo_geovals), intent(in) :: geovals
type(CRTM_Atmosphere_type), intent(inout) :: atm(:)
logical, intent(in), optional :: Is_Active_Sensor
integer, intent(in), optional :: zeroCloudInCRTM(:)
! Optional map from profile index (1..n_Profiles) to the geovals column that
! profile's data actually lives in. Defaults to the identity (column k1 for
! profile k1) when absent, preserving the existing contiguous-range behavior
! for every caller that doesn't pass it.
integer, intent(in), optional :: row_indices(:)
type(crtm_conf), intent(in) :: conf

! Local variables
integer :: k1, jspec, jlevel
integer, allocatable :: src(:)
type(ufo_geoval), pointer :: geoval
type(ufo_geoval), pointer :: geoval_ts, geoval_prs, geoval_prsi, geoval_rh
logical :: apply_rh
character(max_string) :: err_msg
logical :: IsActiveSensor
logical :: seed_cloud
real(kind_real) :: geoval_unit_rescale
real(kind_real) :: co2
real(kind_real) :: seed_reff_min

real(kind_real), allocatable :: geoval_qsat(:,:), airdens(:,:)
real(kind_real), allocatable :: relative_humidity(:,:)
real(kind_real), allocatable :: rain_number(:,:), cloud_ice_number(:,:)
real(kind_real), allocatable :: clouds_mixingratio(:,:)
real(kind_real), allocatable :: cloudmxr_sum(:,:)
integer  :: id_cld(1)

  IsActiveSensor = .false.
  if (present(Is_Active_Sensor)) IsActiveSensor = Is_Active_Sensor

  allocate(src(n_Profiles))
  if (present(row_indices)) then
    src = row_indices
  else
    src = [(k1, k1=1,n_Profiles)]
  end if

  ! Populate the atmosphere structures for CRTM
  ! -------------------------------------------

  call ufo_geovals_get_var(geovals, var_ts, geoval_ts)
  ! Check model levels is consistent in geovals & crtm
  if (geoval_ts%nval /= n_Layers) then
    err_msg = "ufo_crtm_load_atm_data error: layers inconsistent!"
    call fckit_exception % throw(err_msg)
  end if
  call ufo_geovals_get_var(geovals, var_prs,  geoval_prs)
  call ufo_geovals_get_var(geovals, var_prsi, geoval_prsi)
  apply_rh = (conf%Aerosol_Model == "GOCART-GEOS5") .or. (conf%Aerosol_Model == "NAAPS")
  if (apply_rh) call ufo_geovals_get_var(geovals, var_rh, geoval_rh)

  do k1 = 1, n_Profiles
    atm(k1)%Temperature(1:n_Layers) = geoval_ts%vals(:, src(k1))
    atm(k1)%Pressure(1:n_Layers)    = geoval_prs%vals(:, src(k1)) * 0.01   ! to hPa
    atm(k1)%Level_Pressure(:)       = geoval_prsi%vals(:, src(k1)) * 0.01  ! to hPa
    atm(k1)%Climatology = US_STANDARD_ATMOSPHERE
    if (apply_rh) then
      geoval_rh%vals(:, src(k1)) = min(geoval_rh%vals(:, src(k1)), 1.0_kind_real)
      atm(k1)%Relative_Humidity(:) = geoval_rh%vals(:, src(k1))            ! fraction
    end if
  end do

  do jspec = 1, conf%n_Absorbers
    ! O3 Absorber has special treatment for Aerosols
    if (conf%Absorbers(jspec) == var_oz .and. &
      ufo_vars_getindex(geovals%variables, var_oz) < 0 .and. &
      (conf%aerosol_option /= "")) then
      do k1 = 1, n_Profiles
        atm(k1)%Absorber(1:n_Layers, jspec) = ozone_default_value
      end do
    else if (conf%Absorbers(jspec) == var_co2) then
      ! NOTE convert to ppmv
      select case (conf%co2_method_geovars)
        case ("MaunaLoa")
          call co2_mauna_loa_rise(conf%midpoint_julday,co2)
          do k1 = 1, n_Profiles
            atm(k1)%Absorber(1:n_Layers, jspec) = co2
          end do
        case ("EternalConstant")
          do k1 = 1, n_Profiles
            atm(k1)%Absorber(1:n_Layers, jspec) = conf%co2_ppmv_value
          end do
        case ("Background")
          ! this should automatically be filled by Model2GeoVals
          call ufo_geovals_get_var(geovals, conf%Absorbers(jspec), geoval)
          do k1 = 1, n_Profiles
            atm(k1)%Absorber(1:n_Layers, jspec) = co2_rescale_to_ppmv * geoval%vals(:, src(k1))
          end do
        case default
          call fckit_exception % throw("CO2 provided to CRTM is incorrectly set")
      end select
    else
      geoval_unit_rescale = one
      ! NOTE if "water_vapor_mixing_ratio_wrt_dry_air", convert from JEDI's kg/kg to CRTM's g/kg
      if (conf%Absorbers(jspec) == var_mixr) geoval_unit_rescale = kg_to_g
      call ufo_geovals_get_var(geovals, conf%Absorbers(jspec), geoval)
      do k1 = 1, n_Profiles
        atm(k1)%Absorber(1:n_Layers, jspec) = geoval_unit_rescale * geoval%vals(:, src(k1))
      end do
    end if
    do k1 = 1, n_Profiles
      atm(k1)%Absorber_Id(jspec)    = conf%Absorber_Id(jspec)
      atm(k1)%Absorber_Units(jspec) = conf%Absorber_Units(jspec)
    end do
  end do

  do jspec = 1, conf%n_Clouds
    ! cloud species content
    call ufo_geovals_get_var(geovals, conf%Clouds(jspec,1), geoval)
    ! writing the clamped value back to the geoval, presuming the value will be used later.
    geoval%vals(:, src) = max(geoval%vals(:, src), 0.0_kind_real)
    do k1 = 1, n_Profiles
      atm(k1)%Cloud(jspec)%Water_Content = geoval%vals(:, src(k1))
      atm(k1)%Cloud(jspec)%Type = conf%Cloud_Id(jspec)
    end do

    ! effective radius
    if (.not. conf%cal_cloud_reff_in_fov) then
      call ufo_geovals_get_var(geovals, conf%Clouds(jspec,2), geoval)
      geoval%vals(:, src) = max(geoval%vals(:, src), 0.0_kind_real)
      do k1 = 1, n_Profiles
        atm(k1)%Cloud(jspec)%Effective_Radius = geoval%vals(:, src(k1))
      end do
    end if
  end do

  ! When n_Clouds>0, Cloud_Fraction must either be provided as geoval or in conf
  ! If Cloud_Fraction is provided in conf,then use the Cloud_Fraction in conf
  ! If Cloud_Fraction is not provided in conf , then use the Cloud_Fraction in geoval
  ! and make sure the cloud fraction interpolated from background is in the valid range [0.0,1.0]
  if (conf%n_Clouds > 0 .and. (.not. conf%cal_cloud_frac_in_fov)) then
    if ( conf%Cloud_Fraction >= 0.0  ) then
      do k1 = 1, n_Profiles
        atm(k1)%Cloud_Fraction(:) = conf%Cloud_Fraction
      end do
    else if ( ufo_vars_getindex(geovals%variables, var_cldfrac) > 0 ) then
      call ufo_geovals_get_var(geovals, var_cldfrac, geoval)
      geoval%vals(:, src) = min(max(geoval%vals(:, src), 0.0_kind_real), &
                  1.0_kind_real)
      do k1 = 1, n_Profiles
        atm(k1)%Cloud_Fraction(:) = geoval%vals(:, src(k1))
      end do
    end if
  end if
  !
  if (conf%n_Clouds > 0 .and. (conf%cal_cloud_frac_in_fov .or. conf%cal_cloud_reff_in_fov)) then
    ! "saturation_water_vapor_mixing_ratio_wrt_moist_air"
    call ufo_geovals_get_var(geovals, var_qsat, geoval)
    geoval_qsat = geoval%vals(:, src)

    ! "moist_air_density"
    call ufo_geovals_get_var(geovals, var_airdens, geoval)
    airdens = geoval%vals(:, src)

    ! relative humidity from specific humidity, which can be
    ! different from atm(k1)%Relative_Humidity.
    call ufo_geovals_get_var(geovals, var_q, geoval)
    relative_humidity = max(geoval%vals(:, src), qsmall) / geoval_qsat

    call ufo_geovals_get_var(geovals, var_ni, geoval)
    cloud_ice_number = max(geoval%vals(:, src), 0.0_kind_real)

    call ufo_geovals_get_var(geovals, var_nr, geoval)
    rain_number = max(geoval%vals(:, src), 0.0_kind_real)

    ! cloud mixing ratios
    allocate(clouds_mixingratio(n_Layers, n_Profiles))
    allocate(cloudmxr_sum(n_Layers, n_Profiles), source=0.0_kind_real)
    jclouds_loop: do jspec = 1, conf%n_Clouds
      id_cld = findloc(thompson_Cloud_Id, conf%Cloud_Id(jspec))
      if (id_cld(1) > 0) then
        call ufo_geovals_get_var(geovals, trim(thompson_Cloud_var(id_cld(1))), geoval)
        clouds_mixingratio = max(geoval%vals(:, src), 0.0_kind_real)
        if (conf%cal_cloud_frac_in_fov) cloudmxr_sum = cloudmxr_sum + clouds_mixingratio
        if (conf%cal_cloud_reff_in_fov) then
          profile_loop_reff: do k1 = 1, n_Profiles
            ! skip 'call calc_thompson_reff' if clouds are zeroed-out.
            if (cloud_zeroed(k1)) cycle profile_loop_reff
            call calc_thompson_reff(airdens(:,k1),atm(k1)%Temperature,clouds_mixingratio(:,k1), &
                  conf%Cloud_Id(jspec), cloud_ice_number(:,k1), rain_number(:,k1), &
                  n_Layers,atm(k1)%Cloud(jspec)%Effective_Radius)
          end do profile_loop_reff
        end if
      end if
    end do jclouds_loop
    if (conf%cal_cloud_frac_in_fov) then
      profile_loop_cloudfrac: do k1 = 1, n_Profiles
        ! skip 'call calc_thompson_cloudfrac' if clouds are zeroed-out.
        if (cloud_zeroed(k1)) cycle profile_loop_cloudfrac
        call calc_thompson_cloudfrac(n_Layers, conf%flag_deep_conv_mass_flux, &
                atm(k1)%Pressure * 0.1_kind_real, &   ! to kPa
                cloudmxr_sum(:,k1), relative_humidity(:,k1), geoval_qsat(:, k1), &
                atm(k1)%Cloud_Fraction(:))
      end do profile_loop_cloudfrac
    end if
  end if
  !
  ! Set Water_Content, Effective_Radius and Cloud_Fraction = 0 when zeroCloudInCRTM = 1 (true)
  do k1 = 1, n_Profiles
    if (.not. cloud_zeroed(k1)) cycle
    do jspec = 1, conf%n_Clouds
      atm(k1)%Cloud(jspec)%Water_Content = zero
      atm(k1)%Cloud(jspec)%Effective_Radius = zero
    end do
    atm(k1)%Cloud_Fraction(:) = zero
  end do

  if ( (conf%n_Clouds > 0) .and. conf%precip_hydro) then
    do k1 = 1, n_Profiles
      ! set Effective_Radius = zero where Water_Content <= 1.0e-6
      do jspec = 1, conf%n_Clouds
         where(atm(k1)%Cloud(jspec)%Water_Content <= 1.0e-6_kind_real) &
                         atm(k1)%Cloud(jspec)%Effective_Radius = zero
      end do
      ! set minimum cloud_fraction where any species has Water_Content > 1.0e-6
      call set_min_cloud_fraction(k1, 1.0e-6_kind_real, .false., 1.0e-12_kind_real, &
                                  cloudfractionsmall)
    end do
  end if

  if ( (conf%n_Clouds > 0) .and. (.not. IsActiveSensor) .and. conf%Cloud_Seeding ) then
    profile_loop_cs: do k1 = 1, n_Profiles
      ! Do not check and reset these cloud values if clouds are zero-ed out.
      if (cloud_zeroed(k1)) cycle profile_loop_cs
      do jlevel = 1, atm(k1)%n_layers
         ! Check Cloud Content: seed liquid species above -20C and frozen species below 0C,
         ! with a species-dependent minimum effective radius [micron]
         do jspec = 1, conf%n_Clouds
           associate (cld => atm(k1)%Cloud(jspec), tc => atm(k1)%Temperature(jlevel) - t0c)
           select case (cld%Type)
             case (WATER_CLOUD)
               seed_cloud = tc > -20.0_kind_real;  seed_reff_min = 5.001_kind_real
             case (ICE_CLOUD)
               seed_cloud = tc < 0.0_kind_real;    seed_reff_min = 5.001_kind_real
             case (RAIN_CLOUD)
               seed_cloud = tc > -20.0_kind_real;  seed_reff_min = 100.001_kind_real
             case (SNOW_CLOUD)
               seed_cloud = tc < 0.0_kind_real;    seed_reff_min = 50.001_kind_real
             case (GRAUPEL_CLOUD)
               seed_cloud = tc < 0.0_kind_real;    seed_reff_min = 500.001_kind_real
             case (HAIL_CLOUD)
               seed_cloud = tc < 0.0_kind_real;    seed_reff_min = 1000.001_kind_real
             case default
               seed_cloud = .false.
           end select
           if (seed_cloud) then
             cld%Water_Content(jlevel) = max(cld%Water_Content(jlevel), &
                               WATER_CONTENT_THRESHOLD*1.001_kind_real)
             cld%Effective_Radius(jlevel) = max(cld%Effective_Radius(jlevel), seed_reff_min)
           end if
           end associate
         end do
      end do
      ! Check Cloud Fraction: at least 0.001 where any species is seeded
      call set_min_cloud_fraction(k1, WATER_CONTENT_THRESHOLD*1.001_kind_real, .true., &
                                  0.001_kind_real, 0.001_kind_real)
    end do profile_loop_cs
  end if

  if (IsActiveSensor) Atm%Add_Extra_Layers = .false.

  deallocate(src)

contains

  !> True if clouds are zeroed out for profile k (zeroCloudInCRTM(k) == 1).
  logical function cloud_zeroed(k)
    integer, intent(in) :: k
    cloud_zeroed = .false.
    if (present(zeroCloudInCRTM)) cloud_zeroed = (zeroCloudInCRTM(k) == 1)
  end function cloud_zeroed

  !> Set Cloud_Fraction = cf_new on levels of profile k where Cloud_Fraction < cf_below
  !> and any cloud species has Water_Content above wc_min (or >= wc_min if inclusive).
  subroutine set_min_cloud_fraction(k, wc_min, inclusive, cf_below, cf_new)
    integer,         intent(in) :: k
    real(kind_real), intent(in) :: wc_min, cf_below, cf_new
    logical,         intent(in) :: inclusive
    real(kind_real) :: wc_max(atm(k)%n_layers)
    integer :: j, nl
    nl = atm(k)%n_layers
    wc_max = atm(k)%Cloud(1)%Water_Content(1:nl)
    do j = 2, conf%n_Clouds
      wc_max = max(wc_max, atm(k)%Cloud(j)%Water_Content(1:nl))
    end do
    if (inclusive) then
      where (atm(k)%Cloud_Fraction(1:nl) < cf_below &
             .and. wc_max >= wc_min) atm(k)%Cloud_Fraction(1:nl) = cf_new
    else
      where (atm(k)%Cloud_Fraction(1:nl) < cf_below &
             .and. wc_max >  wc_min) atm(k)%Cloud_Fraction(1:nl) = cf_new
    end if
  end subroutine set_min_cloud_fraction

end subroutine ufo_crtm_load_atm_data

! ------------------------------------------------------------------------------

subroutine ufo_crtm_load_sfc_data(n_Profiles, n_Channels, channels, geovals, sfc, chinfo, obss, &
                                  conf, Is_Active_Sensor, Is_Vis_or_UV, zeroCloudInCRTM)

implicit none

integer,                     intent(in)    :: n_Profiles, n_Channels
type(ufo_geovals),           intent(in)    :: geovals
type(CRTM_Surface_type),     intent(inout) :: sfc(:)
type(CRTM_ChannelInfo_type), intent(in)    :: chinfo(:)
type(c_ptr), value,          intent(in)    :: obss
integer(c_int),              intent(in)    :: channels(:)
type(crtm_conf),             intent(in)    :: conf
integer,                     intent(out)   :: zeroCloudInCRTM(n_Profiles)
logical, intent(in) :: Is_Active_Sensor
logical, intent(in) :: Is_Vis_or_UV

integer :: k1, n1
integer :: iLand
logical :: use_mw_vegtyp_soiltyp_data, use_visir_landtyp_data
real(kind_real) :: u(n_Profiles), v(n_Profiles), land_type(n_Profiles)

character(len=MAXVARLEN) :: varname, landtyp_var
character(len=max_string) :: message
character(len=64) :: obsGroupName

real(kind_real), allocatable :: ObsTb(:,:)

  ! Do a quick test to set the group name for ObsValue
  call ufo_crtm_get_var_name(channels(1),varname, Is_Active_Sensor, Is_Vis_or_UV)
  if (obsspace_has(obss, "ObsValue", varname)) then
    obsGroupName = "ObsValue"
  else if (obsspace_has(obss, "DerivedObsValue", varname)) then
    obsGroupName = "DerivedObsValue"
  else
    message = "Group name for observed values is neither ObsValue nor DerivedObsValue"
    call fckit_exception % throw(message)
  end if

  !Pass observation value
  if (.not. Is_Active_Sensor) then
    allocate(ObsTb(n_Profiles, n_Channels))
    do n1 = 1, n_Channels
      call ufo_crtm_get_var_name(channels(n1),varname, Is_Active_Sensor, Is_Vis_or_UV)
      call obsspace_get_db(obss, trim(obsGroupName), varname, ObsTb(:, n1))
    end do
  end if

  do k1 = 1, n_Profiles
    !Pass sensor information
    sfc(k1)%sensordata%sensor_id        = chinfo(1)%sensor_id
    sfc(k1)%sensordata%wmo_sensor_id    = chinfo(1)%wmo_sensor_id
    sfc(k1)%sensordata%wmo_satellite_id = chinfo(1)%wmo_satellite_id
    sfc(k1)%sensordata%sensor_channel   = channels
    if (.not. Is_Active_Sensor) sfc(k1)%sensordata%tb(1:n_Channels) = ObsTb(k1, :)
  end do

  !Water_type
  !** NOTE: need to check how to determine fresh vs sea water types (salinity???)
  sfc(1:n_Profiles)%Water_Type = 1  ! SEA_WATER_TYPE for all SfcOptics

  use_mw_vegtyp_soiltyp_data = any(chinfo%Sensor_Type == MICROWAVE_SENSOR)
  use_visir_landtyp_data = any(chinfo%Sensor_Type == VISIBLE_SENSOR) .or. &
                           any(chinfo%Sensor_Type == INFRARED_SENSOR)

  if (ufo_vars_getindex(geovals%variables, var_sfc_wspeed) > 0 .and. &
      ufo_vars_getindex(geovals%variables, var_sfc_wdir) > 0) then
    ! Directly use model-provided wind speed and direction
    sfc(1:n_Profiles)%Wind_Speed     = sfc_geoval(var_sfc_wspeed)
    sfc(1:n_Profiles)%Wind_Direction = sfc_geoval(var_sfc_wdir)
  else if (ufo_vars_getindex(geovals%variables, var_sfc_u) > 0 .and. &
      ufo_vars_getindex(geovals%variables, var_sfc_v) > 0) then
    ! Convert 2d wind components to speed and direction
    u = sfc_geoval(var_sfc_u)
    v = sfc_geoval(var_sfc_v)
    sfc(1:n_Profiles)%Wind_Speed = sqrt(u**2 + v**2)
    do k1 = 1, n_Profiles
      sfc(k1)%Wind_Direction = uv_to_wdir(u(k1), v(k1))
    end do
  else
    call fckit_exception % throw("ufo_crtm_load_sfc_data error: missing surface wind geovals")
  end if

  sfc(1:n_Profiles)%Water_Coverage    = sfc_geoval(var_sfc_wfrac)
  sfc(1:n_Profiles)%Water_Temperature = sfc_geoval(var_sfc_wtmp)
  sfc(1:n_Profiles)%Ice_Coverage      = sfc_geoval(var_sfc_ifrac)
  sfc(1:n_Profiles)%Ice_Temperature   = sfc_geoval(var_sfc_itmp)
  sfc(1:n_Profiles)%Snow_Coverage     = sfc_geoval(var_sfc_sfrac)
  sfc(1:n_Profiles)%Snow_Temperature  = sfc_geoval(var_sfc_stmp)
  sfc(1:n_Profiles)%Snow_Depth        = sfc_geoval(var_sfc_sdepth)
  sfc(1:n_Profiles)%Land_Coverage     = sfc_geoval(var_sfc_lfrac)

  !Land_Type
  ! + used to lookup land sfc emiss. for IR and VIS
  ! + land sfc emiss. undefined over water/snow/ice
  ! + note that land type classifiation is an option
  if (use_visir_landtyp_data) then
    landtyp_var = ufo_crtm_landtyp_var(conf, "ufo_crtm_load_sfc_data")
    land_type = sfc_geoval(landtyp_var)
    do k1 = 1, n_Profiles
      iLand = int(land_type(k1))
      ! skip water/snow/ice (WSI) categories of the classification
      if (.not.any(iLand == conf%Land_WSI)) sfc(k1)%Land_Type = iLand
    end do
  end if

  sfc(1:n_Profiles)%Land_Temperature    = sfc_geoval(var_sfc_ltmp)
  sfc(1:n_Profiles)%Lai                 = sfc_geoval(var_sfc_lai)
  sfc(1:n_Profiles)%Vegetation_Fraction = sfc_geoval(var_sfc_vegfrac)

  if (use_mw_vegtyp_soiltyp_data) then
    sfc(1:n_Profiles)%Vegetation_Type = nint(sfc_geoval(var_sfc_vegtyp))
    sfc(1:n_Profiles)%Soil_Type       = nint(sfc_geoval(var_sfc_soiltyp))
  end if

  sfc(1:n_Profiles)%Soil_Moisture_Content = sfc_geoval(var_sfc_soilm)
  sfc(1:n_Profiles)%Soil_Temperature      = sfc_geoval(var_sfc_soilt)

  !Sea_Surface_Salinity
  if (conf%salinity_option == "on") sfc(1:n_Profiles)%Salinity = sfc_geoval(var_sfc_sss)

  ! Update Ice_Coverage and Land_Coverage for
  ! Soil_Type or Vegetation_Type corresponding to Glacial land ice
  where (sfc(1:n_Profiles)%Land_Coverage > ZERO .and. &
         (sfc(1:n_Profiles)%Soil_Type == 9 .or. sfc(1:n_Profiles)%Vegetation_Type == 13))
    sfc(1:n_Profiles)%Ice_Coverage  = min(sfc(1:n_Profiles)%Ice_Coverage &
        + sfc(1:n_Profiles)%Land_Coverage, ONE)
    sfc(1:n_Profiles)%Land_Coverage = ZERO
  end where

  if (obsspace_has(obss, "MetaData", "zeroCloudInCRTM")) then
    ! zeroCloudInCRTM is 1 over surface where clear-sky assimilation is conducted.
    call obsspace_get_db(obss, "MetaData", "zeroCloudInCRTM", zeroCloudInCRTM)
  else
    ! Do not zero-out clouds by default.
    zeroCloudInCRTM = 0
  end if

contains

  !> Values of the single-level surface GeoVaL `varname` for all profiles.
  function sfc_geoval(varname) result(vals)
    character(len=*), intent(in) :: varname
    real(kind_real) :: vals(n_Profiles)
    type(ufo_geoval), pointer :: geoval
    call ufo_geovals_get_var(geovals, varname, geoval)
    vals = geoval%vals(1, 1:n_Profiles)
  end function sfc_geoval

end subroutine ufo_crtm_load_sfc_data

! ------------------------------------------------------------------------------

subroutine ufo_crtm_load_geom_data(obss,geo,geo_hf,sensor_id)

use fckit_log_module, only: fckit_log

implicit none
type(c_ptr), value,       intent(in)    :: obss
type(CRTM_Geometry_type), intent(inout) :: geo(:)
type(CRTM_Geometry_type), intent(inout), optional :: geo_hf(:)
character(len=255),       intent(in),    optional :: sensor_id
integer, allocatable :: scan_position(:)
integer :: nlocs
character(len=:), allocatable :: missing_hf  ! channel 10-13 angles not in MetaData

 nlocs = obsspace_get_nlocs(obss)

 geo(:)%Sensor_Zenith_Angle  = abs(metadata("sensorZenithAngle")) ! needs to be absolute value
 geo(:)%Source_Zenith_Angle  = metadata("solarZenithAngle")
 geo(:)%Sensor_Azimuth_Angle = metadata("sensorAzimuthAngle")
 geo(:)%Source_Azimuth_Angle = metadata("solarAzimuthAngle")
 geo(:)%Sensor_Scan_Angle    = metadata("sensorViewAngle") !The Sensor_Scan_Angle is optional

 allocate(scan_position(nlocs))
 call obsspace_get_db(obss, "MetaData", "sensorScanPosition", scan_position)
 geo(:)%Ifov = scan_position(:)

 call clamp_angles(geo)

 ! Read geophysical values for gmi high frequency channels 10-13.
 if (present(sensor_id) .and. present(geo_hf)) then
   if (sensor_id == "gmi_gpm") then
     ! Angles missing from MetaData keep the channel 1-9 value
     geo_hf = geo
     missing_hf = ""
     ! needs to be absolute value
     if (has_hf("sensorZenithAngle1")) geo_hf(:)%Sensor_Zenith_Angle  = abs(metadata("sensorZenithAngle1"))
     if (has_hf("solarZenithAngle1"))  geo_hf(:)%Source_Zenith_Angle  = metadata("solarZenithAngle1")
     if (has_hf("sensorAzimuthAngle1")) geo_hf(:)%Sensor_Azimuth_Angle = metadata("sensorAzimuthAngle1")
     if (has_hf("solarAzimuthAngle1")) geo_hf(:)%Source_Azimuth_Angle = metadata("solarAzimuthAngle1")
     if (has_hf("sensorViewAngle1"))   geo_hf(:)%Sensor_Scan_Angle    = metadata("sensorViewAngle1")
     call clamp_angles(geo_hf)
     if (len(missing_hf) > 0 .and. .not. gmi_hf_geometry_warned) then
       call fckit_log%warning("ufo_crtm_load_geom_data: gmi_gpm MetaData missing:" // missing_hf // &
                              "; channels 10-13 use the channel 1-9 values for these angles")
       gmi_hf_geometry_warned = .true.
     end if
   end if
 end if

contains

  !> Whether MetaData/varname exists; if not, add it to missing_hf.
  logical function has_hf(varname)
    character(len=*), intent(in) :: varname
    has_hf = obsspace_has(obss, "MetaData", varname)
    if (.not. has_hf) missing_hf = missing_hf // " " // varname
  end function has_hf

  !> MetaData/varname for all locations.
  function metadata(varname) result(vals)
    character(len=*), intent(in) :: varname
    real(kind_real) :: vals(nlocs)
    call obsspace_get_db(obss, "MetaData", varname, vals)
  end function metadata

  !> Bring geometry angles into the ranges CRTM accepts.
  subroutine clamp_angles(g)
    type(CRTM_Geometry_type), intent(inout) :: g(:)
    g(:)%Sensor_Zenith_Angle = min(g(:)%Sensor_Zenith_Angle, MAX_SENSOR_ZENITH_ANGLE)
    ! For some microwave instruments the solar and sensor azimuth angles can be
    ! missing  (given a value of 10^11).  Set these to zero to get past CRTM QC.
    where (g(:)%Source_Azimuth_Angle < 0.0_kind_real .or. &
           g(:)%Source_Azimuth_Angle > 360.0_kind_real) &
       g(:)%Source_Azimuth_Angle = 0.0_kind_real
    where (g(:)%Sensor_Azimuth_Angle < 0.0_kind_real .or. &
           g(:)%Sensor_Azimuth_Angle > 360.0_kind_real) &
       g(:)%Sensor_Azimuth_Angle = 0.0_kind_real
    where (abs(g(:)%Source_Zenith_Angle) > 180.0_kind_real) &
       g(:)%Source_Zenith_Angle = 100.0_kind_real
    where (abs(g(:)%Sensor_Scan_Angle) > 80.0_kind_real) &
       g(:)%Sensor_Scan_Angle = 0.0_kind_real
  end subroutine clamp_angles

end subroutine ufo_crtm_load_geom_data

! ------------------------------------------------------------------------------

subroutine ufo_crtm_get_var_name(n,varname, Is_Active_Sensor, Is_Vis_or_UV)

integer, intent(in) :: n
character(len=*), intent(out) :: varname
logical, intent(in) :: Is_Active_Sensor
logical, intent(in) :: Is_Vis_or_UV

character(len=6) :: chan

 write(chan, "(I0)") n
 if (Is_Active_Sensor) then
     varname = "ReflectivityAttenuated_" // chan
 else
     if (Is_Vis_or_UV) then
        varname = "albedo_" // chan
     else
        varname = "brightnessTemperature_" // chan
     end if
 end if

end subroutine ufo_crtm_get_var_name

! -----------------------------------------------------------------------------

!> \brief Determines the wind direction from U and V components
!!
!! \details **uv_to_wdir** Calculates the wind direction, as measured clockwise
!! from north, similar to an azimuth angle.  Takes the eastward and northward
!! wind component magnitudes, respectively, as arguments.
!! Due to the azimuthal convention used here, the inverse equations are:
!! u = w * sin(wdir * deg2rad)
!! v = w * cos(wdir * deg2rad)
!! where w is the wind speed
function uv_to_wdir(u, v) result(wdir)

use ufo_constants_mod, only: zero, one, two, pi, rad2deg

implicit none

real (kind=kind_real), intent(in) :: u !< eastward_wind
real (kind=kind_real), intent(in) :: v !< northward_wind
real (kind=kind_real) :: wdir
real (kind=kind_real) :: windratio, windangle
integer :: iquadrant
real(kind=kind_real),parameter:: windscale = 999999.0_kind_real
real(kind=kind_real),parameter:: windlimit = 0.0001_kind_real
real(kind=kind_real),parameter:: quadcof(4,2) = &
  reshape([zero,  one,  one,  two, &
            one,  -one,  one, -one], [4,2])

  if (u >= zero .and. v >= zero) iquadrant = 1
  if (u >= zero .and. v <  zero) iquadrant = 2
  if (u <  zero .and. v >= zero) iquadrant = 4
  if (u <  zero .and. v <  zero) iquadrant = 3

  if (abs(v) >= windlimit) then
    windratio = u / v
  else
    windratio = zero
    if (abs(u) > windlimit) then
      windratio = windscale * u
    end if
  end if
  windangle = atan(abs(windratio))   ! wind azimuth is in radians
  wdir = ( quadcof(iquadrant, 1) * pi + windangle * quadcof(iquadrant, 2) ) * rad2deg

end function uv_to_wdir

! ------------------------------------------------------------------------------

subroutine ufo_crtm_define_aerosol_model(aerosol_coef_file, aerosol_model)
 ! Possible models are:
 ! CRTM, NAAPS, GOCART-GEOS5, CMAQ
 character(len=*), intent(in) :: aerosol_coef_file
 character(len=*), intent(out) :: aerosol_model

 if (aerosol_coef_file == "AerosolCoeff.nc4" .or. &
     aerosol_coef_file == "AerosolCoeff.bin") then
    aerosol_model = "CRTM"
 else if (aerosol_coef_file == "AerosolCoeff.NAAPS.nc4" .or. &
          aerosol_coef_file == "AerosolCoeff.NAAPS.bin") then
    aerosol_model = "NAAPS"
 else if (aerosol_coef_file == "AerosolCoeff.GOCART-GEOS5.nc4" .or. &
          aerosol_coef_file == "AerosolCoeff.GOCART-GEOS5.bin") then
    aerosol_model = "GOCART-GEOS5"
 else if (aerosol_coef_file == "AerosolCoeff.CMAQ.nc4" .or. &
          aerosol_coef_file == "AerosolCoeff.CMAQ.bin") then
    aerosol_model = "CMAQ"
 else
    call fckit_exception%throw("ufo_crtm_define_aerosol_model: cannot infer the aerosol " // &
                               "model from AerosolCoeff_File " // trim(aerosol_coef_file))
 end if

end subroutine ufo_crtm_define_aerosol_model

! ------------------------------------------------------------------------------

subroutine ufo_crtm_assign_aerosol_names(aerosol_option, var_aerosols)

 character(len=*), intent(in) :: aerosol_option
 character(len=MAXVARLEN), allocatable, intent(out) :: var_aerosols(:)

 character(max_string) :: err_msg

 if (aerosol_option == "aerosols_gocart_default") then
    allocate(var_aerosols(n_aerosols_gocart_default))
    var_aerosols=var_aerosols_gocart_default
 else if (aerosol_option == "aerosols_gocart_gefs") then
    allocate(var_aerosols(n_aerosols_gocart_gefs))
    var_aerosols=var_aerosols_gocart_gefs
 else if (aerosol_option == "aerosols_gocart_ufs") then
    allocate(var_aerosols(n_aerosols_gocart_ufs))
    var_aerosols=var_aerosols_gocart_ufs
 else if (aerosol_option == "aerosols_gocart_geos") then
    allocate(var_aerosols(n_aerosols_gocart_geos))
    var_aerosols=var_aerosols_gocart_geos
 else
    err_msg = "ufo_crtm_assign_aerosol_names: aerosol_option not implemented: " // &
              aerosol_option
    call fckit_exception % throw(err_msg)
  end if

end subroutine ufo_crtm_assign_aerosol_names

! ------------------------------------------------------------------------------

subroutine ufo_crtm_load_aerosol_data(n_profiles, n_layers, geovals, &
                                      conf, var_aerosols, aerosol_model, atm, row_indices)
 use CRTM_aerosolcoeff, only: aeroc
 use fckit_log_module, only: fckit_log

 type(crtm_conf), intent(in)    :: conf
 type(ufo_geovals), intent(in) :: geovals
 type(CRTM_atmosphere_type), intent(inout) :: atm(:)
 type(ufo_geoval), pointer :: geoval

 integer, intent(in) :: n_profiles, n_layers
 ! Optional map from profile index (1..n_profiles) to the geovals column
 ! that profile's data actually lives in. Defaults to the identity when
 ! absent, preserving existing behavior for every caller that omits it.
 integer, intent(in), optional :: row_indices(:)
 integer, allocatable :: src(:)
 integer :: n_aerosols, i, k, m

 real(kind_real), dimension(5), parameter  :: dust_radii = &
      [0.55_kind_real, 1.4_kind_real, 2.4_kind_real, 4.5_kind_real, 8.0_kind_real]
 real(kind_real), dimension(n_layers) :: layer_factors
 real(kind_real), dimension(n_layers, n_profiles) :: rh

 character(len=*), intent(in) :: aerosol_model
 character(len=MAXVARLEN), intent(in) :: var_aerosols(:)
 character(len=MAXVARLEN) :: varname
 character(max_string) :: err_msg

 allocate(src(n_profiles))
 if (present(row_indices)) then
   src = row_indices
 else
   src = [(m, m=1,n_profiles)]
 end if

 varname = var_rh
 call ufo_geovals_get_var(geovals, varname, geoval)
 do m = 1, n_profiles
   rh(1:n_layers,m) = geoval%vals(1:n_layers, src(m))
 end do
 where (rh > 1.0_kind_real) rh=1.0_kind_real

 n_aerosols=size(var_aerosols)

 do m=1,n_profiles

    call ufo_crtm_calculate_aero_layer_factor(atm(m), layer_factors)

    do i=1,n_aerosols

       varname=var_aerosols(i)
       call ufo_geovals_get_var(geovals,varname, geoval)

       atm(m)%aerosol(i)%Concentration(1:n_layers) = &
            max(geoval%vals(:,src(m))*conf%unit_coef*layer_factors, &
                aerosol_concentration_minvalue_layer)

       if (aerosol_model == "CRTM") then

         !Indices for CRTM default LUT
         !DUST_AEROSOL = 1
         !SEASALT_AEROSOL = 2 - 5
         !ORGANIC_CARBON_AEROSOL = 6
         !BLACK_CARBON_AEROSOL = 7
         !SULFATE_AEROSOL = 8
         select case (varname)

         case (var_sulfate)
            atm(m)%aerosol(i)%type  = 8

         case (var_bcphobic)
            atm(m)%aerosol(i)%type  = 7
         case (var_bcphilic)
            atm(m)%aerosol(i)%type  = 7

         case (var_ocphobic)
            atm(m)%aerosol(i)%type  = 6
         case (var_ocphilic)
            atm(m)%aerosol(i)%type  = 6

         case (var_du001)
            atm(m)%aerosol(i)%type  = 1
            atm(m)%aerosol(i)%effective_radius(:)=dust_radii(1)
         case (var_du002)
            atm(m)%aerosol(i)%type  = 1
            atm(m)%aerosol(i)%effective_radius(:)=dust_radii(2)
         case (var_du003)
            atm(m)%aerosol(i)%type  = 1
            atm(m)%aerosol(i)%effective_radius(:)=dust_radii(3)
         case (var_du004)
            atm(m)%aerosol(i)%type  = 1
            atm(m)%aerosol(i)%effective_radius(:)=dust_radii(4)
         case (var_du005)
            atm(m)%aerosol(i)%type  = 1
            atm(m)%aerosol(i)%effective_radius(:)=dust_radii(5)

         case (var_ss001)
            atm(m)%aerosol(i)%type  = 2
         case (var_ss002)
            atm(m)%aerosol(i)%type  = 3
         case (var_ss003)
            atm(m)%aerosol(i)%type  = 4
         case (var_ss004)
            atm(m)%aerosol(i)%type  = 5
         case (var_ss005)
            atm(m)%aerosol(i)%type  = 5

         case default
            ! not in this LUT: type -1 makes CRTM ignore it
            if (m == 1) call fckit_log%warning(trim(varname) // " is not included in " // &
                                               trim(aerosol_model) // " LUT; ignored")
            atm(m)%aerosol(i)%type  = -1
         end select

         select case (varname)

         case (var_sulfate, var_bcphilic, var_ocphilic, &
               var_ss001, var_ss002, var_ss003, var_ss004)
            do k=1,n_layers

               atm(m)%aerosol(i)%effective_radius(k) = &
                  gocart_aerosol_size(atm(m)%aerosol(i)%type, rh(k,m))
            end do

         case (var_bcphobic, var_ocphobic)
            atm(m)%aerosol(i)%effective_radius(:) = AeroC%Reff(1,atm(m)%aerosol(i)%type)

         case default
            ! other species (dust, or ignored with type -1) keep the default radius

         end select


      else if (aerosol_model == "GOCART-GEOS5") then

         ! This is for the NASA GOCART tables, aerosol scheme GOCART-GEOS5 in CRTM
         ! Reff are bin effective radius from NASA tables
         ! 1, 2, 3, 4, 5  = Dust 1, 2, 3, 4, 5;     with Reff 0.64, 1.32, 2.30, 4.17, 7.67
         ! microns
         ! 6, 7, 8, 9, 10 = Sea salt 1, 2, 3, 4, 5; with Reff 0.08, 0.27, 1.07, 2.55, 7.34
         ! microns
         ! 11, 12 = Organic carbon 1, 2;            with Reff 0.09 microns, hydrophobic and
         ! hydrophilic
         ! 13, 14 = Black carbon 1, 2;              with Reff 0.04 microns, hydrophobic and
         ! hydrophilic
         ! 15, 16 = Sulfate 1, 2;                   with Reff 0.16, 0.60 microns
         ! 17, 18, 19 = Nitrate 1, 2, 3;            with Reff 0.16, 2.10, 7.75 microns
         ! 20, 21 = Brown carbon 1,2;              (Pending release, as of Feb2024)

         ! Cheng: In this code block, the biggest difference from the function
         ! assign_gocart_default
         !        is that no rh adjustmenet to Reff is needed. It'll be interpolated within CRTM.
         !        Aerosol optical properties of NASA tables will not be interpolated over Reff
         !        dimension. There is no need to assign effective radius Reff.
         !        Note that the dust aerosol mapping with the defualt CRTM table
         !        was done outside of function assign_gocart_default, because
         !        dust are all assumed hydrophobic that no rh adjustmenet to Reff
         !        was needed.

         select case (varname)

         ! Dust
         case (var_du001)
            atm(m)%aerosol(i)%type  = 1 ! dust bin 1
         case (var_du002)
            atm(m)%aerosol(i)%type  = 2 ! dust bin 2
         case (var_du003)
            atm(m)%aerosol(i)%type  = 3 ! dust bin 3
         case (var_du004)
            atm(m)%aerosol(i)%type  = 4 ! dust bin 4
         case (var_du005)
            atm(m)%aerosol(i)%type  = 5 ! dust bin 5

         ! Sea salt
          case (var_ss001)
            atm(m)%aerosol(i)%type  = 6 ! sea salt bin 1
         case (var_ss002)
            atm(m)%aerosol(i)%type  = 7 ! sea salt bin 2
         case (var_ss003)
            atm(m)%aerosol(i)%type  = 8 ! sea salt bin 3
         case (var_ss004)
            atm(m)%aerosol(i)%type  = 9 ! sea salt bin 4
         case (var_ss005)
            atm(m)%aerosol(i)%type  = 10 ! sea salt bin 5

         ! Organic carbon
         case (var_ocphobic)
            atm(m)%aerosol(i)%type  = 11 ! hydrophobic
         case (var_ocphilic)
            atm(m)%aerosol(i)%type  = 12 ! hydrophilic

         ! Black carbon
         case (var_bcphobic)
            atm(m)%aerosol(i)%type  = 13 ! hydrophobic
         case (var_bcphilic)
            atm(m)%aerosol(i)%type  = 14 ! hydrophilic

         ! Sulfate
         case (var_sulfate)
            atm(m)%aerosol(i)%type  = 15 ! bin 1
         !Jerome: There is only one bin of sulfate so far
         !CASE (var_sulfate)
         !   atm(m)%aerosol(i)%TYPE  = 16 ! bin 2

         ! Nitrate
         case (var_no3an1)
            atm(m)%aerosol(i)%type  = 17 ! bin 1
         case (var_no3an2)
            atm(m)%aerosol(i)%type  = 18 ! bin 2
         case (var_no3an3)
            atm(m)%aerosol(i)%type  = 19 ! bin 3

         ! Brown carbon (pending release, comment out now)
         case (var_brphobic)
            atm(m)%aerosol(i)%type  = 20 ! bin 1
         case (var_brphilic)
            atm(m)%aerosol(i)%type  = 21 ! bin 2

         case default
            ! not in this LUT: type -1 makes CRTM ignore it
            if (m == 1) call fckit_log%warning(trim(varname) // " is not included in " // &
                                               trim(aerosol_model) // " LUT; ignored")
            atm(m)%aerosol(i)%type  = -1
         end select

       else if (aerosol_model == "CMAQ") then
         ! Aerosol scheme CMAQ in CRTM
         ! CMAQ table:
         ! Dust - 1
         ! Soot - 2
         ! Water soluble - 3 (RI of OPAC WASO are used as RI of OC in GOCART)
         ! Sulfate - 4
         ! Sea salt - 5
         ! Water - 6
         ! Insoluble -7
         ! dust-like - 8

         ! Place holder for effective radius variance, set as 1.0
         do k=1,n_layers
            atm(m)%aerosol(i)%effective_variance(k)=1.0_kind_real
         end do

         ! Assign aerosol type
         select case (varname)

           !Dust
           case (var_du001)
              atm(m)%aerosol(i)%type = 1
           case (var_du002)
              atm(m)%aerosol(i)%type = 1
           case (var_du003)
              atm(m)%aerosol(i)%type = 1
           case (var_du004)
              atm(m)%aerosol(i)%type = 1
           case (var_du005)
              atm(m)%aerosol(i)%type = 1

           !Sea Salt
           case (var_ss001)
              atm(m)%aerosol(i)%type = 5
           case (var_ss002)
              atm(m)%aerosol(i)%type = 5
           case (var_ss003)
              atm(m)%aerosol(i)%type = 5
           case (var_ss004)
              atm(m)%aerosol(i)%type = 5
           case (var_ss005)
              atm(m)%aerosol(i)%type = 5

           ! Organic carbon
           case (var_ocphobic)
              atm(m)%aerosol(i)%type = 7
           case (var_ocphilic)
              atm(m)%aerosol(i)%type = 3

           ! Black carbon
           case (var_bcphobic)
              atm(m)%aerosol(i)%type = 7
           case (var_bcphilic)
              atm(m)%aerosol(i)%type = 2

           ! Sulfate
           case (var_sulfate)
              atm(m)%aerosol(i)%type = 4

           ! Nitrate
           case (var_no3an1)
              atm(m)%aerosol(i)%type = 3
           case (var_no3an2)
              atm(m)%aerosol(i)%type = 3
           case (var_no3an3)
              atm(m)%aerosol(i)%type = 3

           !!Brown carbon
           case (var_brphobic)
              atm(m)%aerosol(i)%type = 7
           case (var_brphilic)
              atm(m)%aerosol(i)%type = 3

         case default
            err_msg = trim(varname) // " is not included in " // trim(aerosol_model) // &
                      " LUT: " // trim(conf%aerosol_option) // " not ready in UFO/AODCRTM"
            call fckit_exception % throw(err_msg)

         end select

         select case (varname)

         ! Reff for hydrophilic aerosols
         case (var_sulfate, var_bcphilic, var_ocphilic, &
               var_ss001, var_ss002, var_ss003, var_ss004, &
               var_no3an1, var_no3an2, var_no3an3)
            do k=1,n_layers
               atm(m)%aerosol(i)%effective_radius(k) = &
                  gocart_aerosol_size(atm(m)%aerosol(i)%type, rh(k,m))
            end do
         ! Reff for hydrophobic aerosols
         case (var_du001, var_du002, var_du003, var_du004, var_du005, &
               var_bcphobic, var_ocphobic)
            atm(m)%aerosol(i)%effective_radius(:) = AeroC%Reff(1,atm(m)%aerosol(i)%type)

         case default
            ! other species (dust, or ignored with type -1) keep the default radius

         end select

       else if (aerosol_model == "NAAPS") then
         ! This is for the NRL NAAPS tables, aerosol scheme NAAPS in CRTM
         ! Similar to GOCART-GEOS5 LUT, no Reff needs to be assigned
         ! NAAPS table is recommended for AOD calculation only due to
         ! over-simplified phase function reconstruction.
         ! NAAPS table:
         ! Dust - 1
         ! Smoke - 2
         ! Sea Salt - 3
         ! Anthropogenic and Biogenic Fine Particles - 4

         ! Assign aerosol type
         select case (varname)

         !Dust
         case (var_du001)
            atm(m)%aerosol(i)%type = 1
         case (var_du002)
            atm(m)%aerosol(i)%type = 1
         case (var_du003)
            atm(m)%aerosol(i)%type = 1
         case (var_du004)
            atm(m)%aerosol(i)%type = 1
         case (var_du005)
            atm(m)%aerosol(i)%type = 1

         !Sea Salt
         case (var_ss001)
           atm(m)%aerosol(i)%type  = 3
         case (var_ss002)
           atm(m)%aerosol(i)%type  = 3
         case (var_ss003)
           atm(m)%aerosol(i)%type  = 3
         case (var_ss004)
           atm(m)%aerosol(i)%type  = 3
         case (var_ss005)
           atm(m)%aerosol(i)%type  = 3

         ! Organic carbon
         case (var_ocphobic)
            atm(m)%aerosol(i)%type = 2
         case (var_ocphilic)
            atm(m)%aerosol(i)%type = 2

         ! Black carbon
         case (var_bcphobic)
            atm(m)%aerosol(i)%type = 2
         case (var_bcphilic)
            atm(m)%aerosol(i)%type = 2

         ! Sulfate
         case (var_sulfate)
            atm(m)%aerosol(i)%type = 4

         ! Nitrate
         case (var_no3an1)
            atm(m)%aerosol(i)%type = 4
         case (var_no3an2)
            atm(m)%aerosol(i)%type = 4
         case (var_no3an3)
            atm(m)%aerosol(i)%type = 4

         !Brown carbon
         case (var_brphobic)
            atm(m)%aerosol(i)%type = 2
         case (var_brphilic)
            atm(m)%aerosol(i)%type = 2

         case default
            err_msg = trim(varname) // " is not included in " // trim(aerosol_model) // &
                      " LUT: " // trim(conf%aerosol_option) // " not ready in UFO/AODCRTM"
            call fckit_exception % throw(err_msg)

         end select

       end if
    end do
  end do

  deallocate(src)
end subroutine ufo_crtm_load_aerosol_data

! ------------------------------------------------------------------------------

subroutine calculate_aero_layer_factor_atm_profile(atm, layer_factors)

  type(CRTM_atmosphere_type), intent(in) :: atm
  real(kind_real), intent(out) :: layer_factors(:)

  integer :: k

  do k=1,size(layer_factors)
     ! Convert layer_factors to a dry-air mass basis. delta_p/g is the moist-air
     ! mass, so use the mass fraction 1/(1+mixing_ratio), not the dry partial-pressure
     ! ratio of Cotton eq. (2.4), p_dry=p_total/(1+r_v/r_d*mixing_ratio) (ufo#4450).
     layer_factors(k)=1e-9_kind_real*(atm%Level_Pressure(k)- &
          atm%Level_Pressure(k-1))*100.0_kind_real/grav/ &
          (1.0_kind_real+atm%Absorber(k,1)*1e-3_kind_real)
  end do

end subroutine calculate_aero_layer_factor_atm_profile

! ------------------------------------------------------------------------------

subroutine calculate_aero_layer_factor_atm(atm, layer_factors)

  type(CRTM_atmosphere_type), intent(in) :: atm(:)
  real(kind_real), intent(out) :: layer_factors(:,:)
  integer :: k,m

  do k=1,size(layer_factors,1)
     do m=1,size(layer_factors,2)
        ! Convert layer_factors to a dry-air mass basis. delta_p/g is the moist-air
        ! mass, so use the mass fraction 1/(1+mixing_ratio), not the dry partial-pressure
        ! ratio of Cotton eq. (2.4), p_dry=p_total/(1+r_v/r_d*mixing_ratio) (ufo#4450).
        layer_factors(k,m)=1e-9_kind_real*(atm(m)%Level_Pressure(k)- &
             atm(m)%Level_Pressure(k-1))*100.0_kind_real/grav/ &
             (1.0_kind_real+atm(m)%Absorber(k,1)*1.e-3_kind_real)
     end do
  end do

end subroutine calculate_aero_layer_factor_atm

! ------------------------------------------------------------------------------

function gocart_aerosol_size( itype, rh ) & ! rh input in 0-1
     RESULT(r_eff)   ! in micrometer

  use CRTM_aerosolcoeff, only: aeroc
  implicit none

!
!   modified from a function provided by quanhua liu
!
  integer ,intent(in) :: itype
  real(kind_real)    ,intent(in) :: rh

  integer :: j1,j2,m
  integer :: n_rh
  real(kind_real)    :: h1
  real(kind_real)    :: r_eff

  j2 = 0
  j1 = 1
  n_rh = 35    ! total number of RH in CRTM Default LUTs

  if ( rh <= aeroc%rh(1) ) then
     j1 = 1
  else if ( rh >= aeroc%rh(n_rh) ) then
     j1 = n_rh
  else
     do m = 1, n_rh-1
        if ( rh < aeroc%rh(m+1) .and. rh > aeroc%rh(m) ) then
           j1 = m
           j2 = m+1
           h1 = (rh-aeroc%rh(m))/(aeroc%rh(m+1)-aeroc%rh(m))
           exit
        end if
     end do
  end if

  if ( j2 == 0 ) then
     r_eff = aeroc%reff(j1,itype )
  else
     r_eff = (1.0_kind_real-h1)*aeroc%reff(j1,itype ) + h1*aeroc%reff(j2,itype )
  end if

end function gocart_aerosol_size

! ------------------------------------------------------------------------------

subroutine calc_thompson_reff(rho_air,tsen,qxmr,cloud_type_Id,ni,nr,n_Layers,reff)
 !$$$  subprogram documentation block
 !                .      .    .                                       .
 ! subprogram:    calc_thompson_reff:  calculate effective radius based on Thompson scheme
 !   prgmmr:      Azadeh.Gholoubi
 !
 ! program history log:
 !
 !   2023-01-24
 !
 !   language: f90
 !
 !$$$
 !--------
 use ufo_constants_mod, only: zero, pi, t0c, rhor, rhow, am_s
 implicit none

 ! Declare passed variables
 integer                              ,intent(in   ) :: cloud_type_Id
 integer                              ,intent(in   ) :: n_Layers
 real(kind_real), dimension(n_Layers) ,intent(in   ) :: rho_air   ! [ kg/m3  ]
 real(kind_real), dimension(n_Layers) ,intent(in   ) :: tsen      ! [ K      ]
 real(kind_real), dimension(n_Layers) ,intent(in   ) :: qxmr      ! [ kg/kg  ]
 real(kind_real), dimension(n_Layers) ,intent(in   ) :: ni        ! cloud_ice_number_concentration
 real(kind_real), dimension(n_Layers) ,intent(in   ) :: nr        ! rain_number_concentration
 real(kind_real), dimension(n_Layers) ,intent(inout) :: reff      ! [ micron ]

 ! Declare local variables
 integer(c_int) :: k
 integer(c_int) :: mu_w
 real(kind_real)    :: qx
 real(kind_real)    :: reff_min, reff_max
 real(kind_real):: lam_i,lam_w,lam_r, lam_g,lam_exp, am_r,am_w,am_i,am_g

 ! Parameters
 real(kind_real), parameter :: qmin = 1.0e-12_kind_real          ! [kg/kg ]
 !.. droplet number concentration.
 real(kind_real), parameter :: Nt_c = 100.E6_kind_real  ![m-3]

 ! Parameters for water cloud
 real(kind_real), parameter :: reff_w_min =    2.0_kind_real     ! previous value was 5.0_kind_real
 ! previous value was 10.0_kind_real
 real(kind_real), parameter :: reff_w_max =   25.0_kind_real

 ! Parameters for ice cloud (Heymsfield and McFarquhar 1996)
 real(kind_real), parameter :: rho_i      = 890.0_kind_real      ! [kg/m3 ]
 real(kind_real), parameter :: reff_i_min =    2.5_kind_real     ! previous value was 10_kind_real
 real(kind_real), parameter :: reff_i_max = 250.0_kind_real      ! previous value was 150_kind_real
 real(kind_real), parameter:: mu_i = 0.0_kind_real
 ! minimum number concentration (ccpp-physics)
 real(kind_real), parameter:: ni_min = 1.0e-6_kind_real

 ! Parameters for  rain (Lin 1983)
 ! [micron] ! previous value was 0.0_kind_real
 real(kind_real), parameter :: reff_r_min =       50.0_kind_real
 ! [micron] ! previous value was 10000.0_kind_real
 real(kind_real), parameter :: reff_r_max =   1000.0_kind_real
 real(kind_real), parameter:: mu_r = 0.0_kind_real
 ! Parameters for snow
 ! [micron] ! previous value was 0.0_kind_real
 real(kind_real), parameter :: reff_s_min =       5.0_kind_real
 ! [micron] ! previous value was 10000.0_kind_real
 real(kind_real), parameter :: reff_s_max =   5000.0_kind_real
 ! minimum number concentration (ccpp-physics)
 real(kind_real), parameter:: nr_min = 1.0e-6_kind_real

 !For snow moments conversions  (from Field et al. 2005)
 real(kind_real), dimension(10), parameter:: &
     sa = [ 5.065339_kind_real, -0.062659_kind_real, -3.032362_kind_real,  0.029469_kind_real, &
           -0.000285_kind_real,  0.31255_kind_real,   0.000204_kind_real,  0.003199_kind_real, &
            0.0_kind_real,      -0.015952_kind_real ]
     real(kind_real), dimension(10), parameter:: &
     sb = [ 0.476221_kind_real, -0.015896_kind_real,  0.165977_kind_real,  0.007468_kind_real, &
           -0.000141_kind_real,  0.060366_kind_real,  0.000079_kind_real,  0.000594_kind_real, &
            0.0_kind_real,      -0.003577_kind_real ]
 real(kind_real), parameter :: bm_s = 2.0_kind_real
 real(kind_real), dimension(1), parameter :: cse = [ bm_s + 1.0_kind_real ]
 real(kind_real) :: tc0, smob, smoc, a_, b_, loga_

 ! Parameters for graupel  (Lin 1983)
 ! [kg/m3 ] ! differs from rhog (400) in ufo_constants_mod
 real(kind_real), parameter :: rho_g      =    500.0_kind_real
 ! [micron] ! previous value was 0.0_kind_real
 real(kind_real), parameter :: reff_g_min =    150.0_kind_real
 ! [micron] ! previous value was 10000.0_kind_real
 real(kind_real), parameter :: reff_g_max =   5000.0_kind_real
 real(kind_real), parameter:: mu_g = 0.0_kind_real
 real(kind_real), parameter :: no_exp     =  10.0_kind_real

 !cloud water
 if (cloud_type_Id==WATER_CLOUD) then
    am_w = rhow*pi/6.0_kind_real
    reff_min = reff_w_min
    reff_max = reff_w_max
    do k = n_Layers, 1, -1
       qx = qxmr(k) * rho_air(k)  ! convert mixing ratio (kg/kg) to water content (kg/m3)
       if (qx > qmin) then
          mu_w = max(2, min((nint(1000.E6_kind_real/Nt_c) + 2), 15))
          lam_w=exp(1.0_kind_real / 3.0_kind_real * log ((am_w*Nt_c *gamma(mu_w + 3.0_kind_real &
                  + 1.0_kind_real))/(qx*gamma(mu_w+1.0_kind_real))))

          reff(k) = 0.5_kind_real * ((3.0_kind_real+mu_w)/lam_w)*1.0e6_kind_real
          reff(k) = max(reff_min, min(reff_max, reff(k)))
       else
          reff(k) = zero
       end if
    end do

 ! Cloud Ice
 else if (cloud_type_Id==ICE_CLOUD) then
    am_i = rho_i*pi/6.0_kind_real
    reff_min = reff_i_min
    reff_max = reff_i_max
    do k = n_Layers, 1, -1
       qx = qxmr(k) * rho_air(k)  ! convert mixing ratio (kg/kg) to water content (kg/m3)
       if (qx > qmin .and. ni(k)>ni_min) then
          lam_i=exp(1.0_kind_real / 3.0_kind_real * log((am_i*ni(k) *gamma(mu_i + 3.0_kind_real &
                  + 1.0_kind_real))/(qx*gamma(mu_i+1.0_kind_real))))
          reff(k) = 0.5_kind_real * (3.0_kind_real /lam_i)*1.0e6_kind_real
          reff(k) = max(reff_min, min(reff_max, reff(k)))
       else
          reff(k) = zero
       end if
    end do
 !Rain
 else if (cloud_type_Id==RAIN_CLOUD) then
    am_r = rhor*pi/6.0_kind_real
    reff_min = reff_r_min
    reff_max = reff_r_max
    do k = n_Layers, 1, -1
       qx = qxmr(k) * rho_air(k)  ! convert mixing ratio (kg/kg) to water content (kg/m3)
       if (qx > qmin .and. nr(k)>nr_min) then
          lam_r=exp(1.0_kind_real / 3.0_kind_real * log ((am_r*nr(k) *gamma(mu_r + 3.0_kind_real &
                  + 1.0_kind_real))/(qx*gamma(mu_r + 1.0_kind_real))))
          reff(k) = 0.5_kind_real *(3.0_kind_real/lam_r)*1.0e6_kind_real
          reff(k) = max(reff_min, min(reff_max, reff(k)))
       else
          reff(k) = zero
       end if
    end do

 ! Snow (Field et al. 2005)

 else if (cloud_type_Id==SNOW_CLOUD) then
    reff_min = reff_s_min
    reff_max = reff_s_max
    do k = n_Layers, 1, -1
 !..Calculate bm_s+1 (th) moment, smoc.  Useful for diameter calcs.
       tc0 = min(-0.1_kind_real, tsen(k)-t0c)
       qx = qxmr(k) * rho_air(k)  ! convert mixing ratio (kg/kg) to water content (kg/m3)
       if (qx > qmin) then
          smob =qx/am_s
          loga_ = sa(1) + sa(2)*tc0 + sa(3)*cse(1) &
                + sa(4)*tc0*cse(1) + sa(5)*tc0*tc0 &
                + sa(6)*cse(1)*cse(1) + sa(7)*tc0*tc0*cse(1) &
                + sa(8)*tc0*cse(1)*cse(1) + sa(9)*tc0*tc0*tc0 &
                + sa(10)*cse(1)*cse(1)*cse(1)
          a_ = 10.0_kind_real**loga_
          b_ = sb(1)+ sb(2)*tc0 + sb(3)*cse(1) + sb(4)*tc0*cse(1) &
             + sb(5)*tc0*tc0 + sb(6)*cse(1)*cse(1) &
             + sb(7)*tc0*tc0*cse(1) + sb(8)*tc0*cse(1)*cse(1) &
             + sb(9)*tc0*tc0*tc0 + sb(10)*cse(1)*cse(1)*cse(1)
          smoc = a_ * smob**b_
          reff(k) = max(2.51E-6_kind_real, min(0.5_kind_real*(smoc/smob), 1999.E-6_kind_real))
          reff(k) = max(reff_min, min(reff_max, reff(k)*1.0e6_kind_real))
       else
          reff(k) = zero
       end if
    end do

 ! Graupel
 else if (cloud_type_Id==GRAUPEL_CLOUD) then
    reff_min = reff_g_min
    reff_max = reff_g_max
    am_g = rho_g*pi/6.0_kind_real
    do k = n_Layers, 1, -1
       qx = qxmr(k)*rho_air(k)    ! convert mixing ratio (kg/kg) to water content (kg/m3)
       if (qx > qmin) then
          lam_exp = no_exp* exp(1.0_kind_real/(3.0_kind_real+1.0_kind_real) * &
                  log ((am_g*gamma(3.0_kind_real +1.0_kind_real))/qx))
          lam_g=lam_exp*exp(1/3.0_kind_real*log(gamma(3.0_kind_real+mu_g+1.0_kind_real) &
                  /((3.0_kind_real+mu_g+1.0_kind_real)*(mu_g+1.0_kind_real))))
          reff(k) = 0.5_kind_real *(3.0_kind_real/ lam_g)*1.0e6_kind_real
          reff(k) = max(reff_min, min(reff_max, reff(k)))
       else
          reff(k) = zero
       end if
    end do

 end if

end subroutine calc_thompson_reff

! ------------------------------------------------------------------------------

subroutine calc_thompson_cloudfrac(n_Layers, lmfdeep2, prsl, clwf, rhly, qstl,cldtot)
 !$$$  subprogram documentation block
 !                .      .    .                                       .
 ! subprogram:    calc_thompson_cloudfrac  calculate cloud fractions from cloud using Thompson scheme
 !
 !   prgmmr:      Azadeh.Gholoubi
 !
 ! abstract:
 !
 ! program history log:
 !
 !   2023-01-24
 !
 !   language: f90
 !
 !  ====================  definition of variables  ====================  !
 ! input variables:                                                      !
 !   prsl  (msig) : model layer mean pressure in cb (1000 Pa)            !
 !   qstl  (n_Layers) : layer saturate humidity in gm/gm                 !
 !   rhly  (n_Layers) : layer relative humidity (=qlyr/qstl)             !
 !   n_Layers         : vertical layer/level dimensions                  !
 !   clwf         : sum of the all 5 Hydrometeors                        !
 !   lmfdeep2     : scale-aware mass-flux deep conv scheme flag          !
 ! output variables:                                                     !
 !   cldtot (n_Layers) - layer total cloud fraction                      !
 !  ====================    end of description    =====================  !
 !$$$
 !--------
 use ufo_constants_mod
 implicit none

 !  ---  inputs:
 integer(c_int), intent(in) :: n_Layers
 real(kind_real) , dimension(n_Layers), intent(in) :: prsl, clwf, rhly, qstl
 logical, intent(in) :: lmfdeep2

 !  ---  outputs
 real(kind_real) , dimension(n_Layers), intent(inout) :: cldtot

 !  ---  local variables:
 real(kind_real)  :: clwmin, clwm, clwt, onemrh, value, tem1, tem2
 integer(c_int) :: k

 !> - Compute layer cloud fraction.
 clwmin = 0.0_kind_real
 do k = n_Layers, 2, -1  ! Layer numbers are top-down in UFO.
 ! converting prsl from cb (1000Pa) to bar (100,000 Pa)
    clwt = 1.0e-10_kind_real * (prsl(k)*0.01_kind_real)
    if (clwf(k) > clwt) then
       if(rhly(k) > 0.99_kind_real) then
          cldtot(k) = 1.0_kind_real
       else
          onemrh= max( 1.e-10_kind_real, 1.0_kind_real-rhly(k) )
          clwm  = clwmin / max( 0.01_kind_real, prsl(k)*0.01_kind_real )
          tem1  = min(max((onemrh*qstl(k))**0.49_kind_real,0.0001_kind_real),1.0_kind_real)
          if (lmfdeep2) then
             tem1  = xrc3 / tem1
          else
             tem1  = 100.0_kind_real / tem1
          end if
          value = max( min( tem1*(clwf(k)-clwm), 50.0_kind_real ), 0.0_kind_real )
          tem2  = sqrt( sqrt(rhly(k)) )
          cldtot(k) = max( tem2*(1.0_kind_real-exp(-value)), 0.0_kind_real )
       end if
    else
       cldtot(k) = 0.0_kind_real
    end if
 end do
 where(cldtot < 0.0_kind_real) cldtot = 0.0_kind_real
 where(cldtot > 1.0_kind_real) cldtot = 1.0_kind_real
end subroutine calc_thompson_cloudfrac

! ------------------------------------------------------------------------------

end module ufo_crtm_utils_mod
