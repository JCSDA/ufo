! (C) Copyright 2018 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

!> Fortran module to provide code shared between nonlinear and tlm/adm radiance calculations

MODULE ufo_crtm_utils_mod

use fckit_configuration_module, only: fckit_configuration
use fckit_mpi_module,   only: fckit_mpi_comm
use iso_c_binding
use kinds

use crtm_module

use ufo_geovals_mod, only: ufo_geovals, ufo_geoval, ufo_geovals_get_var
use ufo_vars_mod
use obsspace_mod
use ufo_utils_mod, only: cmp_strings
use CRTM_SpcCoeff, only: CRTM_SpcCoeff_Load, SC

implicit none
private

public crtm_comm_stat_check
public crtm_conf
public crtm_conf_setup
public crtm_conf_delete
public get_var_name
public Load_Atm_Data
public Load_Sfc_Data
public Load_Geom_Data
public ufo_crtm_skip_profiles

PUBLIC Load_Aerosol_Data
public assign_aerosol_names
public define_aerosol_model
public calculate_aero_layer_factor
PUBLIC :: qsmith
PUBLIC :: upper2lower

PUBLIC :: grav,rv_rd, threshold_reflectivity, &
     &aerosol_concentration_minvalue,aerosol_concentration_minvalue_layer

REAL(kind_real), PARAMETER :: &
     &rdgas = 2.8704e+2_kind_real,&
     &rvgas = 4.6150e+2_kind_real,&
     &rv_rd = rvgas/rdgas,&
     &esl = 0.621971831,&
     &zvir =  rv_rd - 1_kind_real,&
     &tice = 273.16_kind_real,&
     &grav = 9.81_kind_real,&
     &aerosol_concentration_minvalue=1.e-16_kind_real,&
     &aerosol_concentration_minvalue_layer=tiny(rdgas),&
     &ozone_default_value=1.e-3_kind_real ! in ppmv in crtm

integer, parameter, public :: max_string=800
real(kind_real), parameter  :: threshold_reflectivity = 100.0

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

 character(len=255), allocatable :: SENSOR_ID(:)
 character(len=255) :: ENDIAN_TYPE
 character(len=255) :: COEFFICIENT_PATH, NC_COEFFICIENT_PATH
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
 real(kind_real) :: unit_coef = 1.0_kind_real
 logical :: Cloud_Seeding = .false.
end type crtm_conf


INTERFACE calculate_aero_layer_factor

   MODULE PROCEDURE calculate_aero_layer_factor_atm_profile,&
        &calculate_aero_layer_factor_atm

END INTERFACE

INTERFACE qsmith

   MODULE PROCEDURE qsmith_atm,qsmith_profiles

END INTERFACE qsmith


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


! Surface Variables

 character(len=MAXVARLEN), parameter :: &
      UFO_Surfaces(7) = &
         [ var_sfc_wtmp, var_sfc_ltmp, var_sfc_itmp, var_sfc_stmp, var_sfc_wspeed, var_sfc_wdir, &
           var_sfc_sss ]

 character(len=MAXVARLEN), parameter :: &
      CRTM_Surfaces(7) = &
         [  character(len=MAXVARLEN):: 'Water_Temperature', 'Land_Temperature', 'Ice_Temperature', &
                                       'Snow_Temperature', 'Wind_Speed', &
                                       'Wind_Direction', 'Salinity' ]

 character(len=MAXVARLEN), parameter :: &
      ValidSurfaceWindGeoVars(2) = [character(len=MAXVARLEN) :: 'vector', 'uv']



contains

! ------------------------------------------------------------------------------

subroutine crtm_conf_setup(conf, f_confOpts, f_confOper, comm)

implicit none
type(crtm_conf),                intent(inout) :: conf
type(fckit_configuration),      intent(in)    :: f_confOpts
type(fckit_configuration),      intent(in)    :: f_confOper
type(fckit_mpi_comm), optional, intent(in)    :: comm

character(*), PARAMETER :: routine_name = 'crtm_conf_setup'
character(len=255) :: IRwaterCoeff, VISwaterCoeff, &
                      IRVISlandCoeff, IRVISsnowCoeff, IRVISiceCoeff, &
                      MWwaterCoeff
integer :: jspec, ivar
character(len=max_string) :: message
character(len=:), allocatable :: str
character(len=:), allocatable :: str_array(:)

CHARACTER(len=MAXVARLEN), ALLOCATABLE :: var_aerosols(:)
logical :: message_flag = .true.

 !Some config needs to come from user
 !-----------------------------------

 !Number of sensors, each call to CRTM will be for a single sensor
 !type (zenith/scan angle will be different)
 conf%n_Sensors = 1

 ! absorbers, clouds and aerosols (should match what model will provide)

 ! Set print rank
 ! --------------
 if (present(comm)) then
   if (comm%rank() .ne. 0) message_flag = .false.
 endif

 ! Absorbers
 !----------
 conf%n_Absorbers = 0
 if (f_confOper%has("Absorbers")) &
   conf%n_Absorbers = conf%n_Absorbers + f_confOper%get_size("Absorbers")

 allocate( conf%Absorbers     ( conf%n_Absorbers ), &
           conf%Absorber_Id   ( conf%n_Absorbers ), &
           conf%Absorber_Units( conf%n_Absorbers ) )

 if (conf%n_Absorbers > 0) then
   call f_confOper%get_or_die("Absorbers",str_array)
   conf%Absorbers(1:conf%n_Absorbers) = str_array
 end if

 ! check for duplications
 do jspec = 2, conf%n_Absorbers
   if ( any(conf%Absorbers(jspec-1) == conf%Absorbers(jspec:conf%n_Absorbers)) ) then
     write(message,*) trim(ROUTINE_NAME),' error: ',trim(conf%Absorbers(jspec)),' is duplicated in Absorbers'
     call abor1_ftn(message)
   end if
 end do

 ! convert from CRTM names to UFO CF names and define Id and Units
 do jspec = 1, conf%n_Absorbers
   ivar = ufo_vars_getindex(CRTM_Absorbers, conf%Absorbers(jspec))
   if (ivar < 1 .or. ivar > size(UFO_Absorbers)) then
     write(message,*) trim(ROUTINE_NAME),' error: ',trim(conf%Absorbers(jspec)),' not supported by UFO_Absorbers'
     call abor1_ftn(message)
   end if
   conf%Absorbers(jspec) = UFO_Absorbers(ivar)
   conf%Absorber_Id(jspec) = CRTM_Absorber_Id(ivar)
   conf%Absorber_Units(jspec) = CRTM_Absorber_Units(ivar)
 end do


 ! Clouds
 !-------
 conf%n_Clouds = 0
 if (f_confOper%has("Clouds")) &
   conf%n_Clouds = f_confOper%get_size("Clouds")
 allocate( conf%Clouds  ( conf%n_Clouds,2), &
           conf%Cloud_Id( conf%n_Clouds ) )
 if (conf%n_Clouds > 0) then
   call f_confOper%get_or_die("Clouds",str_array)
   conf%Clouds(1:conf%n_Clouds,1) = str_array

   if (f_confOper%has("Cloud_Fraction")) then
     call f_confOper%get_or_die("Cloud_Fraction",conf%Cloud_Fraction)
     if ( conf%Cloud_Fraction < 0.0 .or. &
          conf%Cloud_Fraction > 1.0 ) then
       write(message,*) trim(ROUTINE_NAME),' error: must specify ' // &
                        ' 0.0 <= Cloud_Fraction <= 1.0' // &
                        ' or remove Cloud_Fraction from conf' // &
                        ' and provide as a geoval'
       call abor1_ftn(message)
     end if
   else
     message = trim(ROUTINE_NAME) // &
             ': Cloud_Fraction is not provided in conf.' // &
             ' Will request as a geoval.'
     if (message_flag) CALL Display_Message(ROUTINE_NAME, TRIM(message), WARNING )
   end if
   if (f_confOper%has("Cloud_Seeding")) then
     call f_confOper%get_or_die("Cloud_Seeding",conf%Cloud_Seeding)
     if ( conf%Cloud_Seeding ) then
       write(message,*) trim(ROUTINE_NAME),' Cloud Seeding is activated '
     else
       write(message,*) trim(ROUTINE_NAME),' Cloud Seeding is not activated '
     endif
   end if
 end if

 ! check for duplications
 do jspec = 2, conf%n_Clouds
   if ( any(conf%Clouds(jspec-1,1) == conf%Clouds(jspec:conf%n_Clouds,1)) ) then
     write(message,*) trim(ROUTINE_NAME),' error: ',trim(conf%Clouds(jspec,1)), &
                      ' is duplicated in Clouds'
     call abor1_ftn(message)
   end if
 end do

 ! convert from CRTM names to UFO CF names and define Id
 do jspec = 1, conf%n_Clouds
   ivar = ufo_vars_getindex(CRTM_Clouds, conf%Clouds(jspec,1))
   if (ivar < 1 .or. ivar > size(UFO_Clouds)) then
     write(message,*) trim(ROUTINE_NAME),' error: ',trim(conf%Clouds(jspec,1)),' not supported by UFO_Clouds'
     call abor1_ftn(message)
   end if

   conf%Clouds(jspec,1:2) = UFO_Clouds(ivar,1:2)
   conf%Cloud_Id(jspec)   = CRTM_Cloud_Id(ivar)
 end do

 ! Aerosols
 !---------
 IF (f_confOpts%has("AerosolOption")) THEN
    call f_confOpts%get_or_die("AerosolOption",str)
    conf%aerosol_option = str
    conf%aerosol_option = upper2lower(conf%aerosol_option)
    CALL assign_aerosol_names(conf%aerosol_option,var_aerosols)
    conf%n_Aerosols=SIZE(var_aerosols)
 ELSE
    conf%n_Aerosols  = 0
    conf%aerosol_option = ""
 ENDIF
 call f_confOpts%get_or_die("model units coeff", conf%unit_coef)

 ! Surface variables
 !----------
 conf%n_Surfaces = 0
 if (f_confOper%has("Surfaces")) &
   conf%n_Surfaces = conf%n_Surfaces + f_confOper%get_size("Surfaces")

 allocate( conf%Surfaces    ( conf%n_Surfaces ))

 if (conf%n_Surfaces > 0) then
   call f_confOper%get_or_die("Surfaces",str_array)
   conf%Surfaces(1:conf%n_Surfaces) = str_array
 end if

 ! check for duplications
 do jspec = 2, conf%n_Surfaces
   if ( any(conf%Surfaces(jspec-1) == conf%Surfaces(jspec:conf%n_Surfaces)) ) then
     write(message,*) 'crtm_conf_setup error: ',trim(conf%Surfaces(jspec)),' is duplicated in Surfaces'
     call abor1_ftn(message)
   end if
 end do

 ! convert from CRTM names to UFO CF names and define Id and Units
 do jspec = 1, conf%n_Surfaces
   ivar = ufo_vars_getindex(CRTM_Surfaces, conf%Surfaces(jspec))
   if (ivar < 1 .or. ivar > size(UFO_Surfaces)) then
     write(message,*) 'crtm_conf_setup error: ',trim(conf%Surfaces(jspec)),' not supported by UFO_Surfaces'
     call abor1_ftn(message)
   end if
   conf%Surfaces(jspec) = UFO_Surfaces(ivar)

 end do

 ! select between two surface wind geovals options
 ! valid options: vector [default], uv
 if (f_confOper%get('SurfaceWindGeoVars', str)) then
   conf%sfc_wind_geovars = str
 else
   conf%sfc_wind_geovars = 'vector'
 endif
 if (ufo_vars_getindex(ValidSurfaceWindGeoVars, conf%sfc_wind_geovars) < 1) then
   write(message,*) 'crtm_conf_setup error: invalid SurfaceWindGeoVars ',trim(str)
   call abor1_ftn(message)
 end if

 ! Sea_Surface_Salinity
 !---------
 IF (f_confOpts%get("Salinity",str)) THEN
    conf%salinity_option = str
 ELSE
    conf%salinity_option = 'off'
 END IF

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
 if (f_confOpts%has("NC_CoefficientPath")) then
    call f_confOpts%get_or_die("NC_CoefficientPath",str)
    conf%NC_COEFFICIENT_PATH = str
 endif

 ! Cloud coefficient file, model, and format
 conf%Cloud_Model = "CRTM"
 if (f_confOpts%has("Cloud_Model")) then
    call f_confOpts%get_or_die("Cloud_Model",str)
    conf%Cloud_Model = str
 end if

 conf%CloudCoeff_File = "CloudCoeff.bin"
 if (f_confOpts%has("CloudCoeff_File")) then
    call f_confOpts%get_or_die("CloudCoeff_File",str)
    conf%CloudCoeff_File = str
 end if

 conf%CloudCoeff_Format = "Binary"
 if (f_confOpts%has("CloudCoeff_Format")) then
    call f_confOpts%get_or_die("CloudCoeff_Format",str)
    conf%CloudCoeff_Format = str
 end if

 ! Aerosol coefficient file, format, and format
 conf%Aerosol_Model = 'CRTM'
 if (f_confOpts%has("Aerosol_Model")) then
    call f_confOpts%get_or_die("Aerosol_Model",str)
    conf%Aerosol_Model = str
 end if

 conf%AerosolCoeff_File = "AerosolCoeff.bin"
 if (f_confOpts%has("AerosolCoeff_File")) then
    call f_confOpts%get_or_die("AerosolCoeff_File",str)
    conf%AerosolCoeff_File = str
 end if

 conf%AerosolCoeff_Format = "Binary"
 if (f_confOpts%has("AerosolCoeff_Format")) then
    call f_confOpts%get_or_die("AerosolCoeff_Format",str)
    conf%AerosolCoeff_Format = str
 end if

 ! Coefficient file prefixes
 IRwaterCoeff = "Nalli"
 if (f_confOpts%has("IRwaterCoeff")) then
    call f_confOpts%get_or_die("IRwaterCoeff",str)
    IRwaterCoeff = str
 end if

 VISwaterCoeff = "NPOESS"
 if (f_confOpts%has("VISwaterCoeff")) then
    call f_confOpts%get_or_die("VISwaterCoeff",str)
    VISwaterCoeff = str
 end if
 IRVISlandCoeff = "NPOESS"
 if (f_confOpts%has("IRVISlandCoeff")) then
    call f_confOpts%get_or_die("IRVISlandCoeff",str)
    IRVISlandCoeff = str
 end if
 IRVISsnowCoeff = "NPOESS"
 if (f_confOpts%has("IRVISsnowCoeff")) then
    call f_confOpts%get_or_die("IRVISsnowCoeff",str)
    IRVISsnowCoeff = str
 end if
 IRVISiceCoeff = "NPOESS"
 if (f_confOpts%has("IRVISiceCoeff")) then
    call f_confOpts%get_or_die("IRVISiceCoeff",str)
    IRVISiceCoeff = str
 end if
 MWwaterCoeff = "FASTEM6"
 if (f_confOpts%has("MWwaterCoeff")) then
    call f_confOpts%get_or_die("MWwaterCoeff",str)
    MWwaterCoeff = str
 end if

 ! Define water, snow, ice (WSI) categories
 select case (trim(IRVISlandCoeff))
    case ("USGS")
       allocate(conf%Land_WSI(2))
       conf%Land_WSI(1:2) = (/16,24/)
    case ("IGBP")
       allocate(conf%Land_WSI(2))
       conf%Land_WSI(1:2) = (/15,17/)
    case ("NPOESS")
       allocate(conf%Land_WSI(1))
       conf%Land_WSI(1) = -1
    case default
       write(message,*) trim(routine_name), ' error: unknown IR/vis land coeff classification ', &
                        trim(IRVISlandCoeff)
       call abor1_ftn(message)
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
 endif

end subroutine crtm_conf_setup

! -----------------------------------------------------------------------------

subroutine crtm_conf_delete(conf)

implicit none
type(crtm_conf), intent(inout) :: conf

 deallocate(conf%SENSOR_ID)
 deallocate(conf%Land_WSI)
 deallocate(conf%Absorbers)
 deallocate(conf%Absorber_Id)
 deallocate(conf%Absorber_Units)
 deallocate(conf%Clouds)
 deallocate(conf%Cloud_Id)

end subroutine crtm_conf_delete

! ------------------------------------------------------------------------------

subroutine crtm_comm_stat_check(stat, PROGRAM_NAME, message, f_comm)
use fckit_mpi_module,   only: fckit_mpi_comm

implicit none

integer,              intent(in) :: stat
character(*),         intent(in) :: PROGRAM_NAME
character(*),         intent(in) :: message
type(fckit_mpi_comm), intent(in) :: f_comm

character(max_string) :: rank_message

 if ( stat /= SUCCESS ) THEN
    write(rank_message,*) trim(message)," on rank ",f_comm%rank()
    call Display_Message( PROGRAM_NAME, rank_message, FAILURE )
    call abor1_ftn("Abort from "//PROGRAM_NAME)
 end if

end subroutine crtm_comm_stat_check

! ------------------------------------------------------------------------------

subroutine ufo_crtm_skip_profiles(n_Profiles,n_Channels,channels,obss,atm,sfc,  &
                                  Is_Active_Sensor, Is_Vis_or_UV, Options)
! Profiles are skipped when the ObsValue of all channels is missing, if the
! pressure doesn't increase with levels, and if SST is missing for profiles with
! non-zero sea coverage.
! TODO: Use complete QC information
! It would be more comprehensive to use EffectiveQC or EffectiveError. That
! would require those ObsSpace values to be initialized before calls to
! this subroutine within ufo_radiancecrtm_simobs+ufo_radiancecrtm_tlad_settraj.
use missing_values_mod

implicit none
integer,              intent(in)    :: n_Profiles, n_Channels
type(c_ptr), value,   intent(in)    :: obss
integer(c_int),       intent(in)    :: channels(:)
type(CRTM_Atmosphere_type), intent(in) :: atm(:)
type(CRTM_Surface_type),    intent(in) :: sfc(:)
logical, intent(in):: Is_Active_Sensor
logical, intent(in):: Is_Vis_or_UV
type(CRTM_Options_type),    intent(inout) :: Options(:)

integer :: jprofile, jchannel, jlevel
character(len=MAXVARLEN) :: varname
character(len=64) :: obsGroupName
character(len=max_string) :: message
real(kind_real)  :: ObsVal(n_Profiles,n_Channels)

real(c_double)  :: missing_d
real(kind_real) :: missing_r
real(kind_real), parameter  :: lowest_albedo = 0.001

 ! Set missing values
 missing_d = missing_value(missing_d)
 missing_r = missing_value(missing_r)

 ObsVal = missing_d

 ! Do a quick test to set the group name for ObsValue
 call get_var_name(channels(1),varname, Is_Active_Sensor, Is_Vis_or_UV)
 if (obsspace_has(obss, "DerivedObsValue", varname)) then
   obsGroupName = "DerivedObsValue"
 elseif (obsspace_has(obss, "ObsValue", varname)) then
   obsGroupName = "ObsValue"
 else
   write(message,*) 'Group name for observed values is neither ObsValue nor DerivedObsValue'
   call abor1_ftn(message)
 endif

 do jchannel = 1, n_Channels
   call get_var_name(channels(jchannel),varname, Is_Active_Sensor, Is_Vis_or_UV)
   call obsspace_get_db(obss, trim(obsGroupName), varname, ObsVal(:,jchannel))
 enddo

 !Loop over all n_Profiles, i.e. number of locations
 profile_loop: do jprofile = 1, n_Profiles
   ! check whether observations for all channels are missing in the input file
   Options(jprofile)%Skip_Profile = all(ObsVal(jprofile,:) == missing_d)

   ! check for pressure monotonicity
   do jlevel = atm(jprofile)%n_layers, 1, -1
     if ( atm(jprofile)%level_pressure(jlevel) <= atm(jprofile)%level_pressure(jlevel-1) ) then
       Options(jprofile)%Skip_Profile = .TRUE.
       cycle profile_loop
     end if
   end do

   ! check for missing values in water surface temperature when the mask
   ! indicates there is water.
   ! (this can happen with generically coupled atm and ocean components if
   ! the land/sea masks don't match)
   if ((sfc(jprofile)%Water_Temperature == missing_r) .and.   &
       (sfc(jprofile)%Water_Coverage > 0.0) .and. (.not. Is_Vis_or_UV)) then
     Options(jprofile)%Skip_Profile = .TRUE.
   endif

   ! check for all channels in Vis/UV profiles that have ObsValue/albedo
   ! that are below minimum threshold. Skip those.
   if (Is_Vis_or_UV) then
     Options(jprofile)%Skip_Profile = all(ObsVal(jprofile,:) < lowest_albedo)
   endif

   ! check for all channels in active profiles that have ObsValue/Reflectivity
   ! that are beyond threshold. Skip those.
   if (Is_Active_Sensor) then
     ! the second dimension is for channels so if any channel is missing then skip it
     Options(jprofile)%Skip_Profile = any(abs(ObsVal(jprofile,:)) >= threshold_reflectivity)
   endif
 end do profile_loop

end subroutine ufo_crtm_skip_profiles

! ------------------------------------------------------------------------------

SUBROUTINE Load_Atm_Data(n_Profiles, n_Layers, geovals, atm, conf, Is_Active_Sensor)
implicit none

integer, intent(in) :: n_Profiles, n_Layers
type(ufo_geovals), intent(in) :: geovals
type(CRTM_Atmosphere_type), intent(inout) :: atm(:)
logical, intent(in), optional :: Is_Active_Sensor
type(crtm_conf) :: conf

! Local variables
integer :: k1, jspec, jlevel
type(ufo_geoval), pointer :: geoval
character(max_string) :: err_msg
logical :: IsActiveSensor

  if (present(Is_Active_Sensor)) then
     IsActiveSensor = Is_Active_Sensor
  else
     IsActiveSensor = .FALSE.
  endif

  ! Populate the atmosphere structures for CRTM
  ! -------------------------------------------

  call ufo_geovals_get_var(geovals, var_ts, geoval)
  ! Check model levels is consistent in geovals & crtm
  if (geoval%nval /= n_Layers) then
    write(err_msg,*) 'Load_Atm_Data error: layers inconsistent!'
    call abor1_ftn(err_msg)
  endif

  do k1 = 1, n_Profiles
    atm(k1)%Temperature(1:n_Layers) = geoval%vals(:, k1)
  end do

  call ufo_geovals_get_var(geovals, var_prs, geoval)
  do k1 = 1, n_Profiles
    atm(k1)%Pressure(1:n_Layers) = geoval%vals(:, k1) * 0.01  ! to hPa
  end do

  call ufo_geovals_get_var(geovals, var_prsi, geoval)
  do k1 = 1, n_Profiles
    atm(k1)%Level_Pressure(:) = geoval%vals(:, k1) * 0.01     ! to hPa
    atm(k1)%Climatology = US_STANDARD_ATMOSPHERE
  end do

  if ((conf%Aerosol_Model == 'GOCART-GEOS5') .or. &
      (conf%Aerosol_Model == 'NAAPS')) then
    call ufo_geovals_get_var(geovals, var_rh, geoval)
    do k1 = 1, n_Profiles
      WHERE (geoval%vals(:, k1) > 1.0_kind_real) geoval%vals(:, k1) = 1.0_kind_real
      atm(k1)%Relative_Humidity(:) = geoval%vals(:, k1)        ! fraction
      atm(k1)%Climatology = US_STANDARD_ATMOSPHERE
    end do
  endif

  do jspec = 1, conf%n_Absorbers
    ! O3 Absorber has special treatment for Aerosols
    if (cmp_strings(conf%Absorbers(jspec), var_oz) .AND. &
      ufo_vars_getindex(geovals%variables, var_oz) < 0 .AND. &
      (.NOT. cmp_strings(conf%aerosol_option,""))) then
      do k1 = 1, n_Profiles
        atm(k1)%Absorber(1:n_Layers, jspec) = ozone_default_value
      end do
    else
      call ufo_geovals_get_var(geovals, conf%Absorbers(jspec), geoval)
      do k1 = 1, n_Profiles
        atm(k1)%Absorber(1:n_Layers, jspec) = geoval%vals(:, k1)
      end do
    end if
    do k1 = 1, n_Profiles
      atm(k1)%Absorber_Id(jspec) = conf%Absorber_Id(jspec)
      atm(k1)%Absorber_Units(jspec) = conf%Absorber_Units(jspec)
    end do
  end do

  do jspec = 1, conf%n_Clouds
    ! cloud species content
    CALL ufo_geovals_get_var(geovals, conf%Clouds(jspec,1), geoval)
    do k1 = 1, n_Profiles
      atm(k1)%Cloud(jspec)%Water_Content = geoval%vals(:, k1)
      atm(k1)%Cloud(jspec)%Type = conf%Cloud_Id(jspec)
    end do

    ! effective radius
    CALL ufo_geovals_get_var(geovals, conf%Clouds(jspec,2), geoval)
    do k1 = 1, n_Profiles
      atm(k1)%Cloud(jspec)%Effective_Radius = geoval%vals(:, k1)
    end do
  end do

  ! When n_Clouds>0, Cloud_Fraction must either be provided as geoval or in conf
  ! If Cloud_Fraction is provided in conf,then use the Cloud_Fraction in conf
  ! If Cloud_Fraction is not provided in conf , then use the Cloud_Fraction in geoval
  ! and make sure the cloud fraction interpolated from background is in the valid range [0.0,1.0]
  if (conf%n_Clouds > 0) then
    if ( conf%Cloud_Fraction >= 0.0  ) then
      do k1 = 1, n_Profiles
        atm(k1)%Cloud_Fraction(:) = conf%Cloud_Fraction
      end do
    else
      if ( ufo_vars_getindex(geovals%variables, var_cldfrac) > 0 ) then
        CALL ufo_geovals_get_var(geovals, var_cldfrac, geoval)
        do k1 = 1, n_Profiles
          where( geoval%vals(:, k1) < 0.0_kind_real ) geoval%vals(:, k1) = 0.0_kind_real
          where( geoval%vals(:, k1) > 1.0_kind_real ) geoval%vals(:, k1) = 1.0_kind_real
          atm(k1)%Cloud_Fraction(:) =  geoval%vals(:, k1)
        end do
      end if
    end if
  end if


  if ( (conf%n_Clouds > 0) .and. (.NOT. IsActiveSensor) ) then
    if ( conf%Cloud_Seeding ) then 
      do k1 = 1, n_Profiles
        do jlevel = 1, atm(k1)%n_layers
           ! Check Cloud Content
           do jspec = 1, conf%n_Clouds
             if (atm(k1)%Cloud(jspec)%Type == WATER_CLOUD .and. atm(k1)%Temperature(jlevel) - tice > -20.0_kind_real ) then
                 atm(k1)%Cloud(jspec)%Water_Content(jlevel) = max(atm(k1)%Cloud(jspec)%Water_Content(jlevel), WATER_CONTENT_THRESHOLD*1.001_kind_real )
                 atm(k1)%Cloud(jspec)%Effective_Radius(jlevel) = max(atm(k1)%Cloud(jspec)%Effective_Radius(jlevel), 5.001_kind_real)
             end if
             if (atm(k1)%Cloud(jspec)%Type == ICE_CLOUD .and. atm(k1)%Temperature(jlevel) - tice < 0.0_kind_real ) then
                 atm(k1)%Cloud(jspec)%Water_Content(jlevel) = max(atm(k1)%Cloud(jspec)%Water_Content(jlevel), WATER_CONTENT_THRESHOLD*1.001_kind_real )
                 atm(k1)%Cloud(jspec)%Effective_Radius(jlevel) = max(atm(k1)%Cloud(jspec)%Effective_Radius(jlevel), 5.001_kind_real)
             end if
             if (atm(k1)%Cloud(jspec)%Type == RAIN_CLOUD .and. atm(k1)%Temperature(jlevel) - tice > -20.0_kind_real ) then
                 atm(k1)%Cloud(jspec)%Water_Content(jlevel) = max(atm(k1)%Cloud(jspec)%Water_Content(jlevel), WATER_CONTENT_THRESHOLD*1.001_kind_real )
                 atm(k1)%Cloud(jspec)%Effective_Radius(jlevel) = max(atm(k1)%Cloud(jspec)%Effective_Radius(jlevel), 100.001_kind_real)
             end if
             if (atm(k1)%Cloud(jspec)%Type == SNOW_CLOUD .and. atm(k1)%Temperature(jlevel) - tice < 0.0_kind_real ) then
                 atm(k1)%Cloud(jspec)%Water_Content(jlevel) = max(atm(k1)%Cloud(jspec)%Water_Content(jlevel), WATER_CONTENT_THRESHOLD*1.001_kind_real )
                 atm(k1)%Cloud(jspec)%Effective_Radius(jlevel) = max(atm(k1)%Cloud(jspec)%Effective_Radius(jlevel), 50.001_kind_real)
             end if
             if (atm(k1)%Cloud(jspec)%Type == GRAUPEL_CLOUD .and. atm(k1)%Temperature(jlevel) - tice < 0.0_kind_real) then
                 atm(k1)%Cloud(jspec)%Water_Content(jlevel) = max(atm(k1)%Cloud(jspec)%Water_Content(jlevel), WATER_CONTENT_THRESHOLD*1.001_kind_real )
                 atm(k1)%Cloud(jspec)%Effective_Radius(jlevel) = max(atm(k1)%Cloud(jspec)%Effective_Radius(jlevel), 500.001_kind_real)
             end if
             if (atm(k1)%Cloud(jspec)%Type == HAIL_CLOUD .and. atm(k1)%Temperature(jlevel) - tice < 0.0_kind_real) then
                 atm(k1)%Cloud(jspec)%Water_Content(jlevel) = max(atm(k1)%Cloud(jspec)%Water_Content(jlevel), WATER_CONTENT_THRESHOLD*1.001_kind_real )
                 atm(k1)%Cloud(jspec)%Effective_Radius(jlevel) = max(atm(k1)%Cloud(jspec)%Effective_Radius(jlevel), 1000.001_kind_real)
             end if
           end do
           ! Check Cloud Fraction
           do jspec = 1, conf%n_Clouds
             if (atm(k1)%Cloud(jspec)%Water_Content(jlevel) >= WATER_CONTENT_THRESHOLD*1.001_kind_real .and. atm(k1)%Cloud_Fraction(jlevel) < 0.001_kind_real) then
               atm(k1)%Cloud_Fraction(jlevel) = 0.001_kind_real
             end if
           end do
        end do
      end do
    end if
  end if

  if (IsActiveSensor) then
       Atm%Add_Extra_Layers = .FALSE.
  end if

end subroutine Load_Atm_Data

! ------------------------------------------------------------------------------

subroutine Load_Sfc_Data(n_Profiles, n_Channels, channels, geovals, sfc, chinfo, obss, conf,  &
                         Is_Active_Sensor, Is_Vis_or_UV)

implicit none

integer,                     intent(in)    :: n_Profiles, n_Channels
type(ufo_geovals),           intent(in)    :: geovals
type(CRTM_Surface_type),     intent(inout) :: sfc(:)
type(CRTM_ChannelInfo_type), intent(in)    :: chinfo(:)
type(c_ptr), value,          intent(in)    :: obss
integer(c_int),              intent(in)    :: channels(:)
type(crtm_conf),             intent(in)    :: conf
logical, intent(in) :: Is_Active_Sensor
logical, intent(in) :: Is_Vis_or_UV

type(ufo_geoval), pointer :: geoval, u, v
integer :: k1, n1
integer :: iLand
logical :: use_mw_vegtyp_soiltyp_data, use_visir_landtyp_data

character(len=MAXVARLEN) :: varname
character(len=max_string) :: message
character(len=64) :: obsGroupName

real(kind_real), allocatable :: ObsTb(:,:)

  allocate(ObsTb(n_profiles, n_channels))
  ObsTb = 0.0_kind_real

  ! Do a quick test to set the group name for ObsValue
  call get_var_name(channels(1),varname, Is_Active_Sensor, Is_Vis_or_UV)
  if (obsspace_has(obss, "ObsValue", varname)) then
    obsGroupName = "ObsValue"
  elseif (obsspace_has(obss, "DerivedObsValue", varname)) then
    obsGroupName = "DerivedObsValue"
  else
    write(message,*) 'Group name for observed values is neither ObsValue nor DerivedObsValue'
    call abor1_ftn(message)
  endif

  if (.not. Is_Active_Sensor) then
    do n1 = 1, n_Channels
      call get_var_name(channels(n1),varname, Is_Active_Sensor, Is_Vis_or_UV)
      call obsspace_get_db(obss, trim(obsGroupName), varname, ObsTb(:, n1))
    enddo
  end if

  do k1 = 1, n_Profiles
    !Pass sensor information
    sfc(k1)%sensordata%sensor_id        = chinfo(1)%sensor_id
    sfc(k1)%sensordata%wmo_sensor_id    = chinfo(1)%wmo_sensor_id
    sfc(k1)%sensordata%wmo_satellite_id = chinfo(1)%wmo_satellite_id
    sfc(k1)%sensordata%sensor_channel   = channels

    !Pass observation value
    if (.not. Is_Active_Sensor) then
      do n1 = 1, n_channels
         sfc(k1)%sensordata%tb(n1) = ObsTb(k1, n1)
      enddo
    end if

    !Water_type
    !** NOTE: need to check how to determine fresh vs sea water types (salinity???)
    sfc(k1)%Water_Type = 1  ! SEA_WATER_TYPE for all SfcOptics
  end do
  deallocate(ObsTb)

  use_mw_vegtyp_soiltyp_data = any(chinfo%Sensor_Type == MICROWAVE_SENSOR)
  use_visir_landtyp_data = any(chinfo%Sensor_Type == VISIBLE_SENSOR) .or. &
                           any(chinfo%Sensor_Type == INFRARED_SENSOR)

  if (ufo_vars_getindex(geovals%variables, var_sfc_wspeed) > 0 .and. &
      ufo_vars_getindex(geovals%variables, var_sfc_wdir) > 0) then
    ! Directly use model-provided wind speed and direction
    !Wind_Speed
    call ufo_geovals_get_var(geovals, var_sfc_wspeed, geoval)
    do k1 = 1, n_Profiles
      sfc(k1)%Wind_Speed = geoval%vals(1, k1)
    end do

    !Wind_Direction
    call ufo_geovals_get_var(geovals, var_sfc_wdir, geoval)
    do k1 = 1, n_Profiles
      sfc(k1)%Wind_Direction = geoval%vals(1, k1)
    end do
  else if (ufo_vars_getindex(geovals%variables, var_sfc_u) > 0 .and. &
      ufo_vars_getindex(geovals%variables, var_sfc_v) > 0) then
    ! Convert 2d wind components to speed and direction
    call ufo_geovals_get_var(geovals, var_sfc_u, u)
    call ufo_geovals_get_var(geovals, var_sfc_v, v)

    !Wind_Speed
    do k1 = 1, n_Profiles
      sfc(k1)%Wind_Speed = sqrt(u%vals(1, k1)**2 + v%vals(1, k1)**2)
    end do

    !Wind_Direction
    do k1 = 1, n_Profiles
      sfc(k1)%Wind_Direction = uv_to_wdir(u%vals(1, k1), v%vals(1, k1))
    end do
  else
    call abor1_ftn('Load_Sfc_Data error: missing surface wind geovals')
  end if

  !Water_Coverage
  call ufo_geovals_get_var(geovals, var_sfc_wfrac, geoval)
  do k1 = 1, n_Profiles
    sfc(k1)%Water_Coverage = geoval%vals(1, k1)
  end do

  !Water_Temperature
  call ufo_geovals_get_var(geovals, var_sfc_wtmp, geoval)
  do k1 = 1, n_Profiles
    sfc(k1)%Water_Temperature = geoval%vals(1, k1)
  end do

  !Ice_Coverage
  call ufo_geovals_get_var(geovals, var_sfc_ifrac, geoval)
  do k1 = 1, n_Profiles
    sfc(k1)%Ice_Coverage = geoval%vals(1, k1)
  end do

  !Ice_Temperature
  call ufo_geovals_get_var(geovals, var_sfc_itmp, geoval)
  do k1 = 1, n_Profiles
    sfc(k1)%Ice_Temperature = geoval%vals(1, k1)
  end do

  !Snow_Coverage
  call ufo_geovals_get_var(geovals, var_sfc_sfrac, geoval)
  do k1 = 1, n_Profiles
    sfc(k1)%Snow_Coverage = geoval%vals(1, k1)
  end do

  !Snow_Temperature
  call ufo_geovals_get_var(geovals, var_sfc_stmp, geoval)
  do k1 = 1, n_Profiles
    sfc(k1)%Snow_Temperature = geoval%vals(1, k1)
  end do

  !Snow_Depth
  call ufo_geovals_get_var(geovals, var_sfc_sdepth, geoval)
  do k1 = 1, n_Profiles
    sfc(k1)%Snow_Depth = geoval%vals(1, k1)
  end do

  !Land_Coverage
  call ufo_geovals_get_var(geovals, var_sfc_lfrac, geoval)
  do k1 = 1, n_Profiles
    sfc(k1)%Land_Coverage = geoval%vals(1, k1)
  end do

  !Land_Type
  ! + used to lookup land sfc emiss. for IR and VIS
  ! + land sfc emiss. undefined over water/snow/ice
  ! + note that land type classifiation is an option
  if (use_visir_landtyp_data) then
    if (index(conf%IRlandCoeff_File, "NPOESS") == 1) then
      ! this is the default if CRTM not configured otherwise
      call ufo_geovals_get_var(geovals, var_sfc_landtyp_npoess, geoval)
    else if (index(conf%IRlandCoeff_File, "IGBP") == 1) then
      call ufo_geovals_get_var(geovals, var_sfc_landtyp_igbp, geoval)
    else if (index(conf%IRlandCoeff_File, "USGS") == 1) then
      call ufo_geovals_get_var(geovals, var_sfc_landtyp_usgs, geoval)
    else
      write(message,*) "Load_Sfc_Data error: cannot infer land type classification from " &
                       // "IRlandCoeff_File " // trim(conf%IRlandCoeff_File)
      call abor1_ftn(message)
    end if

    do k1 = 1, n_Profiles
      iLand = int(geoval%vals(1, k1))
      if (.not.any(iLand == conf%Land_WSI)) then
        sfc(k1)%Land_Type = iLand
      end if
    end do
  end if

  !Land_Temperature
  call ufo_geovals_get_var(geovals, var_sfc_ltmp, geoval)
  do k1 = 1, n_Profiles
    sfc(k1)%Land_Temperature = geoval%vals(1, k1)
  end do

  !Lai
  call ufo_geovals_get_var(geovals, var_sfc_lai, geoval)
  do k1 = 1, n_Profiles
    sfc(k1)%Lai = geoval%vals(1, k1)
  end do

  !Vegetation_Fraction
  call ufo_geovals_get_var(geovals, var_sfc_vegfrac, geoval)
  do k1 = 1, n_Profiles
    sfc(k1)%Vegetation_Fraction = geoval%vals(1, k1)
  end do

  if (use_mw_vegtyp_soiltyp_data) then
    !Vegetation_Type
    call ufo_geovals_get_var(geovals, var_sfc_vegtyp, geoval)
    do k1 = 1, n_Profiles
      sfc(k1)%Vegetation_Type = int(geoval%vals(1, k1))
    end do

    !Soil_Type
    call ufo_geovals_get_var(geovals, var_sfc_soiltyp, geoval)
    do k1 = 1, n_Profiles
      sfc(k1)%Soil_Type = int(geoval%vals(1, k1))
    end do
  end if

  !Soil_Moisture_Content
  call ufo_geovals_get_var(geovals, var_sfc_soilm, geoval)
  do k1 = 1, n_Profiles
    sfc(k1)%Soil_Moisture_Content = geoval%vals(1, k1)
  end do

  !Soil_Temperature
  call ufo_geovals_get_var(geovals, var_sfc_soilt, geoval)
  do k1 = 1, n_Profiles
    sfc(k1)%Soil_Temperature = geoval%vals(1, k1)
  end do

  !Sea_Surface_Salinity
  if (cmp_strings(conf%salinity_option, "on")) THEN
    call ufo_geovals_get_var(geovals, var_sfc_sss, geoval)
    do k1 = 1, n_Profiles
      sfc(k1)%Salinity = geoval%vals(1, k1)
    end do
  end if

  ! Update Ice_Coverage and Land_Coverage for
  ! Soil_Type or Vegetation_Type corresponding to Glacial land ice
  do k1 = 1, n_Profiles
    if (sfc(k1)%Land_Coverage > ZERO .and. &
       (sfc(k1)%Soil_Type == 9 .or. sfc(k1)%Vegetation_Type == 13)) then
      sfc(k1)%Ice_Coverage = min(sfc(k1)%Ice_Coverage + sfc(k1)%Land_Coverage, ONE)
      sfc(k1)%Land_Coverage = ZERO
    end if
  end do

end subroutine Load_Sfc_Data

! ------------------------------------------------------------------------------

subroutine Load_Geom_Data(obss,geo,geo_hf,sensor_id)

implicit none
type(c_ptr), value,       intent(in)    :: obss
type(CRTM_Geometry_type), intent(inout) :: geo(:)
type(CRTM_Geometry_type), intent(inout), optional :: geo_hf(:)
character(len=255),       intent(in),    optional :: sensor_id
real(kind_real), allocatable :: TmpVar(:)
integer, allocatable :: TmpVar2(:)
integer :: nlocs

 nlocs = obsspace_get_nlocs(obss)
 allocate(TmpVar(nlocs))
 allocate(TmpVar2(nlocs))

 call obsspace_get_db(obss, "MetaData", "sensorZenithAngle", TmpVar)
 geo(:)%Sensor_Zenith_Angle = abs(TmpVar(:)) ! needs to be absolute value
 where (geo(:)%Sensor_Zenith_Angle > 80.0_kind_real) &
    geo(:)%Sensor_Zenith_Angle = 80.0_kind_real

 call obsspace_get_db(obss, "MetaData", "solarZenithAngle", TmpVar)
 geo(:)%Source_Zenith_Angle = TmpVar(:)

 call obsspace_get_db(obss, "MetaData", "sensorAzimuthAngle", TmpVar)
 geo(:)%Sensor_Azimuth_Angle = TmpVar(:)

 call obsspace_get_db(obss, "MetaData", "solarAzimuthAngle", TmpVar)
 geo(:)%Source_Azimuth_Angle = TmpVar(:)

!  For some microwave instruments the solar and sensor azimuth angles can be
!  missing  (given a value of 10^11).  Set these to zero to get past CRTM QC.
 where (geo(:)%Source_Azimuth_Angle < 0.0_kind_real .or. &
        geo(:)%Source_Azimuth_Angle > 360.0_kind_real) &
    geo(:)%Source_Azimuth_Angle = 0.0_kind_real
 where (geo(:)%Sensor_Azimuth_Angle < 0.0_kind_real .or. &
        geo(:)%Sensor_Azimuth_Angle > 360.0_kind_real) &
    geo(:)%Sensor_Azimuth_Angle = 0.0_kind_real

 where (abs(geo(:)%Source_Zenith_Angle) > 180.0_kind_real) &
    geo(:)%Source_Zenith_Angle = 100.0_kind_real

 call obsspace_get_db(obss, "MetaData", "sensorScanPosition", TmpVar2)
 geo(:)%Ifov = TmpVar2(:)

 call obsspace_get_db(obss, "MetaData", "sensorViewAngle", TmpVar) !The Sensor_Scan_Angle is optional
 geo(:)%Sensor_Scan_Angle = TmpVar(:)

 where (abs(geo(:)%Sensor_Scan_Angle) > 80.0_kind_real) &
    geo(:)%Sensor_Scan_Angle = 0.0_kind_real

 ! Read geophysical values for gmi high frequency channels 10-13.
 if (present(sensor_id)) then
   if (cmp_strings(trim(sensor_id),'gmi_gpm')) then
    if ( present(geo_hf) ) then
      geo_hf = geo
      if (obsspace_has(obss, "MetaData", "sensorZenithAngle1")) then
        call obsspace_get_db(obss, "MetaData", "sensorZenithAngle1", TmpVar)
        geo_hf(:)%Sensor_Zenith_Angle = abs(TmpVar(:)) ! needs to be absolute value
      endif
      if (obsspace_has(obss, "MetaData", "solarZenithAngle1")) then
        call obsspace_get_db(obss, "MetaData", "solarZenithAngle1", TmpVar)
        geo_hf(:)%Source_Zenith_Angle = TmpVar(:)
      endif
      if (obsspace_has(obss, "MetaData", "sensorAzimuthAngle1")) then
        call obsspace_get_db(obss, "MetaData", "sensorAzimuthAngle1", TmpVar)
        geo_hf(:)%Sensor_Azimuth_Angle = TmpVar(:)
      endif
      if (obsspace_has(obss, "MetaData", "solarAzimuthAngle1")) then
        call obsspace_get_db(obss, "MetaData", "solarAzimuthAngle1", TmpVar)
        geo_hf(:)%Source_Azimuth_Angle = TmpVar(:)
      endif
      if (obsspace_has(obss, "MetaData", "sensorViewAngle1")) then
        call obsspace_get_db(obss, "MetaData", "sensorViewAngle1", TmpVar)
        geo_hf(:)%Sensor_Scan_Angle = TmpVar(:)
      endif
      ! For some microwave instruments the solar and sensor azimuth angles can be
      ! missing  (given a value of 10^11).  Set these to zero to get past CRTM QC.
      where (geo_hf(:)%Source_Azimuth_Angle < 0.0_kind_real .or. &
            geo_hf(:)%Source_Azimuth_Angle > 360.0_kind_real) &
        geo_hf(:)%Source_Azimuth_Angle = 0.0_kind_real
      where (geo_hf(:)%Sensor_Azimuth_Angle < 0.0_kind_real .or. &
            geo_hf(:)%Sensor_Azimuth_Angle > 360.0_kind_real) &
        geo_hf(:)%Sensor_Azimuth_Angle = 0.0_kind_real
    endif
  endif
 endif

 deallocate(TmpVar)
 deallocate(TmpVar2)

end subroutine Load_Geom_Data

! ------------------------------------------------------------------------------

subroutine get_var_name(n,varname, Is_Active_Sensor, Is_Vis_or_UV)

integer, intent(in) :: n
character(len=*), intent(out) :: varname
logical, intent(in) :: Is_Active_Sensor
logical, intent(in) :: Is_Vis_or_UV

character(len=6) :: chan

 write(chan, '(I0)') n
 if (Is_Active_Sensor) then
     varname = 'ReflectivityAttenuated_' // trim(chan)
 else
     if (Is_Vis_or_UV) then
        varname = 'albedo_' // trim(chan)
     else
        varname = 'brightnessTemperature_' // trim(chan)
     endif
 endif 

end subroutine get_var_name

! -----------------------------------------------------------------------------

!> \brief Determines the wind direction from U and V components
!!
!! \details **uv_to_wdir** Calculates the wind direction, as measured clockwise
!! from north, similar to an azimuth angle.  Takes the eastward and northward
!! wind component magnitudes, respectively, as arguments.
!! Due to the azimuthal convention used here, the inverse equations are:
!! u = w * cos(wdir * deg2rad)
!! v = w * sin(wdir * deg2rad)
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
  reshape((/zero,  one,  one,  two, &
            one,  -one,  one, -one/), (/4,2/))

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
    endif
  endif
  windangle = atan(abs(windratio))   ! wind azimuth is in radians
  wdir = ( quadcof(iquadrant, 1) * pi + windangle * quadcof(iquadrant, 2) ) * rad2deg

end function uv_to_wdir

! -----------------------------------------------------------------------------

   SUBROUTINE define_aerosol_model(aerosol_coef_file, aerosol_model)
    ! Possible models are:
    ! CRTM, NAAPS, GOCART-GEOS5, CMAQ
    CHARACTER(*), INTENT(in) :: aerosol_coef_file
    CHARACTER(*), INTENT(out) :: aerosol_model
    integer :: checkstring

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
    endif

   END SUBROUTINE define_aerosol_model

   SUBROUTINE assign_aerosol_names(aerosol_option, var_aerosols)

    CHARACTER(*), INTENT(in) :: aerosol_option
    CHARACTER(len=MAXVARLEN), ALLOCATABLE, INTENT(out) :: var_aerosols(:)

    CHARACTER(max_string) :: err_msg

    IF (cmp_strings(aerosol_option,"aerosols_gocart_default")) THEN
       ALLOCATE(var_aerosols(n_aerosols_gocart_default))
       var_aerosols=var_aerosols_gocart_default
    ELSEIF (cmp_strings(aerosol_option,"aerosols_gocart_gefs")) THEN
       ALLOCATE(var_aerosols(n_aerosols_gocart_gefs))
       var_aerosols=var_aerosols_gocart_gefs
    ELSEIF (cmp_strings(aerosol_option,"aerosols_gocart_ufs")) THEN
       ALLOCATE(var_aerosols(n_aerosols_gocart_ufs))
       var_aerosols=var_aerosols_gocart_ufs
    ELSEIF (cmp_strings(aerosol_option,"aerosols_gocart_geos")) THEN
       ALLOCATE(var_aerosols(n_aerosols_gocart_geos))
       var_aerosols=var_aerosols_gocart_geos
    ELSE
       WRITE(err_msg,*) 'assign_aerosol_names: aerosol_option not implemented'&
       &//TRIM(aerosol_option)
       call abor1_ftn(err_msg)
     END IF

   END SUBROUTINE assign_aerosol_names

   SUBROUTINE load_aerosol_data(n_profiles, n_layers, geovals,&
     &conf, var_aerosols, aerosol_model, atm)

    USE CRTM_aerosolcoeff, ONLY: aeroc

    TYPE(crtm_conf), INTENT(in)    :: conf
    TYPE(ufo_geovals), INTENT(in) :: geovals
    TYPE(CRTM_atmosphere_type), INTENT(inout) :: atm(:)
    TYPE(ufo_geoval), POINTER :: geoval

    INTEGER, INTENT(in) :: n_profiles, n_layers
    INTEGER :: ivar, n_aerosols, i, k, m

    REAL(kind_real), DIMENSION(5), PARAMETER  :: dust_radii=[&
         &0.55_kind_real,1.4_kind_real,2.4_kind_real,4.5_kind_real,8.0_kind_real]
    REAL(kind_real), DIMENSION(n_layers) :: layer_factors
    REAL(kind_real), DIMENSION(n_layers, n_profiles) :: rh

    CHARACTER(*), PARAMETER :: routine_name = 'Load_Aerosol_Data'
    CHARACTER(*), INTENT(in) :: aerosol_model
    CHARACTER(len=MAXVARLEN) :: var_aerosols(:)
    CHARACTER(len=MAXVARLEN) :: varname
    CHARACTER(max_string) :: err_msg, message


    varname = var_rh
    CALL ufo_geovals_get_var(geovals, varname, geoval)
    rh(1:n_layers,1:n_profiles)=geoval%vals(1:n_layers,1:n_profiles)
    WHERE (rh > 1.0_kind_real) rh=1.0_kind_real

    n_aerosols=SIZE(var_aerosols)

    DO m=1,n_profiles

       CALL calculate_aero_layer_factor(atm(m), layer_factors)

       DO i=1,n_aerosols

          varname=var_aerosols(i)
          CALL ufo_geovals_get_var(geovals,varname, geoval)

          atm(m)%aerosol(i)%Concentration(1:n_layers)=&
               &MAX(geoval%vals(:,m)*conf%unit_coef*layer_factors, &
               &aerosol_concentration_minvalue_layer)


          IF (cmp_strings(TRIM(aerosol_model), "CRTM")) THEN

            !Indices for CRTM default LUT
            !DUST_AEROSOL = 1
            !SEASALT_AEROSOL = 2 - 5
            !ORGANIC_CARBON_AEROSOL = 6
            !BLACK_CARBON_AEROSOL = 7
            !SULFATE_AEROSOL = 8
            SELECT CASE (TRIM(varname))

            CASE (var_sulfate)
               atm(m)%aerosol(i)%TYPE  = 8

            CASE (var_bcphobic)
               atm(m)%aerosol(i)%TYPE  = 7
            CASE (var_bcphilic)
               atm(m)%aerosol(i)%TYPE  = 7

            CASE (var_ocphobic)
               atm(m)%aerosol(i)%TYPE  = 6
            CASE (var_ocphilic)
               atm(m)%aerosol(i)%TYPE  = 6

            CASE (var_du001)
               atm(m)%aerosol(i)%TYPE  = 1
               atm(m)%aerosol(i)%effective_radius(:)=dust_radii(1)
            CASE (var_du002)
               atm(m)%aerosol(i)%TYPE  = 1
               atm(m)%aerosol(i)%effective_radius(:)=dust_radii(2)
            CASE (var_du003)
               atm(m)%aerosol(i)%TYPE  = 1
               atm(m)%aerosol(i)%effective_radius(:)=dust_radii(3)
            CASE (var_du004)
               atm(m)%aerosol(i)%TYPE  = 1
               atm(m)%aerosol(i)%effective_radius(:)=dust_radii(4)
            CASE (var_du005)
               atm(m)%aerosol(i)%TYPE  = 1
               atm(m)%aerosol(i)%effective_radius(:)=dust_radii(5)

            CASE (var_ss001)
               atm(m)%aerosol(i)%TYPE  = 2
            CASE (var_ss002)
               atm(m)%aerosol(i)%TYPE  = 3
            CASE (var_ss003)
               atm(m)%aerosol(i)%TYPE  = 4
            CASE (var_ss004)
               atm(m)%aerosol(i)%TYPE  = 5
            CASE (var_ss005)
               atm(m)%aerosol(i)%TYPE  = 5

            CASE DEFAULT
               write(message,*) 'WARNING!: ', TRIM(varname),&
               &' is not included in ', TRIM(aerosol_model), ' LUT'
               atm(m)%aerosol(i)%TYPE  = -1
            END SELECT

            SELECT CASE (TRIM(varname))

            CASE (var_sulfate, var_bcphilic, var_ocphilic,&
                 &var_ss001, var_ss002, var_ss003, var_ss004)
               DO k=1,n_layers

                  atm(m)%aerosol(i)%effective_radius(k)=&
                   & gocart_aerosol_size(atm(m)%aerosol(i)&
                   &%TYPE, rh(k,m))
               ENDDO

            CASE (var_bcphobic, var_ocphobic)
               atm(m)%aerosol(i)%effective_radius(:)= AeroC&
                &%Reff(1,atm(m)%aerosol(i)%TYPE)

            END SELECT


         ELSEIF (cmp_strings(TRIM(aerosol_model), "GOCART-GEOS5")) THEN

            ! This is for the NASA GOCART tables, aerosol scheme GOCART-GEOS5 in CRTM
            ! Reff are bin effective radius from NASA tables
            ! 1, 2, 3, 4, 5  = Dust 1, 2, 3, 4, 5;     with Reff 0.64, 1.32, 2.30, 4.17, 7.67 microns
            ! 6, 7, 8, 9, 10 = Sea salt 1, 2, 3, 4, 5; with Reff 0.08, 0.27, 1.07, 2.55, 7.34 microns
            ! 11, 12 = Organic carbon 1, 2;            with Reff 0.09 microns, hydrophobic and hydrophilic
            ! 13, 14 = Black carbon 1, 2;              with Reff 0.04 microns, hydrophobic and hydrophilic
            ! 15, 16 = Sulfate 1, 2;                   with Reff 0.16, 0.60 microns
            ! 17, 18, 19 = Nitrate 1, 2, 3;            with Reff 0.16, 2.10, 7.75 microns
            ! 20, 21 = Brown carbon 1,2;              (Pending release, as of Feb2024)

            ! Cheng: In this code block, the biggest difference from the function assign_gocart_default
            !        is that no rh adjustmenet to Reff is needed. It'll be interpolated within CRTM.
            !        Aerosol optical properties of NASA tables will not be interpolated over Reff
            !        dimension. There is no need to assign effective radius Reff.
            !        Note that the dust aerosol mapping with the defualt CRTM table
            !        was done outside of function assign_gocart_default, because
            !        dust are all assumed hydrophobic that no rh adjustmenet to Reff
            !        was needed.

            SELECT CASE (TRIM(varname))

            ! Dust
            CASE (var_du001)
               atm(m)%aerosol(i)%TYPE  = 1 ! dust bin 1
            CASE (var_du002)
               atm(m)%aerosol(i)%TYPE  = 2 ! dust bin 2
            CASE (var_du003)
               atm(m)%aerosol(i)%TYPE  = 3 ! dust bin 3
            CASE (var_du004)
               atm(m)%aerosol(i)%TYPE  = 4 ! dust bin 4
            CASE (var_du005)
               atm(m)%aerosol(i)%TYPE  = 5 ! dust bin 5

            ! Sea salt
             CASE (var_ss001)
               atm(m)%aerosol(i)%TYPE  = 6 ! sea salt bin 1
            CASE (var_ss002)
               atm(m)%aerosol(i)%TYPE  = 7 ! sea salt bin 2
            CASE (var_ss003)
               atm(m)%aerosol(i)%TYPE  = 8 ! sea salt bin 3
            CASE (var_ss004)
               atm(m)%aerosol(i)%TYPE  = 9 ! sea salt bin 4
            CASE (var_ss005)
               atm(m)%aerosol(i)%TYPE  = 10 ! sea salt bin 5

            ! Organic carbon
            CASE (var_ocphobic)
               atm(m)%aerosol(i)%TYPE  = 11 ! hydrophobic
            CASE (var_ocphilic)
               atm(m)%aerosol(i)%TYPE  = 12 ! hydrophilic

            ! Black carbon
            CASE (var_bcphobic)
               atm(m)%aerosol(i)%TYPE  = 13 ! hydrophobic
            CASE (var_bcphilic)
               atm(m)%aerosol(i)%TYPE  = 14 ! hydrophilic

            ! Sulfate
            CASE (var_sulfate)
               atm(m)%aerosol(i)%TYPE  = 15 ! bin 1
            !Jerome: There is only one bin of sulfate so far
            !CASE (var_sulfate)
            !   atm(m)%aerosol(i)%TYPE  = 16 ! bin 2

            ! Nitrate
            CASE (var_no3an1)
               atm(m)%aerosol(i)%TYPE  = 17 ! bin 1
            CASE (var_no3an2)
               atm(m)%aerosol(i)%TYPE  = 18 ! bin 2
            CASE (var_no3an3)
               atm(m)%aerosol(i)%TYPE  = 19 ! bin 3

            ! Brown carbon (pending release, comment out now)
            CASE (var_brphobic)
               atm(m)%aerosol(i)%TYPE  = 20 ! bin 1
            CASE (var_brphilic)
               atm(m)%aerosol(i)%TYPE  = 21 ! bin 2

            CASE DEFAULT
               write(message,*) 'WARNING!: ', TRIM(varname),&
               ' is not included in ', TRIM(aerosol_model), ' LUT'
               atm(m)%aerosol(i)%TYPE  = -1
            END SELECT

          ELSEIF (cmp_strings(trim(aerosol_model), "CMAQ")) THEN
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
            DO k=1,n_layers
               atm(m)%aerosol(i)%effective_variance(k)=1.0_kind_real
            ENDDO

            ! Assign aerosol type
            SELECT CASE (TRIM(varname))

              !Dust
              CASE (var_du001)
                 atm(m)%aerosol(i)%TYPE = 1
              CASE (var_du002)
                 atm(m)%aerosol(i)%TYPE = 1
              CASE (var_du003)
                 atm(m)%aerosol(i)%TYPE = 1
              CASE (var_du004)
                 atm(m)%aerosol(i)%TYPE = 1
              CASE (var_du005)
                 atm(m)%aerosol(i)%TYPE = 1

              !Sea Salt
              CASE (var_ss001)
                 atm(m)%aerosol(i)%TYPE = 5
              CASE (var_ss002)
                 atm(m)%aerosol(i)%TYPE = 5
              CASE (var_ss003)
                 atm(m)%aerosol(i)%TYPE = 5
              CASE (var_ss004)
                 atm(m)%aerosol(i)%TYPE = 5
              CASE (var_ss005)
                 atm(m)%aerosol(i)%TYPE = 5

              ! Organic carbon
              CASE (var_ocphobic)
                 atm(m)%aerosol(i)%TYPE = 7
              CASE (var_ocphilic)
                 atm(m)%aerosol(i)%TYPE = 3

              ! Black carbon
              CASE (var_bcphobic)
                 atm(m)%aerosol(i)%TYPE = 7
              CASE (var_bcphilic)
                 atm(m)%aerosol(i)%TYPE = 2

              ! Sulfate
              CASE (var_sulfate)
                 atm(m)%aerosol(i)%TYPE = 4

              ! Nitrate
              CASE (var_no3an1)
                 atm(m)%aerosol(i)%TYPE = 3
              CASE (var_no3an2)
                 atm(m)%aerosol(i)%TYPE = 3
              CASE (var_no3an3)
                 atm(m)%aerosol(i)%TYPE = 3

              !!Brown carbon
              CASE (var_brphobic)
                 atm(m)%aerosol(i)%TYPE = 7
              CASE (var_brphilic)
                 atm(m)%aerosol(i)%TYPE = 3

            CASE DEFAULT
               write(message,*) 'WARNING!: ', TRIM(varname),&
               ' is not included in ', TRIM(aerosol_model), ' LUT'
               atm(m)%aerosol(i)%TYPE  = -1
            WRITE(err_msg,*) TRIM(conf%aerosol_option)//' not ready in UFO/AODCRTM'
            call abor1_ftn(err_msg)

            END SELECT

            SELECT CASE (TRIM(varname))

            ! Reff for hydrophilic aerosols
            CASE (var_sulfate, var_bcphilic, var_ocphilic,&
                  &var_ss001, var_ss002, var_ss003, var_ss004,&
                  &var_no3an1, var_no3an2, var_no3an3)
               DO k=1,n_layers
                  atm(m)%aerosol(i)%effective_radius(k)=&
                   & gocart_aerosol_size(atm(m)%aerosol(i)&
                   &%TYPE, rh(k,m))
               ENDDO
            ! Reff for hydrophobic aerosols
            CASE (var_du001, var_du002, var_du003, var_du004, var_du005,&
                  &var_bcphobic, var_ocphobic)
               atm(m)%aerosol(i)%effective_radius(:)= AeroC&
                &%Reff(1,atm(m)%aerosol(i)%TYPE)

            END SELECT

          ELSEIF (cmp_strings(TRIM(aerosol_model), "NAAPS")) THEN
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
            SELECT CASE (TRIM(varname))

            !Dust
            CASE (var_du001)
               atm(m)%aerosol(i)%TYPE = 1
            CASE (var_du002)
               atm(m)%aerosol(i)%TYPE = 1
            CASE (var_du003)
               atm(m)%aerosol(i)%TYPE = 1
            CASE (var_du004)
               atm(m)%aerosol(i)%TYPE = 1
            CASE (var_du005)
               atm(m)%aerosol(i)%TYPE = 1

            !Sea Salt
            CASE (var_ss001)
              atm(m)%aerosol(i)%TYPE  = 3
            CASE (var_ss002)
              atm(m)%aerosol(i)%TYPE  = 3
            CASE (var_ss003)
              atm(m)%aerosol(i)%TYPE  = 3
            CASE (var_ss004)
              atm(m)%aerosol(i)%TYPE  = 3
            CASE (var_ss005)
              atm(m)%aerosol(i)%TYPE  = 3

            ! Organic carbon
            CASE (var_ocphobic)
               atm(m)%aerosol(i)%TYPE = 2
            CASE (var_ocphilic)
               atm(m)%aerosol(i)%TYPE = 2

            ! Black carbon
            CASE (var_bcphobic)
               atm(m)%aerosol(i)%TYPE = 2
            CASE (var_bcphilic)
               atm(m)%aerosol(i)%TYPE = 2

            ! Sulfate
            CASE (var_sulfate)
               atm(m)%aerosol(i)%TYPE = 4

            ! Nitrate
            CASE (var_no3an1)
               atm(m)%aerosol(i)%TYPE = 4
            CASE (var_no3an2)
               atm(m)%aerosol(i)%TYPE = 4
            CASE (var_no3an3)
               atm(m)%aerosol(i)%TYPE = 4

            !Brown carbon
            CASE (var_brphobic)
               atm(m)%aerosol(i)%TYPE = 2
            CASE (var_brphilic)
               atm(m)%aerosol(i)%TYPE = 2

            CASE DEFAULT
               write(message,*) 'WARNING!: ', TRIM(varname),&
               ' is not included in ', TRIM(aerosol_model), ' LUT'
               atm(m)%aerosol(i)%TYPE  = -1
            WRITE(err_msg,*) TRIM(conf%aerosol_option)//' not ready in UFO/AODCRTM'
            call abor1_ftn(err_msg)

            END SELECT

          ENDIF
       END DO
     END DO

   END SUBROUTINE load_aerosol_data

   SUBROUTINE calculate_aero_layer_factor_atm_profile(atm, layer_factors)

     TYPE(CRTM_atmosphere_type), INTENT(in) :: atm
     REAL(kind_real), INTENT(out) :: layer_factors(:)

     INTEGER :: k

     DO k=1,SIZE(layer_factors)
        !correct for mixing ratio factor layer_factors
        !being calculated from dry pressure, cotton eq. (2.4)
        !p_dry=p_total/(1+1.61*mixing_ratio)
        layer_factors(k)=1e-9_kind_real*(atm%Level_Pressure(k)-&
             &atm%Level_Pressure(k-1))*100.0_kind_real/grav/&
             &(1.0_kind_real+rv_rd*atm%Absorber(k,1)*1e-3_kind_real)
     ENDDO

   END SUBROUTINE calculate_aero_layer_factor_atm_profile

   SUBROUTINE calculate_aero_layer_factor_atm(atm, layer_factors)

     TYPE(CRTM_atmosphere_type), INTENT(in) :: atm(:)
     REAL(kind_real), INTENT(out) :: layer_factors(:,:)
     INTEGER :: k,m

     DO k=1,SIZE(layer_factors,1)
        DO m=1,SIZE(layer_factors,2)
           !correct for mixing ratio factor layer_factors
           !being calculated from dry pressure, cotton eq. (2.4)
           !p_dry=p_total/(1+1.61*mixing_ratio)
           layer_factors(k,m)=1e-9_kind_real*(atm(m)%Level_Pressure(k)-&
                &atm(m)%Level_Pressure(k-1))*100.0_kind_real/grav/&
                &(1.0_kind_real+rv_rd*atm(m)%Absorber(k,1)*1.e-3_kind_real)
        ENDDO
     ENDDO

   END SUBROUTINE calculate_aero_layer_factor_atm

   FUNCTION gocart_aerosol_size( itype, rh ) & ! rh input in 0-1
        &RESULT(r_eff)   ! in micrometer

     USE CRTM_aerosolcoeff, ONLY: aeroc
     IMPLICIT NONE

!
!   modified from a function provided by quanhua liu
!
     INTEGER ,INTENT(in) :: itype
     REAL(kind_real)    ,INTENT(in) :: rh

     INTEGER :: j1,j2,m
     INTEGER :: n_rh
     REAL(kind_real)    :: h1
     REAL(kind_real)    :: r_eff

     j2 = 0
     j1 = 1
     n_rh = 35    ! total number of RH in CRTM Default LUTs

     IF ( rh <= aeroc%rh(1) ) THEN
        j1 = 1
     ELSE IF ( rh >= aeroc%rh(n_rh) ) THEN
        j1 = n_rh
     ELSE
        DO m = 1, n_rh-1
           IF ( rh < aeroc%rh(m+1) .AND. rh > aeroc%rh(m) ) THEN
              j1 = m
              j2 = m+1
              h1 = (rh-aeroc%rh(m))/(aeroc%rh(m+1)-aeroc%rh(m))
              EXIT
           ENDIF
        ENDDO
     ENDIF

     IF ( j2 == 0 ) THEN
        r_eff = aeroc%reff(j1,itype )
     ELSE
        r_eff = (1.0_kind_real-h1)*aeroc%reff(j1,itype ) + h1*aeroc%reff(j2,itype )
     ENDIF

   END FUNCTION gocart_aerosol_size


   FUNCTION upper2lower(str) RESULT(string)

     IMPLICIT NONE

     CHARACTER(*), INTENT(in) :: str
     CHARACTER(LEN(str))      :: string

     INTEGER :: ic, i

     CHARACTER(26), PARAMETER :: upper = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
     CHARACTER(26), PARAMETER :: lower = 'abcdefghijklmnopqrstuvwxyz'

!   lowcase each letter if it is lowecase
     string = str
     DO i = 1, LEN_TRIM(str)
        ic = INDEX(upper, str(i:i))
        IF (ic > 0) string(i:i) = lower(ic:ic)
     END DO

   END FUNCTION upper2lower

   INTEGER FUNCTION getindex(names,usrname)
     IMPLICIT NONE
     CHARACTER(len=*),INTENT(in) :: names(:)
     CHARACTER(len=*),INTENT(in) :: usrname
     INTEGER i
     getindex=-1
     DO i=1,SIZE(names)
        IF(cmp_strings(usrname, names(i))) THEN
           getindex=i
           EXIT
        ENDIF
     ENDDO
   END FUNCTION getindex

!from fv3

   SUBROUTINE qsmith_atm(atm,rh)

     TYPE(CRTM_atmosphere_type), INTENT(in) :: atm(:)
     REAL(kind_real), INTENT(out),DIMENSION(:,:):: rh

     REAL, ALLOCATABLE :: table(:),des(:)

     REAL es, qs, q
     REAL ap1, eps10
     REAL Tmin
     INTEGER i, k, it, n_layers, n_profiles

     n_layers=SIZE(rh,1)
     n_profiles=SIZE(rh,2)

     Tmin = tice-160.
     eps10  = 10.*esl

     IF( .NOT. ALLOCATED(table) ) CALL qsmith_init(table,des)

     DO k=1,n_layers
        DO i=1,n_profiles
           ap1 = 10.*DIM(atm(i)%Temperature(k), Tmin) + 1.
           ap1 = MIN(2621., ap1)
           it = ap1
           es = table(it) + (ap1-it)*des(it)
           q=atm(i)%Absorber(k,1)*1.e-3/(1.+atm(i)%Absorber(k,1)*1.e-3)
           qs = esl*es*(1.+zvir*q)/(atm(i)%Pressure(k)*100.)
           rh(k,i) = q/qs
        ENDDO
     ENDDO

   END SUBROUTINE qsmith_atm

   SUBROUTINE qsmith_profiles(t,sphum,p,rh)

     REAL(kind_real), DIMENSION(:,:), INTENT(in) :: t,sphum,p
     REAL(kind_real), DIMENSION(:,:), INTENT(out) :: rh

     REAL, ALLOCATABLE :: table(:),des(:)

     REAL es, qs, q
     REAL ap1, eps10
     REAL Tmin
     INTEGER i, k, it, n_layers, n_profiles

     n_layers=SIZE(t,1)
     n_profiles=SIZE(t,2)

     Tmin = tice-160.
     eps10  = 10.*esl

     IF ( .NOT. ALLOCATED(table) ) CALL qsmith_init(table,des)

     DO k=1,n_layers
        DO i=1,n_profiles
           ap1 = 10.*DIM(t(k,i), Tmin) + 1.
           ap1 = MIN(2621., ap1)
           it = ap1
           es = table(it) + (ap1-it)*des(it)
           q=sphum(k,i)
           qs = esl*es*(1.+zvir*q)/p(k,i)
           rh(k,i) = q/qs
        ENDDO
     ENDDO

   END SUBROUTINE qsmith_profiles

   SUBROUTINE qsmith_init(table,des)

     REAL, ALLOCATABLE, INTENT(out) :: table(:),des(:)
     INTEGER, PARAMETER:: length=2621
     INTEGER i

     IF( .NOT. ALLOCATED(table) ) THEN
!                            Generate es table (dT = 0.1 deg. C)

        ALLOCATE ( table(length) )
        ALLOCATE (  des (length) )

        CALL qs_table(length, table)

        DO i=1,length-1
           des(i) = table(i+1) - table(i)
        ENDDO
        des(length) = des(length-1)
     ENDIF

   END SUBROUTINE qsmith_init

   SUBROUTINE qs_table(n,table)
     INTEGER, INTENT(in):: n
     REAL table (n)
     REAL :: dt=0.1
     REAL esbasw, tbasw, tbasi, Tmin, tem, aa, b, c, d, e
     INTEGER i
! Constants
     esbasw = 1013246.0
     tbasw =   373.16
     tbasi =   273.16
     Tmin = tbasi - 160.
!  Compute es over water
!  see smithsonian meteorological tables page 350.
     DO  i=1,n
        tem = Tmin+dt*REAL(i-1)
        aa  = -7.90298*(tbasw/tem-1)
        b   =  5.02808*alog10(tbasw/tem)
        c   = -1.3816e-07*(10**((1-tem/tbasw)*11.344)-1)
        d   =  8.1328e-03*(10**((tbasw/tem-1)*(-3.49149))-1)
        e   =  alog10(esbasw)
        table(i)  = 0.1*10**(aa+b+c+d+e)
     ENDDO

   END SUBROUTINE qs_table

END MODULE ufo_crtm_utils_mod
