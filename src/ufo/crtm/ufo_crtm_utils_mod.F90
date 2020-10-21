! (C) Copyright 2018 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

!> Fortran module to provide code shared between nonlinear and tlm/adm radiance calculations

MODULE ufo_crtm_utils_mod

use fckit_configuration_module, only: fckit_configuration
use iso_c_binding
use kinds

use crtm_module

use ufo_vars_mod
use ufo_geovals_mod, only: ufo_geovals, ufo_geoval, ufo_geovals_get_var
use obsspace_mod

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
public calculate_aero_layer_factor
PUBLIC :: qsmith
PUBLIC :: upper2lower

PUBLIC :: grav,rv_rd,&
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
 character(len=255) :: COEFFICIENT_PATH
 character(len=255) :: &
    IRwaterCoeff_File, IRlandCoeff_File, IRsnowCoeff_File, IRiceCoeff_File, &
    VISwaterCoeff_File, VISlandCoeff_File, VISsnowCoeff_File, VISiceCoeff_File, &
    MWwaterCoeff_File
 integer, allocatable :: Land_WSI(:)
 real(kind_real) :: Cloud_Fraction = -1.0_kind_real
 integer :: inspect
 character(len=MAXVARLEN) :: aerosol_option
 character(len=255) :: salinity_option
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
 character(len=*), parameter :: &
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

 character(len=MAXVARLEN), parameter :: &
      UFO_Clouds(N_VALID_CLOUD_CATEGORIES,2) = &
         reshape( &
            [ var_clw,    var_cli,    var_clr,    var_cls,    var_clg,    var_clh, &
              var_clwefr, var_cliefr, var_clrefr, var_clsefr, var_clgefr, var_clhefr ] &
            , [N_VALID_CLOUD_CATEGORIES,2] )

 ! copy of CLOUD_CATEGORY_NAME defined in CRTM_Cloud_Define
 character(len=*), parameter :: &
      CRTM_Clouds(N_VALID_CLOUD_CATEGORIES) = &
         CLOUD_CATEGORY_NAME(1:N_VALID_CLOUD_CATEGORIES)
 integer, parameter :: &
      CRTM_Cloud_Id(N_VALID_CLOUD_CATEGORIES) = &
         [   WATER_CLOUD, &
               ICE_CLOUD, &
              RAIN_CLOUD, &
              SNOW_CLOUD, &
           GRAUPEL_CLOUD, &
              HAIL_CLOUD  ]

! Surface Variables

 character(len=MAXVARLEN), parameter :: &
      UFO_Surfaces(4) = &
         [ var_sfc_wtmp, var_sfc_wspeed,var_sfc_wdir, var_sfc_sss]

 character(len=MAXVARLEN), parameter :: & 
      CRTM_Surfaces(4) = &
         [  character(len=MAXVARLEN):: 'Water_Temperature', 'Wind_Speed', 'Wind_Direction', 'Salinity' ]

contains

! ------------------------------------------------------------------------------

subroutine crtm_conf_setup(conf, f_confOpts, f_confOper)

implicit none
type(crtm_conf),            intent(inout) :: conf
type(fckit_configuration),  intent(in)    :: f_confOpts
type(fckit_configuration),  intent(in)    :: f_confOper

character(*), PARAMETER :: routine_name = 'crtm_conf_setup'
character(len=255) :: IRwaterCoeff, VISwaterCoeff, &
                      IRVISlandCoeff, IRVISsnowCoeff, IRVISiceCoeff, &
                      MWwaterCoeff
integer :: jspec, ivar
character(len=max_string) :: message
character(len=:), allocatable :: str
character(len=:), allocatable :: str_array(:)

CHARACTER(len=MAXVARLEN), ALLOCATABLE :: var_aerosols(:)


 !Some config needs to come from user
 !-----------------------------------

 !Number of sensors, each call to CRTM will be for a single sensor
 !type (zenith/scan angle will be different)
 conf%n_Sensors = 1

 ! absorbers, clouds and aerosols (should match what model will provide)

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
     CALL Display_Message(ROUTINE_NAME, TRIM(message), WARNING )
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
    case default
       allocate(conf%Land_WSI(1))
       conf%Land_WSI(1) = -1
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

subroutine ufo_crtm_skip_profiles(n_Profiles,n_Channels,channels,obss,Skip_Profiles)
! Profiles are skipped when the ObsValue of all channels is missing.
! TODO: Use complete QC information
! It would be more comprehensive to use EffectiveQC or EffectiveError. That
! would require those ObsSpace values to be initialized before calls to
! this subroutine within ufo_radiancecrtm_simobs+ufo_radiancecrtm_tlad_settraj.
use missing_values_mod

implicit none
integer,              intent(in)    :: n_Profiles, n_Channels
type(c_ptr), value,   intent(in)    :: obss
integer(c_int),       intent(in)    :: channels(:)
logical,              intent(inout) :: Skip_Profiles(:)

integer :: jprofile, jchannel
character(len=MAXVARLEN) :: varname
real(kind_real)  :: ObsVal(n_Profiles,n_Channels)
!real(kind_real) :: EffObsErr(n_Profiles,n_Channels)
!integer         :: EffQC(n_Profiles,n_Channels)

real(c_double) :: missing

 ! Set missing value
 missing = missing_value(missing)

 ObsVal = missing
! EffObsErr = missing
! EffQC = 0

 do jchannel = 1, n_Channels
   call get_var_name(channels(jchannel),varname)
   call obsspace_get_db(obss, "ObsValue", varname, ObsVal(:,jchannel))
!   call obsspace_get_db(obss, "EffectiveError", varname, EffObsErr(:,jchannel))
!   call obsspace_get_db(obss, "EffectiveQC{iter}", varname, EffQC(:,jchannel))
 enddo

 !Loop over all n_Profiles, i.e. number of locations
 do jprofile = 1, n_Profiles
   Skip_Profiles(jprofile) = all(ObsVal(jprofile,:) == missing)
!                       .OR. all(EffObsErr(jprofile,:) == missing) &
!                       .OR. all(EffQC(jprofile,:) /= 0)
 end do

end subroutine ufo_crtm_skip_profiles

! ------------------------------------------------------------------------------

SUBROUTINE Load_Atm_Data(n_Profiles, n_Layers, geovals, atm, conf)

implicit none

integer, intent(in) :: n_Profiles, n_Layers
type(ufo_geovals), intent(in) :: geovals
type(CRTM_Atmosphere_type), intent(inout) :: atm(:)
type(crtm_conf) :: conf

! Local variables
integer :: k1, jspec
type(ufo_geoval), pointer :: geoval
character(max_string) :: err_msg

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

  do jspec = 1, conf%n_Absorbers
    ! O3 Absorber has special treatment for Aerosols
    if ( trim(conf%Absorbers(jspec)) == trim(var_oz) .AND. &
      ufo_vars_getindex(geovals%variables, var_oz) < 0 .AND. &
      TRIM(conf%aerosol_option) /= "" ) then
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
  if (conf%n_Clouds > 0) then
    if ( ufo_vars_getindex(geovals%variables, var_cldfrac) > 0 ) then
      CALL ufo_geovals_get_var(geovals, var_cldfrac, geoval)
      do k1 = 1, n_Profiles
        atm(k1)%Cloud_Fraction(:) = geoval%vals(:, k1)
      end do
    else
      do k1 = 1, n_Profiles
        atm(k1)%Cloud_Fraction(:) = conf%Cloud_Fraction
      end do
    end if
  end if

end subroutine Load_Atm_Data

! ------------------------------------------------------------------------------

subroutine Load_Sfc_Data(n_Profiles, n_Channels, channels, geovals, sfc, chinfo, obss, conf)

implicit none

integer,                     intent(in)    :: n_Profiles, n_Channels
type(ufo_geovals),           intent(in)    :: geovals
type(CRTM_Surface_type),     intent(inout) :: sfc(:)
type(CRTM_ChannelInfo_type), intent(in)    :: chinfo(:)
type(c_ptr), value,          intent(in)    :: obss
integer(c_int),              intent(in)    :: channels(:)
type(crtm_conf),             intent(in)    :: conf

type(ufo_geoval), pointer :: geoval
integer :: k1, n1
integer :: iLand

! Surface type definitions for default SfcOptics definitions
! for IR and VIS, this is the NPOESS reflectivities.
integer, parameter :: TUNDRA_SURFACE_TYPE         = 10  ! NPOESS Land surface type for IR/VIS Land SfcOptics
integer, parameter :: SCRUB_SURFACE_TYPE          =  7  ! NPOESS Land surface type for IR/VIS Land SfcOptics
integer, parameter :: COARSE_SOIL_TYPE            =  1  ! Soil type                for MW land SfcOptics
integer, parameter :: GROUNDCOVER_VEGETATION_TYPE =  7  ! Vegetation type          for MW Land SfcOptics
integer, parameter :: BARE_SOIL_VEGETATION_TYPE   = 11  ! Vegetation type          for MW Land SfcOptics
integer, parameter :: SEA_WATER_TYPE              =  1  ! Water type               for all SfcOptics
integer, parameter :: FRESH_SNOW_TYPE             =  2  ! NPOESS Snow type         for IR/VIS SfcOptics
integer, parameter :: FRESH_ICE_TYPE              =  1  ! NPOESS Ice type          for IR/VIS SfcOptics

character(len=MAXVARLEN) :: varname

real(kind_real), allocatable :: ObsTb(:,:)

  allocate(ObsTb(n_profiles, n_channels))
  ObsTb = 0.0_kind_real

  do n1 = 1, n_Channels
    call get_var_name(channels(n1), varname)
    call obsspace_get_db(obss, "ObsValue", varname, ObsTb(:, n1))
  enddo

  do k1 = 1, n_Profiles
    !Pass sensor information
    sfc(k1)%sensordata%sensor_id        = chinfo(1)%sensor_id
    sfc(k1)%sensordata%wmo_sensor_id    = chinfo(1)%wmo_sensor_id
    sfc(k1)%sensordata%wmo_satellite_id = chinfo(1)%wmo_satellite_id
    sfc(k1)%sensordata%sensor_channel   = channels

    !Pass observation value
    do n1 = 1, n_channels
      sfc(k1)%sensordata%tb(n1) = ObsTb(k1, n1)
    enddo

    !Water_type
    sfc(k1)%Water_Type = SEA_WATER_TYPE    !** NOTE: need to check how to determine fresh vs sea water types (salinity???)
  end do
  deallocate(ObsTb)

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
  call ufo_geovals_get_var(geovals, var_sfc_landtyp, geoval)
  do k1 = 1, n_Profiles
    iLand = int(geoval%vals(1, k1))
    if (.not.any(iLand == conf%Land_WSI)) then
      sfc(k1)%Land_Type = iLand
    end if
  end do

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
  if (TRIM(conf%salinity_option) == "on") THEN
    call ufo_geovals_get_var(geovals, var_sfc_sss, geoval)
    do k1 = 1, n_Profiles
      sfc(k1)%Salinity = geoval%vals(1, k1)
    end do
  end if

end subroutine Load_Sfc_Data

! ------------------------------------------------------------------------------

subroutine Load_Geom_Data(obss,geo)

implicit none
type(c_ptr), value,       intent(in)    :: obss
type(CRTM_Geometry_type), intent(inout) :: geo(:)
real(kind_real), allocatable :: TmpVar(:)
integer :: nlocs

 nlocs = obsspace_get_nlocs(obss)
 allocate(TmpVar(nlocs))

 call obsspace_get_db(obss, "MetaData", "sensor_zenith_angle", TmpVar)
 geo(:)%Sensor_Zenith_Angle = abs(TmpVar(:)) ! needs to be absolute value

 call obsspace_get_db(obss, "MetaData", "solar_zenith_angle", TmpVar)
 geo(:)%Source_Zenith_Angle = TmpVar(:)

 call obsspace_get_db(obss, "MetaData", "sensor_azimuth_angle", TmpVar)
 geo(:)%Sensor_Azimuth_Angle = TmpVar(:)

 call obsspace_get_db(obss, "MetaData", "solar_azimuth_angle", TmpVar)
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

 call obsspace_get_db(obss, "MetaData", "scan_position", TmpVar)
 geo(:)%Ifov = TmpVar(:)

 call obsspace_get_db(obss, "MetaData", "sensor_view_angle", TmpVar) !The Sensor_Scan_Angle is optional
 geo(:)%Sensor_Scan_Angle = TmpVar(:)

 where (abs(geo(:)%Sensor_Scan_Angle) > 80.0_kind_real) &
    geo(:)%Sensor_Scan_Angle = 0.0_kind_real

 deallocate(TmpVar)

end subroutine Load_Geom_Data

! ------------------------------------------------------------------------------

subroutine get_var_name(n,varname)

integer, intent(in) :: n
character(len=*), intent(out) :: varname

character(len=6) :: chan

 write(chan, '(I0)') n
 varname = 'brightness_temperature_' // trim(chan)

end subroutine get_var_name

! -----------------------------------------------------------------------------


SUBROUTINE load_aerosol_data(n_profiles,n_layers,geovals,&
     &aerosol_option,atm)

    USE CRTM_aerosolcoeff, ONLY: aeroc

    INTEGER, INTENT(in) :: n_profiles,n_layers
    TYPE(ufo_geovals), INTENT(in) :: geovals
    TYPE(CRTM_atmosphere_type), INTENT(inout) :: atm(:)

    CHARACTER(*) :: aerosol_option
    CHARACTER(max_string) :: message
    CHARACTER(len=MAXVARLEN) :: varname

    TYPE(ufo_geoval), POINTER :: geoval

    CHARACTER(*), PARAMETER :: routine_name = 'Load_Aerosol_Data'

    REAL(kind_real), DIMENSION(n_layers,n_profiles) :: rh
    INTEGER :: ivar

    IF (TRIM(aerosol_option) == "aerosols_gocart_default") THEN
       varname=var_rh
       CALL ufo_geovals_get_var(geovals, varname, geoval)
       rh(1:n_layers,1:n_profiles)=geoval%vals(1:n_layers,1:n_profiles)
       WHERE (rh > 1_kind_real) rh=1_kind_real
       CALL assign_gocart_default
    ELSEIF (TRIM(aerosol_option) == "aerosols_gocart_merra_2") THEN
       varname=var_rh
       CALL ufo_geovals_get_var(geovals, varname, geoval)
       rh(1:n_layers,1:n_profiles)=geoval%vals(1:n_layers,1:n_profiles)
       WHERE (rh > 1_kind_real) rh=1_kind_real
       CALL assign_gocart_merra_2
    ELSEIF (TRIM(aerosol_option) == "aerosols_other") THEN
       CALL assign_other
    ELSE
       message = 'this aerosol not implemented - check later'
       CALL Display_Message( aerosol_option, message, FAILURE )
       STOP
    ENDIF

  CONTAINS 

    SUBROUTINE assign_gocart_default

      INTEGER, PARAMETER :: ndust_bins=5, nseas_bins=4
      REAL(kind_real), DIMENSION(ndust_bins), PARAMETER  :: dust_radii=[&
           &0.55_kind_real,1.4_kind_real,2.4_kind_real,4.5_kind_real,8.0_kind_real]
      
      INTEGER, DIMENSION(nseas_bins), PARAMETER  :: seas_types=[&
           SEASALT_SSAM_AEROSOL,SEASALT_SSCM1_AEROSOL,SEASALT_SSCM2_AEROSOL,SEASALT_SSCM3_AEROSOL]

      REAL(kind_real), DIMENSION(n_layers) :: layer_factors
      
      INTEGER :: i,k,m

      CHARACTER(len=MAXVARLEN) :: varname

      DO m=1,n_profiles

         CALL calculate_aero_layer_factor(atm(m),layer_factors)
 
         DO i=1,n_aerosols_gocart_default

            varname=var_aerosols_gocart_default(i)
            CALL ufo_geovals_get_var(geovals,varname, geoval)

            atm(m)%aerosol(i)%Concentration(1:n_layers)=&
                 &MAX(geoval%vals(:,m)*layer_factors,aerosol_concentration_minvalue_layer)

            SELECT CASE (TRIM(varname))
            CASE (var_sulfate)
               atm(m)%aerosol(i)%type  = SULFATE_AEROSOL
               DO k=1,n_layers
                  atm(m)%aerosol(i)%effective_radius(k)=&
                       &gocart_aerosol_size(atm(m)%aerosol(i)%type, &
                       &rh(k,m))
               ENDDO

            CASE (var_bcphobic)
               atm(m)%aerosol(i)%type  = BLACK_CARBON_AEROSOL
               atm(m)%aerosol(i)%effective_radius(:)=&
                    &AeroC%Reff(1,atm(m)%aerosol(i)%type)
            CASE (var_bcphilic)
               atm(m)%aerosol(i)%type  = BLACK_CARBON_AEROSOL
               DO k=1,n_layers
                  atm(m)%aerosol(i)%effective_radius(k)=&
                       &gocart_aerosol_size(atm(m)%aerosol(i)%type, &
                       &rh(k,m))
               ENDDO

            CASE (var_ocphobic)
               atm(m)%aerosol(i)%type  = ORGANIC_CARBON_AEROSOL
               atm(m)%aerosol(i)%effective_radius(:)=&
                    &AeroC%Reff(1,atm(m)%aerosol(i)%type)
            CASE (var_ocphilic)
               atm(m)%aerosol(i)%type  = ORGANIC_CARBON_AEROSOL
               DO k=1,n_layers
                  atm(m)%aerosol(i)%effective_radius(k)=&
                       &gocart_aerosol_size(atm(m)%aerosol(i)%type, &
                       &rh(k,m))
               ENDDO

            CASE (var_du001)
               atm(m)%aerosol(i)%type  = DUST_AEROSOL
               atm(m)%aerosol(i)%effective_radius(:)=dust_radii(1)
            CASE (var_du002)
               atm(m)%aerosol(i)%type  = DUST_AEROSOL
               atm(m)%aerosol(i)%effective_radius(:)=dust_radii(2)
            CASE (var_du003)
               atm(m)%aerosol(i)%type  = DUST_AEROSOL
               atm(m)%aerosol(i)%effective_radius(:)=dust_radii(3)
            CASE (var_du004)
               atm(m)%aerosol(i)%type  = DUST_AEROSOL
               atm(m)%aerosol(i)%effective_radius(:)=dust_radii(4)
            CASE (var_du005)
               atm(m)%aerosol(i)%type  = DUST_AEROSOL
               atm(m)%aerosol(i)%effective_radius(:)=dust_radii(5)

            CASE (var_ss001)
               atm(m)%aerosol(i)%type  = seas_types(1)
               DO k=1,n_layers
                  atm(m)%aerosol(i)%effective_radius(k)=&
                       &gocart_aerosol_size(atm(m)%aerosol(i)%type, &
                       &rh(k,m))
               ENDDO
            CASE (var_ss002)
               atm(m)%aerosol(i)%type  = seas_types(2)
               DO k=1,n_layers
                  atm(m)%aerosol(i)%effective_radius(k)=&
                       &gocart_aerosol_size(atm(m)%aerosol(i)%type, &
                       &rh(k,m))
               ENDDO
            CASE (var_ss003)
               atm(m)%aerosol(i)%type  = seas_types(3)
               DO k=1,n_layers
                  atm(m)%aerosol(i)%effective_radius(k)=&
                       &gocart_aerosol_size(atm(m)%aerosol(i)%type, &
                       &rh(k,m))
               ENDDO
            CASE (var_ss004)
               atm(m)%aerosol(i)%type  = seas_types(4)
               DO k=1,n_layers
                  atm(m)%aerosol(i)%effective_radius(k)=&
                       &gocart_aerosol_size(atm(m)%aerosol(i)%type, &
                       &rh(k,m))
               ENDDO
            END SELECT

         ENDDO

      ENDDO


    END SUBROUTINE assign_gocart_default

    SUBROUTINE assign_gocart_merra_2

      message = 'this aerosol not implemented in the CRTM - check later'
      CALL Display_Message( aerosol_option, message, FAILURE )
      STOP

    END SUBROUTINE assign_gocart_merra_2

    SUBROUTINE assign_other

      message = 'this aerosol not implemented - check later'
      CALL Display_Message( aerosol_option, message, FAILURE )
      STOP

    END SUBROUTINE assign_other
    
  END SUBROUTINE load_aerosol_data

  SUBROUTINE assign_aerosol_names(aerosol_option,var_aerosols)
    
    CHARACTER(*), INTENT(in) :: aerosol_option
    CHARACTER(len=MAXVARLEN), ALLOCATABLE, INTENT(out) :: var_aerosols(:)

    CHARACTER(max_string) :: err_msg

    IF (aerosol_option == "aerosols_gocart_default") THEN
       ALLOCATE(var_aerosols(SIZE(var_aerosols_gocart_default)))
       var_aerosols=var_aerosols_gocart_default
    ELSEIF (aerosol_option == "aerosols_gocart_merra_2") THEN
       ALLOCATE(var_aerosols(SIZE(var_aerosols_gocart_merra_2)))
       var_aerosols=var_aerosols_gocart_merra_2
    ELSEIF (aerosol_option == "var_aerosols_other") THEN
       ALLOCATE(var_aerosols(SIZE(var_aerosols_other)))
       var_aerosols=var_aerosols_other
    ELSE
       WRITE(err_msg,*) 'assign_aerosol_names: aerosol_option not implemented '//TRIM(aerosol_option)
        call abor1_ftn(err_msg)
     END IF

   END SUBROUTINE assign_aerosol_names

   SUBROUTINE calculate_aero_layer_factor_atm_profile(atm,layer_factors)

     TYPE(CRTM_atmosphere_type), INTENT(in) :: atm
     REAL(kind_real), INTENT(out) :: layer_factors(:)

     INTEGER :: k

     DO k=1,SIZE(layer_factors)
!correct for mixing ratio factor layer_factors 
!being calculated from dry pressure, cotton eq. (2.4)
!p_dry=p_total/(1+1.61*mixing_ratio)
        layer_factors(k)=1e-9_kind_real*(atm%Level_Pressure(k)-&
             &atm%Level_Pressure(k-1))*100_kind_real/grav/&
             &(1_kind_real+rv_rd*atm%Absorber(k,1)*1e-3_kind_real)
     ENDDO

   END SUBROUTINE calculate_aero_layer_factor_atm_profile

   SUBROUTINE calculate_aero_layer_factor_atm(atm,layer_factors)

     TYPE(CRTM_atmosphere_type), INTENT(in) :: atm(:)
     REAL(kind_real), INTENT(out) :: layer_factors(:,:)

     INTEGER :: k,m

     DO k=1,SIZE(layer_factors,1)
        DO m=1,SIZE(layer_factors,2)
!correct for mixing ratio factor layer_factors 
!being calculated from dry pressure, cotton eq. (2.4)
!p_dry=p_total/(1+1.61*mixing_ratio)
           layer_factors(k,m)=1e-9_kind_real*(atm(m)%Level_Pressure(k)-&
                &atm(m)%Level_Pressure(k-1))*100_kind_real/grav/&
                &(1_kind_real+rv_rd*atm(m)%Absorber(k,1)*1.e-3_kind_real)
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
     REAL(kind_real)    :: h1
     REAL(kind_real)    :: r_eff

     j2 = 0
     j1 = 1
     IF ( rh <= aeroc%rh(1) ) THEN
        j1 = 1
     ELSE IF ( rh >= aeroc%rh(aeroc%n_rh) ) THEN
        j1 = aeroc%n_rh
     ELSE
        DO m = 1, aeroc%n_rh-1
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
        r_eff = (1_kind_real-h1)*aeroc%reff(j1,itype ) + h1*aeroc%reff(j2,itype )
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
        IF(TRIM(usrname)==TRIM(names(i))) THEN
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

