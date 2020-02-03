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

REAL(kind_real), PARAMETER :: &
     &rd = 2.8704e+2_kind_real,&
     &rv = 4.6150e+2_kind_real,&
     &eps_p1 = one+rd/rv,&
     &grav = 9.81_kind_real,&
     &aerosol_concentration_minvalue=1.e-16_kind_real,&
     &aerosol_concentration_minvalue_layer=tiny(rd),& 
     &ozone_default_value=1.e-3_kind_real ! in ppmv in crtm

integer, parameter, public :: max_string=800


!Type for general config
type crtm_conf
 integer :: n_Sensors
 integer :: n_Absorbers
 integer :: n_Clouds
 integer :: n_Aerosols
 character(len=MAXVARLEN), allocatable :: Absorbers(:)
 integer, allocatable :: Absorber_Id(:)
 integer, allocatable :: Absorber_Units(:)
 character(len=MAXVARLEN), allocatable :: Clouds(:,:)
 integer, allocatable :: Cloud_Id(:)

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
 character(len=255) :: aerosol_option
 character(len=255) :: salinity_option
end type crtm_conf

INTERFACE calculate_aero_layer_factor

   MODULE PROCEDURE calculate_aero_layer_factor_atm_profile,&
        &calculate_aero_layer_factor_atm

END INTERFACE

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

contains

! ------------------------------------------------------------------------------

subroutine crtm_conf_setup(conf, f_confOpts, f_confOper)

implicit none
type(crtm_conf),            intent(inout) :: conf
type(fckit_configuration),  intent(in)    :: f_confOpts
type(fckit_configuration),  intent(in)    :: f_confOper

character(*), PARAMETER :: routine_name = 'crtm_conf_setup'
character(len=255) :: IRVISwaterCoeff, IRVISlandCoeff, IRVISsnowCoeff, IRVISiceCoeff, MWwaterCoeff
integer :: jspec, ivar
character(len=max_string) :: message
character(len=:), allocatable :: str
character(kind=c_char,len=MAXVARLEN), allocatable :: char_array(:)
integer(c_size_t),parameter :: csize = MAXVARLEN

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
   call f_confOper%get_or_die("Absorbers",csize,char_array)
   conf%Absorbers(1:conf%n_Absorbers) = char_array
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
   call f_confOper%get_or_die("Clouds",csize,char_array) 
   conf%Clouds(1:conf%n_Clouds,1) = char_array

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
    IF (TRIM(conf%aerosol_option) == "aerosols_gocart_default") THEN
       conf%n_Aerosols=n_aerosols_gocart_default
    ELSEIF (TRIM(conf%aerosol_option) == "aerosols_gocart_esrl") THEN
       conf%n_Aerosols=n_aerosols_gocart_esrl
    ELSEIF (TRIM(conf%aerosol_option) == "aerosols_other") THEN
       conf%n_Aerosols=n_aerosols_other
    ELSE
       conf%n_Aerosols=0
    ENDIF
 ELSE
    conf%n_Aerosols  = 0
    conf%aerosol_option = ""
 ENDIF

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
 IRVISwaterCoeff = "Nalli"
 if (f_confOpts%has("IRVISwaterCoeff")) then
    call f_confOpts%get_or_die("IRVISwaterCoeff",str)
    IRVISwaterCoeff = str
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
 conf%IRwaterCoeff_File = trim(IRVISwaterCoeff)//".IRwater.EmisCoeff.bin"
 conf%IRlandCoeff_File  = trim(IRVISlandCoeff)//".IRland.EmisCoeff.bin"
 conf%IRsnowCoeff_File  = trim(IRVISsnowCoeff)//".IRsnow.EmisCoeff.bin"
 conf%IRiceCoeff_File   = trim(IRVISiceCoeff)//".IRice.EmisCoeff.bin"

 !VIS emissivity coeff files
 conf%VISwaterCoeff_File = trim(IRVISwaterCoeff)//".VISwater.EmisCoeff.bin"
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

SUBROUTINE Load_Atm_Data(N_PROFILES,N_LAYERS,geovals,atm,conf)

implicit none
integer, intent(in) :: N_PROFILES, N_LAYERS
type(ufo_geovals), intent(in) :: geovals
type(CRTM_Atmosphere_type), intent(inout) :: atm(:)
type(crtm_conf) :: conf

! Local variables
integer :: k1,jspec
type(ufo_geoval), pointer :: geoval
character(max_string) :: err_msg


 ! Populate the atmosphere structures for CRTM (atm(k1), for the k1-th profile)
 ! ----------------------------------------------------------------------------
 do k1 = 1,N_PROFILES
    call ufo_geovals_get_var(geovals, var_ts, geoval)

    ! Check model levels is consistent in geovals & crtm
    if (k1 == 1) then
      if (geoval%nval /= n_Layers) then
        write(err_msg,*) 'Load_Atm_Data error: layers inconsistent!'
        call abor1_ftn(err_msg)
      endif
    endif

    atm(k1)%Temperature(1:N_LAYERS) = geoval%vals(:,k1)

    call ufo_geovals_get_var(geovals, var_prs, geoval)
    atm(k1)%Pressure(1:N_LAYERS) = geoval%vals(:,k1) * 0.01  ! to hPa
    call ufo_geovals_get_var(geovals, var_prsi, geoval)
    atm(k1)%Level_Pressure(:) = geoval%vals(:,k1) * 0.01     ! to hPa
    atm(k1)%Climatology         = US_STANDARD_ATMOSPHERE

    do jspec = 1, conf%n_Absorbers
       ! O3 Absorber has special treatment for Aerosols
       if ( trim(conf%Absorbers(jspec)) == trim(var_oz) .AND. &
            ufo_vars_getindex(geovals%variables, var_oz) < 0 .AND. &
            TRIM(conf%aerosol_option) /= "" ) then
          atm(k1)%Absorber(1:N_LAYERS,jspec) = ozone_default_value
       else
          CALL ufo_geovals_get_var(geovals, conf%Absorbers(jspec), geoval)
          atm(k1)%Absorber(1:N_LAYERS,jspec) = geoval%vals(:,k1)
       end if
       atm(k1)%Absorber_Id(jspec)    = conf%Absorber_Id(jspec)
       atm(k1)%Absorber_Units(jspec) = conf%Absorber_Units(jspec)
    end do

    do jspec = 1, conf%n_Clouds
       CALL ufo_geovals_get_var(geovals, conf%Clouds(jspec,1), geoval)
       atm(k1)%Cloud(jspec)%Water_Content = geoval%vals(:,k1)
       CALL ufo_geovals_get_var(geovals, conf%Clouds(jspec,2), geoval)
       atm(k1)%Cloud(jspec)%Effective_Radius = geoval%vals(:,k1)
       atm(k1)%Cloud(jspec)%Type = conf%Cloud_Id(jspec)
    end do
 end do

 ! When n_Clouds>0, Cloud_Fraction must either be provided as geoval or in conf
 if (conf%n_Clouds > 0) then
    if ( ufo_vars_getindex(geovals%variables, var_cldfrac) > 0 ) then
       CALL ufo_geovals_get_var(geovals, var_cldfrac, geoval)
       do k1 = 1,N_PROFILES
          atm(k1)%Cloud_Fraction(:) = geoval%vals(:,k1)
       end do
    else
       do k1 = 1,N_PROFILES
          atm(k1)%Cloud_Fraction(:) = conf%Cloud_Fraction
       end do
    end if
 end if

 end subroutine Load_Atm_Data

! ------------------------------------------------------------------------------

subroutine Load_Sfc_Data(n_Profiles,n_Channels,channels,geovals,sfc,chinfo,obss,conf)

implicit none
integer,                     intent(in)    :: n_Profiles, n_Channels
type(ufo_geovals),           intent(in)    :: geovals
type(CRTM_Surface_type),     intent(inout) :: sfc(:)
type(CRTM_ChannelInfo_type), intent(in)    :: chinfo(:)
type(c_ptr), value,          intent(in)    :: obss
integer(c_int),              intent(in)    :: channels(:)
type(crtm_conf),             intent(in)    :: conf

type(ufo_geoval), pointer :: geoval
integer  :: k1, n1

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

character(len=200) :: varname

real(kind_real), allocatable :: ObsTb(:,:)

 allocate(ObsTb(n_profiles, n_channels))
 ObsTb = 0.0_kind_real

 do n1 = 1,n_Channels
   call get_var_name(channels(n1),varname)
   call obsspace_get_db(obss, "ObsValue", varname, ObsTb(:,n1))
 enddo

 !Loop over all n_Profiles, i.e. number of locations
 do k1 = 1,N_PROFILES

   !Pass sensor information
   sfc(k1)%sensordata%sensor_id        = chinfo(1)%sensor_id
   sfc(k1)%sensordata%wmo_sensor_id    = chinfo(1)%wmo_sensor_id
   sfc(k1)%sensordata%wmo_satellite_id = chinfo(1)%wmo_satellite_id
   sfc(k1)%sensordata%sensor_channel   = channels

   !Pass observation value
   do n1 = 1, n_channels
     sfc(k1)%sensordata%tb(n1) = ObsTb(k1,n1)
   enddo

   !Water_type
   sfc(k1)%Water_Type         = SEA_WATER_TYPE    !** NOTE: need to check how to determine fresh vs sea water types (salinity???)

   !Wind_Speed
   call ufo_geovals_get_var(geovals, var_sfc_wspeed, geoval)
   sfc(k1)%Wind_Speed = geoval%vals(1,k1)

   !Wind_Direction
   call ufo_geovals_get_var(geovals, var_sfc_wdir, geoval)
   sfc(k1)%Wind_Direction = geoval%vals(1,k1)

   !Water_Coverage
   call ufo_geovals_get_var(geovals, var_sfc_wfrac, geoval)
   sfc(k1)%Water_Coverage = geoval%vals(1,k1)

   !Water_Temperature
   call ufo_geovals_get_var(geovals, var_sfc_wtmp, geoval)
   sfc(k1)%Water_Temperature = geoval%vals(1,k1)

   !Ice_Coverage
   call ufo_geovals_get_var(geovals, var_sfc_ifrac, geoval)
   sfc(k1)%Ice_Coverage = geoval%vals(1,k1)

   !Ice_Temperature
   call ufo_geovals_get_var(geovals, var_sfc_itmp, geoval)
   sfc(k1)%Ice_Temperature = geoval%vals(1,k1)

   !Snow_Coverage
   call ufo_geovals_get_var(geovals, var_sfc_sfrac, geoval)
   sfc(k1)%Snow_Coverage      = geoval%vals(1,k1)

   !Snow_Temperature
   call ufo_geovals_get_var(geovals, var_sfc_stmp, geoval)
   sfc(k1)%Snow_Temperature = geoval%vals(1,k1)

   !Snow_Depth
   call ufo_geovals_get_var(geovals, var_sfc_sdepth, geoval)
   sfc(k1)%Snow_Depth = geoval%vals(1,k1)

   !Land_Coverage
   call ufo_geovals_get_var(geovals, var_sfc_lfrac, geoval)
   sfc(k1)%Land_Coverage = geoval%vals(1,k1)

   !Land_Type
   ! + used to lookup land sfc emiss. for IR and VIS
   ! + land sfc emiss. undefined over water/snow/ice
   call ufo_geovals_get_var(geovals, var_sfc_landtyp, geoval)
   if (.not.any(int(geoval%vals(1,k1)) == conf%Land_WSI)) then
      sfc(k1)%Land_Type = int(geoval%vals(1,k1))
   end if

   !Land_Temperature
   call ufo_geovals_get_var(geovals, var_sfc_ltmp, geoval)
   sfc(k1)%Land_Temperature = geoval%vals(1,k1)

   !Lai
   call ufo_geovals_get_var(geovals, var_sfc_lai, geoval)
   sfc(k1)%Lai = geoval%vals(1,k1)

   !Vegetation_Fraction
   call ufo_geovals_get_var(geovals, var_sfc_vegfrac, geoval)
   sfc(k1)%Vegetation_Fraction = geoval%vals(1,k1)

   !Vegetation_Type
   call ufo_geovals_get_var(geovals, var_sfc_vegtyp, geoval)
   sfc(k1)%Vegetation_Type = int(geoval%vals(1,k1))

   !Soil_Type
   call ufo_geovals_get_var(geovals, var_sfc_soiltyp, geoval)
   sfc(k1)%Soil_Type = int(geoval%vals(1,k1))

   !Soil_Moisture_Content
   call ufo_geovals_get_var(geovals, var_sfc_soilm, geoval)
   sfc(k1)%Soil_Moisture_Content = geoval%vals(1,k1)

   !Soil_Temperature
   call ufo_geovals_get_var(geovals, var_sfc_soilt, geoval)
   sfc(k1)%Soil_Temperature = geoval%vals(1,k1)
   
   !Sea_Surface_Salinity
   if (TRIM(conf%salinity_option) == "on") THEN
      call ufo_geovals_get_var(geovals, var_sfc_sss, geoval)
      sfc(k1)%Salinity = geoval%vals(1,k1)
   end if

 end do

 deallocate(ObsTb)

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
       ivar = ufo_vars_getindex(geovals%variables, var_rh)
       IF (ivar < 0) THEN
          message='relative humidity missing as input - will be calculated from tables'
          CALL Display_Message(ROUTINE_NAME, TRIM(message), WARNING )
          CALL qsmith(atm,rh)
       ELSE
          CALL ufo_geovals_get_var(geovals, varname, geoval)
          rh(1:n_layers,1:n_profiles)=geoval%vals(1:n_layers,1:n_profiles)
       ENDIF
       WHERE (rh > 1_kind_real) rh=1_kind_real
       CALL assign_gocart_default
    ELSEIF (TRIM(aerosol_option) == "aerosols_gocart_esrl") THEN
       varname=var_rh
       ivar = ufo_vars_getindex(geovals%variables, var_rh)
       IF (ivar < 0) THEN
          message='relative humidity missing as input - will be calculated from tables'
          CALL Display_Message(ROUTINE_NAME, TRIM(message), WARNING )
          CALL qsmith(atm,rh)
       ELSE
          CALL ufo_geovals_get_var(geovals, varname, geoval)
          rh(1:n_layers,1:n_profiles)=geoval%vals(1:n_layers,1:n_profiles)
       ENDIF
       WHERE (rh > 1_kind_real) rh=1_kind_real
       CALL assign_gocart_esrl
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

      REAL(kind_real), DIMENSION(n_layers) :: ugkg_kgm2
      
      INTEGER :: i,k,m

      DO m=1,n_profiles

         CALL calculate_aero_layer_factor(atm(m),ugkg_kgm2)
 
         DO i=1,n_aerosols_gocart_default

            varname=var_aerosols_gocart_default(i)
            CALL ufo_geovals_get_var(geovals,varname, geoval)

            atm(m)%aerosol(i)%Concentration(1:n_layers)=&
                 &MAX(geoval%vals(:,m)*ugkg_kgm2,aerosol_concentration_minvalue_layer)

            SELECT CASE ( TRIM(varname))
            CASE ('sulf')
               atm(m)%aerosol(i)%type  = SULFATE_AEROSOL
               DO k=1,n_layers
                  atm(m)%aerosol(i)%effective_radius(k)=&
                       &gocart_aerosol_size(atm(m)%aerosol(i)%type, &
                       &rh(k,m))
               ENDDO

            CASE ('bc1')
               atm(m)%aerosol(i)%type  = BLACK_CARBON_AEROSOL
               atm(m)%aerosol(i)%effective_radius(:)=&
                    &AeroC%Reff(1,atm(m)%aerosol(i)%type)
            CASE ('bc2')
               atm(m)%aerosol(i)%type  = BLACK_CARBON_AEROSOL
               DO k=1,n_layers
                  atm(m)%aerosol(i)%effective_radius(k)=&
                       &gocart_aerosol_size(atm(m)%aerosol(i)%type, &
                       &rh(k,m))
               ENDDO

            CASE ('oc1')
               atm(m)%aerosol(i)%type  = ORGANIC_CARBON_AEROSOL
               atm(m)%aerosol(i)%effective_radius(:)=&
                    &AeroC%Reff(1,atm(m)%aerosol(i)%type)
            CASE ('oc2')
               atm(m)%aerosol(i)%type  = ORGANIC_CARBON_AEROSOL
               DO k=1,n_layers
                  atm(m)%aerosol(i)%effective_radius(k)=&
                       &gocart_aerosol_size(atm(m)%aerosol(i)%type, &
                       &rh(k,m))
               ENDDO

            CASE ('dust1')
               atm(m)%aerosol(i)%type  = DUST_AEROSOL
               atm(m)%aerosol(i)%effective_radius(:)=dust_radii(1)
            CASE ('dust2')
               atm(m)%aerosol(i)%type  = DUST_AEROSOL
               atm(m)%aerosol(i)%effective_radius(:)=dust_radii(2)
            CASE ('dust3')
               atm(m)%aerosol(i)%type  = DUST_AEROSOL
               atm(m)%aerosol(i)%effective_radius(:)=dust_radii(3)
            CASE ('dust4')
               atm(m)%aerosol(i)%type  = DUST_AEROSOL
               atm(m)%aerosol(i)%effective_radius(:)=dust_radii(4)
            CASE ('dust5')
               atm(m)%aerosol(i)%type  = DUST_AEROSOL
               atm(m)%aerosol(i)%effective_radius(:)=dust_radii(5)

            CASE ('seas1')
               atm(m)%aerosol(i)%type  = seas_types(1)
               DO k=1,n_layers
                  atm(m)%aerosol(i)%effective_radius(k)=&
                       &gocart_aerosol_size(atm(m)%aerosol(i)%type, &
                       &rh(k,m))
               ENDDO
            CASE ('seas2')
               atm(m)%aerosol(i)%type  = seas_types(2)
               DO k=1,n_layers
                  atm(m)%aerosol(i)%effective_radius(k)=&
                       &gocart_aerosol_size(atm(m)%aerosol(i)%type, &
                       &rh(k,m))
               ENDDO
            CASE ('seas3')
               atm(m)%aerosol(i)%type  = seas_types(3)
               DO k=1,n_layers
                  atm(m)%aerosol(i)%effective_radius(k)=&
                       &gocart_aerosol_size(atm(m)%aerosol(i)%type, &
                       &rh(k,m))
               ENDDO
            CASE ('seas4')
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

    SUBROUTINE assign_gocart_esrl

      INTEGER, PARAMETER :: ndust_bins=5, nseas_bins=4
      REAL(kind_real), DIMENSION(ndust_bins), PARAMETER  :: dust_radii=[&
           &0.55_kind_real,1.4_kind_real,2.4_kind_real,4.5_kind_real,8.0_kind_real]

      REAL(kind_real),PARAMETER  :: p25_radius=0.9_kind_real
!p25_radius <- (0.78*(dust_radii_esrl[1])^3+
!               0.22*(dust_radii_esrl[2])^3)^(1./3.)

      
      INTEGER, DIMENSION(nseas_bins), PARAMETER  :: seas_types=[&
           SEASALT_SSAM_AEROSOL,SEASALT_SSCM1_AEROSOL,SEASALT_SSCM2_AEROSOL,SEASALT_SSCM3_AEROSOL]
      
      REAL(kind_real), DIMENSION(n_layers) :: ugkg_kgm2
      
      INTEGER :: i,k,m
      
      message = 'this aerosol not implemented - check later'
      CALL Display_Message( aerosol_option, message, FAILURE )
      STOP

      DO m=1,n_profiles

         CALL calculate_aero_layer_factor(atm(m),ugkg_kgm2)

         DO i=1,n_aerosols_gocart_esrl
            varname=var_aerosols_gocart_esrl(i)
            CALL ufo_geovals_get_var(geovals,varname, geoval)

            atm(m)%aerosol(i)%Concentration(1:n_layers)=&
                 &MAX(geoval%vals(:,m)*ugkg_kgm2,aerosol_concentration_minvalue_layer)

            SELECT CASE ( TRIM(varname))
            CASE ('sulf')
               atm(m)%aerosol(i)%type  = SULFATE_AEROSOL

               DO k=1,n_layers
                  atm(m)%aerosol(i)%effective_radius(k)=&
                       &gocart_aerosol_size(atm(m)%aerosol(i)%type, &
                       &rh(k,m))
               ENDDO

            CASE ('bc1')
               atm(m)%aerosol(i)%type  = BLACK_CARBON_AEROSOL
               atm(m)%aerosol(i)%effective_radius(:)=&
                    &AeroC%Reff(1,atm(m)%aerosol(i)%type)
            CASE ('bc2')
               atm(m)%aerosol(i)%type  = BLACK_CARBON_AEROSOL
               DO k=1,n_layers
                  atm(m)%aerosol(i)%effective_radius(k)=&
                       &gocart_aerosol_size(atm(m)%aerosol(i)%type, &
                       &rh(k,m))
               ENDDO

            CASE ('oc1')
               atm(m)%aerosol(i)%type  = ORGANIC_CARBON_AEROSOL
               atm(m)%aerosol(i)%effective_radius(:)=&
                    &AeroC%Reff(1,atm(m)%aerosol(i)%type)
            CASE ('oc2')
               atm(m)%aerosol(i)%type  = ORGANIC_CARBON_AEROSOL
               DO k=1,n_layers
                  atm(m)%aerosol(i)%effective_radius(k)=&
                       &gocart_aerosol_size(atm(m)%aerosol(i)%type, &
                       &rh(k,m))
               ENDDO

            CASE ('dust1')
               atm(m)%aerosol(i)%type  = DUST_AEROSOL
               atm(m)%aerosol(i)%effective_radius(:)=dust_radii(1)
            CASE ('dust2')
               atm(m)%aerosol(i)%type  = DUST_AEROSOL
               atm(m)%aerosol(i)%effective_radius(:)=dust_radii(2)
            CASE ('dust3')
               atm(m)%aerosol(i)%type  = DUST_AEROSOL
               atm(m)%aerosol(i)%effective_radius(:)=dust_radii(3)
            CASE ('dust4')
               atm(m)%aerosol(i)%type  = DUST_AEROSOL
               atm(m)%aerosol(i)%effective_radius(:)=dust_radii(4)
            CASE ('dust5')
               atm(m)%aerosol(i)%type  = DUST_AEROSOL
               atm(m)%aerosol(i)%effective_radius(:)=dust_radii(5)

            CASE ('p25')
               atm(m)%aerosol(i)%type  = DUST_AEROSOL
               atm(m)%aerosol(i)%effective_radius(:)=p25_radius

            CASE ('seas1')
               atm(m)%aerosol(i)%type  = seas_types(1)
               DO k=1,n_layers
                  atm(m)%aerosol(i)%effective_radius(k)=&
                       &gocart_aerosol_size(atm(m)%aerosol(i)%type, &
                       &rh(k,m))
               ENDDO
            CASE ('seas2')
               atm(m)%aerosol(i)%type  = seas_types(2)
               DO k=1,n_layers
                  atm(m)%aerosol(i)%effective_radius(k)=&
                       &gocart_aerosol_size(atm(m)%aerosol(i)%type, &
                       &rh(k,m))
               ENDDO
            CASE ('seas3')
               atm(m)%aerosol(i)%type  = seas_types(3)
               DO k=1,n_layers
                  atm(m)%aerosol(i)%effective_radius(k)=&
                       &gocart_aerosol_size(atm(m)%aerosol(i)%type, &
                       &rh(k,m))
               ENDDO
            CASE ('seas4')
               atm(m)%aerosol(i)%type  = seas_types(4)
               DO k=1,n_layers
                  atm(m)%aerosol(i)%effective_radius(k)=&
                       &gocart_aerosol_size(atm(m)%aerosol(i)%type, &
                       &rh(k,m))
               ENDDO
            END SELECT

         ENDDO

      ENDDO

    END SUBROUTINE assign_gocart_esrl

    SUBROUTINE assign_other

      message = 'this aerosol not implemented - check later'
      CALL Display_Message( aerosol_option, message, FAILURE )
      STOP

    END SUBROUTINE assign_other
    
  END SUBROUTINE load_aerosol_data

  SUBROUTINE assign_aerosol_names(aerosol_option,var_aerosols)
    
    CHARACTER(*), INTENT(in) :: aerosol_option
    CHARACTER(*), INTENT(out) :: var_aerosols(:)

    CHARACTER(max_string) :: err_msg

    IF (aerosol_option == "aerosols_gocart_default") THEN
       var_aerosols=var_aerosols_gocart_default
    ELSEIF (aerosol_option == "aerosols_gocart_esrl") THEN
       var_aerosols=var_aerosols_gocart_esrl
    ELSEIF (aerosol_option == "var_aerosols_other") THEN
       var_aerosols=var_aerosols_other
    ELSE
       WRITE(err_msg,*) 'assign_aerosol_names: aerosol_option not implemented '//TRIM(aerosol_option)
        call abor1_ftn(err_msg)
     END IF

   END SUBROUTINE assign_aerosol_names
   
   SUBROUTINE calculate_aero_layer_factor_atm_profile(atm,ugkg_kgm2)

     TYPE(CRTM_atmosphere_type), INTENT(in) :: atm
     REAL(kind_real), INTENT(out) :: ugkg_kgm2(:)

     INTEGER :: k

!rh, ug2kg need to be from top to bottom    
!ug2kg && hPa2Pa
     DO k=1,SIZE(ugkg_kgm2)
!correct for mixing ratio factor ugkg_kgm2 
!being calculated from dry pressure, cotton eq. (2.4)
!p_dry=p_total/(1+1.61*mixing_ratio)
        ugkg_kgm2(k)=1.0e-9_kind_real*(atm%Level_Pressure(k)-&
             &atm%Level_Pressure(k-1))*100_kind_real/grav/&
             &(one+eps_p1*atm%Absorber(k,1)*1.e-3_kind_real)
     ENDDO

   END SUBROUTINE calculate_aero_layer_factor_atm_profile

   SUBROUTINE calculate_aero_layer_factor_atm(atm,ugkg_kgm2)

     TYPE(CRTM_atmosphere_type), INTENT(in) :: atm(:)
     REAL(kind_real), INTENT(out) :: ugkg_kgm2(:,:)

     INTEGER :: k,m

!rh, ug2kg need to be from top to bottom    
!ug2kg && hPa2Pa
     DO k=1,SIZE(ugkg_kgm2,1)
        DO m=1,SIZE(ugkg_kgm2,2)
!correct for mixing ratio factor ugkg_kgm2 
!being calculated from dry pressure, cotton eq. (2.4)
!p_dry=p_total/(1+1.61*mixing_ratio)
           ugkg_kgm2(k,m)=1.0e-9_kind_real*(atm(m)%Level_Pressure(k)-&
                &atm(m)%Level_Pressure(k-1))*100_kind_real/grav/&
                &(one+eps_p1*atm(m)%Absorber(k,1)*1.e-3_kind_real)
        ENDDO
     ENDDO
     
   END SUBROUTINE calculate_aero_layer_factor_atm

   FUNCTION gocart_aerosol_size( itype, rh ) & ! rh input in 0-1
        &RESULT(r_eff)   ! in micrometer

!@mzp: will be modified as in NASA's

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
        r_eff = (one-h1)*aeroc%reff(j1,itype ) + h1*aeroc%reff(j2,itype )
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

   SUBROUTINE qsmith(atm,rh)

     TYPE(CRTM_atmosphere_type), INTENT(in) :: atm(:)
     REAL(kind_real), INTENT(out),DIMENSION(:,:):: rh

! Local:
! input T in deg K; p (Pa)
     REAL, PARAMETER :: rdgas  = rd   !< gas constant for dry air [j/kg/deg]
     REAL, PARAMETER :: rvgas  = rv   !< gas constant for water vapor [J/kg/deg]
     
     REAL, PARAMETER :: esl = 0.621971831
     REAL, PARAMETER :: zvir =  rvgas/rdgas - 1. 
     REAL, PARAMETER :: tice = 273.16
     
     REAL, ALLOCATABLE :: table(:),des(:)

     REAL es, qs
     REAL ap1, eps10
     REAL Tmin
     INTEGER i, k, it, n_layers, n_profiles

     n_layers=SIZE(rh,1)
     n_profiles=SIZE(rh,2)

     Tmin = tice-160.
     eps10  = 10.*esl

     IF( .NOT. ALLOCATED(table) ) CALL  qsmith_init

     DO i=1,n_profiles
        DO k=1,n_layers
           ap1 = 10.*DIM(atm(i)%Temperature(k), Tmin) + 1.
           ap1 = MIN(2621., ap1)
           it = ap1
           es = table(it) + (ap1-it)*des(it)
           qs = esl*es*(1.+zvir*atm(i)%Absorber(k,1)*1.e-3)/(atm(i)%Pressure(k)*100.)
           rh(k,i) = (atm(i)%Absorber(k,1)*1.e-3)/qs
        ENDDO
     ENDDO

   CONTAINS

     SUBROUTINE qsmith_init

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
       REAL:: dt=0.1
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

   END SUBROUTINE qsmith
   
END MODULE ufo_crtm_utils_mod   

