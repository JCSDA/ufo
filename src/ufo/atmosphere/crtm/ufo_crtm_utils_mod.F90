! (C) Copyright 2018 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

!> Fortran module to provide code shared between nonlinear and tlm/adm radiance calculations

module ufo_crtm_utils_mod

use iso_c_binding
use config_mod
use kinds

use crtm_module

use ufo_vars_mod
use ufo_geovals_mod, only: ufo_geovals, ufo_geoval, ufo_geovals_get_var
use ufo_basis_mod, only: ufo_basis
use obsspace_mod

implicit none
private

public crtm_conf
public crtm_conf_setup
public crtm_conf_delete
public Load_Atm_Data
public Load_Sfc_Data
public Load_Geom_Data

PUBLIC Load_Aerosol_Data
PUBLIC check_fwd
public assign_aerosol_names
public calculate_aero_layer_factor

REAL(kind_real), PARAMETER :: &
     &rd = 2.8705e+2_kind_real,&
     &rv = 4.6150e+2_kind_real,&
     &eps_p1 = one+rd/rv,&
     &grav = 9.81_kind_real,&
     &aerosol_concentration_minvalue=1.e-16_kind_real,&
     &aerosol_concentration_minvalue_layer=tiny(rd),& 
     &ozone_default_value=1.e-3_kind_real ! in ppmv in crtm

INTEGER, PARAMETER, public :: min_crtm_n_absorbers = 2

integer, parameter, public :: max_string=800


!Type for general config
type crtm_conf
 integer :: n_Sensors
 integer :: n_Absorbers
 integer :: n_Clouds
 integer :: n_Aerosols
 character(len=255), allocatable :: SENSOR_ID(:)
 character(len=255) :: ENDIAN_TYPE
 character(len=255) :: COEFFICIENT_PATH
 integer :: inspect
 character(len=255) :: aerosol_option
end type crtm_conf

INTERFACE calculate_aero_layer_factor

   MODULE PROCEDURE calculate_aero_layer_factor_atm_profile,&
        &calculate_aero_layer_factor_atm

END INTERFACE

contains

! ------------------------------------------------------------------------------

subroutine crtm_conf_setup(conf, c_conf)

implicit none
type(crtm_conf), intent(inout) :: conf
type(c_ptr),    intent(in)    :: c_conf

 !Some config needs to come from user
 !-----------------------------------

 !Number of sensors, each call to CRTM will be for a single sensor
 !type (zenith/scan angle will be different)
 conf%n_Sensors = 1

 !Number of absorbers, clouds and aerosols (should match what model will provide)
 conf%n_Absorbers = config_get_int(c_conf,"n_Absorbers")
 conf%n_Clouds    = config_get_int(c_conf,"n_Clouds"   )

 IF (config_element_exists(c_conf,"n_Aerosols")) &
      &conf%n_Aerosols  = config_get_int(c_conf,"n_Aerosols" )

 IF (config_element_exists(c_conf,"AerosolOption")) THEN
    conf%aerosol_option = config_get_string(c_conf,LEN(conf%aerosol_option),"AerosolOption")
    conf%aerosol_option = upper2lower(conf%aerosol_option)
    IF (TRIM(conf%aerosol_option) == "aerosols_gocart_nasa") THEN
       conf%n_Aerosols=n_aerosols_gocart_nasa
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

 !Allocate SENSOR_ID
 allocate(conf%SENSOR_ID(conf%n_Sensors))

 !Get sensor ID from config
 conf%SENSOR_ID(conf%n_Sensors) = config_get_string(c_conf,len(conf%SENSOR_ID(conf%n_Sensors)),"Sensor_ID")

 !ENDIAN type
 conf%ENDIAN_TYPE = config_get_string(c_conf,len(conf%ENDIAN_TYPE),"EndianType")

 !Path to coefficient files
 conf%COEFFICIENT_PATH = config_get_string(c_conf,len(conf%COEFFICIENT_PATH),"CoefficientPath")

 conf%inspect = 0
 if (config_element_exists(c_conf,"InspectProfileNumber")) then
   conf%inspect = config_get_int(c_conf,"InspectProfileNumber")
 endif

end subroutine crtm_conf_setup

! -----------------------------------------------------------------------------

subroutine crtm_conf_delete(conf)

implicit none
type(crtm_conf), intent(inout) :: conf

 deallocate(conf%SENSOR_ID)

end subroutine crtm_conf_delete

! ------------------------------------------------------------------------------

SUBROUTINE Load_Atm_Data(N_PROFILES,N_LAYERS,geovals,atm,conf)

implicit none
integer, intent(in) :: N_PROFILES, N_LAYERS
type(ufo_geovals), intent(in) :: geovals
type(CRTM_Atmosphere_type), intent(inout) :: atm(:)
type(crtm_conf) :: conf

! Local variables
integer :: k1,ivar
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
    atm(k1)%Pressure(1:N_LAYERS) = geoval%vals(:,k1)
    call ufo_geovals_get_var(geovals, var_prsi, geoval)
    atm(k1)%Level_Pressure(:) = geoval%vals(:,k1)
    atm(k1)%Climatology         = US_STANDARD_ATMOSPHERE
    atm(k1)%Absorber_Id(1:1)    = (/ H2O_ID /)
    atm(k1)%Absorber_Units(1:1) = (/ MASS_MIXING_RATIO_UNITS /)
    call ufo_geovals_get_var(geovals, var_mixr, geoval)
    atm(k1)%Absorber(1:N_LAYERS,1)       = geoval%vals(:,k1)
    atm(k1)%Absorber_Id(2:2)    = (/ O3_ID /)
    atm(k1)%Absorber_Units(2:2) = (/ VOLUME_MIXING_RATIO_UNITS /)

    ivar = ufo_vars_getindex(geovals%variables, var_oz)
    IF (ivar < 0 .AND. TRIM(conf%aerosol_option) /= "") THEN 
       atm(k1)%Absorber(1:N_LAYERS,2)=ozone_default_value
    ELSE
       CALL ufo_geovals_get_var(geovals, var_oz, geoval)
       atm(k1)%Absorber(1:N_LAYERS,2)       = geoval%vals(:,k1)
    ENDIF

    IF (conf%n_Absorbers > min_crtm_n_absorbers) THEN

       atm(k1)%Absorber_Id(3:3)    = (/ CO2_ID /)
       atm(k1)%Absorber_Units(3:3) = (/ VOLUME_MIXING_RATIO_UNITS /)
       CALL ufo_geovals_get_var(geovals, var_co2, geoval)
       atm(k1)%Absorber(1:N_LAYERS,3)       = geoval%vals(:,k1)
       
    ENDIF

    IF ( conf%n_Clouds >= 1 ) THEN

       atm(k1)%Cloud(1)%Type = WATER_CLOUD
       CALL ufo_geovals_get_var(geovals, var_clw, geoval)
       atm(k1)%Cloud(1)%Water_Content = geoval%vals(:,k1)
       CALL ufo_geovals_get_var(geovals, var_clwefr, geoval)
       atm(k1)%Cloud(1)%Effective_Radius = geoval%vals(:,k1)

!** BTJ added 11/20/2018 for compatibility with CRTM REL 2.3.0+
!** need to map to cloud fraction geoval, if it exists.  For now assume
!** fully filled pixel. 
       atm(k1)%Cloud_Fraction = 1.0_fp  
       
    ENDIF

    IF ( conf%n_Clouds >= 2 ) THEN

       atm(k1)%Cloud(2)%Type = ICE_CLOUD
       CALL ufo_geovals_get_var(geovals, var_cli, geoval)
       atm(k1)%Cloud(2)%Water_Content = geoval%vals(:,k1)
       CALL ufo_geovals_get_var(geovals, var_cliefr, geoval)
       atm(k1)%Cloud(2)%Effective_Radius = geoval%vals(:,k1)

    ENDIF

 end do

 end subroutine Load_Atm_Data

! ------------------------------------------------------------------------------

subroutine Load_Sfc_Data(n_Profiles,n_Layers,N_Channels,channels,geovals,sfc,chinfo,obss)

implicit none
integer,                     intent(in)    :: n_Profiles, n_Layers, N_Channels
type(ufo_geovals),           intent(in)    :: geovals
type(CRTM_Surface_type),     intent(inout) :: sfc(:)
type(CRTM_ChannelInfo_type), intent(in)    :: chinfo(:)
type(c_ptr), value,          intent(in)    :: obss
integer(c_int),              intent(in)    :: channels(:)

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

character(len=100) :: varname_tmplate
character(len=200) :: varname

real(kind_real), allocatable :: ObsTb(:,:)

 varname_tmplate = "brightness_temperature"

 allocate(ObsTb(n_profiles, n_channels))
 ObsTb = 0.0_kind_real
 
 do n1 = 1,n_Channels
   if (any(n1==channels)) then
     !Get the variable name for this channel
     call get_var_name(varname_tmplate,n1,varname)
     call obsspace_get_db(obss, "ObsValue", varname, ObsTb(:,n1))
   endif
 enddo

 !Loop over all n_Profiles, i.e. number of locations
 do k1 = 1,N_PROFILES

   !Pass sensor information
   sfc(k1)%sensordata%sensor_id        = chinfo(1)%sensor_id
   sfc(k1)%sensordata%wmo_sensor_id    = chinfo(1)%wmo_sensor_id
   sfc(k1)%sensordata%wmo_satellite_id = chinfo(1)%wmo_satellite_id
   sfc(k1)%sensordata%sensor_channel   = chinfo(1)%sensor_channel

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

   !Land_Type
   call ufo_geovals_get_var(geovals, var_sfc_landtyp, geoval)
   sfc(k1)%Land_Type = int(geoval%vals(1,k1))

   !Land_Coverage
   call ufo_geovals_get_var(geovals, var_sfc_lfrac, geoval)
   sfc(k1)%Land_Coverage = geoval%vals(1,k1)

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

 call obsspace_get_db(obss, "MetaData", "Sat_Zenith_Angle", TmpVar)
 geo(:)%Sensor_Zenith_Angle = TmpVar(:)

 call obsspace_get_db(obss, "MetaData", "Sol_Zenith_Angle", TmpVar)
 geo(:)%Source_Zenith_Angle = TmpVar(:)

 call obsspace_get_db(obss, "MetaData", "Sat_Azimuth_Angle", TmpVar)
 geo(:)%Sensor_Azimuth_Angle = TmpVar(:)

 call obsspace_get_db(obss, "MetaData", "Sol_Azimuth_Angle", TmpVar)
 geo(:)%Source_Azimuth_Angle = TmpVar(:)

 call obsspace_get_db(obss, "MetaData", "Scan_Position", TmpVar)
 geo(:)%Ifov = TmpVar(:)

 call obsspace_get_db(obss, "MetaData", "Scan_Angle", TmpVar) !The Sensor_Scan_Angle is optional
 geo(:)%Sensor_Scan_Angle = TmpVar(:)

 deallocate(TmpVar)

end subroutine Load_Geom_Data

! ------------------------------------------------------------------------------

subroutine get_var_name(varname_tmplate,n,varname)

character(len=*), intent(in) :: varname_tmplate
integer, intent(in) :: n
character(len=*), intent(out) :: varname

character(len=3) :: chan

 ! pass in varname_tmplate = "brigtness_temperature"
 write(chan, '(I0)') n
 varname = trim(varname_tmplate) // '_' // trim(chan) // '_'

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

    REAL(kind_real), DIMENSION(n_layers,n_profiles) :: rh

    IF (TRIM(aerosol_option) == "aerosols_gocart_nasa") THEN
       varname=var_rh
       CALL ufo_geovals_get_var(geovals, varname, geoval)
       rh(1:n_layers,1:n_profiles)=geoval%vals(1:n_layers,1:n_profiles)
       WHERE (rh > 1_kind_real) rh=1_kind_real
       CALL assign_gocart_nasa
    ELSEIF (TRIM(aerosol_option) == "aerosols_gocart_esrl") THEN
       varname=var_rh
       CALL ufo_geovals_get_var(geovals, varname, geoval)
       rh(1:n_layers,1:n_profiles)=geoval%vals(1:n_layers,1:n_profiles)
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

    SUBROUTINE assign_gocart_nasa

      INTEGER, PARAMETER :: ndust_bins=5, nseas_bins=4
      REAL(kind_real), DIMENSION(ndust_bins), PARAMETER  :: dust_radii=[&
           &0.55_kind_real,1.4_kind_real,2.4_kind_real,4.5_kind_real,8.0_kind_real]
      
      INTEGER, DIMENSION(nseas_bins), PARAMETER  :: seas_types=[&
           SEASALT_SSAM_AEROSOL,SEASALT_SSCM1_AEROSOL,SEASALT_SSCM2_AEROSOL,SEASALT_SSCM3_AEROSOL]

      REAL(kind_real), DIMENSION(n_layers) :: ugkg_kgm2
      
      INTEGER :: i,k,m

      DO m=1,n_profiles

         CALL calculate_aero_layer_factor(atm(m),ugkg_kgm2)
 
         DO i=1,n_aerosols_gocart_nasa

            varname=var_aerosols_gocart_nasa(i)
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

    END SUBROUTINE assign_gocart_nasa

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

      REAL(kind_real), DIMENSION(n_layers) :: ugkg_kgm2
      
      INTEGER :: i,k,m
      
      message = 'this aerosol not implemented - check later'
      CALL Display_Message( aerosol_option, message, FAILURE )
      STOP

    END SUBROUTINE assign_other
    
  END SUBROUTINE load_aerosol_data

  SUBROUTINE assign_aerosol_names(aerosol_option,var_aerosols)
    
    CHARACTER(*), INTENT(in) :: aerosol_option
    CHARACTER(*), INTENT(out) :: var_aerosols(:)

    CHARACTER(max_string) :: err_msg

    IF (aerosol_option == "aerosols_gocart_nasa") THEN
       var_aerosols=var_aerosols_gocart_nasa
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

   SUBROUTINE check_fwd(hofx,obss,n_profiles,n_channels,varname_tmplate,channels)
     
     TYPE(c_ptr), value,       INTENT(in)    :: obss
     INTEGER, INTENT(in) :: n_profiles,n_channels
     CHARACTER(*), INTENT(in) :: varname_tmplate
     REAL(kind_real), DIMENSION(*), INTENT(in) :: hofx
     INTEGER(c_int),           INTENT(in) :: channels(:)  !List of channels to use

     REAL(kind_real), DIMENSION(n_profiles, n_channels) :: &
          &obs, hofxgsi, diff
     REAL(kind_real), DIMENSION(n_channels) :: rmse

     CHARACTER(MAXVARLEN) :: varname,varnamecombo
     CHARACTER(*), PARAMETER :: chofx="@GsiHofXBc"
     CHARACTER(*), PARAMETER :: cobsvalue="@ObsValue"
     
     INTEGER :: i,l,m

     DO l = 1,SIZE(channels)
        CALL get_var_name(varname_tmplate,channels(l),varname)
        varnamecombo=trim(varname)//chofx
        CALL obsspace_get_db(obss, "", varnamecombo, obs(:,channels(l)))
        CALL get_var_name(chofx,channels(l),varname)
        varnamecombo=TRIM(varname)//cobsvalue
        CALL obsspace_get_db(obss, "", varnamecombo, hofxgsi(:,channels(l)))
     ENDDO

     rmse = 0_kind_real

     i = 1

     DO m = 1, n_profiles
        DO l = 1, SIZE(channels)
           diff(m,channels(l)) = hofx(i) - hofxgsi(m,channels(l))
           rmse(channels(l)) = rmse(channels(l)) + diff(m,channels(l))**2
           i = i + 1
        END DO
     ENDDO

     rmse=SQRT(rmse/n_profiles)

     PRINT *,'N_profiles', N_PROFILES
     DO l = 1, SIZE(channels)
        PRINT *, 'Channel: ',channels(l)
        PRINT *, 'Max difference: ', MAXVAL(ABS(diff(:,channels(l))))
        PRINT *, 'RMSE: ', rmse(channels(l))
     ENDDO

   END SUBROUTINE check_fwd

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

   FUNCTION lower2upper(str) RESULT (string)

     IMPLICIT NONE

     CHARACTER(*), INTENT(in) :: str
     CHARACTER(LEN(str))      :: string

     INTEGER :: ic, i

     CHARACTER(26), PARAMETER :: upper = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
     CHARACTER(26), PARAMETER :: lower = 'abcdefghijklmnopqrstuvwxyz'

!   lowcase each letter if it is lowecase
     string = str
     DO i = 1, LEN_TRIM(str)
        ic = INDEX(lower, str(i:i))
        IF (ic > 0) string(i:i) = upper(ic:ic)
     END DO

   END FUNCTION lower2upper

   FUNCTION replace_text(s,text,rep) RESULT(outs) 
     CHARACTER(*)        :: s,text,rep
     CHARACTER(LEN(s)+100) :: outs  ! provide outs with extra 100 char len
     INTEGER             :: i, nt, nr

     outs = s ; nt = LEN_TRIM(text) ; nr = LEN_TRIM(rep)
     DO
        i = INDEX(outs,text(:nt)) ; IF (i == 0) EXIT
        outs = outs(:i-1) // rep(:nr) // outs(i+nt:)
     END DO

   END FUNCTION replace_text

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

 END MODULE ufo_crtm_utils_mod
