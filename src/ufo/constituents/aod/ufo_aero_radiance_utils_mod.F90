! (C) Copyright 2018 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

!> Fortran module to provide code shared between nonlinear and tlm/adm radiance calculations

MODULE ufo_aero_radiance_utils_mod

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

public aero_rad_conf
public aero_rad_conf_setup
public aero_rad_conf_delete
public Load_Atm_Data
public Load_Sfc_Data
public Load_Geom_Data

!@mzp
PUBLIC load_aerosol_data 
INTERFACE load_aerosol_data
MODULE PROCEDURE &
     &load_aerosol_data_gocart_default,&
     &load_aerosol_data_gocart_esrl,&
     &load_aerosol_data_none
END INTERFACE

integer, parameter, public :: max_string=800

!Type for general config
type aero_rad_conf
 integer :: n_Sensors
 integer :: n_Absorbers
 integer :: n_Clouds
 integer :: n_Aerosols
 integer, allocatable :: skiplist(:)
 character(len=255), allocatable :: SENSOR_ID(:)
 character(len=255) :: ENDIAN_TYPE
 character(len=255) :: COEFFICIENT_PATH
end type aero_rad_conf

contains

! ------------------------------------------------------------------------------

SUBROUTINE aero_rad_conf_setup(rc, c_conf)

implicit none
type(aero_rad_conf), intent(inout) :: rc
type(c_ptr),    intent(in)    :: c_conf

character(len=1023) :: SkipChannels
integer :: nskip, i
character(len=100), allocatable :: skiplist_str(:)

 !Some config needs to come from user
 !-----------------------------------

 !Number of sensors, each call to CRTM will be for a single sensor
 !type (zenith/scan angle will be different)
 rc%n_Sensors = 1

 !Number of absorbers, clouds and aerosols (should match what model will provide)
 rc%n_Absorbers = config_get_int(c_conf,"n_Absorbers")
 rc%n_Clouds    = config_get_int(c_conf,"n_Clouds"   )
 rc%n_Aerosols  = config_get_int(c_conf,"n_Aerosols" )

 !Allocate SENSOR_ID
 allocate(rc%SENSOR_ID(rc%n_Sensors))

 !Get sensor ID from config
 rc%SENSOR_ID(rc%n_Sensors) = config_get_string(c_conf,len(rc%SENSOR_ID(rc%n_Sensors)),"Sensor_ID")

 !ENDIAN type
 rc%ENDIAN_TYPE = config_get_string(c_conf,len(rc%ENDIAN_TYPE),"EndianType")

 !Path to coefficient files
 rc%COEFFICIENT_PATH = config_get_string(c_conf,len(rc%COEFFICIENT_PATH),"CoefficientPath")

 !Channels to skip
 if (config_element_exists(c_conf,"SkipChannels")) then
   SkipChannels = config_get_string(c_conf,len(SkipChannels),"SkipChannels")
   nskip = 1 + count(transfer(SkipChannels, 'a', len(SkipChannels)) == ",")
   allocate(skiplist_str(nskip))
   read(SkipChannels,*) skiplist_str
 else
   nskip = 0
 endif
 allocate(rc%skiplist(nskip))
 do i = 1,nskip
   read(skiplist_str(i),*)  rc%skiplist(i)
 enddo

END SUBROUTINE aero_rad_conf_setup

! -----------------------------------------------------------------------------

subroutine aero_rad_conf_delete(rc)

implicit none
type(rad_conf), intent(inout) :: rc

 deallocate(rc%SENSOR_ID)
 deallocate(rc%skiplist)

end subroutine aero_rad_conf_delete

! ------------------------------------------------------------------------------

SUBROUTINE Load_Atm_Data(n_profiles,n_layers,geovals,rc,atm)

implicit none
integer, intent(in) :: n_profiles, n_layers
type(ufo_geovals), intent(in) :: geovals
type(rad_conf), intent(in) :: rc
type(CRTM_Atmosphere_type), intent(inout) :: atm(:)


! Local variables
integer :: k1
type(ufo_geoval), pointer :: geoval
character(MAXVARLEN) :: varname
character(max_string) :: err_msg

 ! Print profile and absorber definitions
 ! --------------------------------------
 do k1 = 1,geovals%nvar
    varname = geovals%variables%fldnames(k1)
    print *, k1, varname
 end do

!@mzp
!var_aerosols - where to set up - declare - from ufo_vars_mod.F90
!possibly from config
!@mzp
 error result when not available for ufo_geovals_get_var call
 e.g.
 CALL ufo_geovals_get_var(geovals, var_oz, geoval)

 ! Populate the atmosphere structures for CRTM (atm(k1), for the k1-th profile)
 ! ----------------------------------------------------------------------------
 do k1 = 1,n_profiles
    call ufo_geovals_get_var(geovals, var_tv, geoval)

    ! Check model levels is consistent in geovals & crtm
    if (k1 == 1) then
      if (geoval%nval /= n_layers) then
        write(err_msg,*) 'Load_Atm_Data error: layers inconsistent!'
        call abor1_ftn(err_msg)
      endif
    endif

    atm(k1)%Temperature(1:n_layers) = geoval%vals(:,k1)

    call ufo_geovals_get_var(geovals, var_prs, geoval)
    atm(k1)%Pressure(1:n_layers) = geoval%vals(:,k1)
    call ufo_geovals_get_var(geovals, var_prsi, geoval)
    atm(k1)%Level_Pressure(:) = geoval%vals(:,k1)
    atm(k1)%Climatology         = US_STANDARD_ATMOSPHERE
    atm(k1)%Absorber_Id(1:1)    = (/ H2O_ID /)
    atm(k1)%Absorber_Units(1:1) = (/ MASS_MIXING_RATIO_UNITS /)
    call ufo_geovals_get_var(geovals, var_mixr, geoval)
    atm(k1)%Absorber(1:n_layers,1)       = geoval%vals(:,k1)
    atm(k1)%Absorber_Id(2:2)    = (/ O3_ID /)
    atm(k1)%Absorber_Units(2:2) = (/ VOLUME_MIXING_RATIO_UNITS /)
    call ufo_geovals_get_var(geovals, var_oz, geoval)

!@mzp
 !if not available  
    IF (error) THEN
       atm(k1)%absorber(1:n_layers,2)=ozone_fill
    ELSE
       atm(k1)%Absorber(1:n_layers,2)       = geoval%vals(:,k1)
    ENDIF

!@mzp
    IF (rc%n_Absorbers >= 3) THEN 
       atm(k1)%Absorber_Id(3:3)    = (/ CO2_ID /)
       atm(k1)%Absorber_Units(3:3) = (/ VOLUME_MIXING_RATIO_UNITS /)
       CALL ufo_geovals_get_var(geovals, var_co2, geoval)
       atm(k1)%Absorber(1:n_layers,3)       = geoval%vals(:,k1)
    ENDIF

!@mzp
    IF (rc%n_Clouds >= 1) THEN
       atm(k1)%Cloud(1)%Type = WATER_CLOUD
       CALL ufo_geovals_get_var(geovals, var_clw, geoval)
       atm(k1)%Cloud(1)%Water_Content = geoval%vals(:,k1)
       CALL ufo_geovals_get_var(geovals, var_clwefr, geoval)
       atm(k1)%Cloud(1)%Effective_Radius = geoval%vals(:,k1)
    ENDIF

    IF (rc%n_Clouds >= 2) THEN
       atm(k1)%Cloud(2)%Type = ICE_CLOUD
       CALL ufo_geovals_get_var(geovals, var_cli, geoval)
       atm(k1)%Cloud(2)%Water_Content = geoval%vals(:,k1)
       CALL ufo_geovals_get_var(geovals, var_cliefr, geoval)
       atm(k1)%Cloud(2)%Effective_Radius = geoval%vals(:,k1)
    ENDIF

    IF (rc%n_Aerosols > 0) THEN
       CALL load_aerosol_data(n_profiles,n_layers,geovals,var_aerosols,atm)       
    ENDIF

 end do

 end subroutine Load_Atm_Data

! ------------------------------------------------------------------------------

subroutine Load_Sfc_Data(n_Profiles,n_Layers,N_Channels,geovals,sfc,chinfo,obss)

implicit none
integer,                     intent(in)    :: n_Profiles, n_Layers, N_Channels
type(ufo_geovals),           intent(in)    :: geovals
type(CRTM_Surface_type),     intent(inout) :: sfc(:)
type(CRTM_ChannelInfo_type), intent(in)    :: chinfo(:)
type(c_ptr), value,          intent(in)    :: obss

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
 
 do n1 = 1,n_Channels
   !Get the variable name for this channel
   call get_var_name(varname_tmplate,n1,varname)
   call obsspace_get_db(obss, "", varname, ObsTb(:,n1))
 enddo

 !Loop over all n_Profiles, i.e. number of locations
 do k1 = 1,n_profiles

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

 call obsspace_get_db(obss, "", "Sat_Zenith_Angle", TmpVar)
 geo(:)%Sensor_Zenith_Angle = TmpVar(:)

 call obsspace_get_db(obss, "", "Sol_Zenith_Angle", TmpVar)
 geo(:)%Source_Zenith_Angle = TmpVar(:)

 call obsspace_get_db(obss, "", "Sat_Azimuth_Angle", TmpVar)
 geo(:)%Sensor_Azimuth_Angle = TmpVar(:)

 call obsspace_get_db(obss, "", "Sol_Azimuth_Angle", TmpVar)
 geo(:)%Source_Azimuth_Angle = TmpVar(:)

 call obsspace_get_db(obss, "", "Scan_Position", TmpVar)
 geo(:)%Ifov = TmpVar(:)

 call obsspace_get_db(obss, "", "Scan_Angle", TmpVar) !The Sensor_Scan_Angle is optional
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

SUBROUTINE load_aerosol_data_gocart_default(n_profiles,n_layers,geovals,&
     &var_aerosols_gocart_default,atm)

    USE CRTM_aerosolcoeff, ONLY: aeroc

    INTEGER, INTENT(in) :: n_profiles,n_layers
    TYPE(ufo_geovals), INTENT(in) :: geovals
    CHARACTER(len=MAXVARLEN), DIMENSION(n_aerosols_gocart_default), &
         &INTENT(in) :: var_aerosols_gocart_default
    TYPE(CRTM_atmosphere_type), INTENT(inout) :: atm(:)

    INTEGER, PARAMETER :: ndust_bins=5, nseas_bins=4
    REAL(kind_real), DIMENSION(ndust_bins), PARAMETER  :: dust_radii=[&
         &0.55_kind_real,1.4_kind_real,2.4_kind_real,4.5_kind_real,8.0_kind_real]
    INTEGER, DIMENSION(nseas_bins), PARAMETER  :: seas_types=[&
         SEASALT_SSAM_AEROSOL,SEASALT_SSCM1_AEROSOL,SEASALT_SSCM2_AEROSOL,    SEASALT_SSCM3_AEROSOL]

    REAL(kind_real), DIMENSION(n_layers) :: ugkg_kgm2,rh

    TYPE(ufo_geoval), POINTER :: geoval
    INTEGER :: nc, nl
    INTEGER :: k1, k2

    INTEGER :: i,k,m

    DO m=1,n_profiles

       CALL calculate_rh_and_depthfactor(atm(m),n_layers,ugkg_kgm2,rh)

       DO i=1,n_aerosols_gocart_default
          varname=var_aerosols_gocart_default(i)
          CALL ufo_geovals_get_var(geovals,varname, geoval)

          atm(m)%aerosol(i)%Concentration(1:n_layers)=&
               &MAX(geoval%vals(:,m)*ugkg_kgm2,aerosol_concentration_minvalue)

          SELECT CASE ( TRIM(varname))
          CASE ('sulf')
             atm(m)%aerosol(i)%type  = SULFATE_AEROSOL
             DO k=1,n_layers
                atm(m)%aerosol(i)%effective_radius(k)=&
                     &gocart_aerosol_size(atm(m)%aerosol(i)%type, &
                     &rh(k))
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
                     &rh(k))
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
                     &rh(k))
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
                     &rh(k))
             ENDDO
          CASE ('seas2')
             atm(m)%aerosol(i)%type  = seas_types(2)
             DO k=1,n_layers
                atm(m)%aerosol(i)%effective_radius(k)=&
                     &gocart_aerosol_size(atm(m)%aerosol(i)%type, &
                     &rh(k))
             ENDDO
          CASE ('seas3')
             atm(m)%aerosol(i)%type  = seas_types(3)
             DO k=1,n_layers
                atm(m)%aerosol(i)%effective_radius(k)=&
                     &gocart_aerosol_size(atm(m)%aerosol(i)%type, &
                     &rh(k))
             ENDDO
          CASE ('seas4')
             atm(m)%aerosol(i)%type  = seas_types(4)
             DO k=1,n_layers
                atm(m)%aerosol(i)%effective_radius(k)=&
                     &gocart_aerosol_size(atm(m)%aerosol(i)%type, &
                     &rh(k))
             ENDDO
          END SELECT

       ENDDO

    ENDDO

    NULLIFY(geoval)

  END SUBROUTINE load_aerosol_data_gocart_default

  SUBROUTINE load_aerosol_data_gocart_esrl(n_profiles,n_layers,geovals,&
       &var_aerosols_gocart_esrl,atm)
    
    USE CRTM_aerosolcoeff, ONLY: aeroc
    
    INTEGER, INTENT(in) :: n_profiles,n_layers
    TYPE(ufo_geovals), INTENT(in) :: geovals
    CHARACTER(len=MAXVARLEN), DIMENSION(n_aerosols_gocart_esrl), &
         &INTENT(in) :: var_aerosols_gocart_esrl
    TYPE(CRTM_atmosphere_type), INTENT(inout) :: atm(:)

    INTEGER, PARAMETER :: ndust_bins=5, nseas_bins=4
    REAL(kind_real), DIMENSION(ndust_bins), PARAMETER  :: dust_radii=[&
         &0.55_kind_real,1.4_kind_real,2.4_kind_real,4.5_kind_real,8.0_kind_real]
    INTEGER, DIMENSION(nseas_bins), PARAMETER  :: seas_types=[&
         SEASALT_SSAM_AEROSOL,SEASALT_SSCM1_AEROSOL,SEASALT_SSCM2_AEROSOL,    SEASALT_SSCM3_AEROSOL]

    REAL(kind_real), DIMENSION(n_layers) :: ugkg_kgm2,rh

    TYPE(ufo_geoval), POINTER :: geoval
    INTEGER :: nc, nl
    INTEGER :: k1, k2

    INTEGER :: i,k,m

    CHARACTER(*), PARAMETER :: PROGRAM_NAME = 'gocart esrl'

    message = 'gocart esrl not implemented'
    CALL Display_Message( PROGRAM_NAME, message, FAILURE )
    STOP

    DO m=1,n_profiles

       CALL calculate_rh_and_depthfactor(atm(m),n_layers,ugkg_kgm2,rh)

       DO i=1,n_aerosols_gocart_esrl
          varname=var_aerosols_gocart_default(i)
          CALL ufo_geovals_get_var(geovals,varname, geoval)

          atm(m)%aerosol(i)%Concentration(1:n_layers)=&
               &MAX(geoval%vals(:,m)*ugkg_kgm2,aerosol_concentration_minvalue)

          SELECT CASE ( TRIM(varname))
          CASE ('sulf')
             atm(m)%aerosol(i)%type  = SULFATE_AEROSOL
!rh needs to be from top to bottom
             DO k=1,n_layers
                atm(m)%aerosol(i)%effective_radius(k)=&
                     &gocart_aerosol_size(atm(m)%aerosol(i)%type, &
                     &rh(k))
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
                     &rh(k))
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
                     &rh(k))
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

!@mzp needs to be calculated between dust1 and dust2
          CASE ('p25')
             atm(m)%aerosol(i)%type  = DUST_AEROSOL
             atm(m)%aerosol(i)%effective_radius(:)=dust_radii(1)

          CASE ('seas1')
             atm(m)%aerosol(i)%type  = seas_types(1)
             DO k=1,n_layers
                atm(m)%aerosol(i)%effective_radius(k)=&
                     &gocart_aerosol_size(atm(m)%aerosol(i)%type, &
                     &rh(k))
             ENDDO
          CASE ('seas2')
             atm(m)%aerosol(i)%type  = seas_types(2)
             DO k=1,n_layers
                atm(m)%aerosol(i)%effective_radius(k)=&
                     &gocart_aerosol_size(atm(m)%aerosol(i)%type, &
                     &rh(k))
             ENDDO
          CASE ('seas3')
             atm(m)%aerosol(i)%type  = seas_types(3)
             DO k=1,n_layers
                atm(m)%aerosol(i)%effective_radius(k)=&
                     &gocart_aerosol_size(atm(m)%aerosol(i)%type, &
                     &rh(k))
             ENDDO
          CASE ('seas4')
             atm(m)%aerosol(i)%type  = seas_types(4)
             DO k=1,n_layers
                atm(m)%aerosol(i)%effective_radius(k)=&
                     &gocart_aerosol_size(atm(m)%aerosol(i)%type, &
                     &rh(k))
             ENDDO
          END SELECT

       ENDDO

    ENDDO

    NULLIFY(geoval)

  END SUBROUTINE load_aerosol_data_gocart_esrl

  SUBROUTINE load_aerosol_data_none(n_profiles,n_layers,geovals,&
       &var_aerosols_none,n_aerosols_none,atm)
    
    USE CRTM_aerosolcoeff, ONLY: aeroc
    
    INTEGER, INTENT(in) :: n_profiles,n_layers
    TYPE(ufo_geovals), INTENT(in) :: geovals
    INTEGER :: n_aerosols_none !=size(var_aerosols_none)
    CHARACTER(len=MAXVARLEN), DIMENSION(n_aerosols_none), &
         &INTENT(in) :: var_aerosols_none
    TYPE(CRTM_atmosphere_type), INTENT(inout) :: atm(:)

    CHARACTER(*), PARAMETER :: PROGRAM_NAME = 'other aerosol'

    CONTINUE

    message = 'this aerosol not implemented - check next week'
    CALL Display_Message( PROGRAM_NAME, message, FAILURE )
    STOP
    
  END SUBROUTINE load_aerosol_data_none
  
  FUNCTION gocart_aerosol_size( itype, eh ) & ! eh input in 0-1
       &RESULT(r_eff)   ! in micrometer

    USE CRTM_aerosolcoeff, ONLY: aeroc
    IMPLICIT NONE

!
!   modified from a function provided by quanhua liu
!
    INTEGER ,INTENT(in) :: itype
    REAL(kind_real)    ,INTENT(in) :: eh

    INTEGER :: j1,j2,k
    REAL(kind_real)    :: h1
    REAL(kind_real)    :: r_eff

    j2 = 0
    IF ( eh <= aeroc%rh(1) ) THEN
       j1 = 1
    ELSE IF ( eh >= aeroc%rh(aeroc%n_rh) ) THEN
       j1 = aeroc%n_rh
    ELSE
       DO k = 1, aeroc%n_rh-1
          IF ( eh < aeroc%rh(k+1) .AND. eh > aeroc%rh(k) ) THEN
             j1 = k
             j2 = k+1
             h1 = (eh-aeroc%rh(k))/(aeroc%rh(k+1)-aeroc%rh(k))
             EXIT
          ENDIF
       ENDDO
    ENDIF

    IF ( j2 == 0 ) THEN
       r_eff = aeroc%reff(j1,itype )
    ELSE
       r_eff = (one-h1)*aeroc%reff(j1,itype ) + h1*aeroc%reff(j2,itype )
    ENDIF

    RETURN

  END FUNCTION gocart_aerosol_size

  SUBROUTINE calculate_rh_and_depthfactor(atm,n_layers,ugkg_kgm2,rh)

    TYPE(CRTM_atmosphere_type), INTENT(in) :: atm
    INTEGER, INTENT(in) :: n_layers
    REAL(kind_real), DIMENSION(n_layers), INTENT(out) :: ugkg_kgm2,rh

    REAL(kind_real), DIMENSION(n_layers) :: tsen,prsl
    INTEGER :: k

!rh, ug2kg need to be from top to bottom    
!ug2kg && hPa2Pa
    DO k=1,n_layers
!correct for mixing ratio factor ugkg_kgm2 
!being calculated from dry pressure, cotton eq. (2.4)
!p_dry=p_total/(1+1.61*mixing_ratio)
       ugkg_kgm2(k)=1.0e-9_kind_real*(atm(m)%level_pressure(k)-&
            &atm(m)%level_pressure(k-1))*100_kind_real/grav/&
            &(one+eps_p1*atm(m)%absorber(k,1)*1.e-3_kind_real)
       prsl(k)=atm(m)%pressure(n_layers-k+1)*0.1_kind_real ! must be in cb for genqsat
       tsen(k)=atm(m)%temperature(n_layers-k+1)
    ENDDO
    
    CALL genqsat(qsat,tsen,prsl,n_layers,ice4qsat)
    
!relative humidity is ratio of specific humidities not mixing ratios
    DO k=1,n_layers
       rh(k)=(atm(m)%absorber(k,1)/(one+atm(m)%absorber(k,1)))*1.e-3_kind_real&
            &qsat(n_layers-k+1)
    ENDDO

    RETURN  
  
  CONTAINS

    SUBROUTINE params()
      
      LOGICAL, PARAMETER :: ice4qsat=.TRUE.

      REAL(kind_real), PARAMETER :: &
           &ttp = 2.7316e+2_kind_real, &
           &psat = 6.1078e+2_kind_real,&
           &rd = 2.8705e+2_kind_real,&
           &rv = 4.6150e+2_kind_real,&
           &cv = 7.1760e+2_kind_real,&
           &cliq = 4.1855e+3_kind_real,&
           &csol = 2.1060e+3_kind_real,&
           &cvap = 1.8460e+3_kind_real,&
           &hvap = 2.5000e+6_kind_real,&
           &hfus = 3.3358e+5_kind_real,&
           &grav = 9.81_kind_real
      
      REAL(kind_real), PARAMETER ::  &
           &tmix = ttp-20_kind_real,&
           &hsub = hvap+hfus,&
           &eps = rd/rv,&
           &eps_p1= one+eps,&
           &omeps=one-eps,&
           &dldt =cvap-cliq,&
           &dldti = cvap-csol,&
           &xa = -(dldt/rv),&
           &xai = -(dldti/rv),&
           &xb = xa+hvap/(rv*ttp),&
           &xbi = xai+hsub/(rv*ttp)
      
    END SUBROUTINE params
       
    SUBROUTINE genqsat(qsat,tsen,prsl,nsig,ice)
      
!   input argument list:
!     tsen      - input sensibile temperature field (nlat,nlon,nsig)
!     prsl      - input layer mean pressure field (nlat,nlon,nsig)
!     nsig      - number of levels                              
!     ice       - logical flag:  t=include ice and ice-water effects,
!                 depending on t, in qsat calcuations.
!                 otherwise, compute qsat with respect to water surface
!
!   output argument list:
!     qsat      - saturation specific humidity (output)
!
! remarks: see modules used
!
! attributes:
!   language: f90
!   machine:  ibm rs/6000 sp
!
      
      IMPLICIT NONE
      
      LOGICAL                               ,INTENT(in   ) :: ice
      REAL(kind_real),DIMENSION(nsig), INTENT(  out) :: qsat
      REAL(kind_real),DIMENSION(nsig),INTENT(in   ) :: tsen,prsl
      INTEGER                       ,INTENT(in   ) :: nsig
      
      INTEGER k
      REAL(kind_real) pw,tdry,tr,es,es2
      REAL(kind_real) w,onep3,esmax
      REAL(kind_real) desidt,deswdt,dwdt,desdt,esi,esw
      REAL(kind_real) :: mint,estmax
      INTEGER :: lmint
      
      onep3 = 1.e3_kind_real
      
      mint=340_kind_real
      lmint=1
      
      DO k=1,nsig
         IF((prsl(k) < 30_kind_real .AND.  &
              prsl(k) > 2_kind_real) .AND.  &
              tsen(k) < mint)THEN
            lmint=k
            mint=tsen(k)
         END IF
      END DO
      
      tdry = mint
      tr = ttp/tdry
      
      IF (tdry >= ttp .OR. .NOT. ice) THEN
         estmax = psat * (tr**xa) * EXP(xb*(one-tr))
      ELSEIF (tdry < tmix) THEN
         estmax = psat * (tr**xai) * EXP(xbi*(one-tr))
      ELSE
         w  = (tdry - tmix) / (ttp - tmix)
         estmax =  w * psat * (tr**xa) * EXP(xb*(one-tr)) &
              + (one-w) * psat * (tr**xai) * EXP(xbi*(one-tr))
      ENDIF
      
      DO k = 1,nsig
         
         tdry = tsen(k)
         tr = ttp/tdry
         IF (tdry >= ttp .OR. .NOT. ice) THEN
            es = psat * (tr**xa) * EXP(xb*(one-tr))
         ELSEIF (tdry < tmix) THEN
            es = psat * (tr**xai) * EXP(xbi*(one-tr))
         ELSE
            esw = psat * (tr**xa) * EXP(xb*(one-tr)) 
            esi = psat * (tr**xai) * EXP(xbi*(one-tr)) 
            w  = (tdry - tmix) / (ttp - tmix)
            es =  w * psat * (tr**xa) * EXP(xb*(one-tr)) &
                 + (one-w) * psat * (tr**xai) * EXP(xbi*(one-tr))
         ENDIF

         pw = onep3*prsl(k)
         esmax = es
         IF(lmint < k)THEN
            esmax=0.1_kind_real*pw
            esmax=MIN(esmax,estmax)
         END IF
         es2=MIN(es,esmax)
         qsat(k) = eps * es2 / (pw - omeps * es2)

      END DO

      RETURN

    END SUBROUTINE genqsat

  END SUBROUTINE calculate_rh_and_depthfactor

END MODULE ufo_radiance_utils_mod
