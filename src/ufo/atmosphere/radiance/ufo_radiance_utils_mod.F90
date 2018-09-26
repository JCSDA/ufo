! (C) Copyright 2018 UCAR
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 

!> Fortran module to provide code shared between nonlinear and tlm/adm radiance calculations

module ufo_radiance_utils_mod

use iso_c_binding
use config_mod
use kinds

use crtm_module

use ioda_obsdb_mod, only: ioda_obsdb, ioda_obsdb_var_to_ovec
use ioda_obs_vectors, only: obs_vector, ioda_obsvec_setup, ioda_obsvec_delete

use ufo_vars_mod
use ufo_geovals_mod, only: ufo_geovals, ufo_geoval, ufo_geovals_get_var
use ufo_basis_mod, only: ufo_basis

implicit none
private

public rad_conf
public rad_conf_setup
public rad_conf_delete
public Load_Atm_Data
public Load_Sfc_Data
public Load_Geom_Data

!Type for general config
type rad_conf
 integer :: n_Sensors
 integer :: n_Absorbers
 integer :: n_Clouds
 integer :: n_Aerosols
 character(len=255), allocatable :: SENSOR_ID(:)
 character(len=255) :: ENDIAN_TYPE
 character(len=255) :: COEFFICIENT_PATH
end type rad_conf

contains

! ------------------------------------------------------------------------------

subroutine rad_conf_setup(rc, c_conf)

implicit none
type(rad_conf), intent(inout) :: rc
type(c_ptr),    intent(in)    :: c_conf


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

end subroutine rad_conf_setup

! ------------------------------------------------------------------------------

subroutine rad_conf_delete(rc)

implicit none
type(rad_conf), intent(inout) :: rc

 deallocate(rc%SENSOR_ID)

end subroutine rad_conf_delete

! ------------------------------------------------------------------------------

subroutine Load_Atm_Data(N_PROFILES,N_LAYERS,geovals,atm)

!Internal subprogam to load some test profile data
implicit none
integer, intent(in) :: N_PROFILES, N_LAYERS
type(ufo_geovals), intent(in) :: geovals   
type(CRTM_Atmosphere_type), intent(inout) :: atm(:)

! Local variables
integer :: k1
type(ufo_geoval), pointer :: geoval
character(MAXVARLEN) :: varname


 ! Print profile and absorber definitions
 ! --------------------------------------
 do k1 = 1,geovals%nvar
    varname = geovals%variables%fldnames(k1)
    print *, k1, varname
 end do

 ! Populate the atmosphere structures for CRTM (atm(k1), for the k1-th profile)
 ! ----------------------------------------------------------------------------
 do k1 = 1,N_PROFILES
    call ufo_geovals_get_var(geovals, var_tv, geoval)
    atm(k1)%Temperature(1:N_LAYERS) = geoval%vals(:,k1) 
    call ufo_geovals_get_var(geovals, var_prs, geoval)
    atm(k1)%Pressure(1:N_LAYERS) = geoval%vals(:,k1) 
    call ufo_geovals_get_var(geovals, var_prsi, geoval)
    atm(k1)%Level_Pressure(0:N_LAYERS) = geoval%vals(:,k1)
    atm(k1)%Climatology         = US_STANDARD_ATMOSPHERE
    atm(k1)%Absorber_Id(1:1)    = (/ H2O_ID /)
    atm(k1)%Absorber_Units(1:1) = (/ MASS_MIXING_RATIO_UNITS /)
    call ufo_geovals_get_var(geovals, var_mixr, geoval)
    atm(k1)%Absorber(1:N_LAYERS,1)       = geoval%vals(:,k1) 
    atm(k1)%Absorber_Id(2:2)    = (/ O3_ID /)
    atm(k1)%Absorber_Units(2:2) = (/ VOLUME_MIXING_RATIO_UNITS /)
    call ufo_geovals_get_var(geovals, var_oz, geoval)
    atm(k1)%Absorber(1:N_LAYERS,2)       = geoval%vals(:,k1) 

    atm(k1)%Absorber_Id(3:3)    = (/ CO2_ID /)
    atm(k1)%Absorber_Units(3:3) = (/ VOLUME_MIXING_RATIO_UNITS /)
    call ufo_geovals_get_var(geovals, var_co2, geoval)
    atm(k1)%Absorber(1:N_LAYERS,3)       = geoval%vals(:,k1)

    atm(k1)%Cloud(1)%Type = WATER_CLOUD
    call ufo_geovals_get_var(geovals, var_clw, geoval)
    atm(k1)%Cloud(1)%Water_Content = geoval%vals(:,k1)
    call ufo_geovals_get_var(geovals, var_clwefr, geoval)
    atm(k1)%Cloud(1)%Effective_Radius = geoval%vals(:,k1)

    atm(k1)%Cloud(2)%Type = ICE_CLOUD
    call ufo_geovals_get_var(geovals, var_cli, geoval)
    atm(k1)%Cloud(2)%Water_Content = geoval%vals(:,k1)
    call ufo_geovals_get_var(geovals, var_cliefr, geoval)
    atm(k1)%Cloud(2)%Effective_Radius = geoval%vals(:,k1)
 end do

 end subroutine Load_Atm_Data

! ------------------------------------------------------------------------------

subroutine Load_Sfc_Data(N_PROFILES,N_LAYERS,geovals,sfc,chinfo)

!Internal subprogam to load some test profile data
implicit none
integer,                     intent(in)    :: N_PROFILES, N_LAYERS
type(ufo_geovals),           intent(in)    :: geovals   
type(CRTM_Surface_type),     intent(inout) :: sfc(:)
type(CRTM_ChannelInfo_type), intent(in)    :: chinfo(:)

type(ufo_geoval), pointer :: geoval
integer  :: k1
      
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
      
 do k1 = 1,N_PROFILES

   sfc(k1)%sensordata%sensor_id        = chinfo(1)%sensor_id
   sfc(k1)%sensordata%wmo_sensor_id    = chinfo(1)%wmo_sensor_id
   sfc(k1)%sensordata%wmo_satellite_id = chinfo(1)%wmo_satellite_id
   sfc(k1)%sensordata%sensor_channel   = chinfo(1)%sensor_channel

   sfc(k1)%Water_Type         = SEA_WATER_TYPE    !** NOTE: need to check how to determine fresh vs sea water types (salinity???)
   call                         ufo_geovals_get_var(geovals, var_sfc_wspeed, geoval)
   sfc(k1)%Wind_Speed         = geoval%vals(1,k1) 
   call                         ufo_geovals_get_var(geovals, var_sfc_wdir, geoval)
   sfc(k1)%Wind_Direction     = geoval%vals(1,k1) 
   call                         ufo_geovals_get_var(geovals, var_sfc_wfrac, geoval)
   sfc(k1)%Water_Coverage     = geoval%vals(1,k1) 
   call                         ufo_geovals_get_var(geovals, var_sfc_wtmp, geoval)
   sfc(k1)%Water_Temperature  = geoval%vals(1,k1) 
   call                         ufo_geovals_get_var(geovals, var_sfc_ifrac, geoval)
   sfc(k1)%Ice_Coverage       = geoval%vals(1,k1) 
   call                         ufo_geovals_get_var(geovals, var_sfc_itmp, geoval)
   sfc(k1)%Ice_Temperature    = geoval%vals(1,k1) 
   call                         ufo_geovals_get_var(geovals, var_sfc_sfrac, geoval)
   sfc(k1)%Snow_Coverage      = geoval%vals(1,k1) 
   call                         ufo_geovals_get_var(geovals, var_sfc_stmp, geoval)
   sfc(k1)%Snow_Temperature   = geoval%vals(1,k1) 
   call                         ufo_geovals_get_var(geovals, var_sfc_sdepth, geoval)
   sfc(k1)%Snow_Depth         = geoval%vals(1,k1)
   call                         ufo_geovals_get_var(geovals, var_sfc_landtyp, geoval)
   sfc(k1)%Land_Type          = geoval%vals(1,k1)    !** NOTE:  is this Land_Type same as CRTM's land type??
   call                         ufo_geovals_get_var(geovals, var_sfc_lfrac, geoval)
   sfc(k1)%Land_Coverage      = geoval%vals(1,k1) 
   call                         ufo_geovals_get_var(geovals, var_sfc_ltmp, geoval)
   sfc(k1)%Land_Temperature   = geoval%vals(1,k1) 
   call                         ufo_geovals_get_var(geovals, var_sfc_lai, geoval)
   sfc(k1)%Lai                = geoval%vals(1,k1) 
   call                         ufo_geovals_get_var(geovals, var_sfc_vegfrac, geoval)
   sfc(k1)%Vegetation_Fraction = geoval%vals(1,k1) 
   call                         ufo_geovals_get_var(geovals, var_sfc_vegtyp, geoval)
   sfc(k1)%Vegetation_Type    = geoval%vals(1,k1) 
   call                         ufo_geovals_get_var(geovals, var_sfc_soiltyp, geoval)
   sfc(k1)%Soil_Type          = geoval%vals(1,k1) 
   call                         ufo_geovals_get_var(geovals, var_sfc_soilm, geoval)
   sfc(k1)%Soil_Moisture_Content = geoval%vals(1,k1) 
   call                         ufo_geovals_get_var(geovals, var_sfc_soilt, geoval)
   sfc(k1)%Soil_Temperature   = geoval%vals(1,k1) 

   ! SRH
   !
   ! This is commented out, instead of deleted, in case we want to recover it in the future.
   ! The new April 15, 2018 00Z observation data contains only 12 channels (omits channels
   ! 7, 8 and 14) since that is what the April 15 GSI run assimilated. Rather than put
   ! in contrived data to fill in all 15 channels, it seemed better to disable it
   ! until we know the correct way to approach this.
   !
   ! do ch = 1, n_channels
   !   sfc(k1)%sensordata%tb(ch) = Radiance_Tbobs(ch, k1)  !** required to match GSI simulated TBs over snow and ice surfaces
   ! enddo

 end do

end subroutine Load_Sfc_Data

! ------------------------------------------------------------------------------

subroutine Load_Geom_Data(obss,geo)

!Internal subprogam to load some test geometry data
!All profiles are given the same value

implicit none
type(ioda_obsdb),         intent(in)    :: obss
type(CRTM_Geometry_type), intent(inout) :: geo(:)

type(obs_vector) :: TmpOvec

 call ioda_obsvec_setup(TmpOvec, obss%nlocs)

 call ioda_obsdb_var_to_ovec(obss, TmpOvec, "Sat_Zenith_Angle")
 geo(:)%Sensor_Zenith_Angle = TmpOvec%values(:)
 call ioda_obsdb_var_to_ovec(obss, TmpOvec, "Sol_Zenith_Angle")
 geo(:)%Source_Zenith_Angle = TmpOvec%values(:)
 call ioda_obsdb_var_to_ovec(obss, TmpOvec, "Sat_Azimuth_Angle")
 geo(:)%Sensor_Azimuth_Angle = TmpOvec%values(:)
 call ioda_obsdb_var_to_ovec(obss, TmpOvec, "Sol_Azimuth_Angle")
 geo(:)%Source_Azimuth_Angle = TmpOvec%values(:)
 call ioda_obsdb_var_to_ovec(obss, TmpOvec, "Scan_Position")
 geo(:)%Ifov = TmpOvec%values(:)
 call ioda_obsdb_var_to_ovec(obss, TmpOvec, "Scan_Angle") !The Sensor_Scan_Angle is optional
 geo(:)%Sensor_Scan_Angle = TmpOvec%values(:)

 call ioda_obsvec_delete(TmpOvec)

end subroutine Load_Geom_Data

! ------------------------------------------------------------------------------

end module ufo_radiance_utils_mod
