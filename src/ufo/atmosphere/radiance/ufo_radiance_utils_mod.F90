! (C) Copyright 2018 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

!> Fortran module to provide code shared between nonlinear and tlm/adm radiance calculations

MODULE ufo_radiance_utils_mod

use iso_c_binding
use config_mod
use kinds

use crtm_module
USE rttov_types, ONLY : rttov_options, rttov_profile, rttov_coefs, &
                        rttov_radiance, rttov_transmission, rttov_emissivity, &
                        rttov_chanprof
USE rttov_const, ONLY : errorstatus_success, deg2rad

use ufo_vars_mod
use ufo_geovals_mod, only: ufo_geovals, ufo_geoval, ufo_geovals_get_var
use ufo_basis_mod, only: ufo_basis
use obsspace_mod

implicit none
private

public rad_conf
public rad_conf_setup
public rad_conf_delete
public Load_Atm_Data
public Load_Sfc_Data
public Load_Geom_Data
public rttov_Load_Atm_Data
public rttov_Load_Geom_Data

integer, parameter, public :: max_string=800

INTEGER, public  :: rttov_errorstatus

!Type for general config
type rad_conf
 integer :: n_Sensors
 integer :: n_Absorbers
 integer :: n_Clouds
 integer :: n_Aerosols
 integer, allocatable :: skiplist(:)
 character(len=255), allocatable :: SENSOR_ID(:)
 character(len=255) :: ENDIAN_TYPE
 character(len=255) :: COEFFICIENT_PATH
 character(len=31)  :: rtmodel
end type rad_conf

TYPE rttov_conf_type
  TYPE(rttov_coefs), ALLOCATABLE :: rttov_coef_array(:)
  TYPE(rttov_options) :: opts
  LOGICAL :: rttov_is_setup = .FALSE.
  CONTAINS 
    PROCEDURE :: set_opts => rttov_set_options
    PROCEDURE :: setup =>    rttov_setup
    
END TYPE rttov_conf_type

TYPE(rttov_conf_type), PUBLIC :: rttov_config

contains

! ------------------------------------------------------------------------------

subroutine rad_conf_setup(rc, c_conf)

implicit none
type(rad_conf), intent(inout) :: rc
type(c_ptr),    intent(in)    :: c_conf

character(len=1023) :: SkipChannels
integer :: nskip, i
character(len=100), allocatable :: skiplist_str(:)

 !Some config needs to come from user
 !-----------------------------------

 !ENDIAN type
 rc%rtmodel = config_get_string(c_conf,len(rc%rtmodel),"RTModel")

 IF (.NOT. (rc%rtmodel == 'RTTOV' .OR. rc%rtmodel == 'CRTM')) THEN
   rc%rtmodel = 'RTTOV'
 ENDIF

!COMMON

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

!CRTM specific
 IF (rc%rtmodel == 'CRTM') THEN
   !ENDIAN type
   rc%ENDIAN_TYPE = config_get_string(c_conf,LEN(rc%ENDIAN_TYPE),"EndianType")
 ENDIF

end subroutine rad_conf_setup

! -----------------------------------------------------------------------------

subroutine rad_conf_delete(rc)

implicit none
type(rad_conf), intent(inout) :: rc

 deallocate(rc%SENSOR_ID)
 deallocate(rc%skiplist)

end subroutine rad_conf_delete

! ------------------------------------------------------------------------------

subroutine Load_Atm_Data(N_PROFILES,N_LAYERS,geovals,atm)

implicit none
integer, intent(in) :: N_PROFILES, N_LAYERS
type(ufo_geovals), intent(in) :: geovals
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

    !** BTJ added 11/20/2018 for compatibility with CRTM REL 2.3.0+
    !** need to map to cloud fraction geoval, if it exists.  For now assume
    !** fully filled pixel. 
    atm(k1)%Cloud_Fraction = 1.0_fp  
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
   call obsspace_get_db(obss, "ObsValue", varname, ObsTb(:,n1))
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

SUBROUTINE RTTOV_Load_Atm_Data(N_PROFILES,N_LAYERS,geovals,obss,profiles)
  
  IMPLICIT NONE
  
  INTEGER, INTENT(in) :: N_PROFILES, N_LAYERS
  TYPE(ufo_geovals), INTENT(in) :: geovals
  TYPE(c_ptr), VALUE,       INTENT(in)    :: obss
  TYPE(rttov_profile), INTENT(inout) :: profiles(:)
  
  ! Local variables
  INTEGER :: k1, nlocs
  INTEGER :: N_LEVELS
  TYPE(ufo_geoval), POINTER :: geoval
  CHARACTER(MAXVARLEN) :: varname
  CHARACTER(max_string) :: err_msg
  
  REAL :: ifrac, sfrac, lfrac
  REAL :: itmp, stmp, ltmp
  REAL :: windsp
  
  REAL(kind_real), ALLOCATABLE :: TmpVar(:)
  REAL, PARAMETER :: q_mixratio_to_ppmv  = 1.60771704e+3_fp ! g/kg -> ppmv
  
  nlocs = obsspace_get_nlocs(obss)
  ALLOCATE(TmpVar(nlocs))
  
  N_LEVELS = N_LAYERS + 1
  
  DO k1 = 1, geovals%nvar
    varname = geovals%variables%fldnames(k1)
    PRINT *, k1, varname
  END DO
  
  !gas_units are ppmv (moist?)
  ! assuming wind direction 0 is N but it could be E?
  
  DO k1 = 1, n_profiles
! gas_units = 1 is mixing_ratio (moist)
! gas_units = 2 is ppmv (moist)
! gas_units is per-profile and cannot be set for individual instruments.

    profiles(k1)%gas_units = 2

    CALL ufo_geovals_get_var(geovals, var_prsi, geoval)
    profiles(k1)%p(1:n_levels) = geoval%vals(:,k1) ! hPa    


    ! ! Check model levels is consistent in geovals & RTTOV
    !   IF (k1 == 1) THEN
    !     IF (geoval%nval /= n_Layers) THEN
    !       WRITE(err_msg,*) 'Load_Atm_Data error: layers inconsistent!'
    !       STOP
    !     ENDIF
    !   ENDIF


    CALL ufo_geovals_get_var(geovals, var_ts, geoval)
    profiles(k1)%t(2:n_levels) = geoval%vals(:,k1)
    profiles(k1)%t(1) = profiles(k1)%t(2)

    CALL ufo_geovals_get_var(geovals, var_mixr, geoval)
    profiles(k1)%q(2:n_levels)       = geoval%vals(:,k1)! * q_mixratio_to_ppmv 
    profiles(k1)%q(1) = profiles(k1)%q(2) 


    IF (ASSOCIATED(profiles(k1)%o3)) THEN
      CALL ufo_geovals_get_var(geovals, var_oz, geoval)
      profiles(k1)%o3(2:n_levels)       = geoval%vals(:,k1) 
      profiles(k1)%o3(1) = profiles(k1)%o3(2) 
    ENDIF


    IF (ASSOCIATED(profiles(k1)%co2)) THEN
      CALL ufo_geovals_get_var(geovals, var_co2, geoval)
      profiles(k1)%co2(2:n_levels)      = geoval%vals(:,k1)
      profiles(k1)%co2(1) = profiles(k1)%co2(2) 
    ENDIF

    IF (ASSOCIATED(profiles(k1)%co2)) THEN
      CALL ufo_geovals_get_var(geovals, var_clw, geoval)
      profiles(k1)%clw(2:n_levels) = geoval%vals(:,k1)
      profiles(k1)%clw(1) = profiles(k1)%clw(2) 
    ENDIF

    ! Near surface
    profiles(k1)%s2m%p = profiles(k1)%p(n_levels)

    !DAR: T2m, q2m is currently defined as temperature at bottom of profile. May not be appropriate.
    profiles(k1)%s2m%t = profiles(k1)%t(n_levels)
    profiles(k1)%s2m%q = profiles(k1)%q(n_levels)! * q_mixratio_to_ppmv 

    !DAR: O3_2m unused
    !        profiles(k1)%s2m%o = profiles(k1)%o3(n_levels)

    ! 10m windspeed - I wonder if this will ultimately be U10 and V10 so won't need to be converted - should expect either.
    CALL ufo_geovals_get_var(geovals, var_sfc_wspeed, geoval)
    windsp = geoval%vals(1,k1)

    CALL ufo_geovals_get_var(geovals, var_sfc_wdir, geoval)       
    profiles(k1)%s2m%u             = windsp * COS(geoval%vals(1,k1) * deg2rad)
    profiles(k1)%s2m%v             = windsp * SIN(geoval%vals(1,k1) * deg2rad)

    !Skin

    !DAR: Assume ocean for now
    profiles(k1)%skin%watertype         = 1 !** NOTE: need to check how to determine fresh vs sea water types (salinity???)

    !DAR: Salinity fixed for now too
    profiles(k1)%skin%salinity          = 35.0_fp

    !DAR: Default fastem parameters. We are not using FASTEM over land so these are unused
    profiles(k1)%skin%fastem            = [3.0_fp, 5.0_fp, 15.0_fp, 0.1_fp, 0.3_fp]

    !Land point or sea point      
    CALL ufo_geovals_get_var(geovals, var_sfc_wfrac, geoval)
    IF (geoval%vals(1,k1) > 0.5) THEN
      profiles(k1)%skin%surftype   = 1 ! sea         

      CALL ufo_geovals_get_var(geovals, var_sfc_wtmp, geoval)
      profiles(k1)%skin%t   = geoval%vals(1,k1)

    ELSE ! land
      profiles(k1)%skin%surftype   = 0 ! land

      !determine land, snow and ice fractions and temperatures to determine average temperature

      CALL ufo_geovals_get_var(geovals, var_sfc_lfrac, geoval)
      lfrac   = geoval%vals(1,k1)
      CALL ufo_geovals_get_var(geovals, var_sfc_sfrac, geoval)
      sfrac   = geoval%vals(1,k1)
      CALL ufo_geovals_get_var(geovals, var_sfc_ifrac, geoval)
      ifrac   = geoval%vals(1,k1)

      CALL ufo_geovals_get_var(geovals, var_sfc_ltmp, geoval)
      ltmp   = geoval%vals(1,k1)
      CALL ufo_geovals_get_var(geovals, var_sfc_stmp, geoval)
      stmp   = geoval%vals(1,k1)
      CALL ufo_geovals_get_var(geovals, var_sfc_itmp, geoval)
      itmp   = geoval%vals(1,k1)

      !Skin temperature is a combination of (i)ce temp, (l)and temp and (s)now temp       
      profiles(k1)%skin%t   = (lfrac * ltmp + sfrac * stmp + ifrac * itmp) / (lfrac + sfrac + ifrac)

    ENDIF

    !DAR: Could/should get emissivity here?
    ! call rttov_get_emissivity()

    CALL obsspace_get_db(obss, "MetaData", "height", TmpVar)
    profiles(:)%elevation = TmpVar(:) / 1000.0 !m -> km for RTTOV

    CALL obsspace_get_db(obss, "MetaData", "latitude", TmpVar)
    profiles(:)%latitude = TmpVar(:)

    CALL obsspace_get_db(obss, "MetaData", "longitude", TmpVar)
    profiles(:)%longitude = TmpVar(:)
  ENDDO

END SUBROUTINE RTTOV_Load_Atm_Data

!
! Internal subprogam to load some test geometry data
!
SUBROUTINE rttov_Load_Geom_Data(obss,profiles)
  ! Satellite viewing geometry
  ! DAR: check it's all within limits

  IMPLICIT NONE

  TYPE(c_ptr), VALUE,       INTENT(in)    :: obss
  TYPE(rttov_profile), INTENT(inout) :: profiles(:)
  REAL(kind_real), ALLOCATABLE :: TmpVar(:)
  INTEGER :: nlocs

  nlocs = obsspace_get_nlocs(obss)
  ALLOCATE(TmpVar(nlocs))

  CALL obsspace_get_db(obss, "", "Sat_Zenith_Angle", TmpVar)
  profiles(:)%zenangle = TmpVar(:)

  CALL obsspace_get_db(obss, "", "Sat_Azimuth_Angle", TmpVar)
  profiles(:)%azangle = TmpVar(:)

  CALL obsspace_get_db(obss, "", "Sol_Zenith_Angle", TmpVar)
  profiles(:)%sunzenangle = TmpVar(:)

  CALL obsspace_get_db(obss, "", "Sol_Azimuth_Angle", TmpVar)
  profiles(:)%sunazangle = TmpVar(:)


  ! DAR-RTTOV doesn't have any code to modify the viewing angle so this will need to be done elsewhere in UFO and passed to RTTOV
  ! call ioda_obsdb_var_to_ovec(UFO_Radiance, TmpOvec, "Scan_Position")
  ! profiles(:)%Ifov = TmpOvec%values(::n_channels)
  ! call ioda_obsdb_var_to_ovec(UFO_Radiance, TmpOvec, "Scan_Angle")
  ! profiles(:)%Sensor_Scan_Angle = TmpOvec%values(::n_channels)

  DEALLOCATE(TmpVar)
END SUBROUTINE Rttov_Load_Geom_Data

SUBROUTINE rttov_set_options(self)
  IMPLICIT NONE
  CLASS(rttov_conf_type), INTENT(INOUT) :: self
  
  self % opts % rt_ir % addsolar            = .FALSE. ! Do not include solar radiation
  self % opts % interpolation % addinterp   = .TRUE.  ! Allow interpolation of input profile
  self % opts % interpolation % interp_mode = 1       ! Set interpolation method
  self % opts % rt_all % addrefrac          = .TRUE.  ! Include refraction in path calc
  self % opts % rt_ir % addclouds           = .FALSE. ! Don't include cloud effects
  self % opts % rt_ir % addaerosl           = .FALSE. ! Don't include aerosol effects
  
  self % opts % rt_ir % ozone_data          = .FALSE. ! Set the relevant flag to .TRUE.
  self % opts % rt_ir % co2_data            = .FALSE. !   when supplying a profile of the
  self % opts % rt_ir % n2o_data            = .FALSE. !   given trace gas (ensure the
  self % opts % rt_ir % ch4_data            = .FALSE. !   coef file supports the gas)
  self % opts % rt_ir % co_data             = .FALSE. !
  self % opts % rt_ir % so2_data            = .FALSE. !
  self % opts % rt_mw % clw_data            = .TRUE.  !
  
  self % opts % config % verbose            = .TRUE.  ! Enable printing of warnings
  self % opts % config % apply_reg_limits   = .TRUE.
  self % opts % config % do_checkinput      = .FALSE.
END SUBROUTINE rttov_set_options
! ------------------------------------------------------------------------------

SUBROUTINE rttov_setup(self, rc, asw)
  CLASS(rttov_conf_type)     :: self
  TYPE(rad_conf), INTENT(in) :: rc
  INTEGER, INTENT(IN) :: asw !allocate switch
  
  CHARACTER(len=255) :: coef_filename
  INTEGER :: i_inst
  
  INCLUDE 'rttov_read_coefs.interface'
  
  rttov_errorstatus = 0

  IF(asw == 1) THEN

    CALL self % set_opts()
    ! --------------------------------------------------------------------------
    ! 2. Read coefficients
    ! --------------------------------------------------------------------------    
    ALLOCATE(self%rttov_coef_array(rc%n_Sensors))
    
    DO i_inst = 1, rc%n_Sensors
      coef_filename = TRIM(rc%COEFFICIENT_PATH) // 'rtcoef_' // TRIM(rc%SENSOR_ID(i_inst)) // '.dat'
      CALL rttov_read_coefs(rttov_errorstatus, &                        !out
                              self%rttov_coef_array(i_inst), & !inout
                              self%opts, &                     !in
                              file_coef=coef_filename)                    !in
        
      IF (rttov_errorstatus /= errorstatus_success) THEN
        WRITE(*,*) 'fatal error reading coefficients'
        !          CALL rttov_exit(errorstatus)
      ELSE
        WRITE(*,*) 'successfully read' // coef_filename
      ENDIF
      
    ENDDO
    
    self%rttov_is_setup =.TRUE.
    
  ELSE !asw == 0
    !Quick and dirty for now
    DEALLOCATE(self%rttov_coef_array)
    self%rttov_is_setup =.FALSE.
  ENDIF
  
END SUBROUTINE rttov_setup

SUBROUTINE get_var_name(varname_tmplate,n,varname)
  
  CHARACTER(len=*), INTENT(in) :: varname_tmplate
  INTEGER, INTENT(in) :: n
  CHARACTER(len=*), INTENT(out) :: varname
  
  CHARACTER(len=3) :: chan
  
  ! pass in varname_tmplate = "brigtness_temperature"
  WRITE(chan, '(I0)') n
  varname = TRIM(varname_tmplate) // '_' // TRIM(chan) // '_'
  
END SUBROUTINE get_var_name

! -----------------------------------------------------------------------------

END MODULE ufo_radiance_utils_mod
