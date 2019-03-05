! (C) Copyright 2018 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

!> Fortran module to provide code shared between nonlinear and tlm/adm radiance calculations

MODULE ufo_radiancerttov_utils_mod

use iso_c_binding
use config_mod
use kinds

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
public load_atm_data_rttov
public load_geom_data_rttov

integer, parameter, public :: max_string=800

INTEGER, public  :: rttov_errorstatus

!Type for general config
type rad_conf
 integer              :: nsensors
 integer, allocatable :: skiplist(:)
 character(len=255), allocatable :: SENSOR_ID(:)
 character(len=255) :: COEFFICIENT_PATH
end type rad_conf

TYPE conf_type_rttov
  TYPE(rttov_coefs), ALLOCATABLE :: rttov_coef_array(:)
  TYPE(rttov_options) :: opts
  LOGICAL :: rttov_is_setup = .FALSE.
  CONTAINS 
    PROCEDURE :: set_opts => set_options_rttov
    PROCEDURE :: setup =>    setup_rttov
  END TYPE conf_type_rttov

TYPE(conf_type_rttov), PUBLIC :: config_rttov

contains

! ------------------------------------------------------------------------------

SUBROUTINE rad_conf_setup(rc, c_conf)

implicit none
type(rad_conf), intent(inout) :: rc
type(c_ptr),    intent(in)    :: c_conf

character(len=1023) :: SkipChannels
integer :: nskip, i
character(len=100), allocatable :: skiplist_str(:)

!Number of sensors, each call to RTTOV will be for a single sensor
!type (zenith/scan angle will be different)
rc%nSensors = 1

!Allocate SENSOR_ID
ALLOCATE(rc%SENSOR_ID(rc%nSensors))

!Get sensor ID from config
rc%SENSOR_ID(rc%nSensors) = config_get_string(c_conf,LEN(rc%SENSOR_ID(rc%nSensors)),"Sensor_ID")

!Path to coefficient files
rc%COEFFICIENT_PATH = config_get_string(c_conf,LEN(rc%COEFFICIENT_PATH),"CoefficientPath")

!Channels to skip
IF (config_element_exists(c_conf,"SkipChannels")) THEN
  SkipChannels = config_get_string(c_conf,LEN(SkipChannels),"SkipChannels")
  nskip = 1 + COUNT(TRANSFER(SkipChannels, 'a', LEN(SkipChannels)) == ",")
  ALLOCATE(skiplist_str(nskip))
  READ(SkipChannels,*) skiplist_str
ELSE
  nskip = 0
ENDIF
ALLOCATE(rc%skiplist(nskip))
DO i = 1,nskip
  READ(skiplist_str(i),*)  rc%skiplist(i)
ENDDO

END SUBROUTINE rad_conf_setup

! -----------------------------------------------------------------------------

SUBROUTINE rad_conf_delete(rc)
  
  IMPLICIT NONE
  TYPE(rad_conf), INTENT(inout) :: rc

  DEALLOCATE(rc%SENSOR_ID)
  DEALLOCATE(rc%skiplist)

END SUBROUTINE rad_conf_delete

! -----------------------------------------------------------------------------

SUBROUTINE load_atm_data_rttov(nprofiles,nlayers,geovals,obss,profiles)
  
  IMPLICIT NONE
  
  INTEGER, INTENT(in) :: nprofiles, nlayers
  TYPE(ufo_geovals), INTENT(in) :: geovals
  TYPE(c_ptr), VALUE,       INTENT(in)    :: obss
  TYPE(rttov_profile), INTENT(inout) :: profiles(:)
  
  ! Local variables
  INTEGER :: k1, nlocs
  INTEGER :: nlevels
  TYPE(ufo_geoval), POINTER :: geoval
  CHARACTER(MAXVARLEN) :: varname
  CHARACTER(max_string) :: err_msg
  
  REAL :: ifrac, sfrac, lfrac
  REAL :: itmp, stmp, ltmp
  REAL :: windsp
  
  REAL(kind_real), ALLOCATABLE :: TmpVar(:)
  REAL, PARAMETER :: q_mixratio_to_ppmv  = 1.60771704e+3 ! g/kg -> ppmv
  
  nlocs = obsspace_get_nlocs(obss)
  ALLOCATE(TmpVar(nlocs))
  
  nlevels = nlayers + 1
  
  DO k1 = 1, geovals%nvar
    varname = geovals%variables%fldnames(k1)
    PRINT *, k1, varname
  END DO
  
  !gas_units are ppmv (moist?)
  ! assuming wind direction 0 is N but it could be E?
  
  DO k1 = 1, nprofiles
! gas_units = 1 is mixing_ratio (moist)
! gas_units = 2 is ppmv (moist)
! gas_units is per-profile and cannot be set for individual instruments.

    profiles(k1)%gas_units = 2

    CALL ufo_geovals_get_var(geovals, var_prsi, geoval)
    profiles(k1)%p(1:nlevels) = geoval%vals(:,k1) ! hPa    


    ! ! Check model levels is consistent in geovals & RTTOV
    !   IF (k1 == 1) THEN
    !     IF (geoval%nval /= nLayers) THEN
    !       WRITE(err_msg,*) 'Load_Atm_Data error: layers inconsistent!'
    !       STOP
    !     ENDIF
    !   ENDIF


    CALL ufo_geovals_get_var(geovals, var_ts, geoval)
    profiles(k1)%t(2:nlevels) = geoval%vals(:,k1)
    profiles(k1)%t(1) = profiles(k1)%t(2)

    CALL ufo_geovals_get_var(geovals, var_mixr, geoval)
    profiles(k1)%q(2:nlevels)       = geoval%vals(:,k1)! * q_mixratio_to_ppmv 
    profiles(k1)%q(1) = profiles(k1)%q(2) 


    IF (ASSOCIATED(profiles(k1)%o3)) THEN
      CALL ufo_geovals_get_var(geovals, var_oz, geoval)
      profiles(k1)%o3(2:nlevels)       = geoval%vals(:,k1) 
      profiles(k1)%o3(1) = profiles(k1)%o3(2) 
    ENDIF


    IF (ASSOCIATED(profiles(k1)%co2)) THEN
      CALL ufo_geovals_get_var(geovals, var_co2, geoval)
      profiles(k1)%co2(2:nlevels)      = geoval%vals(:,k1)
      profiles(k1)%co2(1) = profiles(k1)%co2(2) 
    ENDIF

    IF (ASSOCIATED(profiles(k1)%co2)) THEN
      CALL ufo_geovals_get_var(geovals, var_clw, geoval)
      profiles(k1)%clw(2:nlevels) = geoval%vals(:,k1)
      profiles(k1)%clw(1) = profiles(k1)%clw(2) 
    ENDIF

    ! Near surface
    profiles(k1)%s2m%p = profiles(k1)%p(nlevels)

    !DAR: T2m, q2m is currently defined as temperature at bottom of profile. May not be appropriate.
    profiles(k1)%s2m%t = profiles(k1)%t(nlevels)
    profiles(k1)%s2m%q = profiles(k1)%q(nlevels)! * q_mixratio_to_ppmv 

    !DAR: O3_2m unused
    !        profiles(k1)%s2m%o = profiles(k1)%o3(nlevels)

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
    profiles(k1)%skin%salinity          = 35.0

    !DAR: Default fastem parameters. We are not using FASTEM over land so these are unused
    profiles(k1)%skin%fastem            = [3.0, 5.0, 15.0, 0.1, 0.3]

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

END SUBROUTINE load_atm_data_rttov

!
! Internal subprogam to load some test geometry data
!
SUBROUTINE load_geom_data_rttov(obss,profiles)
  ! Satellite viewing geometry
  ! DAR: check it's all within limits

  IMPLICIT NONE

  TYPE(c_ptr), VALUE,       INTENT(in)    :: obss
  TYPE(rttov_profile), INTENT(inout) :: profiles(:)
  REAL(kind_real), ALLOCATABLE :: TmpVar(:)
  INTEGER :: nlocs

  nlocs = obsspace_get_nlocs(obss)
  ALLOCATE(TmpVar(nlocs))

  CALL obsspace_get_db(obss, "MetaData", "Sat_Zenith_Angle", TmpVar)
  profiles(:)%zenangle = TmpVar(:)

  CALL obsspace_get_db(obss, "MetaData", "Sat_Azimuth_Angle", TmpVar)
  profiles(:)%azangle = TmpVar(:)

  CALL obsspace_get_db(obss, "MetaData", "Sol_Zenith_Angle", TmpVar)
  profiles(:)%sunzenangle = TmpVar(:)

  CALL obsspace_get_db(obss, "MetaData", "Sol_Azimuth_Angle", TmpVar)
  profiles(:)%sunazangle = TmpVar(:)


  ! DAR-RTTOV doesn't have any code to modify the viewing angle so this will need to be done elsewhere in UFO and passed to RTTOV
  ! call ioda_obsdb_var_to_ovec(UFO_Radiance, TmpOvec, "ScanPosition")
  ! profiles(:)%Ifov = TmpOvec%values(::nchannels)
  ! call ioda_obsdb_var_to_ovec(UFO_Radiance, TmpOvec, "ScanAngle")
  ! profiles(:)%Sensor_ScanAngle = TmpOvec%values(::nchannels)

  DEALLOCATE(TmpVar)
END SUBROUTINE load_geom_data_rttov

SUBROUTINE set_options_rttov(self)
  IMPLICIT NONE
  CLASS(conf_type_rttov), INTENT(INOUT) :: self
  
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
END SUBROUTINE set_options_rttov
! ------------------------------------------------------------------------------

SUBROUTINE setup_rttov(self, rc, asw)
  CLASS(conf_type_rttov)     :: self
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
    ALLOCATE(self%rttov_coef_array(rc%nSensors))
    
    DO i_inst = 1, rc%nSensors
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
  
END SUBROUTINE setup_rttov

SUBROUTINE get_var_name(varname_tmplate,n,varname)
  
  CHARACTER(len=*), INTENT(in) :: varname_tmplate
  INTEGER, INTENT(in) :: n
  CHARACTER(len=*), INTENT(out) :: varname
  
  CHARACTER(len=3) :: chan
  
  ! pass in varname_tmplate = "brightness_temperature"
  WRITE(chan, '(I0)') n
  varname = TRIM(varname_tmplate) // '_' // TRIM(chan) // '_'
  
END SUBROUTINE get_var_name

! -----------------------------------------------------------------------------

END MODULE ufo_radiancerttov_utils_mod
