! (C) Copyright 2018 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

!> Fortran module to provide code shared between nonlinear and tlm/adm radiance calculations

module ufo_radiancerttov_utils_mod

  use, intrinsic :: iso_c_binding, only : c_double, c_ptr
  use, intrinsic :: iso_fortran_env, only : stderr=>error_unit, &
                                            stdout=>output_unit

  use datetime_mod, only : datetime, datetime_to_yyyymmddhhmmss
  use fckit_configuration_module, only : fckit_configuration
  use fckit_exception_module, only: fckit_exception
  use fckit_log_module, only : fckit_log
  use kinds, only : kind_real ! from oops
  use missing_values_mod, only : missing_value
  use obsspace_mod, only : obsspace_get_nlocs, obsspace_has, obsspace_get_db, obsspace_put_db, &
    obsspace_get_window

  use rttov_types, only : rttov_options, rttov_profile, rttov_coefs, &
    rttov_radiance, rttov_transmission, rttov_emissivity, rttov_chanprof, &
    rttov_profile_cloud, rttov_options_scatt, rttov_scatt_coef, rttov_scatt_emis_retrieval_type

  use rttov_const, only : gas_id_mixed, gas_id_watervapour, gas_id_ozone, gas_id_wvcont, gas_id_co2, &
    gas_id_n2o, gas_id_co, gas_id_ch4, gas_id_so2, ngases_max, &
    gas_name, mixratio_to_ppmv, &
    errorstatus_success, errorstatus_fatal, &
    interp_rochon, interp_rochon_wfn, ir_scatt_chou, mw_clw_scheme_liebe, mw_clw_scheme_rosenkranz, &
    vis_scatt_dom, &
    rttov_sensor_id => sensor_id, sensor_id_hi, sensor_id_ir, sensor_id_mw, &
    surftype_land, surftype_sea, surftype_seaice, watertype_ocean_water, &
    gas_unit_ppmvdry, gas_unit_specconc, &
    rttov_inst_name => inst_name, rttov_platform_name => platform_name

  use ufo_geovals_mod, only : ufo_geovals, ufo_geoval, ufo_geovals_get_var
  use ufo_utils_mod, only : Ops_SatRad_Qsplit, Ops_Qsat, Ops_QsatWat, cmp_strings, getindex, upper2lower
  use ufo_constants_mod, only : zero, half, one, deg2rad, min_q, min_clw, min_ciw, m_to_km, &
                                g_to_kg, pa_to_hpa, RTTOV_ToA

  use ufo_vars_mod, only : maxvarlen, &
    var_prs, var_ts, var_sfc_t2m, var_sfc_u10, var_sfc_v10, var_ps, &
    var_sfc_q2m, var_sfc_tskin, var_prsi, var_clw, var_cli, var_cldfrac_vol, &
    var_q, var_mixr, var_oz, var_co2, &
    var_radiance, var_tb_clr, var_tb, var_sfc_emiss, var_pmaxlev_weightfunc, var_total_transmit, &
    var_sfc_wdir, var_sfc_wspeed, &
    var_opt_depth, var_lvl_transmit, var_lvl_weightfunc, var_tb_overcast, &
    ufo_vars_getindex

  implicit none
  private

  public rttov_conf
  public rttov_conf_setup
  public rttov_conf_delete
  public rttov_hofxdiags
  public rttov_read_emissivity_from_obsspace

  ! Shared Parameters
  integer, parameter                    :: max_string=800
  integer, parameter                    :: maxvarin = 50

  character(len=maxvarlen), parameter   :: var_surf_type_rttov = "surfaceQualifier"  ! 0 (land), 1 (water), 2 (sea-ice)

  character(len=maxvarlen), dimension(8), public :: varin_default = &
    (/var_prs, var_ts, var_sfc_t2m, &
    var_sfc_u10, var_sfc_v10, var_ps, var_sfc_q2m, &
    var_sfc_tskin /)

  character(len=maxvarlen), dimension(4), public :: varin_scatt = &
    (/var_prsi, var_clw, var_cli, var_cldfrac_vol /)

  ! copy of ABSORBER_ID_NAME defined in rttov_const
  character(len=*), parameter :: &
    RTTOV_Absorbers(ngases_max+2) = &
    [character(len=12) :: gas_name(1:ngases_max),'CLW', &
     'CIW']

  integer, parameter :: &
    RTTOV_Absorber_Id(ngases_max+2) = &
    [gas_id_mixed, &
    gas_id_watervapour, &
    gas_id_ozone,  &
    gas_id_wvcont,  &
    gas_id_co2,  &
    gas_id_n2o,  &
    gas_id_co, &
    gas_id_ch4,  &
    gas_id_so2, 0, 0]

  real(kind_real), parameter :: &
    gas_unit_conv(0:ngases_max) = &
    [1.0_kind_real, & ! 0 index for use with CLW/CIW
    mixratio_to_ppmv(1:ngases_max)] ! N.B. uses ufo-specific array in rttov_const

  character(len=MAXVARLEN), parameter :: null_str = ''

  character(len=MAXVARLEN), parameter :: &
    UFO_Absorbers(ngases_max+2) = &
    [character(len=MAXVARLEN) :: null_str, var_q, var_oz, null_str, var_co2, 'mole_fraction_of_nitrous_oxide_in_air', &
    'mole_fraction_of_carbon_monoxide_in_air', 'mole_fraction_of_methane_in_air', &
    'mole_fraction_of_sulfur_dioxide_in_air', var_clw, var_cli]

  ! Types used
  type, public :: mw_scatt_io

    integer, pointer :: freq_indices(:)
    type(rttov_profile_cloud), pointer :: profiles(:)
    type(rttov_profile_cloud), pointer :: profiles_k(:)
    type(rttov_scatt_emis_retrieval_type) :: emis_retrieval

  contains

    procedure :: set_frequencies => set_freq_indices
      
  end type mw_scatt_io

  type, public :: ufo_rttov_io

    logical, pointer                 :: calcemis(:)     ! Flag to indicate calculation of emissivity within RTTOV
    type(rttov_emissivity), pointer  :: emissivity(:)   ! Input/output surface emissivity
    type(rttov_profile), allocatable :: profiles(:)     ! Input profiles
    type(rttov_profile), allocatable :: profiles_k(:)   ! Output jacobian profiles
    type(rttov_chanprof), pointer    :: chanprof(:)
    type(rttov_transmission)         :: transmission    ! Output transmittances
    type(rttov_radiance)             :: radiance        ! Output radiances

    type(rttov_emissivity), pointer  :: emissivity_k(:) ! Input/output surface emissivity
    type(rttov_transmission)         :: transmission_k  ! Output transmittances
    type(rttov_radiance)             :: radiance_k      ! Output radiances

    type(mw_scatt_io)                :: mw_scatt
    real(kind_real), allocatable     :: ciw(:,:)        ! pointer to either ice from RTTOV-SCATT or diagnosed ice from Qsplit
    real(kind_real), allocatable     :: tc_ozone(:)     ! total column ozone
    logical, allocatable             :: q_profile_reset(:,:)    ! flag to say if the humidity has been reset to a minimum value
    logical, allocatable             :: clw_profile_reset(:,:)  ! flag to say if cloud liquid water has been reset to a minimum value
    logical, allocatable             :: ciw_profile_reset(:,:)  ! flag to say if cloud ice has been reset to a minimum value

    integer, allocatable             :: sensor_index_array(:)

    integer                          :: nchan_inst  ! number of channels being simulated (may be less than full instrument)
    integer                          :: nlocs_total ! nprofiles (including skipped)

  contains

    procedure :: alloc_direct            => ufo_rttov_alloc_direct
    procedure :: alloc_k                 => ufo_rttov_alloc_k
    procedure :: alloc_profiles          => ufo_rttov_alloc_profiles
    procedure :: alloc_profiles_k        => ufo_rttov_alloc_profiles_k
    procedure :: zero_k                  => ufo_rttov_zero_k
    procedure :: init_default_emissivity => ufo_rttov_init_default_emissivity 
    procedure :: setup_rtprof            => ufo_rttov_setup_rtprof
    procedure :: check_rtprof            => ufo_rttov_check_rtprof
    procedure :: print_rtprof            => ufo_rttov_print_rtprof
    procedure :: scale_ozone             => ufo_rttov_scale_ozone
    procedure :: calculate_tc_ozone      => ufo_rttov_calculate_tc_ozone

  end type ufo_rttov_io

  type mw_scatt_conf
    type(rttov_scatt_coef)                :: coef
    type(rttov_options_scatt)             :: opts
    logical                               :: use_totalice
    logical                               :: mmr_snowrain
  end type mw_scatt_conf

  !Type for general config
  type rttov_conf
    integer                               :: nsensors
    integer                               :: ngas

    character(len=MAXVARLEN), allocatable :: Absorbers(:)
    integer, allocatable                  :: Absorber_Id(:)
    real(kind_real)                       :: scale_fac(0:ngases_max)
    logical                               :: RTTOV_GasUnitConv

    integer, allocatable                  :: wmo_id(:)
    character(len=255)                    :: COEFFICIENT_PATH
    character(len=255), allocatable       :: coeffname(:)
    integer,            allocatable       :: instrument_triplet(:,:)
    integer,            allocatable       :: rttov_sensor_type(:)

    type(rttov_coefs),  allocatable       :: rttov_coef_array(:)
    character(len=10)                     :: RTTOV_default_opts
    type(rttov_options)                   :: rttov_opts
    type(mw_scatt_conf)                   :: mw_scatt
    logical                               :: rttov_is_setup

    logical                               :: SatRad_compatibility
    logical                               :: UseRHwaterForQC  ! only used with SatRad compatibility
    logical                               :: UseColdSurfaceCheck  ! to replicate pre-PS45 results
    logical                               :: BoundQToSaturation  ! keep values of q equal or below to saturation
    logical                               :: UseMinimumQ ! apply a lower threshold to the humidity
    logical                               :: UseMinimumClw ! apply a lower threshold to cloud liquid water
    logical                               :: UseMinimumCiw ! apply a lower threshold to cloud ice
    logical                               :: SplitQtotal
    logical                               :: UseQtsplitRain
    logical                               :: RTTOV_profile_checkinput
    logical                               :: Do_MW_Scatt
    logical                               :: UseQCFlagsToSkipHofX

    logical                               :: prof_by_prof
    logical                               :: RTTOV_scale_ref_ozone
    logical                               :: debug       ! debug output

    integer, allocatable                  :: inspect(:)
    integer                               :: nchan_max_sim

    character(len=255)                    :: surface_emissivity_group
    character(len=MAXVARLEN), allocatable :: variablesFromObsSpace(:)

    real(kind_real)                       :: MWScattZeroJacPress ! zero cloud jacobian above threshold pressure level

  contains

    procedure :: set_options => ufo_rttov_set_options
    procedure :: setup =>    ufo_rttov_setup
    procedure :: set_defaults => ufo_rttov_set_defaults

  end type rttov_conf

  type rttov_hofxdiags
    ! ystr_diags contains the name of the rttov variable for ouput e.g.
    ! transmitance or optical depth.  For jacobian output var_tb
    character(len=maxvarlen), allocatable :: ystr_diags(:)
    ! xstr_diags contains the model variable for jacobian output e.g.
    ! var_t or empty if jacobian output is not required
    character(len=maxvarlen), allocatable :: xstr_diags(:)
    integer, allocatable                  :: ch_diags(:)

  contains

    procedure :: reset => ufo_rttov_hofxdiags_delete
    procedure :: parse => ufo_rttov_parse_hofxdiags
    procedure :: populate => ufo_rttov_populate_hofxdiags

  end type rttov_hofxdiags

contains

  ! ------------------------------------------------------------------------------

  subroutine rttov_conf_setup(conf, f_confOpts, f_confOper)
    implicit none

    type(rttov_conf), intent(inout)       :: conf
    type(fckit_configuration), intent(in) :: f_confOpts   ! RTcontrol
    type(fckit_configuration), intent(in) :: f_confOper

    character(*), parameter               :: routine_name = 'rttov_conf_setup'
    character(len=max_string)             :: message
    integer                               :: ivar, jspec
    integer                               :: isensor
    character(len=:), allocatable         :: str
    character(len=:), allocatable         :: tmp_str_array(:)
    character(len=16), allocatable        :: str_array(:)
    character(len=30)                     :: absorber_name

    integer                               :: delim

    logical, allocatable                  :: absorber_mask(:)
    integer                               :: rttov_errorstatus

    include 'rttov_user_options_checkinput.interface'
    include 'rttov_coeffname.interface'

    ! set some defaults
    conf % rttov_is_setup = .false.
    conf % SplitQtotal = .false.

    call f_confOper % get_or_die("Debug",conf % debug)

    ! Absorbers
    !----------

    ! q is mandatory and us expected to be present so is not required to be present in the list of absorbers
    conf%ngas = 1
    absorber_name = "Absorbers"

    ! add additional requested absorbers
    if (f_confOper%has(trim(absorber_name))) &
      conf%ngas = conf%ngas + f_confOper%get_size(trim(absorber_name))

    if (conf%ngas > 1) then
      call f_confOper%get_or_die("Absorbers", tmp_str_array) ! tmp_str_array must be deferred size
      allocate(str_array(size(tmp_str_array)))               ! but str_array must be fixed! 
      str_array(:) = tmp_str_array(:)                              

      allocate(absorber_mask(size(str_array)))
      absorber_mask = .false.

      ! water vapour is mandatory and already included in the RTTOV input list and not required 
      ! in the absorber list as well. It is removed here.
      if ( any(RTTOV_Absorbers(gas_id_watervapour) == str_array)  ) then
        message = trim(routine_name) // trim(RTTOV_Absorbers(gas_id_watervapour)) // &
                  ' is mandatory and not required to be listed in Absorbers'
        call fckit_log%info(message)
        conf%ngas = conf%ngas - 1
      end if

      ! check for duplications
      if (size(str_array) > 1) then
        do jspec = 1, size(str_array) - 1
          if (any(trim(str_array(jspec)) == str_array(jspec+1:)) ) then
            message = trim(routine_name) // trim(str_array(jspec)) // ' is duplicated in Absorbers'
            call fckit_log%info(message)
            conf%ngas = conf%ngas - 1
          else
            absorber_mask(jspec) = .true.
          end if
        enddo
        absorber_mask(size(str_array)) = .true. ! the last element cannot be a duplicate
      else
        absorber_mask(:) = .true.
      end if

    end if

    allocate(conf%Absorbers( conf%ngas ), conf%Absorber_Id( conf%ngas ))

    conf%Absorbers(1) = RTTOV_Absorbers(gas_id_watervapour)
    if (allocated(absorber_mask) .and. conf%ngas > 1) then
      conf%Absorbers(2:conf%ngas) = pack(str_array, absorber_mask)
      if (allocated(tmp_str_array)) deallocate(tmp_str_array)
      if (allocated(str_array)) deallocate(str_array)
      deallocate(absorber_mask)
    end if
    ! convert from RTTOV names to UFO CF names and define Id and Units
    do jspec = 1, conf%ngas
      ivar = ufo_vars_getindex(RTTOV_Absorbers, conf%Absorbers(jspec))

      if (ivar < 1 .or. ivar > size(UFO_Absorbers)) then
        message = trim(routine_name) // ' error: ' // trim(conf%Absorbers(jspec)) // ' not supported by UFO_Absorbers'
        call abor1_ftn(message)
      end if

      conf%Absorbers(jspec) = UFO_Absorbers(ivar)

      conf%Absorber_Id(jspec) = RTTOV_Absorber_Id(ivar)
    end do

    ! set scalar mixing ratio conversion if converting units prior to use in RTTOV
    call f_confOpts % get_or_die("RTTOV_GasUnitConv",conf % RTTOV_GasUnitConv)

    if(conf%RTTOV_GasUnitConv) then 
      conf%scale_fac = gas_unit_conv
    else
      conf%scale_fac = one
    end if

    ! use scaled RTTOV reference profile instead of reading from geovals if Ozone is a required Absorber 
    call f_confOpts % get_or_die("RTTOV_ScaleRefOzone",conf % RTTOV_scale_ref_ozone)

    ! Determine which platform and instrument will be processed by RTTOV. If WMO_ID is present then
    ! it will be used in preference to Platform_Name/Sat_ID and enables the same instrument on
    ! multiple platforms to be processed.
    call f_confOpts % get_or_die("WMO_ID", conf % wmo_id)

    call f_confOpts % get_or_die ("Sat_ID", tmp_str_array) ! tmp_str_array must be deferred size
    allocate(str_array(size(tmp_str_array)))               ! but str_array must be fixed! 
    str_array(:) = tmp_str_array(:)                              

    ! The number of sensors
    conf % nSensors = size(str_array)
    if (.not. (conf % nSensors == size(conf%wmo_id) .or. conf % nSensors == 1)) then
      message = trim(routine_name) // 'Error. Number of Sat_IDs must match WMO_IDs or be ' // &
                                      '1 (process all obs with the same coefficient) '
      call abor1_ftn(message)
    end if

    allocate(conf%instrument_triplet(3,conf%nsensors))

    ! Convert WMO identifier to RTTOV specific platform/satellite id using a lookup table.
    ! An instrument must be added to the case statement in wmo_id_to_rttov_platform if it is to
    ! be processed.
    do isensor=1, conf % nSensors
      delim = index(str_array(isensor), '_')
      if(delim > 0) then
        conf % instrument_triplet(1,isensor) = getindex(rttov_platform_name,upper2lower(str_array(isensor)(1:delim-1)))
        read(str_array(isensor)(delim+1:),*) conf % instrument_triplet(2,isensor)
      else
        message = trim(routine_name) // 'Error. Invalid format for Sat_ID: ' // trim(str_array(isensor))
        call abor1_ftn(message)
      end if
    end do

    call f_confOpts % get_or_die("Instrument_Name",str)   
    ! inst array starts from 0 not 1 so subtract 1
    conf % instrument_triplet(3,:) = getindex(rttov_inst_name, upper2lower(trim(str))) - 1

    allocate(conf%coeffname(conf%nsensors),conf%rttov_sensor_type(conf%nsensors))
    
    ! For now we only support processing with one instrument type but it would be possible to 
    ! extend this to full generic processing
    do isensor = 1, conf%nSensors

      ! get rtcoef name from instrument triplet 
      CALL rttov_coeffname (rttov_errorstatus,  &  ! out
        conf % instrument_triplet(1:3,isensor), &  ! in
        filetype = 'rtcoef',                    &  ! in
        coeffname = conf % coeffname(isensor))     ! out

      !IR=1/MW=2/HI=3/PO=4
      conf % rttov_sensor_type(isensor) = rttov_sensor_id(conf % instrument_triplet(3,isensor))
    enddo

    ! Path to coefficient files
    call f_confOpts % get_or_die("CoefficientPath",str)
    conf % COEFFICIENT_PATH = str

    ! Set interface options - SatRad_compatibility dependent variables
    call f_confOpts % get_or_die("SatRad_compatibility",conf % SatRad_compatibility) ! default = true
    call f_confOpts % get_or_die("UseColdSurfaceCheck",conf % UseColdSurfaceCheck) ! default = false
    call f_confOpts % get_or_die("BoundQToSaturation",conf % BoundQToSaturation) ! default = true
    call f_confOpts % get_or_die("UseRHwaterForQC",conf % UseRHwaterForQC) ! default = true
    call f_confOpts % get_or_die("UseMinimumQ",conf % UseMinimumQ) ! default = true
    call f_confOpts % get_or_die("UseMinimumClw",conf % UseMinimumClw) ! default = true
    call f_confOpts % get_or_die("UseMinimumCiw",conf % UseMinimumCiw) ! default = true
    call f_confOpts % get_or_die("MW_Scatt_zero_jac_press",conf % MWScattZeroJacPress) ! default = 0.0

    ! Settings that control the batching
    call f_confOpts % get_or_die("prof_by_prof",conf % prof_by_prof)
    call f_confOpts % get_or_die("max_channels_per_batch",conf % nchan_max_sim)

    if (conf % nSensors > 1) then
      message = 'Where more than one sensor is processed, fall back to profile-by-profile processing. Setting prof_by_prof to TRUE'
      call fckit_log%info(message)
      conf % prof_by_prof = .true.
    end if

    call f_confOpts % get_or_die("RTTOV_profile_checkinput", conf % RTTOV_profile_checkinput)

    ! Default options (e.g. UKMO_PS45)
    call f_confOpts % get_or_die("RTTOV_default_opts",str)
    conf % RTTOV_default_opts = str
    
    call f_confOpts % get_or_die("Do_MW_Scatt", conf % do_mw_scatt)
    conf % do_mw_scatt = conf % do_mw_scatt .and. any(conf % rttov_sensor_type(:) == sensor_id_mw)

    if (conf % do_mw_scatt .and. .not. conf % prof_by_prof) then
      message = 'RTTOV-SCATT does not support batch processing. Setting prof_by_prof to TRUE'
      call fckit_log%info(message)
      conf % prof_by_prof = .true.
    end if

    ! Decide whether to produce HofX when no channels are active for a profile
    call f_confOpts % get_or_die("UseQCFlagsToSkipHofX", conf % UseQCFlagsToSkipHofX)
    if (conf % UseQCFlagsToSkipHofX .and. .not. conf % prof_by_prof) then
      message = 'UseQCFlagsToSkipHofX does not support batch processing. Setting prof_by_prof to TRUE'
      call fckit_log%info(message)
      conf % prof_by_prof = .true.
    end if

    ! Populate RTTOV operator default options and overrides; read coefficients
    if(.not. conf % rttov_is_setup) then
      call conf % setup(f_confOpts)
    end if

    ! MW specific settings
    if (any(conf % rttov_coef_array(:) % coef % id_sensor == sensor_id_mw)) then
      conf % UseQtsplitRain = .false.
      conf % SplitQtotal = .false.

      if(conf % SatRad_compatibility) then
        conf % SplitQtotal = .true.
        if(.not. conf % Do_MW_Scatt) conf % UseQtsplitRain = .true.
      end if

      ! additional overrides
      conf % rttov_opts % rt_ir % ozone_data = .false.
      conf % rttov_opts % rt_ir % co2_data = .false.
      conf % rttov_opts % rt_ir % n2o_data = .false.
      conf % rttov_opts % rt_ir % ch4_data = .false.
      conf % rttov_opts % rt_ir % so2_data = .false.
    end if
    
    ! Ensure the RTTOV options and coefficients are consistent
    call rttov_user_options_checkinput(rttov_errorstatus, conf % rttov_opts, conf % rttov_coef_array(1))

    if (rttov_errorstatus /= errorstatus_success) then
      message = trim(routine_name) // ': Error in rttov_user_options_checkinput'
      call abor1_ftn(message)
    end if

    ! Default is false; satrad compatibility and mw default is true
    if(f_confOpts % has("QtSplitRain")) then
      call f_confOpts % get_or_die("QtSplitRain", conf % UseQtsplitRain)
    end if

    if (f_confOpts%has("InspectProfileNumber")) then
      call f_confOpts % get_or_die("InspectProfileNumber", conf % inspect)
    else
      allocate(conf % inspect(0))
    end if

    call f_confOpts % get_or_die("surface emissivity group", str)
    conf % surface_emissivity_group = str

    if (f_confOper%has("variables to use from obsspace")) then
      call f_confOper % get_or_die("variables to use from obsspace", tmp_str_array)
      allocate(conf % variablesFromObsSpace(size(tmp_str_array)))
      conf % variablesFromObsSpace(:) = tmp_str_array(:)
    else
      allocate(conf % variablesFromObsSpace(0))
    end if

  end subroutine rttov_conf_setup

  ! -----------------------------------------------------------------------------

  subroutine rttov_conf_delete(conf)

    implicit none
    type(rttov_conf), intent(inout) :: conf
    
    integer                         :: i, rttov_errorstatus

    include 'rttov_dealloc_coefs.interface'

    if (allocated(conf % rttov_coef_array)) then
      do i = 1, size(conf % rttov_coef_array)
        call rttov_dealloc_coefs(rttov_errorstatus, conf % rttov_coef_array(i))
      enddo
      deallocate(conf % rttov_coef_array)
    end if

    ! needed to prevent bugs caused when more than one obs spaces in a yaml file
    if (allocated(conf%Absorber_Id)) deallocate (conf%Absorber_Id)
    if (allocated(conf%Absorbers))   deallocate (conf%Absorbers)

    ! Reset to defaults
    conf % rttov_is_setup = .false.
    conf % SplitQtotal = .false.

  end subroutine rttov_conf_delete

  ! -----------------------------------------------------------------------------

  subroutine ufo_rttov_set_options(self, f_confOpts)
    implicit none

    class(rttov_conf),         intent(inout) :: self
    type(fckit_configuration), intent(in)    :: f_confOpts ! RTcontrol

    include 'rttov_print_opts.interface'
    include 'rttov_print_opts_scatt.interface'

    call self % set_defaults(self%RTTOV_default_opts) ! test, OPS, RTTOV

    !< Switch to enable atmospheric refraction
    if (f_confOpts % has("RTTOV_addrefrac")) then
      call f_confOpts % get_or_die("RTTOV_addrefrac", self % rttov_opts % rt_all % addrefrac)
    end if

    !< Switch for input units in AD/K models
    if (f_confOpts % has("RTTOV_switchrad")) then
      call f_confOpts % get_or_die("RTTOV_switchrad", self % rttov_opts % rt_all % switchrad)
    end if

    !< Switch to enable use of 2m q variable
    if (f_confOpts % has("RTTOV_use_q2m")) then
      call f_confOpts % get_or_die("RTTOV_use_q2m", self % rttov_opts % rt_all % use_q2m)
    end if

    !< Switch for setting Lambertian reflection (IR and MW)
    if (f_confOpts % has("RTTOV_do_lambertian")) then
      call f_confOpts % get_or_die("RTTOV_do_lambertian", self % rttov_opts % rt_all % do_lambertian)
    end if

    !< Switch for fixed/parameterised effective angle for Lambertian option
    if (f_confOpts % has("RTTOV_lambertian_fixed_angle")) then
      call f_confOpts % get_or_die("RTTOV_lambertian_fixed_angle", self % rttov_opts % rt_all % lambertian_fixed_angle)
    end if

    !< Switch to ignore atmospheric curvature
    if (f_confOpts % has("RTTOV_plane_parallel")) then
      call f_confOpts % get_or_die("RTTOV_plane_parallel", self % rttov_opts % rt_all % plane_parallel)
    end if

    !< Linear-in-tau or layer-mean for downwelling radiances
    if (f_confOpts % has("RTTOV_rad_down_lin_tau")) then
      call f_confOpts % get_or_die("RTTOV_rad_down_lin_tau", self % rttov_opts % rt_all % rad_down_lin_tau)
    end if

    !< Switch to apply dtau test in transmit/integrate calculations
    if (f_confOpts % has("RTTOV_dtau_test")) then
      call f_confOpts % get_or_die("RTTOV_dtau_test", self % rttov_opts % rt_all % dtau_test)
    end if

    !< FASTEM version (0-6); 0 => TESSEM2
    if (f_confOpts % has("RTTOV_fastem_version")) then
      call f_confOpts % get_or_die("RTTOV_fastem_version", self % rttov_opts % rt_mw % fastem_version)
    end if

    !< Supply a foam fraction to FASTEM
    if (f_confOpts % has("RTTOV_supply_foam_fraction")) then
      call f_confOpts % get_or_die("RTTOV_supply_foam_fraction", self % rttov_opts % rt_mw % supply_foam_fraction)
    end if

    !< Switch to enable input of cloud liquid water profile
    if (f_confOpts % has("RTTOV_clw_data")) then
      call f_confOpts % get_or_die("RTTOV_clw_data", self % rttov_opts % rt_mw % clw_data)
    end if

    !< MW CLW scheme: 1 => Liebe, 2 => Rosenkranz, 3 => TKC
    if (f_confOpts % has("RTTOV_clw_scheme")) then
      call f_confOpts % get_or_die("RTTOV_clw_scheme", self % rttov_opts % rt_mw % clw_scheme)
    end if

    !< Apply MW CLW calculations on coef/user levels (true/false resp.)
    if (f_confOpts % has("RTTOV_clw_calc_on_coef_lev")) then
      call f_confOpts % get_or_die("RTTOV_clw_calc_on_coef_lev", self % rttov_opts % rt_mw % clw_calc_on_coef_lev)
    end if

    !< Lower pressure limit for MW CLW calculations (hPa)
    if (f_confOpts % has("RTTOV_clw_cloud_top")) then
      call f_confOpts % get_or_die("RTTOV_clw_cloud_top", self % rttov_opts % rt_mw % clw_cloud_top)
    end if

    !< Apply band-correction for Planck radiance and BT calculations
    if (f_confOpts % has("RTTOV_apply_band_correction")) then
      call f_confOpts % get_or_die("RTTOV_apply_band_correction", self % rttov_opts % rt_mw % apply_band_correction)
    end if

    !< Switch to enable RTTOV interpolator
    if (f_confOpts % has("RTTOV_addinterp")) then
      call f_confOpts % get_or_die("RTTOV_addinterp", self % rttov_opts % interpolation % addinterp)
    end if

    !< Interpolation mode (1-5, see user guide)
    if (f_confOpts % has("RTTOV_interp_mode")) then
      call f_confOpts % get_or_die("RTTOV_interp_mode", self % rttov_opts % interpolation % interp_mode)
    end if

    !< Switch to make pressure an active variable in TL/AD/K models
    if (f_confOpts % has("RTTOV_lgradp")) then
      call f_confOpts % get_or_die("RTTOV_lgradp", self % rttov_opts % interpolation % lgradp)
    end if

    !< Switch to assume space boundary at top-most input pressure l
    if (f_confOpts % has("RTTOV_spacetop")) then
      call f_confOpts % get_or_die("RTTOV_spacetop", self % rttov_opts % interpolation % spacetop)
    end if

    !< Switch to extrapolate input profiles using regression limits
    if (f_confOpts % has("RTTOV_reg_limit_extrap")) then
      call f_confOpts % get_or_die("RTTOV_reg_limit_extrap", self % rttov_opts % interpolation % reg_limit_extrap)
    end if

    if (f_confOpts % has("RTTOV_fix_hgpl")) then
      call f_confOpts % get_or_die("RTTOV_fix_hgpl", self % rttov_opts % config % fix_hgpl)
    end if

    if (f_confOpts % has("RTTOV_verbose")) then
      call f_confOpts % get_or_die("RTTOV_verbose", self % rttov_opts % config % verbose)
    end if

    if (f_confOpts % has("RTTOV_do_checkinput")) then
      call f_confOpts % get_or_die("RTTOV_do_checkinput", self % rttov_opts % config % do_checkinput)
    end if

    if (f_confOpts % has("RTTOV_apply_reg_limits")) then
      call f_confOpts % get_or_die("RTTOV_apply_reg_limits", self % rttov_opts % config % apply_reg_limits)
    end if

    !< Solar sea BRDF model (1-2)
    if (f_confOpts % has("RTTOV_solar_sea_brdf_model")) then
      call f_confOpts % get_or_die("RTTOV_solar_sea_brdf_model", self % rttov_opts % rt_ir % solar_sea_brdf_model)
    end if

    !< IR sea emissivity model (1-2)
    if (f_confOpts % has("RTTOV_ir_sea_emis_model")) then
      call f_confOpts % get_or_die("RTTOV_ir_sea_emis_model", self % rttov_opts % rt_ir % ir_sea_emis_model)
    end if

    !< Switch to enable solar simulations
    if (f_confOpts % has("RTTOV_addsolar")) then
      call f_confOpts % get_or_die("RTTOV_addsolar", self % rttov_opts % rt_ir % addsolar)
    end if

    !< Switch to enable Rayleigh single-scattering for VIS/NIR channel
    if (f_confOpts % has("RTTOV_rayleigh_single_scatt")) then
      call f_confOpts % get_or_die("RTTOV_rayleigh_single_scatt", self % rttov_opts % rt_ir % rayleigh_single_scatt)
    end if

    !< Switch to enable NLTE bias correction
    if (f_confOpts % has("RTTOV_do_nlte_correction")) then
      call f_confOpts % get_or_die("RTTOV_do_nlte_correction", self % rttov_opts % rt_ir % do_nlte_correction)
    end if

    !< Switch to enable IR aerosol calculations
    if (f_confOpts % has("RTTOV_addaerosl")) then
      call f_confOpts % get_or_die("RTTOV_addaerosl", self % rttov_opts % rt_ir % addaerosl)
    end if

    !< Switch to supply aerosol optical properties explicitly per channel
    if (f_confOpts % has("RTTOV_user_aer_opt_param")) then
      call f_confOpts % get_or_die("RTTOV_user_aer_opt_param", self % rttov_opts % rt_ir % user_aer_opt_param)
    end if

    !< Switch to enable IR cloudy calculations
    if (f_confOpts % has("RTTOV_addclouds")) then
      call f_confOpts % get_or_die("RTTOV_addclouds", self % rttov_opts % rt_ir % addclouds)
    end if

    !< Switch to supply cloud optical properties explicitly per channel
    if (f_confOpts % has("RTTOV_user_cld_opt_param")) then
      call f_confOpts % get_or_die("RTTOV_user_cld_opt_param", self % rttov_opts % rt_ir % user_cld_opt_param)
    end if

    !< Switch to supply grid-box average cloud concentration or cloud
    if (f_confOpts % has("RTTOV_grid_box_avg_cloud")) then
      call f_confOpts % get_or_die("RTTOV_grid_box_avg_cloud", self % rttov_opts % rt_ir % grid_box_avg_cloud)
    end if

    !!  concentration in cloudy fraction of each layer
    !< Ignore cloud streams with weights lower than this
    if (f_confOpts % has("RTTOV_cldstr_threshold")) then
      call f_confOpts % get_or_die("RTTOV_cldstr_threshold", self % rttov_opts % rt_ir % cldstr_threshold)
    end if

    !< Switch for simplified cloud stream option - USE WITH CAUTION
    if (f_confOpts % has("RTTOV_cldstr_simple")) then
      call f_confOpts % get_or_die("RTTOV_cldstr_simple", self % rttov_opts % rt_ir % cldstr_simple)
    end if

    !< Upper pressure limit for cldstr_simple option (hPa)
    if (f_confOpts % has("RTTOV_cldstr_low_cloud_top")) then
      call f_confOpts % get_or_die("RTTOV_cldstr_low_cloud_top", self % rttov_opts % rt_ir % cldstr_low_cloud_top)
    end if

    !< IR scattering model to use
    if (f_confOpts % has("RTTOV_ir_scatt_model")) then
      call f_confOpts % get_or_die("RTTOV_ir_scatt_model", self % rttov_opts % rt_ir % ir_scatt_model)
    end if

    !< VIS/NIR scattering model to use
    if (f_confOpts % has("RTTOV_vis_scatt_model")) then
      call f_confOpts % get_or_die("RTTOV_vis_scatt_model", self % rttov_opts % rt_ir % vis_scatt_model)
    end if

    !< Number of DOM streams, must be even and not less than 2
    if (f_confOpts % has("RTTOV_dom_nstreams")) then
      call f_confOpts % get_or_die("RTTOV_dom_nstreams", self % rttov_opts % rt_ir % dom_nstreams)
    end if

    !< Convergence criterion for termination of DOM azimuthal loop
    if (f_confOpts % has("RTTOV_dom_accuracy")) then
      call f_confOpts % get_or_die("RTTOV_dom_accuracy", self % rttov_opts % rt_ir % dom_accuracy)
    end if

    !< DOM ignores levels below this optical depth:
    if (f_confOpts % has("RTTOV_dom_opdep_threshold")) then
      call f_confOpts % get_or_die("RTTOV_dom_opdep_threshold", self % rttov_opts % rt_ir % dom_opdep_threshold)
    end if

    !< Switch to enable input of O3 profile
    if (f_confOpts % has("RTTOV_ozone_data")) then
      call f_confOpts % get_or_die("RTTOV_ozone_data", self % rttov_opts % rt_ir % ozone_data)
    end if

    !< Switch to enable input of CO2 profile
    if (f_confOpts % has("RTTOV_co2_data")) then
      call f_confOpts % get_or_die("RTTOV_co2_data", self % rttov_opts % rt_ir % co2_data)
    end if

    !< Switch to enable input of N2O profile
    if (f_confOpts % has("RTTOV_n2o_data")) then
      call f_confOpts % get_or_die("RTTOV_n2o_data", self % rttov_opts % rt_ir % n2o_data)
    end if

    !< Switch to enable input of CO profile
    if (f_confOpts % has("RTTOV_co_data")) then
      call f_confOpts % get_or_die("RTTOV_co_data", self % rttov_opts % rt_ir % co_data)
    end if

    !< Switch to enable input of CH4 profile
    if (f_confOpts % has("RTTOV_ch4_data")) then
      call f_confOpts % get_or_die("RTTOV_ch4_data", self % rttov_opts % rt_ir % ch4_data)
    end if

    !< Switch to enable input of SO2 profile
    if (f_confOpts % has("RTTOV_so2_data")) then
      call f_confOpts % get_or_die("RTTOV_so2_data", self % rttov_opts % rt_ir % so2_data)
    end if

    if (self % do_mw_scatt) then
      !Met Office default is true
      if ( f_ConfOpts % has("MW_Scatt_Use_TotalIce")) then
        call f_confOpts % get_or_die("MW_Scatt_Use_TotalIce", self % mw_scatt % use_totalice)
      end if

      !< Snow and rain input units are: False => kg/m2/s; True => kg/kg
      if ( f_ConfOpts % has("MW_Scatt_MMR_SnowRain")) then
        call f_confOpts % get_or_die("MW_Scatt_MMR_SnowRain", self % mw_scatt % mmr_snowrain)
      end if

      !< RT calculations in radiances instead of BT (default is false, recommended setting is true)
      if ( f_ConfOpts % has("MW_Scatt_lradiance")) then
        call f_confOpts % get_or_die("MW_Scatt_lradiance", self % mw_scatt % opts % lradiance)
      end if

      !< Manually supply effective cloud fraction (default is false)
      if ( f_ConfOpts % has("MW_Scatt_lusercfrac")) then
        call f_confOpts % get_or_die("MW_Scatt_lusercfrac", self % mw_scatt % opts % lusercfrac)
      end if

      !< Threshold for determining if scattering calculations will be performed (default is 0.05)
      if ( f_ConfOpts % has("MW_Scatt_cc_threshold")) then
        call f_confOpts % get_or_die("MW_Scatt_cc_threshold", self % mw_scatt % opts % cc_threshold)
      end if

      !< Hydrometeor TL/AD/K sensitivity includes the indirect effect
      !< through the effective cloud fraction (default is true, but false when PS_configuration is true)
      if ( f_ConfOpts % has("MW_Scatt_hydro_cfrac_tlad")) then
        call f_confOpts % get_or_die("MW_Scatt_hydro_cfrac_tlad", self % mw_scatt % opts % hydro_cfrac_tlad)
      end if

      !< !< Switch for hydrometeor TL/AD sensitivity in layers with zero hydrometeor concentration
      !< (default is false, but true when PS_configuration is true)
      if ( f_ConfOpts % has("MW_Scatt_zero_hydro_tlad")) then
        call f_confOpts % get_or_die("MW_Scatt_zero_hydro_tlad", self % mw_scatt % opts % zero_hydro_tlad)
      end if

    end if

    call rttov_print_opts(self % rttov_opts,lu = stderr)
    if (self % do_MW_scatt) call rttov_print_opts_scatt(self % mw_scatt % opts,lu = stderr)
  end subroutine ufo_rttov_set_options

  ! ------------------------------------------------------------------------------

  subroutine ufo_rttov_setup(self, f_confOpts)
    class(rttov_conf), intent(inout)      :: self
    type(fckit_configuration), intent(in) :: f_confOpts ! RTcontrol

    integer :: i_inst, rttov_errorstatus
    character(len=max_string) :: message

    include 'rttov_read_coefs.interface'
    include 'rttov_read_scattcoeffs.interface'

    rttov_errorstatus = 0

    if (.not. self%rttov_is_setup ) then

      ! --------------------------------------------------------------------------
      ! 1. Setup rttov options
      ! --------------------------------------------------------------------------
      call self % set_options(f_confOpts)

      ! --------------------------------------------------------------------------
      ! 2. Read coefficients
      ! --------------------------------------------------------------------------
      allocate(self % rttov_coef_array(self % nSensors))

      do i_inst = 1, self%nSensors
          
        call rttov_read_coefs(rttov_errorstatus, &       !out
                              self % rttov_coef_array(i_inst), & !inout
                              self % rttov_opts, &           !in
                              instrument = self % instrument_triplet(1:3,i_inst), &
                              path = self % COEFFICIENT_PATH)

        if (rttov_errorstatus /= errorstatus_success) then
          message = 'fatal error reading coefficients: ' // trim(self % coeffname(i_inst))
          call abor1_ftn(message)
        else
          message = 'successfully read RT coefficients: ' // trim(self % coeffname(i_inst))
          call fckit_log%info(message)
        end if
      end do

      ! Read MW scatt coefficients. Only one set of scattering coefficients may be read at the moment which
      ! is not typically an issue because we are only dealing with one instrument and they don't currently 
      ! have different scatt coefs
      if (self % do_mw_scatt) then 
        call rttov_read_scattcoeffs (rttov_errorstatus,        & ! out
                                     self % mw_scatt % opts,    & ! in
                                     self % rttov_coef_array(1), & ! in
                                     self % mw_scatt % coef,          & ! inout
                                     path = self % COEFFICIENT_PATH)

        if (rttov_errorstatus /= errorstatus_success) then
          message = 'fatal error reading compatible MWscatt coefficients for: ' // trim(self % coeffname(1))
          call abor1_ftn(message)
        else
          message = 'successfully read compatible MWscatt coefficients for: ' // trim(self % coeffname(1))
          call fckit_log%info(message)
        end if
      end if
      
      self % rttov_is_setup =.true.
    end if
      
    end subroutine ufo_rttov_setup

  ! ------------------------------------------------------------------------------
  !ufo_rttov_alloc is a wrapper for RTTOV12/13 allocation
  subroutine ufo_rttov_setup_rtprof(self,geovals,obss,conf,ob_info)

    use ufo_rttovonedvarcheck_ob_mod, only : ufo_rttovonedvarcheck_ob

    implicit none

    class(ufo_rttov_io), target,  intent(inout)             :: self
    type(ufo_geovals),            intent(in)                :: geovals
    type(c_ptr), value,           intent(in)                :: obss
    type(rttov_conf),             intent(in)                :: conf
    type(ufo_rttovonedvarcheck_ob), optional, intent(inout) :: ob_info

    ! Local variables
    type(rttov_profile), pointer       :: profiles(:)
    type(rttov_profile_cloud), pointer :: profiles_scatt(:)
    type(ufo_geoval), pointer          :: geoval
    type(datetime), allocatable        :: date_temp(:)

    integer                            :: jspec, isensor, ivar, iprof, ilev
    integer                            :: nlevels
    integer                            :: nprofiles

    character(MAXVARLEN)               :: varname
    character(len=max_string)          :: message

    real(kind_real), allocatable       :: TmpVar(:), windsp(:), p(:), ph(:)
    real(kind_real)                    :: s2m_t(1), s2m_p(1)

    logical                            :: variable_present

    integer                            :: skinTIndexFromObsSpace
    integer                            :: ctpIndexFromObsSpace
    integer                            :: ecaIndexFromObsSpace

    integer                            :: top_level, bottom_level, stride
    integer                            :: level_1000hPa, level_950hpa
    real(kind_real)                    :: NewT

    real(kind_real), allocatable       :: q_temp(:), clw_temp(:), ciw_temp(:), Qtotal(:), qsaturated(:)
    real(kind_real)                    :: scale_fac
    real(c_double)                     :: missing

    integer                            :: year, month, day, hour, minute, second

    integer, allocatable               :: sat_id(:)

    missing = missing_value(missing)

    profiles => self % profiles
    profiles_scatt => self % mw_scatt % profiles

    if(present(ob_info)) then
      self % nlocs_total = 1
    else
      self % nlocs_total = obsspace_get_nlocs(obss)
    end if
    if (self % nlocs_total == 0) return

    nprofiles = min(size(profiles), geovals%nlocs)

    allocate(self % sensor_index_array(nprofiles))
    self % sensor_index_array = -1

    ! Setup variables that need to be read from ObsSpace rather than GeoVaLs
    ! note emissivity is done via the `surface emissivity group` in the yaml
    skinTIndexFromObsSpace = 0 ! initialise
    ctpIndexFromObsSpace = 0 ! initialise
    ecaIndexFromObsSpace = 0 ! initialise
    do ivar = 1, size(conf % variablesFromObsSpace)
      if (index(conf % variablesFromObsSpace(ivar), "skinTemperature") > 0) then
        skinTIndexFromObsSpace = ivar
      else if (index(conf % variablesFromObsSpace(ivar), "pressureAtTopOfCloud") > 0) then
        ctpIndexFromObsSpace = ivar
      else if (index(conf % variablesFromObsSpace(ivar), "cloudAmount") > 0) then
        ecaIndexFromObsSpace = ivar
      else
        message = 'ERROR: ' // trim(conf % variablesFromObsSpace(ivar)) // &
                  ' not setup to be read from ObsSpace => Aborting'
        call fckit_exception % throw(message)
      end if
    end do

    ! store which sensor we will be using to process the observation. 
    ! Primarily for choosing an RTTOV coefficient    
    allocate(sat_id(nprofiles))
    if (present(ob_info)) then 
      sat_id = ob_info % satellite_identifier
    else
      if (obsspace_has(obss, "MetaData", "satelliteIdentifier")) then
        call obsspace_get_db(obss, "MetaData", "satelliteIdentifier", sat_id)
      else
        self % sensor_index_array = 1
      end if
    end if

    do isensor = 1, size(conf % wmo_id)
      where(sat_id == conf % wmo_id(isensor))
        ! If conf%nsensors is 1 then sensor_idx will always be 1.
        self % sensor_index_array = min(isensor, conf%nsensors)
      end where
    end do

    if (allocated(sat_id)) deallocate(sat_id)

    if (present(ob_info)) then
      call datetime_to_yyyymmddhhmmss(ob_info % date, year, month, day, hour, minute, second)
      profiles(1) % date = (/year, month, day/)
      profiles(1) % time = (/hour, minute, second/)
    else
      if (obsspace_has(obss, "MetaData", "dateTime")) then
        allocate(date_temp(obsspace_get_nlocs(obss)))
        call obsspace_get_db(obss, "MetaData", "dateTime", date_temp)
        do iprof=1, nprofiles
          call datetime_to_yyyymmddhhmmss(date_temp(iprof), year, month, day, hour, minute, second)
          profiles(iprof) % date = (/year, month, day/)
          profiles(iprof) % time = (/hour, minute, second/)
        end do
        deallocate(date_temp)
      else
        message = 'Warning: Optional input Date/Time not in database'
        call fckit_log%info(message)
      end if
    end if

    nlevels = size(profiles(1)%p)

    ! Assume that the pressure profile coming from the geovals is increasing in pressure (ToA->surface)...
    top_level = 1
    bottom_level = nlevels
    stride=1

    ! Get pressure levels and check they're up the right way (deprecated)
    if (ufo_vars_getindex(geovals%variables, var_prs) > 0) then
      call ufo_geovals_get_var(geovals, var_prs, geoval)

      allocate(p(nlevels))
      do iprof = 1, nprofiles
        p = geoval%vals(geoval%nval-nlevels+1:geoval%nval, iprof) * Pa_to_hPa ! for RTTOV
        if (iprof == 1) then
          if (p(1) > p(2)) then ! ...but be ready to switch. Assume if one profile is 'upside-down' then they all are
            top_level = nlevels
            bottom_level = 1
            stride = -1
          end if
        end if
        profiles(iprof)%p(top_level:bottom_level:stride) = p(:)
      end do
      deallocate(p)
    end if

! Get temperature
    varname = var_ts
    call ufo_geovals_get_var(geovals, varname, geoval)

    do iprof = 1, nprofiles
      profiles(iprof)%t(top_level:bottom_level:stride) = geoval%vals(:, iprof) ! K
    end do

! Get absorbers
    if (conf % RTTOV_GasUnitConv) then
      !gas_units = 0 is ppmv dry. Conversion will be done prior to use by RTTOV. Currently only scalar conversion performed.
      !this matches satrad conversion
      profiles(1:nprofiles)%gas_units = gas_unit_ppmvdry
    else
      !gas_units = 1 is kg/kg moist. Conversion will be done internally by RTTOV
      profiles(1:nprofiles)%gas_units = gas_unit_specconc
    end if

    do jspec = 1, conf%ngas

      scale_fac = conf%scale_fac(conf%absorber_id(jspec))

      select case (conf%Absorbers(jspec))
      case (var_mixr) ! mixr assumed to be in g/kg
        call ufo_geovals_get_var(geovals, conf%Absorbers(jspec), geoval)
        do iprof = 1, nProfiles
          profiles(iprof)%q(top_level:bottom_level:stride) = geoval%vals(:, iprof) * scale_fac * g_to_kg
        end do
      case (var_q)
        call ufo_geovals_get_var(geovals, conf%Absorbers(jspec), geoval)
        do iprof = 1, nProfiles
          profiles(iprof)%q(top_level:bottom_level:stride) = geoval%vals(:, iprof) * scale_fac
        end do
      case (var_oz)
        if (associated(profiles(1)%o3)) then
          if (.not. allocated(self % tc_ozone)) then 
            allocate(self % tc_ozone(nprofiles))
            self % tc_ozone = missing
          end if
          if (conf % RTTOV_scale_ref_ozone) then
            ! scale ozone by determining TCOzone based on empirically derived coefficients 
            ! and 70 hPa temperature. 
            ! self % tc_ozone populated in calculate_tc_ozone and used in scale_ozone
            if(present(ob_info)) then
              if(any(ob_info % background_ozone == missing)) then
                call self % calculate_tc_ozone(obss)
                call self % scale_ozone(conf)
                ob_info % background_ozone = self % profiles(1) % o3
              else
                self % profiles(1) % o3 = ob_info % background_ozone
              endif
            else ! ob_info not present
              call self % calculate_tc_ozone(obss)
              call self % scale_ozone(conf)
            endif
          else 
            call ufo_geovals_get_var(geovals, conf%Absorbers(jspec), geoval)
            do iprof = 1, nProfiles
              profiles(iprof)%o3(top_level:bottom_level:stride) = geoval%vals(:, iprof) * scale_fac
            end do
            call self % calculate_tc_ozone(obss)
          end if
          
          if(.not. present(ob_info)) then
            call obsspace_put_db(obss, "MetaData", "ozoneTotal", self % tc_ozone)
          end if
        end if
      case (var_co2)
        call ufo_geovals_get_var(geovals, conf%Absorbers(jspec), geoval)
        if (associated(profiles(1)%co2)) then
          do iprof = 1, nProfiles
            profiles(iprof)%co2(top_level:bottom_level:stride) = geoval%vals(:, iprof) * scale_fac 
          end do
        end if
      case (var_clw)
        call ufo_geovals_get_var(geovals, conf%Absorbers(jspec), geoval)
        if (conf % do_mw_scatt) then
          do iprof = 1, nProfiles
            profiles_scatt(iprof)%clw(top_level:bottom_level:stride) = geoval%vals(:, iprof)
          end do
        else if (conf % rttov_opts % rt_mw % clw_data) then
          if (associated(profiles(1)%clw)) then
            do iprof = 1, nprofiles
              profiles(iprof)%clw(top_level:bottom_level:stride) = geoval%vals(:, iprof) ! always kg/kg
            end do
          end if
        end if
      case (var_cli)
        call ufo_geovals_get_var(geovals, conf%Absorbers(jspec), geoval)
        if (conf % do_mw_scatt) then
          if (conf % mw_scatt % use_totalice) then
            do iprof = 1, nProfiles
              profiles_scatt(iprof)%totalice(top_level:bottom_level:stride) = geoval%vals(:, iprof)
            end do
          else
            do iprof = 1, nProfiles
              profiles_scatt(iprof)%ciw(top_level:bottom_level:stride) = geoval%vals(:, iprof)
            end do
          end if
        else
          message = 'ufo_rttov_setup_rtprof: Cloud Ice Water only supported for RTTOV-SCATT'
          call abor1_ftn(message)
        end if
      case default

      end select

    end do

! Get near-surface variables (s2m)

    varname = var_ps
    if (ufo_vars_getindex(geovals%variables, varname) > 0) then
      call ufo_geovals_get_var(geovals, varname, geoval)

      profiles(1:nprofiles)%s2m%p = geoval%vals(1,:) * Pa_to_hPa
    else
      message = 'No near-surface pressure. Using bottom pressure level'
      call fckit_log%info(message)

      do iprof = 1, nprofiles
        profiles(iprof)%s2m%p = profiles(iprof)%p(nlevels)
      enddo
    end if

    varname = var_sfc_t2m ! 2m temperature
    if (ufo_vars_getindex(geovals%variables, varname) > 0) then
      call ufo_geovals_get_var(geovals, varname, geoval) 
      profiles(1:nprofiles)%s2m%t = geoval%vals(1,1:nprofiles)
    else
      message = 'No near-surface temperature. Using bottom temperature level'
      call fckit_log%info(message)
      do iprof = 1, nprofiles
        profiles(iprof)%s2m%t = profiles(iprof)%t(nlevels)
      enddo
    end if

    varname = var_sfc_q2m ! 2m specific humidity
    if (ufo_vars_getindex(geovals%variables, varname) > 0) then
      call ufo_geovals_get_var(geovals, varname, geoval) ! lfric

      profiles(1:nprofiles)%s2m%q = geoval%vals(1,1:nprofiles) * conf%scale_fac(gas_id_watervapour)
    else
      message = 'No near-surface specific humidity. Using bottom q level'
      call fckit_log%info(message)

      do iprof = 1, nprofiles
        profiles(iprof)%s2m%q = profiles(iprof)%q(nlevels)
      enddo
    end if

    varname = var_sfc_u10 ! Eastward-wind in m/s 
    if (ufo_vars_getindex(geovals%variables, varname) > 0) then
      call ufo_geovals_get_var(geovals, varname, geoval)

      profiles(1:nprofiles)%s2m%u = geoval%vals(1,1:nprofiles)
      !assume if eastward then northward too

      varname = var_sfc_v10 ! Northward-wind in m/s 
      call ufo_geovals_get_var(geovals, varname, geoval)

      profiles(1:nprofiles)%s2m%v = geoval%vals(1,1:nprofiles)
    else !! use windspeed and direction instead
      allocate(windsp(nprofiles))
      call ufo_geovals_get_var(geovals, var_sfc_wspeed, geoval)

      windsp(1:nprofiles) = geoval%vals(1,1:nprofiles)

      call ufo_geovals_get_var(geovals, var_sfc_wdir, geoval)

      do iprof = 1, nprofiles
        profiles(iprof)%s2m%u             = windsp(iprof) * cos(geoval%vals(1, iprof) * deg2rad)
        profiles(iprof)%s2m%v             = windsp(iprof) * sin(geoval%vals(1, iprof) * deg2rad)
      enddo
      deallocate(windsp)
    end if

    allocate(TmpVar(self % nlocs_total))

    ! Get Skin (surface) temperature (K)
    if (skinTIndexFromObsSpace > 0) then
      call get_from_obsspace(conf % variablesFromObsSpace(skinTIndexFromObsSpace), &
                             obss, profiles(1:nprofiles)%skin%t)
    else
      varname = var_sfc_tskin
      call ufo_geovals_get_var(geovals, varname, geoval)
      profiles(1:nprofiles)%skin%t = geoval%vals(1,1:nprofiles)
    end if

    ! Setup cloud values
    ! presuure at the top of cloud - ctp
    ! pressureAtTopOfCloud is stored in the ObsSpace in Pa conversion needed as
    ! rttov needs hPa.
    if (ctpIndexFromObsSpace > 0) then
      call get_from_obsspace(conf % variablesFromObsSpace(ctpIndexFromObsSpace), &
                             obss, profiles(1:nprofiles) % ctp)
      profiles(1:nprofiles) % ctp = profiles(1:nprofiles) % ctp * pa_to_hpa
    else
      profiles(1:nprofiles) % ctp = 850.0_kind_real
    end if

    ! cloudAmount - eca
    if (ecaIndexFromObsSpace > 0) then
      call get_from_obsspace(conf % variablesFromObsSpace(ecaIndexFromObsSpace), &
                             obss, profiles(1:nprofiles) % cfraction)
    else
      profiles(1:nprofiles) % cfraction = zero
    end if

!RTTOV-Scatt profile setup

    ! The top half-level is the top of the atmosphere so it should be a nominal
    ! small value. The RTTOV convention is that the bottom half-level is the
    ! surface, but this differs from how the UM defines the bottom half-level.
    ! Let's use the RTTOV convention here.

    if (conf % do_mw_scatt) then  
      if (ufo_vars_getindex(geovals%variables, var_prsi) > 0 ) then
        call ufo_geovals_get_var(geovals, var_prsi, geoval)
        allocate(ph(nlevels+1))
        do iprof = 1, nprofiles
          ph = geoval%vals(geoval%nval-(nlevels+1)+1:geoval%nval, iprof) * Pa_to_hPa
          if (ph(1) > ph(2)) then !upside-down
            profiles_scatt(iprof) % ph(top_level+1:bottom_level:stride) = ph(:)
          else
            profiles_scatt(iprof) % ph(top_level:bottom_level+1:stride) = ph(:)
          end if
          
          profiles_scatt(iprof) % ph(nlevels+1) = profiles(iprof)%s2m%p
          profiles_scatt(iprof) % ph(1) = RTTOV_ToA ! hard-code top of atmosphere 0.0001 hPa
          
        end do
        deallocate(ph)
      else
        message = 'half-level pressures required for RTTOV-SCATT interface'
        call abor1_ftn(message)
      end if
      
      ! cloud fraction from var_cldfrac_vol TODO (IR update)

      ! The input cloud concentrations must be the layer grid-box-average 
      ! concentration (as opposed to the concentration within the cloudy fraction of each layer)
      call ufo_geovals_get_var(geovals, var_cldfrac_vol, geoval)
      do iprof = 1, nprofiles
        profiles_scatt(iprof) % cc(top_level:bottom_level:stride) = &
          geoval%vals(:, iprof)
      enddo
      
      ! The scattering code will use either ciw or totalice depending on the
      ! contents of the coefficient file so we need to assign both.
      ! solid precipitation and rain are currently not supported and are left at the initialised value of 0
      call ufo_geovals_get_var(geovals, var_cli, geoval)
      do iprof = 1, nprofiles
        if (conf % mw_scatt % use_totalice) then
          profiles_scatt(iprof) % totalice(top_level:bottom_level:stride) = &
            geoval%vals(:, iprof)
        else
          profiles_scatt(iprof) % ciw(top_level:bottom_level:stride) = &
            geoval%vals(:, iprof)
        end if
      enddo
    end if

! ---------------------------
! SatRad profile manipulation
! ---------------------------

    if(conf % SatRad_compatibility) then
      allocate(self % q_profile_reset(nProfiles, nlevels))
      self % q_profile_reset(:, :) = .false.
      allocate(self % clw_profile_reset(nProfiles, nlevels))
      self % clw_profile_reset(:, :) = .false.
      allocate(self % ciw_profile_reset(nProfiles, nlevels))
      self % ciw_profile_reset(:, :) = .false.
      do iprof = 1, nProfiles
        !----
        ! Reset low level temperatures over seaice and cold, low land as per Ops_SatRad_SetUpRTprofBg.F90
        ! N.B. I think this should be flagged so it's clear that the background has been modified
        !----
        
        if(profiles(iprof)%skin%surftype /= surftype_sea .and. &
           conf % UseColdSurfaceCheck) then
          if(profiles(iprof)%skin%t < 271.4_kind_real .and. &
            profiles(iprof)%s2m%p  > 950.0_kind_real) then

            level_1000hpa = minloc(abs(profiles(iprof)%p - 1000.0_kind_real),DIM=1)
            level_950hpa = minloc(abs(profiles(iprof)%p - 950.0_kind_real),DIM=1)

            NewT = profiles(iprof)%t(level_950hpa)
            if(profiles(iprof)%s2m%p > 1000.0_kind_real) &
              NewT = max(NewT,profiles(iprof)%t(level_1000hPa))
            NewT = min(NewT, 271.4_kind_real)

            profiles(iprof)%t(level_1000hPa) = max(profiles(iprof)%t(level_1000hPa), NewT)
            profiles(iprof)%s2m%t = max(profiles(iprof)%s2m%t, NewT)
            profiles(iprof)%skin%t = max(profiles(iprof)%skin%t, NewT)
          end if
        end if

        ! -----------------------------------------------
        ! Make sure q does not exceed saturation
        ! -----------------------------------------------
        if(conf % BoundQToSaturation) then
          allocate(qsaturated(nlevels))
          if (conf % UseRHwaterForQC) then
            call Ops_QsatWat (qsaturated(:),    & ! out
                              profiles(iprof) % t(:),  & ! in
                              profiles(iprof) % p(:) / Pa_to_hPa, & ! in convert hPa to Pa
                              nlevels)           ! in
          else
            call Ops_Qsat (qsaturated(:),    & ! out
                           profiles(iprof) % t(:),  & ! in
                           profiles(iprof) % p(:) / Pa_to_hPa, & ! in convert hPa to Pa
                           nlevels)           ! in
          end if

          qsaturated = qsaturated * conf%scale_fac(gas_id_watervapour)

          !qsaturated is assumed to be in kg/kg
          where (profiles(iprof)%q > qsaturated)
            profiles(iprof)%q = qsaturated
          end where
          deallocate(qsaturated)

          ! -----------------------------------------------
          ! Make sure q2m does not exceed saturation
          ! -----------------------------------------------
          allocate(qsaturated(1))
          s2m_t(1) = profiles(iprof)%s2m%t
          s2m_p(1) = profiles(iprof)%s2m%p / Pa_to_hPa
          if (conf % UseRHwaterForQC) then
            call Ops_QsatWat (qsaturated(1:1),  & ! out
                              s2m_t(1:1), & ! in
                              s2m_p(1:1), & ! in
                              1)            ! in
          else
            call Ops_Qsat (qsaturated(1:1),  & ! out
                           s2m_t(1:1), & ! in
                           s2m_p(1:1), & ! in
                           1)            ! in
          end if

          qsaturated(1) = qsaturated(1) * conf%scale_fac(gas_id_watervapour)

          if (profiles(iprof)%s2m%q > qsaturated(1)) profiles(iprof)%s2m%q = qsaturated(1)
          deallocate(qsaturated)
        end if
        
        if(conf % UseMinimumQ) then
          ! Constrain small values to min_q for q profile
          do ilev = 1, nlevels
            if (profiles(iprof) % q(ilev) < min_q * conf % scale_fac(gas_id_watervapour)) then
              profiles(iprof) % q(ilev) = min_q * conf%scale_fac(gas_id_watervapour)
              self % q_profile_reset(iprof, ilev) = .true.
            end if
          end do
          ! Constrain small values to min_q for surface q
          if(profiles(iprof)%s2m%q < min_q * conf%scale_fac(gas_id_watervapour)) profiles(iprof)%s2m%q = min_q * conf%scale_fac(gas_id_watervapour)
        end if

        if(conf % UseMinimumClw) then
          if (conf % do_mw_scatt) then
            ! Constrain small values to min_clw for clw profile
            do ilev = 1, nlevels
              if (profiles_scatt(iprof) % clw(ilev) < min_clw) then
                profiles_scatt(iprof) % clw(ilev) = min_clw
                self % clw_profile_reset(iprof, ilev) = .true.
              end if
            end do
          else if (conf % rttov_opts % rt_mw % clw_data) then
            ! Constrain small values to min_clw for clw profile
            do ilev = 1, nlevels
              if (profiles(iprof) % clw(ilev) < min_clw) then
                profiles(iprof) % clw(ilev) = min_clw
                self % clw_profile_reset(iprof, ilev) = .true.
              end if
            end do
          end if
        end if

        if(conf % UseMinimumCiw .and. conf % do_mw_scatt) then
          ! Cloud Ice Water only supported for RTTOV-SCATT
          if (conf % mw_scatt % use_totalice) then
            ! Constrain small values to min_ciw for totalice profile
            do ilev = 1, nlevels
              if (profiles_scatt(iprof) % totalice(ilev) < min_ciw) then
                profiles_scatt(iprof) % totalice(ilev) = min_ciw
                self % ciw_profile_reset(iprof, ilev) = .true.
              end if
            end do
          else
            ! Constrain small values to min_ciw for ciw profile
            do ilev = 1, nlevels
              if (profiles_scatt(iprof) % ciw(ilev) < min_ciw) then
                profiles_scatt(iprof) % ciw(ilev) = min_ciw
                self % ciw_profile_reset(iprof, ilev) = .true.
              end if
            end do
          end if
        end if

        if(conf % MWScattZeroJacPress > 0.0 .and. conf % do_mw_scatt) then
          ! set cloud jacobians to zero above MWScattZeroJacPress
          do ilev = 1, nlevels
            if (profiles(iprof) % p(ilev) < conf % MWScattZeroJacPress) then
              self % clw_profile_reset(iprof, ilev) = .true.
              self % ciw_profile_reset(iprof, ilev) = .true.
            end if
          end do
        end if

      enddo
    end if

    if (.not. allocated(self % ciw)) allocate(self % ciw(nlevels,nprofiles))

    if(conf % SplitQtotal) then
      allocate(Qtotal(nlevels), q_temp(nlevels), clw_temp(nlevels), ciw_temp(nlevels))
      do iprof = 1, nprofiles

        Qtotal(:) = profiles(iprof) % q(:) / conf%scale_fac(gas_id_watervapour) ! kg/kg
        Qtotal(:) = max(Qtotal(:), min_q)

        if (conf % do_mw_scatt) then
          if (conf % mw_scatt % use_totalice) then          
            Qtotal(:) = Qtotal(:) + &
                        profiles_scatt(iprof) % clw(:) + &
                        profiles_scatt(iprof) % totalice(:)
          else
            Qtotal(:) = Qtotal(:) + &
                        profiles_scatt(iprof) % clw(:) + &
                        profiles_scatt(iprof) % ciw(:)
          end if
        else if (conf % rttov_opts % rt_mw % clw_data) then
          Qtotal(:) = Qtotal(:) + profiles(iprof) % clw(:)
        end if

        ! generate first guess cloud and q based on the qtotal physics
        call Ops_SatRad_Qsplit (1,                     & ! in
          profiles(iprof) % p(:) / Pa_to_hPa,          & ! in convert hPa to Pa
          profiles(iprof) % t(:),                      & ! in
          Qtotal(:),                                   & ! in
          q_temp(:),                                   & ! out
          clw_temp(:),                                 & ! out
          ciw_temp(:),                                 & ! out
          UseQtSplitRain = conf % UseQtSplitRain)

        ! store the profile
        profiles(iprof) % q(:) = q_temp(:) * conf%scale_fac(gas_id_watervapour)

        if (conf % do_mw_scatt) then
          profiles_scatt(iprof) % clw(:) = clw_temp(:)
          if (conf % mw_scatt % use_totalice) then          
            profiles_scatt(iprof) % totalice(:) = ciw_temp(:)
          else
            profiles_scatt(iprof) % ciw(:) = ciw_temp(:)
          end if
        else if (conf % rttov_opts % rt_mw % clw_data) then
          profiles(iprof) % clw(:) = clw_temp(:)
        end if

        self % ciw(:,iprof) = ciw_temp(:)

      enddo

      deallocate(Qtotal, q_temp, clw_temp, ciw_temp)

    else
      ! ensure ciw is available to output in the RTTOV interface
      if (conf % do_mw_scatt) then
        do iprof = 1, nprofiles
          if (conf % mw_scatt % use_totalice) then
            self % ciw(:,iprof) = profiles_scatt(iprof) % totalice(:)
          else
            self % ciw(:,iprof) = profiles_scatt(iprof) % ciw(:)
          endif
        enddo
      end if
    end if

! Set geometry for RTTOV calculation, either from the supplied ob info (1dvar) or the obsspace db
    if(present(ob_info)) then

      profiles(1) % elevation = ob_info % elevation / 1000.0 ! m -> km
      profiles(1) % latitude = ob_info % latitude
      profiles(1) % longitude = ob_info % longitude
      if (ob_info % retrievecloud) then
        profiles(1) % ctp = ob_info % cloudtopp
        profiles(1) % cfraction = ob_info % cloudfrac
      end if

      nprofiles = 1
      nlevels = size(profiles(1) % p)

      profiles(1) % zenangle    = ob_info % sensor_zenith_angle
      profiles(1) % azangle     = ob_info % sensor_azimuth_angle
      profiles(1) % sunzenangle = ob_info % solar_zenith_angle
      profiles(1) % sunazangle  = ob_info % solar_azimuth_angle

      profiles(1) % skin % surftype = ob_info % surface_type

    else

!Set RT profile elevation (ob has priority, otherwise model height from geoval)
      if (obsspace_has(obss, "MetaData", "heightOfSurface")) then
        call obsspace_get_db(obss, "MetaData", "heightOfSurface", TmpVar)
        profiles(1:nprofiles)%elevation = TmpVar(1:nprofiles) * m_to_km !for RTTOV
        message = 'Using MetaData/surface_height for profile elevation'
      else if (ufo_vars_getindex(geovals%variables, "surface_altitude") > 0) then
        call ufo_geovals_get_var(geovals, "surface_altitude", geoval)
        profiles(1:nprofiles)%elevation = geoval%vals(1, 1:nprofiles) * m_to_km
        message = 'Using surf_altitude from GeoVaLs for profile elevation'
      else
        message = 'MetaData elevation not in database'
      end if
      call fckit_log%info(message)

      if (allocated(TmpVar)) deallocate(TmpVar)

!lat/lon
      variable_present = obsspace_has(obss, "MetaData", "latitude")
      if (variable_present) then
        call obsspace_get_db(obss, "MetaData", "latitude", profiles(1:nprofiles)%latitude)
      else
        message = 'Warning: Optional input MetaData/latitude not in database'
        call fckit_log%info(message)
      end if

      variable_present = obsspace_has(obss, "MetaData", "longitude")
      if (variable_present) then
        call obsspace_get_db(obss, "MetaData", "longitude", profiles(1:nprofiles)%longitude)
      else
        message = 'MetaData longitude not in database: check implicit filtering'
        call fckit_log%info(message)
      end if

!Set RTTOV viewing geometry

      ! sensor zenith - RTTOV convention 0-max (coef dependent). Nadir = 0 deg
      variable_present = obsspace_has(obss, "MetaData", "sensorZenithAngle")
      if (variable_present) then
        call obsspace_get_db(obss, "MetaData", "sensorZenithAngle", profiles(1:nprofiles)%zenangle)
      else
        message = 'ERROR: Mandatory input MetaData/sensorZenithAngle not in database. Aborting...'
        call abor1_ftn(message)
      end if

      ! sensor azimuth - convention is 0-360 deg. E=+90
      variable_present = obsspace_has(obss, "MetaData", "sensorAzimuthAngle")
      if (variable_present) then
        call obsspace_get_db(obss, "MetaData", "sensorAzimuthAngle", profiles(1:nprofiles)%azangle)
      else
        message = 'Warning: Optional input MetaData/sensorAzimuthAngle not in database: setting to zero for RTTOV'
        call fckit_log%info(message)
        profiles(1:nprofiles)%azangle = zero
      end if

      ! solar zenith
      variable_present = obsspace_has(obss, "MetaData", "solarZenithAngle")
      if (variable_present) then
        call obsspace_get_db(obss, "MetaData", "solarZenithAngle", profiles(1:nprofiles)%sunzenangle)
      else
        message = 'Warning: Optional input MetaData/solarZenithAngle not in database: setting to zero'
        call fckit_log%info(message)
        profiles(1:nprofiles)%sunzenangle = zero
      end if

      ! solar azimuth
      variable_present = obsspace_has(obss, "MetaData", "solarAzimuthAngle")
      if (variable_present) then
        call obsspace_get_db(obss, "MetaData", "solarAzimuthAngle", profiles(1:nprofiles)%sunazangle)
      else
        message = 'Warning: Optional input MetaData/solarAzimuthAngle not in database: setting to zero for RTTOV'
        call fckit_log%info(message)
        profiles(1:nprofiles)%sunazangle = zero
      end if

      ! RTTOV surface type
      variable_present = obsspace_has(obss, "MetaData", var_surf_type_rttov)
      if (variable_present) then
        call obsspace_get_db(obss, "MetaData", var_surf_type_rttov, profiles(1:nprofiles)%skin%surftype)
      else
        message = 'ERROR: Mandatory input MetaData/surfaceQualifier not in database. Aborting...'
        call abor1_ftn(message)
      end if
    end if

  end subroutine ufo_rttov_setup_rtprof

  subroutine ufo_rttov_check_rtprof(self, conf, iprof, errorstatus)
    implicit none
    
    class(ufo_rttov_io), target,  intent(inout) :: self
    type(rttov_conf),             intent(in)    :: conf
    integer,                      intent(in)    :: iprof
    integer,                      intent(out)   :: errorstatus

    character(10)             :: prof_str
    character(len=max_string) :: message
    integer                   :: sensor_idx

    include 'rttov_print_profile.interface'
    include 'rttov_user_profile_checkinput.interface'

    sensor_idx = self % sensor_index_array(iprof)

    errorstatus = errorstatus_success

    if (sensor_idx < 1) then

      errorstatus = errorstatus_fatal
      write(message, '(A,I4,A,I4,A)') &
      'Bad sensor index (', sensor_idx, ') for profile ', iprof, &
      ' which will not be processed and no further checking will be performed'
      call fckit_log%info(message)
    else
      if(conf % RTTOV_profile_checkinput) then
        call rttov_user_profile_checkinput(errorstatus, &
          conf % rttov_opts, &
          conf % rttov_coef_array(sensor_idx), &
          self % profiles(iprof))

        ! print erroneous profile to stderr
        if(errorstatus /= errorstatus_success .and. conf % debug) then
          write(prof_str,'(I0)') iprof
          self % profiles(iprof) % id = prof_str
          call rttov_print_profile(self % profiles(iprof), lu = stderr)
          write(message, '(A,I0)') 'Error in profile ', iprof
          call fckit_log%info(message)
        end if
        
        if ((conf % rttov_opts % rt_mw % fastem_version >= 3) .and. &
            (conf % rttov_coef_array(sensor_idx) % coef % id_sensor == sensor_id_mw) .and. &
            (self % profiles(iprof) % azangle < zero .or. & 
             self % profiles(iprof) % azangle > 360.0_kind_real)) then
          errorstatus = errorstatus_fatal
          write(message, '(A,I0,A)') 'Bad azimuth angle for requested FASTEM version >= 3 for profile ', &
                                     iprof, ' which will not be processed'
          call fckit_log%info(message)
        endif

        if (self % profiles(iprof) % longitude < -180.0_kind_real .or. & 
          self % profiles(iprof) % longitude > 360.0_kind_real) then
          errorstatus = errorstatus_fatal
          write(message, '(A,I0,A)') 'Bad longitude when using emissivity atlas for profile ', iprof, &
                                     ' which will not be processed'
          call fckit_log%info(message)
        endif
      endif
    endif

  end subroutine ufo_rttov_check_rtprof

  subroutine ufo_rttov_print_rtprof(self, conf, iprof)
    implicit none
    
    class(ufo_rttov_io), target,  intent(inout) :: self
    type(rttov_conf),             intent(in)    :: conf
    integer,                      intent(in)    :: iprof

    character(10) :: prof_str

    include 'rttov_print_profile.interface'
    include 'rttov_print_cld_profile.interface'

    write(*,*) 'profile ', iprof
    write(prof_str,'(i0)') iprof
    self % profiles(iprof) % id = prof_str
 
    call rttov_print_profile(self % profiles(iprof), lu = stdout)
    if (conf % do_MW_scatt) call rttov_print_cld_profile(self % mw_scatt % profiles(iprof), lu = stdout)

  end subroutine ufo_rttov_print_rtprof

  !ufo_rttov_alloc is a wrapper for RTTOV12/13 allocation
  subroutine ufo_rttov_alloc_direct(self, errorstatus, conf, nprofiles, nchannels, nlevels, init, asw)
    implicit none

    class(ufo_rttov_io), target, intent(inout) :: self
    type(rttov_conf),    intent(in)    :: conf

    integer,             intent(out)   :: errorstatus ! Return error status of RTTOV subroutine calls
    integer,             intent(in)    :: asw ! allocate switch
    logical, optional,   intent(in)    :: init ! initialise, default yes
    integer ,            intent(in)    :: nchannels
    integer ,            intent(in)    :: nprofiles
    integer ,            intent(in)    :: nlevels

    logical                            :: init1
    character(len=max_string)          :: message

    include 'rttov_alloc_direct.interface'
    include 'rttov_alloc_emis_ret_terms.interface'

    if (present(init)) then
      init1 = init
    else 
      init1 = .true.
    end if

    ! Allocate structures for rttov_direct
    !profiles are already allocated
    call rttov_alloc_direct(          &
      errorstatus,                    &
      asw,                            &  ! 1 => allocate
      nprofiles,                      &
      nchannels,                      &
      nlevels,                        &
      opts = conf % rttov_opts,           &
      coefs = conf % rttov_coef_array(1), &
      transmission = self % transmission, &
      radiance = self % radiance,         &
      calcemis = self % calcemis,         &
      emissivity = self % emissivity,     &
      init = init1)
    
    if (errorstatus /= errorstatus_success) then
      write(message,'(A, I6)') 'after rttov_alloc_direct error = ', errorstatus
      call abor1_ftn(message)
    end if
    
    if (conf % do_mw_scatt) then
      call rttov_alloc_emis_ret_terms (errorstatus,            & ! in
                                       nchannels, & ! in
                                       self % mw_scatt % emis_retrieval,   & ! in
                                       asw)      ! in
    end if

    if(asw == 1) then
      ! set frequency array indices for mietable coeffs for each channel used by RTTOV-SCATT 
      ! (calls rttov_scatt_setupindex)

      if (conf % do_mw_scatt) then
        allocate(self % mw_scatt % freq_indices(nchannels))
        call self % mw_scatt % set_frequencies(conf % rttov_coef_array(1), nprofiles, nchannels)
      end if

      !additional initialisation not done inside RTTOV
      self % calcemis = .false.
      self % emissivity % emis_in = -one
      self % emissivity % emis_out = -one
    else
      if (conf % do_mw_scatt) then
        if (associated(self % mw_scatt % freq_indices)) deallocate(self % mw_scatt % freq_indices)
      end if
      if (allocated(self % ciw))                      deallocate(self % ciw)
      if (allocated(self % tc_ozone))                 deallocate(self % tc_ozone)
      if (allocated(self % q_profile_reset))          deallocate(self % q_profile_reset)
      if (allocated(self % clw_profile_reset))          deallocate(self % clw_profile_reset)
      if (allocated(self % ciw_profile_reset))          deallocate(self % ciw_profile_reset)
    end if

  end subroutine ufo_rttov_alloc_direct

  !ufo_rttov_alloc is a wrapper for RTTOV12/13 allocation
  subroutine ufo_rttov_alloc_k(self, errorstatus, conf, nprofiles, nchannels, nlevels, init, asw)
    implicit none

    class(ufo_rttov_io), target, intent(inout) :: self
    type(rttov_conf),    intent(in)    :: conf

    integer,             intent(out)   :: errorstatus ! Return error status of RTTOV subroutine calls
    integer,             intent(in)    :: asw ! allocate switch
    logical, optional,   intent(in)    :: init ! initialise, default yes
    integer ,            intent(in)    :: nchannels
    integer ,            intent(in)    :: nprofiles
    integer ,            intent(in)    :: nlevels

    character(len=max_string)          :: message
    logical                            :: init1

    include 'rttov_alloc_k.interface'
   
    if (present(init)) then
      init1 = init
    else
      init1 = .true.
    end if

    call rttov_alloc_k(               &
      errorstatus,                    &
      asw,                            &  ! 1 => allocate
      nprofiles,                      &
      nchannels,                      &
      nlevels,                        &
      opts = conf % rttov_opts,           &
      coefs = conf % rttov_coef_array(1), &
      transmission_k = self % transmission_k, &
      radiance_k = self % radiance_k,         &
      emissivity_k = self % emissivity_k,     &
      init = init1)
    
      if (errorstatus /= errorstatus_success) then
        write(message,'(A, I6)') 'after rttov_alloc_k error = ', errorstatus
        call abor1_ftn(message)
      end if

      if (asw > 0 ) then
   
        ! Inintialize the K-matrix INPUT so that the results are dTb/dx
        ! -------------------------------------------------------------
        
        self % emissivity_k(:) % emis_out = zero
        self % emissivity_k(:) % emis_in = zero
        self % emissivity(:) % emis_out = zero
        self % radiance_k % bt(:) = one
        self % radiance_k % total(:) = one

      end if

      if (errorstatus /= errorstatus_success) then
        write(message,'(A, I6)') 'after rttov_alloc_k error = ', errorstatus
        call abor1_ftn(message)
      end if
   
  end subroutine ufo_rttov_alloc_k

  !ufo_rttov_alloc is a wrapper for RTTOV12/13 allocation
  subroutine ufo_rttov_alloc_profiles(self, errorstatus, conf, nprofiles, nlevels, init, asw)
    implicit none

    class(ufo_rttov_io), target, intent(inout) :: self
    type(rttov_conf),    intent(in)    :: conf

    integer,             intent(out)   :: errorstatus ! Return error status of RTTOV subroutine calls
    integer,             intent(in)    :: asw ! allocate switch
    logical, optional,   intent(in)    :: init ! initialise, default yes
    integer ,            intent(in)    :: nprofiles
    integer ,            intent(in)    :: nlevels
    integer                            :: iprof

    character(len=max_string)          :: message
    logical                            :: init1

    include 'rttov_alloc_prof.interface'
    include 'rttov_alloc_scatt_prof.interface'

    if (present(init)) then
      init1 = init
    else 
      init1 = .true.
    end if

    if (asw == 1) allocate (self % profiles(nprofiles))

    ! Allocate structures for rttov_direct
    call rttov_alloc_prof(          &
      errorstatus,                    &
      nprofiles,                      &
      self % profiles,                &
      nlevels,                        &
      conf % rttov_opts,              &
      asw,                            &
      coefs = conf % rttov_coef_array(1), &
      init=init1)

    if (errorstatus /= errorstatus_success) then
      write(message,'(A, I6)') 'after rttov_alloc_profiles error = ', errorstatus
      call abor1_ftn(message)
    end if
    
    if (conf % do_mw_scatt) then
      ! Allocate the RTTOV-SCATT cloud profiles structure

      if (asw == 1) allocate(self % mw_scatt % profiles(nprofiles))

      call rttov_alloc_scatt_prof(    &
        errorstatus,                    &
        nprofiles,                      &
        self % mw_scatt % profiles(:),  &
        nlevels,                        &
        conf % mw_scatt % use_totalice, &    ! false => separate ciw and snow; true => totalice
        asw,                            &    ! 1 => allocate
        init = .true.,             &
        mmr_snowrain = conf % mw_scatt % mmr_snowrain)  ! snow/rain input units: false => kg/m2/s; true => kg/kg
    end if

    ! Hard code fixed profile inputs which the interface does not currently permit to vary
    if (asw == 1) then

      self % profiles(1:nprofiles) % skin % watertype = watertype_ocean_water  
      self % profiles(1:nprofiles) % s2m % wfetc = 100000.0_kind_real ! wind fetch (m) taken
                                                             ! from users guide
      self % profiles(1:nprofiles) % skin % salinity = 35.0_kind_real

      !Default fastem parameters
      do iprof = 1, size(self % profiles)
        self % profiles(iprof) % skin % fastem = [3.0_kind_real,  &
                                                  5.0_kind_real,  &
                                                  15.0_kind_real, &
                                                  0.1_kind_real,  &
                                                  0.3_kind_real]
      end do

      ! set surftype to invalid so that it must be set explicitly otherwise an error will occur
      self % profiles(:) % skin % surftype = -1
    else
      if (allocated(self % profiles))                 deallocate(self % profiles)
      if (conf % do_mw_scatt) then
        if (associated(self % mw_scatt % profiles))   deallocate(self % mw_scatt % profiles)
      end if
      if (allocated(self % ciw))                      deallocate(self % ciw)
      if (allocated(self % tc_ozone))                 deallocate(self % tc_ozone)
      if (allocated(self % q_profile_reset))          deallocate(self % q_profile_reset)
      if (allocated(self % clw_profile_reset))          deallocate(self % clw_profile_reset)
      if (allocated(self % ciw_profile_reset))          deallocate(self % ciw_profile_reset)
    end if

  end subroutine ufo_rttov_alloc_profiles

  !ufo_rttov_alloc is a wrapper for RTTOV12/13 allocation
  subroutine ufo_rttov_alloc_profiles_k(self, errorstatus, conf, nprofiles, nlevels, init, asw)
    implicit none

    class(ufo_rttov_io), target, intent(inout) :: self
    type(rttov_conf),    intent(in)    :: conf

    integer,             intent(out)   :: errorstatus ! Return error status of RTTOV subroutine calls
    integer,             intent(in)    :: asw ! allocate switch
    logical, optional,   intent(in)    :: init ! initialise, default yes
    integer ,            intent(in)    :: nprofiles
    integer ,            intent(in)    :: nlevels

    character(len=max_string)          :: message
    logical                            :: init1

    include 'rttov_alloc_prof.interface'
    include 'rttov_alloc_scatt_prof.interface'

    if (present(init)) then
      init1 = init
    else 
      init1 = .true.
    end if

    if (asw == 1) then
      allocate (self % profiles_k(nprofiles))
    end if

    call rttov_alloc_prof(          &
      errorstatus,                    &
      nprofiles,                      &
      self % profiles_k,              &
      nlevels,                        &
      conf % rttov_opts,              &
      asw,                            &
      coefs = conf % rttov_coef_array(1), &
      init=init1)

    if (conf % do_mw_scatt) then
      ! Allocate the RTTOV-SCATT cloud profiles structure

      if (asw == 1) allocate(self % mw_scatt % profiles_k(nprofiles))

      call rttov_alloc_scatt_prof(       &
        errorstatus,                     &
        nprofiles,                       &
        self % mw_scatt % profiles_k(:), &
        nlevels,                         &
        conf % mw_scatt % use_totalice,  &    ! false => separate ciw and snow; true => totalice
        asw,                             &    ! 1 => allocate
        init = init1,                    &
        mmr_snowrain = conf % mw_scatt % mmr_snowrain)  ! snow/rain input units: false => kg/m2/s; true => kg/kg
    endif

    if (errorstatus /= errorstatus_success) then
      write(message,'(A, I6)') 'after rttov_alloc_profiles error = ', errorstatus
      call abor1_ftn(message)
    end if
    
    !don't deallocate profiles_k! That's the trajectory.

  end subroutine ufo_rttov_alloc_profiles_k

  subroutine ufo_rttov_zero_k(self, conf, reset_profiles_k)
    implicit none

    class(ufo_rttov_io), target, intent(inout) :: self
    type(rttov_conf), intent(in) :: conf
    logical, intent(in)          :: reset_profiles_k

    include 'rttov_init_prof.interface'
    include 'rttov_init_rad.interface'
    include 'rttov_init_scatt_prof.interface'
    include 'rttov_init_transmission.interface'

    if (reset_profiles_k) then
      call rttov_init_prof(self % profiles_k)
      if(conf % Do_MW_Scatt) then
        call rttov_init_scatt_prof(self % mw_scatt % profiles_k)
      end if
    end if

    call rttov_init_rad(self % radiance_k)
    call rttov_init_transmission(self % transmission_k)

! manually reset emissivity as there's no RTTOV routine to do it
    self % emissivity_k(:) % emis_in = zero
    self % emissivity_k(:) % emis_out = zero
    self % emissivity(:) % emis_out = zero

! Set radiance jacobian scaling (should equal one) 
    self % radiance_k % bt(:) = one
    self % radiance_k % total(:) = one

  end subroutine ufo_rttov_zero_k

  subroutine ufo_rttov_init_default_emissivity(self, conf, prof_list)
    implicit none

    class(ufo_rttov_io), intent(inout) :: self
    type(rttov_conf),    intent(in)    :: conf
    integer,             intent(in)    :: prof_list(:,:)

    integer                            :: prof, all_prof_index, iprof
    integer                            :: start_chan, end_chan

!Emissivity and calcemis are only set for used channels.
!emissivity is already initialised to zero (so RTTOV doesn't complain)
!cycle through the list of good profiles

    start_chan = 0
    end_chan = 0

    do iprof=1, size(prof_list, dim=1)
      if (prof_list(iprof,1) > 0) then
        prof = prof_list(iprof,1)
        all_prof_index = prof_list(iprof,2)
        start_chan = end_chan + 1
        end_chan = end_chan + self % nchan_inst
      else
        cycle
      end if

      if (self % profiles(all_prof_index) % skin % surftype == surftype_sea) then

        ! Calculate by RTTOV
        self % emissivity(start_chan:end_chan) % emis_in = 0.0_kind_real
        self % calcemis(start_chan:end_chan) = .true.

      else

        if (conf % rttov_coef_array(1) % coef % id_sensor == sensor_id_mw) then

          if (self % profiles(all_prof_index) % skin % surftype == surftype_land) then
            self % emissivity(start_chan:end_chan) % emis_in = 0.95_kind_real
          elseif (self % profiles(all_prof_index) % skin % surftype == surftype_seaice) then
            self % emissivity(start_chan:end_chan) % emis_in = 0.92_kind_real
          end if

        elseif ( conf % rttov_coef_array(1) % coef % id_sensor == sensor_id_ir .or. &
                 conf % rttov_coef_array(1) % coef % id_sensor == sensor_id_hi) then
          self % emissivity(start_chan:end_chan) % emis_in = 0.98_kind_real

        end if

        self % calcemis(start_chan:end_chan) = .false.

      end if
    enddo

  end subroutine ufo_rttov_init_default_emissivity

  subroutine ufo_rttov_set_defaults(self, default_opts_set)
    implicit none
    
    class(rttov_conf), intent(inout) :: self
    character(10),     intent(in)    :: default_opts_set

    integer                          :: PS_Number
    logical                          :: PS_configuration
    character(len=max_string)        :: message

    message = 'Setting RTTOV default options to ' // trim(default_opts_set)
    call fckit_log%info(message)

    ! Get PS number if it exists
    if(default_opts_set(1:4) == 'UKMO') then
      PS_configuration = .true.
      read(default_opts_set(8:9),*) PS_Number

      write(message,'(A, I3)') 'Setting RTTOV default options for PS', PS_Number
      call fckit_log%info(message)
    else
      PS_configuration = .false.
      PS_Number = -1
    end if

!Set RTTOV 12.3 defaults explicitly
    self % rttov_opts % config % apply_reg_limits        = .false. !< Switch to restrict input profiles to coef training limits
    self % rttov_opts % config % verbose                 = .true.  !< Switch for verbose output
    self % rttov_opts % config % do_checkinput           = .true.  !< Switch to apply internal profile checking
    self % rttov_opts % config % fix_hgpl                = .false. !< Switch to apply fix to match 2m p with elevation in geometry calculations

    self % rttov_opts % rt_all % addrefrac               = .false. !< Switch to enable atmospheric refraction
    self % rttov_opts % rt_all % switchrad               = .false. !< Switch for input units in AD/K models
    self % rttov_opts % rt_all % use_q2m                 = .true.  !< Switch to enable use of 2m q variable
    self % rttov_opts % rt_all % do_lambertian           = .false. !< Switch for setting Lambertian reflection (IR and MW)
    self % rttov_opts % rt_all % lambertian_fixed_angle  = .true.  !< Switch for fixed/parameterised effective angle for Lambertian option
    self % rttov_opts % rt_all % plane_parallel          = .false. !< Switch to ignore atmospheric curvature           
    self % rttov_opts % rt_all % rad_down_lin_tau        = .true.  !< Linear-in-tau or layer-mean for downwelling radiances    
    self % rttov_opts % rt_all % dtau_test               = .true.  !< Switch to apply dtau test in transmit/integrate calculations

    self % rttov_opts % rt_ir % solar_sea_brdf_model     = 1       !< Solar sea BRDF model (1-2)
    self % rttov_opts % rt_ir % ir_sea_emis_model        = 2       !< IR sea emissivity model (1-2)
    self % rttov_opts % rt_ir % addsolar                 = .false. !< Switch to enable solar simulations
    self % rttov_opts % rt_ir % rayleigh_single_scatt    = .true.  !< Switch to enable Rayleigh single-scattering for VIS/NIR channels
    self % rttov_opts % rt_ir % do_nlte_correction       = .false. !< Switch to enable NLTE bias correction
    self % rttov_opts % rt_ir % addaerosl                = .false. !< Switch to enable IR aerosol calculations
    self % rttov_opts % rt_ir % user_aer_opt_param       = .false. !< Switch to supply aerosol optical properties explicitly per channel
    self % rttov_opts % rt_ir % addclouds                = .false. !< Switch to enable IR cloudy calculations
    self % rttov_opts % rt_ir % user_cld_opt_param       = .false. !< Switch to supply cloud optical properties explicitly per channel
    self % rttov_opts % rt_ir % grid_box_avg_cloud       = .false. !< Switch to supply grid-box average cloud concentration or cloud
    
    self % rttov_opts % rt_ir % cldstr_threshold         = -one           !< Ignore cloud streams with weights lower than this
    self % rttov_opts % rt_ir % cldstr_simple            = .false.        !< Switch for simplified cloud stream option - USE WITH CAUTION
    self % rttov_opts % rt_ir % cldstr_low_cloud_top     = 750._kind_real !< Upper pressure limit for cldstr_simple option (hPa)
    self % rttov_opts % rt_ir % ir_scatt_model           = ir_scatt_chou  !< IR scattering model to use
    self % rttov_opts % rt_ir % vis_scatt_model          = vis_scatt_dom  !< VIS/NIR scattering model to use
    self % rttov_opts % rt_ir % dom_nstreams             = 8              !< Number of DOM streams, must be even and not less than 2
    self % rttov_opts % rt_ir % dom_accuracy             = 0._kind_real   !< Convergence criterion for termination of DOM azimuthal loop
    self % rttov_opts % rt_ir % dom_opdep_threshold      = 0._kind_real   !< DOM ignores levels below this optical depth:
    
    self % rttov_opts % rt_ir % ozone_data               = .false.       !< Switch to enable input of O3 profile
    self % rttov_opts % rt_ir % co2_data                 = .false.       !< Switch to enable input of CO2 profile
    self % rttov_opts % rt_ir % n2o_data                 = .false.       !< Switch to enable input of N2O profile
    self % rttov_opts % rt_ir % co_data                  = .false.       !< Switch to enable input of CO profile
    self % rttov_opts % rt_ir % ch4_data                 = .false.       !< Switch to enable input of CH4 profile
    self % rttov_opts % rt_ir % so2_data                 = .false.       !< Switch to enable input of SO2 profile                      
    
    self % rttov_opts % rt_ir % pc % addpc               = .false. !< Switch to enable PC-RTTOV
    self % rttov_opts % rt_ir % pc % ipcbnd              = -1      !< PC spectral band
    self % rttov_opts % rt_ir % pc % ipcreg              = -1      !< PC predictor channel set
    self % rttov_opts % rt_ir % pc % npcscores           = -1      !< Number of PC scores to compute
    self % rttov_opts % rt_ir % pc % addradrec           = .false. !< Switch for calculation of reconstructed radiances
    
    !> MW-only radiative transfer options
    self % rttov_opts % rt_mw % fastem_version           = 6       !< FASTEM version (0-6); 0 => TESSEM2
    self % rttov_opts % rt_mw % supply_foam_fraction     = .false. !< Supply a foam fraction to FASTEM
    self % rttov_opts % rt_mw % clw_data                 = .false. !< Switch to enable input of cloud liquid water profile
    self % rttov_opts % rt_mw % clw_scheme               = mw_clw_scheme_liebe !< MW CLW scheme: 1 => Liebe, 2 => Rosenkranz, 3 => TKC
    self % rttov_opts % rt_mw % clw_calc_on_coef_lev     = .true.  !< Apply MW CLW calculations on coef/user levels (true/false resp.)
    self % rttov_opts % rt_mw % clw_cloud_top            = 322     !< Lower pressure limit for MW CLW calculations (hPa)
    self % rttov_opts % rt_mw % apply_band_correction    = .true.  !< Apply band-correction for Planck radiance and BT calculations
    
    self % rttov_opts % interpolation % addinterp        = .false.       !< Switch to enable RTTOV interpolator
    self % rttov_opts % interpolation % interp_mode      = interp_rochon !< Interpolation mode 1 (valid options 1-5, see user guide)
    self % rttov_opts % interpolation % lgradp           = .false.       !< Switch to make pressure an active variable in TL/AD/K models
    self % rttov_opts % interpolation % spacetop         = .true.        !< Switch to assume space boundary at top-most input pressure level
    self % rttov_opts % interpolation % reg_limit_extrap = .false.       !< Switch to extrapolate input profiles using regression limits
    
    !> HTFRTC options structure
    self % rttov_opts % htfrtc_opts % htfrtc             = .false. !< Switch to use htfrtc
    self % rttov_opts % htfrtc_opts % n_pc_in            = -1   !< Number of principal components to be used
    self % rttov_opts % htfrtc_opts % reconstruct        = .false. !< Switch to select reconstructed radiances
    self % rttov_opts % htfrtc_opts % simple_cloud       = .false. !< Calculate simple cloud
    self % rttov_opts % htfrtc_opts % overcast           = .false. !< Calculate overcast cloud on all levels
    
    !> RTTOV-Scatt specific options 
    if (self % do_mw_scatt) then
      self % mw_scatt % use_totalice        = .false.        !< False => separate ice and snow; True => total ice      
      self % mw_scatt % mmr_snowrain        = .true.         !< Snow and rain input units are: False => kg/m2/s; True => kg/kg

      !Copy all corresponding RTTOV options to RTTOV_SCATT options structure
      self % mw_scatt % opts % config                       = self % rttov_opts % config
      self % mw_scatt % opts % use_q2m                      = self % rttov_opts % rt_all % use_q2m
      self % mw_scatt % opts % fastem_version               = self % rttov_opts % rt_mw % fastem_version
      self % mw_scatt % opts % supply_foam_fraction         = self % rttov_opts % rt_mw % supply_foam_fraction
      self % mw_scatt % opts % apply_band_correction        = self % rttov_opts % rt_mw % apply_band_correction
      self % mw_scatt % opts % interp_mode                  = self % rttov_opts % interpolation % interp_mode
      self % mw_scatt % opts % reg_limit_extrap             = self % rttov_opts % interpolation % reg_limit_extrap
      self % mw_scatt % opts % lgradp                       = self % rttov_opts % interpolation % lgradp
      self % mw_scatt % opts % dtau_test                    = self % rttov_opts % rt_all % dtau_test
      self % mw_scatt % opts % rad_down_lin_tau             = self % rttov_opts % rt_all % rad_down_lin_tau

      self % mw_scatt % opts % lradiance                    = .false.        !< RT calculations in radiances instead of BT (default is false, recommended setting is true)
      self % mw_scatt % opts % lusercfrac                   = .false.        !< Manually supply effective cloud fraction (default is false)
      self % mw_scatt % opts % cc_threshold                 = 0.05_kind_real !< Threshold for determining if scattering calculations will be performed (default is 0.05)
      self % mw_scatt % opts % hydro_cfrac_tlad             = .true.         !< Switch for hydrometeor TL/AD sensitivity to effective cfrac (default is true).
      self % mw_scatt % opts % zero_hydro_tlad              = .false.        !< Switch for hydrometeor TL/AD sensitivity in layers with zero hydrometeor concentration (default is false).
    end if

    if (PS_configuration) then

      ! Set RTTOV options that different from default and are true for all MetO configurations up to PS45
      if (cmp_strings(default_opts_set(1:4), 'UKMO')) then
        self % rttov_opts % config % verbose                 = .false. ! true if (ProcessMode > VerboseMode .OR. RTTOV_Verbosity > 0)
        self % rttov_opts % config % do_checkinput           = .false. ! we will use the more thorough and verbose user_checkinput
      
        self % rttov_opts % rt_all % switchrad               = .true. 
        self % rttov_opts % rt_all % use_q2m                 = .false.
      
        self % rttov_opts % rt_ir % grid_box_avg_cloud       = .true. ! Assume grid-box average for cloud concentrations
        self % rttov_opts % rt_ir % ozone_data               = .true.  ! Set to true for allocation purposes
      
        self % rttov_opts % rt_mw % clw_data                 = .true.  ! Set to true for allocation purposes
      
        self % rttov_opts % interpolation % addinterp        = .true.  ! Allow interpolation of input profile
        self % rttov_opts % interpolation % interp_mode      = interp_rochon_wfn ! Set interpolation method (4 for all insts at PS44)
      end if
    
      !RTTOV options that are different from RTTOV defaults at and before PS44
      if (PS_Number <= 44) then
        self % rttov_opts % rt_ir % ir_sea_emis_model        = 1 ! Use SSIREM
        
        self % rttov_opts % rt_mw % fastem_version           = 2 ! no support for Fastem-bug so use caution
        
        self % rttov_opts % interpolation % spacetop         = .false. 
      end if
      
      !RTTOV options that are different from RTTOV and PS43 defaults at PS44
      if (PS_Number >= 44) then
        self % rttov_opts % config % apply_reg_limits        = .true.  
        self % rttov_opts % config % fix_hgpl                = .true.  ! This is an RTTOV 13 default
        
        self % rttov_opts % rt_all % dtau_test               = .false. ! This is an RTTOV 13 default
        self % rttov_opts % rt_all % rad_down_lin_tau        = .false. ! This is the recommended setting
      
        self % rttov_opts % rt_mw % clw_calc_on_coef_lev     = .false. ! This is an RTTOV 13 default
      end if
    
      !RTTOV options that are different from RTTOV and PS44 defaults at PS45       
      if (PS_Number == 45) then
        self % rttov_opts % rt_all % addrefrac               = .true. ! This is an RTTOV 13 default
        
        self % rttov_opts % rt_mw % clw_scheme               = mw_clw_scheme_rosenkranz ! This is an RTTOV 13 default
        
        self % rttov_opts % interpolation % reg_limit_extrap = .true. ! This is an RTTOV 13 default
        ! This is not an RTTOV change but the code no longer exists at PS45
        self % UseColdSurfaceCheck = .false. 
      end if

      if (self % do_mw_scatt) then
        self % mw_scatt % use_totalice = .true. ! Met Office default is true
        self % mw_scatt % opts % hydro_cfrac_tlad = .false.
        self % mw_scatt % opts % zero_hydro_tlad  = .true. 

        ! Need resetting following ps configurations above
        self % mw_scatt % opts % config                = self % rttov_opts % config
        self % mw_scatt % opts % use_q2m               = self % rttov_opts % rt_all % use_q2m
        self % mw_scatt % opts % fastem_version        = self % rttov_opts % rt_mw % fastem_version
        self % mw_scatt % opts % supply_foam_fraction  = self % rttov_opts % rt_mw % supply_foam_fraction
        self % mw_scatt % opts % apply_band_correction = self % rttov_opts % rt_mw % apply_band_correction
        self % mw_scatt % opts % interp_mode           = self % rttov_opts % interpolation % interp_mode
        self % mw_scatt % opts % reg_limit_extrap      = self % rttov_opts % interpolation % reg_limit_extrap
        self % mw_scatt % opts % lgradp                = self % rttov_opts % interpolation % lgradp
        self % mw_scatt % opts % dtau_test             = self % rttov_opts % rt_all % dtau_test
        self % mw_scatt % opts % rad_down_lin_tau      = self % rttov_opts % rt_all % rad_down_lin_tau

        ! prior to PS45 there was a bug where these options were not set correctly for RTTOV-SCATT
        ! and were different from RTTOV settings
        if (PS_Number < 45) then
          self % mw_scatt % opts % config % fix_hgpl = .false.
          self % mw_scatt % opts % rad_down_lin_tau  = .true.
          self % mw_scatt % opts % dtau_test         = .true.
        end if
      end if

    end if

  end subroutine ufo_rttov_set_defaults

  subroutine ufo_rttov_populate_hofxdiags(self, RTProf, chanprof, conf, prof_start, hofxdiags)

    class(rttov_hofxdiags), intent(inout) :: self
    type(ufo_rttov_io),   intent(in)    :: RTProf
    type(rttov_chanprof), intent(in)    :: chanprof(:)
    type(rttov_conf),     intent(in)    :: conf
    integer,              intent(in)    :: prof_start
    type(ufo_geovals),    intent(inout) :: hofxdiags    !non-h(x) diagnostics

    integer                       :: jvar, prof, ichan
    integer                       :: coefindex, chan
    integer                       :: nchanprof, nlevels, nprofiles
    integer                       :: rttov_errorstatus
    real(kind_real), allocatable  :: od_level(:), wfunc(:), tstore(:), bt_overcast(:)
    real(kind_real)               :: planck1, planck2, ff_bco, ff_bcs
    real(c_double)                :: missing
    character(len=max_string)     :: message

    include 'rttov_calc_weighting_fn.interface'

    allocate(od_level(size(RTProf % transmission%tau_levels(:,1))))
    allocate(wfunc(size(RTProf % transmission%tau_levels(:,1))))
    allocate(bt_overcast(size(RTProf % radiance % overcast(:,1))))
    allocate(tstore(size(RTProf % radiance % overcast(:,1))))

    missing = missing_value(missing)

    nchanprof = size(chanprof)
    nlevels = size(RTProf % profiles(1) % p)
    nprofiles = RTProf % nlocs_total

    do jvar = 1, hofxdiags%nvar
      if (len(trim(hofxdiags%variables(jvar))) < 1) cycle

      !============================================
      ! Diagnostics used for QC and bias correction
      !============================================

      if (cmp_strings(self % xstr_diags(jvar), "")) then
        ! forward h(x) diags
        select case(trim(self % ystr_diags(jvar)))

          ! variable: optical_thickness_of_atmosphere_layer_CH
          ! variable: transmittances_of_atmosphere_layer_CH
          ! variable: weightingfunction_of_atmosphere_layer_CH
          ! variable: mass_content_of_cloud_ice_in_atmosphere_layer
          ! variable: brightness_temperature_overcast_of_atmosphere_layer_CH
        case (var_opt_depth, var_lvl_transmit, var_lvl_weightfunc, var_cli, var_tb_overcast)

          hofxdiags%geovals(jvar)%nval = nlevels
          if(.not. allocated(hofxdiags%geovals(jvar)%vals)) then
            allocate(hofxdiags%geovals(jvar)%vals(hofxdiags%geovals(jvar)%nval,nprofiles))
            hofxdiags%geovals(jvar)%vals = missing
          end if

          ! get channel/profile
          do ichan = 1, nchanprof
            ! The chanprof contains the channel indexes in the coefficient file.
            ! If cut down coefficients are used (e.g. the coeffs contain channels 6,54,75,...) then
            ! chanprof(ichan)%chan won't contain the instrument channel numbers but the index (1,2,3,...).
            ! The correct instrument channel numbers are stored in ff_ori_chn.
            coefindex = chanprof(ichan)%chan
            if (coefindex < 1) cycle
            chan = conf % rttov_coef_array(1) % coef % ff_ori_chn(coefindex)
            prof = prof_start + chanprof(ichan)%prof - 1

            if(chan == self % ch_diags(jvar)) then
              ! if profile not skipped
              if(cmp_strings(self % ystr_diags(jvar), var_cli)) then
                hofxdiags%geovals(jvar)%vals(:,prof) = RTProf % ciw(:,prof)
              else if(cmp_strings(self % ystr_diags(jvar), var_opt_depth)) then
                od_level(:) = log(RTProf % transmission%tau_levels(:,ichan)) !level->TOA transmittances -> od
                hofxdiags%geovals(jvar)%vals(:,prof) = od_level(1:nlevels-1) - od_level(2:nlevels) ! defined +ve 
              else if (cmp_strings(self % ystr_diags(jvar), var_lvl_transmit)) then
                hofxdiags%geovals(jvar)%vals(:,prof) = RTProf % transmission % tau_levels(1:nlevels-1,ichan) - &
                                                       RTProf % transmission%tau_levels(2:,ichan)
              else if (cmp_strings(self % ystr_diags(jvar), var_lvl_weightfunc)) then
                od_level(:) = log(RTProf % transmission%tau_levels(:,ichan)) !level->TOA transmittances -> od
                call rttov_calc_weighting_fn(rttov_errorstatus, RTProf % profiles(prof)%p, od_level(:), &
                  hofxdiags%geovals(jvar)%vals(:,prof))
              else if (cmp_strings(self % ystr_diags(jvar), var_tb_overcast)) then
                planck1 = conf % rttov_coef_array(1) % coef % planck1(coefindex)
                planck2 = conf % rttov_coef_array(1) % coef % planck2(coefindex)
                ff_bco = conf % rttov_coef_array(1) % coef % ff_bco(coefindex)
                ff_bcs = conf % rttov_coef_array(1) % coef % ff_bcs(coefindex)

                tstore(:) = planck2 / log (1.0 + planck1 / RTProf % radiance % overcast(:, ichan))
                bt_overcast(:) = (tstore(:) - ff_bco) / ff_bcs

                ! The overcast BT is output on layers but the output required is on levels.  To match OPS this
                ! is mapped so the nearest surface level corresponds to the bottom layer.  This leaves the highest
                ! altitude level without a value and isothermal behaviour is assummed hence the copy. e.g.
                ! levels(2:70) = layers(1:69) - not a very good assumption.
                ! levels(1) = layer(1)
                ! where levels are top of the atmosphere to the surface.
                ! MCC - something to look at in the future
                hofxdiags%geovals(jvar)%vals(2:,prof) = bt_overcast(:)
                hofxdiags%geovals(jvar)%vals(1,prof) = bt_overcast(1)
              end if
            end if
          enddo

          ! variable: toa_outgoing_radiance_per_unit_wavenumber_CH [mW / (m^2 sr cm^-1)] (nval=1)
          ! variable: brightness_temperature_assuming_clear_sky_CH
          ! variable: pressure_level_at_peak_of_weightingfunction_CH
          ! variable: toa_total_transmittance_CH
          ! variable: surface_emissivity_CH
        case (var_radiance, var_tb_clr, var_tb, var_sfc_emiss, var_pmaxlev_weightfunc, var_total_transmit)
          ! always returned
          hofxdiags%geovals(jvar)%nval = 1
          if(.not. allocated(hofxdiags%geovals(jvar)%vals)) then
            allocate(hofxdiags%geovals(jvar)%vals(hofxdiags%geovals(jvar)%nval,nprofiles))
            hofxdiags%geovals(jvar)%vals = missing
          end if

          do ichan = 1, nchanprof
            coefindex = chanprof(ichan)%chan
            if (coefindex < 1) cycle
            chan = conf % rttov_coef_array(1) % coef % ff_ori_chn(coefindex)
            prof = prof_start + chanprof(ichan)%prof - 1

            if(chan == self % ch_diags(jvar)) then
              if(cmp_strings(self % ystr_diags(jvar), var_radiance)) then
                hofxdiags%geovals(jvar)%vals(1,prof) = RTProf % radiance % total(ichan)
              else if(cmp_strings(self % ystr_diags(jvar), var_tb_clr)) then
                hofxdiags%geovals(jvar)%vals(1,prof) = RTProf % radiance % bt_clear(ichan)
              else if(cmp_strings(self % ystr_diags(jvar), var_tb)) then
                hofxdiags%geovals(jvar)%vals(1,prof) = RTProf % radiance % bt(ichan)
              else if(cmp_strings(self % ystr_diags(jvar), var_pmaxlev_weightfunc)) then
                call rttov_calc_weighting_fn(rttov_errorstatus, RTProf % profiles(prof)%p, od_level(:), &
                  Wfunc(:))
                hofxdiags%geovals(jvar)%vals(1,prof) = maxloc(Wfunc(:), DIM=1) ! scalar not array(1)
              else if(cmp_strings(self % ystr_diags(jvar), var_total_transmit)) then
                if (conf % do_mw_scatt) then 
                  hofxdiags%geovals(jvar)%vals(1,prof) = RTProf % mw_scatt % emis_retrieval % tau_clr(ichan)
                else
                  hofxdiags%geovals(jvar)%vals(1,prof) = RTProf % transmission % tau_total(ichan)
                end if
              else if(cmp_strings(self % ystr_diags(jvar), var_sfc_emiss)) then
                hofxdiags%geovals(jvar)%vals(1,prof) = RTProf % emissivity(ichan) % emis_out
              end if
            end if
          end do

        case default
          ! not a supported obsdiag but we allocate and initialise here anyway for use later on
          hofxdiags%geovals(jvar)%nval = 1
          if(.not. allocated(hofxdiags%geovals(jvar)%vals)) then
            allocate(hofxdiags%geovals(jvar)%vals(hofxdiags%geovals(jvar)%nval,nprofiles))
            hofxdiags%geovals(jvar)%vals = missing
          end if

          if(conf % debug) then
            write(message,*) 'ufo_radiancerttov_simobs: ObsDiagnostic is unsupported but allocating anyway, ', &
                             trim(hofxdiags%variables(jvar)), shape(hofxdiags%geovals(jvar)%vals)
            call fckit_log%info(message)
          end if

        end select

      else if (cmp_strings(self % ystr_diags(jvar), var_tb)) then
        ! var_tb jacobians
        select case (trim(self % xstr_diags(jvar)))

        case (var_ts,var_mixr,var_q,var_clw,var_cli)

          hofxdiags%geovals(jvar)%nval = nlevels
          if(.not. allocated(hofxdiags%geovals(jvar)%vals)) then
            allocate(hofxdiags%geovals(jvar)%vals(hofxdiags%geovals(jvar)%nval,nprofiles))
            hofxdiags%geovals(jvar)%vals = missing
          end if

          do ichan = 1, nchanprof
            coefindex = chanprof(ichan)%chan
            if (coefindex < 1) cycle
            chan = conf % rttov_coef_array(1) % coef % ff_ori_chn(coefindex)
            prof = prof_start + chanprof(ichan)%prof - 1

            if(chan == self % ch_diags(jvar)) then
              if(self % xstr_diags(jvar) == var_ts) then
                hofxdiags%geovals(jvar)%vals(:,prof) = &
                  RTProf % profiles_k(ichan) % t(:)
              else if(self % xstr_diags(jvar) == var_mixr) then
                hofxdiags%geovals(jvar)%vals(:,prof) = &
                  RTProf % profiles_k(ichan) % q(:) * conf%scale_fac(gas_id_watervapour) / g_to_kg
              else if(self % xstr_diags(jvar) == var_q) then
                hofxdiags%geovals(jvar)%vals(:,prof) = &
                  RTProf % profiles_k(ichan) % q(:) * conf%scale_fac(gas_id_watervapour)
              else if(self % xstr_diags(jvar) == var_clw) then !clw
                if (conf % do_mw_scatt) then
                  hofxdiags%geovals(jvar)%vals(:,prof) = &
                    RTProf % mw_scatt % profiles_k(ichan) % clw(:)
                else if (conf % rttov_opts % rt_mw % clw_data) then
                  hofxdiags%geovals(jvar)%vals(:,prof) = &
                    RTProf % profiles_k(ichan) % clw(:)
                end if
              else if(self % xstr_diags(jvar) == var_cli) then
                if (conf % do_mw_scatt) then
                  if (conf % mw_scatt % use_totalice) then
                    hofxdiags%geovals(jvar)%vals(:,prof) = &
                      RTProf % mw_scatt % profiles_k(ichan) % totalice(:)
                  else
                    hofxdiags%geovals(jvar)%vals(:,prof) = &
                      RTProf % mw_scatt % profiles_k(ichan) % ciw(:)
                  endif
                else
                  hofxdiags%geovals(jvar)%vals(:,prof) = zero
                  if (conf % debug) then
                    message = 'ufo_radiancerttov_simobs: Cloud Ice Water only supported for RTTOV-SCATT'
                    call fckit_log%info(message)
                  end if
                end if
              end if
            end if
          enddo

        case (var_sfc_t2m, var_sfc_tskin, var_sfc_emiss, var_sfc_q2m, var_ps, var_sfc_u10, var_sfc_v10, &
              "cloud_top_pressure", "cloud_fraction")
          hofxdiags%geovals(jvar)%nval = 1
          if(.not. allocated(hofxdiags%geovals(jvar)%vals)) then
            allocate(hofxdiags%geovals(jvar)%vals(hofxdiags%geovals(jvar)%nval,nprofiles))
            hofxdiags%geovals(jvar)%vals = missing
          end if

          do ichan = 1, nchanprof
            coefindex = chanprof(ichan)%chan
            if(coefindex < 1) cycle
            chan = conf % rttov_coef_array(1) % coef % ff_ori_chn(coefindex)
            prof = prof_start + chanprof(ichan)%prof - 1

            if(chan == self % ch_diags(jvar)) then
              if(self % xstr_diags(jvar) == var_sfc_tskin) then
                hofxdiags%geovals(jvar)%vals(1,prof) = &
                  RTProf % profiles_k(ichan) % skin % t
              else if (self % xstr_diags(jvar) == var_sfc_t2m) then
                hofxdiags%geovals(jvar)%vals(1,prof) = &
                  RTProf % profiles_k(ichan) % s2m % t
              else if (self % xstr_diags(jvar) == var_ps) then
                hofxdiags%geovals(jvar)%vals(1,prof) = &
                  RTProf % profiles_k(ichan) % s2m % p
              else if (self % xstr_diags(jvar) == var_sfc_q2m) then
                hofxdiags%geovals(jvar)%vals(1,prof) = &
                  RTProf % profiles_k(ichan) % s2m % q * conf%scale_fac(gas_id_watervapour)
              else if (self % xstr_diags(jvar) == var_sfc_u10) then
                hofxdiags%geovals(jvar)%vals(1,prof) = &
                  RTProf % profiles_k(ichan) % s2m % u
              else if (self % xstr_diags(jvar) == var_sfc_v10) then
                hofxdiags%geovals(jvar)%vals(1,prof) = &
                  RTProf % profiles_k(ichan) % s2m % v
              else if (self % xstr_diags(jvar) == "cloud_top_pressure") then
                hofxdiags%geovals(jvar)%vals(1,prof) = &
                  RTProf % profiles_k(ichan) % ctp
              else if (self % xstr_diags(jvar) == "cloud_fraction") then
                hofxdiags%geovals(jvar)%vals(1,prof) = &
                  RTProf % profiles_k(ichan) % cfraction
              else if (self % xstr_diags(jvar) == var_sfc_emiss) then
                hofxdiags%geovals(jvar)%vals(1,prof) = &
                  RTProf % emissivity_k(ichan) % emis_in
              end if
            end if
          end do

        case default
          if (conf % debug) then
            message = 'ufo_radiancerttov_simobs: Jacobian ObsDiagnostic is unsupported, ' // trim(hofxdiags%variables(jvar))
            call fckit_log%info(message)
          end if  
        end select
      else
        if (conf % debug) then
          message = 'ufo_radiancerttov_simobs: ObsDiagnostic is not recognised, ' // trim(hofxdiags%variables(jvar))
          call fckit_log%info(message)
        end if
      end if

    enddo

    deallocate(od_level, wfunc, tstore, bt_overcast)
  end subroutine ufo_rttov_populate_hofxdiags

  subroutine ufo_rttov_parse_hofxdiags(self, hofxdiags, jacobian_needed)
    implicit none

    class(rttov_hofxdiags), intent(inout)  :: self
    type(ufo_geovals),      intent(in)     :: hofxdiags    !non-h(x) diagnostics
    logical,                intent(out)    :: jacobian_needed
    character(10), parameter  :: jacobianstr = "_jacobian_"

    integer                   :: str_pos(4)
    character(len=maxvarlen)  :: varstr
    integer                   :: jvar
    character(len=max_string) :: err_msg

    jacobian_needed = .false.

     if(hofxdiags%nvar > 0) then
       if (.not. allocated(self % ystr_diags)) allocate (self % ystr_diags(hofxdiags%nvar))
       if (.not. allocated(self % xstr_diags)) allocate (self % xstr_diags(hofxdiags%nvar))
       if (.not. allocated(self % ch_diags)) allocate (self % ch_diags(hofxdiags%nvar))

       self % ch_diags = -9999

       do jvar = 1, hofxdiags%nvar
         varstr = hofxdiags%variables(jvar)

         str_pos(4) = len_trim(varstr)
         if (str_pos(4) < 1) cycle
         str_pos(3) = index(varstr,"_",back=.true.)        !final "_" before channel
         read(varstr(str_pos(3)+1:str_pos(4)),*, err=999) self % ch_diags(jvar)
999      str_pos(1) = index(varstr,jacobianstr) - 1        !position before jacobianstr
         if (str_pos(1) == 0) then
           err_msg = 'parse_hofxdiags: _jacobian_ must be preceded by dependent variable in config: ' // &
                     trim(hofxdiags%variables(jvar))
           call abor1_ftn(err_msg)
         else if (str_pos(1) > 0) then
           !Diagnostic is a Jacobian member (dy/dx)
           self % ystr_diags(jvar) = varstr(1:str_pos(1))
           str_pos(2) = str_pos(1) + len(jacobianstr) + 1 !begin self % xstr_diags
           jacobian_needed = .true.
           str_pos(4) = str_pos(3) - str_pos(2)
           self % xstr_diags(jvar)(1:str_pos(4)) = varstr(str_pos(2):str_pos(3)-1)
           self % xstr_diags(jvar)(str_pos(4)+1:) = ""
         else !null
           !Diagnostic is a dependent variable (y)

           self % xstr_diags(jvar) = ""
           self % ystr_diags(jvar)(1:str_pos(3)-1) = varstr(1:str_pos(3)-1)
           self % ystr_diags(jvar)(str_pos(3):) = ""
           if (self % ch_diags(jvar) < 0) self % ystr_diags(jvar) = varstr
         end if
       end do
     end if

  end subroutine ufo_rttov_parse_hofxdiags

  subroutine ufo_rttov_hofxdiags_delete(self)
    implicit none

    class(rttov_hofxdiags), intent(inout) :: self

    if (allocated(self % ystr_diags)) deallocate (self % ystr_diags)
    if (allocated(self % xstr_diags)) deallocate (self % xstr_diags)
    if (allocated(self % ch_diags)) deallocate (self % ch_diags)
  end subroutine ufo_rttov_hofxdiags_delete

  subroutine set_freq_indices(self, rttov_coeffs, nprofiles, nchannels)
    implicit none
    
    class(mw_scatt_io), intent(inout) :: self  
    type(rttov_coefs), intent(in)     :: rttov_coeffs
    integer, intent(in)               :: nprofiles
    integer, intent(in)               :: nchannels
    
    integer :: nchans_inst
    type(rttov_chanprof), allocatable :: chanprof_dummy(:)

    include 'rttov_scatt_setupindex.interface'

    nchans_inst = rttov_coeffs % coef % fmv_chn

    allocate (chanprof_dummy(nchannels))
    call rttov_scatt_setupindex (nprofiles,            & ! in
                                 nchans_inst,          & ! in
                                 rttov_coeffs,         & ! in
                                 nchannels,            & ! in
                                 chanprof_dummy(:),    & ! out
                                 self % freq_indices(:)) ! out

    deallocate (chanprof_dummy)
  end subroutine set_freq_indices

  subroutine rttov_read_emissivity_from_obsspace(obss, surface_emissivity_group, channels, sfc_emiss)
    implicit none

    type(c_ptr), value, intent(in)  :: obss
    character(len=*), intent(in)    :: surface_emissivity_group
    integer, intent(in)             :: channels(:)
    real(kind_real), intent(out)    :: sfc_emiss(:,:)

    logical                         :: variable_present
    character(len=200)              :: var
    character(len=max_string)       :: message
    integer                         :: ichan

    variable_present = obsspace_has(obss, trim(surface_emissivity_group), trim("emissivity"))
    if (variable_present) then
      do ichan = 1, size(channels)
        ! Read in from the db
        write(var,'(A11,I0)') "emissivity_", channels(ichan)
        call obsspace_get_db(obss, trim(surface_emissivity_group), trim(var), sfc_emiss(ichan,:))
      end do
    else
      message = 'Surface emissivity group provided but not found in the database => aborting'
      call abor1_ftn(message)
    end if

  end subroutine rttov_read_emissivity_from_obsspace

  subroutine ufo_rttov_scale_ozone(self, conf)
    implicit none

    class(ufo_rttov_io), intent(inout) :: self
    type(rttov_conf), intent(in)       :: conf

    integer                            :: inst, i, iprof
    integer                            :: errorstatus
    character(len=max_string)          :: message
    real(c_double)                     :: missing

    include 'rttov_scale_ref_gas_prof.interface'

    missing = missing_value(missing)
    inst = -1

    do i = 1, size(conf % rttov_coef_array)
      if(conf % rttov_coef_array(i) % coef % nozone > 0) inst = i
    end do

    if (inst < 0) then
      message = 'Error: No ozone reference profile in coefficient. Aborting'
      call abor1_ftn(message)
    end if
    
    if (.not. allocated(self % tc_ozone) .or. &
      (allocated(self % tc_ozone) .and. any(self % tc_ozone == missing))) then
      message = 'Error: Total column ozone not calculated. Aborting'
      call abor1_ftn(message)
    end if

    !----
    ! Determine Ozone coefficients
    !----

    do iprof = 1, size(self % profiles)

      !RTTOV will scale the Ozone profile internally
      call rttov_scale_ref_gas_prof(errorstatus, conf % rttov_coef_array(inst), &
        self % profiles(iprof:iprof), &
        o3_col_int_du = self % tc_ozone(iprof), &
        satrad_compatibility = .true., &
        logp = .true.)

      if (errorstatus /= errorstatus_success) then
        message = 'Error in rttov_scale_ref_gas_prof. Aborting'
        call abor1_ftn(message)
      end if
    end do

  end subroutine ufo_rttov_scale_ozone

  subroutine ufo_rttov_calculate_tc_ozone(self,obss)
    implicit none

    class(ufo_rttov_io), target, intent(inout) :: self
    type(c_ptr), value, intent(in)             :: obss

    type(datetime)                     :: win_start, win_end

    type(rttov_profile), pointer       :: profiles(:)
    integer                            :: month, day    ! satrad month for ozone calculation
    integer                            :: dummy, iprof
    integer                            :: p70hpa        ! index of pressure level closest to 70 hPa
    real(kind_real)                    :: t70hpa        ! T at pressure level closest to 70 hPa
    character(len=max_string)          :: message
    integer, parameter                 :: rk = kind_real

    real(kind_real), parameter :: Ozone_c1(12) = &
      (/ -478.1_rk, -485.6_rk, -696.6_rk, -923.1_rk, -859.3_rk, -615.8_rk, &
         -431.1_rk, -449.6_rk, -681.9_rk, -780.4_rk, -528.3_rk, -424.4_rk /)

    real(kind_real), parameter :: Ozone_c2(12) = &
      (/ 3.7322_rk, 3.7532_rk, 4.7797_rk, 5.8444_rk, 5.5484_rk, 4.3689_rk, &
         3.4740_rk, 3.5425_rk, 4.6397_rk, 5.1078_rk, 3.9105_rk, 3.4179_rk /)

    if (.not. allocated(self % tc_ozone)) then
      message = 'Error: Total column ozone not allocated. Aborting'
      call abor1_ftn(message)
    end if

    ! Get the cycle endpoint time which we use to determine the date
    ! as this will always give the OPS equivalent month even in the case that window is shifted back
    ! by one second.
    ! SatRad changeover happens in the 'middle' of the month. This is for
    ! consistency with old code which had updated coefficients supplied on the 3rd
    ! Tuesday of each month in line with operational change practice and we persist it for now
    call obsspace_get_window(obss, win_start, win_end)
    call datetime_to_yyyymmddhhmmss(win_end, dummy, month, day, dummy, dummy, dummy)

    if( day <= 15) then
      month = month - 1
      if (month == 0) month = 12
    end if

    profiles => self % profiles

    do iprof = 1, size(profiles)

      ! calculate t at closest pressure level to 70hPa. This method is less efficent than the binary
      ! search from OPS but it's a very small array. Perhaps we can reimplement this in ufo utils?
      p70hpa = MINLOC(ABS(profiles(iprof) % p - 70.0_rk), dim=1)
      t70hpa = profiles(iprof) % t(p70hpa)

      !For calculating total column ozone from 70 hPa temperature, where the ozone
      !amount in Dobson units is given by
      !o3 = c1 + c2 * t70hPa
      self % tc_ozone(iprof) = Ozone_c1(month) + Ozone_c2(month) * t70hpa
    end do

  end subroutine ufo_rttov_calculate_tc_ozone

  subroutine get_from_obsspace(name, obss, outarray)
    implicit none
    character(len=*), intent(in)   :: name
    type(c_ptr), value, intent(in) :: obss
    real(kind_real), intent(out)   :: outarray(:)

    integer                        :: groupindex
    character(len=MAXVARLEN)       :: groupname, varname
    character(len=max_string)      :: message
    logical                        :: variable_present

    groupindex = index(name, "/")
    groupname = name(1:groupindex-1)
    varname = name(groupindex+1:)
    variable_present = obsspace_has(obss, trim(groupname), trim(varname))
    if (variable_present) then
      call obsspace_get_db(obss, trim(groupname), trim(varname), outarray(:))
    else
      message = 'Requested variable ' // trim(name) // ' not in ObsSpace => Aborting'
      call fckit_exception % throw(message)
    end if

  end subroutine get_from_obsspace

end module ufo_radiancerttov_utils_mod
