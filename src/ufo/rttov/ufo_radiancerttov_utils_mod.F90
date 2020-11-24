! (C) Copyright 2018 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

!> Fortran module to provide code shared between nonlinear and tlm/adm radiance calculations

module ufo_radiancerttov_utils_mod

  use fckit_configuration_module, only: fckit_configuration
  use fckit_log_module, only : fckit_log
  use iso_c_binding
  use kinds
  use missing_values_mod

  use rttov_types, only : rttov_options, rttov_profile, rttov_coefs, &
    rttov_radiance, rttov_transmission, rttov_emissivity, &
    rttov_chanprof

  use rttov_const ! gas_ids

  use ufo_vars_mod
  use ufo_geovals_mod, only: ufo_geovals, ufo_geoval, ufo_geovals_get_var
  use ufo_basis_mod, only: ufo_basis
  use ufo_utils_mod, only: Ops_SatRad_Qsplit
  use obsspace_mod

  implicit none
  private

  public rttov_conf
  public rttov_conf_setup
  public rttov_conf_delete
  public load_atm_data_rttov
  public load_geom_data_rttov
  public parse_hofxdiags
  public populate_hofxdiags
  public ufo_rttov_skip_profiles

  integer, parameter, public            :: max_string=800
  integer, parameter, public            :: MAXVARIN = 50

  !DARFIX this should go somewhere else

  character(len=maxvarlen), public      :: varin_temp(MAXVARIN)
  character(len=max_string), public     :: message

  integer, public                       :: nvars_in
  integer, public                       :: rttov_errorstatus

  ! ystr_diags contains the name of the rttov variable for ouput e.g.
  ! transmitance or optical depth.  For jacobian output var_tb
  character(len=maxvarlen), allocatable :: ystr_diags(:)
  ! xstr_diags contains the model variable for jacobian output e.g.
  ! var_t or empty if jacobian output is not required
  character(len=maxvarlen), allocatable :: xstr_diags(:)
  integer, allocatable                  :: ch_diags(:)

  real(c_double)                        :: missing

  character(len=maxvarlen), dimension(21), target :: varin_default_crtm = &
    (/var_ts, var_prs, var_prsi,                                &
    var_sfc_wfrac, var_sfc_lfrac, var_sfc_ifrac, var_sfc_sfrac, &
    var_sfc_wtmp,  var_sfc_ltmp,  var_sfc_itmp,  var_sfc_stmp,  &
    var_sfc_vegfrac, var_sfc_wspeed, var_sfc_wdir, var_sfc_lai, &
    var_sfc_soilm, var_sfc_soilt, var_sfc_landtyp,              &
    var_sfc_vegtyp, var_sfc_soiltyp, var_sfc_sdepth/)

  !var_ps
  character(len=maxvarlen), dimension(11), target :: varin_default_satrad = &
    (/var_prs, var_ts, var_q, var_sfc_t2m, & 
    var_u, var_v, var_sfc_p2m, var_sfc_q2m, &
    var_sfc_tskin, var_water_type_rttov, var_surf_type_rttov /)

  character(len=maxvarlen), pointer, public :: varin_default(:)

  ! copy of ABSORBER_ID_NAME defined in rttov_const
  character(len=*), parameter :: &
    RTTOV_Absorbers(ngases_max+2) = &
    [gas_name(1:ngases_max),'CLW', &
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
    gas_id_so2, 99, 999]

  character(len=MAXVARLEN), parameter :: null_str = ''

  !DARFIX: need to get correct names (correct units for RTTOV) in ufo_vars_mod
  character(len=MAXVARLEN), parameter :: &
    UFO_Absorbers(ngases_max+2) = &
    [ null_str, var_mixr, var_oz, null_str, var_co2, 'mole_fraction_of_nitrous_oxide_in_air', &
    'mole_fraction_of_carbon_monoxide_in_air', 'mole_fraction_of_methane_in_air', &
    'mole_fraction_of_sulfur_dioxide_in_air', var_clw, var_cli]

  integer, public :: nchan_inst ! number of channels being simulated (may be less than full instrument)
  integer, public :: nchan_sim  ! total number of 'obs' = nprofiles * nchannels
  integer, public :: nlocs_total ! nprofiles (including skipped)
  logical, public :: debug
!Common counters
  integer :: iprof

  type, public :: ufo_rttov_io
    logical, pointer                :: calcemis(:)     ! Flag to indicate calculation of emissivity within RTTOV

    type(rttov_emissivity), pointer :: emissivity(:)   ! Input/output surface emissivity
    type(rttov_profile), pointer    :: profiles(:)     ! Input profiles
    type(rttov_profile), pointer    :: profiles_k(:)   ! Input profiles
    type(rttov_chanprof), pointer   :: chanprof(:)     ! Input profiles
    type(rttov_transmission)        :: transmission    ! Output transmittances
    type(rttov_radiance)            :: radiance        ! Output radiances

    type(rttov_emissivity), pointer :: emissivity_k(:) !Input/output surface emissivity
    type(rttov_transmission)        :: transmission_k  ! Output transmittances
    type(rttov_radiance)            :: radiance_k      ! Output radiances

  contains 

    procedure :: alloc   => alloc_rttov
    procedure :: init_emissivity => ufo_rttov_init_emissivity 

  end type ufo_rttov_io

  !Type for general config
  type rttov_conf
    integer                               :: nsensors
    integer                               :: ngas

    character(len=MAXVARLEN), allocatable :: Absorbers(:)
    integer, allocatable                  :: Absorber_Id(:)

    character(len=255), allocatable       :: SENSOR_ID(:)
    character(len=255)                    :: COEFFICIENT_PATH

    type(rttov_coefs), allocatable        :: rttov_coef_array(:)
    character(len=8)                      :: RTTOV_default_opts = 'RTTOV'
    type(rttov_options)                   :: rttov_opts
    logical                               :: rttov_is_setup = .false.

    logical                               :: SatRad_compatibility = .true.
    logical                               :: qtotal = .true.
    logical                               :: UseQtsplitRain, SplitQtotal = .false. ! true for MW, false otherwise
    logical                               :: RTTOV_profile_checkinput = .false.

    integer                               :: inspect
    integer                               :: nchan_max_sim

  contains

    procedure :: set_options => set_options_rttov
    procedure :: setup =>    setup_rttov
    procedure :: set_defaults => set_defaults_rttov

  end type rttov_conf

contains

  ! ------------------------------------------------------------------------------

  subroutine rttov_conf_setup(conf, f_confOpts, f_confOper)
    implicit none

    type(rttov_conf), intent(inout)       :: conf
    type(fckit_configuration), intent(in) :: f_confOpts ! RTcontrol
    type(fckit_configuration), intent(in) :: f_confOper ! what is this

    character(*), parameter               :: routine_name = 'rttov_conf_setup'
    integer                               :: ivar, jspec
    character(len=:), allocatable         :: str
    character(len=:), allocatable         :: str_array(:)
    logical                               :: varin_satrad = .false.

    !Number of sensors, each call to RTTOV will be for a single sensor
    !type (zenith/scan angle will be different)
    conf % nSensors = 1

    if (f_confOper%has("Debug")) then
      call f_confOper % get_or_die("Debug",debug)
    else
      debug = .false. ! default
    endif

    if (f_confOper%has("GeoVal_type")) then
      call f_confOper%get_or_die("GeoVal_type",str)
      if (trim(str) == 'MetO' .or. trim(str) == 'SatRad') then
        varin_default => varin_default_satrad
        varin_satrad = .true.
      elseif(trim(str) == 'CRTM') then
        varin_default => varin_default_crtm
        varin_satrad = .false.
      else
        write(message,*) trim(ROUTINE_NAME),' error: ',trim(str),' is not a supported GeoVal type'
        call abor1_ftn(message)
      endif
    else
      varin_default => varin_default_satrad
    endif

    ! Absorbers
    !----------
    conf%ngas = 0
    if (f_confOper%has("Absorbers")) &
      conf%ngas = conf%ngas + f_confOper%get_size("Absorbers")

    allocate(conf%Absorbers( conf%ngas ), &
             conf%Absorber_Id( conf%ngas ))

    if (conf%ngas > 0) then
      call f_confOper%get_or_die("Absorbers",str_array)
      conf%Absorbers(1:conf%ngas) = str_array

    end if

    ! check for duplications
    do jspec = 2, conf%ngas
      if ( any(conf%Absorbers(jspec-1) == conf%Absorbers(jspec:conf%ngas)) ) then
        write(message,*) trim(ROUTINE_NAME),' error: ',trim(conf%Absorbers(jspec)),' is duplicated in Absorbers'
        call abor1_ftn(message)
      end if
    end do

    ! convert from CRTM names to UFO CF names and define Id and Units
    do jspec = 1, conf%ngas
      ivar = ufo_vars_getindex(RTTOV_Absorbers, conf%Absorbers(jspec))

      if (ivar < 1 .or. ivar > size(UFO_Absorbers)) then
        write(message,*) trim(ROUTINE_NAME),' error: ',trim(conf%Absorbers(jspec)),' not supported by UFO_Absorbers'
        call abor1_ftn(message)
      end if

      conf%Absorbers(jspec) = UFO_Absorbers(ivar)

      ! DARFIX replace humidity_mixing_ratio with specific_humidity
      ! this is starting to get a bit messy and dangerous
      if(conf%Absorbers(jspec) == var_mixr .and. varin_satrad) then
        conf%Absorbers(jspec) = var_q
      endif
      conf%Absorber_Id(jspec) = RTTOV_Absorber_Id(ivar)
    end do


    ! Allocate SENSOR_ID
    allocate(conf % SENSOR_ID(conf % nSensors))

    ! Get sensor ID from config
    call f_confOpts % get_or_die("Sensor_ID",str)
    conf % SENSOR_ID(conf%nSensors) = str

    ! Path to coefficient files
    call f_confOpts % get_or_die("CoefficientPath",str)
    conf % COEFFICIENT_PATH = str

    if(f_confOpts % has("RTTOV_default_opts")) then
      call f_confOpts % get_or_die("RTTOV_default_opts",str) !test, OPS, RTTOV
      conf % RTTOV_default_opts = str
    endif

    if(f_confOpts % has("SatRad_compatibility")) then
      call f_confOpts % get_or_die("SatRad_compatibility",conf % SatRad_compatibility)
    endif

    if(f_confOpts % has("qtotal")) then
      call f_confOpts % get_or_die("qtotal",conf % qtotal)
    endif

    if(f_confOpts % has("max_channels_per_batch")) then
      call f_confOpts % get_or_die("max_channels_per_batch",conf % nchan_max_sim)
    else
      conf % nchan_max_sim = 10000
    endif

    if( .not. conf % rttov_is_setup) then
      call conf % setup(f_confOpts, asw=1)
    end if

    !DARFIX THIS ONLY WORKS FOR ONE INSTRUMENT
    if (conf % rttov_coef_array(1) % coef % id_sensor == sensor_id_mw) then
      if(conf % rttov_opts % rt_mw % clw_data .and. &
         conf % SatRad_compatibility .and. conf % qtotal) then
        conf % UseQtsplitRain = .true.
        conf % splitQtotal = .true.
      endif

      conf % rttov_opts % rt_ir % ozone_data = .false.
      conf % rttov_opts % rt_ir % co2_data = .false.
      conf % rttov_opts % rt_ir % n2o_data = .false.
      conf % rttov_opts % rt_ir % ch4_data = .false.
      conf % rttov_opts % rt_ir % so2_data = .false.
    else
      conf % UseQtsplitRain = .false.
      conf % splitQtotal = .false.
    endif

    if(f_confOpts % has("RTTOV_profile_checkinput")) then
      call f_confOpts % get_or_die("RTTOV_profile_checkinput",conf % RTTOV_profile_checkinput) !test, OPS, RTTOV
    endif

    conf % inspect = 0
    if (f_confOpts%has("InspectProfileNumber")) then
      call f_confOpts % get_or_die("InspectProfileNumber",conf % inspect)
      conf % RTTOV_profile_checkinput = .true.
    endif

  end subroutine rttov_conf_setup

  ! -----------------------------------------------------------------------------

  subroutine rttov_conf_delete(conf)

    implicit none
    type(rttov_conf), intent(inout) :: conf

    deallocate(conf%SENSOR_ID)
    deallocate(conf%Absorbers)
    deallocate(conf%Absorber_Id)

  end subroutine rttov_conf_delete

  ! -----------------------------------------------------------------------------

  subroutine load_atm_data_rttov(geovals,obss,profiles,prof_start,conf,layer_quantities,ob_info)

    use ufo_constants_mod, only : half, deg2rad, min_q, m_to_km, g_to_kg, pa_to_hpa
    use ufo_rttovonedvarcheck_ob_mod

    implicit none

    type(ufo_geovals),            intent(in)    :: geovals
    type(c_ptr), value,           intent(in)    :: obss
    type(rttov_profile),          intent(inout) :: profiles(:)
    integer,                      intent(in)    :: prof_start
    type(rttov_conf),             intent(in)    :: conf
    logical,                      intent(inout) :: layer_quantities
    type(ufo_rttovonedvarcheck_ob), optional, intent(inout) :: ob_info

    ! Local variables
    integer                      :: jspec
    integer                      :: nlevels
    integer                      :: nprofiles

    type(ufo_geoval), pointer    :: geoval
    character(MAXVARLEN)         :: varname

    real(kind_real)              :: ifrac, sfrac, lfrac, wfrac
    real(kind_real)              :: itmp, stmp, ltmp

    real(kind_real), allocatable :: TmpVar(:), windsp(:), p(:)
    logical                      :: variable_present

    integer                      :: top_level, bottom_level, stride
    real(kind_real)              :: NewT
    integer                      :: level_1000hPa, level_950hpa

    real(kind_real), allocatable :: q_temp(:), clw_temp(:), ciw_temp(:), Qtotal(:)

    if(present(ob_info)) then
      nlocs_total = 1
    else
      nlocs_total = obsspace_get_nlocs(obss)
    end if

    nprofiles = min(size(profiles), geovals%nlocs - prof_start + 1)

    nlevels = size(profiles(1)%p)

    ! Assume that the pressure profile coming from the geovals is increasing in pressure (ToA->surface)...
    top_level = 1
    bottom_level = nlevels
    stride=1

    if (ufo_vars_getindex(geovals%variables, var_prs) > 0) then
      call ufo_geovals_get_var(geovals, var_prs, geoval)

      allocate(p(nlevels))
      do iprof = 1, nprofiles
        p = geoval%vals(geoval%nval-nlevels+1:geoval%nval,prof_start + iprof - 1) * Pa_to_hPa ! for RTTOV
        if (iprof == 1) then
          if (p(1) > p(2)) then ! ...but be ready to switch. Assume if one profile is 'upside-down' then they all are
            top_level = nlevels
            bottom_level = 1
            stride = -1
          endif
        endif
        profiles(iprof)%p(top_level:bottom_level:stride) = p(:)
      end do
      deallocate(p)

      ! DARFIX - it might be that all of the temperatures are on layers so it may be that we need to do interpolation
      ! but it might be better to force this using an option if we know where the geovals are coming from
    else
      call ufo_geovals_get_var(geovals, var_prsi, geoval)
      allocate(p(nlevels-1))
      do iprof = 1, nprofiles
        p = geoval%vals(geoval%nval-(nlevels-1)+1:geoval%nval,prof_start + iprof - 1) * Pa_to_hPa
        if (iprof == 1) then
          if (p(1) > p(2)) then !upside-down
            top_level = nlevels - 1
            bottom_level = 1
            stride = -1
          else
            top_level = 1
            bottom_level = nlevels - 1
            stride = 1       
          endif
        endif
        profiles(iprof)%p(top_level+stride:bottom_level:stride) = half * (p(top_level:bottom_level-stride:stride) + p(top_level+stride:bottom_level:stride))
        profiles(iprof)%p(1) = max( profiles(iprof)%p(2) - half * (profiles(iprof)%p(3) - profiles(iprof)%p(2)),half * profiles(iprof)%p(2))
        profiles(iprof)%p(nlevels) = profiles(iprof)%p(nlevels-1) - half * (profiles(iprof)%p(nlevels-2) - profiles(iprof)%p(nlevels-1))
      end do
      deallocate(p)
    endif

! Get temperature
    varname = var_ts
    call ufo_geovals_get_var(geovals, varname, geoval)

! Check if temperatures are provided on levels as required for RTTOV, otherwise assume that temperatures are layer quantities
! (and all future atmospheric variables) and do some interpolation to prepare for RTTOV. 
! TODO: Put a warning in here that this is happening
    if( size(geoval%vals(:,1)) == nlevels) then
      do iprof = 1, nprofiles
        profiles(iprof)%t(top_level:bottom_level:stride) = geoval%vals(:,prof_start + iprof - 1) ! K
      end do
    else
      layer_quantities = .true.
      bottom_level = nlevels-1

      do iprof = 1, nprofiles
        profiles(iprof)%t(top_level+stride:bottom_level:stride) = half * (profiles(iprof)%t(top_level:bottom_level-stride:stride) + &
          profiles(iprof)%t(top_level+stride:bottom_level:stride))
        profiles(iprof)%t(1) = profiles(iprof)%t(2)
        profiles(iprof)%t(nlevels) = profiles(iprof)%t(nlevels-1) - half * (profiles(iprof)%t(nlevels-2) - profiles(iprof)%t(nlevels-1))
      end do
    endif

! Get absorbers. Assume mass mixing ratio (moist air). The distinction has been made that q will be in kg/kg as per OPS convention
! And mixr will be in g/kg.
! For now all other gases will be in kg/kg and this needs to be handled carefully.

    profiles(1:nprofiles)%gas_units = 1

    do jspec = 1, conf%ngas
      call ufo_geovals_get_var(geovals,conf%Absorbers(jspec) , geoval)

      select case (conf%Absorbers(jspec))
      case (var_mixr)
        do iprof = 1, nProfiles
          profiles(iprof)%q(top_level:bottom_level:stride) = geoval%vals(:,prof_start + iprof - 1) * g_to_kg
        end do
      case (var_q)
        do iprof = 1, nProfiles
          profiles(iprof)%q(top_level:bottom_level:stride) = geoval%vals(:,prof_start + iprof - 1) ! kg/kg
        end do
      case (var_oz)
        if (associated(profiles(1)%o3)) then
          do iprof = 1, nProfiles
            profiles(iprof)%o3(top_level:bottom_level:stride) = geoval%vals(:,prof_start + iprof - 1) ! kg/kg
          end do
        endif
      case (var_co2)
        if (associated(profiles(1)%co2)) then
          do iprof = 1, nProfiles
            profiles(iprof)%co2(top_level:bottom_level:stride) = geoval%vals(:,prof_start + iprof - 1) ! kg/kg
          end do
        endif
      case (var_clw)
        if (associated(profiles(1)%clw)) then
          do iprof = 1, nProfiles
            profiles(iprof)%clw(top_level:bottom_level:stride) = geoval%vals(:,prof_start + iprof - 1) ! kg/kg
          end do
        endif
      case default

      end select

    end do

! Get near-surface variables (s2m)

    varname = var_sfc_p2m
    if (ufo_vars_getindex(geovals%variables, varname) > 0) then
      call ufo_geovals_get_var(geovals, varname, geoval)

      profiles(1:nprofiles)%s2m%p = geoval%vals(1,prof_start:prof_start + nprofiles - 1) * Pa_to_hPa
    else
      write(message,'(A)') 'No near-surface pressure. Using bottom pressure level'
      call fckit_log%info(message)

      do iprof = 1, nprofiles
        profiles(iprof)%s2m%p = profiles(iprof)%p(nlevels)
      enddo
    endif

    varname = var_sfc_t2m ! 2m temperature
    if (ufo_vars_getindex(geovals%variables, varname) > 0) then
      call ufo_geovals_get_var(geovals, varname, geoval) 
      profiles(1:nprofiles)%s2m%t = geoval%vals(1,prof_start:prof_start + nprofiles - 1)
    else
      write(message,'(A)') 'No near-surface temperature. Using bottom temperature level'
      call fckit_log%info(message)

      do iprof = 1, nprofiles
        profiles(iprof)%s2m%t = profiles(iprof)%t(nlevels)
      enddo
    endif

    varname = var_sfc_q2m ! 2m specific humidity (kg/kg)
    if (ufo_vars_getindex(geovals%variables, varname) > 0) then
      call ufo_geovals_get_var(geovals, varname, geoval) ! lfric

      profiles(1:nprofiles)%s2m%q = geoval%vals(1,prof_start:prof_start + nprofiles - 1)
    else
      write(message,'(A)') 'No near-surface specific humidity. Using bottom q level'
      call fckit_log%info(message)

      do iprof = 1, nprofiles
        profiles(iprof)%s2m%q = profiles(iprof)%q(nlevels)
      enddo
    endif

    varname = var_u ! Eastward-wind in m/s 
    if (ufo_vars_getindex(geovals%variables, varname) > 0) then
      call ufo_geovals_get_var(geovals, varname, geoval)

      profiles(1:nprofiles)%s2m%u = geoval%vals(1,prof_start:prof_start + nprofiles - 1)
      !assume if eastward then northward too

      varname = var_v ! Northward-wind in m/s 
      call ufo_geovals_get_var(geovals, varname, geoval)

      profiles(1:nprofiles)%s2m%v = geoval%vals(1,prof_start:prof_start + nprofiles - 1)
    else !! use windspeed and direction instead
      allocate(windsp(nprofiles))
      call ufo_geovals_get_var(geovals, var_sfc_wspeed, geoval)

      windsp(1:nprofiles) = geoval%vals(1,prof_start:prof_start + nprofiles - 1)

      call ufo_geovals_get_var(geovals, var_sfc_wdir, geoval)

      do iprof = 1, nprofiles
        profiles(iprof)%s2m%u             = windsp(iprof) * cos(geoval%vals(1,prof_start + iprof - 1) * deg2rad)
        profiles(iprof)%s2m%v             = windsp(iprof) * sin(geoval%vals(1,prof_start + iprof - 1) * deg2rad)
      enddo
      deallocate(windsp)
    endif

! Get surface type from geoval if it exists else diagnose surface type from water fraction.
! If RTTOV surface type exists then it is expected that the related variables will also be available
! e.g. water type, Tskin
! DARFIX : This will be moved out to a separate subroutine in the next release to support all instruments
!          particularly where we have additional surface metadata for certain instruments

    varname = var_surf_type_rttov ! RTTOV surface type: 0 (land), 1 (water), 2 (sea-ice)
    if (ufo_vars_getindex(geovals%variables, varname) > 0) then
      call ufo_geovals_get_var(geovals, varname, geoval)
      profiles(1:nprofiles)%skin%surftype = int(geoval%vals(1,prof_start:prof_start + nprofiles - 1),kind=jpim)

      varname = var_water_type_rttov ! RTTOV water type: 0 (fresh), 1 (sea)
      call ufo_geovals_get_var(geovals, varname, geoval)
      profiles(1:nprofiles)%skin%watertype = int(geoval%vals(1,prof_start:prof_start + nprofiles - 1),kind=jpim)

      varname = var_sfc_tskin !Skin (surface) temperature (K)
      call ufo_geovals_get_var(geovals, varname, geoval)
      profiles(1:nprofiles)%skin%t = geoval%vals(1,prof_start:prof_start + nprofiles - 1)

      
    else

      ! Try to diagnose RTTOV surface from water fraction
      !DARFIX: this is not consistent with CRTM in any way. Need to find out how they select surface type.

      varname = var_sfc_wfrac
      if (ufo_vars_getindex(geovals%variables, varname) > 0) then
        call ufo_geovals_get_var(geovals, varname, geoval)

        do iprof = 1, nprofiles
          !Land point or sea point
          wfrac = geoval%vals(1,iprof)
          if (wfrac > half) then
            profiles(iprof)%skin%surftype   = surftype_sea
            call ufo_geovals_get_var(geovals, var_sfc_wtmp, geoval)
            profiles(iprof)%skin%t   = geoval%vals(1,prof_start + iprof - 1)
          else
            !maybe it's predominantly land or ice
            !   !determine land, snow and ice fractions and temperatures to determine average temperature
            profiles(iprof)%skin%surftype   = surftype_land ! land

            call ufo_geovals_get_var(geovals, var_sfc_lfrac, geoval) 
            lfrac   = geoval%vals(1,prof_start + iprof - 1)

            call ufo_geovals_get_var(geovals, var_sfc_sfrac, geoval) 
            sfrac   = geoval%vals(1,prof_start + iprof - 1)

            call ufo_geovals_get_var(geovals, var_sfc_ifrac, geoval) 
            ifrac   = geoval%vals(1,prof_start + iprof - 1)
            
            call ufo_geovals_get_var(geovals, var_sfc_ltmp, geoval)
            ltmp   = geoval%vals(1,prof_start + iprof - 1)
            
            call ufo_geovals_get_var(geovals, var_sfc_stmp, geoval)
            stmp   = geoval%vals(1,prof_start + iprof - 1)
            
            call ufo_geovals_get_var(geovals, var_sfc_itmp, geoval)
            itmp   = geoval%vals(1,prof_start + iprof - 1)
            
            !Skin temperature is a combination of (i)ce temp, (l)and temp and (s)now temp
            profiles(iprof)%skin%t   = (lfrac * ltmp + sfrac * stmp + ifrac * itmp) / (lfrac + sfrac + ifrac)
          endif
        end do
      endif
    endif

    !DAR: Salinity fixed for now too
    profiles(1:nprofiles)%skin%salinity = 35.0_kind_real

    !DAR: Default fastem parameters. We are not using FASTEM over land so these are unused
    do iprof = 1,nprofiles
!      profiles(iprof)%skin%fastem            = [3.0, 5.0, 15.0, 0.1, 0.3]
      profiles(iprof)%skin%fastem            = [0,0,0,0,0]
    end do

! ---------------------------
! SatRad profile manipulation
! ---------------------------

    do iprof = 1, nProfiles
      if(conf % SatRad_compatibility) then

        !----
        ! Reset low level temperatures over seaice and cold, low land as per Ops_SatRad_SetUpRTprofBg.F90
        ! N.B. I think this should be flagged so it's clear that the background has been modified
        !----
        
        if(profiles(iprof)%skin%surftype /= surftype_sea) then
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
          endif
        endif

        !min_q fix
        where(profiles(iprof)%q < min_q) profiles(iprof)%q = min_q
        if(profiles(iprof)%s2m%q < min_q) profiles(iprof)%s2m%q = min_q
        
      endif

    enddo

! ----------------------------
! Interpolate layer quantities
! ----------------------------
! Continuation of the profile interpolation if absorber profiles are layer quantities (RTTOV expects levels).
! This is done once any profile manipulation (for q primarily) has been done.

    if(layer_quantities) then
      do iprof=1,nprofiles
        profiles(iprof)%q(top_level+stride:bottom_level:stride) = half * (profiles(iprof)%q(top_level:bottom_level-stride:stride) + &
          profiles(iprof)%q(top_level+stride:bottom_level:stride))
        profiles(iprof)%q(1) = profiles(iprof)%q(2)
        profiles(iprof)%q(nlevels) = profiles(iprof)%q(nlevels-1) - half * (profiles(iprof)%q(nlevels-2) - profiles(iprof)%q(nlevels-1))

        if (associated(profiles(1)%o3)) then
          profiles(iprof)%o3(top_level+stride:bottom_level:stride) = half * &
            (profiles(iprof)%o3(top_level:bottom_level-stride:stride) + profiles(iprof)%o3(top_level+stride:bottom_level:stride))
          profiles(iprof)%o3(1) = profiles(iprof)%o3(2)
          profiles(iprof)%o3(nlevels) = profiles(iprof)%o3(nlevels-1) - &
            half * (profiles(iprof)%o3(nlevels-2) - profiles(iprof)%o3(nlevels-1))
        endif

        if (associated(profiles(1)%clw)) then
          profiles(iprof)%clw(top_level+stride:bottom_level:stride) = half * &
            (profiles(iprof)%clw(top_level:bottom_level-stride:stride) + profiles(iprof)%clw(top_level+stride:bottom_level:stride))
          profiles(iprof)%clw(1) = profiles(iprof)%clw(2)
          profiles(iprof)%clw(nlevels) = profiles(iprof)%clw(nlevels-1) - &
            half * (profiles(iprof)%clw(nlevels-2) - profiles(iprof)%clw(nlevels-1))
        endif
      enddo
    endif


    !non mwscatt only at the moment
    if(conf % SplitQtotal) then

      allocate(Qtotal(nlevels), q_temp(nlevels), clw_temp(nlevels), ciw_temp(nlevels))

      do iprof = 1, nprofiles
        ! compute bg qtotal using q and clw only
        ! currently ice is ignored
        Qtotal(:) = profiles(iprof) % q(:) + profiles(iprof) % clw(:)
        
        ! Make sure there is a minimum humidity
        Qtotal(:) = max(Qtotal(:), Min_q)

        ! Make sure there is a minimum humidity
        Qtotal(:) = max(Qtotal(:), Min_q)

        ! generate first guess cloud and q based on the qtotal physics
        call Ops_SatRad_Qsplit (1,                     & ! in
          profiles(iprof) % p(:) / Pa_to_hPa,          & ! in convert hPa to Pa
          profiles(iprof) % t(:),                      & ! in
          Qtotal(:),                                   & ! in
          q_temp(:),                                   & ! out
          clw_temp(:),                                 & ! out
          ciw_temp(:),                                 & ! out
          UseQtSplitRain = conf % UseQtSplitRain)

        ! store in the profile
        profiles(iprof) % clw(:) = clw_temp(:)
        profiles(iprof) % q(:)   = q_temp(:)

        ! !store non active variable
        ! IF (ALLOCATED(CloudIce)) THEN

        !   CloudIce(toplevel_q_1dvar:) = ciw_temp(:)

        ! ENDIF
      enddo

      deallocate(Qtotal, q_temp, clw_temp, ciw_temp)

    end if

    if(present(ob_info)) then

       write(*,*) "load_atm_data_rttov: getting from ob info"
       profiles(1) % elevation = ob_info % elevation / 1000.0 ! m -> km
       profiles(1) % latitude = ob_info % latitude
       profiles(1) % longitude = ob_info % longitude
       if (ob_info % retrievecloud) then
         profiles(1) % ctp = ob_info % cloudtopp
         profiles(1) % cfraction = ob_info % cloudfrac
       end if

    else

      allocate(TmpVar(nlocs_total))
      variable_present = obsspace_has(obss, "MetaData", "elevation")
      if (variable_present) then
        call obsspace_get_db(obss, "MetaData", "elevation", TmpVar)
        profiles(1:nprofiles)%elevation = TmpVar(prof_start:prof_start + nprofiles - 1) * m_to_km !for RTTOV
      else
        variable_present = obsspace_has(obss, "MetaData", "surface_height")
        if (variable_present) then
          call obsspace_get_db(obss, "MetaData", "surface_height", TmpVar)
          profiles(1:nprofiles)%elevation = TmpVar(prof_start:prof_start + nprofiles - 1) * m_to_km !for RTTOV
        else ! from geoval
          if (ufo_vars_getindex(geovals%variables, 'surface_altitude') > 0) then
            call ufo_geovals_get_var(geovals, 'surface_altitude', geoval)
            profiles(1:nprofiles)%elevation = geoval%vals(1,prof_start:prof_start + nprofiles - 1) * m_to_km
          else
            write(message,'(A)') &
              'MetaData elevation not in database: check implicit filtering'
            call fckit_log%info(message)
          endif
        endif
      end if

      variable_present = obsspace_has(obss, "MetaData", "latitude")
      if (variable_present) then
        call obsspace_get_db(obss, "MetaData", "latitude", TmpVar )
        profiles(1:nprofiles)%latitude = TmpVar(prof_start:prof_start + nprofiles - 1)
      else
        write(message,'(A)') &
          'MetaData latitude not in database: check implicit filtering'
        call fckit_log%info(message)
      end if

      variable_present = obsspace_has(obss, "MetaData", "longitude")
      if (variable_present) then
        call obsspace_get_db(obss, "MetaData", "longitude", TmpVar )
        profiles(1:nprofiles)%longitude = TmpVar(prof_start:prof_start + nprofiles - 1)
      else
        write(message,'(A)') &
          'MetaData longitude not in database: check implicit filtering'
        call fckit_log%info(message)
      end if

    end if

  end subroutine load_atm_data_rttov

  ! Internal subprogam to load some test geometry data
  subroutine load_geom_data_rttov(obss,profiles,prof_start1,ob_info)

    use obsspace_mod, only :  obsspace_get_nlocs, obsspace_get_db
    use ufo_rttovonedvarcheck_ob_mod

    implicit none

    type(c_ptr), value,           intent(in)    :: obss
    type(rttov_profile),          intent(inout) :: profiles(:)
    integer, optional,            intent(in)    :: prof_start1
    type(ufo_rttovonedvarcheck_ob), optional, intent(inout) :: ob_info

    real(kind_real), allocatable                :: TmpVar(:)

    integer                                     :: prof_start
    integer                                     :: nprofiles
    integer                                     :: nlevels

    if(present(prof_start1)) then
      prof_start = prof_start1
    else
      prof_start = 1
    end if

    if(present(ob_info)) then

      nlocs_total = 1
      nprofiles = 1
      nlevels = SIZE(profiles(1) % p)

      profiles(1) % zenangle    = ob_info % sensor_zenith_angle
      profiles(1) % azangle     = ob_info % sensor_azimuth_angle
      profiles(1) % sunzenangle = ob_info % solar_zenith_angle
      profiles(1) % sunazangle  = ob_info % solar_azimuth_angle

    else

      nprofiles = min(size(profiles), nlocs_total - prof_start + 1) !nlocs_total set in read_atm

      allocate(TmpVar(nlocs_total))

      nlevels = size(profiles(1)%p)

      !RTTOV convention 0-max (coef dependent). Nadir = 0 deg
      call obsspace_get_db(obss, "MetaData", "sensor_zenith_angle", TmpVar)
      if (any(TmpVar(prof_start:prof_start + nprofiles - 1) < 0.)) then
        profiles(1:nprofiles)%zenangle = abs(TmpVar(prof_start:prof_start + nprofiles - 1))
        ! DAR need to make a modification to azangle here then?
      else
        profiles(1:nprofiles)%zenangle = abs(TmpVar(prof_start:prof_start + nprofiles - 1))
      endif

      !RTTOV convention is 0-360 deg. E=+90
      call obsspace_get_db(obss, "MetaData", "sensor_azimuth_angle", TmpVar)
      profiles(1:nprofiles)%azangle = TmpVar(prof_start:prof_start + nprofiles - 1)

      call obsspace_get_db(obss, "MetaData", "solar_zenith_angle", TmpVar)
      if (any(TmpVar(prof_start:prof_start + nprofiles - 1) < 0.)) then
        profiles(1:nprofiles)%sunzenangle = abs(TmpVar(prof_start:prof_start + nprofiles - 1))
        ! DAR need to make a modification to sunazangle here then?
      else
        profiles(1:nprofiles)%sunzenangle = abs(TmpVar(prof_start:prof_start + nprofiles - 1))
      endif

      call obsspace_get_db(obss, "MetaData", "solar_azimuth_angle", TmpVar)
      profiles(1:nprofiles)%sunazangle = TmpVar(prof_start:prof_start + nprofiles - 1)

      deallocate(TmpVar)

    end if

  end subroutine load_geom_data_rttov

  ! ------------------------------------------------------------------------------

  subroutine set_options_rttov(self, f_confOpts)
    implicit none

    class(rttov_conf),         intent(inout) :: self
    type(fckit_configuration), intent(in)    :: f_confOpts ! RTcontrol

    include 'rttov_print_opts.interface'

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
      call f_confOpts % get_or_die("RTTOV_", self % rttov_opts % rt_ir % user_cld_opt_param)
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

    ! OPEN(666,file='RTTOV_options.log')
    ! CALL rttov_print_opts(self % rttov_opts,lu=666)
    ! CLOSE(666)
  end subroutine set_options_rttov

  ! ------------------------------------------------------------------------------

  subroutine setup_rttov(self, f_confOpts, asw)
    class(rttov_conf), intent(inout) :: self
    type(fckit_configuration), intent(in) :: f_confOpts ! RTcontrol
    integer, intent(in) :: asw !allocate switch

    character(len=255) :: coef_filename
    character(len=4) :: coef_ext
    integer :: i_inst

    include 'rttov_read_coefs.interface'

    rttov_errorstatus = 0
    coef_ext = '.dat'

    if(asw == 1) then

      ! --------------------------------------------------------------------------
      ! 1. Setup rttov options
      ! --------------------------------------------------------------------------
      call self % set_options(f_confOpts)

      ! --------------------------------------------------------------------------
      ! 2. Read coefficients
      ! --------------------------------------------------------------------------
      allocate(self % rttov_coef_array(self % nSensors))

      do i_inst = 1, self%nSensors
        coef_filename = &
          trim(self % COEFFICIENT_PATH) // 'rtcoef_' // trim(self%SENSOR_ID(i_inst)) // trim(coef_ext)

        call rttov_read_coefs(rttov_errorstatus, &       !out
                              self % rttov_coef_array(i_inst), & !inout
                              self % rttov_opts, &           !in
                              file_coef = coef_filename)     !in

        if (rttov_errorstatus /= errorstatus_success) then
          write(*,*) 'fatal error reading coefficients'
        else
          write(*,*) 'successfully read' // coef_filename
        end if

      end do

      self % rttov_is_setup =.true.
    else !asw == 0
      deallocate(self % rttov_coef_array)
      self%rttov_is_setup =.false.
    end if

  end subroutine setup_rttov

  ! ------------------------------------------------------------------------------

  subroutine get_var_name(n,varname)

    integer, intent(in) :: n
    character(len=*), intent(out) :: varname

    character(len=6) :: chan

    write(chan, '(I0)') n
    varname = 'brightness_temperature_' // trim(chan)

  end subroutine get_var_name

  !ufo_rttov_alloc is a wrapper for RTTOV12/13 allocation
  subroutine alloc_rttov(self, errorstatus, conf, nprofiles, nchannels, nlevels, init, asw)
    implicit none

    class(ufo_rttov_io), intent(inout) :: self
    type(rttov_conf),    intent(in)    :: conf

    integer,             intent(out)   :: errorstatus ! Return error status of RTTOV subroutine calls
    integer,             intent(in)    :: asw ! allocate switch
    logical, optional,   intent(in)    :: init ! initialise, default yes
    integer ,            intent(in)    :: nchannels
    integer ,            intent(in)    :: nprofiles
    integer ,            intent(in)    :: nlevels

    logical                            :: init1
    integer                            :: asw_direct, asw_k

    integer, save                      :: called_from

    include 'rttov_alloc_direct.interface'
    include 'rttov_alloc_k.interface'

    asw_direct = -1
    asw_k = -1

    if(asw == 1) then
      asw_direct = 1
      called_from = 1
    endif
    if(asw == 2) then
      asw_k = 1
      called_from = 2
    endif
    
    if (present(init)) then
      init1 = init
    else 
      init1 = .true.
    endif

    if(asw_direct /= -1 .and. called_from == 1) then
      ! Allocate structures for rttov_direct
      call rttov_alloc_direct(          &
        errorstatus,                    &
        asw_direct,                     &  ! 1 => allocate
        nprofiles,                      &
        nchannels,                      &
        nlevels,                        &
        self % chanprof,                &
        conf % rttov_opts,              &
        self % profiles,                &
        conf % rttov_coef_array(1),     &
        self % transmission,            &
        self % radiance,                &
        calcemis = self % calcemis,     &
        emissivity = self % emissivity, &
        init=init1)

      if (errorstatus /= errorstatus_success) then
        write(message,'(A, I6)') 'after rttov_alloc_direct error = ', errorstatus
        call abor1_ftn(message)
      end if

    endif

    if (asw_k /= -1 .and. called_from == 2) then
      ! Allocate structures for rttov_k
      call rttov_alloc_k(                   &
        errorstatus,                        &
        asw_k,                              &  ! 1 => allocate
        nprofiles,                          &
        nchannels,                          &
        NLEVELS,                            &
        self % chanprof,                    &
        conf % rttov_opts,                  &
        self % profiles,                    &
        self % profiles_k,                  &
        conf % rttov_coef_array(1),         &
        self % transmission,                &
        self % transmission_k,              &
        self % radiance,                    &
        self % radiance_k,                  &
        calcemis = self % calcemis,         &
        emissivity = self % emissivity,     &
        emissivity_k = self % emissivity_k, &
        init=init1)

      if (errorstatus /= errorstatus_success) then
        write(message,'(A, I6)') 'after rttov_alloc_k error = ', errorstatus
        call abor1_ftn(message)
      end if
    endif

    if ( asw /= 0 ) then
      self % calcemis = .false.
      self % emissivity % emis_in = -1.0_kind_real
      self % emissivity % emis_out = -1.0_kind_real
      self % profiles(:) % skin % surftype = -1_jpim
    endif

  end subroutine alloc_rttov

  subroutine ufo_rttov_skip_profiles(nprofiles,nchannels,channels,obss,Skip_Profiles)
    ! Profiles are skipped when the ObsValue of all channels is missing.
    ! TODO: Use complete QC information
    ! It would be more comprehensive to use EffectiveQC or EffectiveError. That
    ! would require those ObsSpace values to be initialized before calls to
    ! this subroutine within ufo_radiancerttov_simobs+ufo_radiancerttov_tlad_settraj.
    implicit none

    integer,              intent(in)    :: nProfiles, nChannels
    type(c_ptr), value,   intent(in)    :: obss
    integer(c_int),       intent(in)    :: channels(:)
    logical,              intent(inout) :: Skip_Profiles(:)

    integer                  :: jprofile, jchannel
    character(len=MAXVARLEN) :: varname
    real(kind_real)          :: ObsVal(nprofiles,nchannels)
    !real(kind_real)         :: EffObsErr(n_Profiles,n_Channels)
    !integer                 :: EffQC(n_Profiles,n_Channels)

    ! Set missing value
    missing = missing_value(missing)

    ObsVal = missing
    ! EffObsErr = missing
    ! EffQC = 0

    do jchannel = 1, nchannels
      call get_var_name(channels(jchannel),varname)
      call obsspace_get_db(obss, "ObsValue", varname, ObsVal(:,jchannel))
      !   call obsspace_get_db(obss, "EffectiveError", varname, EffObsErr(:,jchannel))
      !   call obsspace_get_db(obss, "EffectiveQC{iter}", varname, EffQC(:,jchannel))
    enddo

    !Loop over all n_Profiles, i.e. number of locations
    do jprofile = 1, nprofiles
      Skip_Profiles(jprofile) = all(ObsVal(jprofile,:) == missing) .or. all(ObsVal(jprofile,:) > 400)
      !            .OR. all(EffObsErr(jprofile,:) == missing) &
      !            .OR. all(EffQC(jprofile,:) /= 0)
    end do

  end subroutine ufo_rttov_skip_profiles

  subroutine ufo_rttov_init_emissivity(self, conf)
    class(ufo_rttov_io), intent(inout) :: self
    type(rttov_conf),    intent(in)    :: conf

    integer :: prof, ichan

    self % emissivity(:) % emis_in = -1.0_kind_real ! RMDI?
    self % emissivity(:) % emis_out = -1.0_kind_real ! RMDI?
    self % calcemis(:) = .false.

!Emissivity and calcemis are only set for used channels. So if a profile is skipped then you must not set emis data for the channels that are skipped 

    if ( conf % rttov_coef_array(1) % coef % id_sensor == sensor_id_mw) then
      do ichan = 1, nchan_sim, nchan_inst ! all channels initialised equally
        prof = self % chanprof(ichan)%prof
        if (self % profiles(prof) % skin % surftype == surftype_sea) then

          self % emissivity(ichan:ichan + nchan_inst - 1) % emis_in = 0.0_kind_real
          self % calcemis(ichan:ichan + nchan_inst - 1) = .true.
        else
          !IF ATLAS
            
          !ELSE
          if (self % profiles(prof) % skin % surftype == surftype_land) then
            
            !IF FASTEM (not implemented at MetO) ELSE
            self % emissivity(ichan:ichan + nchan_inst - 1) % emis_in = 0.95_kind_real
          elseif (self % profiles(prof) % skin % surftype == surftype_seaice) then
            
            !IF FASTEM (not implemented at MetO) ELSE
            self % emissivity(ichan:ichan + nchan_inst - 1) % emis_in = 0.92_kind_real
          endif
          !ENDIF !ATLAS
        endif
      enddo
    elseif ( conf % rttov_coef_array(1) % coef % id_sensor == sensor_id_ir .or. &
      conf % rttov_coef_array(1) % coef % id_sensor == sensor_id_hi) then

      do ichan = 1, size (self % chanprof(:)), nchan_inst ! all channels initialised equally
        prof = self % chanprof(ichan)%prof
        if (self % profiles(prof) % skin % surftype == surftype_sea) then
          ! Calculate by SSIREM or IREMIS
          self % emissivity(ichan:ichan + nchan_inst - 1) % emis_in = 0.0_kind_real
          self % calcemis(ichan:ichan + nchan_inst - 1) = .true.
        else
          !IF ATLAS ! CAMEL
          !ELSE
          if (self % profiles(prof) % skin % surftype == surftype_land) then
            self % emissivity(ichan:ichan + nchan_inst - 1) % emis_in = 0.95_kind_real
          elseif (self % profiles(prof) % skin % surftype == surftype_seaice) then
            !IF FASTEM (not implemented at MetO) ELSE
            self % emissivity(ichan:ichan + nchan_inst - 1) % emis_in = 0.92_kind_real
          endif
          !ENDIF !ATLAS
        endif
      enddo
    endif
      
  end subroutine ufo_rttov_init_emissivity

  subroutine set_defaults_rttov(self, default_opts_set)

    class(rttov_conf), intent(inout) :: self
    character(8),      intent(in)    :: default_opts_set

    write(message,'(A, A)') 'Setting RTTOV default options to ', default_opts_set
    call fckit_log%info(message)

    select case (trim(default_opts_set))
    case ('RTTOV') ! RTTOV defaults - needn't set this...
      self % rttov_opts % rt_all % addrefrac               = .false. !< Switch to enable atmospheric refraction
      self % rttov_opts % rt_all % switchrad               = .false. !< Switch for input units in AD/K models
      self % rttov_opts % rt_all % use_q2m                 = .true.  !< Switch to enable use of 2m q variable
      self % rttov_opts % rt_all % do_lambertian           = .false. !< Switch for setting Lambertian reflection (IR and MW)
      self % rttov_opts % rt_all % lambertian_fixed_angle  = .true.  !< Switch for fixed/parameterised effective angle for Lambertian option
      self % rttov_opts % rt_all % plane_parallel          = .false. !< Switch to ignore atmospheric curvature           
      self % rttov_opts % rt_all % rad_down_lin_tau        = .true.  !< Linear-in-tau or layer-mean for downwelling radiances    
      self % rttov_opts % rt_all % dtau_test               = .true.  !< Switch to apply dtau test in transmit/integrate calculations

      self % rttov_opts % rt_ir % solar_sea_brdf_model     = 1     !< Solar sea BRDF model (1-2)
      self % rttov_opts % rt_ir % ir_sea_emis_model        = 2     !< IR sea emissivity model (1-2)
      self % rttov_opts % rt_ir % addsolar                 = .false.  !< Switch to enable solar simulations
      self % rttov_opts % rt_ir % rayleigh_single_scatt    = .true.  !< Switch to enable Rayleigh single-scattering for VIS/NIR channels
      self % rttov_opts % rt_ir % do_nlte_correction       = .false.  !< Switch to enable NLTE bias correction
      self % rttov_opts % rt_ir % addaerosl                = .false.  !< Switch to enable IR aerosol calculations
      self % rttov_opts % rt_ir % user_aer_opt_param       = .false.  !< Switch to supply aerosol optical properties explicitly per channel
      self % rttov_opts % rt_ir % addclouds                = .false.  !< Switch to enable IR cloudy calculations
      self % rttov_opts % rt_ir % user_cld_opt_param       = .false.  !< Switch to supply cloud optical properties explicitly per channel
      self % rttov_opts % rt_ir % grid_box_avg_cloud       = .false.  !< Switch to supply grid-box average cloud concentration or cloud
      !!  concentration in cloudy fraction of each layer
      self % rttov_opts % rt_ir % cldstr_threshold         = -1.0_kind_real !< Ignore cloud streams with weights lower than this
      self % rttov_opts % rt_ir % cldstr_simple            = .false.  !< Switch for simplified cloud stream option - USE WITH CAUTION
      self % rttov_opts % rt_ir % cldstr_low_cloud_top     = 750._kind_real !< Upper pressure limit for cldstr_simple option (hPa)
      self % rttov_opts % rt_ir % ir_scatt_model           = ir_scatt_chou  !< IR scattering model to use
      self % rttov_opts % rt_ir % vis_scatt_model          = vis_scatt_dom  !< VIS/NIR scattering model to use
      self % rttov_opts % rt_ir % dom_nstreams             = 8     !< Number of DOM streams, must be even and not less than 2
      self % rttov_opts % rt_ir % dom_accuracy             = 0._kind_real  !< Convergence criterion for termination of DOM azimuthal loop
      self % rttov_opts % rt_ir % dom_opdep_threshold      = 0._kind_real  !< DOM ignores levels below this optical depth:
      !!  10. is reasonable, not applied if <              = 0
      self % rttov_opts % rt_ir % ozone_data               = .false.       !< Switch to enable input of O3 profile
      self % rttov_opts % rt_ir % co2_data                 = .false.       !< Switch to enable input of CO2 profile
      self % rttov_opts % rt_ir % n2o_data                 = .false.       !< Switch to enable input of N2O profile
      self % rttov_opts % rt_ir % co_data                  = .false.       !< Switch to enable input of CO profile
      self % rttov_opts % rt_ir % ch4_data                 = .false.       !< Switch to enable input of CH4 profile
      self % rttov_opts % rt_ir % so2_data                 = .false.       !< Switch to enable input of SO2 profile                      

      self % rttov_opts % rt_ir % pc % addpc               = .false.  !< Switch to enable PC-RTTOV
      self % rttov_opts % rt_ir % pc % ipcbnd              = -1    !< PC spectral band
      self % rttov_opts % rt_ir % pc % ipcreg              = -1    !< PC predictor channel set
      self % rttov_opts % rt_ir % pc % npcscores           = -1    !< Number of PC scores to compute, if less than 1 npcscores is derived
      !!  from the size of the pccomp%pcscores array
      self % rttov_opts % rt_ir % pc % addradrec           = .false.  !< Switch for calculation of reconstructed radiances

      !> MW-only radiative transfer options
      self % rttov_opts % rt_mw % fastem_version           = 6     !< FASTEM version (0-6); 0 => TESSEM2
      self % rttov_opts % rt_mw % supply_foam_fraction     = .false.  !< Supply a foam fraction to FASTEM
      self % rttov_opts % rt_mw % clw_data                 = .false.  !< Switch to enable input of cloud liquid water profile
      self % rttov_opts % rt_mw % clw_scheme               = 1     !< MW CLW scheme: 1 => Liebe, 2 => Rosenkranz, 3 => TKC
      self % rttov_opts % rt_mw % clw_calc_on_coef_lev     = .true.  !< Apply MW CLW calculations on coef/user levels (true/false resp.)
      self % rttov_opts % rt_mw % clw_cloud_top            = 322    !< Lower pressure limit for MW CLW calculations (hPa)
      self % rttov_opts % rt_mw % apply_band_correction    = .true.  !< Apply band-correction for Planck radiance and BT calculations

      self % rttov_opts % interpolation % addinterp        = .false.    !< Switch to enable RTTOV interpolator
      self % rttov_opts % interpolation % interp_mode      = interp_rochon !< Interpolation mode (1-5, see user guide)
      self % rttov_opts % interpolation % lgradp           = .false.    !< Switch to make pressure an active variable in TL/AD/K models
      self % rttov_opts % interpolation % spacetop         = .true.     !< Switch to assume space boundary at top-most input pressure level
      self % rttov_opts % interpolation % reg_limit_extrap = .false.    !< Switch to extrapolate input profiles using regression limits

      !> HTFRTC options structure
      self % rttov_opts % htfrtc_opts % htfrtc             = .false. !< Switch to use htfrtc
      self % rttov_opts % htfrtc_opts % n_pc_in            = -1   !< Number of principal components to be used
      self % rttov_opts % htfrtc_opts % reconstruct        = .false. !< Switch to select reconstructed radiances
      self % rttov_opts % htfrtc_opts % simple_cloud       = .false. !< Calculate simple cloud
      self % rttov_opts % htfrtc_opts % overcast           = .false. !< Calculate overcast cloud on all levels

    case ('OPS','PS45')
      self % rttov_opts % config % apply_reg_limits        = .true. ! set as true for all instruments
      self % rttov_opts % config % fix_hgpl                = .false. ! true @ PS44
      self % rttov_opts % config % verbose                 = .false. ! true if (ProcessMode > VerboseMode .OR. RTTOV_Verbosity > 0)
      self % rttov_opts % config % do_checkinput           = .false.  ! we will use the more thorough and verbose
      ! rttov_user_profile_checkinput instead to save ourselves an 
      ! RTTOV call where the profile is no good

      self % rttov_opts % rt_all % rad_down_lin_tau        = .true. !FALSE @PS44
      self % rttov_opts % rt_all % dtau_test               = .true. !FALSE @PS44

      if(trim(default_opts_set) == 'OPS') then
        self % rttov_opts % rt_all % addrefrac               = .false. ! At PS44 .TRUE. @ PS45?
      else
        self % rttov_opts % rt_all % addrefrac               = .true. ! At PS44 .TRUE. @ PS45?
      endif

      self % rttov_opts % rt_all % switchrad               = .true. 
      if(trim(default_opts_set) == 'OPS') then
        self % rttov_opts % rt_all % use_q2m                 = .false. ! At PS44 .TRUE. @ PS45?
      endif
      self % rttov_opts % rt_all % do_lambertian           = .false. !
      !  self % rttov_opts % rt_all % plane_parallel       = .FALSE. 

      self % rttov_opts % interpolation % addinterp        = .true. ! Allow interpolation of input profile
      self % rttov_opts % interpolation % interp_mode      = 4 ! Set interpolation method (4 for all insts at PS44)
      !  self % rttov_opts % interpolation % lgradp        = .FALSE. 
      if(trim(default_opts_set) == 'OPS') then
        self % rttov_opts % interpolation % spacetop         = .false. ! At PS44 .TRUE. @ PS45?
        self % rttov_opts % interpolation % reg_limit_extrap = .false. ! At PS44 .TRUE. @ PS45?
      else
        self % rttov_opts % interpolation % reg_limit_extrap = .true. ! At PS44 .TRUE. @ PS45?
      endif

      !  self % rttov_opts % rt_ir % addsolar              = .FALSE.
      !  self % rttov_opts % rt_ir % do_nlte_correction    = .FALSE.
      self % rttov_opts % rt_ir % ir_sea_emis_model        = 1 ! At PS44 2 @ PS45?

      !  self % rttov_opts % rt_ir % user_aer_opt_param    = .FALSE. ! user specifies aerosol scattering optical parameters
      !  self % rttov_opts % rt_ir % user_cld_opt_param    = .FALSE. ! user specifies cloud scattering optical parameters
      !  self % rttov_opts % rt_ir % ir_scatt_model        = 2      ! Chou-scaling
      self % rttov_opts % rt_ir % grid_box_avg_cloud       = .true. ! Assume grid-box average for cloud concentrations

      self % rttov_opts % rt_ir % ozone_data               = .true.  ! There WILL be an ozone profile
      self % rttov_opts % rt_ir % co2_data                 = .false.   ! We do not have profiles
      self % rttov_opts % rt_ir % n2o_data                 = .false.   ! for any other constituents
      self % rttov_opts % rt_ir % ch4_data                 = .false.
      self % rttov_opts % rt_ir % co_data                  = .false.
      self % rttov_opts % rt_ir % so2_data                 = .false.

      ! RTTOV v12 Discrete ordinates method (DOM) not implemented for now 

      self % rttov_opts % rt_mw % clw_calc_on_coef_lev     = .true. ! .FALSE. @PS44
      self % rttov_opts % rt_mw % clw_data                 = .true.  ! Set to true for allocation
      if(trim(default_opts_set) == 'OPS') then
        self % rttov_opts % rt_mw % fastem_version           = 2 ! at PS44 6 at PS45? 
      endif
      !self % rttov_opts % rt_mw % supply_foam_fraction    = RTTOV_Supply_Foam_Fraction

      ! IF (ANY (MwscattSwitch)) THEN
      !  RTTOV12_opts_scatt % config                       = self % rttov_opts % config
      !  RTTOV12_opts_scatt % use_q2m                      = self % rttov_opts % rt_all % use_q2m
      !  RTTOV12_opts_scatt % fastem_version               = self % rttov_opts % rt_mw % fastem_version
      !  RTTOV12_opts_scatt % supply_foam_fraction         = self % rttov_opts % rt_mw % supply_foam_fraction
      !  RTTOV12_opts_scatt % interp_mode                  = self % rttov_opts % interpolation % interp_mode
      !  RTTOV12_opts_scatt % reg_limit_extrap             = self % rttov_opts % interpolation % reg_limit_extrap
      !  RTTOV12_opts_scatt % rttov9_fastem_bug            = self % rttov_opts % compatibility % rttov9_fastem_bug
      ! END IF

    case ('test')
      self % rttov_opts % rt_ir % addsolar                 = .false. ! Do not include solar radiation
      self % rttov_opts % interpolation % addinterp        = .true.  ! Allow interpolation of input profile
      self % rttov_opts % interpolation % interp_mode      = 4       ! Set interpolation method
      self % rttov_opts % interpolation % reg_limit_extrap = .true.  ! reg_limit_extrap
      self % rttov_opts % rt_all % addrefrac               = .true.  ! Include refraction in path calc
      self % rttov_opts % rt_all % switchrad               = .true.  ! Include refraction in path calc
      self % rttov_opts % rt_all % dtau_test               = .false. 
      self % rttov_opts % rt_ir % addclouds                = .false. ! Don't include cloud effects
      self % rttov_opts % rt_ir % addaerosl                = .false. ! Don't include aerosol effects

      self % rttov_opts % rt_ir % ozone_data               = .false. ! Set the relevant flag to .TRUE.
      self % rttov_opts % rt_ir % co2_data                 = .false. !  when supplying a profile of the
      self % rttov_opts % rt_ir % n2o_data                 = .false. !  given trace gas (ensure the
      self % rttov_opts % rt_ir % ch4_data                 = .false. !  coef file supports the gas)
      self % rttov_opts % rt_ir % co_data                  = .false. !
      self % rttov_opts % rt_ir % so2_data                 = .false. !
      self % rttov_opts % rt_mw % clw_data                 = .false. !

      self % rttov_opts % config % verbose                 = .true. ! Enable printing of warnings
      self % rttov_opts % config % apply_reg_limits        = .true.
      self % rttov_opts % config % do_checkinput           = .false.
      self % rttov_opts % config % fix_hgpl                = .false.
    case('default')
      ! RTTOV defaults
    end select

  end subroutine set_defaults_rttov

  subroutine populate_hofxdiags(RTProf, chanprof, hofxdiags)
    use ufo_constants_mod, only : g_to_kg

    type(ufo_rttov_io),   intent(in)    :: RTProf
    type(rttov_chanprof), intent(in)    :: chanprof(:)
    type(ufo_geovals),    intent(inout) :: hofxdiags    !non-h(x) diagnostics

    integer                      :: jvar, chan, prof, ichan 
    integer                      :: nlayers, nchanprof, nlevels, nprofiles
    real(kind_real), allocatable :: od_level(:), wfunc(:)

    include 'rttov_calc_weighting_fn.interface'

    allocate(od_level(size(RTProf % transmission%tau_levels(:,1))))
    allocate(wfunc(size(RTProf % transmission%tau_levels(:,1))))

    missing = missing_value(missing)

    nchanprof = size(chanprof)
    nlevels = size(RTProf % profiles(1) % p)
    nprofiles = maxval(chanprof(:)%prof) !SIZE(RTProf % profiles)

    do jvar = 1, hofxdiags%nvar
      if (len(trim(hofxdiags%variables(jvar))) < 1) cycle

      !============================================
      ! Diagnostics used for QC and bias correction
      !============================================

      if (trim(xstr_diags(jvar)) == "") then
        ! forward h(x) diags
        select case(trim(ystr_diags(jvar)))

          ! variable: optical_thickness_of_atmosphere_layer_CH
          ! variable: transmittances_of_atmosphere_layer_CH        
          ! variable: weightingfunction_of_atmosphere_layer_CH
        case (var_opt_depth, var_lvl_transmit,var_lvl_weightfunc)

          nlayers = nlevels - 1
          hofxdiags%geovals(jvar)%nval = nlevels
          if(.not. allocated(hofxdiags%geovals(jvar)%vals)) &
             allocate(hofxdiags%geovals(jvar)%vals(hofxdiags%geovals(jvar)%nval,nprofiles))
          hofxdiags%geovals(jvar)%vals = missing
          ! get channel/profile
          do ichan = 1, nchanprof
            chan = chanprof(ichan)%chan
            prof = chanprof(ichan)%prof

            if(chan == ch_diags(jvar)) then
              ! if profile not skipped
              if(ystr_diags(jvar) == var_opt_depth) then
                od_level(:) = log(RTProf % transmission%tau_levels(:,chan)) !level->TOA transmittances -> od
                hofxdiags%geovals(jvar)%vals(:,prof) = od_level(1:nlevels-1) - od_level(2:nlevels) ! defined +ve 
              else if (ystr_diags(jvar) == var_lvl_transmit) then
                hofxdiags%geovals(jvar)%vals(:,prof) = RTProf % transmission % tau_levels(1:nlevels-1,chan) - RTProf % transmission%tau_levels(2:,chan)
              else if (ystr_diags(jvar) == var_lvl_weightfunc) then
                od_level(:) = log(RTProf % transmission%tau_levels(:,chan)) !level->TOA transmittances -> od
                call rttov_calc_weighting_fn(rttov_errorstatus, RTProf % profiles(prof)%p, od_level(:), &
                  hofxdiags%geovals(jvar)%vals(:,prof))

              endif
              !endif
            endif
          enddo

          ! variable: toa_outgoing_radiance_per_unit_wavenumber_CH [mW / (m^2 sr cm^-1)] (nval=1)
          ! variable: brightness_temperature_assuming_clear_sky_CH
          ! variable: pressure_level_at_peak_of_weightingfunction_CH
        case (var_radiance, var_tb_clr, var_tb, var_pmaxlev_weightfunc)
          ! always returned
          hofxdiags%geovals(jvar)%nval = 1
          if(.not. allocated(hofxdiags%geovals(jvar)%vals)) &
             allocate(hofxdiags%geovals(jvar)%vals(hofxdiags%geovals(jvar)%nval,nprofiles))
          hofxdiags%geovals(jvar)%vals = missing

          do ichan = 1, nchanprof
            chan = chanprof(ichan)%chan
            prof = chanprof(ichan)%prof

            if(chan == ch_diags(jvar)) then
              if(ystr_diags(jvar) == var_radiance) then
                hofxdiags%geovals(jvar)%vals(1,prof) = RTProf % radiance % total(ichan)
              else if(ystr_diags(jvar) == var_tb_clr) then
                hofxdiags%geovals(jvar)%vals(1,prof) = RTProf % radiance % bt_clear(ichan)
              else if(ystr_diags(jvar) == var_tb) then
                hofxdiags%geovals(jvar)%vals(1,prof) = RTProf % radiance % bt(ichan)
              else if(ystr_diags(jvar) == var_pmaxlev_weightfunc ) then
                call rttov_calc_weighting_fn(rttov_errorstatus, RTProf % profiles(prof)%p, od_level(:), &
                  Wfunc(:))
                hofxdiags%geovals(jvar)%vals(1,prof) = maxloc(Wfunc(:), DIM=1) ! scalar not array(1)
              end if
            endif
          end do

        case default
          write(message,*) 'ufo_radiancerttov_simobs: //&
            & ObsDiagnostic is unsupported, ', &
            & hofxdiags%variables(jvar)
          call fckit_log%info(message)
        end select

      else if (trim(ystr_diags(jvar)) == var_tb) then
        ! var_tb jacobians
        select case (trim(xstr_diags(jvar)))

        case (var_ts,var_mixr,var_q,var_clw,var_cli)
          nlayers = nlevels - 1
          hofxdiags%geovals(jvar)%nval = nlevels
          if(.not. allocated(hofxdiags%geovals(jvar)%vals)) &
            allocate(hofxdiags%geovals(jvar)%vals(hofxdiags%geovals(jvar)%nval,nprofiles))
          hofxdiags%geovals(jvar)%vals = missing

          do ichan = 1, nchanprof
            chan = chanprof(ichan)%chan
            prof = chanprof(ichan)%prof

            if(chan == ch_diags(jvar)) then
              if(xstr_diags(jvar) == var_ts) then
                hofxdiags%geovals(jvar)%vals(:,prof) = &
                  RTProf % profiles_k(ichan) % t(:)
              else if(xstr_diags(jvar) == var_mixr) then
                hofxdiags%geovals(jvar)%vals(:,prof) = &
                  RTProf % profiles_k(ichan) % q(:) / g_to_kg
              else if(xstr_diags(jvar) == var_q) then
                hofxdiags%geovals(jvar)%vals(:,prof) = &
                  RTProf % profiles_k(ichan) % q(:)
              else if(xstr_diags(jvar) == var_clw) then
                hofxdiags%geovals(jvar)%vals(:,prof) = &
                  RTProf % profiles_k(ichan) % clw(:)
              else if(xstr_diags(jvar) == var_cli) then
                ! not in use yet
              endif
            endif
          enddo

        case (var_sfc_t2m, var_sfc_tskin, var_sfc_emiss, var_sfc_q2m, var_sfc_p2m, var_u, var_v)
          hofxdiags%geovals(jvar)%nval = 1
          if(.not. allocated(hofxdiags%geovals(jvar)%vals)) &
            allocate(hofxdiags%geovals(jvar)%vals(hofxdiags%geovals(jvar)%nval,nprofiles))
          hofxdiags%geovals(jvar)%vals = missing

          do ichan = 1, nchanprof
            chan = chanprof(ichan)%chan
            prof = chanprof(ichan)%prof

            if(chan == ch_diags(jvar)) then
              if(xstr_diags(jvar) == var_sfc_tskin) then
                hofxdiags%geovals(jvar)%vals(1,prof) = &
                  RTProf % profiles_k(ichan) % skin % t
              else if (xstr_diags(jvar) == var_sfc_t2m) then
                hofxdiags%geovals(jvar)%vals(1,prof) = &
                  RTProf % profiles_k(ichan) % s2m % t
              else if (xstr_diags(jvar) == var_sfc_p2m) then
                hofxdiags%geovals(jvar)%vals(1,prof) = &
                  RTProf % profiles_k(ichan) % s2m % p
              else if (xstr_diags(jvar) == var_sfc_q2m) then
                hofxdiags%geovals(jvar)%vals(1,prof) = &
                  RTProf % profiles_k(ichan) % s2m % q
              else if (xstr_diags(jvar) == var_u) then
                hofxdiags%geovals(jvar)%vals(1,prof) = &
                  RTProf % profiles_k(ichan) % s2m % u
              else if (xstr_diags(jvar) == var_v) then
                hofxdiags%geovals(jvar)%vals(1,prof) = &
                  RTProf % profiles_k(ichan) % s2m % v
              else if (xstr_diags(jvar) == var_sfc_emiss) then
                hofxdiags%geovals(jvar)%vals(1,prof) = &
                  RTProf % emissivity_k(ichan) % emis_in
              end if
            end if
          end do

        case default
          write(message,*) 'ufo_radiancerttov_simobs: //&
            & ObsDiagnostic is unsupported, ', &
            & hofxdiags%variables(jvar)
          call fckit_log%info(message)
        end select
      else
        write(message,*) 'ufo_radiancerttov_simobs: //&
          & ObsDiagnostic is unsupported, ', &
          & hofxdiags%variables(jvar)
        call fckit_log%info(message)
      end if

    enddo

    deallocate(od_level,wfunc)

  end subroutine populate_hofxdiags

  subroutine parse_hofxdiags(hofxdiags, jacobian_needed)
    implicit none

    type(ufo_geovals),     intent(inout) :: hofxdiags    !non-h(x) diagnostics
    logical,               intent(out)   :: jacobian_needed
    character(10), parameter  :: jacobianstr = "_jacobian_"

    integer                   :: str_pos(4)
    character(len=maxvarlen) :: varstr
    integer                   :: jvar
    character(len=max_string) :: err_msg

    jacobian_needed = .false.

     if(hofxdiags%nvar > 0) then
       if (.not. allocated(ystr_diags)) allocate (ystr_diags(hofxdiags%nvar))
       if (.not. allocated(xstr_diags)) allocate (xstr_diags(hofxdiags%nvar))
       if (.not. allocated(ch_diags)) allocate (ch_diags(hofxdiags%nvar))

       ch_diags = -9999

       do jvar = 1, hofxdiags%nvar
         varstr = hofxdiags%variables(jvar)

         str_pos(4) = len_trim(varstr)
         if (str_pos(4) < 1) cycle
         str_pos(3) = index(varstr,"_",back=.true.)        !final "_" before channel
         read(varstr(str_pos(3)+1:str_pos(4)),*, err=999) ch_diags(jvar)
999      str_pos(1) = index(varstr,jacobianstr) - 1        !position before jacobianstr
         if (str_pos(1) == 0) then
           write(err_msg,*) 'parse_hofxdiags: _jacobian_ must be // &
             & preceded by dependent variable in config: ', &
             & hofxdiags%variables(jvar)
           call abor1_ftn(err_msg)
         else if (str_pos(1) > 0) then
           !Diagnostic is a Jacobian member (dy/dx)
           ystr_diags(jvar) = varstr(1:str_pos(1))
           str_pos(2) = str_pos(1) + len(jacobianstr) + 1 !begin xstr_diags
           jacobian_needed = .true.
           str_pos(4) = str_pos(3) - str_pos(2)
           xstr_diags(jvar)(1:str_pos(4)) = varstr(str_pos(2):str_pos(3)-1)
           xstr_diags(jvar)(str_pos(4)+1:) = ""
         else !null
           !Diagnostic is a dependent variable (y)
           
           xstr_diags(jvar) = ""
           ystr_diags(jvar)(1:str_pos(3)-1) = varstr(1:str_pos(3)-1)
           ystr_diags(jvar)(str_pos(3):) = ""
           if (ch_diags(jvar) < 0) ystr_diags(jvar) = varstr
         end if
       end do
     end if

  end subroutine parse_hofxdiags

end module ufo_radiancerttov_utils_mod
