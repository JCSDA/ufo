
!  (C) Copyright 2024 UCAR
!
!  This software is licensed under the terms of the Apache Licence Version 2.0
!  which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
!

module ufo_vars_mod

implicit none
private
public :: ufo_vars_read, ufo_vars_getindex

integer, parameter, public :: n_aerosols_gocart_default=14,&
     &n_aerosols_gocart_gefs=15,n_aerosols_gocart_ufs=18,&
     &n_aerosols_gocart_geos=18 !soon 21 for geos

integer, parameter, public :: MAXVARLEN=100
character(len=MAXVARLEN), public, parameter :: var_tv   = "virtual_temperature"
character(len=MAXVARLEN), public, parameter :: var_ts   = "air_temperature"
character(len=MAXVARLEN), public, parameter :: var_mixr = "humidity_mixing_ratio" ! g/kg
character(len=MAXVARLEN), public, parameter :: var_q    = "specific_humidity"     ! kg/kg
character(len=MAXVARLEN), public, parameter :: var_u    = "eastward_wind"
character(len=MAXVARLEN), public, parameter :: var_v    = "northward_wind"
character(len=MAXVARLEN), public, parameter :: var_prs  = "air_pressure"
character(len=MAXVARLEN), public, parameter :: var_prsi = "air_pressure_levels"
character(len=MAXVARLEN), public, parameter :: var_prsimo = "air_pressure_levels_minus_one"
character(len=MAXVARLEN), public, parameter :: var_delp   = "air_pressure_thickness"
character(len=MAXVARLEN), public, parameter :: var_ps     = "surface_pressure"
character(len=MAXVARLEN), public, parameter :: var_pmsl   = "surface_pressure_at_mean_sea_level" ! used by MetOffice opsinputs
character(len=MAXVARLEN), public, parameter :: var_z      = "geopotential_height"
character(len=MAXVARLEN), public, parameter :: var_zm     = "geometric_height"
character(len=MAXVARLEN), public, parameter :: var_zi     = "geopotential_height_levels"
character(len=MAXVARLEN), public, parameter :: var_zimo   = "geopotential_height_levels_minus_one"
character(len=MAXVARLEN), public, parameter :: var_sfc_z  = "surface_geopotential_height"
character(len=MAXVARLEN), public, parameter :: var_oz     = "mole_fraction_of_ozone_in_air"
character(len=MAXVARLEN), public, parameter :: var_co2    = "mole_fraction_of_carbon_dioxide_in_air"

! Directly predicted microphysics species mixing ratios (kg/kg) and number concentrations (#/kg)
character(len=MAXVARLEN), public, parameter :: var_qc     = "cloud_liquid_water"  ! liq_wat
character(len=MAXVARLEN), public, parameter :: var_qi     = "cloud_liquid_ice"    ! ice_wat
character(len=MAXVARLEN), public, parameter :: var_qr     = "rain_water"          ! rainwat
character(len=MAXVARLEN), public, parameter :: var_qs     = "snow_water"          ! snowwat
character(len=MAXVARLEN), public, parameter :: var_qg     = "graupel"             ! graupel
character(len=MAXVARLEN), public, parameter :: var_qh     = "hail"                ! hail
character(len=MAXVARLEN), public, parameter :: var_nc     = "cloud_droplet_number_concentration"  ! water_nc
character(len=MAXVARLEN), public, parameter :: var_ni     = "cloud_ice_number_concentration"      ! ice_nc
character(len=MAXVARLEN), public, parameter :: var_nr     = "rain_number_concentration"           ! rain_nc
character(len=MAXVARLEN), public, parameter :: var_ns     = "snow_number_concentration"           ! snow_nc
character(len=MAXVARLEN), public, parameter :: var_ng     = "graupel_number_concentration"        ! graupel_nc
character(len=MAXVARLEN), public, parameter :: var_nh     = "hail_number_concentration"           ! hail_nc
character(len=MAXVARLEN), public, parameter :: var_qvg    = "volume_mixing_ratio_of_graupel_in_air"
character(len=MAXVARLEN), public, parameter :: var_qvh    = "volume_mixing_ratio_of_hail_in_air"

! Derived liquid and ice water paths (kg/m2) computed from mixing ratio times delta-Z times air density.
! The use of "mass_content" names is consistent with the current CCPP standard, but Met Office
! use of these names before passing to RTTOV is actually mixing ratio of each species wrt moist air (kg/kg).
! This should be corrected in the future.
character(len=MAXVARLEN), public, parameter :: var_clw_wp    = "mass_content_of_cloud_liquid_water_in_atmosphere_layer"
character(len=MAXVARLEN), public, parameter :: var_cli_wp    = "mass_content_of_cloud_ice_in_atmosphere_layer"
character(len=MAXVARLEN), public, parameter :: var_clr_wp    = "mass_content_of_rain_in_atmosphere_layer"
character(len=MAXVARLEN), public, parameter :: var_cls_wp    = "mass_content_of_snow_in_atmosphere_layer"
character(len=MAXVARLEN), public, parameter :: var_clg_wp    = "mass_content_of_graupel_in_atmosphere_layer"
character(len=MAXVARLEN), public, parameter :: var_clh_wp    = "mass_content_of_hail_in_atmosphere_layer"

! Microphysics species mixing ratios wrt moist air (kg/kg), for use in RTTOV.
! See comment above: Met Office use of "mass_content" names for this variables is inconsistent with the CCPP standard
character(len=MAXVARLEN), public, parameter :: var_clw    = "mass_content_of_cloud_liquid_water_in_atmosphere_layer"
character(len=MAXVARLEN), public, parameter :: var_cli    = "mass_content_of_cloud_ice_in_atmosphere_layer"

character(len=MAXVARLEN), public, parameter :: var_clwefr = "effective_radius_of_cloud_liquid_water_particle"
character(len=MAXVARLEN), public, parameter :: var_cliefr = "effective_radius_of_cloud_ice_particle"
character(len=MAXVARLEN), public, parameter :: var_clrefr = "effective_radius_of_rain_particle"
character(len=MAXVARLEN), public, parameter :: var_clsefr = "effective_radius_of_snow_particle"
character(len=MAXVARLEN), public, parameter :: var_clgefr = "effective_radius_of_graupel_particle"
character(len=MAXVARLEN), public, parameter :: var_clhefr = "effective_radius_of_hail_particle"
character(len=MAXVARLEN), public, parameter :: var_cldfrac = "cloud_area_fraction_in_atmosphere_layer"
character(len=MAXVARLEN), public, parameter :: var_cldfrac_vol = "cloud_volume_fraction_in_atmosphere_layer"
character(len=MAXVARLEN), public, parameter :: var_sfc_q2m = "specific_humidity_at_two_meters_above_surface" ! (kg/kg)
character(len=MAXVARLEN), public, parameter :: var_sfc_t2m = "surface_temperature" ! (K)
character(len=MAXVARLEN), public, parameter :: var_sfc_tskin = "skin_temperature"  ! (K)
character(len=MAXVARLEN), public, parameter :: var_sfc_wfrac = "water_area_fraction"
character(len=MAXVARLEN), public, parameter :: var_sfc_lfrac = "land_area_fraction"
character(len=MAXVARLEN), public, parameter :: var_sfc_ifrac = "ice_area_fraction"
character(len=MAXVARLEN), public, parameter :: var_sfc_sfrac = "surface_snow_area_fraction"
character(len=MAXVARLEN), public, parameter :: var_sfc_wtmp  = "surface_temperature_where_sea"
character(len=MAXVARLEN), public, parameter :: var_sfc_ltmp  = "surface_temperature_where_land"
character(len=MAXVARLEN), public, parameter :: var_sfc_itmp  = "surface_temperature_where_ice"
character(len=MAXVARLEN), public, parameter :: var_sfc_stmp  = "surface_temperature_where_snow"
character(len=MAXVARLEN), public, parameter :: var_sfc_sdepth  = "surface_snow_thickness"
character(len=MAXVARLEN), public, parameter :: var_sfc_vegfrac = "vegetation_area_fraction"
character(len=MAXVARLEN), public, parameter :: var_sfc_wspeed  = "surface_wind_speed"
character(len=MAXVARLEN), public, parameter :: var_sfc_wdir    = "surface_wind_from_direction"
character(len=MAXVARLEN), public, parameter :: var_sfc_u10     = "uwind_at_10m"
character(len=MAXVARLEN), public, parameter :: var_sfc_v10     = "vwind_at_10m"
character(len=MAXVARLEN), public, parameter :: var_sfc_u       = "surface_eastward_wind"
character(len=MAXVARLEN), public, parameter :: var_sfc_v       = "surface_northward_wind"
character(len=MAXVARLEN), public, parameter :: var_sfc_lai     = "leaf_area_index"
character(len=MAXVARLEN), public, parameter :: var_sfc_soilm   = "volume_fraction_of_condensed_water_in_soil"
character(len=MAXVARLEN), public, parameter :: var_sfc_soilt   = "soil_temperature"
character(len=MAXVARLEN), public, parameter :: var_sfc_landtyp_npoess = "land_type_index_NPOESS"
character(len=MAXVARLEN), public, parameter :: var_sfc_landtyp_igbp   = "land_type_index_IGBP"
character(len=MAXVARLEN), public, parameter :: var_sfc_landtyp_usgs   = "land_type_index_USGS"
character(len=MAXVARLEN), public, parameter :: var_sfc_vegtyp  = "vegetation_type_index"
character(len=MAXVARLEN), public, parameter :: var_sfc_soiltyp = "soil_type"
character(len=MAXVARLEN), public, parameter :: var_geomz       = "height"
character(len=MAXVARLEN), public, parameter :: var_sfc_geomz   = "surface_altitude"
character(len=MAXVARLEN), public, parameter :: var_sfc_rough   = "surface_roughness_length"
character(len=MAXVARLEN), public, parameter :: var_sfc_t       = "surface_temperature"
character(len=MAXVARLEN), public, parameter :: var_sfc_fact10  = "wind_reduction_factor_at_10m"
character(len=MAXVARLEN), public, parameter :: var_observable_domain_mask = "observable_domain_mask"
character(len=MAXVARLEN), public, parameter :: var_sfc_emiss   = "surface_emissivity"
character(len=MAXVARLEN), public, parameter :: var_sfc_sss     = "sea_surface_salinity"
character(len=MAXVARLEN), public, parameter :: var_opt_depth   = "optical_thickness_of_atmosphere_layer"
character(len=MAXVARLEN), public, parameter :: var_radiance    = "toa_outgoing_radiance_per_unit_wavenumber"
character(len=MAXVARLEN), public, parameter :: var_albedo      = "albedo"
character(len=MAXVARLEN), public, parameter :: var_albedo_clr  = "albedo_assuming_clear_sky"
character(len=MAXVARLEN), public, parameter :: var_tb          = "brightness_temperature"
character(len=MAXVARLEN), public, parameter :: var_tb_clr      = "brightness_temperature_assuming_clear_sky"
character(len=MAXVARLEN), public, parameter :: var_tb_overcast = "brightness_temperature_from_atmosphere_layer_to_toa"
character(len=MAXVARLEN), public, parameter :: var_rad_refl      = "radar_reflectivity"
character(len=MAXVARLEN), public, parameter :: var_rad_refl_att  = "radar_reflectivity_attenuated"
character(len=MAXVARLEN), public, parameter :: var_total_transmit= "toa_total_transmittance"
character(len=MAXVARLEN), public, parameter :: var_lvl_transmit= "transmittances_of_atmosphere_layer"
character(len=MAXVARLEN), public, parameter :: var_lvl_weightfunc= "weightingfunction_of_atmosphere_layer"
character(len=MAXVARLEN), public, parameter :: var_pmaxlev_weightfunc= "pressure_level_at_peak_of_weightingfunction"
character(len=MAXVARLEN), public, parameter :: var_tsavg5      = "average_surface_temperature_within_field_of_view"
character(len=MAXVARLEN), public, parameter :: var_sea_fric_vel    = "friction_velocity_over_water"
character(len=MAXVARLEN), public, parameter :: var_obk_length      = "obukhov_length"
character(len=MAXVARLEN), public, parameter :: var_tropprs     = "tropopause_pressure"

character(len=MAXVARLEN), public, parameter :: var_refl        = "equivalent_reflectivity_factor"
character(len=MAXVARLEN), public, parameter :: var_w           = "upward_air_velocity"

character(len=MAXVARLEN), public, parameter :: var_rh          = "relative_humidity" ! dimensionless (0 <= RH <= 1)
character(len=MAXVARLEN), public, parameter :: var_surf_tau = "transmission_at_surface"
character(len=MAXVARLEN), public, parameter :: var_sfc_landmask   = "landmask"       ! 0 (sea), 1 (land)
character(len=MAXVARLEN), public, parameter :: var_sfc_seaicefrac = "seaice_fraction"

character(len=MAXVARLEN), public :: var_seaicefrac      = "sea_ice_category_area_fraction"
character(len=MAXVARLEN), public :: var_seaicethick     = "sea_ice_category_thickness"
character(len=MAXVARLEN), public :: var_seaicesnowthick = "sea_ice_category_snow_thickness"
character(len=MAXVARLEN), public :: var_ocn_chl         = "mass_concentration_of_chlorophyll_in_sea_water"
character(len=MAXVARLEN), public :: var_abs_topo        = "sea_surface_height_above_geoid"
character(len=MAXVARLEN), public :: var_ocn_pot_temp    = "sea_water_potential_temperature"
character(len=MAXVARLEN), public :: var_ocn_con_temp    = "sea_water_conservative_temperature"
character(len=MAXVARLEN), public :: var_ocn_abs_salt    = "sea_water_absolute_salinity"
character(len=MAXVARLEN), public :: var_ocn_pra_salt    = "sea_water_practical_salinity"
character(len=MAXVARLEN), public :: var_ocn_salt        = "sea_water_salinity"
character(len=MAXVARLEN), public :: var_ocn_lay_thick   = "sea_water_cell_thickness"
character(len=MAXVARLEN), public :: var_ocn_depth       = "sea_water_depth"
character(len=MAXVARLEN), public :: var_ocn_sst         = "sea_surface_temperature"
character(len=MAXVARLEN), public :: var_sea_td          = "sea_surface_foundation_temperature"
character(len=MAXVARLEN), public :: var_latent_vap      = "latent_heat_vaporization"
character(len=MAXVARLEN), public :: var_sw_rad          = "net_downwelling_shortwave_radiation"
character(len=MAXVARLEN), public :: var_latent_heat     = "upward_latent_heat_flux_in_air"
character(len=MAXVARLEN), public :: var_sens_heat       = "upward_sensible_heat_flux_in_air"
character(len=MAXVARLEN), public :: var_lw_rad          = "net_downwelling_longwave_radiation"

character(len=MAXVARLEN), public :: var_oz_thick        = "ozone_thickness"
character(len=MAXVARLEN), public :: var_water_vapor     = "water_vapor"
character(len=MAXVARLEN), public :: var_cld_tau         = "cloud_optical_thickness"
character(len=MAXVARLEN), public :: var_cld_lwp         = "cloud_liquid_water_path"
character(len=MAXVARLEN), public :: var_aerosol_tau     = "aerosol_optical_thickness"
character(len=MAXVARLEN), public :: var_scat_albedo     = "single_scattering_albedo"
character(len=MAXVARLEN), public :: var_asym_par        = "asymmetry_parameter"
character(len=MAXVARLEN), public :: var_carb_det        = "Carbon_nitrogen_detritus_concentration"
character(len=MAXVARLEN), public :: var_inorg_carb      = "Particulate_inorganic_carbon"
character(len=MAXVARLEN), public :: var_dis_carb        = "colored_dissolved_organic_carbon"
character(len=MAXVARLEN), public :: var_diatom_conc     = "diatom_concentration"
character(len=MAXVARLEN), public :: var_chloro_conc     = "chlorophyte_concentration"
character(len=MAXVARLEN), public :: var_cyano_conc      = "cyano-bacteria_concentration"
character(len=MAXVARLEN), public :: var_cocco_conc      = "coccolithophore_concentration"
character(len=MAXVARLEN), public :: var_dino_conc       = "dinoflagellate_concentration"
character(len=MAXVARLEN), public :: var_phaeo_conc      = "phaeocystis_concentration"

! GOCART quantities  for AODCRTM
character(len=MAXVARLEN), public, parameter :: var_du001 = "mass_fraction_of_dust001_in_air"
character(len=MAXVARLEN), public, parameter :: var_du002 = "mass_fraction_of_dust002_in_air"
character(len=MAXVARLEN), public, parameter :: var_du003 = "mass_fraction_of_dust003_in_air"
character(len=MAXVARLEN), public, parameter :: var_du004 = "mass_fraction_of_dust004_in_air"
character(len=MAXVARLEN), public, parameter :: var_du005 = "mass_fraction_of_dust005_in_air"
character(len=MAXVARLEN), public, parameter :: var_ss001 = "mass_fraction_of_sea_salt001_in_air"
character(len=MAXVARLEN), public, parameter :: var_ss002 = "mass_fraction_of_sea_salt002_in_air"
character(len=MAXVARLEN), public, parameter :: var_ss003 = "mass_fraction_of_sea_salt003_in_air"
character(len=MAXVARLEN), public, parameter :: var_ss004 = "mass_fraction_of_sea_salt004_in_air"
character(len=MAXVARLEN), public, parameter :: var_ss005 = "mass_fraction_of_sea_salt005_in_air"
character(len=MAXVARLEN), public, parameter :: var_bcphobic = "mass_fraction_of_hydrophobic_black_carbon_in_air"
character(len=MAXVARLEN), public, parameter :: var_bcphilic = "mass_fraction_of_hydrophilic_black_carbon_in_air"
character(len=MAXVARLEN), public, parameter :: var_ocphobic = "mass_fraction_of_hydrophobic_organic_carbon_in_air"
character(len=MAXVARLEN), public, parameter :: var_ocphilic = "mass_fraction_of_hydrophilic_organic_carbon_in_air"
character(len=MAXVARLEN), public, parameter :: var_brphobic = "mass_fraction_of_hydrophobic_brown_carbon_in_air"
character(len=MAXVARLEN), public, parameter :: var_brphilic = "mass_fraction_of_hydrophilic_brown_carbon_in_air"
character(len=MAXVARLEN), public, parameter :: var_sulfate = "mass_fraction_of_sulfate_in_air"
character(len=MAXVARLEN), public, parameter :: var_no3an1 = "mass_fraction_of_nitrate001_in_air"
character(len=MAXVARLEN), public, parameter :: var_no3an2 = "mass_fraction_of_nitrate002_in_air"
character(len=MAXVARLEN), public, parameter :: var_no3an3 = "mass_fraction_of_nitrate003_in_air"

! Values for AODExt
character(len=MAXVARLEN), public, parameter :: var_ext1 = "volume_extinction_in_air_due_to_aerosol_particles_lambda1"
character(len=MAXVARLEN), public, parameter :: var_ext2 = "volume_extinction_in_air_due_to_aerosol_particles_lambda2"
character(len=MAXVARLEN), public, parameter :: var_ext3 = "volume_extinction_in_air_due_to_aerosol_particles_lambda3"
character(len=MAXVARLEN), public, parameter :: var_airdens = "moist_air_density"

! Scaling factors for deriving PM2.5 from CMAQ aerosols in the Aitken (at), accumulation (ac) and coarse (co) modes
character(len=MAXVARLEN), public, parameter :: var_pm25at = "pm25at"
character(len=MAXVARLEN), public, parameter :: var_pm25ac = "pm25ac"
character(len=MAXVARLEN), public, parameter :: var_pm25co = "pm25co"

character(len=MAXVARLEN), dimension(n_aerosols_gocart_default), public, parameter  :: &
     &var_aerosols_gocart_default = [&
     &var_sulfate,&
     &var_bcphobic, var_bcphilic, var_ocphobic, var_ocphilic,&
     &var_du001, var_du002, var_du003, var_du004, var_du005,&
     &var_ss001, var_ss002, var_ss003, var_ss004]

character(len=maxvarlen), dimension(n_aerosols_gocart_gefs), public, parameter :: &
     &var_aerosols_gocart_gefs = [&
     &var_sulfate,&
     &var_bcphobic, var_bcphilic, var_ocphobic, var_ocphilic,&
     &var_du001, var_du002, var_du003, var_du004, var_du005,&
     &var_ss001, var_ss002, var_ss003, var_ss004, var_ss005]

character(len=maxvarlen), dimension(n_aerosols_gocart_ufs), public, parameter :: &
     &var_aerosols_gocart_ufs = [&
     &var_sulfate,&
     &var_bcphobic, var_bcphilic, var_ocphobic, var_ocphilic,& 
     &var_du001, var_du002, var_du003, var_du004, var_du005,&
     &var_ss001, var_ss002, var_ss003, var_ss004, var_ss005,&
     &var_no3an1, var_no3an2, var_no3an3]

character(len=maxvarlen), dimension(n_aerosols_gocart_geos), public, parameter :: &
     &var_aerosols_gocart_geos = [&
     &var_sulfate,& !two bins for sulfate?
     &var_bcphobic, var_bcphilic, var_ocphobic, var_ocphilic,& ! var_brphobic, var_brphilic to be added soon
     &var_du001, var_du002, var_du003, var_du004, var_du005,&
     &var_ss001, var_ss002, var_ss003, var_ss004, var_ss005,&
     &var_no3an1, var_no3an2, var_no3an3]

! ------------------------------------------------------------------------------
contains

subroutine ufo_vars_read(f_vars, vars)
use fckit_configuration_module, only: fckit_configuration
implicit none
type(fckit_configuration), intent(in)                              :: f_vars
character(len=MAXVARLEN), dimension(:), allocatable, intent(inout) :: vars

integer :: nvars
character(len=:), allocatable :: str

if (f_vars%has("nvars")) then
  call f_vars%get_or_die("nvars",nvars)
  if (allocated(vars)) deallocate(vars)
  allocate(vars(nvars))
  call f_vars%get_or_die("variables",str)
  read(str,*) vars
else
  allocate(vars(0))
endif

end subroutine ufo_vars_read

! ------------------------------------------------------------------------------

integer function ufo_vars_getindex(vars, varname)
use ufo_utils_mod, only: cmp_strings
implicit none
character(len=*), intent(in) :: vars(:)
character(len=*), intent(in) :: varname

integer :: ivar

ufo_vars_getindex = -1

do ivar = 1, size(vars)
  if (cmp_strings(vars(ivar), varname)) then
    ufo_vars_getindex = ivar
    exit
  endif
enddo

end function ufo_vars_getindex

! ------------------------------------------------------------------------------

end module ufo_vars_mod
