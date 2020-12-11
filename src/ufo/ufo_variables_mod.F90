!
!  (C) Copyright 2017-2019 UCAR
!
!  This software is licensed under the terms of the Apache Licence Version 2.0
!  which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
!

module ufo_vars_mod

implicit none
private
public :: ufo_vars_read, ufo_vars_getindex

integer, parameter, public :: n_aerosols_gocart_default=14,&
     &n_aerosols_gocart_merra_2=15,n_aerosols_other=1

integer, parameter, public :: MAXVARLEN=100
character(len=MAXVARLEN), public, parameter :: var_tv   = "virtual_temperature"
character(len=MAXVARLEN), public, parameter :: var_ts   = "air_temperature"
character(len=MAXVARLEN), public, parameter :: var_t    = "temperature"
character(len=MAXVARLEN), public, parameter :: var_mixr = "humidity_mixing_ratio" ! g/kg
character(len=MAXVARLEN), public, parameter :: var_q    = "specific_humidity"     ! kg/kg
character(len=MAXVARLEN), public, parameter :: var_u    = "eastward_wind"
character(len=MAXVARLEN), public, parameter :: var_v    = "northward_wind"
character(len=MAXVARLEN), public, parameter :: var_prs  = "air_pressure"
character(len=MAXVARLEN), public, parameter :: var_prsi = "air_pressure_levels"
character(len=MAXVARLEN), public, parameter :: var_delp = "air_pressure_thickness"
character(len=MAXVARLEN), public, parameter :: var_ps   = "surface_pressure"
character(len=MAXVARLEN), public, parameter :: var_z    = "geopotential_height"
character(len=MAXVARLEN), public, parameter :: var_zi   = "geopotential_height_levels"
character(len=MAXVARLEN), public, parameter :: var_sfc_z= "surface_geopotential_height"
character(len=MAXVARLEN), public, parameter :: var_oz   = "mole_fraction_of_ozone_in_air"
character(len=MAXVARLEN), public, parameter :: var_co2  = "mole_fraction_of_carbon_dioxide_in_air"
character(len=MAXVARLEN), public, parameter :: var_clw  = "mass_content_of_cloud_liquid_water_in_atmosphere_layer"
character(len=MAXVARLEN), public, parameter :: var_cli  = "mass_content_of_cloud_ice_in_atmosphere_layer"
character(len=MAXVARLEN), public, parameter :: var_clr  = "mass_content_of_rain_in_atmosphere_layer"
character(len=MAXVARLEN), public, parameter :: var_cls  = "mass_content_of_snow_in_atmosphere_layer"
character(len=MAXVARLEN), public, parameter :: var_clg  = "mass_content_of_graupel_in_atmosphere_layer"
character(len=MAXVARLEN), public, parameter :: var_clh  = "mass_content_of_hail_in_atmosphere_layer"
character(len=MAXVARLEN), public, parameter :: var_clwefr = "effective_radius_of_cloud_liquid_water_particle"
character(len=MAXVARLEN), public, parameter :: var_cliefr = "effective_radius_of_cloud_ice_particle"
character(len=MAXVARLEN), public, parameter :: var_clrefr = "effective_radius_of_rain_particle"
character(len=MAXVARLEN), public, parameter :: var_clsefr = "effective_radius_of_snow_particle"
character(len=MAXVARLEN), public, parameter :: var_clgefr = "effective_radius_of_graupel_particle"
character(len=MAXVARLEN), public, parameter :: var_clhefr = "effective_radius_of_hail_particle"
character(len=MAXVARLEN), public, parameter :: var_cldfrac= "cloud_area_fraction_in_atmosphere_layer"
character(len=MAXVARLEN), public, parameter :: var_sfc_p2m = "air_pressure_at_two_meters_above_surface"      ! (Pa)
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
character(len=MAXVARLEN), public, parameter :: var_sfc_lai     = "leaf_area_index"
character(len=MAXVARLEN), public, parameter :: var_sfc_soilm   = "volume_fraction_of_condensed_water_in_soil"
character(len=MAXVARLEN), public, parameter :: var_sfc_soilt   = "soil_temperature"
character(len=MAXVARLEN), public, parameter :: var_sfc_landtyp = "land_type_index"
character(len=MAXVARLEN), public, parameter :: var_sfc_vegtyp  = "vegetation_type_index"
character(len=MAXVARLEN), public, parameter :: var_sfc_soiltyp = "soil_type"
character(len=MAXVARLEN), public, parameter :: var_geomz       = "height"
character(len=MAXVARLEN), public, parameter :: var_sfc_geomz   = "surface_altitude"
character(len=MAXVARLEN), public, parameter :: var_sfc_rough   = "surface_roughness_length"
character(len=MAXVARLEN), public, parameter :: var_sfc_t       = "surface_temperature"
character(len=MAXVARLEN), public, parameter :: var_sfc_fact10  = "wind_reduction_factor_at_10m"
character(len=MAXVARLEN), public, parameter :: var_sfc_emiss   = "surface_emissivity"
character(len=MAXVARLEN), public, parameter :: var_sfc_sss     = "sea_surface_salinity"
character(len=MAXVARLEN), public, parameter :: var_opt_depth   = "optical_thickness_of_atmosphere_layer"
character(len=MAXVARLEN), public, parameter :: var_radiance    = "toa_outgoing_radiance_per_unit_wavenumber"
character(len=MAXVARLEN), public, parameter :: var_tb          = "brightness_temperature"
character(len=MAXVARLEN), public, parameter :: var_tb_clr      = "brightness_temperature_assuming_clear_sky"
character(len=MAXVARLEN), public, parameter :: var_lvl_transmit= "transmittances_of_atmosphere_layer"
character(len=MAXVARLEN), public, parameter :: var_lvl_weightfunc= "weightingfunction_of_atmosphere_layer"
character(len=MAXVARLEN), public, parameter :: var_pmaxlev_weightfunc= "pressure_level_at_peak_of_weightingfunction"
character(len=MAXVARLEN), public, parameter :: var_tsavg5      = "average_surface_temperature_within_field_of_view"


character(len=MAXVARLEN), public, parameter :: var_refl        = "equivalent_reflectivity_factor"
character(len=MAXVARLEN), public, parameter :: var_w           = "upward_air_velocity"

character(len=MAXVARLEN), public, parameter :: var_rh          = "relative_humidity" ! dimensionless (0 <= RH <= 1)
character(len=MAXVARLEN), public, parameter :: var_water_type_rttov = "water_type"   ! 0 (fresh), 1 (sea)
character(len=MAXVARLEN), public, parameter :: var_surf_type_rttov = "surface_type"  ! 0 (land), 1 (water), 2 (sea-ice)


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
character(len=MAXVARLEN), public :: var_ocn_sst         = "sea_surface_temperature"
character(len=MAXVARLEN), public :: var_sea_td          = "sea_surface_foundation_temperature"
character(len=MAXVARLEN), public :: var_latent_vap      = "latent_heat_vaporization"
character(len=MAXVARLEN), public :: var_sw_rad          = "net_downwelling_shortwave_radiation"
character(len=MAXVARLEN), public :: var_latent_heat     = "upward_latent_heat_flux_in_air"
character(len=MAXVARLEN), public :: var_sens_heat       = "upward_sensible_heat_flux_in_air"
character(len=MAXVARLEN), public :: var_lw_rad          = "net_downwelling_longwave_radiation"
character(len=MAXVARLEN), public :: var_sea_fric_vel    = "friction_velocity_over_water"

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
character(len=MAXVARLEN), public, parameter :: var_sulfate = "mass_fraction_of_sulfate_in_air"
character(len=MAXVARLEN), public, parameter :: var_no3an1 = "mass_fraction_of_nitrate001_in_air"
character(len=MAXVARLEN), public, parameter :: var_no3an2 = "mass_fraction_of_nitrate002_in_air"
character(len=MAXVARLEN), public, parameter :: var_no3an3 = "mass_fraction_of_nitrate003_in_air"

character(len=MAXVARLEN), dimension(n_aerosols_gocart_default), public, parameter  :: &
     &var_aerosols_gocart_default = [&
     &var_sulfate,&
     &var_bcphobic, var_bcphilic, var_ocphobic, var_ocphilic,&
     &var_du001, var_du002, var_du003, var_du004, var_du005,&
     &var_ss001, var_ss002, var_ss003, var_ss004]

character(len=maxvarlen), dimension(n_aerosols_gocart_merra_2), public, parameter :: &
     &var_aerosols_gocart_merra_2 = [&
     &var_sulfate,&
     &var_bcphobic, var_bcphilic, var_ocphobic, var_ocphilic,&
     &var_du001, var_du002, var_du003, var_du004, var_du005,&
     &var_ss001, var_ss002, var_ss003, var_ss004, var_ss005]

character(len=MAXVARLEN), dimension(n_aerosols_other), public, parameter :: &
     &var_aerosols_other = [&
     &"other                                                   "]

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
implicit none
character(len=*), intent(in) :: vars(:)
character(len=*), intent(in) :: varname

integer :: ivar

ufo_vars_getindex = -1

do ivar = 1, size(vars)
  if (trim(vars(ivar)) == trim(varname)) then
    ufo_vars_getindex = ivar
    exit
  endif
enddo

end function ufo_vars_getindex

! ------------------------------------------------------------------------------

end module ufo_vars_mod
