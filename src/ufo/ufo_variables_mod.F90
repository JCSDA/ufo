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

INTEGER, PARAMETER, PUBLIC :: n_aerosols_gocart_default=14,&
     &n_aerosols_gocart_esrl=15,n_aerosols_other=1

integer, parameter, public :: MAXVARLEN=56
character(len=MAXVARLEN), public, parameter :: var_tv   = "virtual_temperature"
character(len=MAXVARLEN), public, parameter :: var_ts   = "air_temperature"
character(len=MAXVARLEN), public, parameter :: var_t    = "temperature"
character(len=MAXVARLEN), public, parameter :: var_prsl = "atmosphere_ln_pressure_coordinate"
character(len=MAXVARLEN), public, parameter :: var_mixr = "humidity_mixing_ratio"
character(len=MAXVARLEN), public, parameter :: var_q    = "specific_humidity"
character(len=MAXVARLEN), public, parameter :: var_u    = "eastward_wind"
character(len=MAXVARLEN), public, parameter :: var_v    = "northward_wind"
character(len=MAXVARLEN), public, parameter :: var_prs  = "air_pressure"
character(len=MAXVARLEN), public, parameter :: var_prsi = "air_pressure_levels"
character(len=MAXVARLEN), public, parameter :: var_ps   = "surface_pressure"
character(len=MAXVARLEN), public, parameter :: var_z    = "geopotential_height"
character(len=MAXVARLEN), public, parameter :: var_zi   = "geopotential_height_levels"
character(len=MAXVARLEN), public, parameter :: var_sfc_z= "surface_geopotential_height"
character(len=MAXVARLEN), public, parameter :: var_oz   = "mass_concentration_of_ozone_in_air"
character(len=MAXVARLEN), public, parameter :: var_co2  = "mass_concentration_of_carbon_dioxide_in_air"
character(len=MAXVARLEN), public, parameter :: var_clw  = "atmosphere_mass_content_of_cloud_liquid_water"
character(len=MAXVARLEN), public, parameter :: var_cli  = "atmosphere_mass_content_of_cloud_ice"
character(len=MAXVARLEN), public, parameter :: var_clwefr = "effective_radius_of_cloud_liquid_water_particle"
character(len=MAXVARLEN), public, parameter :: var_cliefr = "effective_radius_of_cloud_ice_particle"
character(len=MAXVARLEN), public, parameter :: var_sfc_wfrac = "Water_Fraction"
character(len=MAXVARLEN), public, parameter :: var_sfc_lfrac = "Land_Fraction"
character(len=MAXVARLEN), public, parameter :: var_sfc_ifrac = "Ice_Fraction"
character(len=MAXVARLEN), public, parameter :: var_sfc_sfrac = "Snow_Fraction"
character(len=MAXVARLEN), public, parameter :: var_sfc_wtmp  = "Water_Temperature"
character(len=MAXVARLEN), public, parameter :: var_sfc_ltmp  = "Land_Temperature"
character(len=MAXVARLEN), public, parameter :: var_sfc_itmp  = "Ice_Temperature"
character(len=MAXVARLEN), public, parameter :: var_sfc_stmp  = "Snow_Temperature"
character(len=MAXVARLEN), public, parameter :: var_sfc_sdepth  = "Snow_Depth"
character(len=MAXVARLEN), public, parameter :: var_sfc_vegfrac = "Vegetation_Fraction"
character(len=MAXVARLEN), public, parameter :: var_sfc_wspeed  = "Sfc_Wind_Speed"
character(len=MAXVARLEN), public, parameter :: var_sfc_wdir    = "Sfc_Wind_Direction"
character(len=MAXVARLEN), public, parameter :: var_sfc_lai     = "Lai"
character(len=MAXVARLEN), public, parameter :: var_sfc_soilm   = "Soil_Moisture"
character(len=MAXVARLEN), public, parameter :: var_sfc_soilt   = "Soil_Temperature"
character(len=MAXVARLEN), public, parameter :: var_sfc_landtyp = "Land_Type_Index"
character(len=MAXVARLEN), public, parameter :: var_sfc_vegtyp  = "Vegetation_Type"
character(len=MAXVARLEN), public, parameter :: var_sfc_soiltyp = "Soil_Type"
character(len=MAXVARLEN), public, parameter :: var_sfc_rough   = "surface_roughness_length"
character(len=MAXVARLEN), public, parameter :: var_sfc_t   = "surface_temperature"

!@mzp strings have to be same MAXVARLEN length for array constructor
CHARACTER(len=MAXVARLEN), public, parameter :: var_rh          = "relative_humidity"
CHARACTER(len=MAXVARLEN), DIMENSION(n_aerosols_gocart_default), PUBLIC, PARAMETER  :: &
     &var_aerosols_gocart_default = [&
     &"sulf                                                    ",&
     &"bc1                                                     ",&
     &"bc2                                                     ",&
     &"oc1                                                     ",&
     &"oc2                                                     ",&
     &"dust1                                                   ",&
     &"dust2                                                   ",&
     &"dust3                                                   ",&
     &"dust4                                                   ",&
     &"dust5                                                   ",&
     &"seas1                                                   ",&
     &"seas2                                                   ",&
     &"seas3                                                   ",&
     &"seas4                                                   "]
!@mzp var_aerosols_gocart_esrl =[&
!    &var_aerosols_gocart_default,&
!    &"p25                                                     "]
! won't compile
CHARACTER(len=MAXVARLEN), DIMENSION(n_aerosols_gocart_esrl), PUBLIC, PARAMETER :: &
     &var_aerosols_gocart_esrl = [&
     &"sulf                                                    ",&
     &"bc1                                                     ",&
     &"bc2                                                     ",&
     &"oc1                                                     ",&
     &"oc2                                                     ",&
     &"dust1                                                   ",&
     &"dust2                                                   ",&
     &"dust3                                                   ",&
     &"dust4                                                   ",&
     &"dust5                                                   ",&
     &"seas1                                                   ",&
     &"seas2                                                   ",&
     &"seas3                                                   ",&
     &"seas4                                                   ",&
     &"p25                                                     "]

CHARACTER(len=MAXVARLEN), DIMENSION(n_aerosols_other), PUBLIC, PARAMETER :: &
     &var_aerosols_other = [&
     &"other                                                   "]

character(len=MAXVARLEN), public :: var_seaicefrac    = "sea_ice_category_area_fraction"
character(len=MAXVARLEN), public :: var_seaicethick   = "sea_ice_category_thickness"
character(len=MAXVARLEN), public :: var_abs_topo      = "sea_surface_height_above_geoid"
character(len=MAXVARLEN), public :: var_ocn_pot_temp  = "sea_water_potential_temperature"
character(len=MAXVARLEN), public :: var_ocn_con_temp  = "sea_water_conservative_temperature"
character(len=MAXVARLEN), public :: var_ocn_abs_salt  = "sea_water_absolute_salinity"
character(len=MAXVARLEN), public :: var_ocn_pra_salt      = "sea_water_practical_salinity"
character(len=MAXVARLEN), public :: var_ocn_salt      = "sea_water_salinity"
character(len=MAXVARLEN), public :: var_ocn_lay_thick = "sea_water_cell_thickness"
character(len=MAXVARLEN), public :: var_ocn_sst       = "sea_surface_temperature"
character(len=MAXVARLEN), public :: var_sea_td        = "sea_surface_foundation_temperature"
character(len=MAXVARLEN), public :: var_latent_vap    = "latent_heat_vaporization"
character(len=MAXVARLEN), public :: var_sw_rad        = "net_downwelling_shortwave_radiation"
character(len=MAXVARLEN), public :: var_latent_heat   = "upward_latent_heat_flux_in_air"
character(len=MAXVARLEN), public :: var_sens_heat     = "upward_sensible_heat_flux_in_air"
character(len=MAXVARLEN), public :: var_lw_rad        = "net_downwelling_longwave_radiation"
character(len=MAXVARLEN), public :: var_sea_fric_vel  = "friction_velocity_over_water"

! ------------------------------------------------------------------------------
contains

subroutine ufo_vars_read(c_vars, vars)
use iso_c_binding
use config_mod
implicit none
type(c_ptr), intent(in)       :: c_vars
character(len=MAXVARLEN), dimension(:), allocatable, intent(inout) :: vars
character(len=30*MAXVARLEN) :: svars

integer :: nvars

nvars = config_get_int(c_vars, "nvars")

if (allocated(vars)) deallocate(vars)
allocate(vars(nvars))
svars = config_get_string(c_vars,len(svars),"variables")
read(svars,*) vars

end subroutine ufo_vars_read

! ------------------------------------------------------------------------------

integer function ufo_vars_getindex(vars, varname)
implicit none
character(len=MAXVARLEN), intent(in) :: vars(:)
character(MAXVARLEN), intent(in)     :: varname

integer :: ivar

ufo_vars_getindex = -1

do ivar = 1, size(vars)
  if (vars(ivar) == varname) then
    ufo_vars_getindex = ivar
    exit
  endif
enddo

end function ufo_vars_getindex

! ------------------------------------------------------------------------------

end module ufo_vars_mod
