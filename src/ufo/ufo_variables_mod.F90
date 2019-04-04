!
!  (C) Copyright 2017-2019 UCAR
!  
!  This software is licensed under the terms of the Apache Licence Version 2.0
!  which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
!

module ufo_vars_mod

implicit none
private
public :: ufo_vars, ufo_vars_setup, ufo_vars_clone, ufo_vars_delete
public :: ufo_vars_getindex, ufo_vars_nvars, ufo_vars_vnames

INTEGER, PARAMETER, PUBLIC :: n_aerosols_gocart_nasa=14,&
     &n_aerosols_gocart_esrl=15,n_aerosols_other=1

integer, parameter, public :: MAXVARLEN=56
character(len=MAXVARLEN), public, parameter :: var_tv   = "virtual_temperature"
character(len=MAXVARLEN), public, parameter :: var_ts   = "air_temperature"
character(len=MAXVARLEN), public, parameter :: var_t    = "temperature"
character(len=MAXVARLEN), public, parameter :: var_prsl = "atmosphere_ln_pressure_coordinate"
character(len=MAXVARLEN), public, parameter :: var_mixr = "humidity_mixing_ratio"
character(len=MAXVARLEN), public, parameter :: var_q    = "specific_humidity"
character(len=MAXVARLEN), public, parameter :: var_prs  = "air_pressure"
character(len=MAXVARLEN), public, parameter :: var_prsi = "air_pressure_levels"
character(len=MAXVARLEN), public, parameter :: var_z    = "geopotential_height"
character(len=MAXVARLEN), public, parameter :: var_zi   = "geopotential_height_levels"
character(len=MAXVARLEN), public, parameter :: var_sfc_z= "sfc_geopotential_height"
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

!@mzp strings have to be same MAXVARLEN length for array constructor
CHARACTER(len=MAXVARLEN), public, parameter :: var_rh          = "relative_humidity"
CHARACTER(len=MAXVARLEN), DIMENSION(n_aerosols_gocart_nasa), PUBLIC :: var_aerosols_gocart_nasa = [&
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
!    &var_aerosols_gocart_nasa,&
!    &"p25                                                     "]
! won't compile
CHARACTER(len=MAXVARLEN), DIMENSION(n_aerosols_gocart_esrl), PUBLIC :: var_aerosols_gocart_esrl = [&
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

CHARACTER(len=MAXVARLEN), DIMENSION(n_aerosols_other), PUBLIC :: var_aerosols_other = [&
     &"other                                                   "]

character(len=MAXVARLEN), public, parameter :: var_seaicefrac    = "sea_ice_category_area_fraction"
character(len=MAXVARLEN), public, parameter :: var_seaicethick   = "sea_ice_category_thickness"
character(len=MAXVARLEN), public, parameter :: var_abs_topo      = "sea_surface_height_above_geoid"
character(len=MAXVARLEN), public, parameter :: var_ocn_pot_temp  = "sea_water_potential_temperature"
character(len=MAXVARLEN), public, parameter :: var_ocn_con_temp  = "sea_water_conservative_temperature"
character(len=MAXVARLEN), public, parameter :: var_ocn_abs_salt  = "sea_water_absolute_salinity"
character(len=MAXVARLEN), public, parameter :: var_ocn_salt      = "sea_water_practical_salinity"
character(len=MAXVARLEN), public, parameter :: var_ocn_lay_thick = "sea_water_cell_thickness"
character(len=MAXVARLEN), public, parameter :: var_ocn_sst       = "sea_surface_temperature"

! ------------------------------------------------------------------------------

!> Fortran derived type to represent model variables
type :: ufo_vars
  integer :: nv
  character(len=MAXVARLEN), allocatable :: fldnames(:) !< Variable identifiers
end type ufo_vars

! ------------------------------------------------------------------------------
contains
! ------------------------------------------------------------------------------

subroutine ufo_vars_setup(self, c_vars)
use iso_c_binding
use config_mod
implicit none
type(ufo_vars), intent(inout) :: self
type(c_ptr), intent(in)       :: c_vars
character(len=30*MAXVARLEN) :: svars

self%nv = config_get_int(c_vars, "nvars")

allocate(self%fldnames(self%nv))
svars = config_get_string(c_vars,len(svars),"variables")
read(svars,*) self%fldnames

! TODO: a check on whether this var is in the list of defined vars

end subroutine ufo_vars_setup

! ------------------------------------------------------------------------------

subroutine ufo_vars_clone(self, other)
implicit none
type(ufo_vars), intent(in)    :: self
type(ufo_vars), intent(inout) :: other

call ufo_vars_delete(other)
other%nv = self%nv
allocate(other%fldnames(other%nv))
other%fldnames(:) = self%fldnames(:)

end subroutine ufo_vars_clone

! ------------------------------------------------------------------------------

subroutine ufo_vars_delete(self)
implicit none
type(ufo_vars), intent(inout) :: self

if (allocated(self%fldnames)) deallocate(self%fldnames)
self%nv = 0

end subroutine ufo_vars_delete

! ------------------------------------------------------------------------------

integer function ufo_vars_getindex(self, varname)
implicit none
type(ufo_vars), intent(in)       :: self
character(MAXVARLEN), intent(in) :: varname

integer :: ivar

ufo_vars_getindex = -1

do ivar = 1, self%nv
  if (self%fldnames(ivar) == varname) then
    ufo_vars_getindex = ivar
    exit
  endif
enddo

end function ufo_vars_getindex

! ------------------------------------------------------------------------------

integer function ufo_vars_nvars(self) 
implicit none
type(ufo_vars), intent(in) :: self

ufo_vars_nvars = self%nv

end function ufo_vars_nvars

! ------------------------------------------------------------------------------

function ufo_vars_vnames(self) 
implicit none
type(ufo_vars), intent(in) :: self

character(len=MAXVARLEN), dimension(self%nv) :: ufo_vars_vnames

ufo_vars_vnames(1:self%nv) = self%fldnames

end function ufo_vars_vnames

! ------------------------------------------------------------------------------

end module ufo_vars_mod
