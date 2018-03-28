!
!  (C) Copyright 2017 UCAR
!  
!  This software is licensed under the terms of the Apache Licence Version 2.0
!  which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
!

module ufo_vars_mod

implicit none
private
public :: ufo_vars, ufo_vars_setup, ufo_vars_clone, ufo_vars_delete
public :: ufo_vars_getindex, ufo_vars_nvars

integer, parameter, public :: MAXVARLEN=56
character(len=MAXVARLEN), public :: var_tv   = "virtual_temperature"
character(len=MAXVARLEN), public :: var_prsl = "atmosphere_ln_pressure_coordinate"
character(len=MAXVARLEN), public :: var_mixr = "humidity_mixing_ratio"
character(len=MAXVARLEN), public :: var_prs  = "air_pressure"
character(len=MAXVARLEN), public :: var_prsi = "air_pressure_levels"
character(len=MAXVARLEN), public :: var_oz   = "mass_concentration_of_ozone_in_air"
character(len=MAXVARLEN), public :: var_co2  = "mass_concentration_of_carbon_dioxide_in_air"
character(len=MAXVARLEN), public :: var_clw  = "atmosphere_mass_content_of_cloud_liquid_water"
character(len=MAXVARLEN), public :: var_cli  = "atmosphere_mass_content_of_cloud_ice"
character(len=MAXVARLEN), public :: var_clwefr = "effective_radius_of_cloud_liquid_water_particle"
character(len=MAXVARLEN), public :: var_cliefr = "effective_radius_of_cloud_ice_particle"
character(len=MAXVARLEN), public :: var_sfc_wfrac = "Water_Fraction"
character(len=MAXVARLEN), public :: var_sfc_lfrac = "Land_Fraction"
character(len=MAXVARLEN), public :: var_sfc_ifrac = "Ice_Fraction"
character(len=MAXVARLEN), public :: var_sfc_sfrac = "Snow_Fraction"
character(len=MAXVARLEN), public :: var_sfc_wtmp  = "Water_Temperature"
character(len=MAXVARLEN), public :: var_sfc_ltmp  = "Land_Temperature"
character(len=MAXVARLEN), public :: var_sfc_itmp  = "Ice_Temperature"
character(len=MAXVARLEN), public :: var_sfc_stmp  = "Snow_Temperature"
character(len=MAXVARLEN), public :: var_sfc_sdepth  = "Snow_Depth"
character(len=MAXVARLEN), public :: var_sfc_vegfrac = "Vegetation_Fraction"
character(len=MAXVARLEN), public :: var_sfc_wspeed  = "Sfc_Wind_Speed"
character(len=MAXVARLEN), public :: var_sfc_wdir    = "Sfc_Wind_Direction"
character(len=MAXVARLEN), public :: var_sfc_lai     = "Lai"
character(len=MAXVARLEN), public :: var_sfc_soilm   = "Soil_Moisture"
character(len=MAXVARLEN), public :: var_sfc_soilt   = "Soil_Temperature"
character(len=MAXVARLEN), public :: var_sfc_landtyp = "Land_Type_Index"
character(len=MAXVARLEN), public :: var_sfc_vegtyp  = "Vegetation_Type"
character(len=MAXVARLEN), public :: var_sfc_soiltyp = "Soil_Type"
character(len=MAXVARLEN), public :: var_seaicefrac    = "ice_concentration"
character(len=MAXVARLEN), public :: var_stericheight  = "steric_height"
character(len=MAXVARLEN), public :: var_seaicethick   = "ice_thickness"
character(len=MAXVARLEN), public :: var_abs_topo      = "sea_surface_height_above_geoid"
character(len=MAXVARLEN), public :: var_ocn_pot_temp  = "ocean_potential_temperature"
character(len=MAXVARLEN), public :: var_ocn_con_temp  = "ocean_conservative_temperature"
character(len=MAXVARLEN), public :: var_ocn_abs_salt  = "ocean_absolute_salinity"
character(len=MAXVARLEN), public :: var_ocn_salt      = "ocean_salinity"
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
   print *,'--------------------',self%fldnames(ivar),varname
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

end module ufo_vars_mod
