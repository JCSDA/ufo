! (C) Copyright 2017-2018 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

!> Fortran module for atmsfcinterp observation operator

module ufo_atmsfcinterp_mod

  use iso_c_binding
  use config_mod
  use kinds

  use ufo_geovals_mod, only: ufo_geovals, ufo_geoval, ufo_geovals_get_var
  use ufo_basis_mod, only: ufo_basis
  use ufo_vars_mod
  use obsspace_mod

  implicit none
  private

  integer, parameter :: max_string = 800
  real(kind_real), parameter :: grav = 9.80665e+0_kind_real
  !> Fortran derived type for the observation type
  type, extends(ufo_basis), public :: ufo_atmsfcinterp
  private
    integer :: nvars
    real(kind_real) :: magl
    character(len=max_string), public, allocatable :: varin(:)
    character(len=max_string), public, allocatable :: varout(:)
  contains
    procedure :: setup  => ufo_atmsfcinterp_setup
    procedure :: delete => ufo_atmsfcinterp_delete
    procedure :: simobs => ufo_atmsfcinterp_simobs
  end type ufo_atmsfcinterp

contains

! ------------------------------------------------------------------------------
subroutine ufo_atmsfcinterp_setup(self, c_conf)
  implicit none
  class(ufo_atmsfcinterp), intent(inout) :: self
  type(c_ptr),        intent(in)    :: c_conf
  integer :: ii, ivar

  !> Size of variables
  self%nvars = size(config_get_string_vector(c_conf, max_string, "variables"))
  !> Allocate varout: variables in the observation vector
  allocate(self%varout(self%nvars))
  !> Read variable list and store in varout
  self%varout = config_get_string_vector(c_conf, max_string, "variables")
  !> Allocate varin: variables we need from the model
  allocate(self%varin(11))
  !> Set vars_in based on vars_out
  !do ii = 1, self%nvars
  ! self%varin(ii) = self%varout(ii)
  !enddo
  ! TODO: additional variables will need to be added but all of the 'standard'
  ! vars need the same input from the GeoVaLs (T,P,U,V,Q need all for each other)
  !> add geopotential height
  self%varin(1) = "geopotential_height"
  !> need skin temperature for near-surface interpolations
  self%varin(2) = "sfc_skin_temperature"
  !> need surface geopotential height to get difference from phi
  self%varin(3) = "sfc_geopotential_height"
  !> need surface roughness
  self%varid(4) = "surface_roughness"
  !> need surface and atmospheric pressure for potential temperature
  self%varid(5) = "surface_pressure"
  self%varid(6) = "air_pressure"
  self%varid(7) = "air_temperature"
  self%varid(8) = "specific_humidity_mixing_ratio"
  self%varid(9) = "eastward_wind"
  self%varid(10) = "northward_wind"
  self%varid(11) = "landsea_mask"
  print *, 'varin'
  print *, self%varin
  print *, 'varout'
  print *, self%varout

end subroutine ufo_atmsfcinterp_setup

! ------------------------------------------------------------------------------
! TODO: add cleanup of your observation operator (optional)
subroutine ufo_atmsfcinterp_delete(self)
implicit none
class(ufo_atmsfcinterp), intent(inout) :: self

end subroutine ufo_atmsfcinterp_delete

! ------------------------------------------------------------------------------
subroutine ufo_atmsfcinterp_simobs(self, geovals, hofx, obss)
  use atmsfc_mod, only : calc_t2
  implicit none
  class(ufo_atmsfcinterp), intent(in)    :: self
  type(ufo_geovals),  intent(in)    :: geovals
  real(c_double),     intent(inout) :: hofx(:)
  type(c_ptr), value, intent(in)    :: obss
  type(ufo_geoval), pointer :: phi, hgt, profile
  integer :: nlocs, ivar, iobs
  real(kind_real), allocatable :: obselev(:), obshgt(:)
  character(len=MAXVARLEN) :: geovar


  ! for low altitudes, we can assume: phi = z
  ! get geopotential height profile
  call ufo_geovals_get_var(geovals, var_z, phi)
  ! get surface geopotential height
  !call ufo_geovals_get_var(geovals, var_sfc_z, hgt)

  ! get number of obs
  nlocs = obsspace_get_nlocs(obss)

  ! get station elevation from obs
  allocate(obselev(nlocs))
  call obsspace_get_db(obss, "MetaData", "station_elevation",obselev)

  ! get observation height (above sea level)
  allocate(obshgt(nlocs))
  call obsspace_get_db(obss, "MetaData", "height",obshgt)

  do ivar = 1, self%nvars
    ! Get the name of input variable in geovals
    geovar = self%varin(ivar)

    ! Get profile for this variable from geovals
    call ufo_geovals_get_var(geovals, geovar, profile)

    ! Interpolate from geovals to observational location into hofx
    ! Note: hofx holds all variables (varin) for location 1
    ! then all variables for location 2, and so on
    select case (trim(geovar))
      case ("air_temperature")
        ! this subroutine call will mimic that of what is performed in GSI to
        ! compute 2m temperature from other model variables
        call calc_t2
    end select
    !do iobs = 1, nlocs
    !  ! test by just grabbing the lowest model layer value
    !  hofx(ivar + (iobs-1)*self%nvars) = profile%vals(1,iobs) 
    !  if (trim(geovar) == 'surface_pressure') then
    !    hofx(ivar + (iobs-1)*self%nvars) = hofx(ivar +(iobs-1)*self%nvars)*10_kind_real
    !  end if
    !enddo
  enddo




end subroutine ufo_atmsfcinterp_simobs

! ------------------------------------------------------------------------------

end module ufo_atmsfcinterp_mod
