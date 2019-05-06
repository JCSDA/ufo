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
  ! TODO: additional variables may need to be added but all of the 'standard'
  ! vars need the same input from the GeoVaLs (T,P,U,V,Q need all for each other)
  !> add geopotential height
  self%varin(1) = var_z
  !> need skin temperature for near-surface interpolations
  self%varin(2) = var_sfc_t
  !> need surface geopotential height to get difference from phi
  self%varin(3) = var_sfc_z 
  !> need surface roughness
  self%varin(4) = var_sfc_rough 
  !> need surface and atmospheric pressure for potential temperature
  self%varin(5) = var_ps
  self%varin(6) = var_prs
  self%varin(7) = var_ts
  self%varin(8) = var_q
  self%varin(9) = var_u 
  self%varin(10) = var_v 
  self%varin(11) = var_sfc_lfrac

end subroutine ufo_atmsfcinterp_setup

! ------------------------------------------------------------------------------
! TODO: add cleanup of your observation operator (optional)
subroutine ufo_atmsfcinterp_delete(self)
implicit none
class(ufo_atmsfcinterp), intent(inout) :: self

end subroutine ufo_atmsfcinterp_delete

! ------------------------------------------------------------------------------
subroutine ufo_atmsfcinterp_simobs(self, geovals, hofx, obss)
  use atmsfc_mod, only : sfc_wtq_fwd_gsi
  implicit none
  class(ufo_atmsfcinterp), intent(in)    :: self
  type(ufo_geovals),  intent(in)    :: geovals
  real(c_double),     intent(inout) :: hofx(:)
  type(c_ptr), value, intent(in)    :: obss
  type(ufo_geoval), pointer :: phi, hgt, tsfc, roughlen, psfc, prs, &
                               tsen, q, u, v, landmask, &
                               profile
  integer :: nlocs, ivar, iobs
  real(kind_real), allocatable :: obselev(:), obshgt(:)
  real(kind_real) :: outvalue
  character(len=MAXVARLEN) :: geovar


  ! to compute the value near the surface we are going to use
  ! similarity theory which requires a number of near surface parameters
  ! that need to be grabbed from GeoVaLs regardless of the observation type

  ! for low altitudes, we can assume: phi = z
  ! grabbing geopotential height profile and surface elevation,
  ! surface temperature, surface roughness length
  ! surface and profile air pressure
  ! profiles of tsen, q, u, v, and surface land fraction
  call ufo_geovals_get_var(geovals, var_z, phi)
  call ufo_geovals_get_var(geovals, var_sfc_z, hgt)
  call ufo_geovals_get_var(geovals, var_sfc_t, tsfc)
  call ufo_geovals_get_var(geovals, var_sfc_rough, roughlen)
  call ufo_geovals_get_var(geovals, var_ps, psfc)
  call ufo_geovals_get_var(geovals, var_prs, prs)
  call ufo_geovals_get_var(geovals, var_ts, tsen)
  call ufo_geovals_get_var(geovals, var_q, q)
  call ufo_geovals_get_var(geovals, var_u, u)
  call ufo_geovals_get_var(geovals, var_v, v)
  call ufo_geovals_get_var(geovals, var_sfc_lfrac, landmask)
  ! get number of obs
  nlocs = obsspace_get_nlocs(obss)

  ! get station elevation from obs
  allocate(obselev(nlocs))
  call obsspace_get_db(obss, "MetaData", "station_elevation", obselev)

  ! get observation height (above sea level)
  allocate(obshgt(nlocs))
  call obsspace_get_db(obss, "MetaData", "height", obshgt)

  do ivar = 1, self%nvars
    ! Get the name of input variable in geovals
    geovar = self%varin(ivar)

    ! Get profile for this variable from geovals
    call ufo_geovals_get_var(geovals, geovar, profile)

    ! calling a modified version of the sfc_model routine from GSI
    do iobs = 1, nlocs
      call sfc_wtq_fwd_gsi(psfc%vals(1,iobs),tsfc%vals(1,iobs),prs%vals(1,iobs),&
                           tsen%vals(1,iobs),q%vals(1,iobs),u%vals(1,iobs),&
                           v%vals(1,iobs),prs%vals(2,iobs),tsen%vals(2,iobs),&
                           q%vals(2,iobs),phi%vals(1,iobs),roughlen%vals(1,iobs),&
                           landmask%vals(1,iobs),2.0_kind_real,& ! force 2m agl for testing...
                           !landmask%vals(1,iobs),obshgt(iobs)-hgt%vals(1,iobs),&
                           profile%vals(1,iobs),profile%vals(2,iobs),outvalue)
    enddo
  enddo




end subroutine ufo_atmsfcinterp_simobs

! ------------------------------------------------------------------------------

end module ufo_atmsfcinterp_mod
