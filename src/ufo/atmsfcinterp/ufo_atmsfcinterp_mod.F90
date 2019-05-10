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

  integer, parameter :: max_string = 50
  real(kind_real), parameter :: grav = 9.80665e+0_kind_real
  !> Fortran derived type for the observation type
  type, public :: ufo_atmsfcinterp
  private
    integer :: nvars
    real(kind_real) :: magl
    character(len=MAXVARLEN), public, allocatable :: varin(:)
    character(len=MAXVARLEN), public, allocatable :: varout(:)
  contains
    procedure :: setup  => atmsfcinterp_setup_
    procedure :: simobs => atmsfcinterp_simobs_
    !final :: destructor
  end type ufo_atmsfcinterp

contains

! ------------------------------------------------------------------------------
subroutine atmsfcinterp_setup_(self, c_conf)
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
  allocate(self%varin(11+self%nvars))
  do ii = 1, self%nvars
    self%varin(ii+11) = self%varout(ii)
  enddo

  !> add geopotential height
  self%varin(1) = var_z
  !self%varin(self%nvars+1) = var_z
  !> need skin temperature for near-surface interpolations
  !self%varin(self%nvars+2) = var_sfc_t
  self%varin(2) = var_sfc_t
  !> need surface geopotential height to get difference from phi
  !self%varin(self%nvars+3) = var_sfc_z 
  self%varin(3) = var_sfc_z 
  !> need surface roughness
  self%varin(4) = var_sfc_rough 
  !self%varin(self%nvars+4) = var_sfc_rough 
  !> need surface and atmospheric pressure for potential temperature
  !self%varin(self%nvars+5) = var_ps
  !self%varin(self%nvars+6) = var_prs
  !self%varin(self%nvars+7) = var_ts
  !self%varin(self%nvars+8) = var_q
  !self%varin(self%nvars+9) = var_u 
  !self%varin(self%nvars+10) = var_v 
  !self%varin(self%nvars+11) = var_sfc_lfrac
  self%varin(5) = var_ps
  self%varin(6) = var_prs
  self%varin(7) = var_ts
  self%varin(8) = var_q
  self%varin(9) = var_u 
  self%varin(10) = var_v 
  self%varin(11) = var_sfc_lfrac
  print *, self%varin

end subroutine atmsfcinterp_setup_

! ------------------------------------------------------------------------------

subroutine atmsfcinterp_simobs_(self, geovals, obss, nvars, nlocs, hofx)
  use atmsfc_mod, only : sfc_wtq_fwd_gsi
  implicit none
  class(ufo_atmsfcinterp), intent(in)        :: self
  integer, intent(in)                         :: nvars, nlocs
  type(ufo_geovals), intent(in)               :: geovals
  real(c_double),  intent(inout)              :: hofx(nvars, nlocs)
  type(c_ptr), value, intent(in)              :: obss
  type(ufo_geoval), pointer :: phi, hgt, tsfc, roughlen, psfc, prs, &
                               tsen, q, u, v, landmask, &
                               profile
  integer :: ivar, iobs
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

  ! get station elevation from obs
  allocate(obselev(nlocs))
  call obsspace_get_db(obss, "MetaData", "station_elevation", obselev)

  ! get observation height (above sea level)
  allocate(obshgt(nlocs))
  call obsspace_get_db(obss, "MetaData", "height", obshgt)

  do ivar = 1, self%nvars
    ! Get the name of input variable in geovals
    geovar = self%varout(ivar)

    ! Get profile for this variable from geovals
    !call ufo_geovals_get_var(geovals, geovar, profile)

    ! calling a modified version of the sfc_model routine from GSI
    do iobs = 1, nlocs
      call sfc_wtq_fwd_gsi(psfc%vals(1,iobs),tsfc%vals(1,iobs),prs%vals(1,iobs),&
                           tsen%vals(1,iobs),q%vals(1,iobs),u%vals(1,iobs),&
                           v%vals(1,iobs),prs%vals(2,iobs),tsen%vals(2,iobs),&
                           q%vals(2,iobs),phi%vals(1,iobs),roughlen%vals(1,iobs),&
                           landmask%vals(1,iobs),2._kind_real,&
                           !landmask%vals(1,iobs),obshgt(iobs)-hgt%vals(1,iobs),&
                           hofx(ivar,iobs),geovar)
    enddo
  enddo




end subroutine atmsfcinterp_simobs_

! ------------------------------------------------------------------------------

end module ufo_atmsfcinterp_mod
