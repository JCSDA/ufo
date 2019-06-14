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
  use ufo_geovals_mod_c, only: ufo_geovals_registry
  use ufo_basis_mod, only: ufo_basis
  use ufo_vars_mod
  use obsspace_mod

  implicit none
  private

  integer, parameter :: max_string = 50
  !> Fortran derived type for the observation type
  type, public :: ufo_atmsfcinterp
  private
    integer :: nvars
    logical :: use_fact10
    real(kind_real) :: magl
    character(len=MAXVARLEN), public, allocatable :: varin(:)
    character(len=MAXVARLEN), public, allocatable :: varout(:)
  contains
    procedure :: setup  => atmsfcinterp_setup_
    procedure :: simobs => atmsfcinterp_simobs_
    final :: destructor
  end type ufo_atmsfcinterp

contains

! ------------------------------------------------------------------------------
subroutine atmsfcinterp_setup_(self, c_conf, vars)
  implicit none
  class(ufo_atmsfcinterp), intent(inout) :: self
  type(c_ptr),        intent(in)    :: c_conf
  character(len=MAXVARLEN), dimension(:), intent(inout) :: vars
  integer :: ii, ivar, nallvars, istart, fact10tmp
  logical :: fact10nml

  !> Size of variables
  self%nvars = size(vars)
  !> Allocate varout: variables in the observation vector
  allocate(self%varout(self%nvars))
  !> Read variable list and store in varout
  self%varout = vars
  ! check for if we need to look for wind reduction factor
  self%use_fact10 = .false. 
  fact10tmp = config_get_int(c_conf, "use_fact10", 0)
  if (fact10tmp /= 0) then
    self%use_fact10 = .true.
  end if

  !> Allocate varin: variables we need from the model
  if (self%use_fact10) then
    istart = 13
  else
    istart = 12
  end if
  nallvars = self%nvars + istart
  allocate(self%varin(nallvars))
  do ii = 1, self%nvars
    self%varin(ii+istart) = self%varout(ii)
  enddo

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
  self%varin(8) = var_tv
  self%varin(9) = var_q
  self%varin(10) = var_u 
  self%varin(11) = var_v 
  self%varin(12) = var_sfc_lfrac
  if (self%use_fact10)  self%varin(13) = var_sfc_fact10

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
                               tsen, tv, q, u, v, landmask, &
                               profile, rad10
  integer :: ivar, iobs
  real(kind_real), allocatable :: obselev(:), obshgt(:)
  real(kind_real) :: outvalue
  real(kind_real), parameter :: minroughlen = 0.0001_kind_real
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
  call ufo_geovals_get_var(geovals, var_tv, tv)
  call ufo_geovals_get_var(geovals, var_q, q)
  call ufo_geovals_get_var(geovals, var_u, u)
  call ufo_geovals_get_var(geovals, var_v, v)
  call ufo_geovals_get_var(geovals, var_sfc_lfrac, landmask)
  if (self%use_fact10) call ufo_geovals_get_var(geovals, var_sfc_fact10, rad10)

  ! get station elevation from obs
  allocate(obselev(nlocs))
  call obsspace_get_db(obss, "MetaData", "station_elevation", obselev)

  ! get observation height (above sea level)
  allocate(obshgt(nlocs))
  call obsspace_get_db(obss, "MetaData", "height", obshgt)
  
  do iobs = 1, nlocs
    ! minimum roughness length
    z0 = roughlen%vals(1,iobs)
    if (z0 < minroughlen) z0 = minroughlen

    do ivar = 1, self%nvars
      ! Get the name of input variable in geovals
      geovar = self%varout(ivar)
      ! Get profile for this variable from geovals
      call ufo_geovals_get_var(geovals, geovar, profile)

      select case(trim(geovar))
        case("air_temperature", "virtual_temperature")
          ! calling a modified version of the sfc_model routine from GSI
            call sfc_wtq_fwd_gsi(psfc%vals(1,iobs),tsfc%vals(1,iobs),prs%vals(1,iobs),&
                                 tsen%vals(1,iobs),q%vals(1,iobs),u%vals(1,iobs),&
                                 v%vals(1,iobs),prs%vals(2,iobs),tsen%vals(2,iobs),&
                                 q%vals(2,iobs),phi%vals(1,iobs),roughlen%vals(1,iobs),&
                                 landmask%vals(1,iobs),obshgt(iobs)-obselev(iobs),&
                                 hofx(ivar,iobs),geovar)
        case("eastward_wind", "northward_wind")
          if (self%use_fact10) then ! use provided fact10 from model
            hofx(ivar,iobs) = profile%vals(1,iobs) * rad10%vals(1,iobs)
          else ! compute wind reduction factor
            call sfc_wind_fact_gsi(z0, phi%vals(1,iobs), obshgt(iobs)-obselev(iobs), psim, psimz, redfac)
            hofx(ivar,iobs) = profile%vals(1,iobs) * redfac
          end if
      end select
    end do
  end do




end subroutine atmsfcinterp_simobs_

! ------------------------------------------------------------------------------

end module ufo_atmsfcinterp_mod
