! (C) Copyright 2017-2018 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

!> Fortran module for atmsfcinterp observation operator

module ufo_atmsfcinterp_mod

  use oops_variables_mod
  use obs_variables_mod
  use ufo_vars_mod
  use kinds

  implicit none
  private

  !> Fortran derived type for the observation type
  type, public :: ufo_atmsfcinterp
  private
    type(obs_variables), public :: obsvars ! Variables to be simulated
    integer, allocatable, public :: obsvarindices(:) ! Indices of obsvars in the list of all
                                                     ! simulated variables in the ObsSpace.
                                                     ! allocated/deallocated in interface layer
    type(oops_variables), public :: geovars
    logical :: use_fact10
    real(kind_real) :: magl
  contains
    procedure :: setup  => atmsfcinterp_setup_
    procedure :: simobs => atmsfcinterp_simobs_
  end type ufo_atmsfcinterp

contains

! ------------------------------------------------------------------------------
subroutine atmsfcinterp_setup_(self, f_conf)
  use fckit_configuration_module, only: fckit_configuration
  implicit none
  class(ufo_atmsfcinterp), intent(inout) :: self
  type(fckit_configuration), intent(in)  :: f_conf
  integer :: nvars, ivar

  ! if true use the wind reduction factor
  call f_conf%get_or_die("use_fact10", self%use_fact10)

  !> add geopotential height
  call self%geovars%push_back(var_z)
  !> need skin temperature for near-surface interpolations
  call self%geovars%push_back(var_sfc_t)
  !> need surface geopotential height to get difference from phi
  call self%geovars%push_back(var_sfc_z)
  !> need surface roughness
  call self%geovars%push_back(var_sfc_rough)
  !> need surface and atmospheric pressure for potential temperature
  call self%geovars%push_back(var_ps)
  call self%geovars%push_back(var_prs)
  call self%geovars%push_back(var_prsi)
  call self%geovars%push_back(var_ts)
  call self%geovars%push_back(var_tv)
  call self%geovars%push_back(var_q)
  call self%geovars%push_back(var_u)
  call self%geovars%push_back(var_v)
  call self%geovars%push_back(var_sfc_lfrac)
  if (self%use_fact10)  call self%geovars%push_back(var_sfc_fact10)

end subroutine atmsfcinterp_setup_

! ------------------------------------------------------------------------------
subroutine atmsfcinterp_simobs_(self, geovals_in, obss, nvars, nlocs, hofx)
  use atmsfc_mod, only : calc_conv_vel_gsi, sfc_wind_fact_gsi, &
                         calc_psi_vars_gsi
  use thermo_utils_mod, only: calc_theta, gsi_tp_to_qs 
  use ufo_constants_mod, only: grav, rv, rd, rd_over_cp, von_karman
  use ufo_geovals_mod, only: ufo_geovals, ufo_geoval, ufo_geovals_get_var, ufo_geovals_copy, ufo_geovals_reorderzdir
  use ufo_utils_mod, only: cmp_strings
  use obsspace_mod
  use iso_c_binding
  use fckit_log_module, only: fckit_log
  implicit none
  class(ufo_atmsfcinterp), intent(in)        :: self
  integer, intent(in)                         :: nvars, nlocs
  type(ufo_geovals), intent(in)               :: geovals_in
  real(c_double),  intent(inout)              :: hofx(nvars, nlocs)
  type(c_ptr), value, intent(in)              :: obss
  type(ufo_geoval), pointer :: phi, hgt, tsfc, roughlen, psfc, prs, prsi, &
                               tsen, tv, q, u, v, landmask, &
                               profile, rad10
  type(ufo_geovals):: geovals
  integer :: ivar, iobs, iobsvar
  real(kind_real), allocatable :: obselev(:), obshgt(:)
  real(kind_real), parameter :: minroughlen = 1.0e-4_kind_real
  character(len=MAXVARLEN) :: geovar
  character(len=MAXVARLEN) :: var_zdir
  real(kind_real) :: thv1, thv2, th1, thg, thvg, rib, V2, agl, zbot
  real(kind_real) :: gzsoz0, gzzoz0
  real(kind_real) :: redfac, psim, psimz, psih, psihz 
  real(kind_real) :: ttmp1, ttmpg, eg, qg
  real(kind_real) :: z0, zq0, tvsfc 
  real(kind_real), parameter :: fv = rv/rd - 1.0_kind_real
  real(kind_real) :: psit, psitz, ust, psiq, psiqz
  real(kind_real), parameter :: zint0 = 0.01_kind_real ! default roughness over land
  real(kind_real), parameter :: ka = 2.4e-5_kind_real
  character(1024) :: debug_msg

  ! Quickly exit if nlocs is less than one, which happens on many CPUs and
  ! no observations sent to one of them.
  if (nlocs .lt. 1) return

  ! Notes:
  ! (1) Set desired vertical coordinate direction (top2bottom or bottom2top) based
  !     on vertical coodinate variable and reload geovals according to the set
  !     direction
  ! (2) This is done because this observation operator assumes pressure levels
  !     are from bottom to top (bottom2top) with model pressure(1) for surface and
  !     model pressure(nsig+1) for model top

  call ufo_geovals_copy(geovals_in, geovals)  ! dont want to change geovals_in
  var_zdir = var_prsi                         ! vertical coordinate variable
  call ufo_geovals_reorderzdir(geovals, var_zdir, "bottom2top")

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
  call ufo_geovals_get_var(geovals, var_prsi, prsi)
  call ufo_geovals_get_var(geovals, var_ts, tsen)
  call ufo_geovals_get_var(geovals, var_tv, tv)
  call ufo_geovals_get_var(geovals, var_q, q)
  call ufo_geovals_get_var(geovals, var_u, u)
  call ufo_geovals_get_var(geovals, var_v, v)
  call ufo_geovals_get_var(geovals, var_sfc_lfrac, landmask)
  if (self%use_fact10) call ufo_geovals_get_var(geovals, var_sfc_fact10, rad10)

  ! get station elevation from obs
  allocate(obselev(nlocs))
  call obsspace_get_db(obss, "MetaData", "stationElevation", obselev)

  ! get observation height (above sea level)
  allocate(obshgt(nlocs))
  call obsspace_get_db(obss, "MetaData", "height", obshgt)
  
  do iobs = 1, nlocs
    ! minimum roughness length
    z0 = roughlen%vals(1,iobs)
    if (z0 < minroughlen) z0 = minroughlen
    ! roughness length for over water
    zq0 = zint0
    if (landmask%vals(1,iobs) < 0.01) zq0 = z0

    ! get virtual temperature of the ground assuming saturation
    call gsi_tp_to_qs(tsfc%vals(1,iobs), psfc%vals(1,iobs), eg, qg)
    tvsfc = tsfc%vals(1,iobs) * (1.0_kind_real + fv * qg)

    ! get potential temperatures for calculating psi
    call calc_theta(tv%vals(1,iobs), prs%vals(1,iobs), thv1)
    call calc_theta(tv%vals(2,iobs), prs%vals(2,iobs), thv2)
    call calc_theta(tsen%vals(1,iobs), prs%vals(1,iobs), th1)
    call calc_theta(tsfc%vals(1,iobs), psfc%vals(1,iobs), thg)
    call calc_theta(tvsfc, psfc%vals(1,iobs), thvg)

    ! calculate convective velocity
    call calc_conv_vel_gsi(u%vals(1,iobs), v%vals(1,iobs), thvg, thv1, V2)

    ! although there could be first height below sea level, it causes floating-pt-except
    zbot = max(0.1,phi%vals(1,iobs))                ! bottom model level in meters
    agl = max(1.0, (obshgt(iobs)-obselev(iobs)))    ! obs height above ground in meters

    ! calculate bulk richardson number
    rib = (grav * zbot / th1) * (thv1 - thvg) / V2

    gzsoz0 = log(zbot/z0)
    gzzoz0 = log(agl/z0)

    ! calculate parameters regardless of variable
    call calc_psi_vars_gsi(rib, gzsoz0, gzzoz0, thv1, thv2, V2, th1,   &
                           thg, zbot, agl, psim, psih, psimz, psihz)

    do iobsvar = 1, size(self%obsvarindices)
      ! Get the index of the row of hofx to fill
      ivar = self%obsvarindices(iobsvar)

      ! Get the name of input variable in geovals
      geovar = self%geovars%variable(iobsvar)

      ! Get profile for this variable from geovals
      call ufo_geovals_get_var(geovals, geovar, profile)

      select case(trim(geovar))
        case("air_temperature", "virtual_temperature")
          psit = gzsoz0 - psih
          psitz = gzzoz0 - psihz
          if (cmp_strings(geovar, "air_temperature")) then
            ttmp1 = th1
            ttmpg = thg
          else
            ttmp1 = thv1
            ttmpg = thvg
          end if
          hofx(ivar,iobs) = (ttmpg + (ttmp1 - ttmpg)*psitz/psit)*(psfc%vals(1,iobs)/1.0e5_kind_real)**rd_over_cp
        case("eastward_wind", "northward_wind")
          if (self%use_fact10) then ! use provided fact10 from model
            hofx(ivar,iobs) = profile%vals(1,iobs) * rad10%vals(1,iobs)
          else ! compute wind reduction factor
            call sfc_wind_fact_gsi(u%vals(1,iobs), v%vals(1,iobs), tsen%vals(1,iobs), q%vals(1,iobs),&
                                   psfc%vals(1,iobs), prsi%vals(1,iobs), prsi%vals(2,iobs),&
                                   tsfc%vals(1,iobs), z0, landmask%vals(1,iobs), redfac)
            hofx(ivar,iobs) = profile%vals(1,iobs) * redfac
          end if
        case("specific_humidity")
          ust = max(0.01_kind_real, von_karman * sqrt(V2) / (gzsoz0 - psim))
          psiq = log(von_karman*ust*zbot/ka + zbot/zq0) - psih
          psiqz = log(von_karman*ust*agl/ka + agl/zq0) - psihz
          hofx(ivar,iobs) = qg + (q%vals(1,iobs) - qg)*psiqz/psiq
      end select
    end do
  end do

  deallocate(obshgt,obselev)

end subroutine atmsfcinterp_simobs_

! ------------------------------------------------------------------------------

end module ufo_atmsfcinterp_mod
