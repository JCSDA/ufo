! (C) Copyright 2017-2018 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

!> Fortran module to handle radiancecrtm observations

module ufo_crtm_passive_mod

 use CRTM_Module

 use fckit_configuration_module, only: fckit_configuration
 use, intrinsic :: iso_c_binding
 use kinds
 use missing_values_mod

 use obsspace_mod

 use ufo_geovals_mod, only: ufo_geovals, ufo_geoval, ufo_geovals_get_var
 use ufo_vars_mod
 use ufo_crtm_utils_mod

 use ufo_constants_mod, only: deg2rad

 implicit none
 private


 public ufo_crtm_passive_sim
 public ufo_crtm_passive_diag

contains


subroutine ufo_crtm_passive_sim(rts, Options, is_vis_or_uv, nvars, nlocs, n_Profiles, &
                                n_Channels, hofx)

implicit none
integer(c_size_t), intent(in) :: nvars
integer(c_size_t),        intent(in) :: nlocs
integer, intent(in) :: n_Profiles, n_Channels
real(c_double),        intent(inout) :: hofx(nvars, nlocs) !h(x) to return
type(CRTM_RTSolution_type), intent(in) :: rts(:,:)  ! n_channels, n_profiles
type(CRTM_Options_type),  intent(in) :: Options(:)
logical,                  intent(in) :: is_vis_or_uv  ! sensor is visible/UV: hofx is albedo

real(kind=kind_real), parameter :: PI = ACOS(-1.)
real(c_double) :: missing
integer        :: l, m

! Put simulated brightness temperature (or reflectance/albedo) into hofx
! ----------------------------------------------

! Set missing value
missing = missing_value(missing)

!Set to missing, then retrieve non-missing profiles
hofx = missing

#if defined(CRTM_VERSION) && (CRTM_VERSION >= 3)
! For visible or UV, only daytime profiles (solar irradiance > 1) get an albedo;
! a solar zenith angle < 85 deg check is not implemented.
if (is_vis_or_uv) then
   do m = 1, n_Profiles
      if (Options(m)%Skip_Profile) cycle
      if (rts(1,m)%Solar_irradiance > one) then
         do l = 1, n_Channels
            hofx(l,m) = rts(l,m)%Radiance*PI/rts(l,m)%Solar_irradiance  ! Albedo
         end do
      end if
   end do
end if
#endif

if (.not. is_vis_or_uv) then
   do m = 1, n_Profiles
      if (Options(m)%Skip_Profile) cycle
      do l = 1, n_Channels
         hofx(l,m) = rts(l,m)%Brightness_Temperature
      end do
   end do
end if

end subroutine ufo_crtm_passive_sim


subroutine ufo_crtm_passive_diag( &
  rts, rts_K, atm, atm_K, sfc_K, conf, isensor, Options, channels, geovals, &
  obss, nvars, nlocs, n_Profiles, n_Layers, xstr_diags, ystr_diags, ch_diags, &
  hofxdiags, err_stat)
use fckit_mpi_module,   only: fckit_mpi_comm
use fckit_log_module,   only: fckit_log

implicit none

type(CRTM_RTSolution_type), intent(in) :: rts(:,:), rts_K(:,:)
type(CRTM_Atmosphere_type), intent(in) :: atm(:)
type(CRTM_Atmosphere_type),intent(in) :: atm_K(:,:)
type(CRTM_Surface_type),   intent(in) :: sfc_K(:,:)
type(crtm_conf), intent(in) :: conf
integer, intent(in):: channels(:)
type(ufo_geovals),        intent(in) :: geovals      !Inputs from the model
integer(c_size_t),        intent(in) :: nvars, nlocs
integer, intent(in) ::  n_Profiles, n_Layers
type(ufo_geovals),     intent(inout) :: hofxdiags    !non-h(x) diagnostics
type(c_ptr), value,       intent(in) :: obss         !ObsSpace
type(CRTM_Options_type),   intent(in) :: Options(:)
character(len=MAXVARLEN), dimension(:), intent(in) :: &
                          ystr_diags, xstr_diags
integer, intent(in) :: ch_diags(:)
integer, intent(in) :: isensor
integer, intent(out) :: err_stat

! Local Variables
character(max_string) :: err_msg
integer :: jvar, jprofile, jlevel, jchannel, jspec
real(kind=kind_real), parameter :: PI = ACOS(-1.)
real(c_double) :: missing
real(kind_real) :: wfunc_max
real(kind_real) :: Wfunc(n_Layers)
real(kind_real) :: conv_albedo
! gmi_gpm channels 10-13 use the second-feedhorn geometry (MetaData/sensorZenithAngle1)
logical :: gmi_hf
! sensor zenith angle (and its gmi_gpm channel 10-13 version), read on first use
real(kind_real), allocatable :: zenith(:), zenith_hf(:)

! Set missing value
missing = missing_value(missing)

! Put simulated diagnostics into hofxdiags

err_stat = 0
! ----------------------------------------------
do jvar = 1, hofxdiags%nvar
   if (len(trim(hofxdiags%variables(jvar))) < 1) cycle

   if (ch_diags(jvar) > 0) then
      if (count(channels == ch_diags(jvar)) /= 1) then
         write(err_msg,"(a,i0,a,*(1x,i0))") "ufo_crtm_passive_diag: ch_diags(jvar) = ", ch_diags(jvar), &
                                            " must appear exactly once in the h(x) channels:", channels
         call abor1_ftn(err_msg)
      end if
   end if

   ! position of the diagnosed channel in channels (-1 if it is not there)
   jchannel = findloc(channels, ch_diags(jvar), 1)
   if (jchannel == 0) jchannel = -1

   if (allocated(hofxdiags%geovals(jvar)%vals)) then
     deallocate(hofxdiags%geovals(jvar)%vals)
   end if

   gmi_hf = (conf%SENSOR_ID(isensor) == "gmi_gpm" .and. ch_diags(jvar) > 9)
   !============================================
   ! Diagnostics used for QC and bias correction
   !============================================
   if (xstr_diags(jvar) == "") then
      ! forward h(x) diags
      select case(ystr_diags(jvar))
         ! variable: optical_thickness_of_atmosphere_layer by channel and location
         case (var_opt_depth)
            call alloc_diag(n_Layers)
            do jprofile = 1, n_Profiles
               if (Options(jprofile)%Skip_Profile) cycle
               hofxdiags%geovals(jvar)%vals(:,jprofile) = &
                  rts(jchannel,jprofile) % layer_optical_depth(1:n_Layers)
            end do

         ! variable: toa_outgoing_radiance_per_unit_wavenumber by channel and location
         ! [mW / (m^2 sr cm^-1)]
         case (var_radiance)
            call alloc_diag(1)
            do jprofile = 1, n_Profiles
               if (Options(jprofile)%Skip_Profile) cycle
               hofxdiags%geovals(jvar)%vals(1,jprofile) = rts(jchannel,jprofile) % Radiance
            end do

         ! variable: brightness_temperature_assuming_clear_sky by channel and location
         case (var_tb_clr)
            call alloc_diag(1)
            do jprofile = 1, n_Profiles
               if (Options(jprofile)%Skip_Profile) cycle
               ! Note: Using Tb_Clear requires CRTM_Atmosphere_IsFractional(cloud_coverage_flag)
               ! to be true. For CRTM v2.3.0, that happens when
               ! atm(jprofile)%Cloud_Fraction > MIN_COVERAGE_THRESHOLD (1e.-6)
               hofxdiags%geovals(jvar)%vals(1,jprofile) = rts(jchannel,jprofile) % Tb_Clear
            end do

         ! variable: brightness_temperature by channel and location
         case (var_tb)
            call alloc_diag(1)
            do jprofile = 1, n_Profiles
               if (Options(jprofile)%Skip_Profile) cycle
               hofxdiags%geovals(jvar)%vals(1,jprofile) = &
                  rts(jchannel,jprofile) % Brightness_Temperature
            end do

#if defined(CRTM_VERSION) && (CRTM_VERSION >= 3)
         ! variable: albedo by channel and location
         case (var_albedo)
            call alloc_diag(1)
            do jprofile = 1, n_Profiles
               if (Options(jprofile)%Skip_Profile) cycle
               hofxdiags%geovals(jvar)%vals(1,jprofile) = &
                  rts(jchannel,jprofile) % Radiance * PI &
                      / rts(jchannel,jprofile) % Solar_irradiance
            end do

         ! variable: albedo_assuming_clear_sky by channel and location
         case (var_albedo_clr)
            call alloc_diag(1)
            do jprofile = 1, n_Profiles
               if (Options(jprofile)%Skip_Profile) cycle
               hofxdiags%geovals(jvar)%vals(1,jprofile) = &
                  rts(jchannel,jprofile) % R_clear * PI &
                      / rts(jchannel,jprofile) % Solar_irradiance
            end do
#endif

         ! variable: transmittances_of_atmosphere_layer by channel and location
         case (var_lvl_transmit)
            call alloc_diag(n_Layers)
            do jprofile = 1, n_Profiles
               if (Options(jprofile)%Skip_Profile) cycle
               call layer_to_space_transmittance(jprofile, &
                                                 hofxdiags%geovals(jvar)%vals(:,jprofile))
            end do

         ! variable: weightingfunction_of_atmosphere_layer by channel and location
         case (var_lvl_weightfunc)
            call alloc_diag(n_Layers)
            do jprofile = 1, n_Profiles
               if (Options(jprofile)%Skip_Profile) cycle
               call weighting_function(jprofile, hofxdiags%geovals(jvar)%vals(:,jprofile))
            end do

         ! variable: pressure_level_at_peak_of_weightingfunction by channel and location
         case (var_pmaxlev_weightfunc)
            call alloc_diag(1)
            do jprofile = 1, n_Profiles
               if (Options(jprofile)%Skip_Profile) cycle
               call weighting_function(jprofile, Wfunc)
               ! get pressure level at the peak of the weighting function
               wfunc_max = -999.0
               do jlevel = n_Layers-1, 1, -1
                  if (Wfunc(jlevel) > wfunc_max) then
                     wfunc_max = Wfunc(jlevel)
                     hofxdiags%geovals(jvar)%vals(1,jprofile) = jlevel
                  end if
               end do
            end do

         ! variable: surface_emissivity by channel and location
         case (var_sfc_emiss)
            call alloc_diag(1)
            do jprofile = 1, n_Profiles
               if (Options(jprofile)%Skip_Profile) cycle
               hofxdiags%geovals(jvar)%vals(1, jprofile) = &
                  rts(jchannel, jprofile) % surface_emissivity
            end do

         case default
            ! unsupported forward diagnostic: fill with missing
            call alloc_diag(1)
      end select
   else if (ystr_diags(jvar) == var_tb) then
      ! var_tb jacobians
      select case (xstr_diags(jvar))
         ! variable: brightness_temperature_jacobian_air_temperature by channel and location
         case (var_ts)
            call alloc_diag(n_Layers)
            do jprofile = 1, n_Profiles
               if (Options(jprofile)%Skip_Profile) cycle
               hofxdiags%geovals(jvar)%vals(:,jprofile) = &
                  atm_K(jchannel,jprofile) % Temperature(1:n_Layers)
            end do

         ! variable: brightness_temperature_jacobian_<absorber> by channel and location, for
         ! water_vapor_mixing_ratio_wrt_dry_air, mole_fraction_of_carbon_dioxide_in_air and
         ! mole_fraction_of_ozone_in_air (nval==n_Layers --> requires MAXVARLEN=71)
         case (var_mixr, var_co2, var_oz)
            call alloc_diag(n_Layers)
            jspec = ufo_vars_getindex(conf%Absorbers, xstr_diags(jvar))
            do jprofile = 1, n_Profiles
               if (Options(jprofile)%Skip_Profile) cycle
               hofxdiags%geovals(jvar)%vals(:,jprofile) = &
                  atm_K(jchannel,jprofile) % Absorber(1:n_Layers,jspec)
            end do

         ! variable: brightness_temperature_jacobian_mass_content_of_<species>_in_atmosphere_layer
         ! by channel and location, for cloud liquid water, cloud ice, snow, rain, graupel and hail
         case (var_clw_wp, var_cli_wp, var_cls_wp, var_clr_wp, var_clg_wp, var_clh_wp)
            call alloc_diag(n_Layers)
            jspec = ufo_vars_getindex(conf%Clouds(:,1), xstr_diags(jvar))
            do jprofile = 1, n_Profiles
               if (Options(jprofile)%Skip_Profile) cycle
               hofxdiags%geovals(jvar)%vals(:,jprofile) = &
                  atm_K(jchannel,jprofile) % Cloud(jspec) % Water_Content(1:n_Layers)
            end do

         ! variable: brightness_temperature_jacobian_skin_temperature_at_surface by channel and
         ! location
         case (var_sfc_tskin)
            call alloc_diag(1)
            do jprofile = 1, n_Profiles
               if (Options(jprofile)%Skip_Profile) cycle
               hofxdiags%geovals(jvar)%vals(1,jprofile) = &
                  sfc_K(jchannel,jprofile) % water_temperature &
                + sfc_K(jchannel,jprofile) % land_temperature &
                + sfc_K(jchannel,jprofile) % ice_temperature &
                + sfc_K(jchannel,jprofile) % snow_temperature
            end do

         ! variable: brightness_temperature_jacobian_surface_emissivity by channel and location
         ! (nval=1)
         case (var_sfc_emiss)
            call alloc_diag(1)
            do jprofile = 1, n_Profiles
               if (Options(jprofile)%Skip_Profile) cycle
               hofxdiags%geovals(jvar)%vals(1,jprofile) = &
                  rts_K(jchannel,jprofile) % surface_emissivity
            end do

         case default
            err_msg = "ufo_crtm_passive_diag: ObsDiagnostic is unsupported, " // &
                      hofxdiags%variables(jvar)
            call fckit_log%error(err_msg)
            err_stat = 1
      end select
#if defined(CRTM_VERSION) && (CRTM_VERSION >= 3)
   else if (ystr_diags(jvar) == var_albedo) then
      ! var_albedo jacobians
      select case (xstr_diags(jvar))
         ! variable: albedo_jacobian_mass_content_of_<species>_in_atmosphere_layer by channel and
         ! location, for cloud liquid water, cloud ice, snow, rain, graupel and hail
         case (var_clw_wp, var_cli_wp, var_cls_wp, var_clr_wp, var_clg_wp, var_clh_wp)
            call alloc_diag(n_Layers)
            jspec = ufo_vars_getindex(conf%Clouds(:,1), xstr_diags(jvar))
            do jprofile = 1, n_Profiles
               if (Options(jprofile)%Skip_Profile) cycle
               ! Albedo conversion factor
               conv_albedo = PI/rts(jchannel,jprofile)%Solar_irradiance
               hofxdiags%geovals(jvar)%vals(:,jprofile) = conv_albedo *     &
                  atm_K(jchannel,jprofile) % Cloud(jspec) % Water_Content(1:n_Layers)
            end do

         ! variable: albedo_jacobian_surface_emissivity by channel and location
         case (var_sfc_emiss)
            call alloc_diag(1)
            do jprofile = 1, n_Profiles
               if (Options(jprofile)%Skip_Profile) cycle
               ! Albedo conversion factor
               conv_albedo = PI/rts(jchannel,jprofile)%Solar_irradiance
               hofxdiags%geovals(jvar)%vals(1,jprofile) = conv_albedo *     &
                  rts_K(jchannel,jprofile) % surface_emissivity
            end do

         case default
            err_msg = "ufo_crtm_passive_diag: ObsDiagnostic is unsupported, " // &
                      hofxdiags%variables(jvar)
            call fckit_log%error(err_msg)
            err_stat = 1
      end select
#endif
   else
      err_msg = "ufo_crtm_passive_diag: ObsDiagnostic is unsupported, " // &
                hofxdiags%variables(jvar)
      call fckit_log%error(err_msg)
      err_stat = 1
   end if
end do

contains

  !> Allocate hofxdiags variable jvar with nval levels, filled with missing.
  subroutine alloc_diag(nval)
    integer, intent(in) :: nval
    hofxdiags%geovals(jvar)%nval = nval
    allocate(hofxdiags%geovals(jvar)%vals(nval, n_Profiles))
    hofxdiags%geovals(jvar)%vals = missing
  end subroutine alloc_diag

  !> Secant of the sensor zenith angle for profile k (missing if the angle is missing),
  !> using the same angle as the H(x) geometry (ufo_crtm_load_geom_data): absolute
  !> value, at most CRTM's MAX_SENSOR_ZENITH_ANGLE, and for gmi_gpm channels 10-13
  !> sensorZenithAngle1 if present. Angles are read from the ObsSpace once per call.
  real(kind_real) function secant(k)
    integer, intent(in) :: k
    real(kind_real) :: angle
    if (.not. allocated(zenith)) then
      allocate(zenith(n_Profiles))
      call obsspace_get_db(obss, "MetaData", "sensorZenithAngle", zenith)
    end if
    angle = zenith(k)
    if (gmi_hf) then
      if (.not. allocated(zenith_hf)) then
        allocate(zenith_hf(n_Profiles))
        if (obsspace_has(obss, "MetaData", "sensorZenithAngle1")) then
          call obsspace_get_db(obss, "MetaData", "sensorZenithAngle1", zenith_hf)
        else
          zenith_hf = zenith
        end if
      end if
      angle = zenith_hf(k)
    end if
    if (angle == missing) then
      secant = missing
    else
      secant = one/cos(min(abs(angle), MAX_SENSOR_ZENITH_ANGLE)*deg2rad)
    end if
  end function secant

  !> Layer-to-space transmittance of channel jchannel for profile k.
  subroutine layer_to_space_transmittance(k, tao)
    integer,         intent(in)  :: k
    real(kind_real), intent(out) :: tao(n_Layers)
    integer :: l
    real(kind_real) :: total_od, secant_term
    secant_term = secant(k)
    if (secant_term == missing) then
      tao = missing
      return
    end if
    total_od = zero
    do l = 1, n_Layers
      total_od = total_od + rts(jchannel,k) % layer_optical_depth(l)
      tao(l) = exp(-min(limit_exp,total_od*secant_term))
    end do
  end subroutine layer_to_space_transmittance

  !> Weighting function |dTau/dln(p)| of channel jchannel for profile k; the bottom
  !> layer repeats the value of the layer above it.
  subroutine weighting_function(k, wf)
    integer,         intent(in)  :: k
    real(kind_real), intent(out) :: wf(n_Layers)
    real(kind_real) :: tao(n_Layers)
    integer :: l
    call layer_to_space_transmittance(k, tao)
    if (tao(1) == missing) then
      wf = missing
      return
    end if
    do l = n_Layers-1, 1, -1
      wf(l) = abs( (tao(l+1)-tao(l))/ &
                   (log(atm(k)%pressure(l+1))-log(atm(k)%pressure(l))) )
    end do
    wf(n_Layers) = wf(n_Layers-1)
  end subroutine weighting_function

end subroutine ufo_crtm_passive_diag

! ------------------------------------------------------------------------------

end module ufo_crtm_passive_mod
