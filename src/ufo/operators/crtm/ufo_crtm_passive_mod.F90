! (C) Copyright 2017-2018 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

!> Fortran module to handle radiancecrtm observations

module ufo_crtm_passive_mod

 use crtm_module

 use fckit_configuration_module, only: fckit_configuration
 use iso_c_binding
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


subroutine ufo_crtm_passive_sim(rts, Options, nvars, nlocs, n_Profiles, n_Channels, hofx)


 USE crtm_SpcCoeff, ONLY: SC, &
                          SpcCoeff_IsMicrowaveSensor , & 
                          SpcCoeff_IsInfraredSensor  , &
                          SpcCoeff_IsVisibleSensor   , &
                          SpcCoeff_IsUltravioletSensor

implicit none
integer(c_size_t), intent(in) :: nvars
integer(c_size_t),        intent(in) :: nlocs
integer, intent(in) :: n_Profiles, n_Channels
real(c_double),        intent(inout) :: hofx(nvars, nlocs) !h(x) to return
type(CRTM_RTSolution_type), intent(in) :: rts(:,:)  ! n_channels, n_profiles
type(CRTM_Options_type),  intent(in) :: Options(:)

real(kind=kind_real), parameter :: PI = ACOS(-1.)
real(kind=kind_real), parameter :: cos85 = COS(85.)
real(c_double) :: missing
integer        :: l, m
logical        :: is_vis_or_uv = .false.

! Put simulated brightness temperature (or reflectance/albedo) into hofx
! ----------------------------------------------

! Set missing value
missing = missing_value(missing)

!Set to missing, then retrieve non-missing profiles
hofx = missing

IF ( ANY(SpcCoeff_IsVisibleSensor(SC)) .or. ANY(SpcCoeff_IsUltravioletSensor(SC)) ) then
   is_vis_or_uv = .true.
else
   is_vis_or_uv = .false.
end if

#if defined(CRTM_VERSION) && (CRTM_VERSION >= 3)
! For visible or UV, ensure that it is daytime and solar zenith angle is less than 85 deg (not yet).
if (is_vis_or_uv) then
   do m = 1, n_Profiles
      if (.not.Options(m)%Skip_Profile) then
         if (rts(1,m)%Solar_irradiance .gt. 1.0) then    ! .and. rts(1,m)%COS_SUN .gt. cos85) then
            do l = 1, n_Channels
               hofx(l,m) = rts(l,m)%Radiance*PI/rts(l,m)%Solar_irradiance  ! Albedo
            end do
         end if
      end if
   end do
end if
#endif

if (.not. is_vis_or_uv) then
   do m = 1, n_Profiles
      if (.not.Options(m)%Skip_Profile) then
         do l = 1, n_Channels
            hofx(l,m) = rts(l,m)%Brightness_Temperature
         end do
      end if
   end do
end if

end subroutine ufo_crtm_passive_sim


subroutine ufo_crtm_passive_diag(rts, rts_K, atm, atm_K, sfc_K, conf, n_Sensor, Options, channels, geovals, obss, nvars, nlocs, n_Profiles, n_Layers, xstr_diags, ystr_diags, ch_diags, hofxdiags, err_stat)
use fckit_mpi_module,   only: fckit_mpi_comm
use ufo_utils_mod,      only: cmp_strings

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
integer, intent(in) :: n_Sensor
integer, intent(out) :: err_stat

! Local Variables
!character(*), parameter :: PROGRAM_NAME = 'ufo_crtm_passive_diag'
character(max_string) :: err_msg
!integer        :: err_stat, alloc_stat
!integer        :: l, m, n
integer :: jvar, jprofile, jlevel, jchannel, ichannel, jspec
real(kind=kind_real), parameter :: PI = ACOS(-1.)
real(c_double) :: missing
real(kind_real) :: total_od, secant_term, wfunc_max
real(kind_real), allocatable :: TmpVar(:)
real(kind_real), allocatable :: Tao(:)
real(kind_real), allocatable :: Wfunc(:)
real(kind_real) :: conv_albedo
character(len=1) :: angle_hf

! Set missing value
missing = missing_value(missing)

! Put simulated diagnostics into hofxdiags

err_stat = 0
! ----------------------------------------------
do jvar = 1, hofxdiags%nvar
   if (len(trim(hofxdiags%variables(jvar))) < 1) cycle

   if (ch_diags(jvar) > 0) then
      if (size(pack(channels,channels==ch_diags(jvar))) /= 1) then
         write(err_msg,*) 'ufo_radiancecrtm_simobs: mismatch between// &
                           & h(x) channels(', channels,') and// &
                           & ch_diags(jvar) = ', ch_diags(jvar)
         call abor1_ftn(err_msg)
      end if
   end if

   jchannel = -1
   do ichannel = 1, size(channels)
      if (ch_diags(jvar) == channels(ichannel)) then
         jchannel = ichannel
         exit
      end if
   end do

   if (allocated(hofxdiags%geovals(jvar)%vals)) &
      deallocate(hofxdiags%geovals(jvar)%vals)

   angle_hf=achar(0)
   if (cmp_strings(conf%SENSOR_ID(n_Sensor),'gmi_gpm')) then
      if (ch_diags(jvar) > 9) then
         angle_hf="1"
      endif
   endif
   !============================================
   ! Diagnostics used for QC and bias correction
   !============================================
   if (cmp_strings(xstr_diags(jvar), "")) then
      ! forward h(x) diags
      select case(ystr_diags(jvar))
         ! variable: optical_thickness_of_atmosphere_layer_CH
         case (var_opt_depth)
            hofxdiags%geovals(jvar)%nval = n_Layers
            allocate(hofxdiags%geovals(jvar)%vals(hofxdiags%geovals(jvar)%nval,n_Profiles))
            hofxdiags%geovals(jvar)%vals = missing
            do jprofile = 1, n_Profiles
               if (.not.Options(jprofile)%Skip_Profile) then
                  do jlevel = 1, hofxdiags%geovals(jvar)%nval
                     hofxdiags%geovals(jvar)%vals(jlevel,jprofile) = &
                        rts(jchannel,jprofile) % layer_optical_depth(jlevel)
                  end do
               end if
            end do

         ! variable: toa_outgoing_radiance_per_unit_wavenumber_CH [mW / (m^2 sr cm^-1)] (nval=1)
         case (var_radiance)
            hofxdiags%geovals(jvar)%nval = 1
            allocate(hofxdiags%geovals(jvar)%vals(hofxdiags%geovals(jvar)%nval,n_Profiles))
            hofxdiags%geovals(jvar)%vals = missing
            do jprofile = 1, n_Profiles
               if (.not.Options(jprofile)%Skip_Profile) then
                  hofxdiags%geovals(jvar)%vals(1,jprofile) = &
                     rts(jchannel,jprofile) % Radiance
               end if
            end do

         ! variable: brightness_temperature_assuming_clear_sky_CH
         case (var_tb_clr)
            hofxdiags%geovals(jvar)%nval = 1
            allocate(hofxdiags%geovals(jvar)%vals(hofxdiags%geovals(jvar)%nval,n_Profiles))
            hofxdiags%geovals(jvar)%vals = missing
            do jprofile = 1, n_Profiles
               if (.not.Options(jprofile)%Skip_Profile) then
                  ! Note: Using Tb_Clear requires CRTM_Atmosphere_IsFractional(cloud_coverage_flag)
                  ! to be true. For CRTM v2.3.0, that happens when
                  ! atm(jprofile)%Cloud_Fraction > MIN_COVERAGE_THRESHOLD (1e.-6)
                  hofxdiags%geovals(jvar)%vals(1,jprofile) = &
                     rts(jchannel,jprofile) % Tb_Clear
               end if
            end do

         ! variable: brightness_temperature_CH
         case (var_tb)
            hofxdiags%geovals(jvar)%nval = 1
            allocate(hofxdiags%geovals(jvar)%vals(hofxdiags%geovals(jvar)%nval,n_Profiles))
            hofxdiags%geovals(jvar)%vals = missing
            do jprofile = 1, n_Profiles
               if (.not.Options(jprofile)%Skip_Profile) then
                  hofxdiags%geovals(jvar)%vals(1,jprofile) = &
                     rts(jchannel,jprofile) % Brightness_Temperature
               end if
            end do

#if defined(CRTM_VERSION) && (CRTM_VERSION >= 3)
         ! variable: albedo_CH
         case (var_albedo)
            hofxdiags%geovals(jvar)%nval = 1
            allocate(hofxdiags%geovals(jvar)%vals(hofxdiags%geovals(jvar)%nval,n_Profiles))
            hofxdiags%geovals(jvar)%vals = missing
            do jprofile = 1, n_Profiles
               if (.not.Options(jprofile)%Skip_Profile) then
                  hofxdiags%geovals(jvar)%vals(1,jprofile) = &
                     rts(jchannel,jprofile) % Radiance * PI / rts(jchannel,jprofile) % Solar_irradiance 
               end if
            end do

         ! variable: albedo_assuming_clear_sky_CH
         case (var_albedo_clr)
            hofxdiags%geovals(jvar)%nval = 1
            allocate(hofxdiags%geovals(jvar)%vals(hofxdiags%geovals(jvar)%nval,n_Profiles))
            hofxdiags%geovals(jvar)%vals = missing
            do jprofile = 1, n_Profiles
               if (.not.Options(jprofile)%Skip_Profile) then
                  hofxdiags%geovals(jvar)%vals(1,jprofile) = &
                     rts(jchannel,jprofile) % R_clear * PI / rts(jchannel,jprofile) % Solar_irradiance 
               end if
            end do
#endif

         ! variable: transmittances_of_atmosphere_layer_CH
         case (var_lvl_transmit)
            hofxdiags%geovals(jvar)%nval = n_Layers
            allocate(hofxdiags%geovals(jvar)%vals(hofxdiags%geovals(jvar)%nval,n_Profiles))
            hofxdiags%geovals(jvar)%vals = missing
            allocate(TmpVar(n_Profiles))
            !call obsspace_get_db(obss, "MetaData", "sensor_zenith_angle"//angle_hf, TmpVar)
            call obsspace_get_db(obss, "MetaData", "sensorZenithAngle"//angle_hf, TmpVar)
            do jprofile = 1, n_Profiles
               if (.not.Options(jprofile)%Skip_Profile) then
                  secant_term = one/cos(TmpVar(jprofile)*deg2rad)
                  total_od = 0.0
                  do jlevel = 1, n_Layers
                     total_od   = total_od + rts(jchannel,jprofile) % layer_optical_depth(jlevel)
                     hofxdiags%geovals(jvar)%vals(jlevel,jprofile) = &
                        exp(-min(limit_exp,total_od*secant_term))
                  end do
               end if
            end do
            deallocate(TmpVar)

         ! variable: weightingfunction_of_atmosphere_layer_CH
         case (var_lvl_weightfunc)
            hofxdiags%geovals(jvar)%nval = n_Layers
            allocate(hofxdiags%geovals(jvar)%vals(hofxdiags%geovals(jvar)%nval,n_Profiles))
            hofxdiags%geovals(jvar)%vals = missing
            allocate(TmpVar(n_Profiles))
            allocate(Tao(n_Layers))
            call obsspace_get_db(obss, "MetaData", "sensorZenithAngle"//angle_hf, TmpVar)
            do jprofile = 1, n_Profiles
               if (.not.Options(jprofile)%Skip_Profile) then
                  ! get layer-to-space transmittance
                  secant_term = one/cos(TmpVar(jprofile)*deg2rad)
                  total_od = 0.0
                  do jlevel = 1, n_Layers
                     total_od = total_od + rts(jchannel,jprofile) % layer_optical_depth(jlevel)
                     Tao(jlevel) = exp(-min(limit_exp,total_od*secant_term))
                  end do
                  ! get weighting function
                  do jlevel = n_Layers-1, 1, -1
                     hofxdiags%geovals(jvar)%vals(jlevel,jprofile) = &
                        abs( (Tao(jlevel+1)-Tao(jlevel))/ &
                             (log(atm(jprofile)%pressure(jlevel+1))- &
                              log(atm(jprofile)%pressure(jlevel))) )
                  end do
                  hofxdiags%geovals(jvar)%vals(n_Layers,jprofile) = &
                  hofxdiags%geovals(jvar)%vals(n_Layers-1,jprofile)
               end if
            end do
            deallocate(TmpVar)
            deallocate(Tao)

         ! variable: pressure_level_at_peak_of_weightingfunction_CH
         case (var_pmaxlev_weightfunc)
            hofxdiags%geovals(jvar)%nval = 1
            allocate(hofxdiags%geovals(jvar)%vals(hofxdiags%geovals(jvar)%nval,n_Profiles))
            hofxdiags%geovals(jvar)%vals = missing
            allocate(TmpVar(n_Profiles))
            allocate(Tao(n_Layers))
            allocate(Wfunc(n_Layers))
            call obsspace_get_db(obss, "MetaData", "sensorZenithAngle"//angle_hf, TmpVar)
            do jprofile = 1, n_Profiles
               if (.not.Options(jprofile)%Skip_Profile) then
                 ! get layer-to-space transmittance
                  secant_term = one/cos(TmpVar(jprofile)*deg2rad)
                  total_od = 0.0
                  do jlevel = 1, n_Layers
                     total_od = total_od + rts(jchannel,jprofile) % layer_optical_depth(jlevel)
                     Tao(jlevel) = exp(-min(limit_exp,total_od*secant_term))
                  end do
                  ! get weighting function
                  do jlevel = n_Layers-1, 1, -1
                     Wfunc(jlevel) = &
                        abs( (Tao(jlevel+1)-Tao(jlevel))/ &
                             (log(atm(jprofile)%pressure(jlevel+1))- &
                              log(atm(jprofile)%pressure(jlevel))) )
                  end do
                  Wfunc(n_Layers) = Wfunc(n_Layers-1)
                  ! get pressure level at the peak of the weighting function
                  wfunc_max = -999.0
                  do jlevel = n_Layers-1, 1, -1
                     if (Wfunc(jlevel) > wfunc_max) then
                        wfunc_max = Wfunc(jlevel)
                        hofxdiags%geovals(jvar)%vals(1,jprofile) = jlevel
                     endif
                  enddo
               end if
            end do
            deallocate(TmpVar)
            deallocate(Tao)
            deallocate(Wfunc)

         case default
            write(err_msg,*) 'ufo_radiancecrtm_simobs: //&
                              & ObsDiagnostic is unsupported, ', &
                              & hofxdiags%variables(jvar)
            ! call abor1_ftn(err_msg)
            hofxdiags%geovals(jvar)%nval = 1
            allocate(hofxdiags%geovals(jvar)%vals(hofxdiags%geovals(jvar)%nval,n_Profiles))
            hofxdiags%geovals(jvar)%vals = missing
      end select
   else if (ystr_diags(jvar) == var_tb) then
      ! var_tb jacobians
      select case (xstr_diags(jvar))
         ! variable: brightness_temperature_jacobian_air_temperature_CH
         case (var_ts)
            hofxdiags%geovals(jvar)%nval = n_Layers
            allocate(hofxdiags%geovals(jvar)%vals(hofxdiags%geovals(jvar)%nval,n_Profiles))
            hofxdiags%geovals(jvar)%vals = missing
            do jprofile = 1, n_Profiles
               if (.not.Options(jprofile)%Skip_Profile) then
                  do jlevel = 1, hofxdiags%geovals(jvar)%nval
                     hofxdiags%geovals(jvar)%vals(jlevel,jprofile) = &
                        atm_K(jchannel,jprofile) % Temperature(jlevel)
                  end do
               end if
            end do
         ! variable: brightness_temperature_jacobian_humidity_mixing_ratio_CH (nval==n_Layers) --> requires MAXVARLEN=58
         case (var_mixr)
            hofxdiags%geovals(jvar)%nval = n_Layers
            allocate(hofxdiags%geovals(jvar)%vals(hofxdiags%geovals(jvar)%nval,n_Profiles))
            hofxdiags%geovals(jvar)%vals = missing
            jspec = ufo_vars_getindex(conf%Absorbers, var_mixr)
            do jprofile = 1, n_Profiles
               if (.not.Options(jprofile)%Skip_Profile) then
                  do jlevel = 1, hofxdiags%geovals(jvar)%nval
                     hofxdiags%geovals(jvar)%vals(jlevel,jprofile) = &
                        atm_K(jchannel,jprofile) % Absorber(jlevel,jspec)
                  end do
               end if
            end do

         ! variable: brightness_temperature_jacobian_mass_content_of_cloud_liquid_water_in_atmosphere_layer_CH
         case (var_clw_wp)
            hofxdiags%geovals(jvar)%nval = n_Layers
            allocate(hofxdiags%geovals(jvar)%vals(hofxdiags%geovals(jvar)%nval,n_Profiles))
            hofxdiags%geovals(jvar)%vals = missing
            jspec = ufo_vars_getindex(conf%Clouds(:,1), var_clw_wp)
            do jprofile = 1, n_Profiles
               if (.not.Options(jprofile)%Skip_Profile) then
                  do jlevel = 1, hofxdiags%geovals(jvar)%nval
                     hofxdiags%geovals(jvar)%vals(jlevel,jprofile) = &
                        atm_K(jchannel,jprofile) % Cloud(jspec) % Water_Content(jlevel)
                  end do
               end if
            end do

         ! variable: brightness_temperature_jacobian_mass_content_of_cloud_ice_in_atmosphere_layer_CH
         case (var_cli_wp)
            hofxdiags%geovals(jvar)%nval = n_Layers
            allocate(hofxdiags%geovals(jvar)%vals(hofxdiags%geovals(jvar)%nval,n_Profiles))
            hofxdiags%geovals(jvar)%vals = missing
            jspec = ufo_vars_getindex(conf%Clouds(:,1), var_cli_wp)
            do jprofile = 1, n_Profiles
               if (.not.Options(jprofile)%Skip_Profile) then
                  do jlevel = 1, hofxdiags%geovals(jvar)%nval
                     hofxdiags%geovals(jvar)%vals(jlevel,jprofile) = &
                        atm_K(jchannel,jprofile) % Cloud(jspec) % Water_Content(jlevel)
                  end do
               end if
            end do

         ! variable: brightness_temperature_jacobian_mass_content_of_snow_in_atmosphere_layer_CH
         case (var_cls_wp)
            hofxdiags%geovals(jvar)%nval = n_Layers
            allocate(hofxdiags%geovals(jvar)%vals(hofxdiags%geovals(jvar)%nval,n_Profiles))
            hofxdiags%geovals(jvar)%vals = missing
            jspec = ufo_vars_getindex(conf%Clouds(:,1), var_cls_wp)
            do jprofile = 1, n_Profiles
               if (.not.Options(jprofile)%Skip_Profile) then
                  do jlevel = 1, hofxdiags%geovals(jvar)%nval
                     hofxdiags%geovals(jvar)%vals(jlevel,jprofile) = &
                        atm_K(jchannel,jprofile) % Cloud(jspec) % Water_Content(jlevel)
                  end do
               end if
            end do

         ! variable: brightness_temperature_jacobian_mass_content_of_rain_in_atmosphere_layer_CH
         case (var_clr_wp)
            hofxdiags%geovals(jvar)%nval = n_Layers
            allocate(hofxdiags%geovals(jvar)%vals(hofxdiags%geovals(jvar)%nval,n_Profiles))
            hofxdiags%geovals(jvar)%vals = missing
            jspec = ufo_vars_getindex(conf%Clouds(:,1), var_clr_wp)
            do jprofile = 1, n_Profiles
               if (.not.Options(jprofile)%Skip_Profile) then
                  do jlevel = 1, hofxdiags%geovals(jvar)%nval
                     hofxdiags%geovals(jvar)%vals(jlevel,jprofile) = &
                        atm_K(jchannel,jprofile) % Cloud(jspec) % Water_Content(jlevel)
                  end do
               end if
            end do

         ! variable: brightness_temperature_jacobian_mass_content_of_graupel_in_atmosphere_layer_CH
         case (var_clg_wp)
            hofxdiags%geovals(jvar)%nval = n_Layers
            allocate(hofxdiags%geovals(jvar)%vals(hofxdiags%geovals(jvar)%nval,n_Profiles))
            hofxdiags%geovals(jvar)%vals = missing
            jspec = ufo_vars_getindex(conf%Clouds(:,1), var_clg_wp)
            do jprofile = 1, n_Profiles
               if (.not.Options(jprofile)%Skip_Profile) then
                  do jlevel = 1, hofxdiags%geovals(jvar)%nval
                     hofxdiags%geovals(jvar)%vals(jlevel,jprofile) = &
                        atm_K(jchannel,jprofile) % Cloud(jspec) % Water_Content(jlevel)
                  end do
               end if
            end do

         ! variable: brightness_temperature_jacobian_mass_content_of_hail_in_atmosphere_layer_CH
         case (var_clh_wp)
            hofxdiags%geovals(jvar)%nval = n_Layers
            allocate(hofxdiags%geovals(jvar)%vals(hofxdiags%geovals(jvar)%nval,n_Profiles))
            hofxdiags%geovals(jvar)%vals = missing
            jspec = ufo_vars_getindex(conf%Clouds(:,1), var_clh_wp)
            do jprofile = 1, n_Profiles
               if (.not.Options(jprofile)%Skip_Profile) then
                  do jlevel = 1, hofxdiags%geovals(jvar)%nval
                     hofxdiags%geovals(jvar)%vals(jlevel,jprofile) = &
                        atm_K(jchannel,jprofile) % Cloud(jspec) % Water_Content(jlevel)
                  end do
               end if
            end do

         ! variable: brightness_temperature_jacobian_surface_temperature_CH (nval=1)
         case (var_sfc_t)
            hofxdiags%geovals(jvar)%nval = 1
            allocate(hofxdiags%geovals(jvar)%vals(hofxdiags%geovals(jvar)%nval,n_Profiles))
            hofxdiags%geovals(jvar)%vals = missing
            do jprofile = 1, n_Profiles
               if (.not.Options(jprofile)%Skip_Profile) then
                  hofxdiags%geovals(jvar)%vals(1,jprofile) = &
                     sfc_K(jchannel,jprofile) % water_temperature &
                   + sfc_K(jchannel,jprofile) % land_temperature &
                   + sfc_K(jchannel,jprofile) % ice_temperature &
                   + sfc_K(jchannel,jprofile) % snow_temperature
               end if
            end do

         ! variable: brightness_temperature_jacobian_surface_emissivity_CH (nval=1)
         case (var_sfc_emiss)
            hofxdiags%geovals(jvar)%nval = 1
            allocate(hofxdiags%geovals(jvar)%vals(hofxdiags%geovals(jvar)%nval,n_Profiles))
            hofxdiags%geovals(jvar)%vals = missing
            do jprofile = 1, n_Profiles
               if (.not.Options(jprofile)%Skip_Profile) then
                  hofxdiags%geovals(jvar)%vals(1,jprofile) = &
                     rts_K(jchannel,jprofile) % surface_emissivity
               end if
            end do

         case default
            write(err_msg,*) 'ufo_radiancecrtm_simobs: //&
                              & ObsDiagnostic is unsupported, ', &
                              & hofxdiags%variables(jvar)
            !call abor1_ftn(err_msg)
            err_stat = 1
      end select
#if defined(CRTM_VERSION) && (CRTM_VERSION >= 3)
   else if (ystr_diags(jvar) == var_albedo) then
      ! var_albedo jacobians
      select case (xstr_diags(jvar))
         ! variable: albedo_jacobian_mass_content_of_cloud_liquid_water_in_atmosphere_layer_CH
         case (var_clw_wp)
            hofxdiags%geovals(jvar)%nval = n_Layers
            allocate(hofxdiags%geovals(jvar)%vals(hofxdiags%geovals(jvar)%nval,n_Profiles))
            hofxdiags%geovals(jvar)%vals = missing
            jspec = ufo_vars_getindex(conf%Clouds(:,1), var_clw_wp)
            do jprofile = 1, n_Profiles
               if (.not.Options(jprofile)%Skip_Profile) then
                  do jlevel = 1, hofxdiags%geovals(jvar)%nval
                     conv_albedo = PI/rts(jchannel,jprofile)%Solar_irradiance  ! Albedo conversion factor
                     hofxdiags%geovals(jvar)%vals(jlevel,jprofile) = conv_albedo *     &
                        atm_K(jchannel,jprofile) % Cloud(jspec) % Water_Content(jlevel)
                  end do
               end if
            end do

         ! variable: albedo_jacobian_mass_content_of_cloud_ice_in_atmosphere_layer_CH
         case (var_cli_wp)
            hofxdiags%geovals(jvar)%nval = n_Layers
            allocate(hofxdiags%geovals(jvar)%vals(hofxdiags%geovals(jvar)%nval,n_Profiles))
            hofxdiags%geovals(jvar)%vals = missing
            jspec = ufo_vars_getindex(conf%Clouds(:,1), var_cli_wp)
            do jprofile = 1, n_Profiles
               if (.not.Options(jprofile)%Skip_Profile) then
                  do jlevel = 1, hofxdiags%geovals(jvar)%nval
                     conv_albedo = PI/rts(jchannel,jprofile)%Solar_irradiance  ! Albedo conversion factor
                     hofxdiags%geovals(jvar)%vals(jlevel,jprofile) = conv_albedo *     &
                        atm_K(jchannel,jprofile) % Cloud(jspec) % Water_Content(jlevel)
                  end do
               end if
            end do

         ! variable: albedo_jacobian_mass_content_of_snow_in_atmosphere_layer_CH
         case (var_cls_wp)
            hofxdiags%geovals(jvar)%nval = n_Layers
            allocate(hofxdiags%geovals(jvar)%vals(hofxdiags%geovals(jvar)%nval,n_Profiles))
            hofxdiags%geovals(jvar)%vals = missing
            jspec = ufo_vars_getindex(conf%Clouds(:,1), var_cls_wp)
            do jprofile = 1, n_Profiles
               if (.not.Options(jprofile)%Skip_Profile) then
                  do jlevel = 1, hofxdiags%geovals(jvar)%nval
                     conv_albedo = PI/rts(jchannel,jprofile)%Solar_irradiance  ! Albedo conversion factor
                     hofxdiags%geovals(jvar)%vals(jlevel,jprofile) = conv_albedo *     &
                        atm_K(jchannel,jprofile) % Cloud(jspec) % Water_Content(jlevel)
                  end do
               end if
            end do

         ! variable: albedo_jacobian_mass_content_of_rain_in_atmosphere_layer_CH
         case (var_clr_wp)
            hofxdiags%geovals(jvar)%nval = n_Layers
            allocate(hofxdiags%geovals(jvar)%vals(hofxdiags%geovals(jvar)%nval,n_Profiles))
            hofxdiags%geovals(jvar)%vals = missing
            jspec = ufo_vars_getindex(conf%Clouds(:,1), var_clr_wp)
            do jprofile = 1, n_Profiles
               if (.not.Options(jprofile)%Skip_Profile) then
                  do jlevel = 1, hofxdiags%geovals(jvar)%nval
                     conv_albedo = PI/rts(jchannel,jprofile)%Solar_irradiance  ! Albedo conversion factor
                     hofxdiags%geovals(jvar)%vals(jlevel,jprofile) = conv_albedo *     &
                        atm_K(jchannel,jprofile) % Cloud(jspec) % Water_Content(jlevel)
                  end do
               end if
            end do

         ! variable: albedo_jacobian_mass_content_of_graupel_in_atmosphere_layer_CH
         case (var_clg_wp)
            hofxdiags%geovals(jvar)%nval = n_Layers
            allocate(hofxdiags%geovals(jvar)%vals(hofxdiags%geovals(jvar)%nval,n_Profiles))
            hofxdiags%geovals(jvar)%vals = missing
            jspec = ufo_vars_getindex(conf%Clouds(:,1), var_clg_wp)
            do jprofile = 1, n_Profiles
               if (.not.Options(jprofile)%Skip_Profile) then
                  do jlevel = 1, hofxdiags%geovals(jvar)%nval
                     conv_albedo = PI/rts(jchannel,jprofile)%Solar_irradiance  ! Albedo conversion factor
                     hofxdiags%geovals(jvar)%vals(jlevel,jprofile) = conv_albedo *     &
                        atm_K(jchannel,jprofile) % Cloud(jspec) % Water_Content(jlevel)
                  end do
               end if
            end do

         ! variable: albedo_jacobian_mass_content_of_hail_in_atmosphere_layer_CH
         case (var_clh_wp)
            hofxdiags%geovals(jvar)%nval = n_Layers
            allocate(hofxdiags%geovals(jvar)%vals(hofxdiags%geovals(jvar)%nval,n_Profiles))
            hofxdiags%geovals(jvar)%vals = missing
            jspec = ufo_vars_getindex(conf%Clouds(:,1), var_clh_wp)
            do jprofile = 1, n_Profiles
               if (.not.Options(jprofile)%Skip_Profile) then
                  do jlevel = 1, hofxdiags%geovals(jvar)%nval
                     conv_albedo = PI/rts(jchannel,jprofile)%Solar_irradiance  ! Albedo conversion factor
                     hofxdiags%geovals(jvar)%vals(jlevel,jprofile) = conv_albedo *     &
                        atm_K(jchannel,jprofile) % Cloud(jspec) % Water_Content(jlevel)
                  end do
               end if
            end do

         ! variable: albedo_jacobian_surface_emissivity_CH (nval=1)
         case (var_sfc_emiss)
            hofxdiags%geovals(jvar)%nval = 1
            allocate(hofxdiags%geovals(jvar)%vals(hofxdiags%geovals(jvar)%nval,n_Profiles))
            hofxdiags%geovals(jvar)%vals = missing
            do jprofile = 1, n_Profiles
               if (.not.Options(jprofile)%Skip_Profile) then
                  conv_albedo = PI/rts(jchannel,jprofile)%Solar_irradiance  ! Albedo conversion factor
                  hofxdiags%geovals(jvar)%vals(1,jprofile) = conv_albedo *     &
                     rts_K(jchannel,jprofile) % surface_emissivity
               end if
            end do

         case default
            write(err_msg,*) 'ufo_radiancecrtm_simobs: //&
                              & ObsDiagnostic is unsupported, ', &
                              & hofxdiags%variables(jvar)
            !call abor1_ftn(err_msg)
            err_stat = 1
      end select
#endif
   else
      write(err_msg,*) 'ufo_radiancecrtm_simobs: //&
                        & ObsDiagnostic is unsupported, ', &
                        & hofxdiags%variables(jvar)
      !call abor1_ftn(err_msg)
      err_stat = 1
   end if
end do

end subroutine ufo_crtm_passive_diag

! ------------------------------------------------------------------------------

end module ufo_crtm_passive_mod
