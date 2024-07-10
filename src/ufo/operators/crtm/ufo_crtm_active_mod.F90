! (C) Copyright 2017-2023 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

!> Fortran module to handle radiancecrtm active observations

module ufo_crtm_active_mod

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
 use ufo_crtm_passive_mod

 implicit none
 private


 public ufo_crtm_active_sim
 public ufo_crtm_active_diag

contains


subroutine ufo_crtm_active_sim(rts, Options, nvars, nlocs, n_Profiles, n_Channels, hofx, obss)
implicit none
integer(c_size_t),        intent(in) :: nvars
integer(c_size_t),        intent(in) :: nlocs
integer, intent(in) :: n_Profiles, n_Channels
real(c_double),        intent(inout) :: hofx(nvars, nlocs) !h(x) to return
type(CRTM_RTSolution_type), intent(in) :: rts(:,:)  ! n_channels, n_profiles
type(CRTM_Options_type),  intent(in) :: Options(:)
type(c_ptr), value,       intent(in) :: obss         !ObsSpace

real(c_double) :: missing
integer        :: l, m, jlayer
integer, allocatable   :: obs_layer(:)

! allocate observation elevation for reflectivity profiles
allocate( obs_layer( n_Profiles ))

! Put simulated reflectivity into hofx
! ----------------------------------------------


! Set missing value
missing = missing_value(missing)

!Set to missing, then retrieve non-missing profiles
hofx = missing

! Get the level number for reflectivity profiles
call obsspace_get_db(obss, "MetaData", "Layer", obs_layer)

do m = 1, n_Profiles
   if (.not.Options(m)%Skip_Profile) then
      jlayer = obs_layer(m)
      do l = 1, n_Channels
        if (abs(rts(l,m)%Reflectivity_Attenuated(jlayer)) < threshold_reflectivity) then
           hofx(l,m) = rts(l,m)%Reflectivity_Attenuated(jlayer)
        endif
      end do
   end if
end do


end subroutine ufo_crtm_active_sim


subroutine ufo_crtm_active_diag(rts, rts_K, atm, atm_K, sfc_K, conf, n_Sensor, Options,&
             channels, geovals, obss, nvars, nlocs, n_Profiles, n_Layers, xstr_diags, ystr_diags,&
             ch_diags, hofxdiags, err_stat)
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
integer, intent(out)  :: err_stat

! Local Variables
!character(*), parameter :: PROGRAM_NAME = 'ufo_crtm_active_diag'
!character(255) :: message, version
character(max_string) :: err_msg
!integer        :: err_stat, alloc_stat
!integer        :: l, m, n
integer :: jvar, jprofile, jlevel, jchannel, ichannel, jspec, jlayer
real(c_double) :: missing
real(kind_real) :: total_od, secant_term, wfunc_max
real(kind_real), allocatable :: TmpVar(:)
real(kind_real), allocatable :: Tao(:)
real(kind_real), allocatable :: Wfunc(:)
integer, allocatable   :: obs_layer(:)

! allocate observation elevation for reflectivity profiles
allocate( obs_layer( n_Profiles ))

! Set missing value
missing = missing_value(missing)

! Get the level number for reflectivity profiles
call obsspace_get_db(obss, "MetaData", "Layer", obs_layer)

! Put simulated diagnostics into hofxdiags
! We need to call the routines for passive instrument as well when dealing with active obs
call ufo_crtm_passive_diag(rts, rts_K, atm, atm_K, sfc_K, conf, n_Sensor, Options, channels, geovals, obss, nvars, nlocs, n_Profiles, n_Layers, xstr_diags, ystr_diags, ch_diags, hofxdiags, err_stat)

! ----------------------------------------------
do jvar = 1, hofxdiags%nvar
   if (len(trim(hofxdiags%variables(jvar))) < 1) cycle

   if (ch_diags(jvar) > 0) then
      if (size(pack(channels,channels==ch_diags(jvar))) /= 1) then
         write(err_msg,*) 'ufo_crtm_active_diags: mismatch between// &
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

   hofxdiags%geovals(jvar)%nval = n_Layers

   !============================================
   ! Diagnostics used for QC and bias correction
   !============================================
   if (cmp_strings(xstr_diags(jvar), "")) then
      ! forward h(x) diags
      select case(ystr_diags(jvar))
         ! variable: reflectivity_CH Non attenuated
         case (var_rad_refl)
            ! SET n_Layers = 1 => hofxdiags%geovals(jvar)%nval
            hofxdiags%geovals(jvar)%nval = 1 !n_Layers
            allocate(hofxdiags%geovals(jvar)%vals(hofxdiags%geovals(jvar)%nval, n_Profiles))
            hofxdiags%geovals(jvar)%vals = missing
            do jprofile = 1, n_Profiles
               if (.not.Options(jprofile)%Skip_Profile) then
                  jlayer = obs_layer(jprofile)
                  do jlevel = 1, hofxdiags%geovals(jvar)%nval
                     if (abs(rts(jchannel,jprofile)%Reflectivity(jlayer)) < threshold_reflectivity) then 
                          hofxdiags%geovals(jvar)%vals(jlevel,jprofile) = &
                              rts(jchannel,jprofile) % Reflectivity(jlayer)
                     endif
                  end do
               end if
            end do

         ! variable: attenuated_reflectivity_CH
         case (var_rad_refl_att)
            ! SET n_Layers = 1 => hofxdiags%geovals(jvar)%nval
            hofxdiags%geovals(jvar)%nval = 1 !n_Layers
            allocate(hofxdiags%geovals(jvar)%vals(hofxdiags%geovals(jvar)%nval, n_Profiles))
            hofxdiags%geovals(jvar)%vals = missing
            do jprofile = 1, n_Profiles
               if (.not.Options(jprofile)%Skip_Profile) then
                  jlayer = obs_layer(jprofile)
                  do jlevel = 1, hofxdiags%geovals(jvar)%nval
                     if (abs(rts(jchannel,jprofile)%Reflectivity_Attenuated(jlayer)) < threshold_reflectivity) then
                         hofxdiags%geovals(jvar)%vals(jlevel,jprofile) = &
                             rts(jchannel,jprofile) % Reflectivity_Attenuated(jlayer)
                     endif
                  end do
               end if
            end do

         case default
            write(err_msg,*) 'ufo_crtm_active_diags: //&
                              & ObsDiagnostic is unsupported, ', &
                              & hofxdiags%variables(jvar)
            ! SET n_Layers = 1 => hofxdiags%geovals(jvar)%nval
            hofxdiags%geovals(jvar)%nval = 1
            allocate(hofxdiags%geovals(jvar)%vals(hofxdiags%geovals(jvar)%nval, n_Profiles))
            hofxdiags%geovals(jvar)%vals = missing
      end select
   else if ((ystr_diags(jvar) == var_rad_refl) .or. (ystr_diags(jvar) == var_rad_refl_att)) then
      ! var_rad_refl jacobians
      select case (xstr_diags(jvar))
         ! variable: reflectivity_jacobian_mass_content_of_cloud_liquid_water_in_atmosphere_layer_CH
         case (var_clw_wp)
            allocate(hofxdiags%geovals(jvar)%vals(n_Layers, n_Profiles))
            hofxdiags%geovals(jvar)%vals = missing
            jspec = ufo_vars_getindex(conf%Clouds(:,1), var_clw_wp)

            do jprofile = 1, n_Profiles
               if (.not.Options(jprofile)%Skip_Profile) then
                  do jlevel = 1, n_Layers
                     hofxdiags%geovals(jvar)%vals(jlevel,jprofile) = &
                        atm_K(jchannel,jprofile) % Cloud(jspec) % Water_Content(jlevel)
                  end do
               end if
            end do
         ! variable: reflectivity_jacobian_mass_content_of_cloud_ice_in_atmosphere_layer_CH
         case (var_cli_wp)
            allocate(hofxdiags%geovals(jvar)%vals(n_Layers, n_Profiles))
            hofxdiags%geovals(jvar)%vals = missing
            jspec = ufo_vars_getindex(conf%Clouds(:,1), var_cli_wp)

            do jprofile = 1, n_Profiles
               if (.not.Options(jprofile)%Skip_Profile) then
                  do jlevel = 1, n_Layers 
                     hofxdiags%geovals(jvar)%vals(jlevel,jprofile) = &
                        atm_K(jchannel,jprofile) % Cloud(jspec) % Water_Content(jlevel)
                  end do
               end if
            end do
         ! variable: reflectivity_jacobian_mass_content_of_rain_in_atmosphere_layer_CH
         case (var_clr_wp)
            allocate(hofxdiags%geovals(jvar)%vals(n_Layers, n_Profiles))
            hofxdiags%geovals(jvar)%vals = missing
            jspec = ufo_vars_getindex(conf%Clouds(:,1), var_clr_wp)

            do jprofile = 1, n_Profiles
               if (.not.Options(jprofile)%Skip_Profile) then
                  do jlevel = 1, n_Layers 
                     hofxdiags%geovals(jvar)%vals(jlevel,jprofile) = &
                        atm_K(jchannel,jprofile) % Cloud(jspec) % Water_Content(jlevel)
                  end do
               end if
            end do
         ! variable: reflectivity_jacobian_mass_content_of_snow_in_atmosphere_layer_CH
         case (var_cls_wp)
            allocate(hofxdiags%geovals(jvar)%vals(n_Layers, n_Profiles))
            hofxdiags%geovals(jvar)%vals = missing
            jspec = ufo_vars_getindex(conf%Clouds(:,1), var_cls_wp)

            do jprofile = 1, n_Profiles
               if (.not.Options(jprofile)%Skip_Profile) then
                  do jlevel = 1, n_Layers 
                     hofxdiags%geovals(jvar)%vals(jlevel,jprofile) = &
                        atm_K(jchannel,jprofile) % Cloud(jspec) % Water_Content(jlevel)
                  end do
               end if
            end do
         ! variable: reflectivity_jacobian_mass_content_of_graupel_in_atmosphere_layer_CH
         case (var_clg_wp)
            allocate(hofxdiags%geovals(jvar)%vals(n_Layers, n_Profiles))
            hofxdiags%geovals(jvar)%vals = missing
            jspec = ufo_vars_getindex(conf%Clouds(:,1), var_clg_wp)

            do jprofile = 1, n_Profiles
               if (.not.Options(jprofile)%Skip_Profile) then
                  do jlevel = 1, n_Layers 
                     hofxdiags%geovals(jvar)%vals(jlevel,jprofile) = &
                        atm_K(jchannel,jprofile) % Cloud(jspec) % Water_Content(jlevel)
                  end do
               end if
            end do
         ! variable: reflectivity_jacobian_mass_content_of_hail_in_atmosphere_layer_CH
         case (var_clh_wp)
            allocate(hofxdiags%geovals(jvar)%vals(n_Layers, n_Profiles))
            hofxdiags%geovals(jvar)%vals = missing
            jspec = ufo_vars_getindex(conf%Clouds(:,1), var_clh_wp)

            do jprofile = 1, n_Profiles
               if (.not.Options(jprofile)%Skip_Profile) then
                  do jlevel = 1, n_Layers 
                     hofxdiags%geovals(jvar)%vals(jlevel,jprofile) = &
                        atm_K(jchannel,jprofile) % Cloud(jspec) % Water_Content(jlevel)
                  end do
               end if
            end do

         case default
            write(err_msg,*) 'ufo_crtm_active_diags: //&
                              & ObsDiagnostic is unsupported, ', &
                              & hofxdiags%variables(jvar)
            !call abor1_ftn(err_msg)
            err_stat = 1
      end select
   else
      write(err_msg,*) 'ufo_crtm_active_diags: //&
                        & ObsDiagnostic is unsupported, ', &
                        & hofxdiags%variables(jvar)
      !call abor1_ftn(err_msg)
      err_stat = 1
   end if
end do


end subroutine ufo_crtm_active_diag

! ------------------------------------------------------------------------------

end module ufo_crtm_active_mod
