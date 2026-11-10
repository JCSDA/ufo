! (C) Copyright 2017-2023 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

!> Fortran module to handle radiancecrtm active observations

module ufo_crtm_active_mod

 use CRTM_Module

 use fckit_configuration_module, only: fckit_configuration
 use, intrinsic :: iso_c_binding
 use kinds
 use missing_values_mod

 use obsspace_mod

 use ufo_geovals_mod, only: ufo_geovals, ufo_geoval, ufo_geovals_get_var
 use ufo_vars_mod
 use ufo_crtm_utils_mod

 use ufo_constants_mod, only: deg2rad, threshold_reflectivity
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
type(c_ptr), value,       intent(in) :: obss        !ObsSpace

real(c_double) :: missing
integer        :: l, m, jlayer
integer, allocatable   :: modelLayer(:)

! allocate observation elevation for reflectivity profiles
allocate( modelLayer( n_Profiles ) )

! Put simulated reflectivity into hofx
! ----------------------------------------------

! Set missing value
missing = missing_value(missing)

!Set to missing, then retrieve non-missing profiles
hofx = missing

! Get the level number for reflectivity profiles
call obsspace_get_db(obss, "MetaData", "modelLayer", modelLayer)
modelLayer = modelLayer + 1 ! convert from 0 index to 1

do m = 1, n_Profiles
   if (Options(m)%Skip_Profile) cycle
   jlayer = modelLayer(m)
   do l = 1, n_Channels
      if (abs(rts(l,m)%Reflectivity_Attenuated(jlayer)) < threshold_reflectivity) then
         hofx(l,m) = rts(l,m)%Reflectivity_Attenuated(jlayer)
      end if
   end do
end do

end subroutine ufo_crtm_active_sim


subroutine ufo_crtm_active_diag(rts, rts_K, atm, atm_K, sfc_K, conf, isensor, Options,&
             channels, geovals, obss, nvars, nlocs, n_Profiles, n_Layers, xstr_diags, ystr_diags,&
             ch_diags, hofxdiags, err_stat)
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
integer, intent(out)  :: err_stat

! Local Variables
character(max_string) :: err_msg
integer :: jvar, jprofile, jchannel, jspec, jlayer
real(c_double) :: missing
integer, allocatable   :: modelLayer(:)

! allocate observation elevation for reflectivity profiles
allocate( modelLayer( n_Profiles ))

! Set missing value
missing = missing_value(missing)

! Get the level number for reflectivity profiles
call obsspace_get_db(obss, "MetaData", "modelLayer", modelLayer)
modelLayer = modelLayer + 1 ! convert from 0 index to 1

! Put simulated diagnostics into hofxdiags
! We need to call the routines for passive instrument as well when dealing with active obs
call ufo_crtm_passive_diag( &
    rts, rts_K, atm, atm_K, sfc_K, conf, isensor, Options, channels, geovals, &
    obss, nvars, nlocs, n_Profiles, n_Layers, xstr_diags, ystr_diags, ch_diags, &
    hofxdiags, err_stat)

! ----------------------------------------------
do jvar = 1, hofxdiags%nvar
   if (len(trim(hofxdiags%variables(jvar))) < 1) cycle

   if (ch_diags(jvar) > 0) then
      if (count(channels == ch_diags(jvar)) /= 1) then
         write(err_msg,"(a,i0,a,*(1x,i0))") "ufo_crtm_active_diag: ch_diags(jvar) = ", ch_diags(jvar), &
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

   hofxdiags%geovals(jvar)%nval = n_Layers

   !============================================
   ! Diagnostics used for QC and bias correction
   !============================================
   if (xstr_diags(jvar) == "") then
      ! forward h(x) diags (one value per profile, at the observation's model layer)
      select case(ystr_diags(jvar))
         ! variable: reflectivity_CH Non attenuated
         case (var_rad_refl)
            call alloc_diag(1)
            do jprofile = 1, n_Profiles
               if (Options(jprofile)%Skip_Profile) cycle
               jlayer = modelLayer(jprofile)
               if (abs(rts(jchannel,jprofile)%Reflectivity(jlayer)) < &
                   threshold_reflectivity) then
                  hofxdiags%geovals(jvar)%vals(1,jprofile) = &
                     rts(jchannel,jprofile) % Reflectivity(jlayer)
               end if
            end do

         ! variable: attenuated_reflectivity_CH
         case (var_rad_refl_att)
            call alloc_diag(1)
            do jprofile = 1, n_Profiles
               if (Options(jprofile)%Skip_Profile) cycle
               jlayer = modelLayer(jprofile)
               if (abs(rts(jchannel,jprofile)%Reflectivity_Attenuated(jlayer)) < &
                   threshold_reflectivity) then
                  hofxdiags%geovals(jvar)%vals(1,jprofile) = &
                     rts(jchannel,jprofile) % Reflectivity_Attenuated(jlayer)
               end if
            end do

         case default
            ! unsupported forward diagnostic: fill with missing
            call alloc_diag(1)
      end select
   else if ((ystr_diags(jvar) == var_rad_refl) .or. (ystr_diags(jvar) == var_rad_refl_att)) then
      ! var_rad_refl jacobians
      select case (xstr_diags(jvar))
         ! variable: reflectivity_jacobian_mass_content_of_<species>_in_atmosphere_layer_CH,
         ! for cloud liquid water, cloud ice, rain, snow, graupel and hail
         case (var_clw_wp, var_cli_wp, var_clr_wp, var_cls_wp, var_clg_wp, var_clh_wp)
            call alloc_diag(n_Layers)
            jspec = ufo_vars_getindex(conf%Clouds(:,1), xstr_diags(jvar))
            do jprofile = 1, n_Profiles
               if (Options(jprofile)%Skip_Profile) cycle
               hofxdiags%geovals(jvar)%vals(:,jprofile) = &
                  atm_K(jchannel,jprofile) % Cloud(jspec) % Water_Content(1:n_Layers)
            end do

         case default
            err_msg = "ufo_crtm_active_diag: ObsDiagnostic is unsupported, " // &
                      hofxdiags%variables(jvar)
            call fckit_log%error(err_msg)
            err_stat = 1
      end select
   else
      err_msg = "ufo_crtm_active_diag: ObsDiagnostic is unsupported, " // &
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

end subroutine ufo_crtm_active_diag

! ------------------------------------------------------------------------------

end module ufo_crtm_active_mod
