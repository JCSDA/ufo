! (C) Copyright 2022 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

!> Fortran oasim module for observation operator

module ufo_oasim_mod

 use oasim_mod, only: oasim

 use fckit_configuration_module, only: fckit_configuration
 use iso_c_binding
 use kinds

 use ufo_geovals_mod, only: ufo_geovals, ufo_geoval, ufo_geovals_get_var
 use ufo_basis_mod, only: ufo_basis
 use ufo_vars_mod

 use obsspace_mod
 use missing_values_mod


 implicit none
 private

 integer, parameter :: max_string=800

 ! Oasim oject
 type, public :: ufo_oasim
 private
   type(oasim) :: oasim_
 contains
   procedure :: setup  => ufo_oasim_setup
   procedure :: delete => ufo_oasim_delete
   procedure :: simobs => ufo_oasim_simobs

 end type ufo_oasim

contains

! ------------------------------------------------------------------------------
subroutine ufo_oasim_setup(self, f_conf)
class(ufo_oasim),          intent(inout)  :: self
type(fckit_configuration), intent(in)     :: f_conf

character(len=*), parameter :: myname_="ufo_oasim_simobs"
character(max_string) :: err_msg
character(len=:), allocatable :: str

!Path to coefficient files                                                                                                       
call f_conf%get_or_die("CoefficientPath",str)
call self%oasim_%create(str)

end subroutine ufo_oasim_setup

! ------------------------------------------------------------------------------
subroutine ufo_oasim_delete(self)
implicit none
class(ufo_oasim), intent(inout) :: self

call self%oasim_%delete()

end subroutine ufo_oasim_delete

! ------------------------------------------------------------------------------
subroutine ufo_oasim_simobs(self, geovals, hofx, obss)

implicit none

    class(ufo_oasim),          intent(in)    :: self
    type(ufo_geovals),         intent(in)    :: geovals
    real(c_double),            intent(inout) :: hofx(:,:)
    type(c_ptr),  value,       intent(in)    :: obss

    character(len=*), parameter :: myname_="ufo_oasim_simobs"
    character(max_string) :: err_msg

    logical :: is_midnight 
    integer :: km, day_of_year
    real(kind_real) :: cosz, dt
    type(ufo_geoval), pointer :: slp, wspd, ozone, wvapor, rh, cov, cldtau, clwp, cldre
    type(ufo_geoval), pointer :: ta_in, wa_in, asym, dh, cdet, pic, cdc, diatom, chloro, cyano
    type(ufo_geoval), pointer :: cocco, dino, phaeo

    real (kind_real), allocatable :: tirrq(:,:), cdomabsq(:,:), avgq(:,:)
    integer :: obss_nlocs
    integer :: iobs
    real(c_double) :: missing

    ! Set missing flag
    missing = missing_value(missing)
    
    ! check if nobs is consistent in geovals & hofx
    obss_nlocs = obsspace_get_nlocs(obss)

    !nlocs = size(hofx,1)
    if (geovals%nlocs /= size(hofx,1)) then
       write(err_msg,*) myname_, ' error: nobs inconsistent!'
       call abor1_ftn(err_msg)
    endif

    ! check if oasim input variables are in geovals and get them
    call ufo_geovals_get_var(geovals, var_pmsl , slp )
    call ufo_geovals_get_var(geovals, var_sfc_wspeed , wspd ) 
    call ufo_geovals_get_var(geovals, var_oz_thick , ozone )
    call ufo_geovals_get_var(geovals, var_water_vapor , wvapor )
    call ufo_geovals_get_var(geovals, var_rh , rh )
    call ufo_geovals_get_var(geovals, var_cldfrac , cov )
    call ufo_geovals_get_var(geovals, var_cld_tau , cldtau )
    call ufo_geovals_get_var(geovals, var_cld_lwp , clwp )
    call ufo_geovals_get_var(geovals, var_clwefr , cldre )
    call ufo_geovals_get_var(geovals, var_aerosol_tau , ta_in )
    call ufo_geovals_get_var(geovals, var_scat_albedo , wa_in )
    call ufo_geovals_get_var(geovals, var_asym_par , asym )
    call ufo_geovals_get_var(geovals, var_ocn_lay_thick , dh )
    call ufo_geovals_get_var(geovals, var_carb_det , cdet )
    call ufo_geovals_get_var(geovals, var_inorg_carb , pic )
    call ufo_geovals_get_var(geovals, var_dis_carb , cdc )
    call ufo_geovals_get_var(geovals, var_diatom_conc , diatom )
    call ufo_geovals_get_var(geovals, var_chloro_conc , chloro )
    call ufo_geovals_get_var(geovals, var_cyano_conc , cyano )
    call ufo_geovals_get_var(geovals, var_cocco_conc , cocco )
    call ufo_geovals_get_var(geovals, var_phaeo_conc , phaeo )
    call ufo_geovals_get_var(geovals, var_dino_conc , dino )

    ! To Do change this part later to read form obs file
    is_midnight = .FALSE.
    km = size(phaeo%vals,dim=1)
    day_of_year = 90
    cosz = 0.5
    dt = 86400/2.0
    
    allocate(tirrq(obss_nlocs,km))
    allocate(cdomabsq(obss_nlocs,km))
    allocate(avgq(obss_nlocs,km))

    do iobs = 1, obss_nlocs
       ! check if the ocean thickness is positive (valid)
       if (dh%vals(1,iobs) > 0) then

          call self%oasim_%run(km, dt, is_midnight, day_of_year, &
               cosz, slp%vals(1,iobs), wspd%vals(1,iobs), ozone%vals(1,iobs), wvapor%vals(1,iobs), &
               rh%vals(1,iobs), cov%vals(1,iobs), cldtau%vals(1,iobs), clwp%vals(1,iobs), &
               cldre%vals(1,iobs), ta_in%vals(:,iobs), wa_in%vals(:,iobs), asym%vals(:,iobs), &
               dh%vals(:,iobs), cdet%vals(:,iobs), pic%vals(:,iobs), cdc%vals(:,iobs), &
               diatom%vals(:,iobs), chloro%vals(:,iobs), cyano%vals(:,iobs), &
               cocco%vals(:,iobs), dino%vals(:,iobs), phaeo%vals(:,iobs), tirrq(iobs,:), cdomabsq(iobs,:), avgq(iobs,:))
          hofx(iobs,:)=sum(avgq(iobs,:)) 
       endif
    enddo

end subroutine ufo_oasim_simobs

end module ufo_oasim_mod
