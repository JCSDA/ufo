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

call self%create('../../oasim/test/data')

end subroutine ufo_oasim_setup

! ------------------------------------------------------------------------------
subroutine ufo_oasim_delete(self)
implicit none
class(ufo_oasim), intent(inout) :: self

call self%cancel()

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
    call ufo_geovals_get_var(geovals, var_layer_thick , dh )
    call ufo_geovals_get_var(geovals, var_carb_det , cdet )
    call ufo_geovals_get_var(geovals, var_inorg_carb , pic )
    call ufo_geovals_get_var(geovals, var_dis_carb , cdc )
    call ufo_geovals_get_var(geovals, var_diatom_conc , diatom )
    call ufo_geovals_get_var(geovals, var_chloro_conc , chloro )
    call ufo_geovals_get_var(geovals, var_cyano_conc , cyano )
    call ufo_geovals_get_var(geovals, var_cocco_conc , cocco )
    call ufo_geovals_get_var(geovals, var_phaeo_conc , phaeo )
    call ufo_geovals_get_var(geovals, var_dino_conc , dino )

    is_midnight = .FALSE.
    km = 14
    day_of_year = 90
    cosz = 0.5
    dt = 86400/2.0

    do iobs = 1, obss_nlocs

    call self%run(km, dt, is_midnight, day_of_year, &
            cosz, slp%vals(iobs,1), wspd%vals(iobs,1), ozone%vals(iobs,1), wvapor%vals(iobs,1), &
            rh%vals(iobs,1), cov%vals(iobs,1), cldtau%vals(iobs,1), clwp%vals(iobs,1), &
            cldre%vals(iobs,1), ta_in%vals(iobs,:), wa_in%vals(iobs,:), asym%vals(iobs,:), &
            dh%vals(iobs,:), cdet%vals(iobs,:), pic%vals(iobs,:), cdc%vals(iobs,:), &
            diatom%vals(iobs,:), chloro%vals(iobs,:), cyano%vals(iobs,:), &
            cocco%vals(iobs,:), dino%vals(iobs,:), phaeo%vals(iobs,:), tirrq(iobs,:), cdomabsq(iobs,:), avgq(iobs,:))

    enddo
    
end subroutine ufo_oasim_simobs

end module ufo_oasim_mod
