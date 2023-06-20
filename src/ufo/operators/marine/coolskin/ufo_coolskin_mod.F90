! (C) Copyright 2017-2018 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

!> Fortran coolskin module for observation operator

module ufo_coolskin_mod

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

!> Fortran derived type for the observation type
 type, extends(ufo_basis), public :: ufo_coolskin
 private
 contains
   procedure :: setup  => ufo_coolskin_setup
   procedure :: delete => ufo_coolskin_delete
   procedure :: simobs => ufo_coolskin_simobs
 end type ufo_coolskin

contains

! ------------------------------------------------------------------------------
subroutine ufo_coolskin_setup(self, f_conf)
implicit none
class(ufo_coolskin), intent(inout)    :: self
type(fckit_configuration), intent(in) :: f_conf

end subroutine ufo_coolskin_setup

! ------------------------------------------------------------------------------
subroutine ufo_coolskin_delete(self)
implicit none
class(ufo_coolskin), intent(inout) :: self

end subroutine ufo_coolskin_delete

! ------------------------------------------------------------------------------
subroutine ufo_coolskin_simobs(self, geovals, hofx, obss)

use ufo_coolskin_sim_mod
implicit none
    class(ufo_coolskin), intent(in)   :: self
    type(ufo_geovals),  intent(in)    :: geovals
    real(kind_real),     intent(inout) :: hofx(:)
    type(c_ptr), value, intent(in)    :: obss

    character(len=*), parameter :: myname_="ufo_coolskin_simobs"
    character(max_string) :: err_msg
    type(ufo_geoval), pointer :: S_ns,H_I,H_s,R_nl,Td,u
    integer :: obss_nlocs
    integer :: iobs
    real(c_double) :: missing
    real(kind_real) :: dTc

    ! Set missing flag
    missing = missing_value(missing)
    
    ! check if nobs is consistent in geovals & hofx
    obss_nlocs = obsspace_get_nlocs(obss)

    !nlocs = size(hofx,1)
    if (geovals%nlocs /= size(hofx,1)) then
       write(err_msg,*) myname_, ' error: nobs inconsistent!'
       call abor1_ftn(err_msg)
    endif

    ! check if coolskin input variables are in geovals and get them

    call ufo_geovals_get_var(geovals, var_ocn_sst, Td)
    call ufo_geovals_get_var(geovals, var_sw_rad , R_nl )	
    call ufo_geovals_get_var(geovals, var_latent_heat , H_I )
    call ufo_geovals_get_var(geovals, var_sens_heat , H_s )	
    call ufo_geovals_get_var(geovals, var_lw_rad , S_ns )
    call ufo_geovals_get_var(geovals, var_sea_fric_vel , u )
    
    ! simulated obs, hofx(iobs)=Ts 
    do iobs = 1, obss_nlocs
      ! check for missing values
      ! (the atmospheric fields *shouldn't* be masked, so dont
      ! bother checking)
      if (Td%vals(1, iobs) == missing .or. &
          u%vals(1,iobs) == missing) then
        hofx(iobs) = missing
        cycle
      end if

       call ufo_coolskin_sim(hofx(iobs),&
                             dTc,&
                             S_ns%vals(1,iobs),&
                             H_I%vals(1,iobs),&
                             H_s%vals(1,iobs),&
                             R_nl%vals(1,iobs),&
                             Td%vals(1,iobs),&
                             u%vals(1,iobs))
    enddo

    
  end subroutine ufo_coolskin_simobs


!--------------------------------------------------

end module ufo_coolskin_mod
