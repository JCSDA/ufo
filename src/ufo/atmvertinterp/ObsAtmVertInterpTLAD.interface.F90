! (C) Copyright 2017-2018 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

!> Fortran module to handle atmvertinterp observations

module ufo_atmvertinterp_tlad_mod_c

  use fckit_configuration_module, only: fckit_configuration
  use ufo_atmvertinterp_tlad_mod
  use ufo_geovals_mod_c, only: ufo_geovals_registry
  use ufo_geovals_mod,   only: ufo_geovals
  use string_f_c_mod
  implicit none

  private

#define LISTED_TYPE ufo_atmvertinterp_tlad

  !> Linked list interface - defines registry_t type
#include "oops/util/linkedList_i.f"

  !> Global registry
  type(registry_t) :: ufo_atmvertinterp_tlad_registry

  ! ------------------------------------------------------------------------------
contains
  ! ------------------------------------------------------------------------------
  !> Linked list implementation
#include "oops/util/linkedList_c.f"

! ------------------------------------------------------------------------------

subroutine ufo_atmvertinterp_tlad_setup_c(c_key_self, c_conf, c_varconf, c_varlist) bind(c,name='ufo_atmvertinterp_tlad_setup_f90')
use ufo_vars_mod
implicit none
integer(c_int), intent(inout) :: c_key_self
type(c_ptr), intent(in) :: c_conf
type(c_ptr), intent(in) :: c_varconf ! config with variables to be simulated
type(c_ptr), intent(in), value :: c_varlist
character(len=MAXVARLEN), dimension(:), allocatable :: vars

type(ufo_atmvertinterp_tlad), pointer :: self
type(fckit_configuration) :: f_conf, f_varconf

call ufo_atmvertinterp_tlad_registry%setup(c_key_self, self)
f_conf = fckit_configuration(c_conf)
f_varconf = fckit_configuration(c_varconf)

call ufo_vars_read(f_varconf, vars)
call self%setup(f_conf, vars)
deallocate(vars)

!> Update C++ ObsOperator with input variable list
call f_c_push_string_varlist(c_varlist, self%varin)

end subroutine ufo_atmvertinterp_tlad_setup_c

! ------------------------------------------------------------------------------

subroutine ufo_atmvertinterp_tlad_delete_c(c_key_self) bind(c,name='ufo_atmvertinterp_tlad_delete_f90')
implicit none
integer(c_int), intent(inout) :: c_key_self

type(ufo_atmvertinterp_tlad), pointer :: self

call ufo_atmvertinterp_tlad_registry%delete(c_key_self, self)

end subroutine ufo_atmvertinterp_tlad_delete_c

! ------------------------------------------------------------------------------

subroutine ufo_atmvertinterp_tlad_settraj_c(c_key_self, c_key_geovals, c_obsspace) bind(c,name='ufo_atmvertinterp_tlad_settraj_f90')
implicit none
integer(c_int), intent(in) :: c_key_self
integer(c_int), intent(in) :: c_key_geovals
type(c_ptr), value, intent(in) :: c_obsspace

type(ufo_atmvertinterp_tlad), pointer :: self
type(ufo_geovals),            pointer :: geovals

character(len=*), parameter :: myname_="ufo_atmvertinterp_tlad_settraj_c"

call ufo_atmvertinterp_tlad_registry%get(c_key_self, self)
call ufo_geovals_registry%get(c_key_geovals, geovals)

call self%settraj(geovals, c_obsspace)

end subroutine ufo_atmvertinterp_tlad_settraj_c

! ------------------------------------------------------------------------------

subroutine ufo_atmvertinterp_simobs_tl_c(c_key_self, c_key_geovals, c_obsspace, c_nvars, c_nlocs, c_hofx) &
           bind(c,name='ufo_atmvertinterp_simobs_tl_f90')
implicit none
integer(c_int), intent(in) :: c_key_self
integer(c_int), intent(in) :: c_key_geovals
type(c_ptr), value, intent(in) :: c_obsspace
integer(c_int), intent(in) :: c_nvars, c_nlocs
real(c_double), intent(inout) :: c_hofx(c_nvars, c_nlocs)

type(ufo_atmvertinterp_tlad), pointer :: self
type(ufo_geovals),            pointer :: geovals

character(len=*), parameter :: myname_="ufo_atmvertinterp_simobs_tl_c"

call ufo_atmvertinterp_tlad_registry%get(c_key_self, self)
call ufo_geovals_registry%get(c_key_geovals, geovals)

call self%simobs_tl(geovals, c_obsspace, c_nvars, c_nlocs, c_hofx)

end subroutine ufo_atmvertinterp_simobs_tl_c

! ------------------------------------------------------------------------------

subroutine ufo_atmvertinterp_simobs_ad_c(c_key_self, c_key_geovals, c_obsspace, c_nvars, c_nlocs, c_hofx) &
           bind(c,name='ufo_atmvertinterp_simobs_ad_f90')
implicit none
integer(c_int), intent(in) :: c_key_self
integer(c_int), intent(in) :: c_key_geovals
type(c_ptr), value, intent(in) :: c_obsspace
integer(c_int), intent(in) :: c_nvars, c_nlocs
real(c_double), intent(in) :: c_hofx(c_nvars, c_nlocs)

type(ufo_atmvertinterp_tlad), pointer :: self
type(ufo_geovals),            pointer :: geovals

character(len=*), parameter :: myname_="ufo_atmvertinterp_simobs_ad_c"

call ufo_atmvertinterp_tlad_registry%get(c_key_self, self)
call ufo_geovals_registry%get(c_key_geovals, geovals)

call self%simobs_ad(geovals, c_obsspace, c_nvars, c_nlocs, c_hofx)

end subroutine ufo_atmvertinterp_simobs_ad_c

! ------------------------------------------------------------------------------

end module ufo_atmvertinterp_tlad_mod_c
