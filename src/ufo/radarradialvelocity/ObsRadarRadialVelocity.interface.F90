! (C) Copyright 2017-2018 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

!> Fortran radarradialvelocity module for functions on the interface between C++ and Fortran
!  to handle observation operators

module ufo_radarradialvelocity_mod_c

  use fckit_configuration_module, only: fckit_configuration
  use iso_c_binding
  use ufo_radarradialvelocity_mod 
  use string_f_c_mod
  use ufo_geovals_mod,   only: ufo_geovals
  use ufo_geovals_mod_c, only: ufo_geovals_registry
  implicit none
  private

  ! ------------------------------------------------------------------------------
#define LISTED_TYPE ufo_radarradialvelocity

  !> Linked list interface - defines registry_t type
#include "../linkedList_i.f"

  !> Global registry
  type(registry_t) :: ufo_radarradialvelocity_registry

  ! ------------------------------------------------------------------------------

contains

  ! ------------------------------------------------------------------------------
  !> Linked list implementation
#include "../linkedList_c.f"

! ------------------------------------------------------------------------------

subroutine ufo_radarradialvelocity_setup_c(c_key_self, c_conf, c_varconf, c_varlist) bind(c,name='ufo_radarradialvelocity_setup_f90')
use ufo_vars_mod
implicit none
integer(c_int), intent(inout) :: c_key_self
type(c_ptr),    intent(in)    :: c_conf
type(c_ptr), intent(in) :: c_varconf ! config with variables to be simulated
type(c_ptr), intent(in), value :: c_varlist
character(len=MAXVARLEN), dimension(:), allocatable :: vars

type(ufo_radarradialvelocity), pointer :: self
type(fckit_configuration) :: f_conf, f_varconf

call ufo_radarradialvelocity_registry%setup(c_key_self, self)
f_conf = fckit_configuration(c_conf)
f_varconf = fckit_configuration(c_varconf)

call ufo_vars_read(f_varconf, vars)
call self%setup(f_conf, vars)
deallocate(vars)

!> Update C++ ObsOperator with input variable list
call f_c_push_string_varlist(c_varlist, self%varin)

end subroutine ufo_radarradialvelocity_setup_c

! ------------------------------------------------------------------------------

subroutine ufo_radarradialvelocity_delete_c(c_key_self) bind(c,name='ufo_radarradialvelocity_delete_f90')
implicit none
integer(c_int), intent(inout) :: c_key_self
    
type(ufo_radarradialvelocity), pointer :: self

call ufo_radarradialvelocity_registry%delete(c_key_self, self)

end subroutine ufo_radarradialvelocity_delete_c

! ------------------------------------------------------------------------------

subroutine ufo_radarradialvelocity_simobs_c(c_key_self, c_key_geovals, c_obsspace, c_nvars, c_nlocs, &
                                c_hofx) bind(c,name='ufo_radarradialvelocity_simobs_f90')

implicit none
integer(c_int), intent(in) :: c_key_self
integer(c_int), intent(in) :: c_key_geovals
type(c_ptr), value, intent(in) :: c_obsspace
integer(c_int), intent(in)     :: c_nvars, c_nlocs
real(c_double), intent(inout)  :: c_hofx(c_nvars, c_nlocs)

type(ufo_radarradialvelocity), pointer :: self
type(ufo_geovals),       pointer :: geovals

call ufo_radarradialvelocity_registry%get(c_key_self, self)
call ufo_geovals_registry%get(c_key_geovals, geovals)
call self%simobs(geovals, c_obsspace, c_nvars, c_nlocs, c_hofx)

end subroutine ufo_radarradialvelocity_simobs_c

! ------------------------------------------------------------------------------

end module ufo_radarradialvelocity_mod_c
