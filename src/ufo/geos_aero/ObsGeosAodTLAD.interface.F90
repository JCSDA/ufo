! (C) Copyright 2019 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

!> Fortran geosaod module for functions on the interface between C++ and Fortran
!  to handle tl/ad observation operators

module ufo_geosaod_tlad_mod_c

  use iso_c_binding
  use fckit_configuration_module, only: fckit_configuration
  use ufo_geosaod_tlad_mod
  use string_f_c_mod
  use ufo_geovals_mod_c, only: ufo_geovals_registry
  use ufo_geovals_mod,   only: ufo_geovals

  implicit none
  private

#define LISTED_TYPE ufo_geosaod_tlad

  !> Linked list interface - defines registry_t type
#include "../linkedList_i.f"

  !> Global registry
  type(registry_t) :: ufo_geosaod_tlad_registry

contains

  ! ------------------------------------------------------------------------------
  !> Linked list implementation
#include "../linkedList_c.f"

! ------------------------------------------------------------------------------
subroutine ufo_geosaod_tlad_setup_c(c_key_self, c_conf, c_varlist, c_nvars_out) bind(c,name='ufo_geosaod_tlad_setup_f90')

implicit none
integer(c_int), intent(inout) :: c_key_self
type(c_ptr), intent(in)    :: c_conf
integer(c_int), intent(in) :: c_nvars_out
type(c_ptr), intent(in), value :: c_varlist

type(ufo_geosaod_tlad), pointer :: self

type(fckit_configuration)  :: f_conf

f_conf = fckit_configuration(c_conf)


call ufo_geosaod_tlad_registry%setup(c_key_self, self)

call self%setup(f_conf, c_nvars_out)

!> Set vars
call f_c_push_string_varlist(c_varlist, self%varin)

end subroutine ufo_geosaod_tlad_setup_c

! ------------------------------------------------------------------------------

subroutine ufo_geosaod_tlad_delete_c(c_key_self) bind(c,name='ufo_geosaod_tlad_delete_f90')
implicit none
integer(c_int), intent(inout) :: c_key_self
    
type(ufo_geosaod_tlad), pointer :: self

call ufo_geosaod_tlad_registry%get(c_key_self, self)
call self%delete()
call ufo_geosaod_tlad_registry%remove(c_key_self)

end subroutine ufo_geosaod_tlad_delete_c

! ------------------------------------------------------------------------------

subroutine ufo_geosaod_tlad_settraj_c(c_key_self, c_key_geovals, c_obsspace) bind(c,name='ufo_geosaod_tlad_settraj_f90')

implicit none
integer(c_int),     intent(in) :: c_key_self
integer(c_int),     intent(in) :: c_key_geovals
type(c_ptr), value, intent(in) :: c_obsspace

type(ufo_geosaod_tlad), pointer :: self
type(ufo_geovals),      pointer :: geovals

call ufo_geosaod_tlad_registry%get(c_key_self, self)
call ufo_geovals_registry%get(c_key_geovals, geovals)

call self%settraj(geovals, c_obsspace)

end subroutine ufo_geosaod_tlad_settraj_c

! ------------------------------------------------------------------------------

subroutine ufo_geosaod_simobs_tl_c(c_key_self, c_key_geovals, c_obsspace, c_nvars, c_nlocs, c_hofx) bind(c,name='ufo_geosaod_simobs_tl_f90')

implicit none
integer(c_int),     intent(in) :: c_key_self
integer(c_int),     intent(in) :: c_key_geovals
type(c_ptr), value, intent(in) :: c_obsspace
integer(c_int),     intent(in) :: c_nvars, c_nlocs
real(c_double),     intent(inout) :: c_hofx(c_nvars, c_nlocs)

type(ufo_geosaod_tlad), pointer :: self
type(ufo_geovals),      pointer :: geovals

call ufo_geosaod_tlad_registry%get(c_key_self, self)
call ufo_geovals_registry%get(c_key_geovals, geovals)

call self%simobs_tl(geovals, c_obsspace, c_nvars, c_nlocs, c_hofx)

end subroutine ufo_geosaod_simobs_tl_c

! ------------------------------------------------------------------------------

subroutine ufo_geosaod_simobs_ad_c(c_key_self, c_key_geovals, c_obsspace, c_nvars, c_nlocs, c_hofx) bind(c,name='ufo_geosaod_simobs_ad_f90')

implicit none
integer(c_int),     intent(in) :: c_key_self
integer(c_int),     intent(in) :: c_key_geovals
type(c_ptr), value, intent(in) :: c_obsspace
integer(c_int),     intent(in) :: c_nvars, c_nlocs
real(c_double),     intent(in) :: c_hofx(c_nvars, c_nlocs)

type(ufo_geosaod_tlad), pointer :: self
type(ufo_geovals),      pointer :: geovals

call ufo_geosaod_tlad_registry%get(c_key_self, self)
call ufo_geovals_registry%get(c_key_geovals, geovals)

call self%simobs_ad(geovals, c_obsspace, c_nvars, c_nlocs, c_hofx)

end subroutine ufo_geosaod_simobs_ad_c

! ------------------------------------------------------------------------------

end module ufo_geosaod_tlad_mod_c
