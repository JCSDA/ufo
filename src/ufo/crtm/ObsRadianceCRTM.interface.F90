! (C) Copyright 2017-2018 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

!> Fortran module to handle radiancecrtm observations

module ufo_radiancecrtm_mod_c

  use fckit_configuration_module, only: fckit_configuration
  use iso_c_binding
  use ufo_radiancecrtm_mod
  use string_f_c_mod
  use ufo_geovals_mod
  use ufo_geovals_mod_c,   only: ufo_geovals_registry

  implicit none
  private

  ! ------------------------------------------------------------------------------
#define LISTED_TYPE ufo_radiancecrtm

  !> Linked list interface - defines registry_t type
#include "../linkedList_i.f"

  !> Global registry
  type(registry_t) :: ufo_radiancecrtm_registry

  ! ------------------------------------------------------------------------------

contains

  ! ------------------------------------------------------------------------------
  !> Linked list implementation
#include "../linkedList_c.f"

! ------------------------------------------------------------------------------

subroutine ufo_radiancecrtm_setup_c(c_key_self, c_conf, c_nchan, c_channels, c_varlist) &
                                    bind(c,name='ufo_radiancecrtm_setup_f90')
implicit none
integer(c_int), intent(inout) :: c_key_self
type(c_ptr),    intent(in)    :: c_conf
integer(c_int), intent(in) :: c_nchan
integer(c_int), intent(in) :: c_channels(c_nchan)
type(c_ptr), intent(in), value :: c_varlist

type(ufo_radiancecrtm), pointer :: self
type(fckit_configuration) :: f_conf

call ufo_radiancecrtm_registry%setup(c_key_self, self)
f_conf = fckit_configuration(c_conf)

call self%setup(f_conf, c_channels)

!> Update C++ ObsOperator with input variable list
call f_c_push_string_varlist(c_varlist, self%varin)

end subroutine ufo_radiancecrtm_setup_c

! ------------------------------------------------------------------------------

subroutine ufo_radiancecrtm_delete_c(c_key_self) bind(c,name='ufo_radiancecrtm_delete_f90')
implicit none
integer(c_int), intent(inout) :: c_key_self

type(ufo_radiancecrtm), pointer :: self

call ufo_radiancecrtm_registry%get(c_key_self, self)

call self%delete()

call ufo_radiancecrtm_registry%remove(c_key_self)

end subroutine ufo_radiancecrtm_delete_c

! ------------------------------------------------------------------------------
subroutine ufo_radiancecrtm_simobs_c(c_key_self, c_key_geovals, c_obsspace, c_nvars, c_nlocs, &
           c_hofx, c_key_hofxdiags) bind(c,name='ufo_radiancecrtm_simobs_f90')
implicit none
integer(c_int), intent(in) :: c_key_self
integer(c_int), intent(in) :: c_key_geovals
type(c_ptr), value, intent(in) :: c_obsspace
integer(c_int), intent(in) :: c_nvars, c_nlocs
real(c_double), intent(inout) :: c_hofx(c_nvars, c_nlocs)
integer(c_int), intent(in) :: c_key_hofxdiags

type(ufo_radiancecrtm), pointer :: self
type(ufo_geovals),  pointer :: geovals
type(ufo_geovals),  pointer :: hofxdiags

character(len=*), parameter :: myname_="ufo_radiancecrtm_simobs_c"

call ufo_radiancecrtm_registry%get(c_key_self, self)

call ufo_geovals_registry%get(c_key_geovals,geovals)
call ufo_geovals_registry%get(c_key_hofxdiags,hofxdiags)

call self%simobs(geovals, c_obsspace, c_nvars, c_nlocs, c_hofx, hofxdiags)

end subroutine ufo_radiancecrtm_simobs_c

! ------------------------------------------------------------------------------

end module ufo_radiancecrtm_mod_c
