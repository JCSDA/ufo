! (C) Copyright 2017-2018 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

!> Fortran module to handle radiosonde observations

module ufo_radiosonde_mod_c

  use iso_c_binding
  use config_mod
  use ufo_radiosonde_mod
  implicit none

  integer, parameter :: max_string=800

  private

#define LISTED_TYPE ufo_radiosonde

  !> Linked list interface - defines registry_t type
#include "../../linkedList_i.f"

  !> Global registry
  type(registry_t) :: ufo_radiosonde_registry

  ! ------------------------------------------------------------------------------

contains
  ! ------------------------------------------------------------------------------
  !> Linked list implementation
#include "../../linkedList_c.f"

! ------------------------------------------------------------------------------

subroutine ufo_radiosonde_setup_c(c_key_self, c_conf) bind(c,name='ufo_radiosonde_setup_f90')
implicit none
integer(c_int), intent(inout) :: c_key_self
type(c_ptr), intent(in)    :: c_conf

type(ufo_radiosonde), pointer :: self

call ufo_radiosonde_registry%setup(c_key_self, self)

if (config_element_exists(c_conf,"variables")) then
     self%nvars = size(config_get_string_vector(c_conf, max_string, "variables"))
     if (allocated(self%varout)) deallocate(self%varout)
     allocate(self%varout(self%nvars))
     self%varout = config_get_string_vector(c_conf, max_string, "variables")
endif

end subroutine ufo_radiosonde_setup_c

! ------------------------------------------------------------------------------

subroutine ufo_radiosonde_delete_c(c_key_self) bind(c,name='ufo_radiosonde_delete_f90')
implicit none
integer(c_int), intent(inout) :: c_key_self

type(ufo_radiosonde), pointer :: self

call ufo_radiosonde_registry%delete(c_key_self, self)

end subroutine ufo_radiosonde_delete_c

! ------------------------------------------------------------------------------

subroutine ufo_radiosonde_simobs_c(c_key_self, c_key_geovals, c_obsspace, c_nobs, c_hofx, c_bias) bind(c,name='ufo_radiosonde_simobs_f90')

implicit none
integer(c_int), intent(in) :: c_key_self
integer(c_int), intent(in) :: c_key_geovals
type(c_ptr), value, intent(in) :: c_obsspace
integer(c_int), intent(in) :: c_nobs
real(c_double), intent(inout) :: c_hofx(c_nobs)
integer(c_int), intent(in) :: c_bias

type(ufo_radiosonde), pointer :: self

character(len=*), parameter :: myname_="ufo_radiosonde_simobs_c"

call ufo_radiosonde_registry%get(c_key_self, self)
call self%opr_simobs(c_key_geovals, c_obsspace, c_hofx)

end subroutine ufo_radiosonde_simobs_c

! ------------------------------------------------------------------------------

end module ufo_radiosonde_mod_c
