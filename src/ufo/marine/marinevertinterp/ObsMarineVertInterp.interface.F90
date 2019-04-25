! (C) Copyright 2017-2019 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

!> Fortran marinevertinterp module for functions on the interface between C++ and Fortran
!  to handle observation operators

module ufo_marinevertinterp_mod_c

  use iso_c_binding
  use config_mod
  use ufo_marinevertinterp_mod
  use string_f_c_mod
  use ufo_geovals_mod
  use ufo_geovals_mod_c,   only: ufo_geovals_registry
  
  implicit none
  private

  ! ------------------------------------------------------------------------------
#define LISTED_TYPE ufo_marinevertinterp

  !> Linked list interface - defines registry_t type
#include "../../linkedList_i.f"

  !> Global registry
  type(registry_t) :: ufo_marinevertinterp_registry

  ! ------------------------------------------------------------------------------

contains

  ! ------------------------------------------------------------------------------
  !> Linked list implementation
#include "../../linkedList_c.f"

! ------------------------------------------------------------------------------

subroutine ufo_marinevertinterp_setup_c(c_key_self, c_conf, csin, csout, c_str_size) bind(c,name='ufo_marinevertinterp_setup_f90')
implicit none
integer(c_int), intent(inout) :: c_key_self
type(c_ptr),    intent(in)    :: c_conf
integer(c_int), intent(in) :: c_str_size
character(kind=c_char,len=1),intent(inout) :: csin(c_str_size+1),csout(c_str_size+1)

type(ufo_marinevertinterp), pointer :: self

call ufo_marinevertinterp_registry%setup(c_key_self, self)

call self%setup(c_conf)

call f_c_string_vector(self%varout, csout)
call f_c_string_vector(self%varin, csin) 

end subroutine ufo_marinevertinterp_setup_c

! ------------------------------------------------------------------------------

subroutine ufo_marinevertinterp_delete_c(c_key_self) bind(c,name='ufo_marinevertinterp_delete_f90')
implicit none
integer(c_int), intent(inout) :: c_key_self
    
type(ufo_marinevertinterp), pointer :: self

call ufo_marinevertinterp_registry%get(c_key_self, self)

call self%delete()

call ufo_marinevertinterp_registry%remove(c_key_self)

end subroutine ufo_marinevertinterp_delete_c

! ------------------------------------------------------------------------------

subroutine ufo_marinevertinterp_simobs_c(c_key_self, c_key_geovals, c_obsspace, c_nobs, c_hofx, c_bias) bind(c,name='ufo_marinevertinterp_simobs_f90')

implicit none
integer(c_int), intent(in) :: c_key_self
integer(c_int), intent(in) :: c_key_geovals
type(c_ptr), value, intent(in) :: c_obsspace
integer(c_int), intent(in) :: c_nobs
real(c_double), intent(inout) :: c_hofx(c_nobs)
integer(c_int), intent(in) :: c_bias

type(ufo_marinevertinterp), pointer :: self
type(ufo_geovals), pointer :: geovals

call ufo_marinevertinterp_registry%get(c_key_self, self)
call ufo_geovals_registry%get(c_key_geovals,geovals)

call self%simobs(geovals, c_hofx, c_obsspace)

end subroutine ufo_marinevertinterp_simobs_c

! ------------------------------------------------------------------------------

end module ufo_marinevertinterp_mod_c
