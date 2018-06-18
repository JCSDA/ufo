! (C) Copyright 2017 UCAR
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 

!> Fortran module to handle radiance observations

module ufo_radiance_mod_c
  
  use iso_c_binding
  use config_mod
  use ufo_radiance_mod 
  implicit none
  private
  
  ! ------------------------------------------------------------------------------
#define LISTED_TYPE ufo_radiance

  !> Linked list interface - defines registry_t type
#include "../../linkedList_i.f"

  !> Global registry
  type(registry_t) :: ufo_radiance_registry

  ! ------------------------------------------------------------------------------

contains

  ! ------------------------------------------------------------------------------
  !> Linked list implementation
#include "../../linkedList_c.f"

! ------------------------------------------------------------------------------
  
subroutine ufo_radiance_setup_c(c_key_self, c_conf) bind(c,name='ufo_radiance_setup_f90')
implicit none
integer(c_int), intent(inout) :: c_key_self
type(c_ptr), intent(in)    :: c_conf
    
type(ufo_radiance), pointer :: self

call ufo_radiance_registry%setup(c_key_self, self)
    
end subroutine ufo_radiance_setup_c
  
! ------------------------------------------------------------------------------
  
subroutine ufo_radiance_delete_c(c_key_self) bind(c,name='ufo_radiance_delete_f90')
implicit none
integer(c_int), intent(inout) :: c_key_self
    
type(ufo_radiance), pointer :: self

call ufo_radiance_registry%delete(c_key_self, self)
    
end subroutine ufo_radiance_delete_c
  
! ------------------------------------------------------------------------------

subroutine ufo_radiance_eqv_c(c_key_self, c_key_geovals, c_key_obsspace, c_key_hofx, c_bias) bind(c,name='ufo_radiance_eqv_f90')

implicit none
integer(c_int), intent(in) :: c_key_self
integer(c_int), intent(in) :: c_key_geovals
integer(c_int), intent(in) :: c_key_hofx
integer(c_int), intent(in) :: c_key_obsspace
integer(c_int), intent(in) :: c_bias

type(ufo_radiance), pointer :: self

character(len=*), parameter :: myname_="ufo_radiance_eqv_c"

call ufo_radiance_registry%get(c_key_self, self)
call self%opr_eqv(c_key_geovals, c_key_obsspace, c_key_hofx)

end subroutine ufo_radiance_eqv_c
  
end module ufo_radiance_mod_c
