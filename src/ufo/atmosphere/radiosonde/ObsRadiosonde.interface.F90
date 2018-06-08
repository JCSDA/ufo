! (C) Copyright 2017 UCAR
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 

!> Fortran module to handle radiosonde observations

module ufo_radiosonde_mod_c
  
  use iso_c_binding
  use config_mod
  use ufo_radiosonde_mod 
  implicit none
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
    
end subroutine ufo_radiosonde_setup_c
  
! ------------------------------------------------------------------------------
  
subroutine ufo_radiosonde_delete_c(c_key_self) bind(c,name='ufo_radiosonde_delete_f90')
implicit none
integer(c_int), intent(inout) :: c_key_self
    
type(ufo_radiosonde), pointer :: self

call ufo_radiosonde_registry%delete(c_key_self, self)
    
end subroutine ufo_radiosonde_delete_c
  
! ------------------------------------------------------------------------------

subroutine ufo_radiosonde_t_eqv_c(c_key_self, c_key_geovals, c_key_obsspace, c_key_hofx, c_bias) bind(c,name='ufo_radiosonde_t_eqv_f90')

implicit none
integer(c_int), intent(in) :: c_key_self
integer(c_int), intent(in) :: c_key_geovals
integer(c_int), intent(in) :: c_key_hofx
integer(c_int), intent(in) :: c_key_obsspace
integer(c_int), intent(in) :: c_bias

type(ufo_radiosonde), pointer :: self

character(len=*), parameter :: myname_="ufo_radiosonde_t_eqv_c"

call ufo_radiosonde_registry%get(c_key_self, self)
call self%opr_eqv(c_key_geovals, c_key_obsspace, c_key_hofx)

end subroutine ufo_radiosonde_t_eqv_c

end module ufo_radiosonde_mod_c
