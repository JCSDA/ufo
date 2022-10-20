! (C) Copyright 2017 UCAR
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 

!> Fortran module to handle ground-based gnss observations ROPP 1d operator

module ufo_gnssgb_refropp1d_mod_c
  
  use fckit_configuration_module, only: fckit_configuration 
  use iso_c_binding
  use ufo_gnssgb_refropp1d_mod

  implicit none
  private
  
#define LISTED_TYPE ufo_gnssgb_RefROPP1D
  
  !> Linked list interface - defines registry_t type
#include "oops/util/linkedList_i.f"
  
  !> Global registry
  type(registry_t) :: ufo_gnssgb_RefROPP1D_registry
  
  ! ------------------------------------------------------------------------------
contains
  ! ------------------------------------------------------------------------------
  !> Linked list implementation
#include "oops/util/linkedList_c.f"
  
! ------------------------------------------------------------------------------
  
subroutine ufo_gnssgb_refropp1d_setup_c(c_key_self, c_conf) bind(c,name='ufo_gnssgb_refropp1d_setup_f90')
implicit none
integer(c_int), intent(inout)  :: c_key_self
type(c_ptr), value, intent(in) :: c_conf
    
type(ufo_gnssgb_RefROPP1D), pointer :: self
type(fckit_configuration) :: f_conf

call ufo_gnssgb_RefROPP1D_registry%setup(c_key_self, self)
f_conf = fckit_configuration(c_conf)

end subroutine ufo_gnssgb_RefROPP1D_setup_c
  
! ------------------------------------------------------------------------------
  
subroutine ufo_gnssgb_refropp1d_delete_c(c_key_self) bind(c,name='ufo_gnssgb_refropp1d_delete_f90')
implicit none
integer(c_int), intent(inout) :: c_key_self

call ufo_gnssgb_RefROPP1D_registry%remove(c_key_self)

end subroutine ufo_gnssgb_refropp1d_delete_c
  
! ------------------------------------------------------------------------------

subroutine ufo_gnssgb_refropp1d_simobs_c(c_key_self, c_key_geovals, c_obsspace, c_nobs, c_hofx) &
    bind(c,name='ufo_gnssgb_refropp1d_simobs_f90')

implicit none
integer(c_int), intent(in) :: c_key_self
integer(c_int), intent(in) :: c_key_geovals
type(c_ptr), value, intent(in) :: c_obsspace
integer(c_int), intent(in) :: c_nobs
real(c_double), intent(inout) :: c_hofx(c_nobs)

type(ufo_gnssgb_RefROPP1D),     pointer :: self

character(len=*), parameter :: myname_="ufo_gnssgb_refropp1d_simobs_c"
call ufo_gnssgb_RefROPP1D_registry%get(c_key_self, self)
call self%opr_simobs(c_key_geovals, c_obsspace, c_hofx)

end subroutine ufo_gnssgb_refropp1d_simobs_c

! ------------------------------------------------------------------------------

end module ufo_gnssgb_refropp1d_mod_c
