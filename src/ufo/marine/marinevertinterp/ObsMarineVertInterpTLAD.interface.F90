! (C) Copyright 2017-2019 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

!> Fortran marinevertinterp module for functions on the interface between C++ and Fortran
!  to handle tl/ad observation operators

module ufo_marinevertinterp_tlad_mod_c

  use iso_c_binding
  use config_mod
  use ufo_marinevertinterp_tlad_mod
  use string_f_c_mod    
  implicit none
  private

#define LISTED_TYPE ufo_marinevertinterp_tlad

  !> Linked list interface - defines registry_t type
#include "../../linkedList_i.f"

  !> Global registry
  type(registry_t) :: ufo_marinevertinterp_tlad_registry

contains

  ! ------------------------------------------------------------------------------
  !> Linked list implementation
#include "../../linkedList_c.f"

! ------------------------------------------------------------------------------

subroutine ufo_marinevertinterp_tlad_setup_c(c_key_self, c_conf, csin, csout, c_str_size) bind(c,name='ufo_marinevertinterp_tlad_setup_f90')
implicit none
integer(c_int), intent(inout) :: c_key_self
type(c_ptr), intent(in)    :: c_conf
    
type(ufo_marinevertinterp_tlad), pointer :: self
integer(c_int), intent(in) :: c_str_size
character(kind=c_char,len=1),intent(inout) :: csin(c_str_size+1),csout(c_str_size+1)


call ufo_marinevertinterp_tlad_registry%setup(c_key_self, self)

call self%setup(c_conf)

call f_c_string_vector(self%varout, csout)
call f_c_string_vector(self%varin, csin) 

end subroutine ufo_marinevertinterp_tlad_setup_c

! ------------------------------------------------------------------------------

subroutine ufo_marinevertinterp_tlad_delete_c(c_key_self) bind(c,name='ufo_marinevertinterp_tlad_delete_f90')
implicit none
integer(c_int), intent(inout) :: c_key_self
    
type(ufo_marinevertinterp_tlad), pointer :: self

call ufo_marinevertinterp_tlad_registry%get(c_key_self, self)
call self%opr_delete()
call ufo_marinevertinterp_tlad_registry%remove(c_key_self)

end subroutine ufo_marinevertinterp_tlad_delete_c

! ------------------------------------------------------------------------------

subroutine ufo_marinevertinterp_tlad_settraj_c(c_key_self, c_key_geovals, c_obsspace) bind(c,name='ufo_marinevertinterp_tlad_settraj_f90')

implicit none
integer(c_int),     intent(in) :: c_key_self
integer(c_int),     intent(in) :: c_key_geovals
type(c_ptr), value, intent(in) :: c_obsspace

type(ufo_marinevertinterp_tlad), pointer :: self

call ufo_marinevertinterp_tlad_registry%get(c_key_self, self)
call self%opr_settraj(c_key_geovals, c_obsspace)

end subroutine ufo_marinevertinterp_tlad_settraj_c

! ------------------------------------------------------------------------------

subroutine ufo_marinevertinterp_simobs_tl_c(c_key_self, c_key_geovals, c_obsspace, c_nobs, c_hofx) bind(c,name='ufo_marinevertinterp_simobs_tl_f90')

implicit none
integer(c_int), intent(in) :: c_key_self
integer(c_int), intent(in) :: c_key_geovals
type(c_ptr), value, intent(in) :: c_obsspace
integer(c_int), intent(in) :: c_nobs
real(c_double), intent(inout) :: c_hofx(c_nobs)

type(ufo_marinevertinterp_tlad), pointer :: self

call ufo_marinevertinterp_tlad_registry%get(c_key_self, self)
call self%opr_simobs_tl(c_key_geovals, c_obsspace, c_hofx)

end subroutine ufo_marinevertinterp_simobs_tl_c

! ------------------------------------------------------------------------------

subroutine ufo_marinevertinterp_simobs_ad_c(c_key_self, c_key_geovals, c_obsspace, c_nobs, c_hofx) bind(c,name='ufo_marinevertinterp_simobs_ad_f90')

implicit none
integer(c_int), intent(in) :: c_key_self
integer(c_int), intent(in) :: c_key_geovals
type(c_ptr), value, intent(in) :: c_obsspace
integer(c_int), intent(in) :: c_nobs
real(c_double), intent(in) :: c_hofx(c_nobs)

type(ufo_marinevertinterp_tlad), pointer :: self

call ufo_marinevertinterp_tlad_registry%get(c_key_self, self)
call self%opr_simobs_ad(c_key_geovals, c_obsspace, c_hofx)

end subroutine ufo_marinevertinterp_simobs_ad_c

! ------------------------------------------------------------------------------


end module ufo_marinevertinterp_tlad_mod_c
