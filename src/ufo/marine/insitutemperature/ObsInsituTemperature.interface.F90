! (C) Copyright 2017-2018 UCAR
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 

!> Fortran module to handle temperature profile observations

module ufo_insitutemperature_mod_c

  use iso_c_binding
  use config_mod
  use ioda_obs_vectors,   only: obs_vector, ioda_obs_vect_registry
  use ufo_geovals_mod,   only: ufo_geovals
  use ufo_geovals_mod_c, only: ufo_geovals_registry
  use ioda_obs_insitutemperature_mod,   only: ioda_obs_insitutemperature
  use ioda_obs_insitutemperature_mod_c, only: ioda_obs_insitutemperature_registry 
  use ufo_insitutemperature_mod 
  implicit none
  private
  
#define LISTED_TYPE ufo_insitutemperature
  
  !> Linked list interface - defines registry_t type
#include "../../linkedList_i.f"
  
  !> Global registry
  type(registry_t) :: ufo_insitutemperature_registry
  
  ! ------------------------------------------------------------------------------
contains
  ! ------------------------------------------------------------------------------
  !> Linked list implementation
#include "../../linkedList_c.f"
  
! ------------------------------------------------------------------------------
  
subroutine ufo_insitutemperature_setup_c(c_key_self, c_conf) bind(c,name='ufo_insitutemperature_setup_f90')
implicit none
integer(c_int), intent(inout) :: c_key_self
type(c_ptr), intent(in)    :: c_conf
    
type(ufo_insitutemperature), pointer :: self

call ufo_insitutemperature_registry%init()
call ufo_insitutemperature_registry%add(c_key_self)
call ufo_insitutemperature_registry%get(c_key_self, self)
    
end subroutine ufo_insitutemperature_setup_c
  
! ------------------------------------------------------------------------------
  
subroutine ufo_insitutemperature_delete_c(c_key_self) bind(c,name='ufo_insitutemperature_delete_f90')
implicit none
integer(c_int), intent(inout) :: c_key_self
    
type(ufo_insitutemperature), pointer :: self

call ufo_insitutemperature_registry%get(c_key_self, self)
call ufo_insitutemperature_registry%remove(c_key_self)
    
end subroutine ufo_insitutemperature_delete_c
  
! ------------------------------------------------------------------------------

subroutine ufo_insitutemperature_eqv_c(c_key_self, c_key_geovals, c_key_obsspace, c_key_hofx, c_bias) bind(c,name='ufo_insitutemperature_eqv_f90')

implicit none
integer(c_int), intent(in) :: c_key_self
integer(c_int), intent(in) :: c_key_geovals
integer(c_int), intent(in) :: c_key_hofx
integer(c_int), intent(in) :: c_key_obsspace
integer(c_int), intent(in) :: c_bias

type(ufo_insitutemperature), pointer :: self
type(ufo_geovals),    pointer :: geovals
type(obs_vector),     pointer :: hofx
type(ioda_obs_insitutemperature), pointer :: obs_ti
character(len=*), parameter :: myname_="ufo_insitutemperature_eqv_c"

call ufo_insitutemperature_registry%get(c_key_self, self)
call ufo_geovals_registry%get(c_key_geovals,geovals)
call ioda_obs_vect_registry%get(c_key_hofx,hofx)
call ioda_obs_insitutemperature_registry%get(c_key_obsspace,obs_ti)

call ufo_insitutemperature_eqv(self, geovals, hofx, obs_ti)

end subroutine ufo_insitutemperature_eqv_c

end module ufo_insitutemperature_mod_c
