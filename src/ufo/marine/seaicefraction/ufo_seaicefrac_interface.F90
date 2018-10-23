! (C) Copyright 2017-2018 UCAR
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 

!> Fortran module to handle ice concentration observations

module ufo_seaicefrac_mod_c
  
  use iso_c_binding
  use config_mod
  use ioda_obs_vectors,   only: obs_vector, ioda_obs_vect_registry
  use ufo_geovals_mod,   only: ufo_geovals
  use ufo_geovals_mod_c, only: ufo_geovals_registry
  use ioda_obs_seaicefrac_mod,   only: ioda_obs_seaicefrac
  use ioda_obs_seaicefrac_mod_c, only: ioda_obs_seaicefrac_registry 
  use ufo_seaicefrac_mod 
  implicit none
  private
  
#define LISTED_TYPE ufo_seaicefrac
  
  !> Linked list interface - defines registry_t type
#include "../../linkedList_i.f"
  
  !> Global registry
  type(registry_t) :: ufo_seaicefrac_registry
  
  ! ------------------------------------------------------------------------------
contains
  ! ------------------------------------------------------------------------------
  !> Linked list implementation
#include "../../linkedList_c.f"
  
! ------------------------------------------------------------------------------
  
subroutine ufo_seaicefrac_setup_c(c_key_self, c_conf) bind(c,name='ufo_seaicefrac_setup_f90')
implicit none
integer(c_int), intent(inout) :: c_key_self
type(c_ptr), intent(in)    :: c_conf
    
type(ufo_seaicefrac), pointer :: self

call ufo_seaicefrac_registry%init()
call ufo_seaicefrac_registry%add(c_key_self)
call ufo_seaicefrac_registry%get(c_key_self, self)
    
end subroutine ufo_seaicefrac_setup_c
  
! ------------------------------------------------------------------------------
  
subroutine ufo_seaicefrac_delete_c(c_key_self) bind(c,name='ufo_seaicefrac_delete_f90')
implicit none
integer(c_int), intent(inout) :: c_key_self
    
type(ufo_seaicefrac), pointer :: self

call ufo_seaicefrac_registry%get(c_key_self, self)
call ufo_seaicefrac_registry%remove(c_key_self)
    
end subroutine ufo_seaicefrac_delete_c
  
! ------------------------------------------------------------------------------

subroutine ufo_seaicefrac_simobs_c(c_key_self, c_key_geovals, c_obsspace, c_key_hofx, c_bias) bind(c,name='ufo_seaicefrac_simobs_f90')

implicit none
integer(c_int),     intent(in) :: c_key_self
integer(c_int),     intent(in) :: c_key_geovals
integer(c_int),     intent(in) :: c_key_hofx
type(c_ptr), value, intent(in) :: c_obsspace
integer(c_int),     intent(in) :: c_bias

type(ufo_seaicefrac), pointer :: self
type(ufo_geovals),    pointer :: geovals
type(obs_vector),     pointer :: hofx

character(len=*), parameter :: myname_="ufo_seaicefrac_simobs_c"

call ufo_seaicefrac_registry%get(c_key_self, self)
call ufo_geovals_registry%get(c_key_geovals,geovals)
call ioda_obs_vect_registry%get(c_key_hofx,hofx)

call ufo_seaicefrac_simobs(self, geovals, hofx)

end subroutine ufo_seaicefrac_simobs_c
  
end module ufo_seaicefrac_mod_c
