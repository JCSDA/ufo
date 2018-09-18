! (C) Copyright 2017-2018 UCAR
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 

!> Fortran module to handle steric height observations

module ufo_stericheight_mod_c
  
  use iso_c_binding
  use config_mod
  use ioda_obs_vectors,   only: obs_vector, ioda_obs_vect_registry
  use ufo_geovals_mod,   only: ufo_geovals, ufo_geovals_setup
  use ufo_geovals_mod_c, only: ufo_geovals_registry
  use ioda_obs_stericheight_mod,   only: ioda_obs_stericheight
  use ioda_obs_stericheight_mod_c, only: ioda_obs_stericheight_registry 
  use ufo_stericheight_mod
  use ufo_vars_mod
  implicit none
  private
  
#define LISTED_TYPE ufo_stericheight
  
  !> Linked list interface - defines registry_t type
#include "../../linkedList_i.f"
  
  !> Global registry
  type(registry_t) :: ufo_stericheight_registry
  
  ! ------------------------------------------------------------------------------
contains
  ! ------------------------------------------------------------------------------
  !> Linked list implementation
#include "../../linkedList_c.f"
  
! ------------------------------------------------------------------------------
  
subroutine ufo_stericheight_setup_c(c_key_self, c_conf) bind(c,name='ufo_stericheight_setup_f90')
implicit none
integer(c_int), intent(inout) :: c_key_self
type(c_ptr), intent(in)    :: c_conf
    
type(ufo_stericheight), pointer :: self

call ufo_stericheight_registry%init()
call ufo_stericheight_registry%add(c_key_self)
call ufo_stericheight_registry%get(c_key_self, self)
    
end subroutine ufo_stericheight_setup_c
  
! ------------------------------------------------------------------------------
  
subroutine ufo_stericheight_delete_c(c_key_self) bind(c,name='ufo_stericheight_delete_f90')
implicit none
integer(c_int), intent(inout) :: c_key_self
    
type(ufo_stericheight), pointer :: self

call ufo_stericheight_registry%get(c_key_self, self)
call ufo_stericheight_registry%remove(c_key_self)
    
end subroutine ufo_stericheight_delete_c
  
! ------------------------------------------------------------------------------

subroutine ufo_stericheight_simobs_c(c_key_self, c_key_geovals, c_key_obsspace, c_key_hofx, c_bias) bind(c,name='ufo_stericheight_simobs_f90')

implicit none
integer(c_int), intent(in) :: c_key_self
integer(c_int), intent(in) :: c_key_geovals
integer(c_int), intent(in) :: c_key_hofx
integer(c_int), intent(in) :: c_key_obsspace
integer(c_int), intent(in) :: c_bias

type(ufo_stericheight), pointer :: self
type(ufo_geovals),    pointer :: geovals
type(obs_vector),     pointer :: hofx

character(len=*), parameter :: myname_="ufo_stericheight_simobs_c"

print *,myname_

call ufo_stericheight_registry%get(c_key_self, self)
call ufo_geovals_registry%get(c_key_geovals,geovals)
call ioda_obs_vect_registry%get(c_key_hofx,hofx)

call ufo_stericheight_simobs(self, geovals, hofx)

end subroutine ufo_stericheight_simobs_c
  
end module ufo_stericheight_mod_c
