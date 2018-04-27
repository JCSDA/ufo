! (C) Copyright 2017 UCAR
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 

!> Fortran module to handle ice concentration observations

module ufo_seaicefrac_tlad_mod_c
  
  use iso_c_binding
  use config_mod
  use ufo_obs_vectors,   only: obs_vector, ufo_obs_vect_registry
  use ufo_geovals_mod,   only: ufo_geovals
  use ufo_geovals_mod_c, only: ufo_geovals_registry
  use ufo_obs_seaicefrac_mod,   only: ufo_obs_seaicefrac
  use ufo_obs_seaicefrac_mod_c, only: ufo_obs_seaicefrac_registry 
  use ufo_seaicefrac_tlad_mod 
  implicit none
  private
  
#define LISTED_TYPE ufo_seaicefrac_tlad
  
  !> Linked list interface - defines registry_t type
#include "../../linkedList_i.f"
  
  !> Global registry
  type(registry_t) :: ufo_seaicefrac_tlad_registry
  
  ! ------------------------------------------------------------------------------
contains
  ! ------------------------------------------------------------------------------
  !> Linked list implementation
#include "../../linkedList_c.f"
  
! ------------------------------------------------------------------------------
  
subroutine ufo_seaicefrac_tlad_setup_c(c_key_self, c_conf) bind(c,name='ufo_seaicefrac_tlad_setup_f90')
implicit none
integer(c_int), intent(inout) :: c_key_self
type(c_ptr), intent(in)    :: c_conf
    
type(ufo_seaicefrac_tlad), pointer :: self

call ufo_seaicefrac_tlad_registry%init()
call ufo_seaicefrac_tlad_registry%add(c_key_self)
call ufo_seaicefrac_tlad_registry%get(c_key_self, self)
    
end subroutine ufo_seaicefrac_tlad_setup_c
  
! ------------------------------------------------------------------------------
  
subroutine ufo_seaicefrac_tlad_delete_c(c_key_self) bind(c,name='ufo_seaicefrac_tlad_delete_f90')
implicit none
integer(c_int), intent(inout) :: c_key_self
    
type(ufo_seaicefrac_tlad), pointer :: self

call ufo_seaicefrac_tlad_registry%get(c_key_self, self)
call ufo_seaicefrac_tlad_registry%remove(c_key_self)
    
end subroutine ufo_seaicefrac_tlad_delete_c
  
! ------------------------------------------------------------------------------

subroutine ufo_seaicefrac_tlad_settraj_c(c_key_self, c_key_geovals) bind(c,name='ufo_seaicefrac_tlad_settraj_f90')

implicit none
integer(c_int), intent(in) :: c_key_self
integer(c_int), intent(in) :: c_key_geovals

type(ufo_seaicefrac_tlad), pointer :: self
type(ufo_geovals),    pointer :: geovals

character(len=*), parameter :: myname_="ufo_seaicefrac_tlad_settraj_c"

call ufo_seaicefrac_tlad_registry%get(c_key_self, self)
call ufo_geovals_registry%get(c_key_geovals,geovals)

call ufo_seaicefrac_tlad_settraj(self, geovals)

end subroutine ufo_seaicefrac_tlad_settraj_c

! ------------------------------------------------------------------------------

subroutine ufo_seaicefrac_tlad_eqv_tl_c(c_key_self, c_key_geovals, c_key_hofx) bind(c,name='ufo_seaicefrac_tlad_eqv_tl_f90')

implicit none
integer(c_int), intent(in) :: c_key_self
integer(c_int), intent(in) :: c_key_geovals
integer(c_int), intent(in) :: c_key_hofx

type(ufo_seaicefrac_tlad), pointer :: self
type(ufo_geovals),    pointer :: geovals
type(obs_vector),     pointer :: hofx

character(len=*), parameter :: myname_="ufo_seaicefrac_tlad_eqv_tl_c"

call ufo_seaicefrac_tlad_registry%get(c_key_self, self)
call ufo_geovals_registry%get(c_key_geovals,geovals)
call ufo_obs_vect_registry%get(c_key_hofx,hofx)

call ufo_seaicefrac_tlad_eqv_tl(self, geovals, hofx)

end subroutine ufo_seaicefrac_tlad_eqv_tl_c

! ------------------------------------------------------------------------------

subroutine ufo_seaicefrac_tlad_eqv_ad_c(c_key_self, c_key_geovals, c_key_hofx) bind(c,name='ufo_seaicefrac_tlad_eqv_ad_f90')

implicit none
integer(c_int), intent(in) :: c_key_self
integer(c_int), intent(in) :: c_key_geovals
integer(c_int), intent(in) :: c_key_hofx

type(ufo_seaicefrac_tlad), pointer :: self
type(ufo_geovals),    pointer :: geovals
type(obs_vector),     pointer :: hofx


character(len=*), parameter :: myname_="ufo_seaicefrac_tlad_eqv_ad_c"

call ufo_seaicefrac_tlad_registry%get(c_key_self, self)
call ufo_geovals_registry%get(c_key_geovals,geovals)
call ufo_obs_vect_registry%get(c_key_hofx,hofx)

call ufo_seaicefrac_tlad_eqv_ad(self, geovals, hofx)

end subroutine ufo_seaicefrac_tlad_eqv_ad_c

! ------------------------------------------------------------------------------

end module ufo_seaicefrac_tlad_mod_c
