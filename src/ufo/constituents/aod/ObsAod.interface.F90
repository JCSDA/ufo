! (C) Copyright 2017 UCAR
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 

!> Fortran module to handle aod observations

module ufo_aod_mod_c
  
  use iso_c_binding
  use config_mod
  use ufo_obs_vectors,   only: obs_vector, ufo_obs_vect_registry
  use ufo_geovals_mod,   only: ufo_geovals
  use ufo_geovals_mod_c, only: ufo_geovals_registry
  use ufo_obs_aod_mod,   only: ufo_obs_aod
  use ufo_obs_aod_mod_c, only: ufo_obs_aod_registry 
  use ufo_aod_mod 
  implicit none
  private
  
#define LISTED_TYPE ufo_aod
  
  !> Linked list interface - defines registry_t type
#include "../../linkedList_i.f"
  
  !> Global registry
  type(registry_t) :: ufo_aod_registry
  
  ! ------------------------------------------------------------------------------
contains
  ! ------------------------------------------------------------------------------
  !> Linked list implementation
#include "../../linkedList_c.f"
  
! ------------------------------------------------------------------------------
  
subroutine ufo_aod_setup_c(c_key_self, c_conf) bind(c,name='ufo_aod_setup_f90')
implicit none
integer(c_int), intent(inout) :: c_key_self
type(c_ptr), intent(in)    :: c_conf
    
type(ufo_aod), pointer :: self

call ufo_aod_registry%init()
call ufo_aod_registry%add(c_key_self)
call ufo_aod_registry%get(c_key_self, self)
    
end subroutine ufo_aod_setup_c
  
! ------------------------------------------------------------------------------
  
subroutine ufo_aod_delete_c(c_key_self) bind(c,name='ufo_aod_delete_f90')
implicit none
integer(c_int), intent(inout) :: c_key_self
    
type(ufo_aod), pointer :: self

call ufo_aod_registry%get(c_key_self, self)
call ufo_aod_registry%remove(c_key_self)
    
end subroutine ufo_aod_delete_c
  
! ------------------------------------------------------------------------------

subroutine ufo_aod_eqv_c(c_key_self, c_key_geovals, c_key_obsspace, c_key_hofx, c_bias) bind(c,name='ufo_aod_eqv_f90')

implicit none
integer(c_int), intent(in) :: c_key_self
integer(c_int), intent(in) :: c_key_geovals
integer(c_int), intent(in) :: c_key_hofx
integer(c_int), intent(in) :: c_key_obsspace
integer(c_int), intent(in) :: c_bias

type(ufo_aod),     pointer :: self
type(ufo_geovals),        pointer :: geovals
type(obs_vector),         pointer :: hofx
type(ufo_obs_aod), pointer :: obss

character(len=*), parameter :: myname_="ufo_aod_eqv_c"

call ufo_aod_registry%get(c_key_self, self)
call ufo_geovals_registry%get(c_key_geovals,geovals)
call ufo_obs_vect_registry%get(c_key_hofx,hofx)
call ufo_obs_aod_registry%get(c_key_obsspace,obss)

call ufo_aod_eqv(self, geovals, hofx, obss)

end subroutine ufo_aod_eqv_c

! ------------------------------------------------------------------------------

subroutine ufo_aod_settraj_c(c_key_self, c_key_geovals) bind(c,name='ufo_aod_settraj_f90')

implicit none
integer(c_int), intent(in) :: c_key_self
integer(c_int), intent(in) :: c_key_geovals

type(ufo_aod), pointer :: self
type(ufo_geovals),    pointer :: geovals

character(len=*), parameter :: myname_="ufo_aod_settraj_c"

call ufo_aod_registry%get(c_key_self, self)
call ufo_geovals_registry%get(c_key_geovals,geovals)

!call ufo_aod_settraj(self, geovals)

end subroutine ufo_aod_settraj_c

! ------------------------------------------------------------------------------

subroutine ufo_aod_eqv_tl_c(c_key_self, c_key_geovals, c_key_obsspace, c_key_hofx) bind(c,name='ufo_aod_eqv_tl_f90')

implicit none
integer(c_int), intent(in) :: c_key_self
integer(c_int), intent(in) :: c_key_geovals
integer(c_int), intent(in) :: c_key_hofx
integer(c_int), intent(in) :: c_key_obsspace

type(ufo_aod), pointer :: self
type(ufo_geovals),    pointer :: geovals
type(obs_vector),     pointer :: hofx
type(ufo_obs_aod), pointer :: obss

character(len=*), parameter :: myname_="ufo_aod_eqv_tl_c"

call ufo_aod_registry%get(c_key_self, self)
call ufo_geovals_registry%get(c_key_geovals,geovals)
call ufo_obs_vect_registry%get(c_key_hofx,hofx)
call ufo_obs_aod_registry%get(c_key_obsspace,obss)

!call ufo_aod_eqv_tl(self, geovals, hofx, obss)

end subroutine ufo_aod_eqv_tl_c

! ------------------------------------------------------------------------------

subroutine ufo_aod_eqv_ad_c(c_key_self, c_key_geovals, c_key_obsspace, c_key_hofx) bind(c,name='ufo_aod_eqv_ad_f90')

implicit none
integer(c_int), intent(in) :: c_key_self
integer(c_int), intent(in) :: c_key_geovals
integer(c_int), intent(in) :: c_key_hofx
integer(c_int), intent(in) :: c_key_obsspace

type(ufo_aod), pointer :: self
type(ufo_geovals),    pointer :: geovals
type(obs_vector),     pointer :: hofx
type(ufo_obs_aod), pointer :: obss

character(len=*), parameter :: myname_="ufo_aod_eqv_ad_c"

call ufo_aod_registry%get(c_key_self, self)
call ufo_geovals_registry%get(c_key_geovals,geovals)
call ufo_obs_vect_registry%get(c_key_hofx,hofx)
call ufo_obs_aod_registry%get(c_key_obsspace,obss)

!call ufo_aod_eqv_ad(self, geovals, hofx, obss)

end subroutine ufo_aod_eqv_ad_c

! ------------------------------------------------------------------------------

subroutine ufo_obs_get(c_key_self, lcol, c_col, c_key_ovec) bind(c,name='ufo_obsdb_aod_get_f90')  
use  ufo_obs_aod_mod
implicit none
integer(c_int), intent(in) :: c_key_self
integer(c_int), intent(in) :: lcol
character(kind=c_char,len=1), intent(in) :: c_col(lcol+1)
integer(c_int), intent(in) :: c_key_ovec

type(ufo_obs_aod), pointer :: self
type(obs_vector), pointer :: ovec
character(len=lcol) :: col

call ufo_obs_aod_registry%get(c_key_self, self)
call ufo_obs_vect_registry%get(c_key_ovec,ovec)


end subroutine ufo_obs_get
  
end module ufo_aod_mod_c
