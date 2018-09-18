! (C) Copyright 2017-2018 UCAR
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 

!> Fortran module to handle steric height observations

module ufo_stericheight_tlad_mod_c
  
  use iso_c_binding
  use config_mod
  use ioda_obs_vectors,   only: obs_vector, ioda_obs_vect_registry
  use ufo_geovals_mod,   only: ufo_geovals, ufo_geovals_setup
  use ufo_geovals_mod_c, only: ufo_geovals_registry
  use ioda_obs_stericheight_mod,   only: ioda_obs_stericheight
  use ioda_obs_stericheight_mod_c, only: ioda_obs_stericheight_registry 
  use ufo_stericheight_tlad_mod
  use ufo_vars_mod
  implicit none
  private
  
#define LISTED_TYPE ufo_stericheight_tlad
  
  !> Linked list interface - defines registry_t type
#include "../../linkedList_i.f"
  
  !> Global registry
  type(registry_t) :: ufo_stericheight_tlad_registry
  
  ! ------------------------------------------------------------------------------
contains
  ! ------------------------------------------------------------------------------
  !> Linked list implementation
#include "../../linkedList_c.f"
  
! ------------------------------------------------------------------------------
  
subroutine ufo_stericheight_tlad_setup_c(c_key_self, c_conf) bind(c,name='ufo_stericheight_tlad_setup_f90')
implicit none
integer(c_int), intent(inout) :: c_key_self
type(c_ptr), intent(in)    :: c_conf
    
type(ufo_stericheight_tlad), pointer :: self

call ufo_stericheight_tlad_registry%init()
call ufo_stericheight_tlad_registry%add(c_key_self)
call ufo_stericheight_tlad_registry%get(c_key_self, self)
    
end subroutine ufo_stericheight_tlad_setup_c
  
! ------------------------------------------------------------------------------
  
subroutine ufo_stericheight_tlad_delete_c(c_key_self) bind(c,name='ufo_stericheight_tlad_delete_f90')
implicit none
integer(c_int), intent(inout) :: c_key_self
    
type(ufo_stericheight_tlad), pointer :: self

call ufo_stericheight_tlad_registry%get(c_key_self, self)
call ufo_stericheight_tlad_registry%remove(c_key_self)
    
end subroutine ufo_stericheight_tlad_delete_c
  
! ------------------------------------------------------------------------------

subroutine ufo_stericheight_tlad_gettraj(c_key_self, c_nobs, c_vars, c_key_traj) bind(c,name='ufo_stericheight_tlad_gettraj_f90')
use fckit_log_module, only : fckit_log

implicit none
integer(c_int), intent(in) :: c_key_self
integer(c_int), intent(in) :: c_nobs
type(c_ptr), intent(in)       :: c_vars
integer(c_int), intent(inout) :: c_key_traj

type(ufo_stericheight_tlad), pointer :: self
type(ufo_geovals), pointer :: traj
type(ufo_vars) :: cvars
integer nobs

call ufo_stericheight_tlad_registry%get(c_key_self, self)
call ufo_geovals_registry%get(c_key_traj,traj)
call ufo_vars_setup(cvars, c_vars)
!call qg_obsoper_registry%get(c_key_self, self)
!call qg_goms_registry%init()
!call qg_goms_registry%add(c_key_traj)
!call qg_goms_registry%get(c_key_traj,traj)
!call qg_vars_create(vars, c_vars)
!allocate(mobs(c_nobs))
!do jj=1,c_nobs
!  mobs(jj)=jj
!enddo
nobs=c_nobs
call ufo_geovals_setup(traj, cvars, nobs)
!deallocate(mobs)

end subroutine ufo_stericheight_tlad_gettraj

! ------------------------------------------------------------------------------

subroutine ufo_stericheight_tlad_settraj_c(c_key_self, c_key_geovals) bind(c,name='ufo_stericheight_tlad_settraj_f90')

implicit none
integer(c_int), intent(in) :: c_key_self
integer(c_int), intent(in) :: c_key_geovals

type(ufo_stericheight_tlad), pointer :: self
type(ufo_geovals),    pointer :: geovals

character(len=*), parameter :: myname_="ufo_stericheight_tlad_settraj_c"

print *, myname_

call ufo_stericheight_tlad_registry%get(c_key_self, self)
call ufo_geovals_registry%get(c_key_geovals,geovals)

call ufo_stericheight_tlad_settraj(self, geovals)

end subroutine ufo_stericheight_tlad_settraj_c

! ------------------------------------------------------------------------------

subroutine ufo_stericheight_simobs_tl_c(c_key_self, c_key_geovals, c_key_hofx) bind(c,name='ufo_stericheight_simobs_tl_f90')

implicit none
integer(c_int), intent(in) :: c_key_self
integer(c_int), intent(in) :: c_key_geovals
integer(c_int), intent(in) :: c_key_hofx

type(ufo_stericheight_tlad), pointer :: self
type(ufo_geovals),    pointer :: geovals
type(obs_vector),     pointer :: hofx

character(len=*), parameter :: myname_="ufo_stericheight_simobs_tl_c"

call ufo_stericheight_tlad_registry%get(c_key_self, self)
call ufo_geovals_registry%get(c_key_geovals,geovals)
call ioda_obs_vect_registry%get(c_key_hofx,hofx)

call ufo_stericheight_simobs_tl(self, geovals, hofx)

end subroutine ufo_stericheight_simobs_tl_c

! ------------------------------------------------------------------------------

subroutine ufo_stericheight_simobs_ad_c(c_key_self, c_key_geovals, c_key_hofx) bind(c,name='ufo_stericheight_simobs_ad_f90')

implicit none
integer(c_int), intent(in) :: c_key_self
integer(c_int), intent(in) :: c_key_geovals
integer(c_int), intent(in) :: c_key_hofx

type(ufo_stericheight_tlad), pointer :: self
type(ufo_geovals),    pointer :: geovals
type(obs_vector),     pointer :: hofx


character(len=*), parameter :: myname_="ufo_stericheight_simobs_ad_c"

call ufo_stericheight_tlad_registry%get(c_key_self, self)
call ufo_geovals_registry%get(c_key_geovals,geovals)
call ioda_obs_vect_registry%get(c_key_hofx,hofx)

call ufo_stericheight_simobs_ad(self, geovals, hofx)

end subroutine ufo_stericheight_simobs_ad_c

end module ufo_stericheight_tlad_mod_c
