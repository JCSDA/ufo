! (C) Copyright 2017 UCAR
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 

!> Fortran module to handle steric height observations

module ufo_stericheight_mod_c
  
  use iso_c_binding
  use config_mod
  use ufo_obs_data
  use ufo_obs_vectors,   only: obs_vector, ufo_obs_vect_registry
  use ufo_geovals_mod,   only: ufo_geovals, ufo_geovals_setup
  use ufo_geovals_mod_c, only: ufo_geovals_registry
  use ufo_obs_stericheight_mod,   only: ufo_obs_stericheight
  use ufo_obs_stericheight_mod_c, only: ufo_obs_stericheight_registry 
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

subroutine ufo_stericheight_eqv_c(c_key_self, c_key_geovals, c_key_obsspace, c_key_hofx, c_bias) bind(c,name='ufo_stericheight_eqv_f90')

implicit none
integer(c_int), intent(in) :: c_key_self
integer(c_int), intent(in) :: c_key_geovals
integer(c_int), intent(in) :: c_key_hofx
integer(c_int), intent(in) :: c_key_obsspace
integer(c_int), intent(in) :: c_bias

type(ufo_stericheight), pointer :: self
type(ufo_geovals),    pointer :: geovals
type(obs_vector),     pointer :: hofx

character(len=*), parameter :: myname_="ufo_stericheight_eqv_c"

print *,myname_

call ufo_stericheight_registry%get(c_key_self, self)
call ufo_geovals_registry%get(c_key_geovals,geovals)
call ufo_obs_vect_registry%get(c_key_hofx,hofx)

call ufo_stericheight_eqv(self, geovals, hofx)

end subroutine ufo_stericheight_eqv_c

! ------------------------------------------------------------------------------

subroutine ufo_stericheight_gettraj(c_key_self, c_nobs, c_vars, c_key_traj) bind(c,name='ufo_stericheight_gettraj_f90')
use fckit_log_module, only : fckit_log

implicit none
integer(c_int), intent(in) :: c_key_self
integer(c_int), intent(in) :: c_nobs
type(c_ptr), intent(in)       :: c_vars
integer(c_int), intent(inout) :: c_key_traj

type(ufo_stericheight), pointer :: self
type(ufo_geovals), pointer :: traj
type(ufo_vars) :: cvars
integer nobs

call ufo_stericheight_registry%get(c_key_self, self)
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

end subroutine ufo_stericheight_gettraj

! ------------------------------------------------------------------------------

subroutine ufo_stericheight_settraj_c(c_key_self, c_key_geovals) bind(c,name='ufo_stericheight_settraj_f90')

implicit none
integer(c_int), intent(in) :: c_key_self
integer(c_int), intent(in) :: c_key_geovals

type(ufo_stericheight), pointer :: self
type(ufo_geovals),    pointer :: geovals

character(len=*), parameter :: myname_="ufo_stericheight_settraj_c"

print *, myname_

call ufo_stericheight_registry%get(c_key_self, self)
call ufo_geovals_registry%get(c_key_geovals,geovals)

call ufo_stericheight_settraj(self, geovals)

end subroutine ufo_stericheight_settraj_c

! ------------------------------------------------------------------------------

subroutine ufo_stericheight_eqv_tl_c(c_key_self, c_key_geovals, c_key_hofx) bind(c,name='ufo_stericheight_eqv_tl_f90')

implicit none
integer(c_int), intent(in) :: c_key_self
integer(c_int), intent(in) :: c_key_geovals
integer(c_int), intent(in) :: c_key_hofx

type(ufo_stericheight), pointer :: self
type(ufo_geovals),    pointer :: geovals
type(obs_vector),     pointer :: hofx

character(len=*), parameter :: myname_="ufo_stericheight_eqv_tl_c"

call ufo_stericheight_registry%get(c_key_self, self)
call ufo_geovals_registry%get(c_key_geovals,geovals)
call ufo_obs_vect_registry%get(c_key_hofx,hofx)

call ufo_stericheight_eqv_tl(self, geovals, hofx)

end subroutine ufo_stericheight_eqv_tl_c

! ------------------------------------------------------------------------------

subroutine ufo_stericheight_eqv_ad_c(c_key_self, c_key_geovals, c_key_hofx) bind(c,name='ufo_stericheight_eqv_ad_f90')

implicit none
integer(c_int), intent(in) :: c_key_self
integer(c_int), intent(in) :: c_key_geovals
integer(c_int), intent(in) :: c_key_hofx

type(ufo_stericheight), pointer :: self
type(ufo_geovals),    pointer :: geovals
type(obs_vector),     pointer :: hofx


character(len=*), parameter :: myname_="ufo_stericheight_eqv_ad_c"

call ufo_stericheight_registry%get(c_key_self, self)
call ufo_geovals_registry%get(c_key_geovals,geovals)
call ufo_obs_vect_registry%get(c_key_hofx,hofx)

call ufo_stericheight_eqv_ad(self, geovals, hofx)

end subroutine ufo_stericheight_eqv_ad_c

! ------------------------------------------------------------------------------

!subroutine ufo_obs_get(c_key_self, lreq, c_req, lcol, c_col, c_key_ovec) bind(c,name='ufo_obsdbsic_get_f90')
subroutine ufo_obs_get(c_key_self, lcol, c_col, c_key_ovec) bind(c,name='ufo_obsdbsteric_get_f90')  
use  ufo_obs_stericheight_mod
implicit none
integer(c_int), intent(in) :: c_key_self
integer(c_int), intent(in) :: lcol
character(kind=c_char,len=1), intent(in) :: c_col(lcol+1)
integer(c_int), intent(in) :: c_key_ovec

type(ufo_obs_stericheight), pointer :: self
type(obs_vector), pointer :: ovec
character(len=lcol) :: col

call ufo_obs_stericheight_registry%get(c_key_self, self)
call ufo_obs_vect_registry%get(c_key_ovec,ovec)
!call c_f_string(c_req, req)
!call c_f_string(c_col, col)

!call obs_get(self, trim(req), trim(col), ovec)


ovec%nobs = self%nobs
if (c_col(5)//c_col(6)=='rr') then
   ovec%values = 0.1 !self%icefrac_err
!   print *, self%icefrac_err
else
   ovec%values = self%adt
end if


end subroutine ufo_obs_get

  
end module ufo_stericheight_mod_c
