! (C) Copyright 2009-2016 ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
! In applying this licence, ECMWF does not waive the privileges and immunities 
! granted to it by virtue of its status as an intergovernmental organisation nor
! does it submit to any jurisdiction.

!> Fortran module to handle wind speed observations for the QG model
module ufo_wspeed_mod

use iso_c_binding
use config_mod
use duration_mod
use ufo_obs_data
use ufo_obs_vectors
use ufo_obsoper_mod
use ufo_vars_mod
use ufo_locs_mod
use ufo_geovals_mod
use kinds

implicit none
private

! ------------------------------------------------------------------------------
contains
! ------------------------------------------------------------------------------

subroutine c_ufo_wspeed_setup(c_key_self, c_conf) bind(c,name='ufo_wspeed_setup_f90')
implicit none
integer(c_int), intent(inout) :: c_key_self
type(c_ptr), intent(in)    :: c_conf

type(ufo_obsoper), pointer :: self
character(len=1) :: svars(2) = (/"u","v"/)

call ufo_obsoper_registry%init()
call ufo_obsoper_registry%add(c_key_self)
call ufo_obsoper_registry%get(c_key_self, self)

call ufo_oper_setup(self, c_conf, svars, 1)

end subroutine c_ufo_wspeed_setup

! ------------------------------------------------------------------------------

subroutine c_ufo_wspeed_delete(c_key_self) bind(c,name='ufo_wspeed_delete_f90')
implicit none
integer(c_int), intent(inout) :: c_key_self

type(ufo_obsoper), pointer :: self

call ufo_obsoper_registry%get(c_key_self, self)
deallocate(self%varin%fldnames)
call ufo_obsoper_registry%remove(c_key_self)

end subroutine c_ufo_wspeed_delete

! ------------------------------------------------------------------------------
subroutine ufo_wspeed_eqv(c_key_geovals, c_key_hofx, c_bias) bind(c,name='ufo_wspeed_eqv_f90')
implicit none
integer(c_int), intent(in) :: c_key_geovals
integer(c_int), intent(in) :: c_key_hofx
real(c_double), intent(in) :: c_bias
type(ufo_geovals), pointer  :: geovals
type(obs_vect), pointer :: hofx
integer :: io, jo
character(len=250) :: record
real(kind=kind_real) :: zz

if (abs(c_bias) > epsilon(c_bias)) call abor1_ftn ("ufo_wspeed: bias not implemented")

!call ufo_geovals_registry%get(c_key_geovals,geovals)
call ufo_obs_vect_registry%get(c_key_hofx,hofx)


end subroutine ufo_wspeed_eqv
! ------------------------------------------------------------------------------
subroutine ufo_wspeed_equiv_tl(c_key_geovals, c_key_hofx, c_key_traj, c_bias) &
 & bind(c,name='ufo_wspeed_equiv_tl_f90')
implicit none
integer(c_int), intent(in) :: c_key_geovals
integer(c_int), intent(in) :: c_key_hofx
integer(c_int), intent(in) :: c_key_traj
real(c_double), intent(in) :: c_bias
type(ufo_geovals), pointer  :: gom
type(obs_vect), pointer :: hofx
type(ufo_geovals), pointer  :: traj
integer :: io, jo
real(kind=kind_real) :: zz, zu, zv, zt

!call ufo_geovals_registry%get(c_key_geovals,geovals)
!call ufo_obs_vect_registry%get(c_key_hofx,hofx)
!call ufo_geovals_registry%get(c_key_traj,traj)


end subroutine ufo_wspeed_equiv_tl
! ------------------------------------------------------------------------------
subroutine ufo_wspeed_equiv_ad(c_key_gom, c_key_hofx, c_key_traj, c_bias) &
 & bind(c,name='ufo_wspeed_equiv_ad_f90')
implicit none
integer(c_int), intent(in) :: c_key_gom
integer(c_int), intent(in) :: c_key_hofx
integer(c_int), intent(in) :: c_key_traj
real(c_double), intent(inout) :: c_bias
type(ufo_geovals), pointer  :: gom
type(obs_vect), pointer :: hofx
type(ufo_geovals), pointer  :: traj
integer :: io, jo
real(kind=kind_real) :: zz, zu, zv, zt

call ufo_geovals_registry%get(c_key_gom,gom)
call ufo_obs_vect_registry%get(c_key_hofx,hofx)
call ufo_geovals_registry%get(c_key_traj,traj)

end subroutine ufo_wspeed_equiv_ad
! ------------------------------------------------------------------------------
subroutine ufo_wspeed_gettraj(c_key_self, c_nobs, c_key_traj) bind(c,name='ufo_wspeed_gettraj_f90')
use fckit_log_module, only : fckit_log
implicit none
integer(c_int), intent(in) :: c_key_self
integer(c_int), intent(in) :: c_nobs
integer(c_int), intent(inout) :: c_key_traj

end subroutine ufo_wspeed_gettraj
! ------------------------------------------------------------------------------
subroutine ufo_wspeed_settraj(c_key_gom, c_key_traj) bind(c,name='ufo_wspeed_settraj_f90')
use fckit_log_module, only : fckit_log
implicit none
integer(c_int), intent(in) :: c_key_gom
integer(c_int), intent(in) :: c_key_traj


end subroutine ufo_wspeed_settraj
! ------------------------------------------------------------------------------

end module ufo_wspeed_mod
