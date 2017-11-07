! (C) Copyright 2009-2016 ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
! In applying this licence, ECMWF does not waive the privileges and immunities 
! granted to it by virtue of its status as an intergovernmental organisation nor
! does it submit to any jurisdiction.

!> Fortran module to handle wind speed observations for the QG model
module qg_wspeed_mod

use iso_c_binding
use config_mod
use duration_mod
use qg_obs_data
use qg_obs_vectors
use qg_obsoper_mod
use qg_vars_mod
use qg_locs_mod
use qg_goms_mod
use kinds

implicit none
private

! ------------------------------------------------------------------------------
contains
! ------------------------------------------------------------------------------

subroutine c_qg_wspeed_setup(c_key_self, c_conf) bind(c,name='qg_wspeed_setup_f90')
implicit none
integer(c_int), intent(inout) :: c_key_self
type(c_ptr), intent(in)    :: c_conf

type(qg_obsoper), pointer :: self
character(len=1) :: svars(2) = (/"u","v"/)

call qg_obsoper_registry%init()
call qg_obsoper_registry%add(c_key_self)
call qg_obsoper_registry%get(c_key_self, self)

call qg_oper_setup(self, c_conf, svars, 1)

end subroutine c_qg_wspeed_setup

! ------------------------------------------------------------------------------

subroutine c_qg_wspeed_delete(c_key_self) bind(c,name='qg_wspeed_delete_f90')
implicit none
integer(c_int), intent(inout) :: c_key_self

type(qg_obsoper), pointer :: self

call qg_obsoper_registry%get(c_key_self, self)
deallocate(self%varin%fldnames)
call qg_obsoper_registry%remove(c_key_self)

end subroutine c_qg_wspeed_delete

! ------------------------------------------------------------------------------
subroutine qg_wspeed_eqv(c_key_gom, c_key_hofx, c_bias) bind(c,name='qg_wspeed_eqv_f90')
implicit none
integer(c_int), intent(in) :: c_key_gom
integer(c_int), intent(in) :: c_key_hofx
real(c_double), intent(in) :: c_bias
type(qg_goms), pointer  :: gom
type(obs_vect), pointer :: hofx
integer :: io, jo
character(len=250) :: record
real(kind=kind_real) :: zz

if (abs(c_bias) > epsilon(c_bias)) call abor1_ftn ("qg_wspeed: bias not implemented")

call qg_goms_registry%get(c_key_gom,gom)
call qg_obs_vect_registry%get(c_key_hofx,hofx)

do jo=1,gom%nobs
  io=gom%indx(jo)
  zz=gom%values(1,jo)*gom%values(1,jo)+gom%values(2,jo)*gom%values(2,jo)
  hofx%values(1,io)=sqrt(zz)
enddo

end subroutine qg_wspeed_eqv
! ------------------------------------------------------------------------------
subroutine qg_wspeed_equiv_tl(c_key_gom, c_key_hofx, c_key_traj, c_bias) &
 & bind(c,name='qg_wspeed_equiv_tl_f90')
implicit none
integer(c_int), intent(in) :: c_key_gom
integer(c_int), intent(in) :: c_key_hofx
integer(c_int), intent(in) :: c_key_traj
real(c_double), intent(in) :: c_bias
type(qg_goms), pointer  :: gom
type(obs_vect), pointer :: hofx
type(qg_goms), pointer  :: traj
integer :: io, jo
real(kind=kind_real) :: zz, zu, zv, zt

call qg_goms_registry%get(c_key_gom,gom)
call qg_obs_vect_registry%get(c_key_hofx,hofx)
call qg_goms_registry%get(c_key_traj,traj)

do jo=1,gom%nobs
  io=gom%indx(jo)
  zu=traj%values(1,io)
  zv=traj%values(2,io)
  zt=sqrt(zu*zu+zv*zv)
  if (zt>epsilon(zt)) then
    zz=2.0_kind_real*zu*gom%values(1,jo) &
      +2.0_kind_real*zv*gom%values(2,jo)
    hofx%values(1,io)=zz/(2.0_kind_real*zt)
  else
    hofx%values(1,io)=0.0_kind_real
  endif
enddo

end subroutine qg_wspeed_equiv_tl
! ------------------------------------------------------------------------------
subroutine qg_wspeed_equiv_ad(c_key_gom, c_key_hofx, c_key_traj, c_bias) &
 & bind(c,name='qg_wspeed_equiv_ad_f90')
implicit none
integer(c_int), intent(in) :: c_key_gom
integer(c_int), intent(in) :: c_key_hofx
integer(c_int), intent(in) :: c_key_traj
real(c_double), intent(inout) :: c_bias
type(qg_goms), pointer  :: gom
type(obs_vect), pointer :: hofx
type(qg_goms), pointer  :: traj
integer :: io, jo
real(kind=kind_real) :: zz, zu, zv, zt

call qg_goms_registry%get(c_key_gom,gom)
call qg_obs_vect_registry%get(c_key_hofx,hofx)
call qg_goms_registry%get(c_key_traj,traj)

do jo=1,gom%nobs
  io=gom%indx(jo)
  zu=traj%values(1,io)
  zv=traj%values(2,io)
  zt=sqrt(zu*zu+zv*zv)
  if (zt>epsilon(zt)) then
    zz=hofx%values(1,io)/(2.0_kind_real*zt)
    gom%values(1,jo)=2.0_kind_real*zu*zz
    gom%values(2,jo)=2.0_kind_real*zv*zz
  else
    gom%values(1,jo)=0.0_kind_real
    gom%values(2,jo)=0.0_kind_real
  endif
enddo

end subroutine qg_wspeed_equiv_ad
! ------------------------------------------------------------------------------
subroutine qg_wspeed_gettraj(c_key_self, c_nobs, c_key_traj) bind(c,name='qg_wspeed_gettraj_f90')
use fckit_log_module, only : fckit_log
implicit none
integer(c_int), intent(in) :: c_key_self
integer(c_int), intent(in) :: c_nobs
integer(c_int), intent(inout) :: c_key_traj

type(qg_obsoper), pointer :: self
type(qg_goms), pointer :: traj
integer, allocatable :: mobs(:)
integer :: jj

call qg_obsoper_registry%get(c_key_self, self)

call qg_goms_registry%init()
call qg_goms_registry%add(c_key_traj)
call qg_goms_registry%get(c_key_traj,traj)
allocate(mobs(c_nobs))
do jj=1,c_nobs
  mobs(jj)=jj
enddo
call gom_setup(traj, self%varin, mobs)
deallocate(mobs)

end subroutine qg_wspeed_gettraj
! ------------------------------------------------------------------------------
subroutine qg_wspeed_settraj(c_key_gom, c_key_traj) bind(c,name='qg_wspeed_settraj_f90')
use fckit_log_module, only : fckit_log
implicit none
integer(c_int), intent(in) :: c_key_gom
integer(c_int), intent(in) :: c_key_traj

type(qg_goms), pointer :: gom
type(qg_goms), pointer :: traj
integer :: jo, io

call qg_goms_registry%get(c_key_gom,gom)
call qg_goms_registry%get(c_key_traj,traj)

do jo=1,gom%nobs
  io=gom%indx(jo)
  traj%values(1,io)=gom%values(1,jo)
  traj%values(2,io)=gom%values(2,jo)
enddo

end subroutine qg_wspeed_settraj
! ------------------------------------------------------------------------------

end module qg_wspeed_mod
