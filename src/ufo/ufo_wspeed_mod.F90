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
  use ufo_obs_vectors
  use ufo_vars_mod
  use ufo_locs_mod
  use ufo_geovals_mod
  use kinds
  
  implicit none
  private
  
  ! ------------------------------------------------------------------------------
  
  !> Fortran derived type for stream function observations for the QG model
  type :: ufo_obsoper
    integer :: nothing_here
  end type ufo_obsoper
  
#define LISTED_TYPE ufo_obsoper
  
  !> Linked list interface - defines registry_t type
#include "linkedList_i.f"
  
  !> Global registry
  type(registry_t) :: ufo_wspeed_registry
  
  ! ------------------------------------------------------------------------------
contains
  ! ------------------------------------------------------------------------------
  !> Linked list implementation
#include "linkedList_c.f"
  
  ! ------------------------------------------------------------------------------
  
  subroutine c_ufo_wspeed_setup(c_key_self, c_conf) bind(c,name='ufo_wspeed_setup_f90')
    implicit none
    integer(c_int), intent(inout) :: c_key_self
    type(c_ptr), intent(in)    :: c_conf
    
    type(ufo_obsoper), pointer :: self
    
    call ufo_wspeed_registry%init()
    call ufo_wspeed_registry%add(c_key_self)
    call ufo_wspeed_registry%get(c_key_self, self)
    
  end subroutine c_ufo_wspeed_setup
  
  ! ------------------------------------------------------------------------------
  
  subroutine c_ufo_wspeed_delete(c_key_self) bind(c,name='ufo_wspeed_delete_f90')
    implicit none
    integer(c_int), intent(inout) :: c_key_self
    
    type(ufo_obsoper), pointer :: self
    
    call ufo_wspeed_registry%get(c_key_self, self)
    call ufo_wspeed_registry%remove(c_key_self)
    
  end subroutine c_ufo_wspeed_delete
  
  ! ------------------------------------------------------------------------------
  
  subroutine ufo_wspeed_eqv(c_key_geovals, c_key_hofx, c_bias) bind(c,name='ufo_wspeed_eqv_f90')
    implicit none
    integer(c_int), intent(in) :: c_key_geovals
    integer(c_int), intent(in) :: c_key_hofx
    real(c_double), intent(in) :: c_bias
    type(ufo_geovals), pointer  :: geovals
    type(obs_vector), pointer :: hofx
    call ufo_obs_vect_registry%get(c_key_hofx,hofx)
    hofx%values(:) = 1.0
  end subroutine ufo_wspeed_eqv
  
  ! ------------------------------------------------------------------------------
  
  subroutine c_ufo_wspeed_inputs(c_key_self, c_key_vars) bind(c,name='ufo_wspeed_inputs_f90')
    implicit none
    integer(c_int), intent(in)    :: c_key_self
    integer(c_int), intent(inout) :: c_key_vars
    
    type(ufo_obsoper), pointer :: self
    type(ufo_vars), pointer :: vars
    
    call ufo_wspeed_registry%get(c_key_self, self)
    call ufo_vars_registry%init()
    call ufo_vars_registry%add(c_key_vars)
    call ufo_vars_registry%get(c_key_vars, vars)
    
  end subroutine c_ufo_wspeed_inputs
  
  ! ------------------------------------------------------------------------------
  subroutine ufo_wspeed_equiv_tl(c_key_geovals, c_key_hofx, c_key_traj, c_bias) &
       & bind(c,name='ufo_wspeed_equiv_tl_f90')
    implicit none
    integer(c_int), intent(in) :: c_key_geovals
    integer(c_int), intent(in) :: c_key_hofx
    integer(c_int), intent(in) :: c_key_traj
    real(c_double), intent(in) :: c_bias
  end subroutine ufo_wspeed_equiv_tl
  ! ------------------------------------------------------------------------------
  subroutine ufo_wspeed_equiv_ad(c_key_gom, c_key_hofx, c_key_traj, c_bias) &
       & bind(c,name='ufo_wspeed_equiv_ad_f90')
    implicit none
    integer(c_int), intent(in) :: c_key_gom
    integer(c_int), intent(in) :: c_key_hofx
    integer(c_int), intent(in) :: c_key_traj
    real(c_double), intent(inout) :: c_bias
  end subroutine ufo_wspeed_equiv_ad
  ! ------------------------------------------------------------------------------
  
end module ufo_wspeed_mod
