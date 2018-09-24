! (C) Copyright 2017-2018 UCAR
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 

!> Fortran example module for functions on the interface between C++ and Fortran
!  to handle tl/ad observation operators

! TODO: replace "example" with your_observation_operator_name through the file

module ufo_example_tlad_mod_c
  
  use iso_c_binding
  use config_mod
  use ioda_obs_vectors,   only: obs_vector, ioda_obs_vect_registry
  use ufo_geovals_mod,   only: ufo_geovals
  use ufo_geovals_mod_c, only: ufo_geovals_registry
  use ioda_obsdb_mod,   only: ioda_obsdb
  use ioda_obsdb_mod_c, only: ioda_obsdb_registry 
  use ufo_example_tlad_mod 
  implicit none
  private
  
#define LISTED_TYPE ufo_example_tlad
  
  !> Linked list interface - defines registry_t type
#include "../../linkedList_i.f"
  
  !> Global registry
  type(registry_t) :: ufo_example_tlad_registry
  
  ! ------------------------------------------------------------------------------
contains
  ! ------------------------------------------------------------------------------
  !> Linked list implementation
#include "../../linkedList_c.f"
  
! ------------------------------------------------------------------------------
  
subroutine ufo_example_tlad_setup_c(c_key_self, c_conf) bind(c,name='ufo_example_tlad_setup_f90')
implicit none
integer(c_int), intent(inout) :: c_key_self
type(c_ptr), intent(in)    :: c_conf
    
type(ufo_example_tlad), pointer :: self

call ufo_example_tlad_registry%init()
call ufo_example_tlad_registry%add(c_key_self)
call ufo_example_tlad_registry%get(c_key_self, self)

! TODO: add call to your Fortran routine to setup the tl/ad observation operator, if needed
!       (defined in ufo_<your_obs_operator_name>_tlad_mod.F90)
! example: call ufo_example_tlad_setup(self)
    
end subroutine ufo_example_tlad_setup_c
  
! ------------------------------------------------------------------------------
  
subroutine ufo_example_tlad_delete_c(c_key_self) bind(c,name='ufo_example_tlad_delete_f90')
implicit none
integer(c_int), intent(inout) :: c_key_self
    
type(ufo_example_tlad), pointer :: self

call ufo_example_tlad_registry%get(c_key_self, self)

! TODO: replace with the call to your Fortran routine to destruct the tl/ad observation operator
!       (defined in ufo_<your_obs_operator_name>_tlad_mod.F90)
call ufo_example_tlad_delete(self)

call ufo_example_tlad_registry%remove(c_key_self)
    
end subroutine ufo_example_tlad_delete_c
  
! ------------------------------------------------------------------------------

subroutine ufo_example_tlad_settraj_c(c_key_self, c_key_geovals, c_key_obsspace) bind(c,name='ufo_example_tlad_settraj_f90')

implicit none
integer(c_int), intent(in) :: c_key_self
integer(c_int), intent(in) :: c_key_geovals
integer(c_int), intent(in) :: c_key_obsspace

type(ufo_example_tlad), pointer :: self
type(ufo_geovals),      pointer :: geovals
type(ioda_obsdb),       pointer :: obss

character(len=*), parameter :: myname_="ufo_example_tlad_settraj_c"

call ufo_example_tlad_registry%get(c_key_self, self)
call ufo_geovals_registry%get(c_key_geovals,geovals)
call ioda_obsdb_registry%get(c_key_obsspace,obss)

! TODO: replace with the call to your Fortran routine for setting tl/ad trajectory
!       (defined in ufo_<your_obs_operator_name>_tlad_mod.F90)
call ufo_example_tlad_settraj(self, geovals, obss)

end subroutine ufo_example_tlad_settraj_c

! ------------------------------------------------------------------------------

subroutine ufo_example_simobs_tl_c(c_key_self, c_key_geovals, c_key_obsspace, c_key_hofx) bind(c,name='ufo_example_simobs_tl_f90')

implicit none
integer(c_int), intent(in) :: c_key_self
integer(c_int), intent(in) :: c_key_geovals
integer(c_int), intent(in) :: c_key_hofx
integer(c_int), intent(in) :: c_key_obsspace

type(ufo_example_tlad), pointer :: self
type(ufo_geovals),      pointer :: geovals
type(obs_vector),       pointer :: hofx
type(ioda_obsdb),       pointer :: obss

character(len=*), parameter :: myname_="ufo_example_simobs_tl_c"

call ufo_example_tlad_registry%get(c_key_self, self)
call ufo_geovals_registry%get(c_key_geovals,geovals)
call ioda_obs_vect_registry%get(c_key_hofx,hofx)
call ioda_obsdb_registry%get(c_key_obsspace,obss)

! TODO: replace with the call to your Fortran routine for tl obs operator
!       (defined in ufo_<your_obs_operator_name>_tlad_mod.F90)
call ufo_example_simobs_tl(self, geovals, hofx, obss)

end subroutine ufo_example_simobs_tl_c

! ------------------------------------------------------------------------------

subroutine ufo_example_simobs_ad_c(c_key_self, c_key_geovals, c_key_obsspace, c_key_hofx) bind(c,name='ufo_example_simobs_ad_f90')

implicit none
integer(c_int), intent(in) :: c_key_self
integer(c_int), intent(in) :: c_key_geovals
integer(c_int), intent(in) :: c_key_hofx
integer(c_int), intent(in) :: c_key_obsspace

type(ufo_example_tlad), pointer :: self
type(ufo_geovals),      pointer :: geovals
type(obs_vector),       pointer :: hofx
type(ioda_obsdb),       pointer :: obss

character(len=*), parameter :: myname_="ufo_example_simobs_ad_c"

call ufo_example_tlad_registry%get(c_key_self, self)
call ufo_geovals_registry%get(c_key_geovals,geovals)
call ioda_obs_vect_registry%get(c_key_hofx,hofx)
call ioda_obsdb_registry%get(c_key_obsspace,obss)

! TODO: replace with the call to your Fortran routine for ad obs operator
!       (defined in ufo_<your_obs_operator_name>_tlad_mod.F90)
call ufo_example_simobs_ad(self, geovals, hofx, obss)

end subroutine ufo_example_simobs_ad_c

  
end module ufo_example_tlad_mod_c
