! (C) Copyright 2017 UCAR
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 

!> Fortran example module for functions on the interface between C++ and Fortran
!  to handle observation operators

! TODO: replace "example" with your_observation_operator_name through the file

module ufo_example_mod_c
  
  use iso_c_binding
  use config_mod
  use ioda_obs_vectors,   only: obs_vector, ioda_obs_vect_registry
  use ufo_geovals_mod,   only: ufo_geovals
  use ufo_geovals_mod_c, only: ufo_geovals_registry
  use ioda_obsdb_mod,   only: ioda_obsdb
  use ioda_obsdb_mod_c, only: ioda_obsdb_registry 
  use ufo_example_mod 
  implicit none
  private
  
#define LISTED_TYPE ufo_example
  !> Linked list interface - defines registry_t type
#include "../../linkedList_i.f"
  
  !> Global registry
  type(registry_t) :: ufo_example_registry
  
  ! ------------------------------------------------------------------------------
contains
  ! ------------------------------------------------------------------------------
  !> Linked list implementation
#include "../../linkedList_c.f"
  
! ------------------------------------------------------------------------------
  
subroutine ufo_example_setup_c(c_key_self, c_conf) bind(c,name='ufo_example_setup_f90')
implicit none
integer(c_int), intent(inout) :: c_key_self
type(c_ptr), intent(in)    :: c_conf
    
type(ufo_example), pointer :: self

call ufo_example_registry%init()
call ufo_example_registry%add(c_key_self)
call ufo_example_registry%get(c_key_self, self)

! TODO: add call to your Fortran routine to setup the observation operator, if needed
!       (defined in ufo_<your_obs_operator_name>_mod.F90)
! example: call ufo_example_setup(self)
    
end subroutine ufo_example_setup_c
  
! ------------------------------------------------------------------------------
  
subroutine ufo_example_delete_c(c_key_self) bind(c,name='ufo_example_delete_f90')
implicit none
integer(c_int), intent(inout) :: c_key_self
    
type(ufo_example), pointer :: self

call ufo_example_registry%get(c_key_self, self)

! TODO: add call to your Fortran routine to destruct the observation operator, if needed
!       (defined in ufo_<your_obs_operator_name>_mod.F90)
! example: call ufo_example_delete(self)

call ufo_example_registry%remove(c_key_self)
    
end subroutine ufo_example_delete_c
  
! ------------------------------------------------------------------------------

subroutine ufo_example_eqv_c(c_key_self, c_key_geovals, c_key_obsspace, c_key_hofx, c_bias) bind(c,name='ufo_example_eqv_f90')

implicit none
integer(c_int), intent(in) :: c_key_self
integer(c_int), intent(in) :: c_key_geovals
integer(c_int), intent(in) :: c_key_hofx
integer(c_int), intent(in) :: c_key_obsspace
integer(c_int), intent(in) :: c_bias

type(ufo_example), pointer :: self
type(ufo_geovals), pointer :: geovals
type(obs_vector),  pointer :: hofx
type(ioda_obsdb),  pointer :: obss

character(len=*), parameter :: myname_="ufo_example_eqv_c"

call ufo_example_registry%get(c_key_self, self)
call ufo_geovals_registry%get(c_key_geovals,geovals)
call ioda_obs_vect_registry%get(c_key_hofx,hofx)
call ioda_obsdb_registry%get(c_key_obsspace,obss)

! TODO: replace with the call to your Fortran routine for observation operator
!       (defined in ufo_<your_obs_operator_name>_mod.F90)
call ufo_example_eqv(self, geovals, hofx, obss)

end subroutine ufo_example_eqv_c
  
end module ufo_example_mod_c
