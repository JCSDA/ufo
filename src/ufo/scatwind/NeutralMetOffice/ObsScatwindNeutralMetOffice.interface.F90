!-------------------------------------------------------------------------------
! (C) British Crown Copyright 2020 Met Office
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
!-------------------------------------------------------------------------------

!> Fortran module to handle scatwind observations - Met Office neutral operator

module ufo_scatwind_neutralmetoffice_mod_c
  
  use iso_c_binding
  use ufo_scatwind_neutralmetoffice_mod
  use ufo_geovals_mod,    only: ufo_geovals
  use ufo_geovals_mod_c,  only: ufo_geovals_registry

  implicit none
  private
  
#define LISTED_TYPE ufo_scatwind_neutralmetoffice
  
  !> Linked list interface - defines registry_t type
#include "oops/util/linkedList_i.f"
  
  !> Global registry
  type(registry_t) :: ufo_scatwind_neutralmetoffice_registry
  
  ! ------------------------------------------------------------------------------
contains
  ! ------------------------------------------------------------------------------
  !> Linked list implementation
#include "oops/util/linkedList_c.f"
  
! ------------------------------------------------------------------------------
  
subroutine ufo_scatwind_neutralmetoffice_setup_c(c_key_self, c_obsvars, &
                                                 c_geovars) bind(c,name='ufo_scatwind_neutralmetoffice_setup_f90')
use oops_variables_mod
implicit none
integer(c_int), intent(inout)  :: c_key_self
type(c_ptr), intent(in), value :: c_obsvars !< variables to be simulated
type(c_ptr), intent(in), value :: c_geovars !< variables requested from the model

type(ufo_scatwind_NeutralMetOffice), pointer :: self

call ufo_scatwind_neutralmetoffice_registry%setup(c_key_self, self)

self%obsvars = oops_variables(c_obsvars)
self%geovars = oops_variables(c_geovars)

call self%setup()

end subroutine ufo_scatwind_neutralmetoffice_setup_c
  
! ------------------------------------------------------------------------------
  
subroutine ufo_scatwind_neutralmetoffice_delete_c(c_key_self) bind(c,name='ufo_scatwind_neutralmetoffice_delete_f90')
implicit none
integer(c_int), intent(inout) :: c_key_self
    
type(ufo_scatwind_NeutralMetOffice), pointer :: self

call ufo_scatwind_NeutralMetOffice_registry%delete(c_key_self,self)

end subroutine ufo_scatwind_neutralmetoffice_delete_c
  
! ------------------------------------------------------------------------------

subroutine ufo_scatwind_neutralmetoffice_simobs_c(c_key_self, c_key_geovals, &
                                                  c_obsspace, c_nvars, c_nlocs, &
                                                  c_hofx) bind(c,name='ufo_scatwind_neutralmetoffice_simobs_f90')

implicit none
integer(c_int), intent(in) :: c_key_self
integer(c_int), intent(in) :: c_key_geovals
type(c_ptr), value, intent(in) :: c_obsspace
integer(c_int), intent(in)     :: c_nvars, c_nlocs
real(c_double), intent(inout)  :: c_hofx(c_nvars, c_nlocs)

type(ufo_scatwind_NeutralMetOffice),     pointer :: self
type(ufo_geovals),                       pointer :: geovals
character(len=*), parameter :: myname_="ufo_scatwind_neutralmetoffice_simobs_c"

call ufo_scatwind_NeutralMetOffice_registry%get(c_key_self, self)
call ufo_geovals_registry%get(c_key_geovals, geovals)

call self%simobs(geovals, c_obsspace, c_nvars, c_nlocs, c_hofx)

end subroutine ufo_scatwind_neutralmetoffice_simobs_c

! ------------------------------------------------------------------------------

end module ufo_scatwind_neutralmetoffice_mod_c
