! (C) Copyright 2017-2018 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

!> Fortran seaicethickness module for functions on the interface between C++ and Fortran
!  to handle observation operators

module ufo_seaicethickness_mod_c

  use fckit_configuration_module, only: fckit_configuration 
  use iso_c_binding
  use ufo_seaicethickness_mod 
  implicit none
  private

  ! ------------------------------------------------------------------------------
#define LISTED_TYPE ufo_seaicethickness

  !> Linked list interface - defines registry_t type
#include "oops/util/linkedList_i.f"

  !> Global registry
  type(registry_t) :: ufo_seaicethickness_registry

  ! ------------------------------------------------------------------------------

contains

  ! ------------------------------------------------------------------------------
  !> Linked list implementation
#include "oops/util/linkedList_c.f"

! ------------------------------------------------------------------------------

subroutine ufo_seaicethickness_setup_c(c_key_self, c_conf, c_obsvars) bind(c,name='ufo_seaicethickness_setup_f90')
use obs_variables_mod
implicit none
integer(c_int), intent(inout)  :: c_key_self
type(c_ptr), value, intent(in) :: c_conf
type(c_ptr), value, intent(in) :: c_obsvars ! variables to be simulated                                                        

type(ufo_seaicethickness), pointer :: self
type(fckit_configuration) :: f_conf

call ufo_seaicethickness_registry%setup(c_key_self, self)
f_conf = fckit_configuration(c_conf)

self%obsvars = obs_variables(c_obsvars)

call self%setup(f_conf)

end subroutine ufo_seaicethickness_setup_c

! ------------------------------------------------------------------------------

subroutine ufo_seaicethickness_delete_c(c_key_self) bind(c,name='ufo_seaicethickness_delete_f90')
implicit none
integer(c_int), intent(inout) :: c_key_self
    
type(ufo_seaicethickness), pointer :: self

call ufo_seaicethickness_registry%get(c_key_self, self)

call self%delete()

call ufo_seaicethickness_registry%remove(c_key_self)

end subroutine ufo_seaicethickness_delete_c

! ------------------------------------------------------------------------------

subroutine ufo_seaicethickness_simobs_c(c_key_self, c_key_geovals, c_obsspace, c_nobs, c_hofx) &
    bind(c,name='ufo_seaicethickness_simobs_f90')

implicit none
integer(c_int), intent(in) :: c_key_self
integer(c_int), intent(in) :: c_key_geovals
type(c_ptr), value, intent(in) :: c_obsspace
integer(c_int), intent(in) :: c_nobs
real(c_double), intent(inout) :: c_hofx(c_nobs)

type(ufo_seaicethickness), pointer :: self

call ufo_seaicethickness_registry%get(c_key_self, self)
call self%opr_simobs(c_key_geovals, c_obsspace, c_hofx)

end subroutine ufo_seaicethickness_simobs_c

! ------------------------------------------------------------------------------

end module ufo_seaicethickness_mod_c
