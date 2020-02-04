
! (C) Copyright 2017-2018 UCAR
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 

module ufo_onedvarfortran_mod_c

use iso_c_binding
use ufo_onedvarfortran_mod
use ufo_onedvarfortran_utils_mod, only: max_string_length
use ufo_geovals_mod
use ufo_geovals_mod_c,   only: ufo_geovals_registry

implicit none
private

#define LISTED_TYPE ufo_onedvarfortran

!> Linked list interface - defines registry_t type
#include "oops/util/linkedList_i.f"

!> Global registry
type(registry_t) :: ufo_onedvarfortran_registry

! ------------------------------------------------------------------------------
contains
! ------------------------------------------------------------------------------
!> Linked list implementation
#include "oops/util/linkedList_c.f"
! ------------------------------------------------------------------------------

subroutine ufo_onedvarfortran_create_c(c_self, c_obspace, c_conf, c_nchan, c_channels) bind(c,name='ufo_onedvarfortran_create_f90')
implicit none
integer(c_int), intent(inout)  :: c_self
type(c_ptr), value, intent(in) :: c_obspace
type(c_ptr), value, intent(in) :: c_conf
integer(c_int), intent(in) :: c_nchan
integer(c_int), intent(in) :: c_channels(c_nchan)

type(ufo_onedvarfortran), pointer :: self

call ufo_onedvarfortran_registry%setup(c_self, self)
call ufo_onedvarfortran_create(self, c_obspace, c_conf, c_channels)

end subroutine ufo_onedvarfortran_create_c

! ------------------------------------------------------------------------------

subroutine ufo_onedvarfortran_delete_c(c_self) bind(c,name='ufo_onedvarfortran_delete_f90')
implicit none
integer(c_int), intent(inout) :: c_self

type(ufo_onedvarfortran), pointer :: self

call ufo_onedvarfortran_registry%get(c_self, self)
call ufo_onedvarfortran_delete(self)
call ufo_onedvarfortran_registry%delete(c_self, self)

end subroutine ufo_onedvarfortran_delete_c

! ------------------------------------------------------------------------------

subroutine ufo_onedvarfortran_prior_c(c_self, c_geovals) bind(c,name='ufo_onedvarfortran_prior_f90')
implicit none
integer(c_int), intent(in) :: c_self
integer(c_int), intent(in) :: c_geovals

type(ufo_onedvarfortran), pointer :: self
type(ufo_geovals), pointer :: geovals

call ufo_onedvarfortran_registry%get(c_self, self)
call ufo_geovals_registry%get(c_geovals, geovals)

call ufo_onedvarfortran_prior(self, geovals)

end subroutine ufo_onedvarfortran_prior_c

! ------------------------------------------------------------------------------

subroutine ufo_onedvarfortran_post_c(c_self, c_nvar, c_nloc, c_hofx, c_var, c_geovals, c_conf) bind(c,name='ufo_onedvarfortran_post_f90')
implicit none
integer(c_int), intent(in) :: c_self
integer(c_int), intent(in) :: c_nvar
integer(c_int), intent(in) :: c_nloc
real(c_double), intent(in) :: c_hofx(c_nvar, c_nloc)
type(c_ptr), value, intent(in) :: c_var
integer(c_int), intent(in) :: c_geovals
type(c_ptr), value, intent(in) :: c_conf

type(ufo_onedvarfortran), pointer :: self
type(ufo_geovals), pointer :: geovals
character(len=max_string_length), allocatable :: vars(:)

allocate(vars(c_nvar))
vars = config_get_string_vector(c_var, max_string_length, "variables")

call ufo_onedvarfortran_registry%get(c_self, self)
call ufo_geovals_registry%get(c_geovals, geovals)

call ufo_onedvarfortran_post(self, c_hofx, vars, geovals, c_conf)

end subroutine ufo_onedvarfortran_post_c

! ------------------------------------------------------------------------------

end module ufo_onedvarfortran_mod_c

