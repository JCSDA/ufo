
! (C) Copyright 2017-2020 Met Office
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 

module ufo_rttovonedvarcheck_mod_c

use fckit_configuration_module, only: fckit_configuration
use iso_c_binding
use oops_variables_mod
use ufo_geovals_mod
use ufo_geovals_mod_c,   only: ufo_geovals_registry
use ufo_rttovonedvarcheck_mod

implicit none
private

#define LISTED_TYPE ufo_rttovonedvarcheck

!> Linked list interface - defines registry_t type
#include "oops/util/linkedList_i.f"

!> Global registry
type(registry_t) :: ufo_rttovonedvarcheck_registry

! ------------------------------------------------------------------------------
contains
! ------------------------------------------------------------------------------
!> Linked list implementation
#include "oops/util/linkedList_c.f"
! ------------------------------------------------------------------------------

subroutine ufo_rttovonedvarcheck_create_c(c_self, c_obspace, c_conf, c_nchan, c_channels) &
                                            bind(c,name='ufo_rttovonedvarcheck_create_f90')
implicit none
integer(c_int), intent(inout)  :: c_self
type(c_ptr), value, intent(in) :: c_obspace
type(c_ptr), value, intent(in) :: c_conf
integer(c_int), intent(in) :: c_nchan
integer(c_int), intent(in) :: c_channels(c_nchan)

type(ufo_rttovonedvarcheck), pointer :: self
type(fckit_configuration) :: f_conf

call ufo_rttovonedvarcheck_registry%setup(c_self, self)
f_conf = fckit_configuration(c_conf)

call ufo_rttovonedvarcheck_create(self, c_obspace, f_conf, c_channels)

end subroutine ufo_rttovonedvarcheck_create_c

! ------------------------------------------------------------------------------

subroutine ufo_rttovonedvarcheck_delete_c(c_self) &
                      bind(c,name='ufo_rttovonedvarcheck_delete_f90')
implicit none
integer(c_int), intent(inout) :: c_self

type(ufo_rttovonedvarcheck), pointer :: self

call ufo_rttovonedvarcheck_registry%get(c_self, self)
call ufo_rttovonedvarcheck_delete(self)
call ufo_rttovonedvarcheck_registry%delete(c_self, self)

end subroutine ufo_rttovonedvarcheck_delete_c

! ------------------------------------------------------------------------------

subroutine ufo_rttovonedvarcheck_apply_c(c_self, c_vars, c_geovals, c_nobs, c_apply) &
               bind(c,name='ufo_rttovonedvarcheck_apply_f90')
implicit none
integer(c_int), intent(in)     :: c_self
type(c_ptr), value, intent(in) :: c_vars     !< List of variables
integer(c_int), intent(in)     :: c_geovals
integer(c_int), intent(in)     :: c_nobs
character(c_char), intent(in)  :: c_apply(c_nobs)

type(ufo_rttovonedvarcheck), pointer :: self
type(oops_variables)                 :: vars
type(ufo_geovals), pointer           :: geovals
integer                              :: ii
logical                              :: apply(c_nobs)

call ufo_rttovonedvarcheck_registry%get(c_self, self)
call ufo_geovals_registry%get(c_geovals, geovals)

vars = oops_variables(c_vars)

! Convert character to logical for passing to Fortran
apply(:) = .false.
where (c_apply == 'T')
  apply = .true.
end where

call ufo_rttovonedvarcheck_apply(self, vars, geovals, apply)

end subroutine ufo_rttovonedvarcheck_apply_c

! ------------------------------------------------------------------------------

end module ufo_rttovonedvarcheck_mod_c

