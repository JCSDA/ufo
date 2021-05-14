
! (C) Copyright 2017-2020 Met Office
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 

module ufo_gnssroonedvarcheck_mod_c

use fckit_configuration_module, only: fckit_configuration
use iso_c_binding
use oops_variables_mod
use ufo_geovals_mod
use ufo_geovals_mod_c,   only: ufo_geovals_registry
use ufo_gnssroonedvarcheck_mod

implicit none
private

#define LISTED_TYPE ufo_gnssroonedvarcheck

!> Linked list interface - defines registry_t type
#include "oops/util/linkedList_i.f"

!> Global registry
type(registry_t) :: ufo_gnssroonedvarcheck_registry

! ------------------------------------------------------------------------------------------------
contains
! ------------------------------------------------------------------------------------------------
!> Linked list implementation
#include "oops/util/linkedList_c.f"
! ------------------------------------------------------------------------------------------------

subroutine ufo_gnssroonedvarcheck_create_c(c_self, c_obspace, c_conf, c_onedvarflag) &
                        bind(c,name='ufo_gnssroonedvarcheck_create_f90')

!> \brief Interface to the Fortran create method
!!
!! \author Met Office
!!
!! \date 09/06/2020: Created
!!
implicit none
integer(c_int), intent(inout)  :: c_self     !< self - inout
type(c_ptr), value, intent(in) :: c_obspace  !< obsspace - input
type(c_ptr), value, intent(in) :: c_conf     !< yaml configuration - input
integer(c_int), intent(in) :: c_onedvarflag  !< flag for qc manager logging - input

type(ufo_gnssroonedvarcheck), pointer :: self
type(fckit_configuration) :: f_conf

call ufo_gnssroonedvarcheck_registry%setup(c_self, self)
f_conf = fckit_configuration(c_conf)

call ufo_gnssroonedvarcheck_create(self, c_obspace, f_conf, c_onedvarflag)

end subroutine ufo_gnssroonedvarcheck_create_c

! ------------------------------------------------------------------------------------------------

subroutine ufo_gnssroonedvarcheck_delete_c(c_self) &
                      bind(c,name='ufo_gnssroonedvarcheck_delete_f90')

!> \brief Interface to the Fortran delete method
!!
!! \author Met Office
!!
!! \date 09/06/2020: Created
!!

implicit none
integer(c_int), intent(inout) :: c_self !< self - input

type(ufo_gnssroonedvarcheck), pointer :: self

call ufo_gnssroonedvarcheck_registry%get(c_self, self)
call ufo_gnssroonedvarcheck_delete(self)
call ufo_gnssroonedvarcheck_registry%delete(c_self, self)

end subroutine ufo_gnssroonedvarcheck_delete_c

! ------------------------------------------------------------------------------------------------

subroutine ufo_gnssroonedvarcheck_apply_c(c_self, c_geovals, c_nobs, c_apply) &
               bind(c,name='ufo_gnssroonedvarcheck_apply_f90')

!> \brief Interface to filter apply method
!!
!! \author Met Office
!!
!! \date 09/06/2020: Created
!!

implicit none
integer(c_int), intent(in)     :: c_self          !< self - input
integer(c_int), intent(in)     :: c_geovals       !< Geovals - input
integer(c_int), intent(in)     :: c_nobs          !< number of observations - input
character(c_char), intent(in)  :: c_apply(c_nobs) !< apply flag (converted to logical) - input

type(ufo_gnssroonedvarcheck), pointer :: self
type(oops_variables)                 :: vars
type(oops_variables)                 :: retrieval_vars
type(ufo_geovals), pointer           :: geovals
integer                              :: ii
logical                              :: apply(c_nobs)

call ufo_gnssroonedvarcheck_registry%get(c_self, self)
call ufo_geovals_registry%get(c_geovals, geovals)

! Convert character to logical for passing to Fortran
apply(:) = .false.
where (c_apply == 'T')
  apply = .true.
end where

call ufo_gnssroonedvarcheck_apply(self, geovals, apply)

end subroutine ufo_gnssroonedvarcheck_apply_c

! ------------------------------------------------------------------------------------------------

end module ufo_gnssroonedvarcheck_mod_c

