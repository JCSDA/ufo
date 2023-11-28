
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

subroutine ufo_gnssroonedvarcheck_create_c(c_self, &
                                           c_obspace, &
                                           filename_length, &
                                           input_filename, &
                                           capsupersat, &
                                           cost_funct_test, &
                                           Delta_ct2, &
                                           Delta_factor, &
                                           min_temp_grad, &
                                           n_iteration_test, &
                                           OB_test, &
                                           pseudo_ops, &
                                           vert_interp_ops, &
                                           y_test, &
                                           c_onedvarflag, &
                                           nchans, &
                                           chanList, &
                                           noSuperCheck) &
                        bind(c,name='ufo_gnssroonedvarcheck_create_f90')

!> \brief Interface to the Fortran create method
!!
!! \author Met Office
!!
!! \date 09/06/2020: Created
!!
use string_f_c_mod

implicit none

integer(c_int), intent(inout)             :: c_self            !< self - inout
type(c_ptr), value, intent(in)            :: c_obspace         !< obsspace - input
integer(c_int), intent(in)                :: filename_length   !< Length of the filename string
character(kind=c_char, len=1), intent(in) :: input_filename(filename_length+1)  !< B-matrix filename
logical(c_bool), intent(in)               :: capsupersat       !< Whether to remove super-saturation (wrt ice?)
real(c_float), intent(in)                 :: cost_funct_test   !< Threshold value for the cost function convergence test
real(c_float), intent(in)                 :: Delta_ct2         !< Threshold used in calculating convergence
real(c_float), intent(in)                 :: Delta_factor      !< Threshold used in calculating convergence
real(c_float), intent(in)                 :: min_temp_grad     !< The minimum vertical temperature gradient allowed
integer(c_int), intent(in)                :: n_iteration_test  !< Maximum number of iterations in the 1DVar
real(c_float), intent(in)                 :: OB_test           !< Threshold for the O-B throughout the profile
logical(c_bool), intent(in)               :: pseudo_ops        !< Whether to use pseudo levels in forward operator
logical(c_bool), intent(in)               :: vert_interp_ops   !< Whether to use ln(p) or exner in vertical interpolation
real(c_float), intent(in)                 :: y_test            !< Threshold on distance between observed and solution bending angles
integer(c_int), intent(in)                :: c_onedvarflag     !< flag for qc manager logging - input
integer(c_int), intent(in)                :: nchans            !< Number of channels (levels) to be used
integer(c_int), intent(in)                :: chanList(nchans)  !< List of channels to use
logical(c_bool), intent(in)               :: noSuperCheck      !< Whether to avoid using super-refraction check in operator

character(len=filename_length) :: bmatrix_filename  ! Location of the B-matrix file
integer                        :: ifname            ! Loop variable for filename
integer(c_int), allocatable    :: localChanList(:)  ! Allocated list of channels (even if nchans=0)

type(ufo_gnssroonedvarcheck), pointer :: self

call ufo_gnssroonedvarcheck_registry%setup(c_self, self)

! copy over the char* into a Fortran character
call c_f_string(input_filename, bmatrix_filename)

if (nchans == 0) then
  allocate(localChanList(1))
  localChanList(1) = 0
else
  allocate(localChanList(nchans))
  localChanList(1:nchans) = chanList(1:nchans)
end if

call ufo_gnssroonedvarcheck_create(self, &
                                   c_obspace, &
                                   bmatrix_filename, &
                                   capsupersat, &
                                   cost_funct_test, &
                                   Delta_ct2, &
                                   Delta_factor, &
                                   min_temp_grad, &
                                   n_iteration_test, &
                                   OB_test, &
                                   pseudo_ops, &
                                   vert_interp_ops, &
                                   y_test, &
                                   c_onedvarflag, &
                                   localChanList, &
                                   noSuperCheck)

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
call ufo_gnssroonedvarcheck_registry%remove(c_self)

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

