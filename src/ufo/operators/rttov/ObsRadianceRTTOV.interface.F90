! (C) Copyright 2017-2018 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

!> Fortran module to handle radiancerttov observations

module ufo_radiancerttov_mod_c

  use fckit_configuration_module, only: fckit_configuration
  use iso_c_binding
  use ufo_radiancerttov_mod 
  use ufo_geovals_mod
  use ufo_geovals_mod_c,   only: ufo_geovals_registry

  implicit none
  private

  ! ------------------------------------------------------------------------------
#define LISTED_TYPE ufo_radiancerttov

  !> Linked list interface - defines registry_t type
#include "oops/util/linkedList_i.f"

  !> Global registry
  type(registry_t) :: ufo_radiancerttov_registry

  ! ------------------------------------------------------------------------------

contains

  ! ------------------------------------------------------------------------------
  !> Linked list implementation
#include "oops/util/linkedList_c.f"

! ------------------------------------------------------------------------------

subroutine ufo_radiancerttov_setup_c(c_key_self, c_conf, c_nchan, c_channels, c_varlist, &
                                    c_qc_size, c_qc_passed) &
                                    bind(c,name='ufo_radiancerttov_setup_f90')
use oops_variables_mod
implicit none
integer(c_int), intent(inout)  :: c_key_self
type(c_ptr), value, intent(in) :: c_conf
integer(c_int), intent(in) :: c_nchan
integer(c_int), intent(in) :: c_channels(c_nchan)
type(c_ptr), intent(in), value :: c_varlist
integer(c_int), intent(in) :: c_qc_size              !< number of flags calssified as passed - in
integer(c_int), intent(in) :: c_qc_passed(c_qc_size) !< flags for passivated obs. - in

type(oops_variables) :: oops_vars
type(ufo_radiancerttov), pointer :: self
type(fckit_configuration) :: f_conf

call ufo_radiancerttov_registry%setup(c_key_self, self)
f_conf = fckit_configuration(c_conf)

call self%setup(f_conf, c_channels, c_qc_passed)

!> Update C++ ObsOperator with input variable list
oops_vars = oops_variables(c_varlist)
call oops_vars%push_back( self%varin )

end subroutine ufo_radiancerttov_setup_c

! ------------------------------------------------------------------------------

subroutine ufo_radiancerttov_delete_c(c_key_self) bind(c,name='ufo_radiancerttov_delete_f90')
implicit none
integer(c_int), intent(inout) :: c_key_self
    
type(ufo_radiancerttov), pointer :: self

call ufo_radiancerttov_registry%get(c_key_self, self)

call self%delete()

call ufo_radiancerttov_registry%remove(c_key_self)

end subroutine ufo_radiancerttov_delete_c

! ------------------------------------------------------------------------------
subroutine ufo_radiancerttov_simobs_c(c_key_self, c_key_geovals, c_obsspace, c_nvars, c_nlocs, &
           c_hofx, c_key_hofxdiags, c_qc_flags) bind(c,name='ufo_radiancerttov_simobs_f90')
implicit none
integer(c_int), intent(in) :: c_key_self
integer(c_int), intent(in) :: c_key_geovals
type(c_ptr), value, intent(in) :: c_obsspace
integer(c_int), intent(in) :: c_nvars, c_nlocs
real(c_double), intent(inout) :: c_hofx(c_nvars, c_nlocs)
integer(c_int), intent(in) :: c_key_hofxdiags
type(c_ptr),value, intent(in) :: c_qc_flags

type(ufo_radiancerttov), pointer :: self
type(ufo_geovals),  pointer :: geovals
type(ufo_geovals),  pointer :: hofxdiags

character(len=*), parameter :: myname_="ufo_radiancerttov_simobs_c"

call ufo_radiancerttov_registry%get(c_key_self, self)

call ufo_geovals_registry%get(c_key_geovals,geovals)
call ufo_geovals_registry%get(c_key_hofxdiags,hofxdiags)

call self%simobs(geovals, c_obsspace, c_nvars, c_nlocs, c_hofx, &
                 hofxdiags, qcf_p=c_qc_flags)

end subroutine ufo_radiancerttov_simobs_c

! ------------------------------------------------------------------------------

end module ufo_radiancerttov_mod_c
