! (C) Copyright 2017-2018 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

!> Fortran module for Fortran-C++ interface functions for CRTM tl/ad obs operator

module ufo_radiancecrtm_tlad_mod_c

  use fckit_configuration_module, only: fckit_configuration
  use fckit_mpi_module,   only: fckit_mpi_comm
  use iso_c_binding
  use ufo_radiancecrtm_tlad_mod
  use ufo_geovals_mod
  use ufo_geovals_mod_c,   only: ufo_geovals_registry

  implicit none
  private

#define LISTED_TYPE ufo_radiancecrtm_tlad

  !> Linked list interface - defines registry_t type
#include "oops/util/linkedList_i.f"

  !> Global registry
  type(registry_t) :: ufo_radiancecrtm_tlad_registry

contains

  ! ------------------------------------------------------------------------------
  !> Linked list implementation
#include "oops/util/linkedList_c.f"

! ------------------------------------------------------------------------------

subroutine ufo_radiancecrtm_tlad_setup_c(c_key_self, c_conf, c_nchan, c_channels, c_varlist, &
                                         c_comm) &
                                    bind(c,name='ufo_radiancecrtm_tlad_setup_f90')
use oops_variables_mod
implicit none
integer(c_int), intent(inout)  :: c_key_self
type(c_ptr), value, intent(in) :: c_conf
integer(c_int), intent(in) :: c_nchan
integer(c_int), intent(in) :: c_channels(c_nchan)
type(c_ptr), intent(in), value :: c_varlist
type(c_ptr), value, intent(in) :: c_comm

type(oops_variables) :: oops_vars
type(ufo_radiancecrtm_tlad), pointer :: self
type(fckit_configuration) :: f_conf
type(fckit_mpi_comm) :: f_comm

call ufo_radiancecrtm_tlad_registry%setup(c_key_self, self)
f_conf = fckit_configuration(c_conf)

f_comm = fckit_mpi_comm(c_comm)

call self%setup(f_conf, c_channels, f_comm)

!> Update C++ ObsOperator with input variable list
oops_vars = oops_variables(c_varlist)
call oops_vars%push_back( self%varin )

end subroutine ufo_radiancecrtm_tlad_setup_c

! ------------------------------------------------------------------------------

subroutine ufo_radiancecrtm_tlad_delete_c(c_key_self) bind(c,name='ufo_radiancecrtm_tlad_delete_f90')
implicit none
integer(c_int), intent(inout) :: c_key_self

type(ufo_radiancecrtm_tlad), pointer :: self

call ufo_radiancecrtm_tlad_registry%get(c_key_self, self)
call self%delete()
call ufo_radiancecrtm_tlad_registry%remove(c_key_self)

end subroutine ufo_radiancecrtm_tlad_delete_c

! ------------------------------------------------------------------------------

subroutine ufo_radiancecrtm_tlad_settraj_c(c_key_self, c_key_geovals, c_obsspace, c_key_hofxdiags) &
                                       bind(c,name='ufo_radiancecrtm_tlad_settraj_f90')

implicit none
integer(c_int),     intent(in) :: c_key_self
integer(c_int),     intent(in) :: c_key_geovals
type(c_ptr), value, intent(in) :: c_obsspace
integer(c_int),     intent(in) :: c_key_hofxdiags

type(ufo_radiancecrtm_tlad), pointer :: self
type(ufo_geovals),       pointer :: geovals
type(ufo_geovals),       pointer :: hofxdiags

character(len=*), parameter :: myname_="ufo_radiancecrtm_tlad_settraj_c"

call ufo_radiancecrtm_tlad_registry%get(c_key_self, self)
call ufo_geovals_registry%get(c_key_geovals,geovals)
call ufo_geovals_registry%get(c_key_hofxdiags,hofxdiags)

call self%settraj(geovals, c_obsspace, hofxdiags)

end subroutine ufo_radiancecrtm_tlad_settraj_c

! ------------------------------------------------------------------------------

subroutine ufo_radiancecrtm_simobs_tl_c(c_key_self, c_key_geovals, c_obsspace, c_nvars, c_nlocs, &
	c_hofx, c_qc_flags) bind(c,name='ufo_radiancecrtm_simobs_tl_f90')
implicit none
integer(c_int),     intent(in)    :: c_key_self
integer(c_int),     intent(in)    :: c_key_geovals
type(c_ptr), value, intent(in)    :: c_obsspace
integer(c_int),     intent(in)    :: c_nvars, c_nlocs
real(c_double),     intent(inout) :: c_hofx(c_nvars, c_nlocs)
type(c_ptr), value, intent(in) :: c_qc_flags
type(ufo_radiancecrtm_tlad), pointer :: self
type(ufo_geovals),       pointer :: geovals

character(len=*), parameter :: myname_="ufo_radiancecrtm_simobs_tl_c"

call ufo_radiancecrtm_tlad_registry%get(c_key_self, self)
call ufo_geovals_registry%get(c_key_geovals,geovals)

call self%simobs_tl(geovals, c_obsspace, c_nvars, c_nlocs, c_hofx, c_qc_flags)

end subroutine ufo_radiancecrtm_simobs_tl_c

! ------------------------------------------------------------------------------

subroutine ufo_radiancecrtm_simobs_ad_c(c_key_self, c_key_geovals, c_obsspace, c_nvars, c_nlocs, &
	c_hofx, c_qc_flags) bind(c,name='ufo_radiancecrtm_simobs_ad_f90')

implicit none
integer(c_int),     intent(in) :: c_key_self
integer(c_int),     intent(in) :: c_key_geovals
type(c_ptr), value, intent(in) :: c_obsspace
integer(c_int),     intent(in) :: c_nvars, c_nlocs
real(c_double),     intent(in) :: c_hofx(c_nvars, c_nlocs)
type(c_ptr), value, intent(in) :: c_qc_flags
type(ufo_radiancecrtm_tlad), pointer :: self
type(ufo_geovals),       pointer :: geovals

character(len=*), parameter :: myname_="ufo_radiancecrtm_simobs_ad_c"

call ufo_radiancecrtm_tlad_registry%get(c_key_self, self)
call ufo_geovals_registry%get(c_key_geovals,geovals)

call self%simobs_ad(geovals, c_obsspace, c_nvars, c_nlocs, c_hofx, c_qc_flags)

end subroutine ufo_radiancecrtm_simobs_ad_c

! ------------------------------------------------------------------------------

end module ufo_radiancecrtm_tlad_mod_c
