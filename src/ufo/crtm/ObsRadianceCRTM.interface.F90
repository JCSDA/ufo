! (C) Copyright 2017-2018 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

!> Fortran module to handle radiancecrtm observations

module ufo_radiancecrtm_mod_c

  use iso_c_binding
  use config_mod
  use ufo_radiancecrtm_mod
  use ufo_geovals_mod
  use ufo_geovals_mod_c,   only: ufo_geovals_registry

  implicit none
  private

  ! ------------------------------------------------------------------------------
#define LISTED_TYPE ufo_radiancecrtm

  !> Linked list interface - defines registry_t type
#include "../linkedList_i.f"

  !> Global registry
  type(registry_t) :: ufo_radiancecrtm_registry

  ! ------------------------------------------------------------------------------

contains

  ! ------------------------------------------------------------------------------
  !> Linked list implementation
#include "../linkedList_c.f"

! ------------------------------------------------------------------------------

subroutine ufo_radiancecrtm_setup_c(c_key_self, c_conf) bind(c,name='ufo_radiancecrtm_setup_f90')
implicit none
integer(c_int), intent(inout) :: c_key_self
type(c_ptr),    intent(in)    :: c_conf

type(ufo_radiancecrtm), pointer :: self

call ufo_radiancecrtm_registry%setup(c_key_self, self)

call self%setup(c_conf)

end subroutine ufo_radiancecrtm_setup_c

! ------------------------------------------------------------------------------

subroutine ufo_radiancecrtm_delete_c(c_key_self) bind(c,name='ufo_radiancecrtm_delete_f90')
implicit none
integer(c_int), intent(inout) :: c_key_self

type(ufo_radiancecrtm), pointer :: self

call ufo_radiancecrtm_registry%get(c_key_self, self)

call self%delete()

call ufo_radiancecrtm_registry%remove(c_key_self)

end subroutine ufo_radiancecrtm_delete_c

! ------------------------------------------------------------------------------

subroutine ufo_radiancecrtm_simobs_c(c_key_self, c_key_geovals, c_obsspace, c_nobs, c_hofx, c_bias, &
                                 c_nchan, c_channels) bind(c,name='ufo_radiancecrtm_simobs_f90')

implicit none
integer(c_int), intent(in) :: c_key_self
integer(c_int), intent(in) :: c_key_geovals
type(c_ptr), value, intent(in) :: c_obsspace
integer(c_int), intent(in) :: c_nobs
real(c_double), intent(inout) :: c_hofx(c_nobs)
integer(c_int), intent(in) :: c_bias
integer(c_int), intent(in) :: c_nchan
integer(c_int), intent(in) :: c_channels(c_nchan)


type(ufo_radiancecrtm), pointer :: self
type(ufo_geovals),  pointer :: geovals

character(len=*), parameter :: myname_="ufo_radiancecrtm_simobs_c"

call ufo_radiancecrtm_registry%get(c_key_self, self)

call ufo_geovals_registry%get(c_key_geovals,geovals)

call self%simobs(geovals, c_hofx, c_obsspace, c_channels)

end subroutine ufo_radiancecrtm_simobs_c

! ------------------------------------------------------------------------------

end module ufo_radiancecrtm_mod_c
