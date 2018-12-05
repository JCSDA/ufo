! (C) Copyright 2017-2018 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

!> Fortran module to handle aircraft observations

module ufo_aircraft_mod_c

  use iso_c_binding
  use config_mod
  use ufo_aircraft_mod
  implicit none

  integer, parameter :: max_string=800

  private

#define LISTED_TYPE ufo_aircraft

  !> Linked list interface - defines registry_t type
#include "../../linkedList_i.f"

  !> Global registry
  type(registry_t) :: ufo_aircraft_registry

  ! ------------------------------------------------------------------------------
contains
  ! ------------------------------------------------------------------------------
  !> Linked list implementation
#include "../../linkedList_c.f"

! ------------------------------------------------------------------------------

subroutine ufo_aircraft_setup_c(c_key_self, c_conf) bind(c,name='ufo_aircraft_setup_f90')
implicit none
integer(c_int), intent(inout) :: c_key_self
type(c_ptr), intent(in)    :: c_conf

type(ufo_aircraft), pointer :: self

call ufo_aircraft_registry%setup(c_key_self, self)

call self%setup(c_conf)

end subroutine ufo_aircraft_setup_c

! ------------------------------------------------------------------------------

subroutine ufo_aircraft_delete_c(c_key_self) bind(c,name='ufo_aircraft_delete_f90')
implicit none
integer(c_int), intent(inout) :: c_key_self

type(ufo_aircraft), pointer :: self

call ufo_aircraft_registry%delete(c_key_self, self)

end subroutine ufo_aircraft_delete_c

! ------------------------------------------------------------------------------

subroutine ufo_aircraft_simobs_c(c_key_self, c_key_geovals, c_obsspace, c_nobs, c_hofx, c_bias) bind(c,name='ufo_aircraft_simobs_f90')

implicit none
integer(c_int), intent(in) :: c_key_self
integer(c_int), intent(in) :: c_key_geovals
type(c_ptr), value, intent(in) :: c_obsspace
integer(c_int), intent(in) :: c_nobs
real(c_double), intent(inout) :: c_hofx(c_nobs)
integer(c_int), intent(in) :: c_bias

type(ufo_aircraft), pointer :: self

character(len=*), parameter :: myname_="ufo_aircraft_simobs_c"

call ufo_aircraft_registry%get(c_key_self, self)

call self%opr_simobs(c_key_geovals, c_obsspace, c_hofx)

end subroutine ufo_aircraft_simobs_c

! ------------------------------------------------------------------------------

subroutine ufo_aircraft_locateobs_c(c_key_self, c_obsspace, c_t1, c_t2, c_locs) bind(c,name='ufo_aircraft_locateobs_f90')

implicit none
integer(c_int), intent(in)     :: c_key_self
type(c_ptr), value, intent(in) :: c_obsspace
type(c_ptr), intent(in)        :: c_t1, c_t2
integer(c_int), intent(inout)  :: c_locs

type(ufo_aircraft), pointer :: self

character(len=*), parameter :: myname_="ufo_aircraft_locateobs_c"

call ufo_aircraft_registry%get(c_key_self, self)
call self%opr_locateobs(c_obsspace, c_t1, c_t2, c_locs)

end subroutine ufo_aircraft_locateobs_c

! ------------------------------------------------------------------------------

end module ufo_aircraft_mod_c
