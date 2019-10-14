! (C) Copyright 2017-2018 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

!> Fortran timeoper module for functions on the interface between C++ and Fortran
!  to handle observation operators

module ufo_timeoper_mod_c

use iso_c_binding
use fckit_configuration_module, only: fckit_configuration
use ufo_locs_mod
use ufo_locs_mod_c
use ufo_timeoper_mod
use ufo_timeoper_locs_mod

use ufo_geovals_mod,   only: ufo_geovals
use ufo_geovals_mod_c, only: ufo_geovals_registry

implicit none

public :: ufo_timeoper_registry

private

  ! ------------------------------------------------------------------------------
#define LISTED_TYPE ufo_timeoper

!> Linked list interface - defines registry_t type
#include "oops/util/linkedList_i.f"

!> Global registry
type(registry_t) :: ufo_timeoper_registry

! ------------------------------------------------------------------------------

contains

  ! ----------------------------------------------------------------------------
  !> Linked list implementation
#include "oops/util/linkedList_c.f"

! ------------------------------------------------------------------------------

subroutine ufo_timeoper_setup_c(c_key_self, c_conf, c_size, &
  c_tstencil) bind(c,name='ufo_timeoper_setup_f90')

implicit none
integer(c_int), intent(inout) :: c_key_self
type(c_ptr),    intent(in)    :: c_conf
integer(c_int), intent(in)    :: c_size  ! obsspace vector length
integer(c_int), intent(in)    :: c_tstencil

type(ufo_timeoper), pointer :: self
type(fckit_configuration) :: f_conf

f_conf = fckit_configuration(c_conf)
call ufo_timeoper_registry%setup(c_key_self, self)
call self%setup(f_conf, c_size, c_tstencil)

end subroutine ufo_timeoper_setup_c


! ------------------------------------------------------------------------------
subroutine ufo_timeoper_set_timeweight_c(c_key_self, c_conf, c_obsspace, &
  c_t0, c_t3, c_st, c_ts) bind(c,name='ufo_timeoper_set_timeweight_f90')
use datetime_mod
use duration_mod
implicit none
integer(c_int),     intent(inout)  :: c_key_self  ! operator key
type(c_ptr),        intent(in)     :: c_conf
type(c_ptr), value, intent(in)     :: c_obsspace
type(c_ptr),        intent(in)     :: c_t0, c_t3
type(c_ptr),        intent(in)     :: c_st        ! stateTime
type(c_ptr),        intent(in)     :: c_ts        ! windowSub

type(ufo_timeoper),        pointer :: self

type(datetime) :: t0, t3, st
type(duration) :: ts
type(fckit_configuration) :: f_conf

f_conf = fckit_configuration(c_conf)
call c_f_datetime(c_t0, t0)
call c_f_datetime(c_t3, t3)
call c_f_datetime(c_st, st)
call c_f_duration(c_ts, ts)

call ufo_timeoper_registry%get(c_key_self, self)
call self%set_timeweight(f_conf, c_obsspace, &
                         t0, t3, st, ts)

end subroutine ufo_timeoper_set_timeweight_c

! ------------------------------------------------------------------------------

subroutine ufo_timeoper_delete_c(c_key_self) &
  bind(c,name='ufo_timeoper_delete_f90')
implicit none
integer(c_int), intent(inout) :: c_key_self

type(ufo_timeoper), pointer :: self

call ufo_timeoper_registry%delete(c_key_self, self)

end subroutine ufo_timeoper_delete_c

! ------------------------------------------------------------------------------

subroutine ufo_timeoper_simobs_c(c_key_self, c_key_geovals, &
  c_obsspace) bind(c,name='ufo_timeoper_simobs_f90')

implicit none
integer(c_int), intent(in) :: c_key_self
integer(c_int), intent(inout) :: c_key_geovals
type(c_ptr), value, intent(in) :: c_obsspace

type(ufo_timeoper), pointer :: self
type(ufo_geovals),    pointer :: geovals

call ufo_timeoper_registry%get(c_key_self, self)
call ufo_geovals_registry%get(c_key_geovals, geovals)
call self%simobs(geovals, c_obsspace)

end subroutine ufo_timeoper_simobs_c

! ------------------------------------------------------------------------------
subroutine ufo_timeoper_locs_init_c(c_key_self, c_key_locs, c_obsspace, &
  c_t0, c_t1, c_t2, c_t3, c_st) bind(c,name='ufo_timeoper_locs_init_f90')
use datetime_mod
use duration_mod
implicit none
integer(c_int),     intent(in)     :: c_key_self  ! operator key
integer(c_int),     intent(inout)  :: c_key_locs  ! location key
type(c_ptr), value, intent(in)     :: c_obsspace
type(c_ptr),        intent(in)     :: c_t0, c_t1, c_t2, c_t3, c_st

type(ufo_locs),            pointer :: locs
type(ufo_timeoper),        pointer :: self

type(datetime) :: t0, t1, t2, t3, st

call c_f_datetime(c_t0, t0)
call c_f_datetime(c_t1, t1)
call c_f_datetime(c_t2, t2)
call c_f_datetime(c_t3, t3)
call c_f_datetime(c_st, st)


call ufo_locs_registry%get(c_key_locs, locs)
call ufo_timeoper_registry%get(c_key_self, self)
call ufo_timeoper_locs_init(self, locs, c_obsspace, &
                            t0, t1, t2, t3, st)

end subroutine ufo_timeoper_locs_init_c


end module ufo_timeoper_mod_c
