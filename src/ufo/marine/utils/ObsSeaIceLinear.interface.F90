! (C) Copyright 2017-2019 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

!> Fortran seaicelinear module for functions on the interface between C++ and Fortran
!  to handle linearized observation operators

module ufo_seaicelinear_mod_c

  use fckit_configuration_module, only: fckit_configuration
  use iso_c_binding
  use ufo_geovals_mod
  use ufo_geovals_mod_c,   only: ufo_geovals_registry
  use ufo_vars_mod

  implicit none

  private
  integer, parameter :: max_string=800
  type, public :: ufo_seaicelinear
     integer :: ncat = -1      !< number of ice categories
  end type ufo_seaicelinear

#define LISTED_TYPE ufo_seaicelinear

  !> Linked list interface - defines registry_t type
#include "oops/util/linkedList_i.f"

  !> Global registry
  type(registry_t) :: ufo_seaicelinear_registry

contains

  ! ------------------------------------------------------------------------------
  !> Linked list implementation
#include "oops/util/linkedList_c.f"

! ------------------------------------------------------------------------------

subroutine ufo_seaicelinear_setup_c(c_key_self, c_conf) bind(c,name='ufo_seaicelinear_setup_f90')
integer(c_int), intent(inout)  :: c_key_self
type(c_ptr), value, intent(in) :: c_conf

type(ufo_seaicelinear), pointer :: self
type(fckit_configuration) :: f_conf

call ufo_seaicelinear_registry%setup(c_key_self, self)
f_conf = fckit_configuration(c_conf)
!call self%setup(f_conf)

end subroutine ufo_seaicelinear_setup_c

! ------------------------------------------------------------------------------

subroutine ufo_seaicelinear_delete_c(c_key_self) bind(c,name='ufo_seaicelinear_delete_f90')
integer(c_int), intent(inout) :: c_key_self

type(ufo_seaicelinear), pointer :: self

call ufo_seaicelinear_registry%get(c_key_self, self)
call ufo_seaicelinear_registry%remove(c_key_self)

end subroutine ufo_seaicelinear_delete_c

! ------------------------------------------------------------------------------

subroutine ufo_seaicelinear_settraj_c(c_key_self, c_key_geovals, c_obsspace)&
  bind(c,name='ufo_seaicelinear_settraj_f90')
integer(c_int),     intent(in) :: c_key_self
integer(c_int),     intent(in) :: c_key_geovals
type(c_ptr), value, intent(in) :: c_obsspace

type(ufo_seaicelinear), pointer :: self
type(ufo_geovals),    pointer :: geovals
type(ufo_geoval), pointer :: geoval

call ufo_seaicelinear_registry%get(c_key_self, self)
call ufo_geovals_registry%get(c_key_geovals,geovals)

call ufo_geovals_get_var(geovals, var_seaicefrac, geoval)
self%ncat = geoval%nval

end subroutine ufo_seaicelinear_settraj_c

! ------------------------------------------------------------------------------

subroutine ufo_seaicelinear_alloc_ad_c(c_key_self, c_key_geovals, c_obsspace, c_nobs, c_hofx)&
  bind(c,name='ufo_seaicelinear_alloc_ad_f90')
integer(c_int), intent(in) :: c_key_self
integer(c_int), intent(in) :: c_key_geovals
type(c_ptr), value, intent(in) :: c_obsspace
integer(c_int), intent(in) :: c_nobs
real(c_double), intent(in) :: c_hofx(c_nobs)

type(ufo_seaicelinear), pointer :: self
type(ufo_geovals),    pointer :: geovals
type(ufo_geoval), pointer :: geoval
character(max_string) :: err_msg

call ufo_seaicelinear_registry%get(c_key_self, self)
call ufo_geovals_registry%get(c_key_geovals,geovals)

! check if nlocs is consistent in geovals & hofx
if (geovals%nlocs /= size(c_hofx,1)) then
  write(err_msg,*) ' error: nlocs inconsistent!'
  call abor1_ftn(err_msg)
endif

if (.not. geovals%linit ) geovals%linit=.true.

! check if sea ice fraction variables is in geovals and get it
call ufo_geovals_get_var(geovals, var_seaicefrac, geoval)

if (.not.(allocated(geoval%vals))) then
   if (self%ncat < 1) then
     write(err_msg,*)' unknown number of categories'
     call abor1_ftn(err_msg)
   endif
   allocate(geoval%vals(self%ncat,size(c_hofx,1)))
end if
end subroutine ufo_seaicelinear_alloc_ad_c

! ------------------------------------------------------------------------------

end module ufo_seaicelinear_mod_c
