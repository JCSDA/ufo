!
! (C) Copyright 2017 UCAR
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
!
module ufo_geovals_mod_c

use iso_c_binding
use ufo_geovals_mod
use ufo_vars_mod
use kinds

implicit none
private
! ------------------------------------------------------------------------------

contains

! ------------------------------------------------------------------------------

subroutine ufo_geovals_create_c(c_key_self) bind(c,name='ufo_geovals_create_f90')

implicit none
integer(c_int), intent(inout) :: c_key_self

type(ufo_geovals), pointer :: self

call ufo_geovals_registry%init()
call ufo_geovals_registry%add(c_key_self)
call ufo_geovals_registry%get(c_key_self, self)

self%lalloc = .false.
self%linit  = .false.

end subroutine ufo_geovals_create_c

! ------------------------------------------------------------------------------

subroutine ufo_geovals_delete_c(c_key_self) bind(c,name='ufo_geovals_delete_f90')

implicit none
integer(c_int), intent(inout) :: c_key_self

type(ufo_geovals), pointer :: self

call ufo_geovals_registry%get(c_key_self, self)

call ufo_geovals_delete(self)

call ufo_geovals_registry%remove(c_key_self)

end subroutine ufo_geovals_delete_c

! ------------------------------------------------------------------------------

subroutine ufo_geovals_zero_c(c_key_self) bind(c,name='ufo_geovals_zero_f90')
implicit none
integer(c_int), intent(in) :: c_key_self
type(ufo_geovals), pointer :: self

call ufo_geovals_registry%get(c_key_self, self)

call ufo_geovals_zero(self)

end subroutine ufo_geovals_zero_c

! ------------------------------------------------------------------------------

subroutine ufo_geovals_random_c(c_key_self) bind(c,name='ufo_geovals_random_f90')
implicit none
integer(c_int), intent(in) :: c_key_self
type(ufo_geovals), pointer :: self

call ufo_geovals_registry%get(c_key_self, self)

call ufo_geovals_random(self)

end subroutine ufo_geovals_random_c

! ------------------------------------------------------------------------------

subroutine ufo_geovals_dotprod_c(c_key_self, c_key_other, prod) bind(c,name='ufo_geovals_dotprod_f90')
implicit none
integer(c_int), intent(in) :: c_key_self, c_key_other
real(c_double), intent(inout) :: prod
type(ufo_geovals), pointer :: self, other

call ufo_geovals_registry%get(c_key_self, self)
call ufo_geovals_registry%get(c_key_other, other)

call ufo_geovals_dotprod(self, other, prod)

end subroutine ufo_geovals_dotprod_c

! ------------------------------------------------------------------------------

subroutine ufo_geovals_minmaxavg_c(c_key_self, kobs, pmin, pmax, prms) bind(c,name='ufo_geovals_minmaxavg_f90')
implicit none
integer(c_int), intent(in) :: c_key_self
integer(c_int), intent(inout) :: kobs
real(c_double), intent(inout) :: pmin, pmax, prms
type(ufo_geovals), pointer :: self

call ufo_geovals_registry%get(c_key_self, self)

call ufo_geovals_minmaxavg(self, kobs, pmin, pmax, prms)

end subroutine ufo_geovals_minmaxavg_c

! ------------------------------------------------------------------------------

subroutine ufo_geovals_read_file_c(c_key_self, c_conf) bind(c,name='ufo_geovals_read_file_f90')
use config_mod

implicit none
integer(c_int), intent(in) :: c_key_self
type(c_ptr), intent(in)    :: c_conf

type(ufo_geovals), pointer :: self

character(128) :: filename

call ufo_geovals_registry%get(c_key_self, self)

! read filename for config
filename = config_get_string(c_conf,len(filename),"filename")

call ufo_geovals_read_netcdf(self, filename)

end subroutine ufo_geovals_read_file_c

! ------------------------------------------------------------------------------

subroutine ufo_geovals_write_file_c(c_key_self, c_conf) bind(c,name='ufo_geovals_write_file_f90')
use config_mod
implicit none
integer(c_int), intent(in) :: c_key_self
type(c_ptr), intent(in) :: c_conf
type(ufo_geovals), pointer :: self

call ufo_geovals_registry%get(c_key_self, self)

end subroutine ufo_geovals_write_file_c

! ------------------------------------------------------------------------------

end module ufo_geovals_mod_c
