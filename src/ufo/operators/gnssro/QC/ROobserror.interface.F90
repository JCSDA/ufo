!
! (C) Copyright 2017-2018 UCAR
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
!
module ufo_roobserror_mod_c

use fckit_configuration_module, only: fckit_configuration
use iso_c_binding
use ufo_roobserror_mod
use kinds

use ufo_geovals_mod
use ufo_geovals_mod_c,   only: ufo_geovals_registry

implicit none
private

#define LISTED_TYPE ufo_roobserror

!> Linked list interface - defines registry_t type
#include "oops/util/linkedList_i.f"

!> Global registry
type(registry_t) :: ufo_roobserror_registry

! ------------------------------------------------------------------------------
contains
! ------------------------------------------------------------------------------
!> Linked list implementation
#include "oops/util/linkedList_c.f"
! ------------------------------------------------------------------------------

subroutine ufo_roobserror_create_c(c_self, c_obspace, c_conf, c_filtervar) bind(c,name='ufo_roobserror_create_f90')
use oops_variables_mod
implicit none
integer(c_int), intent(inout)  :: c_self
type(c_ptr), value, intent(in) :: c_obspace
type(c_ptr), value, intent(in) :: c_conf
type(c_ptr), value, intent(in) :: c_filtervar
type(ufo_roobserror), pointer  :: self
type(fckit_configuration)      :: f_conf

call ufo_roobserror_registry%setup(c_self, self)
f_conf = fckit_configuration(c_conf)

self%obsvar   = oops_variables(c_filtervar)
self%variable = self%obsvar%variable(1)

call ufo_roobserror_create(self, c_obspace, f_conf)

end subroutine ufo_roobserror_create_c

! ------------------------------------------------------------------------------

subroutine ufo_roobserror_delete_c(c_self) bind(c,name='ufo_roobserror_delete_f90')
implicit none
integer(c_int), intent(inout) :: c_self

type(ufo_roobserror), pointer :: self

call ufo_roobserror_registry%get(c_self, self)
call ufo_roobserror_delete(self)
call ufo_roobserror_registry%remove(c_self)

end subroutine ufo_roobserror_delete_c

! ------------------------------------------------------------------------------

subroutine ufo_roobserror_prior_c(c_self, air_nobs, air_nlevs, air_temperature, &
                                  height_nobs, height_nlevs, geopotential_height) bind(c,name='ufo_roobserror_prior_f90')

implicit none

integer(c_int), intent(in)    :: c_self      ! The object containing configuration info
integer(c_int), intent(in)   :: air_nlevs    ! The number of vertical levels in air_temperature
integer(c_int), intent(in)   :: air_nobs     ! The number of observations in air_temperature
real(c_float), intent(in)   :: air_temperature(1:air_nobs,1:air_nlevs)
                                             ! The geovals with air temperature
integer(c_int), intent(in)   :: height_nlevs ! The number of vertical levels in geopotential_height
integer(c_int), intent(in)   :: height_nobs  ! The number of observations in geopotential_height
real(c_float), intent(in)   :: geopotential_height(1:height_nobs,1:height_nlevs)
                                             ! The geovals with geopotential height

type(ufo_roobserror), pointer :: self
character(len=200)            :: ErrorMessage   ! Error message to be output

call ufo_roobserror_registry%get(c_self, self)

if (height_nlevs /= air_nlevs) then
  write(ErrorMessage, '(2A,2I7)') "air_temperature and geopotential_height must have the same number of levels", &
    " I received: ", air_nlevs, height_nlevs
  call abor1_ftn(ErrorMessage)
end if

if (height_nobs /= air_nobs) then
  write(ErrorMessage, '(2A,2I7)') "air_temperature and geopotential_height must have the same number of locations", &
    " I received: ", air_nobs, height_nobs
  call abor1_ftn(ErrorMessage)
end if

call ufo_roobserror_prior(self, air_nobs, air_nlevs, air_temperature, &
    geopotential_height)

end subroutine ufo_roobserror_prior_c

! ------------------------------------------------------------------------------

end module ufo_roobserror_mod_c
