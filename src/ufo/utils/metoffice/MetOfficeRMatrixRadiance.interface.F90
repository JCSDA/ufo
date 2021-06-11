!-------------------------------------------------------------------------------
! (C) Crown Copyright 2021 Met Office
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
!-------------------------------------------------------------------------------

module ufo_metoffice_rmatrixradiance_mod_c

use fckit_configuration_module, only: fckit_configuration
use iso_c_binding
use kinds
use ufo_metoffice_rmatrixradiance_mod

implicit none

private

#define LISTED_TYPE ufo_metoffice_rmatrixradiance

!> Linked list interface - defines registry_t type
#include "oops/util/linkedList_i.f"

!> Global registry
type(registry_t) :: ufo_metoffice_rmatrixradiance_registry

contains

!> Linked list implementation
#include "oops/util/linkedList_c.f"

!-------------------------------------------------------------------------------
subroutine ufo_metoffice_rmatrixradiance_setup_c(c_self, c_conf, nchans, wmoid, rtype) &
           bind(c, name='ufo_metoffice_rmatrixradiance_setup_f90')

implicit none
integer(c_int), intent(inout)            :: c_self
type(c_ptr), value, intent(in)           :: c_conf
integer(c_size_t), intent(inout)         :: nchans
integer(c_size_t), intent(inout)         :: wmoid
integer(c_size_t), intent(inout)         :: rtype

type(ufo_metoffice_rmatrixradiance), pointer :: self
type(fckit_configuration)                  :: f_conf
character(len=:), allocatable              :: str
character(len=200)                         :: filepath

! Interface and setup
call ufo_metoffice_rmatrixradiance_registry % setup(c_self, self)

! Get filepath from configuration
f_conf = fckit_configuration(c_conf)
call f_conf % get_or_die("RMatrix", str)
filepath = str

! Call Fortran
call self % setup(trim(filepath))

! Pass back dimension and wmoid
nchans = self % nchans
wmoid = self % wmo_id
rtype = self % rtype

end subroutine ufo_metoffice_rmatrixradiance_setup_c

!-------------------------------------------------------------------------------
subroutine ufo_metoffice_rmatrixradiance_delete_c(c_self) &
           bind(c, name='ufo_metoffice_rmatrixradiance_delete_f90')

implicit none
integer(c_int), intent(inout) :: c_self

! Interface and setup
type(ufo_metoffice_rmatrixradiance), pointer :: self
call ufo_metoffice_rmatrixradiance_registry % get(c_self, self)

! Delete
call self % delete()
call ufo_metoffice_rmatrixradiance_registry % delete(c_self, self)

end subroutine ufo_metoffice_rmatrixradiance_delete_c

!-------------------------------------------------------------------------------
subroutine ufo_metoffice_rmatrixradiance_getelements_c(c_self, nchans, channels, obs_error) &
                                bind(C, name='ufo_metoffice_rmatrixradiance_getelements_f90')

implicit none
integer(c_int), intent(inout)    :: c_self
integer(c_size_t), intent(in)    :: nchans
integer(c_int), intent(inout)    :: channels(nchans)
real(c_float), intent(inout)     :: obs_error(nchans)

type(ufo_metoffice_rmatrixradiance), pointer :: self

call ufo_metoffice_rmatrixradiance_registry%get(c_self, self)

channels = self % channels(1:nchans)
obs_error = self % errors(1:nchans)

end subroutine ufo_metoffice_rmatrixradiance_getelements_c

!-------------------------------------------------------------------------------
subroutine ufo_metoffice_rmatrixradiance_print_c(c_self) &
           bind(c, name='ufo_metoffice_rmatrixradiance_print_f90')

implicit none
integer(c_int), intent(inout) :: c_self

! Interface and setup
type(ufo_metoffice_rmatrixradiance), pointer :: self
call ufo_metoffice_rmatrixradiance_registry % get(c_self, self)

! Print information about object
call self % print()

end subroutine ufo_metoffice_rmatrixradiance_print_c

!-------------------------------------------------------------------------------
end module ufo_metoffice_rmatrixradiance_mod_c
