! (C) Copyright 2018 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

!> Fortran example module for for functions on the interface between C++ and Fortran
!  to handle observation space

! TODO: replace example with your_observation_space_name through the file

module ufo_obs_example_mod_c

use iso_c_binding
use string_f_c_mod
use config_mod
use datetime_mod
use duration_mod
use ufo_geovals_mod
use ufo_geovals_mod_c, only : ufo_geovals_registry
use ufo_locs_mod
use ufo_locs_mod_c, only : ufo_locs_registry
use ufo_obs_vectors
use ufo_vars_mod
use ufo_obs_example_mod
use fckit_log_module, only : fckit_log
use kinds

implicit none
private

public :: ufo_obs_example_registry

! ------------------------------------------------------------------------------
integer, parameter :: max_string=800
! ------------------------------------------------------------------------------

#define LISTED_TYPE ufo_obs_example

!> Linked list interface - defines registry_t type
#include "../../linkedList_i.f"

!> Global registry
type(registry_t) :: ufo_obs_example_registry

! ------------------------------------------------------------------------------
contains
! ------------------------------------------------------------------------------
!> Linked list implementation
#include "../../linkedList_c.f"

! ------------------------------------------------------------------------------

subroutine ufo_obsdb_example_setup_c(c_key_self, c_conf) bind(c,name='ufo_obsdb_example_setup_f90')
implicit none
integer(c_int), intent(inout) :: c_key_self
type(c_ptr), intent(in)       :: c_conf !< configuration

type(ufo_obs_example), pointer :: self
character(len=max_string)  :: fin
character(len=max_string)  :: MyObsType
character(len=255) :: record

if (config_element_exists(c_conf,"ObsData.ObsDataIn")) then
  fin  = config_get_string(c_conf,max_string,"ObsData.ObsDataIn.obsfile")
else
  fin  = ""
endif
call fckit_log%info(record)

call ufo_obs_example_registry%init()
call ufo_obs_example_registry%add(c_key_self)
call ufo_obs_example_registry%get(c_key_self, self)
if (trim(fin) /= "") then
! TODO: replace with the call to your Fortran routine for reading observations
!       (defined in ufo_obs_<your_obs_space_name>_mod.F90)
  call ufo_obs_example_read(fin, self)
endif

end subroutine ufo_obsdb_example_setup_c

! ------------------------------------------------------------------------------

subroutine ufo_obsdb_example_getlocations_c(c_key_self, c_t1, c_t2, c_key_locs) bind(c,name='ufo_obsdb_example_getlocations_f90')
implicit none
integer(c_int), intent(in)    :: c_key_self
type(c_ptr), intent(in)       :: c_t1, c_t2
integer(c_int), intent(inout) :: c_key_locs

type(ufo_obs_example), pointer :: self
type(datetime) :: t1, t2
type(ufo_locs), pointer :: locs

call ufo_obs_example_registry%get(c_key_self, self)
call c_f_datetime(c_t1, t1)
call c_f_datetime(c_t2, t2)

call ufo_locs_registry%init()
call ufo_locs_registry%add(c_key_locs)
call ufo_locs_registry%get(c_key_locs,locs)

! TODO: replace with the call to your Fortran routine for getting locations of obs
!       (defined in ufo_obs_<your_obs_space_name>_mod.F90)
call ufo_obs_example_getlocs(self, locs)

end subroutine ufo_obsdb_example_getlocations_c

! ------------------------------------------------------------------------------

subroutine ufo_obsdb_example_generate_c(c_key_self, c_conf, c_t1, c_t2) bind(c,name='ufo_obsdb_example_generate_f90')
implicit none
integer(c_int), intent(inout) :: c_key_self
type(c_ptr), intent(in)       :: c_conf !< configuration
type(c_ptr), intent(in)       :: c_t1, c_t2

type(ufo_obs_example), pointer :: self
type(datetime) :: t1, t2
integer :: nobs

call ufo_obs_example_registry%get(c_key_self, self)
call c_f_datetime(c_t1, t1)
call c_f_datetime(c_t2, t2)

nobs = config_get_int(c_conf, "nobs")

! TODO: replace with the call to your Fortran routine for generating random observations
!       (defined in ufo_obs_<your_obs_space_name>_mod.F90)
call ufo_obs_example_generate(self, nobs)

end subroutine ufo_obsdb_example_generate_c

! ------------------------------------------------------------------------------

subroutine ufo_obsdb_example_nobs_c(c_key_self, kobs) bind(c,name='ufo_obsdb_example_nobs_f90')
implicit none
integer(c_int), intent(in) :: c_key_self
integer(c_int), intent(inout) :: kobs
type(ufo_obs_example), pointer :: self

call ufo_obs_example_registry%get(c_key_self, self)
kobs = self%nobs

end subroutine ufo_obsdb_example_nobs_c

! ------------------------------------------------------------------------------

subroutine ufo_obsdb_example_delete_c(c_key_self) bind(c,name='ufo_obsdb_example_delete_f90')
implicit none
integer(c_int), intent(inout) :: c_key_self
type(ufo_obs_example), pointer :: self

call ufo_obs_example_registry%get(c_key_self, self)

! TODO: replace with the call to your Fortran routine to destruct the obsspace
!       (defined in ufo_obs_<your_obs_space_name>_mod.F90)
call ufo_obs_example_delete(self)
call ufo_obs_example_registry%remove(c_key_self)

end subroutine ufo_obsdb_example_delete_c

! ------------------------------------------------------------------------------

end module ufo_obs_example_mod_c
