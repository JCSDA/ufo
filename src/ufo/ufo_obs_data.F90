! (C) Copyright 2009-2016 ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
! In applying this licence, ECMWF does not waive the privileges and immunities 
! granted to it by virtue of its status as an intergovernmental organisation nor
! does it submit to any jurisdiction.

!> Handle observations for the QG model

module ufo_obs_data

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
use ufo_obs_data_mod
use fckit_log_module, only : fckit_log
use kinds

implicit none
private

public :: obs_setup, obs_delete, obs_get, obs_put, max_string
public :: ufo_obs_data_registry

! ------------------------------------------------------------------------------
integer, parameter :: max_string=800
! ------------------------------------------------------------------------------

#define LISTED_TYPE obs_data

!> Linked list interface - defines registry_t type
#include "linkedList_i.f"

!> Global registry
type(registry_t) :: ufo_obs_data_registry

! ------------------------------------------------------------------------------
contains
! ------------------------------------------------------------------------------
!> Linked list implementation
#include "linkedList_c.f"

! ------------------------------------------------------------------------------

subroutine obs_setup(fin, fout, obtype, self)
implicit none
type(obs_data), intent(inout) :: self
character(len=*), intent(in) :: fin, fout
character(len=*), intent(in) :: obtype

self%filein =fin
self%fileout=fout

call self%SetupBasis(self%filein,obtype,self%nobs,self%nlocs)
call fckit_log%debug("TRACE: ufo_obs_data:obs_setup: done")

end subroutine obs_setup

! ------------------------------------------------------------------------------

subroutine obs_delete(self)
implicit none
type(obs_data), intent(inout) :: self

if (self%fileout/="") call obs_write(self)
call self%Delete()

end subroutine obs_delete

! ------------------------------------------------------------------------------

subroutine obs_get(self, col, ovec)
implicit none
type(obs_data), intent(in) :: self
character(len=*), intent(in) :: col
type(obs_vector), intent(inout) :: ovec

end subroutine obs_get

! ------------------------------------------------------------------------------

subroutine obs_put(self, col, ovec)
implicit none
type(obs_data), intent(inout) :: self
character(len=*), intent(in) :: col
type(obs_vector), intent(in) :: ovec

end subroutine obs_put

! ------------------------------------------------------------------------------

subroutine ufo_obsdb_getlocations_c(c_key_self, c_t1, c_t2, c_key_locs) bind(c,name='ufo_obsdb_getlocations_f90')
implicit none
integer(c_int), intent(in) :: c_key_self
type(c_ptr), intent(in) :: c_t1, c_t2
integer(c_int), intent(inout) :: c_key_locs

type(obs_data), pointer :: self
type(datetime) :: t1, t2
type(ufo_locs), pointer :: locs

call ufo_obs_data_registry%get(c_key_self, self)
call c_f_datetime(c_t1, t1)
call c_f_datetime(c_t2, t2)

call ufo_locs_registry%init()
call ufo_locs_registry%add(c_key_locs)
call ufo_locs_registry%get(c_key_locs,locs)

call self%Obspoint%GetLocs(self%nlocs, locs)
!print *, 'getlocs: ', locs%nlocs, locs%lat(1:10)

end subroutine ufo_obsdb_getlocations_c

! ------------------------------------------------------------------------------

subroutine ufo_obsdb_getgeovals_c(c_key_self, c_vars, c_t1, c_t2, c_key_geovals) bind(c,name='ufo_obsdb_getgeovals_f90')
implicit none
integer(c_int), intent(in) :: c_key_self
type(c_ptr), intent(in)    :: c_vars
type(c_ptr), intent(in) :: c_t1, c_t2
integer(c_int), intent(inout) :: c_key_geovals

type(obs_data), pointer :: self
type(datetime) :: t1, t2
type(ufo_geovals), pointer :: geovals
integer :: nvar
character(len=max_string) :: svars
character(len=MAXVARLEN), allocatable :: cvars(:)
type(ufo_vars) :: vars

call ufo_obs_data_registry%get(c_key_self, self)
call c_f_datetime(c_t1, t1)
call c_f_datetime(c_t2, t2)

nvar = config_get_int(c_vars, "nvars")
allocate(cvars(nvar))
svars = config_get_string(c_vars,len(svars),"variables")
read(svars,*) cvars
call ufo_vars_setup(vars, cvars)

call ufo_geovals_registry%init()
call ufo_geovals_registry%add(c_key_geovals)
call ufo_geovals_registry%get(c_key_geovals,geovals)

call ufo_geovals_init(geovals)
call ufo_geovals_setup(geovals, vars, self%nlocs)

end subroutine ufo_obsdb_getgeovals_c

! ------------------------------------------------------------------------------

subroutine ufo_obsdb_nobs_c(c_key_self, kobs) bind(c,name='ufo_obsdb_nobs_f90')
implicit none
integer(c_int), intent(in) :: c_key_self
integer(c_int), intent(inout) :: kobs
type(obs_data), pointer :: self
call ufo_obs_data_registry%get(c_key_self, self)
kobs = self%nobs
end subroutine ufo_obsdb_nobs_c

! ------------------------------------------------------------------------------

subroutine obs_time_get(self, col, t1, t2, ovec)
implicit none
type(obs_data), intent(in)    :: self
character(len=*), intent(in)  :: col
type(datetime), intent(in)    :: t1, t2
type(obs_vector), intent(inout) :: ovec

end subroutine obs_time_get

! ------------------------------------------------------------------------------
!  Private
! ------------------------------------------------------------------------------

subroutine obs_write(self)
implicit none
type(obs_data), intent(in) :: self
integer :: iout, icol, jc, jo
real(kind=kind_real), allocatable :: ztmp(:)
character(len=20) :: stime

end subroutine obs_write

! ------------------------------------------------------------------------------

subroutine ufo_obsdb_setup_c(c_key_self, c_conf) bind(c,name='ufo_obsdb_setup_f90')
implicit none
integer(c_int), intent(inout) :: c_key_self
type(c_ptr), intent(in)    :: c_conf !< configuration

type(obs_data), pointer :: self
character(len=max_string) :: fin, fout
character(len=max_string+30) :: record
character(len=max_string) :: MyObsType

if (config_element_exists(c_conf,"ObsData.ObsDataIn")) then
  fin  = config_get_string(c_conf,max_string,"ObsData.ObsDataIn.obsfile")
else
  fin  = ""
endif
MyObsType = trim(config_get_string(c_conf,max_string,"ObsType"))
write(record,*)'ufo_obsdb_setup_c: file in =',trim(fin)
call fckit_log%info(record)

!fout = config_get_string(c_conf,max_string,"ObsData.ObsDataOut.obsfile")
!write(record,*)'ufo_obsdb_setup_c: file out=',trim(fout)
!call fckit_log%info(record)
fout = ""

call ufo_obs_data_registry%init()
call ufo_obs_data_registry%add(c_key_self)
call ufo_obs_data_registry%get(c_key_self, self)
call obs_setup(trim(fin), trim(fout), MyObsType, self)

end subroutine ufo_obsdb_setup_c

! ------------------------------------------------------------------------------

subroutine ufo_obsdb_delete_c(c_key_self) bind(c,name='ufo_obsdb_delete_f90')
implicit none
integer(c_int), intent(inout) :: c_key_self
type(obs_data), pointer :: self

call ufo_obs_data_registry%get(c_key_self, self)
call obs_delete(self)
call ufo_obs_data_registry%remove(c_key_self)

end subroutine ufo_obsdb_delete_c

! ------------------------------------------------------------------------------

subroutine ufo_obsdb_get_c(c_key_self, lcol, c_col, c_key_ovec) bind(c,name='ufo_obsdb_get_f90')
implicit none
integer(c_int), intent(in) :: c_key_self
integer(c_int), intent(in) :: lcol
character(kind=c_char,len=1), intent(in) :: c_col(lcol+1)
integer(c_int), intent(in) :: c_key_ovec

type(obs_data), pointer :: self
type(obs_vector), pointer :: ovec
character(len=lcol) :: col

call ufo_obs_data_registry%get(c_key_self, self)
call ufo_obs_vect_registry%get(c_key_ovec,ovec)
call c_f_string(c_col, col)

call obs_get(self, trim(col), ovec)

end subroutine ufo_obsdb_get_c

! ------------------------------------------------------------------------------

subroutine ufo_obsdb_put_c(c_key_self, lcol, c_col, c_key_ovec) bind(c,name='ufo_obsdb_put_f90')
implicit none
integer(c_int), intent(in) :: c_key_self
integer(c_int), intent(in) :: lcol
character(kind=c_char,len=1), intent(in) :: c_col(lcol+1)
integer(c_int), intent(in) :: c_key_ovec

type(obs_data), pointer :: self
type(obs_vector), pointer :: ovec
character(len=lcol) :: col

call ufo_obs_data_registry%get(c_key_self, self)
call ufo_obs_vect_registry%get(c_key_ovec,ovec)
call c_f_string(c_col, col)

call obs_put(self, trim(col), ovec)

end subroutine ufo_obsdb_put_c

! ------------------------------------------------------------------------------

end module ufo_obs_data
