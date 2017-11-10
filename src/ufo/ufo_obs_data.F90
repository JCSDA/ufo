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
use ufo_locs_mod
use ufo_obs_vectors
use ufo_vars_mod
use fckit_log_module, only : fckit_log
use kinds

implicit none
private

public :: obs_data, obs_setup, obs_delete, obs_get, obs_put, obs_count, max_string
public :: obs_data_registry

! ------------------------------------------------------------------------------
integer, parameter :: max_string=800
! ------------------------------------------------------------------------------

!> A type to represent observation data
type obs_data
  integer :: ngrp
  character(len=max_string) :: filein, fileout
  type(group_data), pointer :: grphead => null()
end type obs_data

#define LISTED_TYPE obs_data

!> Linked list interface - defines registry_t type
#include "linkedList_i.f"

!> Global registry
type(registry_t) :: obs_data_registry

! ------------------------------------------------------------------------------

!> A type to represent a linked list of observation group data
type group_data
  character(len=50) :: grpname
  type(group_data), pointer :: next => null()
  integer :: nobs
  integer, allocatable :: seqnos(:)
  type(datetime), allocatable :: times(:)
  type(column_data), pointer :: colhead => null()
end type group_data

! ------------------------------------------------------------------------------

!> A type to represent a linked list of observation columns
type column_data
  character(len=50) :: colname
  type(column_data), pointer :: next => null()
  integer :: ncol
  real(kind=kind_real), allocatable :: values(:,:)
end type column_data

! ------------------------------------------------------------------------------

!> Fortran generic
interface obs_count
  module procedure obs_count_time, obs_count_all, obs_count_indx
end interface obs_count

! ------------------------------------------------------------------------------
contains
! ------------------------------------------------------------------------------
!> Linked list implementation
#include "linkedList_c.f"

! ------------------------------------------------------------------------------

subroutine obs_setup(fin, fout, self)
implicit none
type(obs_data), intent(inout) :: self
character(len=*), intent(in) :: fin, fout

self%ngrp=0
self%filein =fin
self%fileout=fout

if (self%filein/="") call obs_read(self)
call fckit_log%debug("TRACE: ufo_obs_data:obs_setup: done")

end subroutine obs_setup

! ------------------------------------------------------------------------------

subroutine obs_delete(self)
implicit none
type(obs_data), intent(inout) :: self
type(group_data), pointer :: jgrp
type(column_data), pointer :: jcol
integer :: jo

if (self%fileout/="") call obs_write(self)

do while (associated(self%grphead))
  jgrp=>self%grphead
  self%grphead=>jgrp%next
  do jo=1,jgrp%nobs
    call datetime_delete(jgrp%times(jo))
  enddo
  deallocate(jgrp%times)
  do while (associated(jgrp%colhead))
    jcol=>jgrp%colhead
    jgrp%colhead=>jcol%next
    deallocate(jcol%values)
    deallocate(jcol)
  enddo
  deallocate(jgrp)
enddo

end subroutine obs_delete

! ------------------------------------------------------------------------------

subroutine obs_get(self, req, col, ovec)
implicit none
type(obs_data), intent(in) :: self
character(len=*), intent(in) :: req, col
type(obs_vect), intent(inout) :: ovec

end subroutine obs_get

! ------------------------------------------------------------------------------

subroutine obs_put(self, req, col, ovec)
implicit none
type(obs_data), intent(inout) :: self
character(len=*), intent(in) :: req, col
type(obs_vect), intent(in) :: ovec

end subroutine obs_put

! ------------------------------------------------------------------------------

subroutine ufo_obsdb_locations_c(c_key_self, lreq, c_req, c_t1, c_t2, c_key_locs) bind(c,name='ufo_obsdb_locations_f90')
implicit none
integer(c_int), intent(in) :: c_key_self
integer(c_int), intent(in) :: lreq
character(kind=c_char,len=1), intent(in) :: c_req(lreq+1)
type(c_ptr), intent(in) :: c_t1, c_t2
integer(c_int), intent(inout) :: c_key_locs

type(obs_data), pointer :: self
character(len=lreq) :: req
type(datetime) :: t1, t2
type(ufo_locs), pointer :: locs
type(obs_vect) :: ovec
character(len=8) :: col="Location"

call obs_data_registry%get(c_key_self, self)
call c_f_string(c_req, req)
call c_f_datetime(c_t1, t1)
call c_f_datetime(c_t2, t2)

call obs_time_get(self, req, col, t1, t2, ovec)

call ufo_locs_registry%init()
call ufo_locs_registry%add(c_key_locs)
call ufo_locs_registry%get(c_key_locs,locs)
     
call ufo_loc_setup(locs, ovec)

deallocate(ovec%values)

end subroutine ufo_obsdb_locations_c

! ------------------------------------------------------------------------------

subroutine ufo_obsdb_getgeovals_c(c_key_self, c_key_vars, c_t1, c_t2, c_key_geovals) bind(c,name='ufo_obsdb_getgeovals_f90')
implicit none
integer(c_int), intent(in) :: c_key_self
integer(c_int), intent(in) :: c_key_vars
type(c_ptr), intent(in) :: c_t1, c_t2
integer(c_int), intent(inout) :: c_key_geovals

type(obs_data), pointer :: self
type(ufo_vars), pointer :: vars
type(datetime) :: t1, t2
type(ufo_geovals), pointer :: geovals

integer :: nobs

call obs_data_registry%get(c_key_self, self)
!call c_f_string(c_req, req)
call ufo_vars_registry%get(c_key_vars, vars)
call c_f_datetime(c_t1, t1)
call c_f_datetime(c_t2, t2)

!call obs_count(self, req, t1, t2, nobs)
nobs=1
!call obs_count(self, req, t1, t2, mobs)

allocate(geovals)
call ufo_geovals_registry%init()
call ufo_geovals_registry%add(c_key_geovals)
call ufo_geovals_registry%get(c_key_geovals,geovals)

geovals%lalloc = .false. ! very bad! should just call init that adds to registry 
geovals%linit  = .false. ! and initalizes!!!

call ufo_geovals_setup(geovals, vars, nobs)

!deallocate(mobs)

end subroutine ufo_obsdb_getgeovals_c

! ------------------------------------------------------------------------------

subroutine ufo_obsdb_generate_c(c_key_self, lreq, c_req, c_conf, c_bgn, c_step, ktimes, kobs) bind(c,name='ufo_obsdb_generate_f90')
implicit none
integer(c_int), intent(in) :: c_key_self
integer(c_int), intent(in) :: lreq
character(kind=c_char,len=1), intent(in) :: c_req(lreq+1)
type(c_ptr), intent(in)    :: c_conf
type(c_ptr), intent(in)    :: c_bgn
type(c_ptr), intent(in)    :: c_step
integer(c_int), intent(in)  :: ktimes
integer(c_int), intent(inout) :: kobs

end subroutine ufo_obsdb_generate_c

! ------------------------------------------------------------------------------

subroutine ufo_obsdb_nobs_c(c_key_self, lreq, c_req, kobs) bind(c,name='ufo_obsdb_nobs_f90')
implicit none
integer(c_int), intent(in) :: c_key_self
integer(c_int), intent(in) :: lreq
character(kind=c_char,len=1), intent(in) :: c_req(lreq+1)
integer(c_int), intent(inout) :: kobs

end subroutine ufo_obsdb_nobs_c

! ------------------------------------------------------------------------------

subroutine obs_time_get(self, req, col, t1, t2, ovec)
implicit none
type(obs_data), intent(in)    :: self
character(len=*), intent(in)  :: req, col
type(datetime), intent(in)    :: t1, t2
type(obs_vect), intent(inout) :: ovec

end subroutine obs_time_get

! ------------------------------------------------------------------------------

subroutine obs_count_time(self, req, t1, t2, kobs)
implicit none
type(obs_data), intent(in)   :: self
character(len=*), intent(in) :: req
type(datetime), intent(in)   :: t1, t2
integer, intent(inout)       :: kobs

end subroutine obs_count_time

! ------------------------------------------------------------------------------

subroutine obs_count_indx(self, req, t1, t2, kobs)
implicit none
type(obs_data), intent(in)   :: self
character(len=*), intent(in) :: req
type(datetime), intent(in)   :: t1, t2
integer, intent(inout)       :: kobs(:)

end subroutine obs_count_indx

! ------------------------------------------------------------------------------

subroutine obs_count_all(self, req, kobs)
implicit none
type(obs_data), intent(in) :: self
character(len=*), intent(in) :: req
integer, intent(inout) :: kobs

end subroutine obs_count_all

! ------------------------------------------------------------------------------

subroutine obs_create(self, req, times, locs)
implicit none
type(obs_data), intent(inout) :: self
character(len=*), intent(in) :: req
type(datetime), intent(in) :: times(:)
type(obs_vect), intent(in) :: locs

end subroutine obs_create

! ------------------------------------------------------------------------------
!  Private
! ------------------------------------------------------------------------------

subroutine obs_read(self)
implicit none
type(obs_data), intent(inout) :: self
integer :: iin, icol, jo, jc, jg, ncol
type(group_data), pointer :: jgrp
type(column_data), pointer :: jcol
real(kind=kind_real), allocatable :: ztmp(:)
character(len=20) :: stime
character(len=max_string+50) :: record

end subroutine obs_read

! ------------------------------------------------------------------------------

subroutine obs_write(self)
implicit none
type(obs_data), intent(in) :: self
integer :: iout, icol, jc, jo
type(group_data), pointer :: jgrp
type(column_data), pointer :: jcol
real(kind=kind_real), allocatable :: ztmp(:)
character(len=20) :: stime

end subroutine obs_write

! ------------------------------------------------------------------------------

subroutine findgroup(self,req,find)
type(obs_data), intent(in) :: self
character(len=*), intent(in) :: req
type(group_data), pointer, intent(inout) :: find

find=>self%grphead
do while (associated(find))
  if (find%grpname==req) exit
  find=>find%next
enddo

end subroutine findgroup

! ------------------------------------------------------------------------------

subroutine findcolumn(grp,col,find)
type(group_data), intent(in) :: grp
character(len=*), intent(in) :: col
type(column_data), pointer, intent(inout) :: find

find=>grp%colhead
do while (associated(find))
  if (find%colname==col) exit
  find=>find%next
enddo

end subroutine findcolumn

! ------------------------------------------------------------------------------
subroutine ufo_obsdb_setup_c(c_key_self, c_conf) bind(c,name='ufo_obsdb_setup_f90')
implicit none
integer(c_int), intent(inout) :: c_key_self
type(c_ptr), intent(in)    :: c_conf !< configuration

type(obs_data), pointer :: self
character(len=max_string) :: fin, fout
character(len=max_string+30) :: record

if (config_element_exists(c_conf,"ObsData.ObsDataIn")) then
  fin  = config_get_string(c_conf,max_string,"ObsData.ObsDataIn.obsfile")
else
  fin  = ""
endif
write(record,*)'ufo_obsdb_setup_c: file in =',trim(fin)
call fckit_log%info(record)

fout = config_get_string(c_conf,max_string,"ObsData.ObsDataOut.obsfile")
write(record,*)'ufo_obsdb_setup_c: file out=',trim(fout)
call fckit_log%info(record)

call obs_data_registry%init()
call obs_data_registry%add(c_key_self)
call obs_data_registry%get(c_key_self, self)
call obs_setup(trim(fin), trim(fout), self)

end subroutine ufo_obsdb_setup_c

! ------------------------------------------------------------------------------

subroutine ufo_obsdb_delete_c(c_key_self) bind(c,name='ufo_obsdb_delete_f90')
implicit none
integer(c_int), intent(inout) :: c_key_self
type(obs_data), pointer :: self

call obs_data_registry%get(c_key_self, self)
call obs_delete(self)
call obs_data_registry%remove(c_key_self)

end subroutine ufo_obsdb_delete_c

! ------------------------------------------------------------------------------

subroutine ufo_obsdb_get_c(c_key_self, lreq, c_req, lcol, c_col, c_key_ovec) bind(c,name='ufo_obsdb_get_f90')
implicit none
integer(c_int), intent(in) :: c_key_self
integer(c_int), intent(in) :: lreq, lcol
character(kind=c_char,len=1), intent(in) :: c_req(lreq+1), c_col(lcol+1)
integer(c_int), intent(in) :: c_key_ovec

type(obs_data), pointer :: self
type(obs_vect), pointer :: ovec
character(len=lreq) :: req
character(len=lcol) :: col

call obs_data_registry%get(c_key_self, self)
call ufo_obs_vect_registry%get(c_key_ovec,ovec)
call c_f_string(c_req, req)
call c_f_string(c_col, col)

call obs_get(self, trim(req), trim(col), ovec)

end subroutine ufo_obsdb_get_c

! ------------------------------------------------------------------------------

subroutine ufo_obsdb_put_c(c_key_self, lreq, c_req, lcol, c_col, c_key_ovec) bind(c,name='ufo_obsdb_put_f90')
implicit none
integer(c_int), intent(in) :: c_key_self
integer(c_int), intent(in) :: lreq, lcol
character(kind=c_char,len=1), intent(in) :: c_req(lreq+1), c_col(lcol+1)
integer(c_int), intent(in) :: c_key_ovec

type(obs_data), pointer :: self
type(obs_vect), pointer :: ovec
character(len=lreq) :: req
character(len=lcol) :: col

call obs_data_registry%get(c_key_self, self)
call ufo_obs_vect_registry%get(c_key_ovec,ovec)
call c_f_string(c_req, req)
call c_f_string(c_col, col)

call obs_put(self, trim(req), trim(col), ovec)

end subroutine ufo_obsdb_put_c

! ------------------------------------------------------------------------------

end module ufo_obs_data
