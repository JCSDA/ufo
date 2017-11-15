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
use read_diag, only: set_radiag,&
                     diag_header_fix_list,&
                     diag_header_chan_list,&
                     diag_data_name_list,&
                     read_radiag_header,&
                     set_netcdf_read, &
                     open_radiag, &
                     close_radiag
use read_diag, only: read_radiag_data,&
                     diag_data_fix_list,&
                     diag_data_extra_list,&
                     diag_data_chan_list
use nc_diag_read_mod, only: nc_diag_read_init, nc_diag_read_close
use kinds

implicit none
private

public :: obs_data, obs_setup, obs_delete, obs_get, obs_put, max_string
public :: obs_data_registry

! ------------------------------------------------------------------------------
integer, parameter :: max_string=800
! ------------------------------------------------------------------------------

!> A type to represent observation data
type obs_data
  integer :: nobs
  character(len=max_string) :: filein, fileout
  type(diag_header_fix_list )              ::  header_fix
  type(diag_header_chan_list),allocatable  ::  header_chan(:)
  type(diag_data_name_list)                ::  header_name

  type(diag_data_fix_list)                 ::  datafix
  type(diag_data_chan_list)  ,allocatable  ::  datachan(:)
  type(diag_data_extra_list) ,allocatable  ::  dataextra(:,:)
end type obs_data

#define LISTED_TYPE obs_data

!> Linked list interface - defines registry_t type
#include "linkedList_i.f"

!> Global registry
type(registry_t) :: obs_data_registry

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

self%filein =fin
self%fileout=fout

! ugly fix for conventional for now
if (self%filein == 'Data/amsua_n19_wprofiles.nc4') then
  call obs_read(self)
  self%nobs = 15
elseif (self%filein == 'Data/diag_t_01_wprofiles.nc4') then
  self%nobs = 915
elseif (self%filein == 'Data/diag_q_01_wprofiles.nc4') then
  self%nobs = 73651
elseif (self%filein == 'Data/diag_uv_01_wprofiles.nc4') then
  self%nobs = 110119
elseif (self%filein == 'Data/diag_ps_01_wprofiles.nc4') then
  self%nobs = 84283
else
  print *, 'Error: dont know how to read ', trim(self%filein)
endif

call fckit_log%debug("TRACE: ufo_obs_data:obs_setup: done")

end subroutine obs_setup

! ------------------------------------------------------------------------------

subroutine obs_delete(self)
implicit none
type(obs_data), intent(inout) :: self

if (self%fileout/="") call obs_write(self)

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

subroutine ufo_obsdb_locations_c(c_key_self, c_t1, c_t2, c_key_locs) bind(c,name='ufo_obsdb_locations_f90')
implicit none
integer(c_int), intent(in) :: c_key_self
type(c_ptr), intent(in) :: c_t1, c_t2
integer(c_int), intent(inout) :: c_key_locs

type(obs_data), pointer :: self
type(datetime) :: t1, t2
type(ufo_locs), pointer :: locs
type(obs_vector) :: ovec
character(len=8) :: col="Location"

call obs_data_registry%get(c_key_self, self)
call c_f_datetime(c_t1, t1)
call c_f_datetime(c_t2, t2)

call obs_time_get(self, col, t1, t2, ovec)

call ufo_locs_registry%init()
call ufo_locs_registry%add(c_key_locs)
call ufo_locs_registry%get(c_key_locs,locs)
     
call ufo_loc_setup(locs, ovec)

!diag_data_fix_list%lon
!diag_data_fix_list%lat
!diag_data_fix_list%obstime

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

call obs_data_registry%get(c_key_self, self)
call ufo_vars_registry%get(c_key_vars, vars)
call c_f_datetime(c_t1, t1)
call c_f_datetime(c_t2, t2)

call ufo_geovals_registry%init()
call ufo_geovals_registry%add(c_key_geovals)
call ufo_geovals_registry%get(c_key_geovals,geovals)

call ufo_geovals_init(geovals)
call ufo_geovals_setup(geovals, vars, self%nobs)

end subroutine ufo_obsdb_getgeovals_c

! ------------------------------------------------------------------------------

subroutine ufo_obsdb_nobs_c(c_key_self, kobs) bind(c,name='ufo_obsdb_nobs_f90')
implicit none
integer(c_int), intent(in) :: c_key_self
integer(c_int), intent(inout) :: kobs
type(obs_data), pointer :: self
call obs_data_registry%get(c_key_self, self)
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

subroutine obs_read(self)
use ncd_kinds, only: i_kind
implicit none
character(len=*),parameter :: myname_ ="ufo_obs_data:obs_read"
type(obs_data), intent(inout) :: self
integer(i_kind) :: ier
integer(i_kind) :: luin=0
integer(i_kind) :: npred = 7   
integer(i_kind) :: iversion=30303
logical :: lverbose  = .true.  ! control verbose
logical :: retrieval = .false. ! true when dealing with SST retrievals
character(len=max_string) :: ncfname


ncfname = self%filein
call set_netcdf_read(.true.)
call open_radiag(ncfname,luin)
call set_radiag("version",iversion,ier)

call read_radiag_header(luin,npred,retrieval,self%header_fix,self%header_chan,self%header_name,ier,lverbose)

print*, myname_, ': Found this many channels: ', self%header_fix%nchan
print*, myname_, ': Observation type in file: ', self%header_fix%obstype
print*, myname_, ': Date of input file:       ', self%header_fix%idate

self%nobs=0
do while (ier .ge. 0)
   call read_radiag_data ( luin, self%header_fix, .false., self%datafix, self%datachan, &
                           self%dataextra, ier )

   if (ier .lt. 0) cycle
   self%nobs = self%nobs + 1
enddo
print *, myname_, ' Total number of observations in file: ', self%nobs
call close_radiag(ncfname, luin)
end subroutine obs_read

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

if (config_element_exists(c_conf,"ObsData.ObsDataIn")) then
  fin  = config_get_string(c_conf,max_string,"ObsData.ObsDataIn.obsfile")
else
  fin  = ""
endif

write(record,*)'ufo_obsdb_setup_c: file in =',trim(fin)
call fckit_log%info(record)

!fout = config_get_string(c_conf,max_string,"ObsData.ObsDataOut.obsfile")
!write(record,*)'ufo_obsdb_setup_c: file out=',trim(fout)
!call fckit_log%info(record)
fout = ""

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

subroutine ufo_obsdb_get_c(c_key_self, lcol, c_col, c_key_ovec) bind(c,name='ufo_obsdb_get_f90')
implicit none
integer(c_int), intent(in) :: c_key_self
integer(c_int), intent(in) :: lcol
character(kind=c_char,len=1), intent(in) :: c_col(lcol+1)
integer(c_int), intent(in) :: c_key_ovec

type(obs_data), pointer :: self
type(obs_vector), pointer :: ovec
character(len=lcol) :: col

call obs_data_registry%get(c_key_self, self)
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

call obs_data_registry%get(c_key_self, self)
call ufo_obs_vect_registry%get(c_key_ovec,ovec)
call c_f_string(c_col, col)

call obs_put(self, trim(col), ovec)

end subroutine ufo_obsdb_put_c

! ------------------------------------------------------------------------------

end module ufo_obs_data
