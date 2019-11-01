!
! (C) Copyright 2017-2018 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
!
module ufo_geovals_mod_c

use fckit_configuration_module, only: fckit_configuration
use fckit_mpi_module, only: fckit_mpi_comm
use iso_c_binding
use ufo_geovals_mod
use ufo_locs_mod
use ufo_locs_mod_c, only : ufo_locs_registry
use kinds

implicit none

public :: ufo_geovals_registry

private
integer, parameter :: max_string=800

#define LISTED_TYPE ufo_geovals

!> Linked list interface - defines registry_t type
#include "oops/util/linkedList_i.f"

!> Global registry
type(registry_t) :: ufo_geovals_registry

! ------------------------------------------------------------------------------
contains
! ------------------------------------------------------------------------------
!> Linked list implementation
#include "oops/util/linkedList_c.f"
! ------------------------------------------------------------------------------
!> Setup GeoVaLs (store nlocs, variables; don't do allocation yet)
subroutine ufo_geovals_setup_c(c_key_self, c_nlocs, c_vars) bind(c,name='ufo_geovals_setup_f90')
use oops_variables_mod
implicit none
integer(c_int), intent(inout)  :: c_key_self
integer(c_int), intent(in)     :: c_nlocs
type(c_ptr), value, intent(in) :: c_vars

type(ufo_geovals), pointer :: self
type(oops_variables) :: vars

call ufo_geovals_registry%init()
call ufo_geovals_registry%add(c_key_self)
call ufo_geovals_registry%get(c_key_self, self)

vars = oops_variables(c_vars)
call ufo_geovals_setup(self, vars, c_nlocs)

end subroutine ufo_geovals_setup_c

! ------------------------------------------------------------------------------
!> Copy one GeoVaLs object into another

subroutine ufo_geovals_copy_c(c_key_self, c_key_other) bind(c,name='ufo_geovals_copy_f90')
implicit none
integer(c_int), intent(in) :: c_key_self
integer(c_int), intent(in) :: c_key_other
type(ufo_geovals), pointer :: self
type(ufo_geovals), pointer :: other

call ufo_geovals_registry%get(c_key_self, self)
call ufo_geovals_registry%get(c_key_other, other)

call ufo_geovals_copy(self, other)

end subroutine ufo_geovals_copy_c

! ------------------------------------------------------------------------------
!> Analytic init

subroutine ufo_geovals_analytic_init_c(c_key_self, c_key_locs, c_conf) bind(c,name='ufo_geovals_analytic_init_f90')

implicit none
integer(c_int), intent(in) :: c_key_self
integer(c_int), intent(in) :: c_key_locs
type(c_ptr), intent(in)    :: c_conf

type(ufo_geovals), pointer :: self
type(ufo_locs), pointer :: locs
character(len=30) :: ic
character(len=:), allocatable :: str
type(fckit_configuration) :: f_conf

call ufo_geovals_registry%get(c_key_self, self)
call ufo_locs_registry%get(c_key_locs,locs)

f_conf = fckit_configuration(c_conf)
call f_conf%get_or_die("analytic_init",str)
ic = str
call ufo_geovals_analytic_init(self,locs,ic)

end subroutine ufo_geovals_analytic_init_c

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

subroutine ufo_geovals_abs_c(c_key_self) bind(c,name='ufo_geovals_abs_f90')
implicit none
integer(c_int), intent(in) :: c_key_self
type(ufo_geovals), pointer :: self

call ufo_geovals_registry%get(c_key_self, self)

call ufo_geovals_abs(self)

end subroutine ufo_geovals_abs_c

! ------------------------------------------------------------------------------

subroutine ufo_geovals_rms_c(c_key_self,vrms) bind(c,name='ufo_geovals_rms_f90')
implicit none
integer(c_int), intent(in) :: c_key_self
real(c_double), intent(inout) :: vrms
type(ufo_geovals), pointer :: self

call ufo_geovals_registry%get(c_key_self, self)

call ufo_geovals_rms(self,vrms)

end subroutine ufo_geovals_rms_c

! ------------------------------------------------------------------------------

subroutine ufo_geovals_random_c(c_key_self) bind(c,name='ufo_geovals_random_f90')
implicit none
integer(c_int), intent(in) :: c_key_self
type(ufo_geovals), pointer :: self

call ufo_geovals_registry%get(c_key_self, self)

call ufo_geovals_random(self)

end subroutine ufo_geovals_random_c

! ------------------------------------------------------------------------------

subroutine ufo_geovals_scalmult_c(c_key_self, zz) bind(c,name='ufo_geovals_scalmult_f90')
implicit none
integer(c_int), intent(in) :: c_key_self
real(c_double), intent(in) :: zz
type(ufo_geovals), pointer :: self

call ufo_geovals_registry%get(c_key_self, self)

call ufo_geovals_scalmult(self, zz)

end subroutine ufo_geovals_scalmult_c

! ------------------------------------------------------------------------------

subroutine ufo_geovals_assign_c(c_key_self, c_key_rhs) bind(c,name='ufo_geovals_assign_f90')
implicit none
integer(c_int), intent(in) :: c_key_self
integer(c_int), intent(in) :: c_key_rhs
type(ufo_geovals), pointer :: self
type(ufo_geovals), pointer :: rhs

call ufo_geovals_registry%get(c_key_self, self)
call ufo_geovals_registry%get(c_key_rhs, rhs)

call ufo_geovals_assign(self, rhs)

end subroutine ufo_geovals_assign_c

! ------------------------------------------------------------------------------

subroutine ufo_geovals_add_c(c_key_self, c_key_other) bind(c,name='ufo_geovals_add_f90')
implicit none
integer(c_int), intent(in) :: c_key_self
integer(c_int), intent(in) :: c_key_other
type(ufo_geovals), pointer :: self
type(ufo_geovals), pointer :: other

call ufo_geovals_registry%get(c_key_self, self)
call ufo_geovals_registry%get(c_key_other, other)

call ufo_geovals_add(self, other)

end subroutine ufo_geovals_add_c

! ------------------------------------------------------------------------------

subroutine ufo_geovals_diff_c(c_key_self, c_key_other) bind(c,name='ufo_geovals_diff_f90')
implicit none
integer(c_int), intent(in) :: c_key_self
integer(c_int), intent(in) :: c_key_other
type(ufo_geovals), pointer :: self
type(ufo_geovals), pointer :: other

call ufo_geovals_registry%get(c_key_self, self)
call ufo_geovals_registry%get(c_key_other, other)

call ufo_geovals_diff(self, other)

end subroutine ufo_geovals_diff_c

! ------------------------------------------------------------------------------

subroutine ufo_geovals_schurmult_c(c_key_self, c_key_other) bind(c,name='ufo_geovals_schurmult_f90')
implicit none
integer(c_int), intent(in) :: c_key_self
integer(c_int), intent(in) :: c_key_other
type(ufo_geovals), pointer :: self
type(ufo_geovals), pointer :: other

call ufo_geovals_registry%get(c_key_self, self)
call ufo_geovals_registry%get(c_key_other, other)

call ufo_geovals_schurmult(self, other)

end subroutine ufo_geovals_schurmult_c

! ------------------------------------------------------------------------------

subroutine ufo_geovals_normalize_c(c_key_self, c_key_other) bind(c,name='ufo_geovals_normalize_f90')
implicit none
integer(c_int), intent(in) :: c_key_self
integer(c_int), intent(in) :: c_key_other
type(ufo_geovals), pointer :: self
type(ufo_geovals), pointer :: other

call ufo_geovals_registry%get(c_key_self, self)
call ufo_geovals_registry%get(c_key_other, other)

call ufo_geovals_normalize(self, other)

end subroutine ufo_geovals_normalize_c

! ------------------------------------------------------------------------------

subroutine ufo_geovals_dotprod_c(c_key_self, c_key_other, prod, lcname, cname) bind(c,name='ufo_geovals_dotprod_f90')
use string_f_c_mod
implicit none
integer(c_int), intent(in) :: c_key_self, c_key_other
real(c_double), intent(inout) :: prod
integer(c_int),intent(in) :: lcname                        !< Communicator name length
character(kind=c_char,len=1),intent(in) :: cname(lcname+1) !< Communicator name

type(ufo_geovals), pointer :: self, other
type(fckit_mpi_comm) :: f_comm
character(len=lcname) :: name

call ufo_geovals_registry%get(c_key_self, self)
call ufo_geovals_registry%get(c_key_other, other)

call c_f_string(cname, name)
f_comm = fckit_mpi_comm(name)


call ufo_geovals_dotprod(self, other, prod, f_comm)

end subroutine ufo_geovals_dotprod_c

! ------------------------------------------------------------------------------

subroutine ufo_geovals_minmaxavg_c(c_key_self, kobs, kvar, pmin, pmax, prms) bind(c,name='ufo_geovals_minmaxavg_f90')
implicit none
integer(c_int), intent(in) :: c_key_self
integer(c_int), intent(inout) :: kobs
integer(c_int), intent(in) :: kvar
real(c_double), intent(inout) :: pmin, pmax, prms
type(ufo_geovals), pointer :: self

call ufo_geovals_registry%get(c_key_self, self)

call ufo_geovals_minmaxavg(self, kobs, kvar, pmin, pmax, prms)

end subroutine ufo_geovals_minmaxavg_c

! ------------------------------------------------------------------------------

subroutine ufo_geovals_nlocs_c(c_key_self, kobs) bind(c, name='ufo_geovals_nlocs_f90')
implicit none
integer(c_int),    intent(in) :: c_key_self
integer(c_size_t), intent(inout) :: kobs
type(ufo_geovals), pointer :: self

call ufo_geovals_registry%get(c_key_self, self)
kobs = self%nlocs

end subroutine ufo_geovals_nlocs_c

! ------------------------------------------------------------------------------

subroutine ufo_geovals_nlevs_c(c_key_self, lvar, c_var, nlevs) bind(c, name='ufo_geovals_nlevs_f90')
use ufo_vars_mod, only: MAXVARLEN
use string_f_c_mod
implicit none
integer(c_int), intent(in) :: c_key_self
integer(c_int), intent(in) :: lvar
character(kind=c_char, len=1), intent(in) :: c_var(lvar+1)
integer(c_int), intent(out) :: nlevs

type(ufo_geoval), pointer :: geoval
character(len=MAXVARLEN) :: varname
type(ufo_geovals), pointer :: self

call c_f_string(c_var, varname)
call ufo_geovals_registry%get(c_key_self, self)

call ufo_geovals_get_var(self, varname, geoval)

nlevs = size(geoval%vals,1)

end subroutine ufo_geovals_nlevs_c

! ------------------------------------------------------------------------------

subroutine ufo_geovals_get2d_c(c_key_self, lvar, c_var, nlocs, values) bind(c, name='ufo_geovals_get2d_f90')
use ufo_vars_mod, only: MAXVARLEN
use string_f_c_mod
implicit none
integer(c_int), intent(in) :: c_key_self
integer(c_int), intent(in) :: lvar
character(kind=c_char, len=1), intent(in) :: c_var(lvar+1)
integer(c_int), intent(in) :: nlocs
real(c_float), intent(inout) :: values(nlocs)

character(max_string) :: err_msg
type(ufo_geoval), pointer :: geoval
character(len=MAXVARLEN) :: varname
type(ufo_geovals), pointer :: self

call c_f_string(c_var, varname)
call ufo_geovals_registry%get(c_key_self, self)

call ufo_geovals_get_var(self, varname, geoval)

if (size(geoval%vals,1) /= 1) then
  write(err_msg,*)'ufo_geovals_get2d_f90',trim(varname),'is not a 2D var:',size(geoval%vals,1), ' levels'
  call abor1_ftn(err_msg)
endif
if (nlocs /= size(geoval%vals,2)) then
  write(err_msg,*)'ufo_geovals_get2d_f90',trim(varname),'error locs number:',nlocs,size(geoval%vals,2)
  call abor1_ftn(err_msg)
endif

values(:) = geoval%vals(1,:)

end subroutine ufo_geovals_get2d_c

! ------------------------------------------------------------------------------

subroutine ufo_geovals_get_c(c_key_self, lvar, c_var, lev, nlocs, values) bind(c, name='ufo_geovals_get_f90')
use ufo_vars_mod, only: MAXVARLEN
use string_f_c_mod
implicit none
integer(c_int), intent(in) :: c_key_self
integer(c_int), intent(in) :: lvar
character(kind=c_char, len=1), intent(in) :: c_var(lvar+1)
integer(c_int), intent(in) :: lev
integer(c_int), intent(in) :: nlocs
real(c_float), intent(inout) :: values(nlocs)

character(max_string) :: err_msg
type(ufo_geoval), pointer :: geoval
character(len=MAXVARLEN) :: varname
type(ufo_geovals), pointer :: self

call c_f_string(c_var, varname)
call ufo_geovals_registry%get(c_key_self, self)

call ufo_geovals_get_var(self, varname, geoval)

if (lev<1 .or. lev>size(geoval%vals,1)) then
  write(err_msg,*)'ufo_geovals_get_f90',trim(varname),'level out of range:',lev,size(geoval%vals,1)
  call abor1_ftn(err_msg)
endif
if (nlocs /= size(geoval%vals,2)) then
  write(err_msg,*)'ufo_geovals_get_f90',trim(varname),'error locs number:',nlocs,size(geoval%vals,2)
  call abor1_ftn(err_msg)
endif

values(:) = geoval%vals(lev,:)

end subroutine ufo_geovals_get_c

! ------------------------------------------------------------------------------

subroutine ufo_geovals_maxloc_c(c_key_self, mxval, iloc, ivar) bind(c,name='ufo_geovals_maxloc_f90')
implicit none
integer(c_int), intent(in) :: c_key_self
real(c_double), intent(inout) :: mxval
integer(c_int), intent(inout) :: iloc, ivar
type(ufo_geovals), pointer :: self

call ufo_geovals_registry%get(c_key_self, self)

call ufo_geovals_maxloc(self, mxval, iloc, ivar)

end subroutine ufo_geovals_maxloc_c

! ------------------------------------------------------------------------------

subroutine ufo_geovals_read_file_c(c_key_self, c_conf, c_obspace, c_vars) bind(c,name='ufo_geovals_read_file_f90')
use oops_variables_mod
use datetime_mod

implicit none
integer(c_int), intent(inout)  :: c_key_self
type(c_ptr), intent(in)        :: c_conf
type(c_ptr), value, intent(in) :: c_obspace
type(c_ptr), value, intent(in) :: c_vars

type(ufo_geovals), pointer :: self
character(max_string)      :: filename
integer :: loc_multiplier
character(len=:), allocatable :: str
type(fckit_configuration) :: f_conf
type(oops_variables)      :: vars

call ufo_geovals_registry%init()
call ufo_geovals_registry%add(c_key_self)
call ufo_geovals_registry%get(c_key_self, self)

! read filename for config
f_conf = fckit_configuration(c_conf)
call f_conf%get_or_die("filename",str)
filename = str

if (f_conf%has("loc_multiplier")) then
  call f_conf%get_or_die("loc_multiplier", loc_multiplier)
else
  loc_multiplier = 1
endif

vars = oops_variables(c_vars)
! read geovals
call ufo_geovals_read_netcdf(self, filename, loc_multiplier, c_obspace, vars)

end subroutine ufo_geovals_read_file_c

! ------------------------------------------------------------------------------

subroutine ufo_geovals_write_file_c(c_key_self, c_conf, lcname, cname) bind(c,name='ufo_geovals_write_file_f90')
use string_f_c_mod
implicit none
integer(c_int), intent(in) :: c_key_self
type(c_ptr), intent(in) :: c_conf
integer(c_int),intent(in) :: lcname                        !< Communicator name length
character(kind=c_char,len=1),intent(in) :: cname(lcname+1) !< Communicator name

type(ufo_geovals), pointer :: self
character(max_string) :: fout, filename

type(fckit_mpi_comm)      :: comm
character(len=10)         :: cproc
integer                   :: ppos
character(len=:), allocatable :: str
type(fckit_configuration) :: f_conf
character(len=lcname) :: name

! read filename for config
f_conf = fckit_configuration(c_conf)
call f_conf%get_or_die("filename",str)
filename = str

! get the process rank number
call c_f_string(cname, name)
comm= fckit_mpi_comm(name)

write(cproc,fmt='(i4.4)') comm%rank()

! Find the left-most dot in the file name, and use that to pick off the file name
! and file extension.
ppos = scan(trim(filename), '.', BACK=.true.)
if (ppos > 0) then
 ! found a file extension
 fout = filename(1:ppos-1) // '_' // trim(adjustl(cproc)) // trim(filename(ppos:))
else
 ! no file extension
 fout = trim(filename) // '_' // trim(adjustl(cproc))
endif

call ufo_geovals_registry%get(c_key_self, self)
call ufo_geovals_write_netcdf(self, fout)

end subroutine ufo_geovals_write_file_c

! ------------------------------------------------------------------------------

end module ufo_geovals_mod_c
