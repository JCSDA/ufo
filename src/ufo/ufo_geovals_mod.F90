!
! (C) Copyright 2017 UCAR
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
!
module ufo_geovals_mod

use iso_c_binding
use ufo_vars_mod
use kinds

implicit none
private
public :: ufo_geovals, ufo_geovals_setup
public :: ufo_geovals_registry

! ------------------------------------------------------------------------------

!> type to hold interpolated field for one variable, one observation
type :: ufo_geoval
  real, allocatable :: vals(:)   !< values (vertical profile or single value for now)
  integer :: nval                !< number of values in vals array
end type ufo_geoval

!> type to hold interpolated fields required by the obs operators
type :: ufo_geovals
  integer :: nobs                !< number of observations
  integer :: nvar                !< number of variables (supposed to be the
                                 !  same for same obs operator

  type(ufo_geoval), allocatable :: geovals(:,:)  !< array of interpolated
                                                 !  vertical profiles (nvar, nobs)

  type(ufo_vars) :: variables    !< variables list

  logical :: lalloc              !< .true. if type was initialized and allocated
                                 !  (only geovals are allocated, not the arrays
                                 !   inside of the ufo_geoval type)
  logical :: linit               !< .true. if all the ufo_geoval arrays inside geovals
                                 !  were allocated and have data
end type ufo_geovals

#define LISTED_TYPE ufo_geovals

!> Linked list interface - defines registry_t type
#include "linkedList_i.f"

!> Global registry
type(registry_t) :: ufo_geovals_registry

! ------------------------------------------------------------------------------
contains
! ------------------------------------------------------------------------------
!> Linked list implementation
#include "linkedList_c.f"

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
subroutine ufo_geovals_setup(self, vars, nobs)
implicit none
type(ufo_geovals), intent(inout) :: self
type(ufo_vars), intent(in) :: vars
integer, intent(in) :: nobs

call ufo_geovals_delete(self)

self%nobs = nobs
self%nvar = vars%nv
self%variables = vars

allocate(self%geovals(self%nvar,self%nobs))
self%lalloc = .true.

end subroutine ufo_geovals_setup

! ------------------------------------------------------------------------------

subroutine ufo_geovals_delete(self)
implicit none
type(ufo_geovals), intent(inout) :: self

integer :: i, j

if (self%linit) then
  do i = 1, self%nvar
    do j = 1, self%nobs
      deallocate(self%geovals(i,j)%vals)
    enddo
  enddo
  self%linit = .false.
endif
if (self%lalloc) then
  deallocate(self%geovals)
  self%lalloc = .false.
  self%nvar = 0
  self%nobs = 0
endif

end subroutine ufo_geovals_delete

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

logical function ufo_geovals_get_var(self, iobs, varname, geoval)
implicit none
type(ufo_geovals), intent(in)    :: self
integer, intent(in)              :: iobs
character(MAXVARLEN), intent(in) :: varname
type(ufo_geoval), intent(out)    :: geoval

integer :: ivar

if (.not. self%lalloc .or. .not. self%linit) then
  call abor1_ftn("ufo_geovals_get_var: geovals not allocated")
endif

ivar = ufo_vars_getindex(self%variables, varname)
if (ivar < 0) then
  ufo_geovals_get_var = .false.
else
  ufo_geovals_get_var = .true.
  geoval = self%geovals(ivar, iobs)
endif

end function ufo_geovals_get_var

! ------------------------------------------------------------------------------

subroutine ufo_geovals_zero_c(c_key_self) bind(c,name='ufo_geovals_zero_f90')
implicit none
integer(c_int), intent(in) :: c_key_self
type(ufo_geovals), pointer :: self
integer :: i, j

call ufo_geovals_registry%get(c_key_self, self)

if (.not. self%lalloc) then
  call abor1_ftn("ufo_geovals_zero: geovals not allocated")
endif
if (.not. self%linit) then
  ! TODO: abort! for now just allocating 1
  do i = 1, self%nvar
    do j = 1, self%nobs
      self%geovals(i,j)%nval = 1
      allocate(self%geovals(i,j)%vals(1))
    enddo
  enddo
  self%linit = .true.
endif
do i = 1, self%nvar
  do j = 1, self%nobs
    self%geovals(i,j)%vals = 0.0_kind_real
  enddo
enddo
end subroutine ufo_geovals_zero_c

! ------------------------------------------------------------------------------

subroutine ufo_geovals_random_c(c_key_self) bind(c,name='ufo_geovals_random_f90')
use random_vectors_mod
implicit none
integer(c_int), intent(in) :: c_key_self
type(ufo_geovals), pointer :: self
integer :: i, j

call ufo_geovals_registry%get(c_key_self, self)

if (.not. self%lalloc) then
  call abor1_ftn("ufo_geovals_random: geovals not allocated")
endif
if (.not. self%linit) then
  ! TODO: abort! for now just allocating 1
  do i = 1, self%nvar
    do j = 1, self%nobs
      self%geovals(i,j)%nval = 1
      allocate(self%geovals(i,j)%vals(1))
    enddo
  enddo
  self%linit = .true.
endif
do i = 1, self%nvar
  do j = 1, self%nobs
    self%geovals(i,j)%vals = 1.0_kind_real
  enddo
enddo

end subroutine ufo_geovals_random_c

! ------------------------------------------------------------------------------

subroutine ufo_geovals_dotprod_c(c_key_self, c_key_other, prod) bind(c,name='ufo_geovals_dotprod_f90')
implicit none
integer(c_int), intent(in) :: c_key_self, c_key_other
real(c_double), intent(inout) :: prod
type(ufo_geovals), pointer :: self, other
integer :: jo

call ufo_geovals_registry%get(c_key_self, self)
call ufo_geovals_registry%get(c_key_other, other)

if (.not. self%lalloc .or. .not. self%linit) then
  call abor1_ftn("ufo_geovals_dotprod: geovals not allocated")
endif

prod=0.0_kind_real
do jo=1,self%nobs
  prod=prod+self%geovals(1,jo)%vals(1)*other%geovals(1,jo)%vals(1)
enddo

end subroutine ufo_geovals_dotprod_c

! ------------------------------------------------------------------------------

subroutine ufo_geovals_minmaxavg_c(c_key_self, kobs, pmin, pmax, prms) bind(c,name='ufo_geovals_minmaxavg_f90')
implicit none
integer(c_int), intent(in) :: c_key_self
integer(c_int), intent(inout) :: kobs
real(c_double), intent(inout) :: pmin, pmax, prms
type(ufo_geovals), pointer :: self

call ufo_geovals_registry%get(c_key_self, self)

kobs = self%nobs
pmin=0. !minval(self%values(:,:))
pmax=0. !maxval(self%values(:,:))
prms=0. !sqrt(sum(self%values(:,:)**2)/real(self%nobs*self%nvar,kind_real))

end subroutine ufo_geovals_minmaxavg_c

! ------------------------------------------------------------------------------

subroutine ufo_geovals_read_file_c(c_key_self, c_conf) bind(c,name='ufo_geovals_read_file_f90')
use config_mod
use fckit_log_module, only : fckit_log
use nc_diag_read_mod, only: nc_diag_read_get_var,  nc_diag_read_get_global_attr
use nc_diag_read_mod, only: nc_diag_read_get_dim
use nc_diag_read_mod, only: nc_diag_read_init, nc_diag_read_close

implicit none
integer(c_int), intent(in) :: c_key_self
type(c_ptr), intent(in)    :: c_conf

type(ufo_geovals), pointer :: self
type(ufo_vars) :: vars

integer :: nvar_prof, nvar_surf_real, nvar_surf_int
integer :: it, iwv, ipr, iprl, ioz

character(128) :: filename

integer :: iunit
integer :: nobs, nsig, nsig_plus_one

real(8), allocatable :: field(:,:)
real(8), allocatable :: field1d(:)
integer, allocatable :: field1di(:)

integer :: iobs, ivar, nval
type(ufo_geoval) :: geoval
character(MAXVARLEN) :: varname
logical :: lfound

call ufo_geovals_registry%get(c_key_self, self)

! variables hardcoded for the CRTM
nvar_prof = 5
it = 1; iwv = 2; ipr = 3; iprl = 4; ioz = 5 ! indices of vars
nvar_surf_real = 9;
nvar_surf_int = 1;

! allocate and fill in variables
self%nvar = nvar_prof + nvar_surf_real + nvar_surf_int

vars%nv = nvar_prof + nvar_surf_real + nvar_surf_int
allocate(vars%fldnames(vars%nv))
vars%fldnames(it)   = 'Temperature'
vars%fldnames(iwv)  = 'Water vapor'
vars%fldnames(ipr)  = 'Pressure'
vars%fldnames(iprl) = 'Level pressure'
vars%fldnames(ioz)  = 'Ozone'

vars%fldnames(nvar_prof+1) = 'Water_Fraction'
vars%fldnames(nvar_prof+2) = 'Land_Fraction'
vars%fldnames(nvar_prof+3) = 'Ice_Fraction'
vars%fldnames(nvar_prof+4) = 'Snow_Fraction'
vars%fldnames(nvar_prof+5) = 'Water_Temperature'
vars%fldnames(nvar_prof+6) = 'Land_Temperature'
vars%fldnames(nvar_prof+7) = 'Ice_Temperature'
vars%fldnames(nvar_prof+8) = 'Snow_Temperature'
vars%fldnames(nvar_prof+9) = 'Vegetation_Fraction'
vars%fldnames(nvar_prof+nvar_surf_real+1) = 'Land_Type_Index' ! int!!!

! read filename for config
filename = config_get_string(c_conf,len(filename),"filename")

call nc_diag_read_init(filename, iunit)
nobs = nc_diag_read_get_dim(iunit,'nobs')
nsig = nc_diag_read_get_dim(iunit,'nsig')
nsig_plus_one = nc_diag_read_get_dim(iunit,'nsig_plus_one')

! allocate geovals structure TODO: need to clean it first!
call ufo_geovals_setup(self, vars, nobs)

! read temperature
nval = nsig; ivar = it
allocate(field(nval, nobs))
call nc_diag_read_get_var(iunit, 'tvp', field)
do iobs = 1, nobs
  self%geovals(ivar,iobs)%nval = nval
  allocate(self%geovals(ivar,iobs)%vals(nval))
  self%geovals(ivar,iobs)%vals = field(:,iobs)
enddo
deallocate(field)

! read water vapor (humidity)
nval = nsig; ivar = iwv
allocate(field(nval, nobs))
call nc_diag_read_get_var(iunit, 'qvp', field)
do iobs = 1, nobs
  self%geovals(ivar,iobs)%nval = nval
  allocate(self%geovals(ivar,iobs)%vals(nval))
  self%geovals(ivar,iobs)%vals = 1000.*field(:,iobs)/(1.-field(:,iobs))
enddo
deallocate(field)

! read pressure
nval = nsig; ivar = ipr
allocate(field(nval, nobs))
call nc_diag_read_get_var(iunit, 'prsltmp', field)
do iobs = 1, nobs
  self%geovals(ivar,iobs)%nval = nval
  allocate(self%geovals(ivar,iobs)%vals(nval))
  self%geovals(ivar,iobs)%vals = 10.*field(:,iobs)
enddo
deallocate(field)

! read level pressure
nval = nsig_plus_one; ivar = iprl
allocate(field(nval, nobs))
call nc_diag_read_get_var(iunit, 'prsitmp', field)
do iobs = 1, nobs
  self%geovals(ivar,iobs)%nval = nval
  allocate(self%geovals(ivar,iobs)%vals(nval))
  self%geovals(ivar,iobs)%vals = 10.*field(:,iobs)
enddo
deallocate(field)

! read ozone
nval = nsig; ivar = ioz
allocate(field(nval, nobs))
call nc_diag_read_get_var(iunit, 'poz', field)
do iobs = 1, nobs
  self%geovals(ivar,iobs)%nval = nval
  allocate(self%geovals(ivar,iobs)%vals(nval))
  self%geovals(ivar,iobs)%vals = field(:,iobs)
enddo
deallocate(field)

! read surface stuff
nval = 1
allocate(field1d(nobs))
do ivar = nvar_prof+1, nvar_prof+nvar_surf_real
  call nc_diag_read_get_var(iunit, self%variables%fldnames(ivar), field1d)
  do iobs = 1, nobs
    self%geovals(ivar,iobs)%nval = nval
    allocate(self%geovals(ivar,iobs)%vals(nval))
    self%geovals(ivar,iobs)%vals(1) = field1d(iobs)
  enddo
enddo
deallocate(field1d)

! read surface stuff (integer)
nval = 1
allocate(field1di(nobs))
do ivar = nvar_prof+nvar_surf_real+1, nvar_prof+nvar_surf_real+nvar_surf_int
  call nc_diag_read_get_var(iunit, self%variables%fldnames(ivar), field1di)
  do iobs = 1, nobs
    self%geovals(ivar,iobs)%nval = nval
    allocate(self%geovals(ivar,iobs)%vals(nval))
    self%geovals(ivar,iobs)%vals(1) = field1di(iobs)
  enddo
enddo
deallocate(field1di)


self%linit = .true.

call nc_diag_read_close(filename)

! Example of getting a variable below:
!varname = 'Ozone'
!lfound =  ufo_geovals_get_var(self, 1, varname, geoval)
!if (lfound) then
!  print *, 'geoval test: ', trim(varname), geoval%nval, geoval%vals
!else
!  print *, 'geoval test: ', trim(varname), ' doesnt exist'
!endif

end subroutine ufo_geovals_read_file_c

! ------------------------------------------------------------------------------

subroutine ufo_geovals_write_file_c(c_key_self, c_conf) bind(c,name='ufo_geovals_write_file_f90')
use config_mod
use fckit_log_module, only : fckit_log
implicit none
integer(c_int), intent(in) :: c_key_self
type(c_ptr), intent(in) :: c_conf
type(ufo_geovals), pointer :: self

call ufo_geovals_registry%get(c_key_self, self)

end subroutine ufo_geovals_write_file_c

! ------------------------------------------------------------------------------

end module ufo_geovals_mod
