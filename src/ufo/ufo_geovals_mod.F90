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
public :: ufo_geovals, geovals_setup
public :: ufo_geovals_registry

! ------------------------------------------------------------------------------

!> type to hold interpolated field for one variable, one observation
type :: ufo_geoval
  real, allocatable :: vals(:)   !< values (vertical profile or single value for now)
  integer :: nvals               !< number of values in vals array
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

subroutine geovals_setup(self, vars, kobs)
implicit none
type(ufo_geovals), intent(inout) :: self
type(ufo_vars), intent(in) :: vars
integer, intent(in) :: kobs(:)

self%nobs = size(kobs)
self%nvar = vars%nv
self%variables = vars

allocate(self%geovals(self%nvar,self%nobs))
self%lalloc = .true.

end subroutine geovals_setup

! ------------------------------------------------------------------------------

subroutine ufo_geovals_delete_c(c_key_self) bind(c,name='ufo_geovals_delete_f90')

implicit none
integer(c_int), intent(inout) :: c_key_self

type(ufo_geovals), pointer :: self
integer :: i, j


call ufo_geovals_registry%get(c_key_self, self)
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
endif
call ufo_geovals_registry%remove(c_key_self)

end subroutine ufo_geovals_delete_c

! ------------------------------------------------------------------------------

subroutine ufo_geovals_get_var(self, iobs, varname, geoval)
implicit none
type(ufo_geovals), intent(in)    :: self
integer, intent(in)              :: iobs
character(MAXVARLEN), intent(in) :: varname
type(ufo_geoval), intent(out)    :: geoval

integer :: ivar

ivar = ufo_vars_getindex(self%variables, varname)
if (ivar < 0) then
  ! abort!, no such variable
endif
geoval = self%geovals(ivar, iobs)


end subroutine ufo_geovals_get_var

! ------------------------------------------------------------------------------

subroutine ufo_geovals_zero_c(c_key_self) bind(c,name='ufo_geovals_zero_f90')
implicit none
integer(c_int), intent(in) :: c_key_self
type(ufo_geovals), pointer :: self
integer :: i, j

call ufo_geovals_registry%get(c_key_self, self)

if (.not. self%linit) then
  ! TODO: abort! for now just allocating 1
  do i = 1, self%nvar
    do j = 1, self%nobs
      self%geovals(i,j)%nvals = 1
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

!call random_vector(self%values(:,:))

if (.not. self%linit) then
  ! TODO: abort! for now just allocating 1
  do i = 1, self%nvar
    do j = 1, self%nobs
      self%geovals(i,j)%nvals = 1
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
implicit none
integer(c_int), intent(in) :: c_key_self
type(c_ptr), intent(in)    :: c_conf

type(ufo_geovals), pointer :: self
character(128) :: filename
integer :: i, j
integer, allocatable :: nlen(:)
type(ufo_geoval) :: geoval
character(MAXVARLEN) :: varname

call ufo_geovals_registry%get(c_key_self, self)

! read filename for config
filename = config_get_string(c_conf,len(filename),"filename")

call ufo_vars_readconfig(self%variables, c_conf)

self%nobs=1
self%nvar=self%variables%nv

allocate(self%geovals(self%nvar,self%nobs))

self%lalloc = .true.

allocate(nlen(self%nvar))
nlen(:) = 1
do i = 1, self%nvar
  do j = 1, self%nobs
    self%geovals(i,j)%nvals = nlen(i)
    allocate(self%geovals(i,j)%vals(nlen(i)))
    self%geovals(i,j)%vals(:) = i
  enddo
enddo
deallocate(nlen)

self%linit = .true.

!varname = 'u'
!call ufo_geovals_get_var(self, 1, varname, geoval)
!print *, 'geoval test: ', geoval%nvals, geoval%vals
!varname = 'prse'
!call ufo_geovals_get_var(self, 1, varname, geoval)
!print *, 'geoval test: ', geoval%nvals, geoval%vals
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
