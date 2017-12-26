!
! (C) Copyright 2017 UCAR
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
!
module ufo_geovals_mod

use ufo_vars_mod
use kinds

implicit none
private
integer, parameter :: max_string=800

public :: ufo_geovals, ufo_geoval, ufo_geovals_get_var
public :: ufo_geovals_init, ufo_geovals_setup, ufo_geovals_delete, ufo_geovals_print
public :: ufo_geovals_zero, ufo_geovals_random, ufo_geovals_dotprod
public :: ufo_geovals_minmaxavg
public :: ufo_geovals_read_netcdf

! ------------------------------------------------------------------------------

!> type to hold interpolated field for one variable, one observation
type :: ufo_geoval
  real(kind_real), allocatable :: vals(:,:) !< values (nval, nobs)
  integer :: nval                !< number of values in profile
  integer :: nobs                !< number of observations
end type ufo_geoval

!> type to hold interpolated fields required by the obs operators
type :: ufo_geovals
  integer :: nobs                !< number of observations
  integer :: nvar                !< number of variables (supposed to be the
                                 !  same for same obs operator

  type(ufo_geoval), allocatable :: geovals(:)  !< array of interpolated
                                               !  vertical profiles for all obs (nvar)

  type(ufo_vars) :: variables    !< variables list

  logical :: lalloc              !< .true. if type was initialized and allocated
                                 !  (only geovals are allocated, not the arrays
                                 !   inside of the ufo_geoval type)
  logical :: linit               !< .true. if all the ufo_geoval arrays inside geovals
                                 !  were allocated and have data
end type ufo_geovals


! ------------------------------------------------------------------------------
contains

! ------------------------------------------------------------------------------

subroutine ufo_geovals_init(self)
implicit none
type(ufo_geovals), intent(inout) :: self

self%lalloc = .false.
self%linit  = .false.
self%nvar = 0
self%nobs = 0

end subroutine ufo_geovals_init

! ------------------------------------------------------------------------------

subroutine ufo_geovals_setup(self, vars, nobs)
implicit none
type(ufo_geovals), intent(inout) :: self
type(ufo_vars), intent(in) :: vars
integer, intent(in) :: nobs

integer :: ivar

call ufo_geovals_delete(self)

self%nobs = nobs
self%nvar = vars%nv
call ufo_vars_clone(vars, self%variables) 

allocate(self%geovals(self%nvar))
do ivar = 1, self%nvar
  self%geovals(ivar)%nobs = nobs
  self%geovals(ivar)%nval = 0
enddo
self%lalloc = .true.

end subroutine ufo_geovals_setup

! ------------------------------------------------------------------------------

subroutine ufo_geovals_delete(self)
implicit none
type(ufo_geovals), intent(inout) :: self

integer :: ivar

if (self%linit) then
  do ivar = 1, self%nvar
    deallocate(self%geovals(ivar)%vals)
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

logical function ufo_geovals_get_var(self, varname, geoval)
implicit none
type(ufo_geovals), intent(in)    :: self
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
  geoval = self%geovals(ivar)
endif

end function ufo_geovals_get_var

! ------------------------------------------------------------------------------

subroutine ufo_geovals_zero(self) 
implicit none
type(ufo_geovals), intent(inout) :: self
integer :: ivar

if (.not. self%lalloc) then
  call abor1_ftn("ufo_geovals_zero: geovals not allocated")
endif
if (.not. self%linit) then
  ! TODO: abort! for now just allocating 1
  do ivar = 1, self%nvar
    self%geovals(ivar)%nval = 1
    allocate(self%geovals(ivar)%vals(1,self%nobs))
  enddo
  self%linit = .true.
endif
do ivar = 1, self%nvar
  self%geovals(ivar)%vals = 0.0
enddo

end subroutine ufo_geovals_zero

! ------------------------------------------------------------------------------

subroutine ufo_geovals_random(self) 
implicit none
type(ufo_geovals), intent(inout) :: self
integer :: ivar

if (.not. self%lalloc) then
  call abor1_ftn("ufo_geovals_random: geovals not allocated")
endif
if (.not. self%linit) then
  ! TODO: abort! for now just allocating 1
  do ivar = 1, self%nvar
    self%geovals(ivar)%nval = 1
    allocate(self%geovals(ivar)%vals(1,self%nobs))
  enddo
  self%linit = .true.
endif
do ivar = 1, self%nvar
  self%geovals(ivar)%vals = 1.0
enddo

end subroutine ufo_geovals_random

! ------------------------------------------------------------------------------

subroutine ufo_geovals_dotprod(self, other, prod) 
implicit none
real(kind_real), intent(inout) :: prod
type(ufo_geovals), intent(in) :: self, other
integer :: jo

if (.not. self%lalloc .or. .not. self%linit) then
  call abor1_ftn("ufo_geovals_dotprod: geovals not allocated")
endif

if (.not. other%lalloc .or. .not. other%linit) then
  call abor1_ftn("ufo_geovals_dotprod: geovals not allocated")
endif

! just something to put in (dot product of the 1st var and 1st element in the profile
prod=0.0
do jo=1,self%nobs
  prod=prod+self%geovals(1)%vals(1,jo)*other%geovals(1)%vals(1,jo)
enddo
prod = prod / real(self%nobs,kind_real)

end subroutine ufo_geovals_dotprod

! ------------------------------------------------------------------------------

subroutine ufo_geovals_minmaxavg(self, kobs, pmin, pmax, prms) 
implicit none
integer, intent(inout) :: kobs
real(kind_real), intent(inout) :: pmin, pmax, prms
type(ufo_geovals), intent(in) :: self

kobs = self%nobs
pmin=0. !minval(self%values(:,:))
pmax=0. !maxval(self%values(:,:))
prms=0. !sqrt(sum(self%values(:,:)**2)/real(self%nobs*self%nvar,kind_real))

end subroutine ufo_geovals_minmaxavg

! ------------------------------------------------------------------------------

subroutine ufo_geovals_read_netcdf(self, filename, vars)
use netcdf, only: NF90_DOUBLE, NF90_INT
use nc_diag_read_mod, only: nc_diag_read_get_var
use nc_diag_read_mod, only: nc_diag_read_get_dim
use nc_diag_read_mod, only: nc_diag_read_get_var_dims, nc_diag_read_check_var
use nc_diag_read_mod, only: nc_diag_read_get_var_type
use nc_diag_read_mod, only: nc_diag_read_init, nc_diag_read_close
implicit none
type(ufo_geovals), intent(inout)  :: self
character(max_string), intent(in) :: filename
type(ufo_vars), intent(in)        :: vars

integer :: iunit, ivar, nobs, nval
integer :: nvardim, vartype
integer, allocatable, dimension(:) :: vardims

real(kind_real), allocatable :: fieldr2d(:,:), fieldr1d(:)
integer, allocatable :: fieldi2d(:,:), fieldi1d(:)

! open netcdf file and read dimensions
call nc_diag_read_init(filename, iunit)
nobs = nc_diag_read_get_dim(iunit,'nobs')

! allocate geovals structure
call ufo_geovals_setup(self, vars, nobs)

do ivar = 1, vars%nv
  if (.not. nc_diag_read_check_var(iunit, vars%fldnames(ivar))) then
     print *, 'ERROR reading geovals: var ', trim(vars%fldnames(ivar)), ' doesnt exist'
     stop 1
  endif
  !> get dimensions of variable
  if (allocated(vardims)) deallocate(vardims)
  call nc_diag_read_get_var_dims(iunit, vars%fldnames(ivar), nvardim, vardims)
  !> get variable type
  vartype = nc_diag_read_get_var_type(iunit, vars%fldnames(ivar))
  !> read 1d vars (only double precision and integer for now)
  if (nvardim == 1) then
    if (vardims(1) /= nobs) then
       print *, 'ERROR reading geovals: var dim /= nobs'
       stop 2
    endif
    nval = 1

    !> allocate geoval for this variable
    self%geovals(ivar)%nval = nval
    allocate(self%geovals(ivar)%vals(nval,nobs))

    if (vartype == NF90_DOUBLE) then
       allocate(fieldr1d(vardims(1)))
       call nc_diag_read_get_var(iunit, vars%fldnames(ivar), fieldr1d)  
       self%geovals(ivar)%vals(1,:) = fieldr1d
       deallocate(fieldr1d)
    elseif (vartype == NF90_INT) then
       allocate(fieldi1d(vardims(1)))
       call nc_diag_read_get_var(iunit, vars%fldnames(ivar), fieldi1d)
       self%geovals(ivar)%vals(1,:) = fieldi1d
       deallocate(fieldi1d)
    else
       print *, 'ERROR reading geovals: only can read double and int'
       stop 3
    endif
  !> read 2d vars (only double precision and integer for now)
  elseif (nvardim == 2) then
    if (vardims(2) /= nobs) then
       print *, 'ERROR reading geovals: var dim /= nobs'
       stop 2
    endif
    nval = vardims(1)

    !> allocate geoval for this variable
    self%geovals(ivar)%nval = nval
    allocate(self%geovals(ivar)%vals(nval,nobs))

    if (vartype == NF90_DOUBLE) then
       allocate(fieldr2d(vardims(1), vardims(2)))
       call nc_diag_read_get_var(iunit, vars%fldnames(ivar), fieldr2d)
       self%geovals(ivar)%vals = fieldr2d
       deallocate(fieldr2d)
    elseif (vartype == NF90_INT) then
       allocate(fieldi2d(vardims(1), vardims(2)))
       call nc_diag_read_get_var(iunit, vars%fldnames(ivar), fieldi2d)
       self%geovals(ivar)%vals = fieldi2d
       deallocate(fieldi2d)
    else
       print *, 'ERROR reading geovals: only can read double and int'
       stop 3
    endif
  !> only 1d & 2d vars
  else
    print *, 'ERROR reading geovals: only can read 1d and 2d fields'
    stop 4
  endif
enddo

self%linit = .true.

call nc_diag_read_close(filename)

end subroutine ufo_geovals_read_netcdf

! ------------------------------------------------------------------------------

subroutine ufo_geovals_print(self, iobs)
implicit none
type(ufo_geovals), intent(in) :: self
integer, intent(in) :: iobs

type(ufo_geoval) :: geoval
character(MAXVARLEN) :: varname
logical :: lfound


integer :: ivar

do ivar = 1, self%nvar
  varname = self%variables%fldnames(ivar)
  lfound =  ufo_geovals_get_var(self, varname, geoval)
  if (lfound) then
    print *, 'geoval test: ', trim(varname), geoval%nval, geoval%vals(:,iobs)
  else
    print *, 'geoval test: ', trim(varname), ' doesnt exist'
  endif
enddo

end subroutine ufo_geovals_print

! ------------------------------------------------------------------------------

end module ufo_geovals_mod
