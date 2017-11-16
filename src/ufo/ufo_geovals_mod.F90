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
public :: ufo_geovals, ufo_geoval, ufo_geovals_get_var
public :: ufo_geovals_registry
public :: ufo_geovals_init, ufo_geovals_setup, ufo_geovals_delete, ufo_geovals_print
public :: ufo_geovals_zero, ufo_geovals_random, ufo_geovals_dotprod
public :: ufo_geovals_minmaxavg
public :: ufo_geovals_read_t_netcdf, ufo_geovals_read_raob_t_netcdf
public :: ufo_geovals_read_q_netcdf
public :: ufo_geovals_read_uv_netcdf, ufo_geovals_read_ps_netcdf
public :: ufo_geovals_read_rad_netcdf

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

call ufo_geovals_delete(self)

self%nobs = nobs
self%nvar = vars%nv
call ufo_vars_clone(vars, self%variables) 

allocate(self%geovals(self%nvar))
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
integer(c_int), intent(inout) :: kobs
real(kind_real), intent(inout) :: pmin, pmax, prms
type(ufo_geovals), intent(in) :: self

kobs = self%nobs
pmin=0. !minval(self%values(:,:))
pmax=0. !maxval(self%values(:,:))
prms=0. !sqrt(sum(self%values(:,:)**2)/real(self%nobs*self%nvar,kind_real))

end subroutine ufo_geovals_minmaxavg

! ------------------------------------------------------------------------------

subroutine ufo_geovals_read_prof_netcdf(self, filename, vars, varsfile)
use nc_diag_read_mod, only: nc_diag_read_get_var
use nc_diag_read_mod, only: nc_diag_read_get_dim
use nc_diag_read_mod, only: nc_diag_read_init, nc_diag_read_close

implicit none
type(ufo_geovals), intent(inout) :: self
character(128), intent(in)       :: filename
type(ufo_vars), intent(in)       :: vars, varsfile

integer :: iunit
integer :: nobs, nsig

real(kind_real), allocatable :: field(:,:)

integer :: ivar, nval

! open netcdf file and read dimensions
call nc_diag_read_init(filename, iunit)
nobs = nc_diag_read_get_dim(iunit,'nobs')
nsig = nc_diag_read_get_dim(iunit,'nsig')

! allocate geovals structure
call ufo_geovals_setup(self, vars, nobs)

nval = nsig
allocate(field(nval, nobs))
do ivar = 1, vars%nv
  call nc_diag_read_get_var(iunit, varsfile%fldnames(ivar), field)
  self%geovals(ivar)%nval = nval
  allocate(self%geovals(ivar)%vals(nval,nobs))
  self%geovals(ivar)%vals = field
enddo
deallocate(field)
self%linit = .true.

call nc_diag_read_close(filename)

end subroutine ufo_geovals_read_prof_netcdf

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

subroutine ufo_geovals_read_t_netcdf(self, filename)
implicit none
type(ufo_geovals), intent(inout) :: self
character(128), intent(in)       :: filename

type(ufo_vars) :: vars, varsfile

integer :: nvar_prof

! variables hardcoded for the temperature
nvar_prof = 6

! allocate and fill in variables
vars%nv = nvar_prof; varsfile%nv = nvar_prof
allocate(vars%fldnames(vars%nv), varsfile%fldnames(varsfile%nv))
vars%fldnames(1) = 'Virtual temperature';  varsfile%fldnames(1) = 'tvtmp'
vars%fldnames(2) = 'Specific humidity';    varsfile%fldnames(2) = 'qtmp'
vars%fldnames(3) = 'U-wind';               varsfile%fldnames(3) = 'utmp'
vars%fldnames(4) = 'V-wind';               varsfile%fldnames(4) = 'vtmp'
vars%fldnames(5) = 'LogPressure';          varsfile%fldnames(5) = 'prsltmp'
vars%fldnames(6) = 'Geopotential height';  varsfile%fldnames(6) = 'hsges'

call ufo_geovals_read_prof_netcdf(self, filename, vars, varsfile)

!call ufo_geovals_print(self, 1)

end subroutine ufo_geovals_read_t_netcdf

! ------------------------------------------------------------------------------
subroutine ufo_geovals_read_raob_t_netcdf(self, filename)
use nc_diag_read_mod, only: nc_diag_read_get_var
use nc_diag_read_mod, only: nc_diag_read_get_dim
use nc_diag_read_mod, only: nc_diag_read_init, nc_diag_read_close

implicit none
type(ufo_geovals), intent(inout) :: self
character(128), intent(in)       :: filename

type(ufo_vars) :: vars, varsfile

integer :: iunit
integer :: nobs, nobs_raob, nsig
integer :: iobs, iobs_raob

real(kind_real), allocatable :: field(:,:)
integer, allocatable         :: obstype(:)
real(kind_real), allocatable :: tvflag(:)

integer :: ivar, nval

integer :: nvar_prof

integer, parameter :: raobtype = 120

! variables hardcoded for the temperature
nvar_prof = 6

! allocate and fill in variables
vars%nv = nvar_prof; varsfile%nv = nvar_prof
allocate(vars%fldnames(vars%nv), varsfile%fldnames(varsfile%nv))
vars%fldnames(1) = 'Virtual temperature';  varsfile%fldnames(1) = 'tvtmp'
vars%fldnames(2) = 'Specific humidity';    varsfile%fldnames(2) = 'qtmp'
vars%fldnames(3) = 'U-wind';               varsfile%fldnames(3) = 'utmp'
vars%fldnames(4) = 'V-wind';               varsfile%fldnames(4) = 'vtmp'
vars%fldnames(5) = 'LogPressure';          varsfile%fldnames(5) = 'prsltmp'
vars%fldnames(6) = 'Geopotential height';  varsfile%fldnames(6) = 'hsges'


! open netcdf file and read dimensions
call nc_diag_read_init(filename, iunit)
nobs = nc_diag_read_get_dim(iunit,'nobs')
nsig = nc_diag_read_get_dim(iunit,'nsig')

allocate(obstype(nobs), tvflag(nobs))
call nc_diag_read_get_var(iunit, "Observation_Type", obstype)
call nc_diag_read_get_var(iunit, "Setup_QC_Mark", tvflag)

nobs_raob = count(obstype == raobtype .and. tvflag == 0)
! allocate geovals structure
call ufo_geovals_setup(self, vars, nobs_raob)

nval = nsig
allocate(field(nval, nobs))
do ivar = 1, vars%nv
  call nc_diag_read_get_var(iunit, varsfile%fldnames(ivar), field)
  self%geovals(ivar)%nval = nval
  allocate(self%geovals(ivar)%vals(nval,nobs_raob))
  iobs_raob = 1
  do iobs = 1, nobs
    if (obstype(iobs) == raobtype .and. tvflag(iobs) == 0) then
      self%geovals(ivar)%vals(:,iobs_raob) = field(:,iobs)
      iobs_raob = iobs_raob + 1
    endif
  enddo
enddo
deallocate(field)
self%linit = .true.

call nc_diag_read_close(filename)

!call ufo_geovals_print(self, 1)

end subroutine ufo_geovals_read_raob_t_netcdf

! ------------------------------------------------------------------------------

subroutine ufo_geovals_read_uv_netcdf(self, filename)
implicit none
type(ufo_geovals), intent(inout) :: self
character(128), intent(in)       :: filename

type(ufo_vars) :: vars, varsfile

integer :: nvar_prof

nvar_prof = 5

! allocate and fill in variables
vars%nv = nvar_prof; varsfile%nv = nvar_prof
allocate(vars%fldnames(vars%nv), varsfile%fldnames(varsfile%nv))
vars%fldnames(1) = 'U-wind';               varsfile%fldnames(1) = 'uges'
vars%fldnames(2) = 'V-wind';               varsfile%fldnames(2) = 'vges'
vars%fldnames(3) = 'LogPressure';          varsfile%fldnames(3) = 'prsltmp'
vars%fldnames(4) = 'Geopotential height';  varsfile%fldnames(4) = 'zges'
vars%fldnames(5) = 'Virtual temperature';  varsfile%fldnames(5) = 'tges'

call ufo_geovals_read_prof_netcdf(self, filename, vars, varsfile)

!call ufo_geovals_print(self, 1)

end subroutine ufo_geovals_read_uv_netcdf

! ------------------------------------------------------------------------------
subroutine ufo_geovals_read_q_netcdf(self, filename)
implicit none
type(ufo_geovals), intent(inout) :: self
character(128), intent(in)       :: filename

type(ufo_vars) :: vars, varsfile

integer :: nvar_prof

! variables hardcoded for the temperature
nvar_prof = 3

! allocate and fill in variables
vars%nv = nvar_prof; varsfile%nv = nvar_prof
allocate(vars%fldnames(vars%nv), varsfile%fldnames(varsfile%nv))
vars%fldnames(1) = 'Specific humidity';    varsfile%fldnames(1) = 'qtmp'
vars%fldnames(2) = 'Saturation';           varsfile%fldnames(2) = 'qgtmp'
vars%fldnames(3) = 'LogPressure';          varsfile%fldnames(3) = 'prsltmp'

call ufo_geovals_read_prof_netcdf(self, filename, vars, varsfile)

end subroutine ufo_geovals_read_q_netcdf

! ------------------------------------------------------------------------------

subroutine ufo_geovals_read_ps_netcdf(self, filename)
use nc_diag_read_mod, only: nc_diag_read_get_var
use nc_diag_read_mod, only: nc_diag_read_get_dim
use nc_diag_read_mod, only: nc_diag_read_init, nc_diag_read_close

implicit none
type(ufo_geovals), intent(inout) :: self
character(128), intent(in)       :: filename

type(ufo_vars) :: vars, varsfile

integer :: nvar_prof, nvar_surf

integer :: iunit
integer :: nobs, nsig

real(kind_real), allocatable :: field(:,:)
real(kind_real), allocatable :: field1d(:)

integer :: ivar, nval

! variables hardcoded for the surface pressure
nvar_prof = 1; nvar_surf = 2

! allocate and fill in variables
vars%nv = nvar_prof+nvar_surf; varsfile%nv = nvar_prof+nvar_surf
allocate(vars%fldnames(vars%nv), varsfile%fldnames(varsfile%nv))
vars%fldnames(1) = 'LogPressure';       varsfile%fldnames(1) = 'prsltmp'

vars%fldnames(2) = 'Surface pressure';  varsfile%fldnames(2) = 'psges'
vars%fldnames(3) = 'Surface height';    varsfile%fldnames(3) = 'zsges'

! open netcdf file and read dimensions
call nc_diag_read_init(filename, iunit)
nobs = nc_diag_read_get_dim(iunit,'nobs')
nsig = nc_diag_read_get_dim(iunit,'nsig')

! allocate geovals structure
call ufo_geovals_setup(self, vars, nobs)

nval = nsig
allocate(field(nval, nobs))
do ivar = 1, nvar_prof
  call nc_diag_read_get_var(iunit, varsfile%fldnames(ivar), field)
  self%geovals(ivar)%nval = nval
  allocate(self%geovals(ivar)%vals(nval,nobs))
  self%geovals(ivar)%vals = field
enddo
deallocate(field)

nval = 1
allocate(field1d(nobs))
do ivar = nvar_prof+1, nvar_prof+nvar_surf
  call nc_diag_read_get_var(iunit, varsfile%fldnames(ivar), field1d)
  self%geovals(ivar)%nval = nval
  allocate(self%geovals(ivar)%vals(nval,nobs))
  self%geovals(ivar)%vals(1,:) = field1d(:)
enddo
deallocate(field1d)

self%linit = .true.

call nc_diag_read_close(filename)

!call ufo_geovals_print(self, 1)

end subroutine ufo_geovals_read_ps_netcdf

! ------------------------------------------------------------------------------

subroutine ufo_geovals_read_rad_netcdf(self, filename)
use nc_diag_read_mod, only: nc_diag_read_get_var
use nc_diag_read_mod, only: nc_diag_read_get_dim
use nc_diag_read_mod, only: nc_diag_read_init, nc_diag_read_close

implicit none
type(ufo_geovals), intent(inout) :: self
character(128), intent(in)       :: filename

type(ufo_vars) :: vars

integer :: nvar_prof, nvar_surf_real, nvar_surf_int
integer :: it, iwv, ipr, iprl, ioz, icl1, icl2

integer :: iunit
integer :: nobs, nchans, nobs_all
integer :: nsig, nsig_plus_one

real(kind_real), allocatable :: field(:,:)
real(kind_real), allocatable :: field1d(:)
integer, allocatable :: field1di(:)

integer :: ivar, iobs, nval

! variables hardcoded for the CRTM  !** Note: we'll need to revisit these in the future -BTJ 11.15.2017
nvar_prof = 7
it = 1; iwv = 2; ipr = 3; iprl = 4; ioz = 5; icl1 = 6; icl2 = 7 ! indices of vars
nvar_surf_real = 14;
nvar_surf_int = 3;

! allocate and fill in variables
self%nvar = nvar_prof + nvar_surf_real + nvar_surf_int

vars%nv = nvar_prof + nvar_surf_real + nvar_surf_int
allocate(vars%fldnames(vars%nv))
vars%fldnames(it)   = 'Temperature'
vars%fldnames(iwv)  = 'Water vapor'
vars%fldnames(ipr)  = 'Pressure'
vars%fldnames(iprl) = 'Level pressure'
vars%fldnames(ioz)  = 'Ozone'
vars%fldnames(icl1) = 'Cloud liquid'
vars%fldnames(icl2) = 'Cloud ice'

vars%fldnames(nvar_prof+1) = 'Water_Fraction'
vars%fldnames(nvar_prof+2) = 'Land_Fraction'
vars%fldnames(nvar_prof+3) = 'Ice_Fraction'
vars%fldnames(nvar_prof+4) = 'Snow_Fraction'
vars%fldnames(nvar_prof+5) = 'Water_Temperature'
vars%fldnames(nvar_prof+6) = 'Land_Temperature'
vars%fldnames(nvar_prof+7) = 'Ice_Temperature'
vars%fldnames(nvar_prof+8) = 'Snow_Temperature'
vars%fldnames(nvar_prof+9) = 'Vegetation_Fraction'
vars%fldnames(nvar_prof+10) = 'Sfc_Wind_Speed'
vars%fldnames(nvar_prof+11) = 'Sfc_Wind_Direction'
vars%fldnames(nvar_prof+12) = "Lai"
vars%fldnames(nvar_prof+13) = "Soil_Moisture"
vars%fldnames(nvar_prof+14) = "Soil_Temperature"
vars%fldnames(nvar_prof+nvar_surf_real+1) = 'Land_Type_Index' ! int!!!
vars%fldnames(nvar_prof+nvar_surf_real+2) = "Vegetation_Type"
vars%fldnames(nvar_prof+nvar_surf_real+3) = "Soil_Type"

! open netcdf file and read dimensions
call nc_diag_read_init(filename, iunit)
nobs_all = nc_diag_read_get_dim(iunit,'nobs')
nchans = nc_diag_read_get_dim(iunit, 'nchans')

nobs = nobs_all / nchans

nsig = nc_diag_read_get_dim(iunit,'nsig')
nsig_plus_one = nc_diag_read_get_dim(iunit,'nsig_plus_one')

! allocate geovals structure
call ufo_geovals_setup(self, vars, nobs)

! read temperature
nval = nsig; ivar = it
allocate(field(nval, nobs_all))
call nc_diag_read_get_var(iunit, 'tvp', field)
self%geovals(ivar)%nval = nval
allocate(self%geovals(ivar)%vals(nval,nobs))
do iobs = 1, nobs
  self%geovals(ivar)%vals(:,iobs) = field(:,iobs*nchans)
enddo
deallocate(field)

! read water vapor (humidity)
nval = nsig; ivar = iwv
allocate(field(nval, nobs_all))
call nc_diag_read_get_var(iunit, 'qvp', field)
self%geovals(ivar)%nval = nval
allocate(self%geovals(ivar)%vals(nval,nobs))
do iobs = 1, nobs
  self%geovals(ivar)%vals(:,iobs) = 1000.*field(:,iobs*nchans) / (1.-field(:,iobs*nchans))
enddo
deallocate(field)

! read pressure
nval = nsig; ivar = ipr
allocate(field(nval, nobs_all))
call nc_diag_read_get_var(iunit, 'prsltmp', field)
self%geovals(ivar)%nval = nval
allocate(self%geovals(ivar)%vals(nval,nobs))
do iobs = 1, nobs
  self%geovals(ivar)%vals(:,iobs) = 10.*field(:,iobs*nchans)
enddo
deallocate(field)

! read level pressure
nval = nsig_plus_one; ivar = iprl
allocate(field(nval, nobs_all))
call nc_diag_read_get_var(iunit, 'prsitmp', field)
self%geovals(ivar)%nval = nval
allocate(self%geovals(ivar)%vals(nval,nobs))
do iobs = 1, nobs
  self%geovals(ivar)%vals(:,iobs) = 10.*field(:,iobs*nchans)
enddo
deallocate(field)

! read ozone
nval = nsig; ivar = ioz
allocate(field(nval, nobs_all))
call nc_diag_read_get_var(iunit, 'poz', field)
self%geovals(ivar)%nval = nval
allocate(self%geovals(ivar)%vals(nval,nobs))
do iobs = 1, nobs
  self%geovals(ivar)%vals(:,iobs) = field(:,iobs*nchans)
enddo
deallocate(field)

! read liquid cloud
nval = nsig; ivar = icl1
allocate(field(nval, nobs_all))
call nc_diag_read_get_var(iunit, 'cloud1', field)
self%geovals(ivar)%nval = nval
allocate(self%geovals(ivar)%vals(nval,nobs))
do iobs = 1, nobs
  self%geovals(ivar)%vals(:,iobs) = field(:,iobs*nchans)
enddo
deallocate(field)

! read ice cloud
nval = nsig; ivar = icl2
allocate(field(nval, nobs_all))
call nc_diag_read_get_var(iunit, 'cloud2', field)
self%geovals(ivar)%nval = nval
allocate(self%geovals(ivar)%vals(nval,nobs))
do iobs = 1, nobs
  self%geovals(ivar)%vals(:,iobs) = field(:,iobs*nchans)
enddo
deallocate(field)

! read surface stuff
nval = 1
allocate(field1d(nobs_all))
do ivar = nvar_prof+1, nvar_prof+nvar_surf_real
  call nc_diag_read_get_var(iunit, self%variables%fldnames(ivar), field1d)
  self%geovals(ivar)%nval = nval
  allocate(self%geovals(ivar)%vals(nval,nobs))
  do iobs = 1, nobs
    self%geovals(ivar)%vals(1,iobs) = field1d(iobs*nchans)
  enddo
enddo
deallocate(field1d)

! read surface stuff (integer)
nval = 1
allocate(field1di(nobs_all))
do ivar = nvar_prof+nvar_surf_real+1, nvar_prof+nvar_surf_real+nvar_surf_int
  call nc_diag_read_get_var(iunit, self%variables%fldnames(ivar), field1di)
  self%geovals(ivar)%nval = nval
  allocate(self%geovals(ivar)%vals(nval,nobs))
  do iobs = 1, nobs
    self%geovals(ivar)%vals(1,iobs) = field1di(iobs*nchans)
  enddo
enddo
deallocate(field1di)

self%linit = .true.

call nc_diag_read_close(filename)

end subroutine ufo_geovals_read_rad_netcdf

! ------------------------------------------------------------------------------

end module ufo_geovals_mod
