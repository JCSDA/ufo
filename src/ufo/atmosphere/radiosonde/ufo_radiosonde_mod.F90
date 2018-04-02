! (C) Copyright 2017 UCAR
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 

!> Fortran module to handle radiosonde observations

module ufo_radiosonde_mod
  
  use ufo_obs_radiosonde_mod
  use ufo_obs_vectors
  use ufo_vars_mod
  use ufo_locs_mod
  use ufo_geovals_mod
  use kinds
  
  implicit none
  public :: ufo_radiosonde
  public :: ufo_radiosonde_t_eqv
  public :: ufo_radiosonde_settraj
  public :: ufo_radiosonde_t_eqv_tl
  public :: ufo_radiosonde_t_eqv_ad
  private
  integer, parameter :: max_string=800

!> Fortran derived type for radiosonde_t trajectory
type :: ufo_radiosonde
   type(ufo_geoval) :: prsl   !< as a vertical coordinate
   logical :: ltraj = .false. !< trajectory set?
end type ufo_radiosonde

! ------------------------------------------------------------------------------

contains
    
! ------------------------------------------------------------------------------

subroutine ufo_radiosonde_t_eqv(self, geovals, hofx, obss)

implicit none
type(ufo_radiosonde), intent(in)     :: self
type(ufo_geovals), intent(in)        :: geovals
type(obs_vector),  intent(inout)     :: hofx
type(ufo_obs_radiosonde), intent(in) :: obss

character(len=*), parameter :: myname_="ufo_radiosonde_t_eqv"
character(max_string) :: err_msg

real(kind_real), allocatable :: omf(:), obs(:)

integer :: iobs
real(kind_real) :: z, dz
real(kind_real), allocatable :: pressure(:)
type(ufo_geoval), pointer :: prsl, tv

integer, save :: run = 0

run = run + 1

print *, myname_, ' nobs: ', geovals%nobs, hofx%nobs

! check if nobs is consistent in geovals & hofx
if (geovals%nobs /= hofx%nobs) then
  write(err_msg,*) myname_, ' error: nobs inconsistent!'
  call abor1_ftn(err_msg)
endif

! check if prsl variable is in geovals and get it
if (.not. ufo_geovals_get_var(geovals, var_prsl, prsl)) then
  write(err_msg,*) myname_, trim(var_prsl), ' doesnt exist'
  call abor1_ftn(err_msg)
endif

! check if tv variable is in geovals and get it
if (.not. ufo_geovals_get_var(geovals, var_tv, tv)) then
  write(err_msg,*) myname_, trim(var_tv), ' doesnt exist'
  call abor1_ftn(err_msg)
endif

! observation of pressure (for vertical interpolation)
allocate(pressure(geovals%nobs))
pressure = obss%mass(:)%pressure

! obs operator
do iobs = 1, hofx%nobs
  z = log(pressure(iobs)/10.)
  dz = interp_weight(z, prsl%vals(:,iobs), prsl%nval)
  hofx%values(iobs) = vert_interp(tv%vals(:,iobs), tv%nval, dz)
  write(102,*)hofx%values(iobs)
enddo

deallocate(pressure)

end subroutine ufo_radiosonde_t_eqv

! ------------------------------------------------------------------------------

subroutine ufo_radiosonde_settraj(self, geovals)
implicit none
type(ufo_radiosonde), intent(inout) :: self
type(ufo_geovals), intent(in)       :: geovals

character(len=*), parameter :: myname_="ufo_radiosonde_settraj"
character(max_string) :: err_msg

type(ufo_geoval), pointer :: prsl, tv

!Check if radiosondes in geovals and get it
if (.not. ufo_geovals_get_var(geovals, var_prsl, prsl)) then
  write(err_msg,*) myname_, trim(var_prsl), ' doesnt exist'
  call abor1_ftn(err_msg)
endif

self%prsl  = prsl
self%ltraj = .true.

end subroutine ufo_radiosonde_settraj

! ------------------------------------------------------------------------------

subroutine ufo_radiosonde_t_eqv_tl(self, geovals, hofx, obss)
implicit none
type(ufo_radiosonde), intent(in)     :: self
type(ufo_geovals),    intent(in)     :: geovals
type(obs_vector),     intent(inout)  :: hofx
type(ufo_obs_radiosonde), intent(in) :: obss

character(len=*), parameter :: myname_="ufo_radiosonde_t_eqv_tl"
character(max_string) :: err_msg

real(kind_real), allocatable :: pres(:)

integer :: iobs
real(kind_real) :: z, dz
real(kind_real), allocatable :: pressure(:)
type(ufo_geoval), pointer :: tv_d

print *, myname_, ' nobs: ', geovals%nobs, hofx%nobs

! check if trajectory was set
if (.not. self%ltraj) then
  write(err_msg,*) myname_, ' trajectory wasnt set!'
  call abor1_ftn(err_msg)
endif

! check if nobs is consistent in geovals & hofx
if (geovals%nobs /= hofx%nobs) then
  write(err_msg,*) myname_, ' error: nobs inconsistent!'
  call abor1_ftn(err_msg)
endif

! check if tv variable is in geovals and get it
if (.not. ufo_geovals_get_var(geovals, var_tv, tv_d)) then
  write(err_msg,*) myname_, trim(var_tv), ' doesnt exist'
  call abor1_ftn(err_msg)
endif

! observation of pressure (for vertical interpolation)
allocate(pressure(geovals%nobs))
pressure = obss%mass(:)%pressure

! tangent linear obs operator
do iobs = 1, hofx%nobs
  z = log(pressure(iobs)/10.)
  dz = interp_weight(z, self%prsl%vals(:,iobs), self%prsl%nval)
  hofx%values(iobs) = vert_interp_tl(tv_d%vals(:,iobs), tv_d%nval, dz)
enddo

deallocate(pressure)

end subroutine ufo_radiosonde_t_eqv_tl

! ------------------------------------------------------------------------------

subroutine ufo_radiosonde_t_eqv_ad(self, geovals, hofx, obss)
implicit none
type(ufo_radiosonde), intent(in)     :: self
type(ufo_geovals),    intent(in)     :: geovals
type(obs_vector),     intent(inout)  :: hofx
type(ufo_obs_radiosonde), intent(in) :: obss

character(len=*), parameter :: myname_="ufo_radiosonde_t_eqv_ad"
character(max_string) :: err_msg

integer :: iobs
real(kind_real) :: z, dz
real(kind_real), allocatable :: pressure(:)
type(ufo_geoval), pointer :: tv_d, prsl_d

print *, myname_, ' nobs: ', geovals%nobs, hofx%nobs

! check if trajectory was set
if (.not. self%ltraj) then
  write(err_msg,*) myname_, ' trajectory wasnt set!'
  call abor1_ftn(err_msg)
endif

! check if nobs is consistent in geovals & hofx
if (geovals%nobs /= hofx%nobs) then
  write(err_msg,*) myname_, ' error: nobs inconsistent!'
  call abor1_ftn(err_msg)
endif

! check if tv variable is in geovals and get it
if (.not. ufo_geovals_get_var(geovals, var_prsl, prsl_d)) then
  write(err_msg,*) myname_, trim(var_prsl), ' doesnt exist'
  call abor1_ftn(err_msg)
endif

! check if tv variable is in geovals and get it
if (.not. ufo_geovals_get_var(geovals, var_tv, tv_d)) then
  write(err_msg,*) myname_, trim(var_tv), ' doesnt exist'
  call abor1_ftn(err_msg)
endif

! observation of pressure (for vertical interpolation)
allocate(pressure(geovals%nobs))
pressure = obss%mass(:)%pressure

! adjoint obs operator
tv_d%vals = 0.0
prsl_d%vals = 0.0
do iobs = 1, hofx%nobs
  z = log(pressure(iobs)/10.)
  dz = interp_weight(z, self%prsl%vals(:,iobs), self%prsl%nval)
  call vert_interp_ad(tv_d%vals(:,iobs), tv_d%nval, dz, hofx%values(iobs))
enddo

deallocate(pressure)

end subroutine ufo_radiosonde_t_eqv_ad

! ------------------------------------------------------------------------------
  
real(kind_real) function interp_weight(d, x, nx) 
implicit none

integer, intent(in) :: nx
real(kind_real), intent(in)    :: x(nx)
real(kind_real), intent(in)    :: d

integer ::  ix  

! Case in which x is in decreasing order
if(d>=x(1)) then
  ix = 1 
else
  ix = 1 
  do while (d < x(ix))
    ix = ix + 1 
    if (ix == nx) exit
  enddo
  ix = ix - 1 
endif
interp_weight = real(ix,kind_real) + (d-x(ix)) / (x(ix+1)-x(ix))

end function interp_weight

! ------------------------------------------------------------------------------

real(kind_real) function vert_interp(f, nsig, dz) 
implicit none

integer :: nsig
real(kind_real), intent(in)  :: f(nsig)
real(kind_real), intent(in)  :: dz

integer :: iz, izp 
real(kind_real) :: delz, delzp

iz=int(dz)
iz=max(1,min(iz,nsig))
izp=min(iz+1,nsig)

delz=dz-float(iz)
delz=max(0.,min(delz,1.))
delzp=1.-delz

vert_interp = f(iz)*delzp + f(izp)*delz

end function vert_interp

! ------------------------------------------------------------------------------

real(kind_real) function vert_interp_tl(f_tl, nsig, dz) 
implicit none

integer :: nsig
real(kind_real), intent(in)  :: f_tl(nsig)
real(kind_real), intent(in)  :: dz

integer :: iz, izp 
real(kind_real) :: delz, delzp

iz=int(dz)
iz=max(1,min(iz,nsig))
izp=min(iz+1,nsig)

delz=dz-float(iz)
delz=max(0.,min(delz,1.))
delzp=1.-delz

vert_interp_tl = f_tl(iz)*delzp + f_tl(izp)*delz

end function vert_interp_tl

! ------------------------------------------------------------------------------

subroutine vert_interp_ad(f_ad, nsig, dz, hofx_ad) 
implicit none

integer        , intent(in)    :: nsig
real(kind_real), intent(inout) :: f_ad(nsig)
real(kind_real), intent(in)    :: dz
real(kind_real), intent(in)    :: hofx_ad

integer :: iz, izp 
real(kind_real) :: delz, delzp

iz=int(dz)
iz=max(1,min(iz,nsig))
izp=min(iz+1,nsig)

delz=dz-float(iz)
delz=max(0.,min(delz,1.))
delzp=1.-delz

f_ad(iz)  = f_ad(iz)  + hofx_ad * delzp
f_ad(izp) = f_ad(izp) + hofx_ad * delz

end subroutine vert_interp_ad

end module ufo_radiosonde_mod
