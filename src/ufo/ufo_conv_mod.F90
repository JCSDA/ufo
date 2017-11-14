! (C) Copyright 2017 UCAR
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 

!> Fortran module to handle conventional observations

module ufo_conv_mod
  
  use iso_c_binding
  use config_mod
  use ufo_obs_data
  use ufo_obs_vectors
  use ufo_vars_mod
  use ufo_locs_mod
  use ufo_geovals_mod
  use kinds
  
  implicit none
  private
  
  ! ------------------------------------------------------------------------------
  
  !> Fortran derived type for stream function observations for the QG model
  type :: ufo_obsoper
     integer :: nothing_here_yet
  end type ufo_obsoper
  
#define LISTED_TYPE ufo_obsoper
  
  !> Linked list interface - defines registry_t type
#include "linkedList_i.f"
  
  !> Global registry
  type(registry_t) :: ufo_conv_registry
  
  ! ------------------------------------------------------------------------------
contains
  ! ------------------------------------------------------------------------------
  !> Linked list implementation
#include "linkedList_c.f"
  
! ------------------------------------------------------------------------------
  
subroutine ufo_conv_setup_c(c_key_self, c_conf) bind(c,name='ufo_conv_setup_f90')
implicit none
integer(c_int), intent(inout) :: c_key_self
type(c_ptr), intent(in)    :: c_conf
    
type(ufo_obsoper), pointer :: self

call ufo_conv_registry%init()
call ufo_conv_registry%add(c_key_self)
call ufo_conv_registry%get(c_key_self, self)
    
end subroutine ufo_conv_setup_c
  
! ------------------------------------------------------------------------------
  
subroutine ufo_conv_delete_c(c_key_self) bind(c,name='ufo_conv_delete_f90')
implicit none
integer(c_int), intent(inout) :: c_key_self
    
type(ufo_obsoper), pointer :: self

call ufo_conv_registry%get(c_key_self, self)
call ufo_conv_registry%remove(c_key_self)
    
end subroutine ufo_conv_delete_c

! ------------------------------------------------------------------------------
real function interp_weight(d, x, nx) 
implicit none

integer, intent(in) :: nx
real, intent(in)    :: x(nx)
real, intent(in)    :: d

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
interp_weight = real(ix) + (d-x(ix)) / (x(ix+1)-x(ix))

end function interp_weight

! ------------------------------------------------------------------------------
real function vert_interp(f, nsig, dz) 
implicit none

integer :: nsig
real, intent(in)  :: f(nsig)
real, intent(in)  :: dz

integer :: iz, izp 
real    :: delz, delzp

iz=int(dz)
iz=max(1,min(iz,nsig))
izp=min(iz+1,nsig)

delz=dz-float(iz)
delz=max(0.,min(delz,1.))
delzp=1.-delz

vert_interp = f(iz)*delzp + f(izp)*delz

end function vert_interp
  
! ------------------------------------------------------------------------------
subroutine ufo_conv_q_eqv_c(c_key_geovals, c_key_hofx, c_bias) bind(c,name='ufo_conv_q_eqv_f90')
use nc_diag_read_mod, only: nc_diag_read_get_var
use nc_diag_read_mod, only: nc_diag_read_get_dim
use nc_diag_read_mod, only: nc_diag_read_init, nc_diag_read_close

implicit none
integer(c_int), intent(in) :: c_key_geovals
integer(c_int), intent(in) :: c_key_hofx
real(c_double), intent(in) :: c_bias
type(ufo_geovals), pointer  :: geovals
type(obs_vector), pointer :: hofx

character(128) :: filename
integer :: iunit

real(8), allocatable :: pres(:), omf(:), obs(:)
integer, allocatable :: obstype(:)

integer :: iobs, nobs
real :: z, dz

logical :: lfound
type(ufo_geoval) :: geoval
character(MAXVARLEN) :: varname

! Get pointers to geovals and hofx
call ufo_geovals_registry%get(c_key_geovals,geovals)
call ufo_obs_vect_registry%get(c_key_hofx,hofx)

! open netcdf file and read some stuff (it should be in the obs_data)
filename='Data/diag_q_01_wprofiles.nc4'
call nc_diag_read_init(filename, iunit)
nobs = nc_diag_read_get_dim(iunit,'nobs')
allocate(pres(nobs))
call nc_diag_read_get_var(iunit, "Pressure", pres)
allocate(obstype(nobs))
call nc_diag_read_get_var(iunit, "Observation_Type", obstype)
allocate(obs(nobs), omf(nobs))
call nc_diag_read_get_var(iunit, "Observation", obs)
call nc_diag_read_get_var(iunit, "Obs_Minus_Forecast_unadjusted", omf)
call nc_diag_read_close(filename)

if (nobs /= geovals%nobs) then
  print *, 'convq: error: nobs inconsistent!'
endif

print *, 'ufoconvq: nobs ', nobs, geovals%nobs, hofx%nobs
do iobs = 1, nobs
  varname = 'LogPressure'
  lfound =  ufo_geovals_get_var(geovals, iobs, varname, geoval)
  if (lfound) then
    z = log(pres(iobs)/10.) 
    dz = interp_weight(z, geoval%vals, geoval%nval)
    ! hardcoded for ships, buoys (?)
    if((obstype(iobs) > 179 .and. obstype(iobs) < 186) .or. obstype(iobs) == 199) dz=1.
    varname = 'Specific humidity'
    lfound = ufo_geovals_get_var(geovals, iobs, varname, geoval)
    if (lfound) then
      hofx%values(iobs) = vert_interp(geoval%vals, geoval%nval, dz)
!      !print *, 'convq test: interpolated q: ', hofx%values(iobs)
!      !print *, 'convq test: from gsi: ', obs(iobs) - omf(iobs)
    else
      print *, 'convq test: ', trim(varname), ' doesnt exist'
    endif
  else
    print *, 'convq test: ', trim(varname), ' doesnt exist'
  endif
enddo
print *, 'conv q test: max diff: ', maxval(abs(hofx%values-(obs-omf))/abs(hofx%values))

deallocate(obstype, obs, omf, pres)

end subroutine ufo_conv_q_eqv_c

! ------------------------------------------------------------------------------
subroutine ufo_conv_t_eqv_c(c_key_geovals, c_key_hofx, c_bias) bind(c,name='ufo_conv_t_eqv_f90')
use nc_diag_read_mod, only: nc_diag_read_get_var
use nc_diag_read_mod, only: nc_diag_read_get_dim
use nc_diag_read_mod, only: nc_diag_read_init, nc_diag_read_close

implicit none
integer(c_int), intent(in) :: c_key_geovals
integer(c_int), intent(in) :: c_key_hofx
real(c_double), intent(in) :: c_bias
type(ufo_geovals), pointer  :: geovals
type(obs_vector), pointer :: hofx

character(128) :: filename
integer :: iunit

real(8), allocatable :: pres(:), omf(:), obs(:)
integer, allocatable :: obstype(:)

integer :: iobs, nobs
real :: z, dz

logical :: lfound
type(ufo_geoval) :: geoval
character(MAXVARLEN) :: varname

! Get pointers to geovals and hofx
call ufo_geovals_registry%get(c_key_geovals,geovals)
call ufo_obs_vect_registry%get(c_key_hofx,hofx)

! open netcdf file and read some stuff (it should be in the obs_data)
filename='Data/diag_t_01_wprofiles.nc4'
call nc_diag_read_init(filename, iunit)
nobs = nc_diag_read_get_dim(iunit,'nobs')
allocate(pres(nobs))
call nc_diag_read_get_var(iunit, "Pressure", pres)
allocate(obstype(nobs))
call nc_diag_read_get_var(iunit, "Observation_Type", obstype)
allocate(obs(nobs), omf(nobs))
call nc_diag_read_get_var(iunit, "Observation", obs)
call nc_diag_read_get_var(iunit, "Obs_Minus_Forecast_unadjusted", omf)
call nc_diag_read_close(filename)

if (nobs /= geovals%nobs) then
  print *, 'convt: error: nobs inconsistent!'
endif

print *, 'ufoconvt: nobs ', nobs, geovals%nobs, hofx%nobs
do iobs = 1, nobs
  varname = 'LogPressure'
  lfound =  ufo_geovals_get_var(geovals, iobs, varname, geoval)
  if (lfound) then
    z = log(pres(iobs)/10.)
    dz = interp_weight(z, geoval%vals, geoval%nval)
    ! hardcoded for ships, buoys (?)
    if((obstype(iobs) > 179 .and. obstype(iobs) < 186) .or. obstype(iobs) == 199) dz=1.
    varname = 'Virtual temperature'
    lfound = ufo_geovals_get_var(geovals, iobs, varname, geoval)
    if (lfound) then
      hofx%values(iobs) = vert_interp(geoval%vals, geoval%nval, dz)
!      print *, 'convt test: interpolated q: ', hofx%values(iobs)
!      print *, 'convt test: from gsi: ', obs(iobs) - omf(iobs)
    else
      print *, 'convt test: ', trim(varname), ' doesnt exist'
    endif
  else
    print *, 'convt test: ', trim(varname), ' doesnt exist'
  endif
enddo
print *, 'conv t test: max diff: ', maxval(abs(hofx%values-(obs-omf))/abs(hofx%values))

deallocate(obstype, obs, omf, pres)

end subroutine ufo_conv_t_eqv_c
  
! ------------------------------------------------------------------------------
subroutine ufo_conv_u_eqv_c(c_key_geovals, c_key_hofx, c_bias) bind(c,name='ufo_conv_u_eqv_f90')
use nc_diag_read_mod, only: nc_diag_read_get_var
use nc_diag_read_mod, only: nc_diag_read_get_dim
use nc_diag_read_mod, only: nc_diag_read_init, nc_diag_read_close

implicit none
integer(c_int), intent(in) :: c_key_geovals
integer(c_int), intent(in) :: c_key_hofx
real(c_double), intent(in) :: c_bias
type(ufo_geovals), pointer  :: geovals
type(obs_vector), pointer :: hofx

character(128) :: filename
integer :: iunit

real(8), allocatable :: pres(:), omf(:), obs(:)
integer, allocatable :: obstype(:)

integer :: iobs, nobs
real :: z, dz

logical :: lfound
type(ufo_geoval) :: geoval
character(MAXVARLEN) :: varname

! Get pointers to geovals and hofx
call ufo_geovals_registry%get(c_key_geovals,geovals)
call ufo_obs_vect_registry%get(c_key_hofx,hofx)

! open netcdf file and read some stuff (it should be in the obs_data)
filename='Data/diag_uv_01_wprofiles.nc4'
call nc_diag_read_init(filename, iunit)
nobs = nc_diag_read_get_dim(iunit,'nobs')
allocate(pres(nobs))
call nc_diag_read_get_var(iunit, "Pressure", pres)
allocate(obstype(nobs))
call nc_diag_read_get_var(iunit, "Observation_Type", obstype)
allocate(obs(nobs), omf(nobs))
call nc_diag_read_get_var(iunit, "u_Observation", obs)
call nc_diag_read_get_var(iunit, "u_Obs_Minus_Forecast_unadjusted", omf)
call nc_diag_read_close(filename)

if (nobs /= geovals%nobs) then
  print *, 'conv u: error: nobs inconsistent!'
endif

print *, 'ufo conv u: nobs ', nobs, geovals%nobs, hofx%nobs
do iobs = 1, nobs
  varname = 'LogPressure'
  lfound =  ufo_geovals_get_var(geovals, iobs, varname, geoval)
  if (lfound) then
    z = log(pres(iobs)/10.)
    dz = interp_weight(z, geoval%vals, geoval%nval)
    ! hardcoded for ships, buoys (?)
!    if((obstype(iobs) > 179 .and. obstype(iobs) < 186) .or. obstype(iobs) == 199) dz=1.
    varname = 'U-wind'
    lfound = ufo_geovals_get_var(geovals, iobs, varname, geoval)
    if (lfound) then
      hofx%values(iobs) = vert_interp(geoval%vals, geoval%nval, dz)
!      print *, 'conv u test: interpolated u: ', hofx%values(iobs)
!      print *, 'conv u test: from gsi: ', obs(iobs) - omf(iobs)
    else
      print *, 'conv u test: ', trim(varname), ' doesnt exist'
    endif
  else
    print *, 'conv u test: ', trim(varname), ' doesnt exist'
  endif
enddo
print *, 'conv u test: max diff: ', maxval(abs(hofx%values-(obs-omf))/abs(hofx%values))

deallocate(obstype, obs, omf, pres)

end subroutine ufo_conv_u_eqv_c

! ------------------------------------------------------------------------------
subroutine ufo_conv_ps_eqv_c(c_key_geovals, c_key_hofx, c_bias) bind(c,name='ufo_conv_ps_eqv_f90')
use nc_diag_read_mod, only: nc_diag_read_get_var
use nc_diag_read_mod, only: nc_diag_read_get_dim
use nc_diag_read_mod, only: nc_diag_read_init, nc_diag_read_close

implicit none
integer(c_int), intent(in) :: c_key_geovals
integer(c_int), intent(in) :: c_key_hofx
real(c_double), intent(in) :: c_bias
type(ufo_geovals), pointer  :: geovals
type(obs_vector), pointer :: hofx

character(128) :: filename
integer :: iunit

real(8), allocatable :: omf(:), obs(:)

integer :: iobs, nobs
real :: z, dz

logical :: lfound
type(ufo_geoval) :: geoval
character(MAXVARLEN) :: varname

! Get pointers to geovals and hofx
call ufo_geovals_registry%get(c_key_geovals,geovals)
call ufo_obs_vect_registry%get(c_key_hofx,hofx)

! open netcdf file and read some stuff (it should be in the obs_data)
filename='Data/diag_ps_01_wprofiles.nc4'
call nc_diag_read_init(filename, iunit)
nobs = nc_diag_read_get_dim(iunit,'nobs')
allocate(obs(nobs), omf(nobs))
call nc_diag_read_get_var(iunit, "Observation", obs)
call nc_diag_read_get_var(iunit, "Obs_Minus_Forecast_unadjusted", omf)
call nc_diag_read_close(filename)

if (nobs /= geovals%nobs) then
  print *, 'conv ps: error: nobs inconsistent!'
endif

print *, 'ufo conv ps: nobs ', nobs, geovals%nobs, hofx%nobs
do iobs = 1, nobs
  varname = 'LogSurface pressure'
  lfound =  ufo_geovals_get_var(geovals, iobs, varname, geoval)
  if (lfound) then
    hofx%values(iobs) = geoval%vals(1)
!    print *, 'conv ps test: interpolated ps: ', hofx%values(iobs)
!    print *, 'conv ps test: from gsi: ', obs(iobs) - omf(iobs)
  else
    print *, 'conv ps test: ', trim(varname), ' doesnt exist'
  endif
enddo
print *, 'conv ps test: max diff: ', maxval(abs(hofx%values-(obs-omf))/abs(hofx%values))

deallocate(obs, omf)

end subroutine ufo_conv_ps_eqv_c

  
! ------------------------------------------------------------------------------
  
subroutine ufo_conv_inputs_c(c_key_self, c_key_vars) bind(c,name='ufo_conv_inputs_f90')
implicit none
integer(c_int), intent(in)    :: c_key_self
integer(c_int), intent(inout) :: c_key_vars
    
type(ufo_obsoper), pointer :: self
type(ufo_vars), pointer :: vars
    
call ufo_conv_registry%get(c_key_self, self)
call ufo_vars_registry%init()
call ufo_vars_registry%add(c_key_vars)
call ufo_vars_registry%get(c_key_vars, vars)
    
end subroutine ufo_conv_inputs_c
  
  
end module ufo_conv_mod
