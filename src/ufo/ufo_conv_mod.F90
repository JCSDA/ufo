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

real(8), allocatable :: tvflag(:), pres(:), omf(:), obs(:)
real(8), allocatable :: pres_raob(:), omf_raob(:), obs_raob(:)
integer, allocatable :: obstype(:)

integer :: iobs, nobs, iobs_raob, nobs_raob
real :: z, dz
real(8) :: rmse
integer, parameter :: raobtype = 120
logical :: lfound
type(ufo_geoval) :: geoval_pr, geoval_tv
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
allocate(obstype(nobs),tvflag(nobs))
call nc_diag_read_get_var(iunit, "Observation_Type", obstype)
call nc_diag_read_get_var(iunit, "Setup_QC_Mark", tvflag)
allocate(obs(nobs), omf(nobs))
call nc_diag_read_get_var(iunit, "Observation", obs)
call nc_diag_read_get_var(iunit, "Obs_Minus_Forecast_unadjusted", omf)
call nc_diag_read_close(filename)

nobs_raob = count(obstype == raobtype .and. tvflag == 0)
allocate(pres_raob(nobs_raob), obs_raob(nobs_raob))
allocate(omf_raob(nobs_raob))
iobs_raob = 1
do iobs = 1, nobs
  if (obstype(iobs) == raobtype .and. tvflag(iobs) == 0) then
    !print *, iobs, obstype(iobs), tvflag(iobs)
    pres_raob(iobs_raob) = pres(iobs)
    obs_raob(iobs_raob) = obs(iobs)
    omf_raob(iobs_raob) = omf(iobs)
    iobs_raob = iobs_raob + 1
  endif
enddo

rmse = 0.
do iobs = 1, nobs_raob
  rmse = rmse + (obs_raob(iobs)-omf_raob(iobs))*(obs_raob(iobs)-omf_raob(iobs))
enddo
print *, 'rmse=', sqrt(rmse/real(nobs_raob, 8))

if (nobs_raob /= geovals%nobs) then
  print *, 'convt: error: nobs inconsistent!'
endif

varname = 'LogPressure'
lfound =  ufo_geovals_get_var(geovals, varname, geoval_pr)
if (lfound) then
  varname = 'Virtual temperature'
  lfound = ufo_geovals_get_var(geovals, varname, geoval_tv)
  if (lfound) then
    do iobs = 1, nobs_raob
      z = log(pres_raob(iobs)/10.)
      dz = interp_weight(z, geoval_pr%vals(:,iobs), geoval_pr%nval)
      hofx%values(iobs) = vert_interp(geoval_tv%vals(:,iobs), geoval_tv%nval, dz)
    enddo
  else
    print *, 'convt test: ', trim(varname), ' doesnt exist'
  endif
else
  print *, 'convt test: ', trim(varname), ' doesnt exist'
endif
print *, 'conv t test: max diff: ', maxval(abs(hofx%values-(obs_raob-omf_raob))/abs(hofx%values))

deallocate(obstype, obs, omf, pres, tvflag)
deallocate(obs_raob, omf_raob, pres_raob)

end subroutine ufo_conv_t_eqv_c

  
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
