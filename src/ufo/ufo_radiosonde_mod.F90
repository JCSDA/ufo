! (C) Copyright 2017 UCAR
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 

!> Fortran module to handle radiosondeentional observations

module ufo_radiosonde_mod
  
  use ufo_obs_data
  use ufo_obs_data_mod
  use ufo_obs_vectors
  use ufo_vars_mod
  use ufo_locs_mod
  use ufo_geovals_mod
  use kinds
  
  implicit none
  public :: ufo_radiosonde_t_eqv
  private

  ! ------------------------------------------------------------------------------
contains
  
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
subroutine ufo_radiosonde_t_eqv(geovals, obss, hofx)
use raobDiag_mod, only: RaobDiag
use ufo_obs_data_mod, only: Radiosonde

implicit none
type(ufo_geovals), intent(in)    :: geovals
type(obs_data),    intent(inout) :: obss
type(obs_vector),  intent(inout) :: hofx

character(len=*), parameter :: myname_="ufo_radiosonde_t_eqv"
integer :: iunit

real(kind_real), allocatable :: omf(:), obs(:)
real(kind_real), pointer :: pres(:)

integer :: iobs, nobs
real(kind_real) :: z, dz
real(kind_real) :: rmse

logical :: lfound
type(ufo_geoval) :: geoval_pr, geoval_tv
character(MAXVARLEN) :: varname

! Get observations from obs-structure
nobs = obss%nobs
allocate(obs(nobs), omf(nobs))
obss%Obspoint => Radiosonde
pres=>obss%lev
obs=Radiosonde%mass(:)%Observation
omf=Radiosonde%mass(:)%Obs_Minus_Forecast_unadjusted

print *, myname_, ' nobs: ', nobs, geovals%nobs, hofx%nobs

rmse = 0.
do iobs = 1, nobs
  rmse = rmse + (obs(iobs)-omf(iobs))*(obs(iobs)-omf(iobs))
enddo
print *, 'rmse=', sqrt(rmse/real(nobs, kind_real))

if (nobs /= geovals%nobs) then
  print *, myname_, ' error: nobs inconsistent!'
endif

varname = 'LogPressure'
lfound =  ufo_geovals_get_var(geovals, varname, geoval_pr)
if (lfound) then
  varname = 'Virtual temperature'
  lfound = ufo_geovals_get_var(geovals, varname, geoval_tv)
  if (lfound) then
    do iobs = 1, nobs
      z = log(pres(iobs)/10.)
      dz = interp_weight(z, geoval_pr%vals(:,iobs), geoval_pr%nval)
      hofx%values(iobs) = vert_interp(geoval_tv%vals(:,iobs), geoval_tv%nval, dz)
    enddo
  else
    print *, myname_, trim(varname), ' doesnt exist'
  endif
else
  print *, myname_, trim(varname), ' doesnt exist'
endif
print *, myname_, ' radiosonde t test: max diff: ', maxval(abs(hofx%values-(obs-omf))/abs(hofx%values))

nullify(pres)
deallocate(obs, omf)

end subroutine ufo_radiosonde_t_eqv

end module ufo_radiosonde_mod
