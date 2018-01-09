! (C) Copyright 2017 UCAR
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 

!> Fortran module to handle ice concentration observations

module ufo_seaicefrac_mod
  
  use ufo_obs_seaicefrac_mod
  use ufo_obs_vectors
  use ufo_vars_mod
  use ufo_locs_mod
  use ufo_geovals_mod
  use kinds
  
  implicit none
  public :: ufo_seaicefrac_eqv
  private

  ! ------------------------------------------------------------------------------
contains
  
  
! ------------------------------------------------------------------------------
subroutine ufo_seaicefrac_eqv(geovals, obss, hofx)
implicit none
type(ufo_geovals), intent(in)    :: geovals
type(ufo_obs_seaicefrac), intent(inout) :: obss
type(obs_vector),  intent(inout) :: hofx

character(len=*), parameter :: myname_="ufo_seaicefrac_eqv"
integer :: iunit

real(kind_real), allocatable :: omf(:), obs(:)

integer :: iobs, nobs
real(kind_real)  :: rmse

type(ufo_geoval) :: geoval

! Get observations from obs-structure
nobs = obss%nobs
print *, myname_, ' nobs: ', nobs, geovals%nobs, hofx%nobs

if (nobs /= geovals%nobs .or. nobs /= hofx%nobs) then
  print *, myname_, ' error: nobs inconsistent!'
  stop 6
endif

if (.not. ufo_geovals_get_var(geovals, var_seaicefrac, geoval)) then
  print *, myname_, trim(var_seaicefrac), ' doesnt exist'
  stop 5
endif

do iobs = 1, nobs
  hofx%values(iobs) = sum(geoval%vals(:,iobs))
enddo

!rmse = 0.
!do iobs = 1, nobs
!  rmse = rmse + (obs(iobs)-omf(iobs))*(obs(iobs)-omf(iobs))
!enddo
!print *, 'rmse=', sqrt(rmse/real(nobs, kind_real))

!print *, myname_, ' seaicefrac t test: max diff: ', maxval(abs(hofx%values-(obs-omf))/abs(hofx%values))


end subroutine ufo_seaicefrac_eqv

end module ufo_seaicefrac_mod
