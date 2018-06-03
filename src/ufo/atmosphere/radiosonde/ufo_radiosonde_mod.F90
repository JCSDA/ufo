! (C) Copyright 2017 UCAR
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 

!> Fortran module to handle radiosonde observations

module ufo_radiosonde_mod
  
  use ioda_obsdb_mod
  use ioda_obs_vectors
  use ufo_vars_mod
  use ioda_locs_mod
  use ufo_geovals_mod
  use kinds
  use vert_interp_mod

  implicit none
  public :: ufo_radiosonde
  public :: ufo_radiosonde_t_eqv
  private
  integer, parameter :: max_string=800

!> Fortran derived type for radiosonde_t trajectory
type :: ufo_radiosonde
end type ufo_radiosonde

! ------------------------------------------------------------------------------

contains
    
! ------------------------------------------------------------------------------

subroutine ufo_radiosonde_t_eqv(self, geovals, hofx, obss)

implicit none
type(ufo_radiosonde), intent(in)  :: self
type(ufo_geovals), intent(in)     :: geovals
type(obs_vector),  intent(inout)  :: hofx
type(ioda_obsdb),  intent(in)     :: obss

character(len=*), parameter :: myname_="ufo_radiosonde_t_eqv"
character(max_string) :: err_msg

integer :: iobs
real(kind_real) :: wf
integer :: wi
type(obs_vector) :: pressure
type(ufo_geoval), pointer :: prsl, tv

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
call ioda_obsvec_setup(pressure, obss%nobs)
call ioda_obsdb_var_to_ovec(obss, pressure, "Pressure")

! obs operator
do iobs = 1, hofx%nobs
  call vert_interp_weights(prsl%nval,log(pressure%values(iobs)/10.),prsl%vals(:,iobs),wi,wf)
  call vert_interp_apply(tv%nval, tv%vals(:,iobs), hofx%values(iobs), wi, wf)
enddo
call ioda_obsvec_delete(pressure)

end subroutine ufo_radiosonde_t_eqv

! ------------------------------------------------------------------------------

end module ufo_radiosonde_mod
