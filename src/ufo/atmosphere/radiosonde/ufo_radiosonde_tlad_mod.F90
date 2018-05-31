! (C) Copyright 2017 UCAR
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 

!> Fortran module to handle radiosonde observations

module ufo_radiosonde_tlad_mod
  
  use ioda_obsdb_mod
  use ioda_obs_vectors
  use ufo_vars_mod
  use ioda_locs_mod
  use ufo_geovals_mod
  use kinds
  use vert_interp_mod

  implicit none
  public :: ufo_radiosonde_tlad
  public :: ufo_radiosonde_tlad_settraj
  public :: ufo_radiosonde_tlad_t_eqv_tl
  public :: ufo_radiosonde_tlad_t_eqv_ad
  public :: ufo_radiosonde_tlad_delete
  private
  integer, parameter :: max_string=800

!> Fortran derived type for radiosonde_t trajectory
type :: ufo_radiosonde_tlad
   integer :: nval, nobs
   logical :: ltraj = .false. !< trajectory set?
   real(kind_real), allocatable :: wf(:)
   integer, allocatable :: wi(:)
end type ufo_radiosonde_tlad

! ------------------------------------------------------------------------------

contains
    
! ------------------------------------------------------------------------------

subroutine ufo_radiosonde_tlad_settraj(self, geovals, obss)
implicit none
type(ufo_radiosonde_tlad), intent(inout) :: self
type(ufo_geovals),         intent(in)    :: geovals
type(ioda_obsdb),          intent(in)    :: obss

character(len=*), parameter :: myname_="ufo_radiosonde_tlad_settraj"
character(max_string) :: err_msg

type(obs_vector) :: pressure
type(ufo_geoval), pointer :: prsl
integer :: iobs

!Check if radiosondes in geovals and get it
if (.not. ufo_geovals_get_var(geovals, var_prsl, prsl)) then
  write(err_msg,*) myname_, trim(var_prsl), ' doesnt exist'
  call abor1_ftn(err_msg)
endif

!Make sure nothing already allocated
call ufo_radiosonde_tlad_delete(self)

!Keep copy of dimensions
self%nobs = prsl%nobs
self%nval = prsl%nval

allocate(self%wi(self%nobs))
allocate(self%wf(self%nobs))

! observation of pressure (for vertical interpolation)
call ioda_obsdb_var_to_ovec(obss, pressure, "Pressure")

! compute interpolation weights
do iobs = 1, self%nobs
  call vert_interp_weights(self%nval,log(pressure%values(iobs)/10.),prsl%vals(:,iobs),self%wi(iobs),self%wf(iobs))
enddo

self%ltraj = .true.

end subroutine ufo_radiosonde_tlad_settraj

! ------------------------------------------------------------------------------

subroutine ufo_radiosonde_tlad_t_eqv_tl(self, geovals, hofx, obss)
implicit none
type(ufo_radiosonde_tlad), intent(in)     :: self
type(ufo_geovals),         intent(in)     :: geovals
type(obs_vector),          intent(inout)  :: hofx
type(ioda_obsdb),          intent(in)     :: obss

character(len=*), parameter :: myname_="ufo_radiosonde_tlad_t_eqv_tl"
character(max_string) :: err_msg

integer :: iobs
type(ufo_geoval), pointer :: tv_d

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

! tangent linear obs operator (linear)
do iobs = 1, hofx%nobs
  call vert_interp_apply_tl(tv_d%nval, tv_d%vals(:,iobs), hofx%values(iobs), self%wi(iobs), self%wf(iobs))
enddo

end subroutine ufo_radiosonde_tlad_t_eqv_tl

! ------------------------------------------------------------------------------

subroutine ufo_radiosonde_tlad_t_eqv_ad(self, geovals, hofx, obss)
implicit none
type(ufo_radiosonde_tlad), intent(in)     :: self
type(ufo_geovals),         intent(inout)  :: geovals
type(obs_vector),          intent(in)     :: hofx
type(ioda_obsdb),          intent(in)     :: obss

character(len=*), parameter :: myname_="ufo_radiosonde_tlad_t_eqv_ad"
character(max_string) :: err_msg

integer :: iobs
type(ufo_geoval), pointer :: tv_d, prsl_d

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

! allocate if not yet allocated
if (.not. allocated(tv_d%vals)) then
   tv_d%nobs = self%nobs
   tv_d%nval = self%nval
   allocate(tv_d%vals(tv_d%nval,tv_d%nobs))
endif
if (.not. allocated(prsl_d%vals)) then
   prsl_d%nobs = self%nobs
   prsl_d%nval = self%nval
   allocate(prsl_d%vals(prsl_d%nval,prsl_d%nobs))
endif
if (.not. geovals%linit ) geovals%linit=.true.

! adjoint obs operator
tv_d%vals = 0.0_kind_real
prsl_d%vals = 0.0_kind_real
do iobs = 1, hofx%nobs
  call vert_interp_apply_ad(tv_d%nval, tv_d%vals(:,iobs), hofx%values(iobs), self%wi(iobs), self%wf(iobs))
enddo

end subroutine ufo_radiosonde_tlad_t_eqv_ad

! ------------------------------------------------------------------------------

subroutine ufo_radiosonde_tlad_delete(self)
implicit none
type(ufo_radiosonde_tlad), intent(inout) :: self

character(len=*), parameter :: myname_="ufo_radiosonde_tlad_delete"

self%nval = 0
self%nobs = 0
if (allocated(self%wi)) deallocate(self%wi)
if (allocated(self%wf)) deallocate(self%wf)
self%ltraj = .false.

end subroutine ufo_radiosonde_tlad_delete

! ------------------------------------------------------------------------------

end module ufo_radiosonde_tlad_mod
