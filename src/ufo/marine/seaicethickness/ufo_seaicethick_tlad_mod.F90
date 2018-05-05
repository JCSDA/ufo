! (C) Copyright 2017 UCAR
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 

!> Fortran module to handle ice concentration observations

module ufo_seaicethick_tlad_mod
  
use ufo_obs_seaicethick_mod
use ufo_obs_vectors
use ufo_vars_mod
use ufo_locs_mod
use ufo_geovals_mod
use kinds
  
implicit none
public :: ufo_seaicethick_tlad
public :: ufo_seaicethick_tlad_delete
public :: ufo_seaicethick_tlad_settraj
public :: ufo_seaicethick_tlad_eqv_tl
public :: ufo_seaicethick_tlad_eqv_ad
private
integer, parameter :: max_string=800

!> Fortran derived type for sea ice fraction observation operator
type :: ufo_seaicethick_tlad
   type(ufo_geoval) :: icethick !< ice thickness (traj)
   type(ufo_geoval) :: icefrac  !< ice fraction  (traj)
   logical :: ltraj = .false.   !< trajectory set?
end type ufo_seaicethick_tlad


! ------------------------------------------------------------------------------

contains
 
! ------------------------------------------------------------------------------

subroutine ufo_seaicethick_tlad_delete(self)
implicit none
type(ufo_seaicethick_tlad), intent(inout) :: self

self%ltraj = .false.

end subroutine ufo_seaicethick_tlad_delete

! ------------------------------------------------------------------------------

subroutine ufo_seaicethick_tlad_settraj(self, geovals)
implicit none
type(ufo_seaicethick_tlad), intent(inout) :: self
type(ufo_geovals), intent(in)       :: geovals

character(len=*), parameter :: myname_="ufo_seaicethick_tlad_settraj"
character(max_string) :: err_msg

type(ufo_geoval), pointer :: icethick, icefrac

! check if sea ice thickness variables is in geovals and get it
if (.not. ufo_geovals_get_var(geovals, var_seaicethick, icethick)) then
  write(err_msg,*) myname_, trim(var_seaicethick), ' doesnt exist'
  call abor1_ftn(err_msg)
endif

! check if sea ice fraction variables is in geovals and get it
if (.not. ufo_geovals_get_var(geovals, var_seaicefrac, icefrac)) then
  write(err_msg,*) myname_, trim(var_seaicefrac), ' doesnt exist'
  call abor1_ftn(err_msg)
endif

self%icethick = icethick
self%icefrac  = icefrac
self%ltraj    = .true.

end subroutine ufo_seaicethick_tlad_settraj


! ------------------------------------------------------------------------------

subroutine ufo_seaicethick_tlad_eqv_tl(self, geovals, hofx)
implicit none
type(ufo_seaicethick_tlad), intent(in) :: self
type(ufo_geovals), intent(in)    :: geovals
type(obs_vector),  intent(inout) :: hofx

character(len=*), parameter :: myname_="ufo_seaicethick_tlad_eqv_tl"
character(max_string) :: err_msg

integer :: iobs, icat, ncat
type(ufo_geoval), pointer :: icethick_d, icefrac_d

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

! check if sea ice fraction variable is in geovals and get it
if (.not. ufo_geovals_get_var(geovals, var_seaicefrac, icefrac_d)) then
  write(err_msg,*) myname_, trim(var_seaicefrac), ' doesnt exist'
  call abor1_ftn(err_msg)
endif

! check if sea ice thickness variable is in geovals and get it
if (.not. ufo_geovals_get_var(geovals, var_seaicethick, icethick_d)) then
  write(err_msg,*) myname_, trim(var_seaicethick), ' doesnt exist'
  call abor1_ftn(err_msg)
endif

! sea ice thickness obs operator
ncat = icefrac_d%nval
hofx%values = 0.0
do iobs = 1, hofx%nobs
   do icat = 1, ncat
     hofx%values(iobs) = hofx%values(iobs) +                                         &
                         self%icefrac%vals(icat,iobs) * icethick_d%vals(icat,iobs) / 905.0 + &
                         icefrac_d%vals(icat,iobs) * self%icethick%vals(icat,iobs) /905.0
   enddo
enddo

end subroutine ufo_seaicethick_tlad_eqv_tl

! ------------------------------------------------------------------------------

subroutine ufo_seaicethick_tlad_eqv_ad(self, geovals, hofx)
implicit none
type(ufo_seaicethick_tlad), intent(in) :: self
type(ufo_geovals), intent(inout) :: geovals
type(obs_vector),  intent(inout)    :: hofx

character(len=*), parameter :: myname_="ufo_seaicethick_tlad_eqv_ad"
character(max_string) :: err_msg

integer :: iobs, icat, ncat
type(ufo_geoval), pointer :: icefrac_d, icethick_d

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

! check if sea ice fraction variable is in geovals and get it
if (.not. ufo_geovals_get_var(geovals, var_seaicefrac, icefrac_d)) then
  write(err_msg,*) myname_, trim(var_seaicefrac), ' doesnt exist'
  call abor1_ftn(err_msg)
endif

! check if sea ice thickness variable is in geovals and get it
if (.not. ufo_geovals_get_var(geovals, var_seaicethick, icethick_d)) then
  write(err_msg,*) myname_, trim(var_seaicethick), ' doesnt exist'
  call abor1_ftn(err_msg)
endif

ncat = self%icethick%nval
if (.not.(allocated(icefrac_d%vals) .or. .not. allocated(icethick_d%vals))) then
   if (ncat < 1) then
     write(err_msg,*) myname_, ' unknown number of categories'
     call abor1_ftn(err_msg)
   endif
   if (.not. allocated(icefrac_d%vals))  allocate(icefrac_d%vals(ncat,hofx%nobs))
   if (.not. allocated(icethick_d%vals)) allocate(icethick_d%vals(ncat, hofx%nobs))
end if

! backward sea ice thickness obs operator

if (.not. allocated(icefrac_d%vals))  allocate(icefrac_d%vals(ncat,hofx%nobs))
if (.not. allocated(icethick_d%vals)) allocate(icethick_d%vals(ncat, hofx%nobs))
icethick_d%vals = 0.0
icefrac_d%vals = 0.0
do iobs = 1, hofx%nobs
   do icat = 1, ncat
      icefrac_d%vals(icat,iobs)  = icefrac_d%vals(icat,iobs) + self%icethick%vals(icat,iobs) * hofx%values(iobs) / 905.0
      icethick_d%vals(icat,iobs) = icethick_d%vals(icat,iobs) + self%icefrac%vals(icat,iobs) * hofx%values(iobs) / 905.0
   enddo
enddo

end subroutine ufo_seaicethick_tlad_eqv_ad

end module ufo_seaicethick_tlad_mod
