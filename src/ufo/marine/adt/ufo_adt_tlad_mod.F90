! (C) Copyright 2017 UCAR
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 

!> Fortran module to handle adt observations

module ufo_adt_tlad_mod
  
use ioda_obs_adt_mod
use ioda_obs_vectors
use ufo_vars_mod
use ioda_locs_mod
use ufo_geovals_mod
use kinds
  
implicit none
public :: ufo_adt_tlad
public :: ufo_adt_tlad_delete
public :: ufo_adt_tlad_settraj
public :: ufo_adt_tlad_eqv_tl
public :: ufo_adt_tlad_eqv_ad
private
integer, parameter :: max_string=800

!> Fortran derived type for adt observation operator
type :: ufo_adt_tlad
   type(ufo_geoval) :: adt !< adt (traj)
   logical :: ltraj = .false.   !< trajectory set?
end type ufo_adt_tlad


! ------------------------------------------------------------------------------

contains
 
! ------------------------------------------------------------------------------

subroutine ufo_adt_tlad_delete(self)
implicit none
type(ufo_adt_tlad), intent(inout) :: self

self%ltraj = .false.

end subroutine ufo_adt_tlad_delete

! ------------------------------------------------------------------------------

subroutine ufo_adt_tlad_settraj(self, geovals)
implicit none
type(ufo_adt_tlad), intent(inout) :: self
type(ufo_geovals), intent(in)       :: geovals

character(len=*), parameter :: myname_="ufo_adt_tlad_settraj"
character(max_string) :: err_msg

type(ufo_geoval), pointer :: adt 

! check if adt variables is in geovals and get it
if (.not. ufo_geovals_get_var(geovals, var_abs_topo, adt)) then
  write(err_msg,*) myname_, trim(var_abs_topo), ' doesnt exist'
  call abor1_ftn(err_msg)
endif

self%adt = adt 
self%ltraj    = .true.

end subroutine ufo_adt_tlad_settraj


! ------------------------------------------------------------------------------

subroutine ufo_adt_tlad_eqv_tl(self, geovals, hofx)
implicit none
type(ufo_adt_tlad), intent(in) :: self
type(ufo_geovals), intent(in)    :: geovals
type(obs_vector),  intent(inout) :: hofx

character(len=*), parameter :: myname_="ufo_adt_tlad_eqv_tl"
character(max_string) :: err_msg

integer :: iobs, icat, ncat
type(ufo_geoval), pointer :: adt_d 

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

! check if adt variable is in geovals and get it
if (.not. ufo_geovals_get_var(geovals, var_abs_topo, adt_d)) then
  write(err_msg,*) myname_, trim(var_abs_topo), ' doesnt exist'
  call abor1_ftn(err_msg)
endif

! adt obs operator
hofx%values = 0.0
do iobs = 1, hofx%nobs
     hofx%values(iobs) = adt_d%vals(1,iobs)
enddo

print *,'tl hx done!'

end subroutine ufo_adt_tlad_eqv_tl

! ------------------------------------------------------------------------------

subroutine ufo_adt_tlad_eqv_ad(self, geovals, hofx)
implicit none
type(ufo_adt_tlad), intent(in) :: self
type(ufo_geovals), intent(inout) :: geovals
type(obs_vector),  intent(inout)    :: hofx

character(len=*), parameter :: myname_="ufo_adt_tlad_eqv_ad"
character(max_string) :: err_msg

integer :: iobs, icat, ncat
type(ufo_geoval), pointer :: adt_d

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

if (.not. geovals%linit ) geovals%linit=.true.

! check if adt variable is in geovals and get it
if (.not. ufo_geovals_get_var(geovals, var_abs_topo, adt_d)) then
  write(err_msg,*) myname_, trim(var_abs_topo), ' doesnt exist'
  call abor1_ftn(err_msg)
endif

! backward adt obs operator

if (.not. allocated(adt_d%vals))  allocate(adt_d%vals(1,hofx%nobs))
adt_d%vals = 0.0
do iobs = 1, hofx%nobs
      adt_d%vals(1,iobs) =  adt_d%vals(1,iobs) + hofx%values(iobs) 
enddo

print *,'ad hx done!'

end subroutine ufo_adt_tlad_eqv_ad

end module ufo_adt_tlad_mod
