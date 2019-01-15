! (C) Copyright 2017-2018 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

!> Fortran module for seaicethickness tl/ad observation operator

module ufo_seaicethickness_tlad_mod

 use iso_c_binding
 use config_mod
 use kinds

 use ufo_geovals_mod, only: ufo_geovals, ufo_geoval, ufo_geovals_get_var
 use ufo_basis_tlad_mod, only: ufo_basis_tlad
 use ufo_vars_mod
 use obsspace_mod

 implicit none
 private

 integer, parameter :: max_string=800

 !> Fortran derived type for the tl/ad observation operator
 type, extends(ufo_basis_tlad), public :: ufo_seaicethickness_tlad
 private
  type(ufo_geoval) :: icethick !< ice thickness (traj)
  type(ufo_geoval) :: icefrac  !< ice fraction  (traj)
 contains
  procedure :: setup  => ufo_seaicethickness_tlad_setup
  procedure :: delete  => ufo_seaicethickness_tlad_delete
  procedure :: settraj => ufo_seaicethickness_tlad_settraj
  procedure :: simobs_tl  => ufo_seaicethickness_simobs_tl
  procedure :: simobs_ad  => ufo_seaicethickness_simobs_ad
 end type ufo_seaicethickness_tlad

contains

! ------------------------------------------------------------------------------
subroutine ufo_seaicethickness_tlad_setup(self, c_conf)
implicit none
class(ufo_seaicethickness_tlad), intent(inout) :: self
type(c_ptr),              intent(in)    :: c_conf

end subroutine ufo_seaicethickness_tlad_setup

! ------------------------------------------------------------------------------
subroutine ufo_seaicethickness_tlad_delete(self)
implicit none
class(ufo_seaicethickness_tlad), intent(inout) :: self

self%ltraj = .false.

end subroutine ufo_seaicethickness_tlad_delete

! ------------------------------------------------------------------------------
subroutine ufo_seaicethickness_tlad_settraj(self, geovals, obss)
implicit none
class(ufo_seaicethickness_tlad), intent(inout) :: self
type(ufo_geovals),       intent(in)    :: geovals
type(c_ptr), value,      intent(in)    :: obss

character(len=*), parameter :: myname_="ufo_seaicethick_tlad_settraj"
character(max_string) :: err_msg

type(ufo_geoval), pointer :: icethick, icefrac

! check if sea ice thickness variables is in geovals and get it
call ufo_geovals_get_var(geovals, var_seaicethick, icethick)

! check if sea ice fraction variables is in geovals and get it
call ufo_geovals_get_var(geovals, var_seaicefrac, icefrac)

self%icethick = icethick
self%icefrac  = icefrac
self%ltraj    = .true.

end subroutine ufo_seaicethickness_tlad_settraj

! ------------------------------------------------------------------------------
subroutine ufo_seaicethickness_simobs_tl(self, geovals, hofx, obss)
implicit none
class(ufo_seaicethickness_tlad), intent(in)    :: self
type(ufo_geovals),       intent(in)    :: geovals
real(c_double),          intent(inout) :: hofx(:)
type(c_ptr), value,      intent(in)    :: obss

character(len=*), parameter :: myname_="ufo_seaicethick_simobs_tl"
character(max_string) :: err_msg

integer :: iobs, icat, ncat
type(ufo_geoval), pointer :: icethick_d, icefrac_d

print *, myname_, ' nobs: ', geovals%nobs, size(hofx,1)

! check if trajectory was set
if (.not. self%ltraj) then
  write(err_msg,*) myname_, ' trajectory wasnt set!'
  call abor1_ftn(err_msg)
endif

! check if nobs is consistent in geovals & hofx
if (geovals%nobs /= size(hofx,1)) then
  write(err_msg,*) myname_, ' error: nobs inconsistent!'
  call abor1_ftn(err_msg)
endif

! check if sea ice fraction variable is in geovals and get it
call ufo_geovals_get_var(geovals, var_seaicefrac, icefrac_d)

! check if sea ice thickness variable is in geovals and get it
call ufo_geovals_get_var(geovals, var_seaicethick, icethick_d)

! sea ice thickness obs operator
ncat = icefrac_d%nval
hofx = 0.0
do iobs = 1, size(hofx,1)
   do icat = 1, ncat
     hofx(iobs) = hofx(iobs) +                                         &
                  self%icefrac%vals(icat,iobs) * icethick_d%vals(icat,iobs) / 905.0 + &
                  icefrac_d%vals(icat,iobs) * self%icethick%vals(icat,iobs) /905.0
   enddo
enddo


end subroutine ufo_seaicethickness_simobs_tl

! ------------------------------------------------------------------------------
subroutine ufo_seaicethickness_simobs_ad(self, geovals, hofx, obss)
implicit none
class(ufo_seaicethickness_tlad), intent(in)    :: self
type(ufo_geovals),       intent(inout) :: geovals
real(c_double),          intent(in)    :: hofx(:)
type(c_ptr), value,      intent(in)    :: obss


character(len=*), parameter :: myname_="ufo_seaicethick_simobs_ad"
character(max_string) :: err_msg

integer :: iobs, icat, ncat
type(ufo_geoval), pointer :: icefrac_d, icethick_d

! check if trajectory was set
if (.not. self%ltraj) then
  write(err_msg,*) myname_, ' trajectory wasnt set!'
  call abor1_ftn(err_msg)
endif

! check if nobs is consistent in geovals & hofx
if (geovals%nobs /= size(hofx,1)) then
  write(err_msg,*) myname_, ' error: nobs inconsistent!'
  call abor1_ftn(err_msg)
endif

if (.not. geovals%linit ) geovals%linit=.true.

! check if sea ice fraction variable is in geovals and get it
call ufo_geovals_get_var(geovals, var_seaicefrac, icefrac_d)

! check if sea ice thickness variable is in geovals and get it
call ufo_geovals_get_var(geovals, var_seaicethick, icethick_d)

ncat = self%icethick%nval
if (.not.(allocated(icefrac_d%vals) .or. .not. allocated(icethick_d%vals))) then
   if (ncat < 1) then
     write(err_msg,*) myname_, ' unknown number of categories'
     call abor1_ftn(err_msg)
   endif
   if (.not. allocated(icefrac_d%vals))  allocate(icefrac_d%vals(ncat,size(hofx,1)))
   if (.not. allocated(icethick_d%vals)) allocate(icethick_d%vals(ncat, size(hofx,1)))
end if

!print *, 'in ad: hofx=', hofx

! backward sea ice thickness obs operator

print *,'ncat=',ncat
if (.not. allocated(icefrac_d%vals))  allocate(icefrac_d%vals(ncat,size(hofx,1)))
if (.not. allocated(icethick_d%vals)) allocate(icethick_d%vals(ncat, size(hofx,1)))

icethick_d%vals = 0.0
icefrac_d%vals = 0.0
do iobs = 1, size(hofx,1)
   do icat = 1, ncat
      icefrac_d%vals(icat,iobs)  = icefrac_d%vals(icat,iobs) + self%icethick%vals(icat,iobs) * hofx(iobs) / 905.0
      icethick_d%vals(icat,iobs) = icethick_d%vals(icat,iobs) + self%icefrac%vals(icat,iobs) * hofx(iobs) / 905.0
      !print *, 'in ad: thick=', icethick_d%vals(:,iobs)
      !print *, 'in ad: frac=', icefrac_d%vals(:,iobs)
   enddo
enddo

end subroutine ufo_seaicethickness_simobs_ad

! ------------------------------------------------------------------------------

end module ufo_seaicethickness_tlad_mod
