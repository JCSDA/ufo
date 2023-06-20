! (C) Copyright 2017-2018 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

!> Fortran coolskin module for tl/ad observation operator

module ufo_coolskin_tlad_mod

 use fckit_configuration_module, only: fckit_configuration 
 use iso_c_binding
 use kinds
 
 use ufo_geovals_mod, only: ufo_geovals, ufo_geoval, ufo_geovals_get_var
 use ufo_basis_tlad_mod, only: ufo_basis_tlad
 use ufo_vars_mod
 use obsspace_mod
 use missing_values_mod

 implicit none
 private

 integer, parameter :: max_string=800

 !> Fortran derived type for the tl/ad observation operator
 type, extends(ufo_basis_tlad), public :: ufo_coolskin_tlad
 private
  integer          :: nlocs           !< Local number of obs
  real(c_double)   :: r_miss_val      !< Missing value flag  
  real (kind=kind_real), allocatable :: jac(:,:)   !< Jacobian     [6*nobs]
 contains
  procedure :: setup  => ufo_coolskin_tlad_setup
  procedure :: delete  => ufo_coolskin_tlad_delete
  procedure :: settraj => ufo_coolskin_tlad_settraj
  procedure :: simobs_tl  => ufo_coolskin_simobs_tl
  procedure :: simobs_ad  => ufo_coolskin_simobs_ad
 end type ufo_coolskin_tlad

contains

! ------------------------------------------------------------------------------
subroutine ufo_coolskin_tlad_setup(self, f_conf)
class(ufo_coolskin_tlad), intent(inout) :: self
type(fckit_configuration), intent(in)   :: f_conf

! Set missing flag
self%r_miss_val = missing_value(self%r_miss_val)

end subroutine ufo_coolskin_tlad_setup

! ------------------------------------------------------------------------------
subroutine ufo_coolskin_tlad_delete(self)
class(ufo_coolskin_tlad), intent(inout) :: self

if (allocated(self%jac)) deallocate(self%jac)
self%ltraj = .false.

end subroutine ufo_coolskin_tlad_delete

! ------------------------------------------------------------------------------
subroutine ufo_coolskin_tlad_settraj(self, geovals, obss)
use ufo_coolskin_sim_mod

implicit none
class(ufo_coolskin_tlad), intent(inout) :: self
type(ufo_geovals),        intent(in)    :: geovals
type(c_ptr), value,       intent(in)    :: obss

character(len=*), parameter :: myname_="ufo_coolskin_tlad_settraj"
type(ufo_geoval), pointer :: S_ns,H_I,H_s,R_nl,Td,u

integer        :: iobs
self%nlocs = obsspace_get_nlocs(obss)

! Get trajectory geovals (state at which the Jacobian is computed)
call ufo_geovals_get_var(geovals, var_ocn_sst, Td)
call ufo_geovals_get_var(geovals, var_sw_rad , R_nl )
call ufo_geovals_get_var(geovals, var_latent_heat , H_I )
call ufo_geovals_get_var(geovals, var_sens_heat , H_s )
call ufo_geovals_get_var(geovals, var_lw_rad , S_ns )
call ufo_geovals_get_var(geovals, var_sea_fric_vel , u )

self%ltraj  = .true.

! Initialize and compute traj dependent Jacobian
allocate(self%jac(6,self%nlocs))
do iobs = 1, self%nlocs
   ! check for missing values
   ! (the atmospheric fields *shouldn't* be masked, so dont bother checking)
   if (Td%vals(1, iobs) == self%r_miss_val .or. &
       u%vals(1,iobs) == self%r_miss_val) then
     self%jac(:,iobs) = 0.0
     cycle
   end if

   call ufo_coolskin_jac(self%jac(:,iobs),&
                         S_ns%vals(1,iobs),&
                         H_I%vals(1,iobs),&
                         H_s%vals(1,iobs),&
                         R_nl%vals(1,iobs),&
                         Td%vals(1,iobs),&
                         u%vals(1,iobs))
enddo

end subroutine ufo_coolskin_tlad_settraj

! ------------------------------------------------------------------------------
subroutine ufo_coolskin_simobs_tl(self, geovals, hofx, obss)

implicit none
class(ufo_coolskin_tlad), intent(in)   :: self
type(ufo_geovals),       intent(in)    :: geovals
real(c_double),          intent(inout) :: hofx(:)
type(c_ptr), value,      intent(in)    :: obss

character(len=*), parameter :: myname_="ufo_coolskin_simobs_tl"
character(max_string) :: err_msg
integer :: iobs, nobs
type(ufo_geoval), pointer :: S_ns,H_I,H_s,R_nl,Td,u


! Get geovals increments
call ufo_geovals_get_var(geovals, var_ocn_sst, Td)
call ufo_geovals_get_var(geovals, var_sw_rad , R_nl )
call ufo_geovals_get_var(geovals, var_latent_heat , H_I )
call ufo_geovals_get_var(geovals, var_sens_heat , H_s )
call ufo_geovals_get_var(geovals, var_lw_rad , S_ns )
call ufo_geovals_get_var(geovals, var_sea_fric_vel , u )

! check if trajectory was set
if (.not. self%ltraj) then
  write(err_msg,*) myname_, ' trajectory wasnt set!'
  call abor1_ftn(err_msg)
endif

! check if nobs is consistent in geovals & hofx
nobs = self%nlocs

if (geovals%nlocs /= nobs) then
  write(err_msg,*) myname_, ' error: nobs inconsistent!'
  call abor1_ftn(err_msg)
endif


! Perturbation coolskin obs operator
hofx = 0.0
do iobs = 1, self%nlocs
   ! check for missing values
   ! (the atmospheric fields *shouldn't* be masked, so dont bother checking)
   if (Td%vals(1,iobs) == self%r_miss_val .or. &
        u%vals(1,iobs) == self%r_miss_val) then
          hofx(iobs) = self%r_miss_val
      cycle
    end if

   hofx(iobs) = self%jac(1,iobs)*S_ns%vals(1,iobs) + &
                self%jac(2,iobs)*H_I%vals(1,iobs) + &
                self%jac(3,iobs)*H_s%vals(1,iobs) + &
                self%jac(4,iobs)*R_nl%vals(1,iobs) + &
                self%jac(5,iobs)*Td%vals(1,iobs) + &
                self%jac(6,iobs)*u%vals(1,iobs)
enddo

end subroutine ufo_coolskin_simobs_tl

! ------------------------------------------------------------------------------
subroutine ufo_coolskin_simobs_ad(self, geovals, hofx, obss)

implicit none
class(ufo_coolskin_tlad), intent(in)   :: self
type(ufo_geovals),       intent(inout) :: geovals
real(c_double),          intent(in)    :: hofx(:)
type(c_ptr), value,      intent(in)    :: obss

character(len=*), parameter :: myname_="ufo_coolskin_simobs_ad"
character(max_string) :: err_msg

integer :: iobs, nobs

type(ufo_geoval), pointer :: S_ns, H_I, H_s, R_nl, Td, u


! check if trajectory was set
if (.not. self%ltraj) then
  write(err_msg,*) myname_, ' trajectory wasnt set!'
  call abor1_ftn(err_msg)
endif

! check if nobs is consistent in geovals & hofx
nobs = self%nlocs
if (geovals%nlocs /= nobs) then
  write(err_msg,*) myname_, ' error: nobs inconsistent!'
  call abor1_ftn(err_msg)
endif

if (.not. geovals%linit ) geovals%linit=.true.

! Associate geovals pointers
call ufo_geovals_get_var(geovals, var_ocn_sst, Td)
call ufo_geovals_get_var(geovals, var_sw_rad , R_nl )
call ufo_geovals_get_var(geovals, var_latent_heat , H_I )
call ufo_geovals_get_var(geovals, var_sens_heat , H_s )
call ufo_geovals_get_var(geovals, var_lw_rad , S_ns )
call ufo_geovals_get_var(geovals, var_sea_fric_vel , u )

! Apply adjoint obs operator
do iobs = 1, nobs
   if (hofx(iobs)/=self%r_miss_val) then
      S_ns%vals(1,iobs) = S_ns%vals(1,iobs)  + self%jac(1,iobs)*hofx(iobs)
      H_I%vals(1,iobs)  = H_I%vals(1,iobs)   + self%jac(2,iobs)*hofx(iobs)
      H_s%vals(1,iobs)  = H_s%vals(1,iobs)   + self%jac(3,iobs)*hofx(iobs)
      R_nl%vals(1,iobs) = R_nl%vals(1,iobs)  + self%jac(4,iobs)*hofx(iobs)
      Td%vals(1,iobs)   = Td%vals(1,iobs)    + self%jac(5,iobs)*hofx(iobs)
      u%vals(1,iobs)    = u%vals(1,iobs)     + self%jac(6,iobs)*hofx(iobs)
   end if
enddo

end subroutine ufo_coolskin_simobs_ad

!--------------------------------------------------

end module ufo_coolskin_tlad_mod
