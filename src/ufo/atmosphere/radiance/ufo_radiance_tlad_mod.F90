! (C) Copyright 2017 UCAR
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 

!> Fortran module to handle radiance observations

module ufo_radiance_tlad_mod
  
  use ufo_obs_radiance_mod
  use ufo_obs_vectors
  use ufo_vars_mod
  use ufo_locs_mod
  use ufo_geovals_mod
  use kinds  

  use crtm_module

  implicit none
  public :: ufo_radiance_tlad
  public :: ufo_radiance_tlad_delete
  public :: ufo_radiance_tlad_settraj
  public :: ufo_radiance_tlad_eqv_tl
  public :: ufo_radiance_tlad_eqv_ad
  private
  integer, parameter :: max_string=800

!> Fortran derived type for radiance trajectory
type :: ufo_radiance_tlad
   logical :: ltraj = .false. !< trajectory set?
end type ufo_radiance_tlad
  
! ------------------------------------------------------------------------------

contains
  
! ------------------------------------------------------------------------------

subroutine ufo_radiance_tlad_delete(self)
implicit none
type(ufo_radiance_tlad), intent(inout)  :: self

self%ltraj = .false.

end subroutine ufo_radiance_tlad_delete

! ------------------------------------------------------------------------------
 
subroutine ufo_radiance_tlad_settraj(self, geovals)
implicit none
type(ufo_radiance_tlad), intent(inout) :: self
type(ufo_geovals), intent(in)       :: geovals

character(len=*), parameter :: myname_="ufo_radiance_tlad_settraj"
character(max_string) :: err_msg

! Nothing here yet

self%ltraj = .false. !.true.

end subroutine ufo_radiance_tlad_settraj

! ------------------------------------------------------------------------------

subroutine ufo_radiance_tlad_eqv_tl(self, geovals, hofx, obss)
implicit none
type(ufo_radiance_tlad), intent(in)     :: self
type(ufo_geovals),    intent(in)     :: geovals
type(obs_vector),     intent(inout)  :: hofx
type(ufo_obs_radiance), intent(in) :: obss

character(len=*), parameter :: myname_="ufo_radiance_tlad_eqv_tl"
character(max_string) :: err_msg

! Nothing here yet

end subroutine ufo_radiance_tlad_eqv_tl

! ------------------------------------------------------------------------------

subroutine ufo_radiance_tlad_eqv_ad(self, geovals, hofx, obss)
implicit none
type(ufo_radiance_tlad), intent(in)     :: self
type(ufo_geovals),    intent(in)     :: geovals
type(obs_vector),     intent(inout)  :: hofx
type(ufo_obs_radiance), intent(in) :: obss

character(len=*), parameter :: myname_="ufo_radiance_tlad_eqv_ad"
character(max_string) :: err_msg

! Nothing here yet

end subroutine ufo_radiance_tlad_eqv_ad

! ------------------------------------------------------------------------------
 
end module ufo_radiance_tlad_mod
