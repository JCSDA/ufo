! (C) Copyright 2017-2018 UCAR
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 

!> Fortran module to handle radiance observations

module ufo_radiance_tlad_mod
  
  use ioda_obsdb_mod
  use ioda_obs_vectors
  use ufo_vars_mod
  use ioda_locs_mod
  use ufo_geovals_mod
  use kinds  
  use ufo_basis_tlad_mod, only: ufo_basis_tlad

  use crtm_module

  implicit none

  public :: ufo_radiance_tlad
  private
  integer, parameter :: max_string=800

  !> Fortran derived type for radiance trajectory
  type, extends(ufo_basis_tlad) :: ufo_radiance_tlad
  contains
    procedure :: delete  => ufo_radiance_tlad_delete
    procedure :: settraj => ufo_radiance_tlad_settraj 
    procedure :: eqv_tl  => ufo_radiance_tlad_eqv_tl
    procedure :: eqv_ad  => ufo_radiance_tlad_eqv_ad
  end type ufo_radiance_tlad
  
contains

! ------------------------------------------------------------------------------

subroutine ufo_radiance_tlad_delete(self)
  implicit none
  class(ufo_radiance_tlad), intent(inout)  :: self

  self%ltraj = .false.

end subroutine ufo_radiance_tlad_delete

! ------------------------------------------------------------------------------
 
subroutine ufo_radiance_tlad_settraj(self, geovals, obss)
  implicit none
  class(ufo_radiance_tlad), intent(inout) :: self
  type(ufo_geovals),     intent(in)    :: geovals
  type(ioda_obsdb),      intent(in)    :: obss

  character(len=*), parameter :: myname_="ufo_radiance_tlad_settraj"
  character(max_string) :: err_msg

  ! Nothing here yet

  self%ltraj = .false. !.true.

end subroutine ufo_radiance_tlad_settraj

! ------------------------------------------------------------------------------

subroutine ufo_radiance_tlad_eqv_tl(self, geovals, hofx, obss)
  implicit none
  class(ufo_radiance_tlad), intent(in)     :: self
  type(ufo_geovals),     intent(in)     :: geovals
  type(obs_vector),      intent(inout)  :: hofx
  type(ioda_obsdb),      intent(in)     :: obss

  character(len=*), parameter :: myname_="ufo_radiance_tlad_eqv_tl"
  character(max_string) :: err_msg

  ! Nothing here yet

end subroutine ufo_radiance_tlad_eqv_tl

! ------------------------------------------------------------------------------

subroutine ufo_radiance_tlad_eqv_ad(self, geovals, hofx, obss)
  implicit none
  class(ufo_radiance_tlad), intent(in) :: self
  type(ufo_geovals),     intent(inout)    :: geovals
  type(obs_vector),      intent(in)       :: hofx
  type(ioda_obsdb),      intent(in)       :: obss

  character(len=*), parameter :: myname_="ufo_radiance_tlad_eqv_ad"
  character(max_string) :: err_msg

  ! Nothing here yet

end subroutine ufo_radiance_tlad_eqv_ad

! ------------------------------------------------------------------------------
 
end module ufo_radiance_tlad_mod
