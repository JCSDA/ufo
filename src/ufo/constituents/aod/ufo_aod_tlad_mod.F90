! (C) Copyright 2017-2018 UCAR
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 

!> Fortran module to handle aod observations

MODULE ufo_aod_tlad_mod
  use iso_c_binding
  use ufo_vars_mod
  use ufo_locs_mod
  use ufo_geovals_mod
  use kinds  
  use ufo_basis_tlad_mod, only: ufo_basis_tlad
  implicit none
  
  public :: ufo_aod_tlad
    
  !> Fortran derived type for aod trajectory
  type, extends(ufo_basis_tlad) :: ufo_aod_tlad
  contains
    procedure :: delete  => ufo_aod_tlad_delete
    procedure :: settraj => ufo_aod_tlad_settraj 
    procedure :: simobs_tl  => ufo_aod_simobs_tl
    procedure :: simobs_ad  => ufo_aod_simobs_ad
  end type ufo_aod_tlad

contains

! ------------------------------------------------------------------------------

  subroutine ufo_aod_tlad_delete(self)
    implicit none
    class(ufo_aod_tlad), intent(inout)  :: self

    ! Nothing here yet

  end subroutine ufo_aod_tlad_delete

! ------------------------------------------------------------------------------

  subroutine ufo_aod_tlad_settraj(self, geovals, obss)
    implicit none
    class(ufo_aod_tlad), intent(inout) :: self
    type(ufo_geovals),   intent(in)    :: geovals
    type(c_ptr), value,  intent(in)    :: obss

    character(len=*), parameter :: myname_="ufo_aod_tlad_settraj"

    ! Nothing here yet

  end subroutine ufo_aod_tlad_settraj

! ------------------------------------------------------------------------------

  subroutine ufo_aod_simobs_tl(self, geovals, hofx, obss)
    implicit none
    class(ufo_aod_tlad), intent(in)     :: self
    type(ufo_geovals),   intent(in)     :: geovals
    real(kind_real),     intent(inout)  :: hofx(:)
    type(c_ptr), value,  intent(in)     :: obss

    character(len=*), parameter :: myname_="ufo_aod_simobs_tl"

    ! Nothing here yet

  end subroutine ufo_aod_simobs_tl

! ------------------------------------------------------------------------------

  subroutine ufo_aod_simobs_ad(self, geovals, hofx, obss)
    implicit none
    class(ufo_aod_tlad), intent(in)    :: self
    type(ufo_geovals),   intent(inout) :: geovals
    real(kind_real),     intent(in)    :: hofx(:)
    type(c_ptr), value,  intent(in)    :: obss

    character(len=*), parameter :: myname_="ufo_aod_simobs_tl"

    ! Nothing here yet

  end subroutine ufo_aod_simobs_ad

! ------------------------------------------------------------------------------

END MODULE ufo_aod_tlad_mod
