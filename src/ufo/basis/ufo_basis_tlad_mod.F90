! (C) Copyright 2017 UCAR
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 

module ufo_basis_tlad_mod

  use iso_c_binding
  use ufo_geovals_mod
  use ufo_geovals_mod_c,   only: ufo_geovals_registry
  use ioda_obs_vectors
  use ioda_obsdb_mod, only: ioda_obsdb
  use ioda_obsdb_mod_c, only: ioda_obsdb_registry

  type, abstract :: ufo_basis_tlad
    private
      logical, public :: ltraj = .false. !< trajectory set?
  contains
    procedure, non_overridable    :: opr_delete  => opr_delete_
    procedure, non_overridable    :: opr_settraj => opr_settraj_
    procedure, non_overridable    :: opr_eqv_tl  => opr_eqv_tl_
    procedure, non_overridable    :: opr_eqv_ad  => opr_eqv_ad_
    procedure(delete_),  deferred :: delete
    procedure(settraj_), deferred :: settraj
    procedure(eqv_tl_),  deferred :: eqv_tl
    procedure(eqv_ad_),  deferred :: eqv_ad
  end type ufo_basis_tlad

  abstract interface
  
  ! ------------------------------------------------------------------------------

    subroutine delete_(self)
      import ufo_basis_tlad
      implicit none
      class(ufo_basis_tlad), intent(inout) :: self
    end subroutine

! ------------------------------------------------------------------------------

    subroutine settraj_(self, geovals, obss)
      import ufo_basis_tlad, ufo_geovals, ioda_obsdb
      implicit none
      class(ufo_basis_tlad), intent(inout) :: self
      type(ufo_geovals),    intent(in) :: geovals
      type(ioda_obsdb),     intent(in) :: obss
    end subroutine

! ------------------------------------------------------------------------------

    subroutine eqv_tl_(self, geovals, hofx, obss)
      import ufo_basis_tlad, ufo_geovals, obs_vector, ioda_obsdb
      implicit none
      class(ufo_basis_tlad), intent(in) :: self
      type(ufo_geovals),     intent(in)    :: geovals
      type(obs_vector),      intent(inout) :: hofx
      type(ioda_obsdb),      intent(in)    :: obss
    end subroutine

! ------------------------------------------------------------------------------

    subroutine eqv_ad_(self, geovals, hofx, obss)
      import ufo_basis_tlad, ufo_geovals, obs_vector, ioda_obsdb
      implicit none
      class(ufo_basis_tlad), intent(in)    :: self
      type(ufo_geovals),     intent(inout)    :: geovals
      type(obs_vector),      intent(in)    :: hofx
      type(ioda_obsdb),      intent(in)    :: obss
    end subroutine

  ! ------------------------------------------------------------------------------

  end interface

contains
! ------------------------------------------------------------------------------

    subroutine opr_delete_(self)
      implicit none
      class(ufo_basis_tlad), intent(inout)  :: self
      
      call self%delete()
    
    end subroutine opr_delete_
    
! ------------------------------------------------------------------------------
    
    subroutine opr_settraj_(self, c_key_geovals, c_key_obsspace)
      implicit none
    
      class(ufo_basis_tlad), intent(inout) :: self
      integer(c_int), intent(in) :: c_key_geovals
      integer(c_int), intent(in) :: c_key_obsspace
      
      type(ufo_geovals),    pointer :: geovals
      type(ioda_obsdb),     pointer :: obss
    
      call ufo_geovals_registry%get(c_key_geovals,geovals)
      call ioda_obsdb_registry%get(c_key_obsspace,obss)
      
      call self%settraj(geovals, obss)
      
    end subroutine opr_settraj_
    
! ------------------------------------------------------------------------------
    
    subroutine opr_eqv_tl_(self, c_key_geovals, c_key_obsspace, c_key_hofx)
      implicit none
    
      class(ufo_basis_tlad), intent(in) :: self
      integer(c_int), intent(in) :: c_key_geovals
      integer(c_int), intent(in) :: c_key_hofx
      integer(c_int), intent(in) :: c_key_obsspace
      
      type(ufo_geovals),    pointer :: geovals
      type(obs_vector),     pointer :: hofx
      type(ioda_obsdb),     pointer :: obss
      
      call ufo_geovals_registry%get(c_key_geovals,geovals)
      call ioda_obs_vect_registry%get(c_key_hofx,hofx)
      call ioda_obsdb_registry%get(c_key_obsspace,obss)
      
      call self%eqv_tl(geovals, hofx, obss)
      
    end subroutine opr_eqv_tl_
    
! ------------------------------------------------------------------------------
    
    subroutine opr_eqv_ad_(self, c_key_geovals, c_key_obsspace, c_key_hofx)
      implicit none
    
      class(ufo_basis_tlad), intent(in) :: self
      integer(c_int), intent(in) :: c_key_geovals
      integer(c_int), intent(in) :: c_key_hofx
      integer(c_int), intent(in) :: c_key_obsspace
      
      type(ufo_geovals),    pointer :: geovals
      type(obs_vector),     pointer :: hofx
      type(ioda_obsdb),     pointer :: obss
    
      call ufo_geovals_registry%get(c_key_geovals,geovals)
      call ioda_obs_vect_registry%get(c_key_hofx,hofx)
      call ioda_obsdb_registry%get(c_key_obsspace,obss)
      
      call self%eqv_ad(geovals, hofx, obss)
      
    end subroutine opr_eqv_ad_

! ------------------------------------------------------------------------------

end module ufo_basis_tlad_mod
