! (C) Copyright 2017-2018 UCAR
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 

module ufo_basis_mod

  use iso_c_binding
  use ufo_geovals_mod
  use ufo_geovals_mod_c,   only: ufo_geovals_registry
  use ioda_obs_vectors
  use ioda_obsdb_mod, only: ioda_obsdb
  use ioda_obsdb_mod_c, only: ioda_obsdb_registry

  type, abstract :: ufo_basis
    private
  contains
    procedure, non_overridable :: opr_eqv  => opr_eqv_
    procedure(eqv_), deferred  :: eqv
  end type ufo_basis

  abstract interface

  ! ------------------------------------------------------------------------------

    subroutine eqv_(self, geovals, hofx, obss)
      import ufo_basis, ufo_geovals, obs_vector, ioda_obsdb
      implicit none
      class(ufo_basis),         intent(in)    :: self
      type(ufo_geovals),        intent(in)    :: geovals
      type(obs_vector),         intent(inout) :: hofx
      type(ioda_obsdb), target, intent(in)    :: obss
    end subroutine
  
  ! ------------------------------------------------------------------------------

  end interface
  
contains

! ------------------------------------------------------------------------------
    
    subroutine opr_eqv_(self, c_key_geovals, c_key_obsspace, c_key_hofx)
      implicit none
    
      class(ufo_basis), intent(in)    :: self
      integer(c_int), intent(in) :: c_key_geovals
      integer(c_int), intent(in) :: c_key_hofx
      integer(c_int), intent(in) :: c_key_obsspace
    
      type(ufo_geovals), pointer :: geovals
      type(obs_vector),  pointer :: hofx
      type(ioda_obsdb),  pointer :: obss
    
      call ufo_geovals_registry%get(c_key_geovals,geovals)
      call ioda_obs_vect_registry%get(c_key_hofx,hofx)
      call ioda_obsdb_registry%get(c_key_obsspace,obss)
    
      call self%eqv(geovals, hofx, obss)
    
    end subroutine opr_eqv_
    
! ------------------------------------------------------------------------------

end module ufo_basis_mod
