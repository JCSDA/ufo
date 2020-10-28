! (C) Copyright 2017-2018 UCAR
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 

module ufo_basis_mod

  use iso_c_binding
  use ufo_geovals_mod
  use ufo_geovals_mod_c,   only: ufo_geovals_registry

  type, abstract :: ufo_basis
    private
  contains
    procedure, non_overridable :: opr_simobs  => opr_simobs_
    procedure(simobs_), deferred  :: simobs
  end type ufo_basis

  abstract interface

  ! ------------------------------------------------------------------------------

    subroutine simobs_(self, geovals, hofx, obss)
      use iso_c_binding
      import ufo_basis, ufo_geovals
      implicit none
      class(ufo_basis),         intent(in)    :: self
      type(ufo_geovals),        intent(in)    :: geovals
      real(c_double),           intent(inout) :: hofx(:)
      type(c_ptr), value,       intent(in)    :: obss
    end subroutine simobs_
  
  ! ------------------------------------------------------------------------------

  end interface
  
contains

! ------------------------------------------------------------------------------
    
  subroutine opr_simobs_(self, c_key_geovals, c_obsspace, c_hofx)
    implicit none
  
    class(ufo_basis), intent(in)    :: self
    integer(c_int), intent(in) :: c_key_geovals
    type(c_ptr), value, intent(in) :: c_obsspace
    real(c_double), intent(inout) :: c_hofx(:)
  
    type(ufo_geovals), pointer :: geovals
  
    call ufo_geovals_registry%get(c_key_geovals,geovals)
  
    call self%simobs(geovals, c_hofx, c_obsspace)
  
  end subroutine opr_simobs_
    
! ------------------------------------------------------------------------------

end module ufo_basis_mod
