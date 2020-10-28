! (C) Copyright 2017-2018 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

module ufo_basis_tlad_mod

  use iso_c_binding
  use ufo_geovals_mod
  use ufo_geovals_mod_c,   only: ufo_geovals_registry

  type, abstract :: ufo_basis_tlad
    private
      logical, public :: ltraj = .false. !< trajectory set?
  contains
    procedure, non_overridable    :: opr_delete  => opr_delete_
    procedure, non_overridable    :: opr_settraj => opr_settraj_
    procedure, non_overridable    :: opr_simobs_tl  => opr_simobs_tl_
    procedure, non_overridable    :: opr_simobs_ad  => opr_simobs_ad_
    procedure(delete_),  deferred :: delete
    procedure(settraj_), deferred :: settraj
    procedure(simobs_tl_),  deferred :: simobs_tl
    procedure(simobs_ad_),  deferred :: simobs_ad
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
      import ufo_basis_tlad, ufo_geovals, c_ptr
      implicit none
      class(ufo_basis_tlad), intent(inout) :: self
      type(ufo_geovals),    intent(in) :: geovals
      type(c_ptr), value,   intent(in) :: obss
    end subroutine

! ------------------------------------------------------------------------------

    subroutine simobs_tl_(self, geovals, hofx, obss)
      use iso_c_binding
      import ufo_basis_tlad, ufo_geovals
      implicit none
      class(ufo_basis_tlad), intent(in) :: self
      type(ufo_geovals),     intent(in)    :: geovals
      real(c_double),        intent(inout) :: hofx(:)
      type(c_ptr), value,    intent(in)    :: obss
    end subroutine

! ------------------------------------------------------------------------------

    subroutine simobs_ad_(self, geovals, hofx, obss)
      use iso_c_binding
      import ufo_basis_tlad, ufo_geovals
      implicit none
      class(ufo_basis_tlad), intent(in)    :: self
      type(ufo_geovals),     intent(inout) :: geovals
      real(c_double),        intent(in)    :: hofx(:)
      type(c_ptr), value,    intent(in)    :: obss
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

    subroutine opr_settraj_(self, c_key_geovals, c_obsspace)
      implicit none

      class(ufo_basis_tlad), intent(inout) :: self
      integer(c_int), intent(in) :: c_key_geovals
      type(c_ptr), value, intent(in) :: c_obsspace

      type(ufo_geovals),    pointer :: geovals

      call ufo_geovals_registry%get(c_key_geovals,geovals)

      call self%settraj(geovals, c_obsspace)
    end subroutine opr_settraj_

! ------------------------------------------------------------------------------

    subroutine opr_simobs_tl_(self, c_key_geovals, c_obsspace, c_hofx)
      implicit none

      class(ufo_basis_tlad), intent(in) :: self
      integer(c_int), intent(in) :: c_key_geovals
      real(c_double), intent(inout) :: c_hofx(:)
      type(c_ptr), value, intent(in) :: c_obsspace

      type(ufo_geovals),    pointer :: geovals

      call ufo_geovals_registry%get(c_key_geovals,geovals)

      call self%simobs_tl(geovals, c_hofx, c_obsspace)
    end subroutine opr_simobs_tl_

! ------------------------------------------------------------------------------

    subroutine opr_simobs_ad_(self, c_key_geovals, c_obsspace, c_hofx)
      implicit none

      class(ufo_basis_tlad), intent(in) :: self
      integer(c_int), intent(in) :: c_key_geovals
      real(c_double), intent(in) :: c_hofx(:)
      type(c_ptr), value, intent(in) :: c_obsspace

      type(ufo_geovals),    pointer :: geovals

      call ufo_geovals_registry%get(c_key_geovals,geovals)

      call self%simobs_ad(geovals, c_hofx, c_obsspace)
    end subroutine opr_simobs_ad_

! ------------------------------------------------------------------------------

end module ufo_basis_tlad_mod
