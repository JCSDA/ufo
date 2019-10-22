! (C) Copyright 2017-2018 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

!> Fortran module for identity tl/ad observation operator

module ufo_identity_tlad_mod

 use oops_variables_mod
 use ufo_vars_mod
 implicit none

 ! ------------------------------------------------------------------------------

 !> Fortran derived type for the tl/ad observation operator
 type, public :: ufo_identity_tlad
 private
   type(oops_variables), public :: obsvars
   type(oops_variables), public :: geovars
 contains
   procedure :: setup  => identity_tlad_setup_
   procedure :: settraj => identity_tlad_settraj_
   procedure :: simobs_tl  => identity_simobs_tl_
   procedure :: simobs_ad  => identity_simobs_ad_
 end type ufo_identity_tlad

contains

! ------------------------------------------------------------------------------
subroutine identity_tlad_setup_(self)
  implicit none
  class(ufo_identity_tlad), intent(inout) :: self

  integer :: ivar, nvars

  !> copy observed variables to variables requested from the model
  nvars = self%obsvars%nvars()
  do ivar = 1, nvars
    call self%geovars%push_back(self%obsvars%variable(ivar))
  enddo

end subroutine identity_tlad_setup_

!------------------------------------------------------------------------------
subroutine identity_tlad_settraj_(self, geovals, obss)
  use iso_c_binding
  use ufo_geovals_mod, only: ufo_geovals
  implicit none
  class(ufo_identity_tlad), intent(inout) :: self
  type(ufo_geovals),       intent(in)    :: geovals
  type(c_ptr), value,      intent(in)    :: obss

! since observation operator is linear, don't care about trajectory itself

end subroutine identity_tlad_settraj_

! ------------------------------------------------------------------------------
subroutine identity_simobs_tl_(self, geovals, hofx, obss)
  use iso_c_binding
  use kinds
  use ufo_geovals_mod, only: &
   ufo_geovals,             &
   ufo_geoval,              &
   ufo_geovals_get_var
  use obsspace_mod
  implicit none
  class(ufo_identity_tlad), intent(in)    :: self
  type(ufo_geovals),  intent(in)     :: geovals
  real(c_double),     intent(inout)  :: hofx(:)
  type(c_ptr), value, intent(in)     :: obss

  integer :: iobs, ivar, nlocs, nvars
  type(ufo_geoval), pointer :: point
  character(len=MAXVARLEN) :: geovar

  !> Get the observation vertical coordinates
  nlocs = obsspace_get_nlocs(obss)
  nvars = self%obsvars%nvars()
  do ivar = 1, nvars
    !> Get the name of input variable in geovals
    geovar = self%geovars%variable(ivar)

    !> Get profile for this variable from geovals
    call ufo_geovals_get_var(geovals, geovar, point)

    !> Here we just apply a identity hofx
    do iobs = 1, nlocs
      hofx(ivar + (iobs-1)*nvars) = point%vals(1,iobs)
    enddo
  enddo

end subroutine identity_simobs_tl_

! ------------------------------------------------------------------------------
subroutine identity_simobs_ad_(self, geovals, hofx, obss)
  use iso_c_binding
  use kinds
  use ufo_geovals_mod, only: &
   ufo_geovals,             &
   ufo_geoval,              &
   ufo_geovals_get_var
  use obsspace_mod
  use missing_values_mod
  implicit none
  class(ufo_identity_tlad), intent(in)   :: self
  type(ufo_geovals),       intent(inout) :: geovals
  real(c_double),          intent(in)    :: hofx(:)
  type(c_ptr), value,      intent(in)    :: obss

  integer :: iobs, ivar, nlocs, nvars
  type(ufo_geoval), pointer :: point
  character(len=MAXVARLEN)  :: geovar
  real(c_double) :: missing

  !> Set missing value
  missing = missing_value(missing)

  if (.not. geovals%linit ) geovals%linit=.true.

  nlocs = obsspace_get_nlocs(obss)
  nvars = self%obsvars%nvars()
  do ivar = 1, nvars
    !> Get the name of input variable in geovals
    geovar = self%geovars%variable(ivar)

    !> Get profile for this variable from geovals
    call ufo_geovals_get_var(geovals, geovar, point)

    if (.not.(allocated(point%vals))) then
      point%nval=1
      allocate(point%vals(1,size(hofx,1)))
      point%vals = 0.0
    end if

    ! backward obs operator
    do iobs = 1, nlocs
      if (hofx(ivar + (iobs-1)*nvars) /= missing) then
        point%vals(1,iobs) = hofx(ivar + (iobs-1)*nvars)
      endif
    enddo
  enddo

end subroutine identity_simobs_ad_

end module ufo_identity_tlad_mod
