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
   type(oops_variables), public :: obsvars ! Variables to be simulated
   integer, allocatable, public :: obsvarindices(:) ! Indices of obsvars in the list of all
                                                    ! simulated variables in the ObsSpace
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

  !> copy simulated variables to variables requested from the model
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
subroutine identity_simobs_tl_(self, geovals, obss, nvars, nlocs, hofx)
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
  type(c_ptr), value, intent(in)     :: obss
  integer,             intent(in)    :: nvars, nlocs
  real(c_double),     intent(inout)  :: hofx(nvars, nlocs)

  integer :: iobs, iobsvar, ivar
  type(ufo_geoval), pointer :: point
  character(len=MAXVARLEN) :: geovar

  do iobsvar = 1, size(self%obsvarindices)
    ! Get the index of the row of hofx to fill
    ivar = self%obsvarindices(iobsvar)

    !> Get the name of input variable in geovals
    geovar = self%geovars%variable(iobsvar)

    !> Get profile for this variable from geovals
    call ufo_geovals_get_var(geovals, geovar, point)

    !> Here we just apply a identity hofx
    do iobs = 1, nlocs
      hofx(ivar, iobs) = point%vals(1, iobs)
    enddo
  enddo

end subroutine identity_simobs_tl_

! ------------------------------------------------------------------------------
subroutine identity_simobs_ad_(self, geovals, obss, nvars, nlocs, hofx)
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
  type(c_ptr), value,      intent(in)    :: obss
  integer,                 intent(in)    :: nvars, nlocs
  real(c_double),          intent(in)    :: hofx(nvars, nlocs)

  integer :: iobs, iobsvar, ivar
  type(ufo_geoval), pointer :: point
  character(len=MAXVARLEN)  :: geovar
  real(c_double) :: missing

  !> Set missing value
  missing = missing_value(missing)

  if (.not. geovals%linit ) geovals%linit=.true.

  do iobsvar = 1, size(self%obsvarindices)
    ! Get the index of the row of hofx to fill
    ivar = self%obsvarindices(iobsvar)

    !> Get the name of input variable in geovals
    geovar = self%geovars%variable(iobsvar)

    !> Get profile for this variable from geovals
    call ufo_geovals_get_var(geovals, geovar, point)

    ! backward obs operator
    do iobs = 1, nlocs
      if (hofx(ivar, iobs) /= missing) then
        point%vals(1, iobs) = hofx(ivar, iobs)
      endif
    enddo
  enddo

end subroutine identity_simobs_ad_

end module ufo_identity_tlad_mod
