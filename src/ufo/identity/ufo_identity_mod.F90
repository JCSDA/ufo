! (C) Copyright 2017-2018 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

! Fortran module for identity observation operator
!---------------------------------------------------------------------------------------------------
module ufo_identity_mod

 use oops_variables_mod
 use ufo_vars_mod

! Fortran derived type for the observation type
!---------------------------------------------------------------------------------------------------
 type, public :: ufo_identity
   type(oops_variables), public :: obsvars
   type(oops_variables), public :: geovars
 contains
   procedure :: setup  => identity_setup_
   procedure :: simobs => identity_simobs_
 end type ufo_identity

contains

! ------------------------------------------------------------------------------
subroutine identity_setup_(self)
implicit none
class(ufo_identity), intent(inout) :: self

integer :: nvars, ivar

! set variables requested from the model (same as simulated):
nvars = self%obsvars%nvars()
do ivar = 1, nvars
  call self%geovars%push_back(self%obsvars%variable(ivar))
enddo

end subroutine identity_setup_


! ------------------------------------------------------------------------------
subroutine identity_simobs_(self, geovals, obss, nvars, nlocs, hofx)
  use ufo_geovals_mod
  use obsspace_mod
  use iso_c_binding
  implicit none
  class(ufo_identity), intent(in)    :: self
  type(ufo_geovals),  intent(in)     :: geovals
  integer,            intent(in)     :: nvars, nlocs
  real(c_double),     intent(inout)  :: hofx(nvars, nlocs)
  type(c_ptr), value, intent(in)     :: obss

  integer :: iobs, ivar
  type(ufo_geoval), pointer :: point
  character(len=MAXVARLEN) :: geovar

  do ivar = 1, nvars
    !> Get the name of input variable in geovals
    geovar = self%geovars%variable(ivar)

    !> Get profile for this variable from geovals
    call ufo_geovals_get_var(geovals, geovar, point)

    !> Here we just apply a identity hofx
    do iobs = 1, nlocs
      hofx(ivar,iobs) = point%vals(1,iobs)
    enddo
  enddo

end subroutine identity_simobs_

end module ufo_identity_mod
