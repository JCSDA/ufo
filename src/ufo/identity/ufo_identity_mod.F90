! (C) Copyright 2017-2018 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

! Fortran module for identity observation operator
!---------------------------------------------------------------------------------------------------
module ufo_identity_mod

 use iso_c_binding
 use kinds
 use ufo_vars_mod
 use ufo_geovals_mod

 use obsspace_mod

 integer, parameter :: max_string=800

! Fortran derived type for the observation type
!---------------------------------------------------------------------------------------------------
 type, public :: ufo_identity
 private
    integer, public :: nvars
    character(len=max_string), public, allocatable :: varin(:)
    character(len=max_string), public, allocatable :: varout(:)
 contains
   procedure :: setup  => identity_setup_
   procedure :: simobs => identity_simobs_
   final :: destructor
 end type ufo_identity

contains

! ------------------------------------------------------------------------------
subroutine identity_setup_(self, c_conf)
   use config_mod
   implicit none
   class(ufo_identity), intent(inout) :: self
   type(c_ptr),        intent(in)     :: c_conf

   integer :: ii

  !> Size of variables
  self%nvars = size(config_get_string_vector(c_conf, max_string, "variables"))

  !> Allocate varout: variables in the observation vector
  allocate(self%varin(self%nvars))

  !> Read variable list and store in varout
  self%varin = config_get_string_vector(c_conf, max_string, "variables")

  !> -----------------------------------------------------------------------------
  !> Allocate varin: variables we need from the model
  !> need additional slot to hold vertical coord.

  allocate(self%varout(self%nvars))

  !> Set vars_in based on vars_out
  do ii = 1, self%nvars
    self%varout(ii) = self%varin(ii)
  enddo

end subroutine identity_setup_


! ------------------------------------------------------------------------------
subroutine identity_simobs_(self, geovals, obss, nvars, nlocs, hofx)
  implicit none
  class(ufo_identity), intent(in)    :: self
  type(ufo_geovals),  intent(in)     :: geovals
  integer,            intent(in)     :: nvars, nlocs
  real(c_double),     intent(inout)  :: hofx(nvars, nlocs)
  type(c_ptr), value, intent(in)     :: obss

  integer :: iobs, ivar
  type(ufo_geoval), pointer :: point
  character(len=MAXVARLEN) :: geovar

  do ivar = 1, self%nvars
    !> Get the name of input variable in geovals
    geovar = self%varin(ivar)

    !> Get profile for this variable from geovals
    call ufo_geovals_get_var(geovals, geovar, point)

    !> Here we just apply a identity hofx
    do iobs = 1, nlocs
      hofx(ivar,iobs) = point%vals(1,iobs)
    enddo
  enddo

end subroutine identity_simobs_


! ------------------------------------------------------------------------------
subroutine destructor(self)
  type(ufo_identity), intent(inout) :: self
  if (allocated(self%varout)) deallocate(self%varout)
  if (allocated(self%varin)) deallocate(self%varin)
end subroutine destructor


end module ufo_identity_mod
