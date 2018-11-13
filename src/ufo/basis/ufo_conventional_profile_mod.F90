! (C) Copyright 2017-2018 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

module ufo_conventional_profile_mod

  use iso_c_binding
  use kinds
  use ufo_vars_mod
  use ufo_geovals_mod
  use ufo_geovals_mod_c,   only: ufo_geovals_registry
  use vert_interp_mod
  use ufo_basis_mod, only: ufo_basis
  use obsspace_mod

  integer, parameter :: max_string=800

! ------------------------------------------------------------------------------

  type, extends(ufo_basis) :: ufo_conventional_profile
   private
     integer, public :: nvars
     character(len=max_string), public, allocatable :: varin(:)
     character(len=max_string), public, allocatable :: varout(:)
  contains
    procedure :: simobs    => conventional_profile_simobs_
    final :: destructor
  end type ufo_conventional_profile

! ------------------------------------------------------------------------------
contains
! ------------------------------------------------------------------------------

subroutine conventional_profile_simobs_(self, geovals, hofx, obss)

  implicit none
  class(ufo_conventional_profile), intent(in) :: self
  type(ufo_geovals), intent(in)               :: geovals
  real(c_double),  intent(inout)              :: hofx(:)
  type(c_ptr), value, intent(in)              :: obss

  integer :: iobs, ivar, nlocs
  real(kind_real), dimension(:), allocatable :: obspressure
  type(ufo_geoval), pointer :: presprofile, profile
  real(kind_real), allocatable :: wf(:)
  integer, allocatable :: wi(:)
  character(len=MAXVARLEN) :: geovar

  ! Get pressure profiles from geovals
  call ufo_geovals_get_var(geovals, var_prsl, presprofile)

  ! Get the observation vertical coordinates
  nlocs = obsspace_get_nlocs(obss)
  allocate(obspressure(nlocs))
  call obsspace_get_db(obss, "MetaData", "air_pressure", obspressure)

  ! Allocate arrays for interpolation weights
  allocate(wi(nlocs))
  allocate(wf(nlocs))

  ! Calculate the interpolation weights
  do iobs = 1, nlocs
    call vert_interp_weights(presprofile%nval, log(obspressure(iobs)/10.), &
                             presprofile%vals(:,iobs), wi(iobs), wf(iobs))
  enddo

  ivar = 1

  ! Get the name of input variable in geovals
  geovar = self%varout(ivar)
  if (trim(geovar) == "air_temperature") geovar = "virtual_temperature"

  ! Get profile for this variable from geovals
  call ufo_geovals_get_var(geovals, geovar, profile)

  ! Interpolate from geovals to observational location into hofx
  do iobs = 1, nlocs
    call vert_interp_apply(profile%nval, profile%vals(:,iobs), &
                             & hofx(ivar+(iobs-1)*self%nvars), &
                             & wi(iobs), wf(iobs))
  enddo

  ! Cleanup memory
  deallocate(obspressure)
  deallocate(wi)
  deallocate(wf)

end subroutine conventional_profile_simobs_

! ------------------------------------------------------------------------------

subroutine  destructor(self)
  type(ufo_conventional_profile), intent(inout) :: self
  if (allocated(self%varout)) deallocate(self%varout)
  if (allocated(self%varin)) deallocate(self%varin)
end subroutine destructor

! ------------------------------------------------------------------------------

end module ufo_conventional_profile_mod
