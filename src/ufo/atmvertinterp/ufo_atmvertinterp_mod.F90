! (C) Copyright 2017-2018 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

module ufo_atmvertinterp_mod

  use iso_c_binding
  use kinds
  use ufo_vars_mod
  use ufo_geovals_mod
  use ufo_geovals_mod_c,   only: ufo_geovals_registry
  use vert_interp_mod
  use obsspace_mod

  integer, parameter :: max_string=800

! ------------------------------------------------------------------------------

  type, public :: ufo_atmvertinterp
   private
     integer :: nvars  ! number of variables to be interpolated
     character(len=max_string), public, allocatable :: varin(:)    ! size nvars+1 (+1 for log pressure)
     character(len=max_string), public, allocatable :: varout(:)   ! size nvars
   contains
     procedure :: setup  => atmvertinterp_setup_
     procedure :: simobs => atmvertinterp_simobs_
     final :: destructor
  end type ufo_atmvertinterp

! ------------------------------------------------------------------------------
contains
! ------------------------------------------------------------------------------

subroutine atmvertinterp_setup_(self, c_conf)
  use config_mod
  implicit none
  class(ufo_atmvertinterp), intent(inout) :: self
  type(c_ptr), intent(in)    :: c_conf

  integer :: ii

  !> Size of variables
  self%nvars = size(config_get_string_vector(c_conf, max_string, "variables"))
  !> Allocate varout: variables in the observation vector
  allocate(self%varout(self%nvars))
  !> Read variable list and store in varout
  self%varout = config_get_string_vector(c_conf, max_string, "variables")
  !> Allocate varin: variables we need from the model
  !  need additional slot to hold vertical coord.
  allocate(self%varin(self%nvars+1))
  !> Set vars_in based on vars_out
  do ii = 1, self%nvars
    self%varin(ii) = self%varout(ii)
  enddo
  !> Put log pressure to the varin (vars from the model) list
  self%varin(self%nvars+1) = "atmosphere_ln_pressure_coordinate"

end subroutine atmvertinterp_setup_

! ------------------------------------------------------------------------------

subroutine atmvertinterp_simobs_(self, geovals, obss, nvars, nlocs, hofx)

  implicit none
  class(ufo_atmvertinterp), intent(in)        :: self
  integer, intent(in)                         :: nvars, nlocs
  type(ufo_geovals), intent(in)               :: geovals
  real(c_double),  intent(inout)              :: hofx(nvars, nlocs)
  type(c_ptr), value, intent(in)              :: obss

  integer :: iobs, ivar
  real(kind_real), dimension(:), allocatable :: obspressure
  type(ufo_geoval), pointer :: presprofile, profile
  real(kind_real), allocatable :: wf(:)
  integer, allocatable :: wi(:)
  character(len=MAXVARLEN) :: geovar

  ! Get pressure profiles from geovals
  call ufo_geovals_get_var(geovals, var_prsl, presprofile)

  ! Get the observation vertical coordinates
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

  do ivar = 1, self%nvars
    ! Get the name of input variable in geovals
    geovar = self%varin(ivar)

    ! Get profile for this variable from geovals
    call ufo_geovals_get_var(geovals, geovar, profile)

    ! Interpolate from geovals to observational location into hofx
    do iobs = 1, nlocs
      call vert_interp_apply(profile%nval, profile%vals(:,iobs), &
                             & hofx(ivar,iobs), wi(iobs), wf(iobs))
    enddo
  enddo
  ! Cleanup memory
  deallocate(obspressure)
  deallocate(wi)
  deallocate(wf)

end subroutine atmvertinterp_simobs_

! ------------------------------------------------------------------------------

subroutine destructor(self)
  type(ufo_atmvertinterp), intent(inout) :: self
  if (allocated(self%varout)) deallocate(self%varout)
  if (allocated(self%varin)) deallocate(self%varin)
end subroutine destructor

! ------------------------------------------------------------------------------

end module ufo_atmvertinterp_mod
