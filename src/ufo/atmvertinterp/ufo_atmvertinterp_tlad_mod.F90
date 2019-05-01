! (C) Copyright 2017-2018 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

module ufo_atmvertinterp_tlad_mod

  use iso_c_binding
  use kinds
  use ufo_vars_mod
  use ufo_geovals_mod
  use ufo_geovals_mod_c,   only: ufo_geovals_registry
  use vert_interp_mod
  use obsspace_mod
  use missing_values_mod

  integer, parameter :: max_string=800

! ------------------------------------------------------------------------------

  type, public :: ufo_atmvertinterp_tlad
  private
     integer :: nvars
     character(len=max_string), public, allocatable :: varin(:)
     integer :: nval, nlocs
     real(kind_real), allocatable :: wf(:)
     integer, allocatable :: wi(:)
  contains
    procedure :: setup => atmvertinterp_tlad_setup_
    procedure :: cleanup => atmvertinterp_tlad_cleanup_
    procedure :: settraj => atmvertinterp_tlad_settraj_
    procedure :: simobs_tl => atmvertinterp_simobs_tl_
    procedure :: simobs_ad => atmvertinterp_simobs_ad_
    final :: destructor
  end type ufo_atmvertinterp_tlad

! ------------------------------------------------------------------------------
contains
! ------------------------------------------------------------------------------

subroutine atmvertinterp_tlad_setup_(self, c_conf)
  use config_mod
  implicit none
  class(ufo_atmvertinterp_tlad), intent(inout) :: self
  type(c_ptr), intent(in)    :: c_conf

  integer :: ii

  !> Size of variables
  self%nvars = size(config_get_string_vector(c_conf, max_string, "variables"))
  !> Allocate varin
  allocate(self%varin(self%nvars))
  !> Read variable list and store in varin
  self%varin = config_get_string_vector(c_conf, max_string, "variables")

end subroutine atmvertinterp_tlad_setup_

! ------------------------------------------------------------------------------

subroutine atmvertinterp_tlad_settraj_(self, geovals, obss)
  implicit none
  class(ufo_atmvertinterp_tlad), intent(inout) :: self
  type(ufo_geovals),         intent(in)    :: geovals
  type(c_ptr), value,        intent(in)    :: obss

  real(kind_real), allocatable :: obspressure(:)
  type(ufo_geoval), pointer :: presprofile
  integer :: iobs

  ! Make sure nothing already allocated
  call self%cleanup()

  ! Get pressure profiles from geovals
  call ufo_geovals_get_var(geovals, var_prsl, presprofile)
  self%nval = presprofile%nval

  ! Get the observation vertical coordinates
  self%nlocs = obsspace_get_nlocs(obss)
  allocate(obspressure(self%nlocs))
  call obsspace_get_db(obss, "MetaData", "air_pressure", obspressure)

  ! Allocate arrays for interpolation weights
  allocate(self%wi(self%nlocs))
  allocate(self%wf(self%nlocs))

  ! Calculate the interpolation weights
  do iobs = 1, self%nlocs
    call vert_interp_weights(presprofile%nval, log(obspressure(iobs)/10.), &
                             presprofile%vals(:,iobs), self%wi(iobs), self%wf(iobs))
  enddo

  ! Cleanup memory
  deallocate(obspressure)
end subroutine atmvertinterp_tlad_settraj_

! ------------------------------------------------------------------------------

subroutine atmvertinterp_simobs_tl_(self, geovals, obss, nvars, nlocs, hofx)
  implicit none
  class(ufo_atmvertinterp_tlad), intent(in) :: self
  type(ufo_geovals),         intent(in) :: geovals
  integer,                   intent(in) :: nvars, nlocs
  real(c_double),         intent(inout) :: hofx(nvars, nlocs)
  type(c_ptr), value,        intent(in) :: obss
  
  integer :: iobs, ivar
  type(ufo_geoval), pointer :: profile
  character(len=MAXVARLEN) :: geovar

  do ivar = 1, self%nvars
    ! Get the name of input variable in geovals
    geovar = self%varin(ivar)

    ! Get profile for this variable from geovals
    call ufo_geovals_get_var(geovals, geovar, profile)

    ! Interpolate from geovals to observational location into hofx
    do iobs = 1, self%nlocs
      call vert_interp_apply_tl(profile%nval, profile%vals(:,iobs), &
                                & hofx(ivar,iobs), self%wi(iobs), self%wf(iobs))
    enddo
  enddo
end subroutine atmvertinterp_simobs_tl_

! ------------------------------------------------------------------------------

subroutine atmvertinterp_simobs_ad_(self, geovals, obss, nvars, nlocs, hofx)
  implicit none
  class(ufo_atmvertinterp_tlad), intent(in) :: self
  type(ufo_geovals),         intent(inout) :: geovals
  integer,                   intent(in)    :: nvars, nlocs
  real(c_double),            intent(in)    :: hofx(nvars, nlocs)
  type(c_ptr), value,        intent(in)    :: obss
  
  integer :: iobs, ivar
  type(ufo_geoval), pointer :: profile
  character(len=MAXVARLEN) :: geovar
  real(c_double) :: missing

  missing = missing_value(missing)

  do ivar = 1, self%nvars
    ! Get the name of input variable in geovals
    geovar = self%varin(ivar)

    ! Get pointer to profile for this variable in geovals
    call ufo_geovals_get_var(geovals, geovar, profile)
      
    ! Allocate geovals profile if not yet allocated
    if (.not. allocated(profile%vals)) then
       profile%nlocs = self%nlocs
       profile%nval  = self%nval
       allocate(profile%vals(profile%nval, profile%nlocs))
       profile%vals(:,:) = 0.0_kind_real
    endif
    if (.not. geovals%linit ) geovals%linit=.true.

    ! Adjoint of interpolate, from hofx into geovals
    do iobs = 1, self%nlocs
      if (hofx(ivar,iobs) /= missing) then
        call vert_interp_apply_ad(profile%nval, profile%vals(:,iobs), &
                                & hofx(ivar,iobs), self%wi(iobs), self%wf(iobs))
      endif
    enddo
  enddo
end subroutine atmvertinterp_simobs_ad_

! ------------------------------------------------------------------------------

subroutine atmvertinterp_tlad_cleanup_(self)
  implicit none
  class(ufo_atmvertinterp_tlad), intent(inout) :: self
  self%nval = 0
  self%nlocs = 0
  if (allocated(self%wi)) deallocate(self%wi)
  if (allocated(self%wf)) deallocate(self%wf)
end subroutine atmvertinterp_tlad_cleanup_

! ------------------------------------------------------------------------------

subroutine  destructor(self)
  type(ufo_atmvertinterp_tlad), intent(inout)  :: self

  call self%cleanup()
  self%nvars = 0
  if (allocated(self%varin)) deallocate(self%varin)

end subroutine destructor

! ------------------------------------------------------------------------------

end module ufo_atmvertinterp_tlad_mod
