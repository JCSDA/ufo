! (C) Copyright 2017-2018 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

module ufo_conventional_profile_tlad_mod

  use iso_c_binding
  use kinds
  use ufo_vars_mod
  use ufo_geovals_mod
  use ufo_geovals_mod_c,   only: ufo_geovals_registry
  use vert_interp_mod
  use ufo_basis_tlad_mod, only: ufo_basis_tlad
  use obsspace_mod

  integer, parameter :: max_string=800

  type, extends(ufo_basis_tlad) :: ufo_conventional_profile_tlad
   private
     integer :: nval, nlocs
     real(kind_real), allocatable :: wf(:)
     integer, allocatable :: wi(:)
     integer, public :: nvars
     character(len=max_string), public, allocatable :: varin(:)
     character(len=max_string), public, allocatable :: varout(:)
  contains
    procedure :: delete => conventional_profile_tlad_delete_
    procedure :: settraj => conventional_profile_tlad_settraj_
    procedure :: simobs_tl => conventional_profile_simobs_tl_
    procedure :: simobs_ad => conventional_profile_simobs_ad_
    final :: destructor
  end type ufo_conventional_profile_tlad
contains

! ------------------------------------------------------------------------------

subroutine conventional_profile_tlad_settraj_(self, geovals, obss)
  implicit none
  class(ufo_conventional_profile_tlad), intent(inout) :: self
  type(ufo_geovals),         intent(in)    :: geovals
  type(c_ptr), value,        intent(in)    :: obss

  real(kind_real), allocatable :: obspressure(:)
  type(ufo_geoval), pointer :: presprofile
  integer :: iobs

  ! Make sure nothing already allocated
  call self%delete()

  ! Get pressure profiles from geovals
  call ufo_geovals_get_var(geovals, var_prsl, presprofile)

  self%nval = presprofile%nval
  self%nlocs = obsspace_get_nlocs(obss)

  ! Get the vertical coordinate and its dimension for this variable
  allocate(obspressure(self%nlocs))
  call obsspace_get_db(obss, "MetaData", "air_pressure", obspressure)

  ! Allocate arrays for vertical coordinate (pressure) and interpolation weights
  allocate(self%wi(self%nlocs))
  allocate(self%wf(self%nlocs))

  ! Calculate the interpolation weights
  do iobs = 1, self%nlocs
    call vert_interp_weights(presprofile%nval, log(obspressure(iobs)/10.), &
                             presprofile%vals(:,iobs), self%wi(iobs), self%wf(iobs))
  enddo

  self%ltraj = .true.
  ! Cleanup memory
  deallocate(obspressure)
end subroutine conventional_profile_tlad_settraj_

! ------------------------------------------------------------------------------

subroutine conventional_profile_simobs_tl_(self, geovals, hofx, obss)
  implicit none
  class(ufo_conventional_profile_tlad), intent(in) :: self
  type(ufo_geovals),         intent(in) :: geovals
  real(c_double),         intent(inout) :: hofx(:)
  type(c_ptr), value,        intent(in) :: obss
  
  character(len=*), parameter :: myname_="ufo_conventional_profile_simobs_tl"
  character(max_string) :: err_msg

  integer :: iobs, ivar

  type(ufo_geoval), pointer :: profile
  character(len=MAXVARLEN) :: geovar

  ! check if trajectory was set
  if (.not. self%ltraj) then
    write(err_msg,*) myname_, ' trajectory wasnt set!'
    call abor1_ftn(err_msg)
  endif

  do ivar = 1, self%nvars
    ! Get the name of input variable in geovals
    geovar = self%varout(ivar)
    if (trim(geovar) == "air_temperature") geovar = "virtual_temperature"

    ! Get profile for this variable in geovals
    call ufo_geovals_get_var(geovals, geovar, profile)

    ! Interpolate from geovals to observational location into hofx
    do iobs = 1, self%nlocs
      call vert_interp_apply_tl(profile%nval, profile%vals(:,iobs), &
                              & hofx(ivar+(iobs-1)*self%nvars), &
                              & self%wi(iobs), self%wf(iobs))
    enddo
  enddo

end subroutine conventional_profile_simobs_tl_

! ------------------------------------------------------------------------------

subroutine conventional_profile_simobs_ad_(self, geovals, hofx, obss)
  implicit none
  class(ufo_conventional_profile_tlad), intent(in)     :: self
  type(ufo_geovals),         intent(inout)  :: geovals
  real(c_double),            intent(in)     :: hofx(:)
  type(c_ptr), value,        intent(in)     :: obss
  
  character(len=*), parameter :: myname_="ufo_conventional_profile_simobs_ad"
  character(max_string) :: err_msg
  integer :: iobs, ivar
  real(c_double) :: missing_value

  type(ufo_geoval), pointer :: profile
  character(len=MAXVARLEN) :: geovar

  if (.not. self%ltraj) then
    write(err_msg,*) myname_, ' trajectory wasnt set!'
    call abor1_ftn(err_msg)
  endif

  missing_value = obspace_missing_value()

  do ivar = 1, self%nvars
    ! Get the name of input variable in geovals
    geovar = self%varout(ivar)
    if (trim(geovar) == "air_temperature") geovar = "virtual_temperature"

    ! Get profile for this variable in geovals
    call ufo_geovals_get_var(geovals, geovar, profile)
      
    ! Allocate geovals profile if not yet allocated
    if (.not. allocated(profile%vals)) then
       profile%nobs = self%nlocs
       profile%nval = self%nval
       allocate(profile%vals(profile%nval, profile%nobs))
       profile%vals(:,:) = 0.0_kind_real
    endif
    if (.not. geovals%linit ) geovals%linit=.true.

    ! Interpolate from geovals to observational location into hofx
    do iobs = 1, self%nlocs
      if (hofx(ivar+(iobs-1)*self%nvars) /= missing_value) then
        call vert_interp_apply_ad(profile%nval, profile%vals(:,iobs), &
                                & hofx(ivar+(iobs-1)*self%nvars), &
                                & self%wi(iobs), self%wf(iobs))
      endif
    enddo
  enddo

end subroutine conventional_profile_simobs_ad_

! ------------------------------------------------------------------------------

subroutine conventional_profile_tlad_delete_(self)
  implicit none
  class(ufo_conventional_profile_tlad), intent(inout) :: self
  self%nval = 0
  self%ltraj = .false.
! Only deleting trajectory
  if (allocated(self%wi)) deallocate(self%wi)
  if (allocated(self%wf)) deallocate(self%wf)
end subroutine conventional_profile_tlad_delete_

! ------------------------------------------------------------------------------

subroutine  destructor(self)
  type(ufo_conventional_profile_tlad), intent(inout)  :: self
  self%nval = 0
  self%ltraj = .false.
  if (allocated(self%varin)) deallocate(self%varin)
  if (allocated(self%varout)) deallocate(self%varout)
  if (allocated(self%wi)) deallocate(self%wi)
  if (allocated(self%wf)) deallocate(self%wf)
end subroutine destructor

! ------------------------------------------------------------------------------

end module ufo_conventional_profile_tlad_mod
