! (C) Copyright 2017-2019 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

module ufo_atmvertinterp_tlad_mod

  use oops_variables_mod
  use ufo_vars_mod
  use ufo_geovals_mod
  use vert_interp_mod
  use missing_values_mod


! ------------------------------------------------------------------------------

  type, public :: ufo_atmvertinterp_tlad
  private
    type(oops_variables), public :: obsvars
    type(oops_variables), public :: geovars
    integer :: nval, nlocs
    real(kind_real), allocatable :: wf(:)
    integer, allocatable :: wi(:)
    character(len=MAXVARLEN), public :: v_coord ! GeoVaL to use to interpolate in vertical
    logical, public :: use_ln ! if T, use ln(v_coord) not v_coord
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

subroutine atmvertinterp_tlad_setup_(self, grid_conf)
  use fckit_configuration_module, only: fckit_configuration
  implicit none
  class(ufo_atmvertinterp_tlad), intent(inout) :: self
  type(fckit_configuration), intent(in) :: grid_conf

  character(kind=c_char,len=:), allocatable :: coord_name
  integer :: ivar, nvars

  !> Fill in variables requested from the model
  nvars = self%obsvars%nvars()
  do ivar = 1, nvars
    call self%geovars%push_back(self%obsvars%variable(ivar))
  enddo
  !> grab what vertical coordinate/variable to use from the config
  self%use_ln = .false.

  if( grid_conf%has("VertCoord") ) then
      call grid_conf%get_or_die("VertCoord",coord_name)
      self%v_coord = coord_name
      if( trim(self%v_coord) .eq. var_prs ) self%use_ln = .true.
  else  ! default
      self%v_coord = var_prs
      self%use_ln  = .true.
  endif

end subroutine atmvertinterp_tlad_setup_

! ------------------------------------------------------------------------------

subroutine atmvertinterp_tlad_settraj_(self, geovals, obss)
  use obsspace_mod
  implicit none
  class(ufo_atmvertinterp_tlad), intent(inout) :: self
  type(ufo_geovals),         intent(in)    :: geovals
  type(c_ptr), value,        intent(in)    :: obss

  real(kind_real), allocatable :: obsvcoord(:)
  type(ufo_geoval), pointer :: vcoordprofile
  integer :: iobs
  real(kind_real), allocatable :: tmp(:)
  real(kind_real) :: tmp2

  ! Make sure nothing already allocated
  call self%cleanup()

  ! Get pressure profiles from geovals
  call ufo_geovals_get_var(geovals, self%v_coord, vcoordprofile)
  self%nval = vcoordprofile%nval

  ! Get the observation vertical coordinates
  self%nlocs = obsspace_get_nlocs(obss)
  allocate(obsvcoord(self%nlocs))
  call obsspace_get_db(obss, "MetaData", self%v_coord, obsvcoord)

  ! Allocate arrays for interpolation weights
  allocate(self%wi(self%nlocs))
  allocate(self%wf(self%nlocs))

  ! Calculate the interpolation weights
  allocate(tmp(vcoordprofile%nval))
  do iobs = 1, self%nlocs
    if (self%use_ln) then
      tmp = log(vcoordprofile%vals(:,iobs))
      tmp2 = log(obsvcoord(iobs))
    else
      tmp = vcoordprofile%vals(:,iobs)
      tmp2 = obsvcoord(iobs)
    end if
    call vert_interp_weights(vcoordprofile%nval, tmp2, tmp, self%wi(iobs), self%wf(iobs))
  enddo

  ! Cleanup memory
  deallocate(obsvcoord)
  deallocate(tmp)

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

  do ivar = 1, nvars
    ! Get the name of input variable in geovals
    geovar = self%geovars%variable(ivar)

    ! Get profile for this variable from geovals
    call ufo_geovals_get_var(geovals, geovar, profile)

    ! Interpolate from geovals to observational location into hofx
    do iobs = 1, nlocs
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

  do ivar = 1, nvars
    ! Get the name of input variable in geovals
    geovar = self%geovars%variable(ivar)

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

end subroutine destructor

! ------------------------------------------------------------------------------

end module ufo_atmvertinterp_tlad_mod
