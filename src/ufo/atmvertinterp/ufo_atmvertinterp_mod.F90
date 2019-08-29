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
  use fckit_configuration_module, only: fckit_configuration

  integer, parameter :: max_string=800

! ------------------------------------------------------------------------------

  type, public :: ufo_atmvertinterp
   private
     integer :: nvars  ! number of variables to be interpolated
     character(len=max_string), public, allocatable :: varin(:)    ! size nvars+1 (+1 for log pressure)
     character(len=max_string), public, allocatable :: varout(:)   ! size nvars
     character(len=MAXVARLEN), public :: v_coord ! GeoVaL to use to interpolate in vertical
     logical, public :: use_ln ! if T, use ln(v_coord) not v_coord
   contains
     procedure :: setup  => atmvertinterp_setup_
     procedure :: simobs => atmvertinterp_simobs_
     final :: destructor
  end type ufo_atmvertinterp

! ------------------------------------------------------------------------------
contains
! ------------------------------------------------------------------------------

subroutine atmvertinterp_setup_(self, grid_conf, vars)
  implicit none
  class(ufo_atmvertinterp), intent(inout) :: self
  character(len=MAXVARLEN), dimension(:), intent(inout) :: vars
  type(fckit_configuration) :: grid_conf
  character(kind=c_char,len=:), allocatable :: coord_name

  !> Size of variables
  self%nvars = size(vars)
  !> Allocate varout: variables in the observation vector
  allocate(self%varout(self%nvars))
  !> Read variable list and store in varout
  self%varout = vars
  !> Allocate varin: variables we need from the model
  !  need additional slot to hold vertical coord.
  allocate(self%varin(self%nvars+1))
  !> Set vars_in based on vars_out
  self%varin(1:self%nvars) = self%varout(1:self%nvars)

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

  self%varin(self%nvars+1) = self%v_coord

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
  real(kind_real), dimension(:), allocatable :: obsvcoord
  type(ufo_geoval), pointer :: vcoordprofile, profile
  real(kind_real), allocatable :: wf(:)
  integer, allocatable :: wi(:)
  character(len=MAXVARLEN) :: geovar

  real(kind_real), allocatable :: tmp(:)
  real(kind_real) :: tmp2

  ! Get pressure profiles from geovals
  call ufo_geovals_get_var(geovals, self%v_coord, vcoordprofile)

  ! Get the observation vertical coordinates
  allocate(obsvcoord(nlocs))
  call obsspace_get_db(obss, "MetaData", self%v_coord, obsvcoord)

  ! Allocate arrays for interpolation weights
  allocate(wi(nlocs))
  allocate(wf(nlocs))

  ! Calculate the interpolation weights
  allocate(tmp(vcoordprofile%nval))
  do iobs = 1, nlocs
    if (self%use_ln) then
      tmp = log(vcoordprofile%vals(:,iobs))
      tmp2 = log(obsvcoord(iobs))
    else
      tmp = vcoordprofile%vals(:,iobs)
      tmp2 = obsvcoord(iobs)
    end if
    call vert_interp_weights(vcoordprofile%nval, tmp2, tmp, wi(iobs), wf(iobs))
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
  deallocate(obsvcoord)
  deallocate(wi)
  deallocate(wf)

  deallocate(tmp)

end subroutine atmvertinterp_simobs_

! ------------------------------------------------------------------------------

subroutine destructor(self)
  type(ufo_atmvertinterp), intent(inout) :: self
  if (allocated(self%varout)) deallocate(self%varout)
  if (allocated(self%varin)) deallocate(self%varin)
end subroutine destructor

! ------------------------------------------------------------------------------

end module ufo_atmvertinterp_mod
