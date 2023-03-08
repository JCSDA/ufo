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
  use ufo_interp_param_mod

! ------------------------------------------------------------------------------

  type, public :: ufo_atmvertinterp_tlad
  private
    type(oops_variables), public :: obsvars ! Variables to be simulated
    integer, allocatable, public :: obsvarindices(:) ! Indices of obsvars in the list of all
                                                     ! simulated variables in the ObsSpace.
                                                     ! allocated/deallocated at interface layer
    logical, public :: use_constant_vcoord ! if T, use constant vertical coordinate specified
                                           ! in configuration instead of geoval
    real, allocatable, public :: const_v_coord(:) ! if use_constant_vcoord, holds values of
                                                  ! constant vertical coordinate
    type(oops_variables), public :: geovars
    integer :: nval, nlocs
    real(kind_real), allocatable :: wf(:)
    integer, allocatable :: wi(:)
    character(len=MAXVARLEN), public :: v_coord ! GeoVaL to use to interpolate in vertical
    character(len=MAXVARLEN), public :: o_v_coord ! Observation vertical coordinate
    character(len=MAXVARLEN), public :: o_v_group ! Observation vertical coordinate group
    character(len=MAXVARLEN), public :: interp_method ! Vertical interpolation method

    integer, public :: selected_interp
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
  character(kind=c_char,len=:), allocatable :: coord_group
  character(kind=c_char,len=:), allocatable :: interp_method
  integer :: ivar, nlevs

  !> grab what vertical coordinate/variable to use from the config
  !> check if constant vertical coordinate is provided
  self%use_constant_vcoord = grid_conf%has("constant vertical coordinate values")
  if (self%use_constant_vcoord) then
    nlevs = grid_conf%get_size("constant vertical coordinate values")
    allocate(self%const_v_coord(nlevs))
    call grid_conf%get_or_die("constant vertical coordinate values", self%const_v_coord)
  !> if constant values aren't provided, get geoval name for vertical coordinate
  else
    call grid_conf%get_or_die("vertical coordinate",coord_name)
    self%v_coord = coord_name
  endif

  call grid_conf%get_or_die("interpolation method",interp_method)
  self%interp_method = interp_method

  !> Linear interpolation is used by default.
  self%selected_interp = LINEAR_INTERP
  if(trim(self%interp_method) == "linear") then
    self%selected_interp = LINEAR_INTERP
  else if(trim(self%interp_method) == "log-linear") then
    self%selected_interp = LOG_LINEAR_INTERP
  else if(trim(self%interp_method) == "nearest-neighbor") then
    self%selected_interp = NEAREST_NEIGHBOR_INTERP
  else
    !> the method is automatic
    if (trim(self%interp_method) == "automatic") then
       !> Log-linear interpolation is used when v_coord is pressure
       if ((trim(self%v_coord) .eq. var_prs) .or. &
           (trim(self%v_coord) .eq. var_prsi) .or. &
           (trim(self%v_coord) .eq. var_prsimo)) then
         self%selected_interp = LOG_LINEAR_INTERP
       !> Nearest-Neighbor is used when const vertical coordinate used.
       else if (self%use_constant_vcoord) then
         self%selected_interp = NEAREST_NEIGHBOR_INTERP
       endif
    endif
  endif

  !> Determine observation vertical coordinate.
  !  Use the model vertical coordinate unless the option
  !  'observation vertical coordinate' is specified.
  if ( grid_conf%has("observation vertical coordinate") ) then
     call grid_conf%get_or_die("observation vertical coordinate",coord_name)
     self%o_v_coord = coord_name
  else
     self%o_v_coord = self%v_coord
  endif

  !> Determine observation vertical coordinate group.
  !  Use MetaData unless the option
  !  'observation vertical coordinate' is specified.
  if ( grid_conf%has("observation vertical coordinate group") ) then
    call grid_conf%get_or_die("observation vertical coordinate group",coord_group)
    self%o_v_group = coord_group
  else
    self%o_v_group = "MetaData"
  endif

end subroutine atmvertinterp_tlad_setup_

! ------------------------------------------------------------------------------

subroutine atmvertinterp_tlad_settraj_(self, geovals, obss)
  use missing_values_mod
  use obsspace_mod
  implicit none
  class(ufo_atmvertinterp_tlad), intent(inout) :: self
  type(ufo_geovals),         intent(in)    :: geovals
  type(c_ptr), value,        intent(in)    :: obss

  real(kind_real), allocatable :: obsvcoord(:)
  type(ufo_geoval), pointer :: vcoordprofile
  integer :: ilev, iobs, nlevs
  real(kind_real), allocatable :: tmp(:)
  real(kind_real) :: tmp2
  real(kind_real) :: missing

  ! Make sure nothing already allocated
  call self%cleanup()

  ! Get pressure profiles from geovals
  call ufo_geovals_get_var(geovals, self%v_coord, vcoordprofile)
  self%nval = vcoordprofile%nval

  ! Get the observation vertical coordinates
  self%nlocs = obsspace_get_nlocs(obss)
  allocate(obsvcoord(self%nlocs))
  call obsspace_get_db(obss, self%o_v_group, self%o_v_coord, obsvcoord)

  ! Set missing value
  if (self%nlocs > 0) then
     missing = missing_value(obsvcoord(1))
  end if

  ! Allocate arrays for interpolation weights
  allocate(self%wi(self%nlocs))
  allocate(self%wf(self%nlocs))

  ! Calculate the interpolation weights
  if (self%use_constant_vcoord) then
    nlevs = size(self%const_v_coord)
    allocate(tmp(nlevs))
    tmp = self%const_v_coord
    if (self%selected_interp == LOG_LINEAR_INTERP) then
      do ilev = 1, nlevs
        tmp(ilev) = log(tmp(ilev))
      enddo
    endif
  else
    nlevs = vcoordprofile%nval
    allocate(tmp(vcoordprofile%nval))
  endif

  do iobs = 1, self%nlocs
    if (.not. self%use_constant_vcoord) then
      if (self%selected_interp == LOG_LINEAR_INTERP) then
        ! the lines below are computing a "missing value safe" log, that passes missing value inputs
        ! through to the output. the simpler "tmp = log(rhs)" produces NaN for missing value inputs.
        do ilev = 1, vcoordprofile%nval
          if (vcoordprofile%vals(ilev,iobs) /= missing) then
            tmp(ilev) = log(vcoordprofile%vals(ilev,iobs))
          else
            tmp(ilev) = missing
          end if
        end do
      else
        tmp = vcoordprofile%vals(:,iobs)
      endif
    endif

    if (self%selected_interp == LOG_LINEAR_INTERP) then
      if (obsvcoord(iobs) /= missing) then
         tmp2 = log(obsvcoord(iobs))
      else
         tmp2 = missing
      end if
    else
      tmp2 = obsvcoord(iobs)
    end if

    if (self%selected_interp == NEAREST_NEIGHBOR_INTERP) then
      call nearestneighbor_interp_index(nlevs, tmp2, tmp, self%wi(iobs))
    else
      call vert_interp_weights(nlevs, tmp2, tmp, self%wi(iobs), self%wf(iobs))
    end if
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

  integer :: iobs, iobsvar, ivar
  type(ufo_geoval), pointer :: profile
  character(len=MAXVARLEN) :: geovar

  do iobsvar = 1, size(self%obsvarindices)
    ! Get the index of the row of hofx to fill
    ivar = self%obsvarindices(iobsvar)

    ! Get the name of input variable in geovals
    geovar = self%geovars%variable(iobsvar)

    ! Get profile for this variable from geovals
    call ufo_geovals_get_var(geovals, geovar, profile)

    ! Interpolate from geovals to observational location into hofx
    if (self%selected_interp == NEAREST_NEIGHBOR_INTERP) then
      do iobs = 1, nlocs
        call nearestneighbor_interp_apply_tl(profile%nval, profile%vals(:,iobs), &
                                           & hofx(ivar,iobs), self%wi(iobs))
      enddo
    else
      do iobs = 1, nlocs
        call vert_interp_apply_tl(profile%nval, profile%vals(:,iobs), &
                                & hofx(ivar,iobs), self%wi(iobs), self%wf(iobs))
      enddo
    end if
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

  integer :: iobs, iobsvar, ivar
  type(ufo_geoval), pointer :: profile
  character(len=MAXVARLEN) :: geovar
  real(c_double) :: missing

  missing = missing_value(missing)

  do iobsvar = 1, size(self%obsvarindices)
    ! Get the index of the row of hofx to fill
    ivar = self%obsvarindices(iobsvar)

    ! Get the name of input variable in geovals
    geovar = self%geovars%variable(iobsvar)

    ! Get pointer to profile for this variable in geovals
    call ufo_geovals_get_var(geovals, geovar, profile)

    ! Adjoint of interpolate, from hofx into geovals
    ! Interpolate from geovals to observational location into hofx
    if (self%selected_interp == NEAREST_NEIGHBOR_INTERP) then
      do iobs = 1, nlocs
        if (hofx(ivar,iobs) /= missing) then
          call nearestneighbor_interp_apply_ad(profile%nval, profile%vals(:,iobs), &
                                             & hofx(ivar,iobs), self%wi(iobs))
        endif
      enddo
    else
      do iobs = 1, nlocs
        if (hofx(ivar,iobs) /= missing) then
          call vert_interp_apply_ad(profile%nval, profile%vals(:,iobs), &
                                  & hofx(ivar,iobs), self%wi(iobs), self%wf(iobs))
        endif
      enddo
    end if
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
