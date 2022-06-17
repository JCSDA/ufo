! (C) Copyright 2017-2019 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

module ufo_atmvertinterp_mod

use oops_variables_mod
use ufo_vars_mod
! ------------------------------------------------------------------------------

  type, public :: ufo_atmvertinterp
     type(oops_variables), public :: geovars
     type(oops_variables), public :: obsvars ! Variables to be simulated
     integer, allocatable, public :: obsvarindices(:) ! Indices of obsvars in the list of all
                                                      ! simulated variables in the ObsSpace
     character(len=MAXVARLEN), public :: v_coord ! GeoVaL to use to interpolate in vertical
     character(len=MAXVARLEN), public :: o_v_coord ! Observation vertical coordinate
     character(len=MAXVARLEN), public :: o_v_group ! Observation vertical coordinate group
     character(len=MAXVARLEN), public :: interp_method ! Vertical interpolation method

     logical, public :: use_ln ! if T, use ln(v_coord) not v_coord
   contains
     procedure :: setup  => atmvertinterp_setup_
     procedure :: simobs => atmvertinterp_simobs_
  end type ufo_atmvertinterp

! ------------------------------------------------------------------------------
contains
! ------------------------------------------------------------------------------

subroutine atmvertinterp_setup_(self, grid_conf)
  use iso_c_binding
  use fckit_configuration_module, only: fckit_configuration
  implicit none
  class(ufo_atmvertinterp), intent(inout) :: self
  type(fckit_configuration), intent(in)   :: grid_conf

  character(kind=c_char,len=:), allocatable :: coord_name
  character(kind=c_char,len=:), allocatable :: coord_group
  character(kind=c_char,len=:), allocatable :: interp_method
  integer :: ivar, nvars

  !> Size of variables
  nvars = self%obsvars%nvars()
  !> Fill in geovars: variables we need from the model
  !  need additional slot to hold vertical coord.
  do ivar = 1, nvars
    call self%geovars%push_back(self%obsvars%variable(ivar))
  enddo
  !> grab what vertical coordinate/variable to use from the config
  call grid_conf%get_or_die("vertical coordinate",coord_name)
  self%v_coord = coord_name

  call grid_conf%get_or_die("interpolation method",interp_method)
  self%interp_method = interp_method

  !> Linear interpolation is used by default.
  self%use_ln = .false.
  !> Log-linear interpolation is used either if it is explicitly requested
  !  or the method is automatically determined based on the vertical coordinate used.
  if ((trim(self%interp_method) == "automatic" .and. &
       ((trim(self%v_coord) .eq. var_prs) .or. &
        (trim(self%v_coord) .eq. var_prsi) .or. &
        (trim(self%v_coord) .eq. var_prsimo))) .or. &
      (trim(self%interp_method) == "log-linear")) then
     self%use_ln = .true.
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

  call self%geovars%push_back(self%v_coord)

end subroutine atmvertinterp_setup_

! ------------------------------------------------------------------------------

subroutine atmvertinterp_simobs_(self, geovals, obss, nvars, nlocs, hofx)
  use kinds
  use missing_values_mod
  use obsspace_mod
  use vert_interp_mod
  use ufo_geovals_mod
  implicit none
  class(ufo_atmvertinterp), intent(in)        :: self
  integer, intent(in)                         :: nvars, nlocs
  type(ufo_geovals), intent(in)               :: geovals
  real(c_double),  intent(inout)              :: hofx(nvars, nlocs)
  type(c_ptr), value, intent(in)              :: obss

  integer :: iobs, ivar, iobsvar
  real(kind_real), dimension(:), allocatable :: obsvcoord
  type(ufo_geoval), pointer :: vcoordprofile, profile
  real(kind_real), allocatable :: wf(:)
  integer, allocatable :: wi(:)
  character(len=MAXVARLEN) :: geovar

  real(kind_real), allocatable :: tmp(:)
  real(kind_real) :: tmp2
  real(kind_real) :: missing

  ! Get pressure profiles from geovals
  call ufo_geovals_get_var(geovals, self%v_coord, vcoordprofile)

  ! Get the observation vertical coordinates
  allocate(obsvcoord(nlocs))
  call obsspace_get_db(obss, self%o_v_group, self%o_v_coord, obsvcoord)

  ! Set missing value
  if (nlocs > 0) then
     missing = missing_value(obsvcoord(1))
  end if

  ! Allocate arrays for interpolation weights
  allocate(wi(nlocs))
  allocate(wf(nlocs))

  ! Calculate the interpolation weights
  allocate(tmp(vcoordprofile%nval))
  do iobs = 1, nlocs
    if (self%use_ln) then
      tmp = log(vcoordprofile%vals(:,iobs))
      if (obsvcoord(iobs) /= missing) then
         tmp2 = log(obsvcoord(iobs))
      else
         tmp2 = missing
      end if
    else
      tmp = vcoordprofile%vals(:,iobs)
      tmp2 = obsvcoord(iobs)
    end if
    call vert_interp_weights(vcoordprofile%nval, tmp2, tmp, wi(iobs), wf(iobs))
  enddo

  do iobsvar = 1, size(self%obsvarindices)
    ! Get the index of the row of hofx to fill
    ivar = self%obsvarindices(iobsvar)

    ! Get the name of input variable in geovals
    geovar = self%geovars%variable(iobsvar)

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

end module ufo_atmvertinterp_mod