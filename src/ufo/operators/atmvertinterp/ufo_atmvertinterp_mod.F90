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
                                                      ! simulated variables in the ObsSpace.
                                                      ! allocated/deallocated at interface layer
     character(len=MAXVARLEN), public :: v_coord ! GeoVaL to use to interpolate in vertical
     character(len=MAXVARLEN), public :: o_v_coord ! Observation vertical coordinate
     character(len=MAXVARLEN), public :: o_v_group ! Observation vertical coordinate group
     character(len=MAXVARLEN), public :: interp_method ! Vertical interpolation method

     logical, public :: use_ln ! if T, use ln(v_coord) not v_coord
     logical, public :: use_fact10 ! Apply scaling factor to winds below lowest model level
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

  !> Apply scaling to winds below lowest model level
  self%use_fact10 = .false.
  if ( grid_conf%has("apply near surface wind scaling") ) then
    call grid_conf%get_or_die("apply near surface wind scaling", self%use_fact10)
  endif
  if (self%use_fact10) call self%geovars%push_back("wind_reduction_factor_at_10m")

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

  integer :: ilev, iobs, ivar, iobsvar
  real(kind_real), dimension(:), allocatable :: obsvcoord
  type(ufo_geoval), pointer :: vcoordprofile, profile, fact10
  real(kind_real), allocatable :: wf(:)
  integer, allocatable :: wi(:)
  character(len=MAXVARLEN) :: geovar

  real(kind_real), allocatable :: tmp(:)
  real(kind_real) :: tmp2
  real(kind_real) :: missing

  real(kind_real), allocatable :: wind_scaling_factor(:)

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

  ! If scaling the wind record the observation pressure and get GeoVaLs
  if (self%use_fact10) then
    allocate(wind_scaling_factor(nlocs))
    wind_scaling_factor = 1.0_kind_real
    call ufo_geovals_get_var(geovals, "wind_reduction_factor_at_10m", fact10)
  end if

  ! Calculate the interpolation weights
  allocate(tmp(vcoordprofile%nval))
  do iobs = 1, nlocs
    if (self%use_ln) then
      ! the lines below are computing a "missing value safe" log, that passes missing value inputs
      ! through to the output. the simpler "tmp = log(rhs)" produces NaN for missing value inputs.
      do ilev = 1, vcoordprofile%nval
        if (vcoordprofile%vals(ilev,iobs) /= missing) then
          tmp(ilev) = log(vcoordprofile%vals(ilev,iobs))
        else
          tmp(ilev) = missing
        end if
      end do
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

    ! Set scaling factor
    if (self%use_fact10) then
      if (tmp2 >= tmp(1)) wind_scaling_factor(iobs) = fact10%vals(1,iobs)
    end if
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

  ! Apply a scaling to winds below lowest model level
  if (self%use_fact10) then
    ! Loop over the variables
    do iobsvar = 1, size(self%obsvarindices)
      ! Check that this is a typical wind variable
      if ((trim(self%obsvars%variable(iobsvar)) == 'eastward_wind') .or. &
          (trim(self%obsvars%variable(iobsvar)) == 'northward_wind')) then
        ! Get the index of the row of hofx to fill
        ivar = self%obsvarindices(iobsvar)
        ! Loop over the observations
        do iobs = 1, nlocs
          ! Apply wind scaling
          hofx(ivar,iobs) = hofx(ivar,iobs) * wind_scaling_factor(iobs)
        enddo
      end if
    enddo
  endif

  ! Cleanup memory
  deallocate(obsvcoord)
  deallocate(wi)
  deallocate(wf)

  deallocate(tmp)

  if (allocated(wind_scaling_factor)) deallocate(wind_scaling_factor)

end subroutine atmvertinterp_simobs_

! ------------------------------------------------------------------------------

end module ufo_atmvertinterp_mod
