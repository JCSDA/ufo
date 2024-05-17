! (C) Copyright 2017-2023 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

module ufo_vertinterp_tlad_mod

  use oops_variables_mod
  use obs_variables_mod
  use ufo_vars_mod
  use ufo_geovals_mod
  use vert_interp_mod
  use missing_values_mod
  use ufo_interp_param_mod

! ------------------------------------------------------------------------------

  type, public :: ufo_vertinterp_tlad
  private
    type(obs_variables), public :: obsvars ! Variables to be simulated
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

    logical, public :: hofx_scaling ! Apply scaling factor to hofx
    character(len=MAXVARLEN), public :: hofx_scaling_field
    character(len=MAXVARLEN), public :: hofx_scaling_field_group
    real(kind_real), allocatable :: scaling_field(:)

    ! Backup coordinate/method for interpolation
    logical :: use_backup_coordinate
    character(len=MAXVARLEN), public :: o_v_coord_backup     ! Obs vertical coordinate (backup)
    character(len=MAXVARLEN), public :: o_v_group_backup     ! Obs vertical coord group (backup)
    character(len=MAXVARLEN), public :: v_coord_backup       ! GeoVaL vert coordinate (backup)
    character(len=MAXVARLEN), public :: interp_method_backup ! Interpolation method (backup)
    integer, public :: selected_interp_backup

  contains
    procedure :: setup => vertinterp_tlad_setup_
    procedure :: cleanup => vertinterp_tlad_cleanup_
    procedure :: settraj => vertinterp_tlad_settraj_
    procedure :: simobs_tl => vertinterp_simobs_tl_
    procedure :: simobs_ad => vertinterp_simobs_ad_
    final :: destructor
  end type ufo_vertinterp_tlad

! ------------------------------------------------------------------------------
contains
! ------------------------------------------------------------------------------

subroutine vertinterp_tlad_setup_(self, grid_conf)
  use fckit_configuration_module, only: fckit_configuration
  implicit none
  class(ufo_vertinterp_tlad), intent(inout) :: self
  type(fckit_configuration), intent(in) :: grid_conf

  character(kind=c_char,len=:), allocatable :: coord_name
  character(kind=c_char,len=:), allocatable :: coord_group
  character(kind=c_char,len=:), allocatable :: interp_method
  character(kind=c_char,len=:), allocatable :: hofx_scaling_field
  character(kind=c_char,len=:), allocatable :: hofx_scaling_field_group
  integer :: ivar, nlevs
  character(len=MAXVARLEN) :: interp_method_backup

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

  !> Scale hofx by an incoming field. Can come from GeoVaLs or ObsSpace
  self%hofx_scaling = .false.
  if ( grid_conf%has("hofx scaling field") ) then
    self%hofx_scaling = .true.
    ! Get field name
    call grid_conf%get_or_die("hofx scaling field", hofx_scaling_field)
    self%hofx_scaling_field = hofx_scaling_field
    ! Get field name group
    self%hofx_scaling_field_group = "GeoVaLs"
    if ( grid_conf%has("hofx scaling field group") ) then
      call grid_conf%get_or_die("hofx scaling field group", hofx_scaling_field_group)
      self%hofx_scaling_field_group = hofx_scaling_field_group
    endif
    ! If the group is GeoVaLs then push back the variable name
    if (trim(self%hofx_scaling_field_group) == "GeoVaLs") then
      call self%geovars%push_back(trim(self%hofx_scaling_field))
    endif
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

  !> Look to see if the user wants to use a backup coordinate for the interpolation
  self%use_backup_coordinate = .false.
  if ( grid_conf%has("observation vertical coordinate backup") ) then

    ! Set flag to true
    self%use_backup_coordinate = .true.

    !> Use of a backup coordinate is not tested with also using self%use_constant_vcoord
    if (self%use_constant_vcoord) &
      call abor1_ftn('Requesting a backup coordinate in the vertical interpolation, but ' // &
                     'also using a constant vertical coordinate is not supported.')

    !> Get the name of the backup coordinate
    call grid_conf%get_or_die("observation vertical coordinate backup", coord_name)
    self%o_v_coord_backup = coord_name

    ! Set others to defaults
     self%o_v_group_backup     = self%o_v_group
     self%v_coord_backup       = self%v_coord
     self%interp_method_backup = self%interp_method

    !> Get group backup
    if ( grid_conf%has("observation vertical coordinate group backup") ) then
      call grid_conf%get_or_die("observation vertical coordinate group backup", coord_group)
      self%o_v_group_backup = coord_group
    endif

    !> Get model backgup coodinate
    if ( grid_conf%has("vertical coordinate backup") ) then
      call grid_conf%get_or_die("vertical coordinate backup", coord_name)
      self%v_coord_backup = coord_name
    endif

    !> Get interpolation method backup
    call grid_conf%get_or_die("interpolation method backup", interp_method)
    interp_method_backup = interp_method

    !> Linear interpolation is used by default.
    self%selected_interp_backup = LINEAR_INTERP
    if (trim(interp_method_backup) == "linear") then
      self%selected_interp_backup = LINEAR_INTERP
    else if (trim(interp_method_backup) == "log-linear") then
      self%selected_interp_backup = LOG_LINEAR_INTERP
    else if (trim(interp_method_backup) == "nearest-neighbor") then
      self%selected_interp_backup = NEAREST_NEIGHBOR_INTERP
    else
      !> the method is automatic
      if (trim(interp_method_backup) == "automatic") then
         !> Log-linear interpolation is used when v_coord is pressure
         if ((trim(self%v_coord_backup) .eq. var_prs) .or. &
             (trim(self%v_coord_backup) .eq. var_prsi) .or. &
             (trim(self%v_coord_backup) .eq. var_prsimo)) then
           self%selected_interp_backup = LOG_LINEAR_INTERP
         endif
      endif
    endif

    !> Assert that if nearest neighbor is chosen for the regular interpolation, then it is also
    !  chosen for the backup interpolation
    if ((self%selected_interp == NEAREST_NEIGHBOR_INTERP .and. &
         self%selected_interp_backup .ne. NEAREST_NEIGHBOR_INTERP) .or. &
        (self%selected_interp .ne. NEAREST_NEIGHBOR_INTERP .and. &
         self%selected_interp_backup == NEAREST_NEIGHBOR_INTERP)) &
      call abor1_ftn('If the regular interpolation method is nearest neighbor, then the ' // &
                     'backup interpolation method must also be nearest neighbor (and vice versa).')

  endif

end subroutine vertinterp_tlad_setup_

! ------------------------------------------------------------------------------

subroutine vertinterp_tlad_settraj_(self, geovals, obss)
  use missing_values_mod
  use obsspace_mod
  implicit none
  class(ufo_vertinterp_tlad), intent(inout) :: self
  type(ufo_geovals),         intent(in)    :: geovals
  type(c_ptr), value,        intent(in)    :: obss

  real(kind_real), allocatable :: obsvcoord(:), obsvcoord_backup(:)
  type(ufo_geoval), pointer :: vcoordprofile, vcoordprofile_backup
  integer :: ilev, iobs, nlevs
  real(kind_real), allocatable :: tmp(:)
  real(kind_real) :: tmp2
  real(kind_real) :: missing
  integer, allocatable :: selected_interp(:)
  ! Scaling by field
  type(ufo_geoval), pointer :: scaling_field_gval

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

  ! Turn selected interpolation into an array
  allocate(selected_interp(self%nlocs))
  selected_interp = self%selected_interp

  ! If using a backup coordinate for the interpolation, get the backup coordinate
  if (self%use_backup_coordinate) then

    ! Get the backup observation vertical coordinates
    allocate(obsvcoord_backup(self%nlocs))
    call obsspace_get_db(obss, self%o_v_group_backup, self%o_v_coord_backup, obsvcoord_backup)

    ! Get the backup coorindate from the model
    call ufo_geovals_get_var(geovals, self%v_coord_backup, vcoordprofile_backup)

    ! Loop over observations and use backup if necessary
    do iobs = 1, self%nlocs
      if (obsvcoord(iobs) == missing) then
        ! Use backup coordinate for this observation
        obsvcoord(iobs) = obsvcoord_backup(iobs)

        ! Use backup interpolation method for this observation
        selected_interp(iobs) = self%selected_interp_backup

        ! Use backup coordinate profile for this observation
        vcoordprofile%vals(:, iobs) = vcoordprofile_backup%vals(:, iobs)
      endif
    enddo
  endif

  do iobs = 1, self%nlocs
    if (.not. self%use_constant_vcoord) then
      if (selected_interp(iobs) == LOG_LINEAR_INTERP) then
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

    if (selected_interp(iobs) == LOG_LINEAR_INTERP) then
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

  ! Scaling to hofx
  if (self%hofx_scaling) then

    ! Get the scaling factor
    allocate(self%scaling_field(self%nlocs))
    if (trim(self%hofx_scaling_field_group) == "GeoVaLs") then
      call ufo_geovals_get_var(geovals, self%hofx_scaling_field, scaling_field_gval)
      self%scaling_field(:) = scaling_field_gval%vals(0, :)
    else
      call obsspace_get_db(obss, self%hofx_scaling_field_group, self%hofx_scaling_field, &
                           self%scaling_field)
    endif
  endif

  ! Cleanup memory
  deallocate(obsvcoord)
  deallocate(tmp)

end subroutine vertinterp_tlad_settraj_

! ------------------------------------------------------------------------------

subroutine vertinterp_simobs_tl_(self, geovals, obss, nvars, nlocs, hofx)
  implicit none
  class(ufo_vertinterp_tlad), intent(in)    :: self
  type(ufo_geovals),          intent(in)    :: geovals
  integer,                    intent(in)    :: nvars
  integer,                    intent(in)    :: nlocs
  real(c_double),             intent(inout) :: hofx(nvars, nlocs)
  type(c_ptr), value,         intent(in)    :: obss

  integer :: iobs, iobsvar, ivar
  real(c_double) :: missing
  type(ufo_geoval), pointer :: profile
  character(len=MAXVARLEN) :: geovar

  missing = missing_value(missing)

  ! Loop over the variables
  do iobsvar = 1, size(self%obsvarindices)

    ! Get the index of the row of hofx to fill
    ivar = self%obsvarindices(iobsvar)

    ! Get pointer to profile for this variable in geovals
    geovar = self%geovars%variable(iobsvar)
    call ufo_geovals_get_var(geovals, geovar, profile)

    ! Interpolate from geovals to observational location into hofx
    if (self%selected_interp == NEAREST_NEIGHBOR_INTERP) then
      do iobs = 1, nlocs
        if (hofx(ivar,iobs) /= missing) then
          call nearestneighbor_interp_apply_tl(profile%nval, profile%vals(:,iobs), &
                                             & hofx(ivar,iobs), self%wi(iobs))
        end if
      enddo
    else
      do iobs = 1, nlocs
        if (hofx(ivar,iobs) /= missing) then
          call vert_interp_apply_tl(profile%nval, profile%vals(:,iobs), &
                                  & hofx(ivar,iobs), self%wi(iobs), self%wf(iobs))
        end if
      enddo
    end if

    ! Scaling to hofx
    if (self%hofx_scaling) then
      do iobs = 1, nlocs
        if (hofx(ivar,iobs) /= missing) then
          hofx(ivar,iobs) = hofx(ivar,iobs) * self%scaling_field(iobs)
        endif
      enddo
    endif

  enddo

end subroutine vertinterp_simobs_tl_

! ------------------------------------------------------------------------------

subroutine vertinterp_simobs_ad_(self, geovals, obss, nvars, nlocs, hofx_in)
  implicit none
  class(ufo_vertinterp_tlad), intent(in)    :: self
  type(ufo_geovals),          intent(inout) :: geovals
  integer,                    intent(in)    :: nvars
  integer,                    intent(in)    :: nlocs
  real(c_double),             intent(in)    :: hofx_in(nvars, nlocs)
  type(c_ptr), value,         intent(in)    :: obss

  integer :: iobs, iobsvar, ivar
  type(ufo_geoval), pointer :: profile
  character(len=MAXVARLEN) :: geovar
  real(c_double) :: missing
  real(c_double) :: hofx(nvars, nlocs)

  missing = missing_value(missing)

  ! Copy hofx_in to hofx, in case it needs to be scaled
  hofx = hofx_in

  ! Loop over the variables
  do iobsvar = 1, size(self%obsvarindices)

    ! Get the index of the row of hofx to fill
    ivar = self%obsvarindices(iobsvar)

    ! Get pointer to profile for this variable in geovals
    geovar = self%geovars%variable(iobsvar)
    call ufo_geovals_get_var(geovals, geovar, profile)

    ! Scaling to hofx
    if (self%hofx_scaling) then
      do iobs = 1, nlocs
        if (hofx(ivar,iobs) /= missing) then
          hofx(ivar,iobs) = hofx(ivar,iobs) * self%scaling_field(iobs)
        endif
      enddo
    endif

    ! Adjoint of interpolatation
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

end subroutine vertinterp_simobs_ad_

! ------------------------------------------------------------------------------

subroutine vertinterp_tlad_cleanup_(self)
  implicit none
  class(ufo_vertinterp_tlad), intent(inout) :: self
  self%nval = 0
  self%nlocs = 0
  if (allocated(self%wi)) deallocate(self%wi)
  if (allocated(self%wf)) deallocate(self%wf)
  if (allocated(self%scaling_field)) deallocate(self%scaling_field)
end subroutine vertinterp_tlad_cleanup_

! ------------------------------------------------------------------------------

subroutine  destructor(self)
  type(ufo_vertinterp_tlad), intent(inout)  :: self

  call self%cleanup()

end subroutine destructor

! ------------------------------------------------------------------------------

end module ufo_vertinterp_tlad_mod
