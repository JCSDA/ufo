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
    procedure :: setup  => conventional_profile_setup_
    procedure :: simobs => conventional_profile_simobs_
    procedure :: locateobs => conventional_profile_locateobs_

    final :: destructor
  end type ufo_conventional_profile

! ------------------------------------------------------------------------------
contains
! ------------------------------------------------------------------------------

subroutine conventional_profile_setup_(self, c_conf)
  use config_mod
  implicit none
  class(ufo_conventional_profile), intent(inout) :: self
  type(c_ptr), intent(in)    :: c_conf

  integer :: ii

  if (config_element_exists(c_conf,"variables")) then
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
       if (trim(self%varout(ii)) .eq. "air_temperature") then
         self%varin(ii) = "virtual_temperature"
       else
         self%varin(ii) = self%varout(ii)
       endif
    enddo
    !> Put log pressure to the varin (vars from the model) list
    self%varin(self%nvars+1) = "atmosphere_ln_pressure_coordinate"
  endif

end subroutine conventional_profile_setup_

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

  do ivar = 1, self%nvars
    ! Get the name of input variable in geovals
    geovar = self%varin(ivar)

    ! Get profile for this variable from geovals
    call ufo_geovals_get_var(geovals, geovar, profile)

    ! Interpolate from geovals to observational location into hofx
    ! Note: hofx holds all variables (varin) for location 1
    ! then all variables for location 2, and so on
    do iobs = 1, nlocs
      call vert_interp_apply(profile%nval, profile%vals(:,iobs), &
                             & hofx(ivar + (iobs-1)*self%nvars), wi(iobs), wf(iobs))
    enddo
  enddo
  ! Cleanup memory
  deallocate(obspressure)
  deallocate(wi)
  deallocate(wf)

end subroutine conventional_profile_simobs_

! ------------------------------------------------------------------------------

subroutine conventional_profile_locateobs_(self, obss, t1, t2, locs)
  use datetime_mod
  use twindow_utils_mod
  use fckit_log_module, only : fckit_log
  use ufo_locs_mod, only: ufo_locs, ufo_locs_setup

  implicit none

  class(ufo_conventional_profile), intent(in) :: self
  type(c_ptr), value, intent(in)              :: obss
  type(datetime), intent(in)                  :: t1, t2
  type(ufo_locs), intent(inout)               :: locs

  integer :: nlocs
  type(datetime) :: refdate

  character(len=*),parameter:: &
     myname = "conventional_profile_locateobs_"
  character(len=255) :: record
  integer :: i
  integer :: tw_nlocs
  integer, dimension(:), allocatable :: tw_indx
  real(kind_real), dimension(:), allocatable :: time, lon, lat

  ! Local copies pre binning
  nlocs = obsspace_get_nlocs(obss)
  refdate = obsspace_get_refdate(obss)

  allocate(time(nlocs), lon(nlocs), lat(nlocs))

  !!Each operator may have its own way to derive time, lon, lat from MetaData
  !!BEGIN THIS PART CAN BE UNIQUE FOR SOME OBS OPERATORS
  call obsspace_get_db(obss, "MetaData", "time", time)

  ! Generate the timing window indices
  allocate(tw_indx(nlocs))
  call gen_twindow_index(refdate, t1, t2, nlocs, time, tw_indx, tw_nlocs)

  call obsspace_get_db(obss, "MetaData", "longitude", lon)
  call obsspace_get_db(obss, "MetaData", "latitude", lat)
  !!END THIS PART CAN BE UNIQUE FOR SOME OBS OPERATORS

  !Setup ufo locations
  call ufo_locs_setup(locs, tw_nlocs)
  do i = 1, tw_nlocs
    locs%lon(i)  = lon(tw_indx(i))
    locs%lat(i)  = lat(tw_indx(i))
    locs%time(i) = time(tw_indx(i))
  enddo
  locs%indx = tw_indx(1:tw_nlocs)

  deallocate(time, lon, lat, tw_indx)


  write(record,*) myname,': allocated/assigned obs locations'
  call fckit_log%info(record)

end subroutine conventional_profile_locateobs_

! ------------------------------------------------------------------------------

subroutine  destructor(self)
  type(ufo_conventional_profile), intent(inout) :: self
  if (allocated(self%varout)) deallocate(self%varout)
  if (allocated(self%varin)) deallocate(self%varin)
end subroutine destructor

! ------------------------------------------------------------------------------

end module ufo_conventional_profile_mod
