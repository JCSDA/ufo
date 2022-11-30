! (C) Copyright 2021 Met Office UK
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

module ufo_backgrounderrorvertinterp_mod

use iso_c_binding,      only: c_ptr
use oops_variables_mod, only: oops_variables
use ufo_geovals_mod,    only: ufo_geovals

contains

!> For each obs diagnostic called <var>_background_error, where <var> belongs to the set of variable
!> names @p obsvars, fill this diagnostic with estimates of the background error of variable <var>
!> at observation locations.
subroutine ufo_backgrounderrorvertinterp_fillobsdiags(obs_vcoord_name, obs_vcoord_group, vcoord_name, &
                                                      geovals, obsspace, nlocs, obsvars, obsdiags)
  use kinds,              only: kind_real
  use missing_values_mod, only: missing_value
  use obsspace_mod,       only: obsspace_get_db
  use ufo_vars_mod,       only: MAXVARLEN, var_prs, var_prsi, var_prsimo
  use ufo_geovals_mod,    only: ufo_geoval, ufo_geovals, ufo_geovals_get_var
  use vert_interp_mod,    only: vert_interp_weights, vert_interp_apply
  implicit none

  ! Name of the variable with vertical coordinates of observations
  character(len=*), intent(in)     :: obs_vcoord_name
  ! Group of the variable with vertical coordinates of observations
  character(len=*), intent(in)     :: obs_vcoord_group
  ! Name of the GeoVaL with the vertical coordinate levels to use for
  ! interpolation of background errors
  character(len=*), intent(in)     :: vcoord_name
  type(ufo_geovals), intent(in)    :: geovals
  type(c_ptr), value, intent(in)   :: obsspace
  integer, intent(in)              :: nlocs
  type(oops_variables), intent(in) :: obsvars
  type(ufo_geovals), intent(inout) :: obsdiags

  logical                          :: use_ln
  integer                          :: iobs, ivar
  real(kind_real)                  :: obs_vcoord(nlocs)
  type(ufo_geoval), pointer        :: vcoord_profile, background_error_profile
  real(kind_real)                  :: wf(nlocs)
  integer                          :: wi(nlocs)
  character(len=MAXVARLEN)         :: varstr
  integer                          :: lenvarstr
  real(kind_real), allocatable     :: interp_nodes(:)
  real(kind_real)                  :: interp_point
  real(kind_real)                  :: missing

  character(len=*), parameter      :: suffix = "_background_error"

  ! Get vertical coordinate profiles from geovals
  call ufo_geovals_get_var(geovals, vcoord_name, vcoord_profile)

  ! Get the observation vertical coordinates
  call obsspace_get_db(obsspace, obs_vcoord_group, obs_vcoord_name, obs_vcoord)

  ! Set missing value
  if (nlocs > 0) then
     missing = missing_value(obs_vcoord(1))
  end if

  ! Use logarithmic interpolation if the vertical coordinate is background_error_air_pressure
  use_ln = (vcoord_name .eq. "background_error_" // var_prs)

  ! Calculate the interpolation weights
  allocate(interp_nodes(vcoord_profile%nval))
  do iobs = 1, nlocs
    if (use_ln) then
      interp_nodes = log(vcoord_profile%vals(:,iobs))
      if (obs_vcoord(iobs) /= missing) then
         interp_point = log(obs_vcoord(iobs))
      else
         interp_point = missing
      end if
    else
      interp_nodes = vcoord_profile%vals(:,iobs)
      interp_point = obs_vcoord(iobs)
    end if
    call vert_interp_weights(vcoord_profile%nval, interp_point, interp_nodes, &
                             wi(iobs), wf(iobs))
  enddo

  do ivar = 1, obsdiags%nvar
    varstr = obsdiags%variables(ivar)    
    lenvarstr = len_trim(varstr)

    ! We need to fill this diagnostic if:
    ! (a) its name is long enough to be of the form `<var>_background_error`;
    if (lenvarstr <= len(suffix)) cycle
    ! (b) its name actually *is* of the form `<var>_background_error`;
    if (varstr(lenvarstr - len(suffix)+1:lenvarstr) /= suffix) cycle
    ! (c) <var> belongs to the list obsvars.
    if (.not. obsvars%has(varstr(:lenvarstr - len(suffix)))) cycle

    ! All tests passed -- fill the diagnostic.

    ! Get background error profile from geovals
    call ufo_geovals_get_var(geovals, varstr, background_error_profile)

    ! Allocate the background error diagnostic.
    if (allocated(obsdiags%geovals(ivar)%vals)) deallocate(obsdiags%geovals(ivar)%vals)
    obsdiags%geovals(ivar)%nval = 1
    allocate(obsdiags%geovals(ivar)%vals(obsdiags%geovals(ivar)%nval, nlocs))

    ! Interpolate the profile at observation location into the obsdiag
    do iobs = 1, nlocs
      call vert_interp_apply(background_error_profile%nval, &
                             background_error_profile%vals(:,iobs), &
                             obsdiags%geovals(ivar)%vals(1,iobs), &
                             wi(iobs), wf(iobs))
    enddo
  enddo

  ! Free memory
  deallocate(interp_nodes)
end subroutine ufo_backgrounderrorvertinterp_fillobsdiags

end module ufo_backgrounderrorvertinterp_mod
