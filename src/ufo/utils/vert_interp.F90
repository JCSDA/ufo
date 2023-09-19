! (C) Copyright 2018 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

!> Fortran module to perform linear interpolation

module vert_interp_mod

use, intrinsic :: iso_c_binding
use kinds, only: kind_real
use missing_values_mod
use oops_variables_mod

implicit none
public

contains

! ------------------------------------------------------------------------------

subroutine vert_interp_weights(nlev,obl,vec,wi,wf)

integer,         intent(in ) :: nlev       !Number of model levels
real(kind_real), intent(in ) :: obl        !Observation location
real(kind_real), intent(in ) :: vec(nlev)  !Structured vector of grid points
integer,         intent(out) :: wi         !Index for interpolation
real(kind_real), intent(out) :: wf         !Weight for interpolation

integer         :: k
real(kind_real) :: missing

missing = missing_value(obl)

! If the observation is missing then set both the index and weight to missing.
if (obl == missing) then
   wi = missing_value(nlev)
   wf = missing
   return
end if

if (vec(1) < vec(nlev)) then !Pressure increases with index

  if (obl < vec(1)) then
     wi = 1
     wf = 1.0
  elseif (obl > vec(nlev)) then
     wi = nlev - 1
     wf = 0.0
  else
     do k = 1,nlev-1
        if (obl >= vec(k) .and. obl <= vec(k+1)) then
           wi = k
        endif
     enddo
     wf = (vec(wi+1) - obl)/(vec(wi+1) - vec(wi))
  endif

else !Pressure decreases with index

  if (obl > vec(1)) then
     wi = 1
     wf = 1.0
  elseif (obl < vec(nlev)) then
     wi = nlev - 1
     wf = 0.0
  else
     do k = 1,nlev-1
        if (obl >= vec(k+1) .and. obl <= vec(k)) then
           wi = k
        endif
     enddo
     wf = (vec(wi+1) - obl)/(vec(wi+1) - vec(wi))
  endif

endif

end subroutine vert_interp_weights

! ------------------------------------------------------------------------------

subroutine vert_interp_apply(nlev, fvec, f, wi, wf)

integer,         intent(in ) :: nlev        !Number of model levels
real(kind_real), intent(in ) :: fvec(nlev)  !Field at grid points
integer,         intent(in ) :: wi          !Index for interpolation
real(kind_real), intent(in ) :: wf          !Weight for interpolation
real(kind_real), intent(out) :: f           !Output at obs location using linear interp

if (wi == missing_value(nlev)) then
  f = missing_value(f)
elseif (fvec(wi) == missing_value(f) .or. fvec(wi+1) == missing_value(f)) then
  f = missing_value(f)
else
  f = fvec(wi)*wf + fvec(wi+1)*(1.0-wf)
endif

end subroutine vert_interp_apply

! ------------------------------------------------------------------------------

subroutine vert_interp_apply_tl(nlev, fvec_tl, f_tl, wi, wf)

integer,         intent(in)  :: nlev
real(kind_real), intent(in)  :: fvec_tl(nlev)
integer,         intent(in)  :: wi
real(kind_real), intent(in)  :: wf
real(kind_real), intent(out) :: f_tl

if (wi == missing_value(nlev)) then
  f_tl = missing_value(f_tl)
elseif (fvec_tl(wi) == missing_value(f_tl) .or. fvec_tl(wi+1) == missing_value(f_tl)) then
  f_tl = missing_value(f_tl)
else
  f_tl = fvec_tl(wi)*wf + fvec_tl(wi+1)*(1.0_kind_real-wf)
endif

end subroutine vert_interp_apply_tl

! ------------------------------------------------------------------------------

subroutine vert_interp_apply_ad(nlev, fvec_ad, f_ad, wi, wf)

integer,         intent(in)    :: nlev
real(kind_real), intent(inout) :: fvec_ad(nlev)
integer,         intent(in)    :: wi
real(kind_real), intent(in)    :: wf
real(kind_real), intent(in)    :: f_ad
real(kind_real) :: missing

missing = missing_value(missing)

! Do not modify the adjoint if the weight index is missing.
! This occurs when the observed vertical coordinate is missing.
if (wi == missing_value(nlev)) return

if (fvec_ad(wi) == missing .or. f_ad == missing) then
  fvec_ad(wi  ) = 0.0_kind_real
else
  fvec_ad(wi  ) = fvec_ad(wi  ) + f_ad*wf
endif
if (fvec_ad(wi+1) == missing .or. f_ad == missing) then
  fvec_ad(wi+1) = 0.0_kind_real
else
  fvec_ad(wi+1) = fvec_ad(wi+1) + f_ad*(1.0_kind_real-wf)
endif

end subroutine vert_interp_apply_ad

! ------------------------------------------------------------------------------

subroutine nearestneighbor_interp_index(nlev,obl,vec,idx)

integer,         intent(in ) :: nlev       !Number of model levels
real(kind_real), intent(in ) :: obl        !Observation location
real(kind_real), intent(in ) :: vec(nlev)  !Structured vector of grid points
integer,         intent(out) :: idx        !Index for interpolation

integer         :: k
real(kind_real) :: missing

missing = missing_value(obl)

! If the observation is missing then set both the index and weight to missing.
if (obl == missing) then
   idx= missing_value(nlev)
   return
end if

if (vec(1) < vec(nlev)) then !increases with index

  if (obl <= vec(1)) then
     idx = 1
  else
     idx = nlev
     do k = 2, nlev
        if (obl <= vec(k)) then
           if((obl - vec(k-1)) <= (vec(k)-obl)) then
             idx = k - 1
           else
             idx = k
           endif
           exit
        endif
     enddo
  endif

else !decreases with index

  if (obl >= vec(1)) then
     idx = 1
  else
     idx = nlev
     do k = 2, nlev
        if (obl >= vec(k)) then
           if((vec(k-1)-obl) <= (obl-vec(k))) then
             idx = k - 1
           else
             idx = k
           endif
           exit
        endif
     enddo
  endif

endif

end subroutine nearestneighbor_interp_index

! ------------------------------------------------------------------------------

subroutine nearestneighbor_interp_apply(nlev, fvec, f, idx)

integer,         intent(in ) :: nlev        !Number of model levels
real(kind_real), intent(in ) :: fvec(nlev)  !Field at grid points
integer,         intent(in ) :: idx         !Index for interpolation
real(kind_real), intent(out) :: f           !Output at obs location using linear interp

if (idx == missing_value(nlev)) then
  f = missing_value(f)
else
  f = fvec(idx)
endif

end subroutine nearestneighbor_interp_apply

! ------------------------------------------------------------------------------

subroutine nearestneighbor_interp_apply_tl(nlev, fvec_tl, f_tl, idx)

integer,         intent(in)  :: nlev
real(kind_real), intent(in)  :: fvec_tl(nlev)
integer,         intent(in)  :: idx
real(kind_real), intent(out) :: f_tl

if (idx== missing_value(nlev)) then
  f_tl = missing_value(f_tl)
else if (fvec_tl(idx) == missing_value(f_tl)) then
  f_tl = missing_value(f_tl)
else
  f_tl = fvec_tl(idx)
endif

end subroutine nearestneighbor_interp_apply_tl

! ------------------------------------------------------------------------------

subroutine nearestneighbor_interp_apply_ad(nlev, fvec_ad, f_ad, idx)

integer,         intent(in)    :: nlev
real(kind_real), intent(inout) :: fvec_ad(nlev)
integer,         intent(in)    :: idx
real(kind_real), intent(in)    :: f_ad
real(kind_real) :: missing

missing = missing_value(missing)

! Do not modify the adjoint if the weight index is missing.
! This occurs when the observed vertical coordinate is missing.
if (idx == missing_value(nlev)) return

if (fvec_ad(idx) == missing .or. f_ad == missing) then
  fvec_ad(idx) = 0.0_kind_real
else
  fvec_ad(idx) = fvec_ad(idx) + f_ad
endif

end subroutine nearestneighbor_interp_apply_ad

! ------------------------------------------------------------------------------

subroutine check_adjustment_function(adjustment_function, observation_vertical_coordinate, &
                                     geovars)

  ! Arguments
  character(len=*),               intent(in)    :: adjustment_function
  character(len=*),               intent(in)    :: observation_vertical_coordinate
  type(oops_variables), optional, intent(inout) :: geovars

  ! Locals
  integer :: i
  character(255), dimension(1) :: valid_adjustment_functions

  ! Populate list
  valid_adjustment_functions(1) = 'subtract scaled station elevation'

  ! Check if the input string is in the list
  do i = 1, size(valid_adjustment_functions)

    if (trim(adjustment_function) == trim(valid_adjustment_functions(i))) then
      ! If the adjustment function is subtract scaled station elevation check that the
      ! observation vertical coordinate is station elevation
      if (trim(adjustment_function) == "subtract scaled station elevation") then
        ! Add surface_geometric_height to the list of geovars
        if (present(geovars)) call geovars%push_back('surface_geometric_height')
        ! Check that height is the observation vertical coordinate
        if (trim(observation_vertical_coordinate) .ne. "height") &
          call abor1_ftn('Observation vertical coordinate adjustment using ' // &
                         'subtract scaled station elevation is only supported for height ' // &
                         'coordinates.')
      endif
      return
    end if
  end do

  call abor1_ftn("Invalid adjustment function: " // trim(adjustment_function))

end subroutine check_adjustment_function

! ------------------------------------------------------------------------------

subroutine adjust_obs_coordinate_subtract_scaled_station_elevation(nlocs, geovals, obss, obsvcoord)

use obsspace_mod
use ufo_geovals_mod

! Arguments
integer,            intent(in)    :: nlocs
type(ufo_geovals),  intent(in)    :: geovals
type(c_ptr), value, intent(in)    :: obss
real(kind_real),    intent(inout) :: obsvcoord(nlocs)

! Locals
integer :: iobs
real(kind_real) :: missing, factz, height_above_station
real(kind_real), allocatable :: stationElevation(:)
type(ufo_geoval), pointer :: sgh

! Missing values
missing = missing_value(obsvcoord(1))

! Allocate and get the station elevation
allocate(stationElevation(nlocs))
call obsspace_get_db(obss, 'MetaData', 'stationElevation', stationElevation)

! This function needs the surface geometric height from the geovals
call ufo_geovals_get_var(geovals, "surface_geometric_height", sgh)

! Subtract off combination of surface station elevation and station elevation depending on how close
! to surface
do iobs = 1, nlocs

    ! If the observation coordinate is missing, skip it
    if (obsvcoord(iobs) == missing) continue

    ! Height above the station
    height_above_station = obsvcoord(iobs) - stationElevation(iobs)

    ! Set a scaling factor for the difference between the surface height and the station elevation]
    ! F is a function that approximately scales between 0 (at the station) and 1 (at 1000m above the
    ! station)
    if (height_above_station > 1000.0) then
      factz = 1.0
    elseif (height_above_station > 10.0) then
      ! This function seems strange but is consistent with the old code in GSI. Should
      ! be tested in the future to see whether the formula can be improved.
      factz = height_above_station/990.0  ! - 1.0/99.0  ! Should there be a constant -1/99 here?
    else
      factz = 0.0
    endif

    ! Subtract off the scaled station elevation
    ! When f is zero (near the station) subtract off the station elevation
    ! When f is one  (far above the station) subtract off the surface height
    obsvcoord(iobs) = obsvcoord(iobs) -  factz        * sgh%vals(1, iobs) + &
                                        (factz - 1.0) * stationElevation(iobs)

enddo

end subroutine adjust_obs_coordinate_subtract_scaled_station_elevation

! ------------------------------------------------------------------------------

end module vert_interp_mod
