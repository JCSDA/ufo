! (C) Copyright 2021 Met Office UK
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

module ufo_backgrounderroridentity_mod

use iso_c_binding,      only: c_ptr
use oops_variables_mod, only: oops_variables
use ufo_geovals_mod,    only: ufo_geovals

contains

!> For each obs diagnostic called <var>_background_error, where <var> belongs to the set of variable
!> names @p obsvars, fill this diagnostic with estimates of the background error of variable <var>
!> at observation locations.
subroutine ufo_backgrounderroridentity_fillobsdiags(geovals, nlocs, obsvars, obsdiags)
  use kinds,           only: kind_real
  use ufo_geovals_mod, only: ufo_geoval, ufo_geovals, ufo_geovals_get_var
  use ufo_vars_mod,    only: MAXVARLEN
  implicit none

  type(ufo_geovals), intent(in)    :: geovals
  integer, intent(in)              :: nlocs
  type(oops_variables), intent(in) :: obsvars
  type(ufo_geovals), intent(inout) :: obsdiags

  type(ufo_geoval), pointer        :: background_error
  integer                          :: ivar
  character(len=MAXVARLEN)         :: varstr
  integer                          :: lenvarstr

  character(len=*), parameter      :: suffix = "_background_error"

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

    ! Get the background error geoval.
    call ufo_geovals_get_var(geovals, varstr, background_error)

    ! Allocate the background error diagnostic.
    if (allocated(obsdiags%geovals(ivar)%vals)) deallocate(obsdiags%geovals(ivar)%vals)
    obsdiags%geovals(ivar)%nval = 1
    allocate(obsdiags%geovals(ivar)%vals(obsdiags%geovals(ivar)%nval, nlocs))

    ! Copy the geoval to the diagnostic.
    obsdiags%geovals(ivar)%vals(1, 1:nlocs) = background_error%vals(1, 1:nlocs)
  enddo
end subroutine ufo_backgrounderroridentity_fillobsdiags

end module ufo_backgrounderroridentity_mod
