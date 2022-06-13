! (C) Copyright 2021 Met Office UK
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

module ufo_backgrounderroridentity_mod_c

use iso_c_binding
use ufo_backgrounderroridentity_mod, only: ufo_backgrounderroridentity_fillobsdiags
implicit none

contains

! ------------------------------------------------------------------------------

subroutine ufo_backgrounderroridentity_fillobsdiags_c(c_key_geovals, c_nlocs, &
                                                      c_obsvars, c_key_obsdiags) &
  bind(c, name='ufo_backgrounderroridentity_fillobsdiags_f90')

  use oops_variables_mod, only: oops_variables
  use ufo_geovals_mod,    only: ufo_geovals
  use ufo_geovals_mod_c,  only: ufo_geovals_registry
  implicit none

  integer(c_int), intent(in)     :: c_key_geovals
  integer(c_int), intent(in)     :: c_nlocs
  type(c_ptr), value, intent(in) :: c_obsvars
  integer(c_int), intent(in)     :: c_key_obsdiags

  type(ufo_geovals), pointer :: geovals
  type(oops_variables)       :: obsvars
  type(ufo_geovals), pointer :: obsdiags

  call ufo_geovals_registry%get(c_key_geovals, geovals)
  obsvars = oops_variables(c_obsvars)
  call ufo_geovals_registry%get(c_key_obsdiags, obsdiags)

  call ufo_backgrounderroridentity_fillobsdiags(geovals, c_nlocs, obsvars, obsdiags)

end subroutine ufo_backgrounderroridentity_fillobsdiags_c

end module ufo_backgrounderroridentity_mod_c
