! (C) Copyright 2021 Met Office UK
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

module ufo_backgrounderrorvertinterp_mod_c

use iso_c_binding
use ufo_backgrounderrorvertinterp_mod, only: ufo_backgrounderrorvertinterp_fillobsdiags
implicit none

contains

! ------------------------------------------------------------------------------

subroutine ufo_backgrounderrorvertinterp_fillobsdiags_c(len_obs_vcoord, c_obs_vcoord, &
                                                        len_obs_vgroup, c_obs_vgroup, &
                                                        len_vcoord, c_vcoord, &
                                                        c_key_geovals, c_obsspace, c_nlocs, &
                                                        c_obsvars, c_key_obsdiags) &
  bind(c, name='ufo_backgrounderrorvertinterp_fillobsdiags_f90')

  use string_f_c_mod,     only: c_f_string
  use oops_variables_mod, only: oops_variables
  use ufo_geovals_mod,    only: ufo_geovals
  use ufo_geovals_mod_c,  only: ufo_geovals_registry
  use ufo_vars_mod,       only: MAXVARLEN
  implicit none

  integer(c_int), intent(in) :: len_obs_vcoord
  character(kind=c_char, len=1), intent(in) :: c_obs_vcoord(len_obs_vcoord + 1)
  integer(c_int), intent(in) :: len_obs_vgroup
  character(kind=c_char, len=1), intent(in) :: c_obs_vgroup(len_obs_vgroup + 1)
  integer(c_int), intent(in) :: len_vcoord
  character(kind=c_char, len=1), intent(in) :: c_vcoord(len_vcoord + 1)
  integer(c_int), intent(in) :: c_key_geovals
  type(c_ptr), value, intent(in) :: c_obsspace
  integer(c_int), intent(in) :: c_nlocs
  type(c_ptr), value, intent(in) :: c_obsvars
  integer(c_int), intent(in) :: c_key_obsdiags

  character(len=MAXVARLEN) :: obs_vcoord
  character(len=MAXVARLEN) :: obs_vgroup
  character(len=MAXVARLEN) :: vcoord
  type(ufo_geovals), pointer :: geovals
  type(oops_variables)       :: obsvars
  type(ufo_geovals), pointer :: obsdiags

  call c_f_string(c_obs_vcoord, obs_vcoord)
  call c_f_string(c_obs_vgroup, obs_vgroup)
  call c_f_string(c_vcoord, vcoord)
  call ufo_geovals_registry%get(c_key_geovals, geovals)
  obsvars = oops_variables(c_obsvars)
  call ufo_geovals_registry%get(c_key_obsdiags, obsdiags)

  call ufo_backgrounderrorvertinterp_fillobsdiags(obs_vcoord, obs_vgroup, vcoord, &
                                                  geovals, c_obsspace, c_nlocs, obsvars, obsdiags)

end subroutine ufo_backgrounderrorvertinterp_fillobsdiags_c

end module ufo_backgrounderrorvertinterp_mod_c
