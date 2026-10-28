! (C) Copyright 2021 Met Office
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

!> Fortran module to handle ground based GNSS observations operator at the Met Office

module ufo_groundgnss_metoffice_mod_c

  use fckit_configuration_module, only: fckit_configuration
  use, intrinsic :: iso_c_binding
  use ufo_groundgnss_metoffice_mod

  implicit none
  private

#define LISTED_TYPE ufo_groundgnss_MetOffice

  !> Linked list interface - defines registry_t type
#include "oops/util/linkedList_i.f"

  !> Global registry
  type(registry_t) :: ufo_groundgnss_MetOffice_registry

  ! ------------------------------------------------------------------------------
contains
  ! ------------------------------------------------------------------------------
  !> Linked list implementation
#include "oops/util/linkedList_c.f"

! ------------------------------------------------------------------------------

subroutine ufo_groundgnss_metoffice_setup_c(c_key_self, &
                                            c_vertInterpOPS, &
                                            c_pseudoLevels, &
                                            c_minTempGrad, &
                                            c_dryRefractivityConstant, &
                                            c_wetRefractivityConstant, &
                                            c_numFOV, &
                                            c_fovRadius) bind(c,name="ufo_groundgnss_metoffice_setup_f90")
implicit none
integer(c_int), intent(inout)  :: c_key_self
logical(c_bool), intent(in) :: c_vertInterpOPS
logical(c_bool), intent(in) :: c_pseudoLevels
real(c_float), intent(in) :: c_minTempGrad
real(c_float), intent(in) :: c_dryRefractivityConstant
real(c_float), intent(in) :: c_wetRefractivityConstant
integer(c_int), intent(in) :: c_numFOV
real(c_float), intent(in) :: c_fovRadius

type(ufo_groundgnss_MetOffice), pointer :: self

call ufo_groundgnss_MetOffice_registry%setup(c_key_self, self)

call self%setup(c_vertInterpOPS, c_pseudoLevels, c_minTempGrad, &
                c_dryRefractivityConstant, c_wetRefractivityConstant, &
                c_numFOV, c_fovRadius)

end subroutine ufo_groundgnss_MetOffice_setup_c

! ------------------------------------------------------------------------------

subroutine ufo_groundgnss_metoffice_delete_c(c_key_self) bind(c,name="ufo_groundgnss_metoffice_delete_f90")
implicit none
integer(c_int), intent(inout) :: c_key_self

call ufo_groundgnss_MetOffice_registry%remove(c_key_self)

end subroutine ufo_groundgnss_metoffice_delete_c

! ------------------------------------------------------------------------------

subroutine ufo_groundgnss_metoffice_simobs_c(c_key_self, c_key_geovals, c_obsspace, &
                                             c_nobs, c_hofx) bind(c,name="ufo_groundgnss_metoffice_simobs_f90")

implicit none
integer(c_int), intent(in) :: c_key_self
integer(c_int), intent(in) :: c_key_geovals
type(c_ptr), value, intent(in) :: c_obsspace
integer(c_int), intent(in) :: c_nobs
real(c_double), intent(inout) :: c_hofx(c_nobs)

type(ufo_groundgnss_MetOffice),     pointer :: self

character(len=*), parameter :: myname_="ufo_groundgnss_metoffice_simobs_c"
call ufo_groundgnss_MetOffice_registry%get(c_key_self, self)
call self%opr_simobs(c_key_geovals, c_obsspace, c_hofx)

end subroutine ufo_groundgnss_metoffice_simobs_c

! ------------------------------------------------------------------------------

subroutine ufo_groundgnss_metoffice_locs_init_c(c_key_self, c_obsspace, c_npaths, &
                                                c_lons, c_lats) &
    bind(c,name="ufo_groundgnss_metoffice_locs_init_f90")
implicit none
integer(c_int),     intent(in)    :: c_key_self
type(c_ptr), value, intent(in)    :: c_obsspace
integer(c_int),     intent(in)    :: c_npaths
real(c_float),      intent(inout) :: c_lons(c_npaths)
real(c_float),      intent(inout) :: c_lats(c_npaths)

type(ufo_groundgnss_MetOffice), pointer :: self

call ufo_groundgnss_MetOffice_registry%get(c_key_self, self)
call ufo_groundgnss_metoffice_locs_init(self, c_obsspace, c_npaths, c_lons, c_lats)

end subroutine ufo_groundgnss_metoffice_locs_init_c

! ------------------------------------------------------------------------------

end module ufo_groundgnss_metoffice_mod_c
