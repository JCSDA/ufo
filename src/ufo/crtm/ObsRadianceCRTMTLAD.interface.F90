! (C) Copyright 2017-2018 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

!> Fortran module for Fortran-C++ interface functions for CRTM tl/ad obs operator

module ufo_radiancecrtm_tlad_mod_c

  use iso_c_binding
  use config_mod
  use ufo_radiancecrtm_tlad_mod
  use string_f_c_mod
  use ufo_geovals_mod
  use ufo_geovals_mod_c,   only: ufo_geovals_registry

  implicit none
  private

#define LISTED_TYPE ufo_radiancecrtm_tlad

  !> Linked list interface - defines registry_t type
#include "../linkedList_i.f"

  !> Global registry
  type(registry_t) :: ufo_radiancecrtm_tlad_registry

contains

  ! ------------------------------------------------------------------------------
  !> Linked list implementation
#include "../linkedList_c.f"

! ------------------------------------------------------------------------------

subroutine ufo_radiancecrtm_tlad_setup_c(c_key_self, c_conf, c_nchan, c_channels, csin, &
                                    c_str_size) bind(c,name='ufo_radiancecrtm_tlad_setup_f90')
implicit none
integer(c_int), intent(inout) :: c_key_self
type(c_ptr),    intent(in)    :: c_conf
integer(c_int), intent(in) :: c_nchan
integer(c_int), intent(in) :: c_channels(c_nchan)
integer(c_int), intent(in) :: c_str_size
character(kind=c_char,len=1),intent(inout) :: csin(c_str_size+1)

type(ufo_radiancecrtm_tlad), pointer :: self

call ufo_radiancecrtm_tlad_registry%setup(c_key_self, self)

call self%setup(c_conf, c_channels)

call f_c_string_vector(self%varin, csin)

end subroutine ufo_radiancecrtm_tlad_setup_c

! ------------------------------------------------------------------------------

subroutine ufo_radiancecrtm_tlad_delete_c(c_key_self) bind(c,name='ufo_radiancecrtm_tlad_delete_f90')
implicit none
integer(c_int), intent(inout) :: c_key_self

type(ufo_radiancecrtm_tlad), pointer :: self

call ufo_radiancecrtm_tlad_registry%get(c_key_self, self)
call self%delete()
call ufo_radiancecrtm_tlad_registry%remove(c_key_self)

end subroutine ufo_radiancecrtm_tlad_delete_c

! ------------------------------------------------------------------------------

subroutine ufo_radiancecrtm_tlad_settraj_c(c_key_self, c_key_geovals, c_obsspace) &
                                       bind(c,name='ufo_radiancecrtm_tlad_settraj_f90')

implicit none
integer(c_int),     intent(in) :: c_key_self
integer(c_int),     intent(in) :: c_key_geovals
type(c_ptr), value, intent(in) :: c_obsspace

type(ufo_radiancecrtm_tlad), pointer :: self
type(ufo_geovals),       pointer :: geovals

character(len=*), parameter :: myname_="ufo_radiancecrtm_tlad_settraj_c"

call ufo_radiancecrtm_tlad_registry%get(c_key_self, self)
call ufo_geovals_registry%get(c_key_geovals,geovals)

call self%settraj(geovals, c_obsspace)

end subroutine ufo_radiancecrtm_tlad_settraj_c

! ------------------------------------------------------------------------------

subroutine ufo_radiancecrtm_simobs_tl_c(c_key_self, c_key_geovals, c_obsspace, c_nvars, c_nlocs, c_hofx) &
                                    bind(c,name='ufo_radiancecrtm_simobs_tl_f90')

implicit none
integer(c_int),     intent(in)    :: c_key_self
integer(c_int),     intent(in)    :: c_key_geovals
type(c_ptr), value, intent(in)    :: c_obsspace
integer(c_int),     intent(in)    :: c_nvars, c_nlocs
real(c_double),     intent(inout) :: c_hofx(c_nvars, c_nlocs)

type(ufo_radiancecrtm_tlad), pointer :: self
type(ufo_geovals),       pointer :: geovals

character(len=*), parameter :: myname_="ufo_radiancecrtm_simobs_tl_c"

call ufo_radiancecrtm_tlad_registry%get(c_key_self, self)
call ufo_geovals_registry%get(c_key_geovals,geovals)

call self%simobs_tl(geovals, c_obsspace, c_nvars, c_nlocs, c_hofx)

end subroutine ufo_radiancecrtm_simobs_tl_c

! ------------------------------------------------------------------------------

subroutine ufo_radiancecrtm_simobs_ad_c(c_key_self, c_key_geovals, c_obsspace, c_nvars, c_nlocs, c_hofx) &
                                    bind(c,name='ufo_radiancecrtm_simobs_ad_f90')

implicit none
integer(c_int),     intent(in) :: c_key_self
integer(c_int),     intent(in) :: c_key_geovals
type(c_ptr), value, intent(in) :: c_obsspace
integer(c_int),     intent(in) :: c_nvars, c_nlocs
real(c_double),     intent(in) :: c_hofx(c_nvars, c_nlocs)

type(ufo_radiancecrtm_tlad), pointer :: self
type(ufo_geovals),       pointer :: geovals

character(len=*), parameter :: myname_="ufo_radiancecrtm_simobs_ad_c"

call ufo_radiancecrtm_tlad_registry%get(c_key_self, self)
call ufo_geovals_registry%get(c_key_geovals,geovals)

call self%simobs_ad(geovals, c_obsspace, c_nvars, c_nlocs, c_hofx)

end subroutine ufo_radiancecrtm_simobs_ad_c

! ------------------------------------------------------------------------------

end module ufo_radiancecrtm_tlad_mod_c
