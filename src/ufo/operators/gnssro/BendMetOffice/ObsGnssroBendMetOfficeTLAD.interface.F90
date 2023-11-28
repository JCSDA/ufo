! (C) Copyright 2017 UCAR
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 

!> Fortran module to handle gnssro observations

module ufo_gnssro_bendmetoffice_tlad_mod_c
  
  use iso_c_binding, only: c_int, c_bool, c_float, c_double, c_ptr
  use ufo_gnssro_bendmetoffice_tlad_mod
  implicit none
  private
  
#define LISTED_TYPE ufo_gnssro_bendmetoffice_tlad

  !> Linked list interface - defines registry_t type
#include "oops/util/linkedList_i.f"

  !> Global registry
  type(registry_t) :: ufo_gnssro_bendmetoffice_tlad_registry

  ! ------------------------------------------------------------------------------
contains
  ! ------------------------------------------------------------------------------
  !> Linked list implementation
#include "oops/util/linkedList_c.f"
  
! ------------------------------------------------------------------------------

subroutine ufo_gnssro_bendmetoffice_tlad_setup_c(c_key_self, &
                                                 vert_interp_ops, &
                                                 pseudo_ops, &
                                                 min_temp_grad, &
                                                 nchans, &
                                                 chanList, &
                                                 noSuperCheck) bind(c,name='ufo_gnssro_bendmetoffice_tlad_setup_f90')
implicit none
integer(c_int), intent(inout) :: c_key_self        !< Reference to this object
logical(c_bool), intent(in)   :: vert_interp_ops   !< Whether to do vertical interpolation using ln(p)
logical(c_bool), intent(in)   :: pseudo_ops        !< Whether to use pseudo-levels
real(c_float), intent(in)     :: min_temp_grad     !< Minimum temperature gradient
integer(c_int), intent(in)    :: nchans            !< Number of channels (levels) to be used
integer(c_int), intent(in)    :: chanList(nchans)  !< List of channels to use
logical(c_bool), intent(in)   :: noSuperCheck      !< Whether to avoid using super-refraction check in operator

integer(c_int)                :: noChans(1)        !< Channel list when no channels are used

type(ufo_gnssro_bendmetoffice_tlad), pointer :: self

call ufo_gnssro_bendmetoffice_tlad_registry%setup(c_key_self, self)

if (nchans == 0) then
  noChans(1) = 0
  call self%setup(vert_interp_ops, pseudo_ops, min_temp_grad, noChans, noSuperCheck)
else
  call self%setup(vert_interp_ops, pseudo_ops, min_temp_grad, chanList, noSuperCheck)
end if

end subroutine ufo_gnssro_bendmetoffice_tlad_setup_c
  
! ------------------------------------------------------------------------------
  
subroutine ufo_gnssro_bendmetoffice_tlad_delete_c(c_key_self) bind(c,name='ufo_gnssro_bendmetoffice_tlad_delete_f90')
implicit none
integer(c_int), intent(inout) :: c_key_self
    
type(ufo_gnssro_bendmetoffice_tlad), pointer :: self

call ufo_gnssro_bendmetoffice_tlad_registry%get(c_key_self, self)
call self%opr_delete()
call ufo_gnssro_bendmetoffice_tlad_registry%remove(c_key_self)
    
end subroutine ufo_gnssro_bendmetoffice_tlad_delete_c
  
! ------------------------------------------------------------------------------

subroutine ufo_gnssro_bendmetoffice_tlad_settraj_c(c_key_self, c_key_geovals, c_obsspace) bind(c,name='ufo_gnssro_bendmetoffice_tlad_settraj_f90')

implicit none
integer(c_int),     intent(in) :: c_key_self
integer(c_int),     intent(in) :: c_key_geovals
type(c_ptr), value, intent(in) :: c_obsspace

type(ufo_gnssro_bendmetoffice_tlad), pointer :: self

character(len=*), parameter :: myname_="ufo_gnssro_bendmetoffice_tlad_settraj_c"

call ufo_gnssro_bendmetoffice_tlad_registry%get(c_key_self, self)
call self%opr_settraj(c_key_geovals, c_obsspace)

end subroutine ufo_gnssro_bendmetoffice_tlad_settraj_c

! ------------------------------------------------------------------------------

subroutine ufo_gnssro_bendmetoffice_simobs_tl_c(c_key_self, c_key_geovals, c_obsspace, c_nobs, c_hofx) bind(c,name='ufo_gnssro_bendmetoffice_simobs_tl_f90')

implicit none
integer(c_int), intent(in) :: c_key_self
integer(c_int), intent(in) :: c_key_geovals
type(c_ptr), value, intent(in) :: c_obsspace
integer(c_int), intent(in) :: c_nobs
real(c_double), intent(inout) :: c_hofx(c_nobs)

type(ufo_gnssro_bendmetoffice_tlad), pointer :: self

character(len=*), parameter :: myname_="ufo_gnssro_bendmetoffice_simobs_tl_c"

call ufo_gnssro_bendmetoffice_tlad_registry%get(c_key_self, self)
call self%opr_simobs_tl(c_key_geovals, c_obsspace, c_hofx)

end subroutine ufo_gnssro_bendmetoffice_simobs_tl_c

! ------------------------------------------------------------------------------

subroutine ufo_gnssro_bendmetoffice_simobs_ad_c(c_key_self, c_key_geovals, c_obsspace, c_nobs, c_hofx) bind(c,name='ufo_gnssro_bendmetoffice_simobs_ad_f90')

implicit none
integer(c_int), intent(in) :: c_key_self
integer(c_int), intent(in) :: c_key_geovals
type(c_ptr), value, intent(in) :: c_obsspace
integer(c_int), intent(in) :: c_nobs
real(c_double), intent(in) :: c_hofx(c_nobs)

type(ufo_gnssro_bendmetoffice_tlad), pointer :: self

character(len=*), parameter :: myname_="ufo_gnssro_bendmetoffice_simobs_ad_c"

call ufo_gnssro_bendmetoffice_tlad_registry%get(c_key_self, self)
call self%opr_simobs_ad(c_key_geovals, c_obsspace, c_hofx)

end subroutine ufo_gnssro_bendmetoffice_simobs_ad_c

! ------------------------------------------------------------------------------
  
end module ufo_gnssro_bendmetoffice_tlad_mod_c
