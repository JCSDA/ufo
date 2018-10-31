! (C) Copyright 2017 UCAR
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 

!> Fortran module to handle ice concentration observations

module ufo_seasurfacetemp_tlad_mod_c
  
  use iso_c_binding
  use config_mod
  use ufo_geovals_mod,   only: ufo_geovals
  use ufo_geovals_mod_c, only: ufo_geovals_registry
  use ufo_seasurfacetemp_tlad_mod 
  implicit none
  private
  
#define LISTED_TYPE ufo_seasurfacetemp_tlad
  
  !> Linked list interface - defines registry_t type
#include "../../linkedList_i.f"
  
  !> Global registry
  type(registry_t) :: ufo_seasurfacetemp_tlad_registry
  
  ! ------------------------------------------------------------------------------
contains
  ! ------------------------------------------------------------------------------
  !> Linked list implementation
#include "../../linkedList_c.f"
  
! ------------------------------------------------------------------------------
  
subroutine ufo_seasurfacetemp_tlad_setup_c(c_key_self, c_conf) bind(c,name='ufo_seasurfacetemp_tlad_setup_f90')
implicit none
integer(c_int), intent(inout) :: c_key_self
type(c_ptr), intent(in)    :: c_conf
    
type(ufo_seasurfacetemp_tlad), pointer :: self

call ufo_seasurfacetemp_tlad_registry%init()
call ufo_seasurfacetemp_tlad_registry%add(c_key_self)
call ufo_seasurfacetemp_tlad_registry%get(c_key_self, self)
    
end subroutine ufo_seasurfacetemp_tlad_setup_c
  
! ------------------------------------------------------------------------------
  
subroutine ufo_seasurfacetemp_tlad_delete_c(c_key_self) bind(c,name='ufo_seasurfacetemp_tlad_delete_f90')
implicit none
integer(c_int), intent(inout) :: c_key_self
    
type(ufo_seasurfacetemp_tlad), pointer :: self

call ufo_seasurfacetemp_tlad_registry%get(c_key_self, self)
call ufo_seasurfacetemp_tlad_registry%remove(c_key_self)
    
end subroutine ufo_seasurfacetemp_tlad_delete_c
  
! ------------------------------------------------------------------------------

subroutine ufo_seasurfacetemp_tlad_settraj_c(c_key_self, c_key_geovals) bind(c,name='ufo_seasurfacetemp_tlad_settraj_f90')

implicit none
integer(c_int), intent(in) :: c_key_self
integer(c_int), intent(in) :: c_key_geovals

type(ufo_seasurfacetemp_tlad), pointer :: self
type(ufo_geovals),    pointer :: geovals

character(len=*), parameter :: myname_="ufo_seasurfacetemp_tlad_settraj_c"

call ufo_seasurfacetemp_tlad_registry%get(c_key_self, self)
call ufo_geovals_registry%get(c_key_geovals,geovals)

call ufo_seasurfacetemp_tlad_settraj(self, geovals)

end subroutine ufo_seasurfacetemp_tlad_settraj_c

! ------------------------------------------------------------------------------

subroutine ufo_seasurfacetemp_simobs_tl_c(c_key_self, c_key_geovals, c_nobs, c_hofx) bind(c,name='ufo_seasurfacetemp_simobs_tl_f90')

implicit none
integer(c_int), intent(in) :: c_key_self
integer(c_int), intent(in) :: c_key_geovals
integer(c_int),     intent(in) :: c_nobs
real(c_double),     intent(inout) :: c_hofx(c_nobs)

type(ufo_seasurfacetemp_tlad), pointer :: self
type(ufo_geovals),    pointer :: geovals

character(len=*), parameter :: myname_="ufo_seasurfacetemp_simobs_tl_c"

call ufo_seasurfacetemp_tlad_registry%get(c_key_self, self)
call ufo_geovals_registry%get(c_key_geovals,geovals)

call ufo_seasurfacetemp_simobs_tl(self, geovals, c_hofx)

end subroutine ufo_seasurfacetemp_simobs_tl_c

! ------------------------------------------------------------------------------

subroutine ufo_seasurfacetemp_simobs_ad_c(c_key_self, c_key_geovals, c_nobs, c_hofx) bind(c,name='ufo_seasurfacetemp_simobs_ad_f90')

implicit none
integer(c_int), intent(in) :: c_key_self
integer(c_int), intent(in) :: c_key_geovals
integer(c_int),     intent(in) :: c_nobs
real(c_double),     intent(inout) :: c_hofx(c_nobs)

type(ufo_seasurfacetemp_tlad), pointer :: self
type(ufo_geovals),    pointer :: geovals

character(len=*), parameter :: myname_="ufo_seasurfacetemp_simobs_ad_c"

call ufo_seasurfacetemp_tlad_registry%get(c_key_self, self)
call ufo_geovals_registry%get(c_key_geovals,geovals)

call ufo_seasurfacetemp_simobs_ad(self, geovals, c_hofx)

end subroutine ufo_seasurfacetemp_simobs_ad_c

! ------------------------------------------------------------------------------

end module ufo_seasurfacetemp_tlad_mod_c
