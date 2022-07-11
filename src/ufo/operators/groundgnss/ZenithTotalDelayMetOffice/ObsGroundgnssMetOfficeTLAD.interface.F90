! (C) Copyright 2021 Met Office
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 

!> Fortran module to handle ground gnss observations

module ufo_groundgnss_metoffice_tlad_mod_c

  use fckit_configuration_module, only: fckit_configuration
  use ufo_groundgnss_metoffice_tlad_mod
  implicit none
  private
  
#define LISTED_TYPE ufo_groundgnss_metoffice_tlad

  !> Linked list interface - defines registry_t type
#include "oops/util/linkedList_i.f"

  !> Global registry
  type(registry_t) :: ufo_groundgnss_metoffice_tlad_registry

  ! ------------------------------------------------------------------------------
contains
  ! ------------------------------------------------------------------------------
  !> Linked list implementation
#include "oops/util/linkedList_c.f"
  
! ------------------------------------------------------------------------------
  
subroutine ufo_groundgnss_metoffice_tlad_setup_c(c_key_self, c_conf) bind(c,name='ufo_groundgnss_metoffice_tlad_setup_f90')
implicit none
integer(c_int), intent(inout)  :: c_key_self
type(c_ptr), value, intent(in) :: c_conf
    
type(ufo_groundgnss_metoffice_tlad), pointer :: self
type(fckit_configuration)                    :: f_conf

call ufo_groundgnss_metoffice_tlad_registry%setup(c_key_self, self)
f_conf = fckit_configuration(c_conf)

call self%setup(f_conf)

end subroutine ufo_groundgnss_metoffice_tlad_setup_c
  
! ------------------------------------------------------------------------------
  
subroutine ufo_groundgnss_metoffice_tlad_delete_c(c_key_self) bind(c,name='ufo_groundgnss_metoffice_tlad_delete_f90')
implicit none
integer(c_int), intent(inout) :: c_key_self
    
type(ufo_groundgnss_metoffice_tlad), pointer :: self

call ufo_groundgnss_metoffice_tlad_registry%get(c_key_self, self)
call self%opr_delete()
call ufo_groundgnss_metoffice_tlad_registry%remove(c_key_self)
    
end subroutine ufo_groundgnss_metoffice_tlad_delete_c
  
! ------------------------------------------------------------------------------

subroutine ufo_groundgnss_metoffice_tlad_settraj_c(c_key_self, c_key_geovals, c_obsspace) bind(c,name='ufo_groundgnss_metoffice_tlad_settraj_f90')

implicit none
integer(c_int),     intent(in) :: c_key_self
integer(c_int),     intent(in) :: c_key_geovals
type(c_ptr), value, intent(in) :: c_obsspace

type(ufo_groundgnss_metoffice_tlad), pointer :: self

character(len=*), parameter :: myname_="ufo_groundgnss_metoffice_tlad_settraj_c"

call ufo_groundgnss_metoffice_tlad_registry%get(c_key_self, self)
call self%opr_settraj(c_key_geovals, c_obsspace)

end subroutine ufo_groundgnss_metoffice_tlad_settraj_c

! ------------------------------------------------------------------------------

subroutine ufo_groundgnss_metoffice_simobs_tl_c(c_key_self, c_key_geovals, c_obsspace, c_nobs, c_hofx) bind(c,name='ufo_groundgnss_metoffice_simobs_tl_f90')

implicit none
integer(c_int), intent(in) :: c_key_self
integer(c_int), intent(in) :: c_key_geovals
type(c_ptr), value, intent(in) :: c_obsspace
integer(c_int), intent(in) :: c_nobs
real(c_double), intent(inout) :: c_hofx(c_nobs)

type(ufo_groundgnss_metoffice_tlad), pointer :: self

character(len=*), parameter :: myname_="ufo_groundgnss_metoffice_simobs_tl_c"

call ufo_groundgnss_metoffice_tlad_registry%get(c_key_self, self)
call self%opr_simobs_tl(c_key_geovals, c_obsspace, c_hofx)

end subroutine ufo_groundgnss_metoffice_simobs_tl_c

! ------------------------------------------------------------------------------

subroutine ufo_groundgnss_metoffice_simobs_ad_c(c_key_self, c_key_geovals, c_obsspace, c_nobs, c_hofx) bind(c,name='ufo_groundgnss_metoffice_simobs_ad_f90')

implicit none
integer(c_int), intent(in) :: c_key_self
integer(c_int), intent(in) :: c_key_geovals
type(c_ptr), value, intent(in) :: c_obsspace
integer(c_int), intent(in) :: c_nobs
real(c_double), intent(in) :: c_hofx(c_nobs)

type(ufo_groundgnss_metoffice_tlad), pointer :: self

character(len=*), parameter :: myname_="ufo_groundgnss_metoffice_simobs_ad_c"

call ufo_groundgnss_metoffice_tlad_registry%get(c_key_self, self)
call self%opr_simobs_ad(c_key_geovals, c_obsspace, c_hofx)

end subroutine ufo_groundgnss_metoffice_simobs_ad_c

! ------------------------------------------------------------------------------
  
end module ufo_groundgnss_metoffice_tlad_mod_c
