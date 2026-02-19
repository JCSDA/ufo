! (C) British Crown Copyright 2021 Met Office
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 

!> Fortran module to handle gnssro observations - refractivity Met Office 1d operator

module ufo_gnssro_refmetoffice_mod_c
  
  use logger_mod, only: oops_log
  use iso_c_binding
  use ufo_gnssro_refmetoffice_mod
  use ufo_geovals_mod
  use ufo_geovals_mod_c,   only: ufo_geovals_registry

  implicit none
  private
  
#define LISTED_TYPE ufo_gnssro_RefMetOffice
  
  !> Linked list interface - defines registry_t type
#include "oops/util/linkedList_i.f"
  
  !> Global registry
  type(registry_t) :: ufo_gnssro_RefMetOffice_registry
  
  ! ------------------------------------------------------------------------------
contains
  ! ------------------------------------------------------------------------------
  !> Linked list implementation
#include "oops/util/linkedList_c.f"
  
! ------------------------------------------------------------------------------
  
subroutine ufo_gnssro_refmetoffice_setup_c(c_key_self, &
                                           vert_interp_ops, &
                                           pseudo_ops, &
                                           min_temp_grad, &
                                           dryRefractivityConstant, &
                                           wetRefractivityConstant) bind(c,name='ufo_gnssro_refmetoffice_setup_f90')
implicit none
integer(c_int), intent(inout) :: c_key_self
logical(c_bool), intent(in) :: vert_interp_ops
logical(c_bool), intent(in) :: pseudo_ops
real(c_float), intent(in) :: min_temp_grad
real(c_float), intent(in) :: dryRefractivityConstant
real(c_float), intent(in) :: wetRefractivityConstant

type(ufo_gnssro_RefMetOffice), pointer :: self

call ufo_gnssro_refmetoffice_registry%setup(c_key_self, self)

call self%setup(vert_interp_ops, pseudo_ops, min_temp_grad, dryRefractivityConstant, &
                wetRefractivityConstant)

end subroutine ufo_gnssro_refmetoffice_setup_c
  
! ------------------------------------------------------------------------------
  
subroutine ufo_gnssro_refmetoffice_delete_c(c_key_self) bind(c,name='ufo_gnssro_refmetoffice_delete_f90')
implicit none
integer(c_int), intent(inout) :: c_key_self

call ufo_gnssro_RefMetOffice_registry%remove(c_key_self)

end subroutine ufo_gnssro_refmetoffice_delete_c
  
! ------------------------------------------------------------------------------

subroutine ufo_gnssro_refmetoffice_simobs_c(c_key_self, c_key_geovals, c_obsspace, &
                                            c_nobs, c_hofx, c_key_obs_diags) &
                                            bind(c,name='ufo_gnssro_refmetoffice_simobs_f90')

implicit none
integer(c_int), intent(in)     :: c_key_self       ! Key giving pointer to self object
integer(c_int), intent(in)     :: c_key_geovals    ! Key giving pointer to geovals object
type(c_ptr), value, intent(in) :: c_obsspace       ! Pointer to obs-space object
integer(c_int), intent(in)     :: c_nobs           ! Number of observations
real(c_double), intent(inout)  :: c_hofx(c_nobs)   ! Array of calculated H(x) object
integer(c_int), intent(in)     :: c_key_obs_diags  ! Key giving pointer to obs diagnostics object

type(ufo_gnssro_RefMetOffice), pointer  :: self            ! Self object
type(ufo_geovals),  pointer             :: obs_diags       ! Observations diagnostics
type(ufo_geovals), pointer              :: geovals         ! Geovals object
character(len=*), parameter             :: myname_="ufo_gnssro_refmetoffice_simobs_c"
character(len=200)                      :: output_message  ! Message to be output

write(output_message, *) 'Beginning interface', c_key_obs_diags, c_key_geovals, c_key_self
call oops_log % trace(output_message)

call ufo_gnssro_RefMetOffice_registry % get(c_key_self, self)
call ufo_geovals_registry % get(c_key_obs_diags, obs_diags)
call ufo_geovals_registry % get(c_key_geovals, geovals)

call self%simobs(geovals, c_obsspace, c_hofx, obs_diags)

write(output_message, *) 'Finishing interface'
call oops_log % trace(output_message)

end subroutine ufo_gnssro_refmetoffice_simobs_c

! ------------------------------------------------------------------------------

end module ufo_gnssro_refmetoffice_mod_c
