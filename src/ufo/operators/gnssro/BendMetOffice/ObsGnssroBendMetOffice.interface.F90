! (C) Copyright 2017 UCAR
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 

!> Fortran module to handle gnssro observations-bending angle Met Office 1d operator

module ufo_gnssro_bendmetoffice_mod_c
  
  use fckit_log_module,  only : fckit_log
  use iso_c_binding
  use ufo_gnssro_bendmetoffice_mod
  use ufo_geovals_mod
  use ufo_geovals_mod_c,   only: ufo_geovals_registry

  implicit none
  private
  
#define LISTED_TYPE ufo_gnssro_BendMetOffice
  
  !> Linked list interface - defines registry_t type
#include "oops/util/linkedList_i.f"
  
  !> Global registry
  type(registry_t) :: ufo_gnssro_BendMetOffice_registry
  
  ! ------------------------------------------------------------------------------
contains
  ! ------------------------------------------------------------------------------
  !> Linked list implementation
#include "oops/util/linkedList_c.f"
  
! ------------------------------------------------------------------------------

subroutine ufo_gnssro_bendmetoffice_setup_c(c_key_self, &
                                            vert_interp_ops, &
                                            pseudo_ops, &
                                            min_temp_grad, &
                                            nchans, &
                                            chanList, &
                                            noSuperCheck) bind(c,name='ufo_gnssro_bendmetoffice_setup_f90')
implicit none
integer(c_int), intent(inout) :: c_key_self        !< Reference to this object
logical(c_bool), intent(in)   :: vert_interp_ops   !< Whether to do vertical interpolation using ln(p)
logical(c_bool), intent(in)   :: pseudo_ops        !< Whether to use pseudo-levels
real(c_float), intent(in)     :: min_temp_grad     !< Minimum temperature gradient
integer(c_int), intent(in)    :: nchans            !< Number of channels (levels) to be used
integer(c_int), intent(in)    :: chanList(nchans)  !< List of channels to use
logical(c_bool), intent(in)   :: noSuperCheck      !< Whether to avoid using super-refraction check in operator

integer(c_int)                :: noChans(1)        !< Channel list when no channels are used

type(ufo_gnssro_BendMetOffice), pointer :: self

call ufo_gnssro_bendmetoffice_registry%setup(c_key_self, self)

if (nchans == 0) then
  noChans(1) = 0
  call self%setup(vert_interp_ops, pseudo_ops, min_temp_grad, noChans, noSuperCheck)
else
  call self%setup(vert_interp_ops, pseudo_ops, min_temp_grad, chanList, noSuperCheck)
end if

end subroutine ufo_gnssro_bendmetoffice_setup_c
  
! ------------------------------------------------------------------------------
  
subroutine ufo_gnssro_bendmetoffice_delete_c(c_key_self) bind(c,name='ufo_gnssro_bendmetoffice_delete_f90')
implicit none
integer(c_int), intent(inout) :: c_key_self

call ufo_gnssro_BendMetOffice_registry%remove(c_key_self)

end subroutine ufo_gnssro_bendmetoffice_delete_c
  
! ------------------------------------------------------------------------------

subroutine ufo_gnssro_bendmetoffice_simobs_c(c_key_self, c_key_geovals, c_obsspace, &
                                             c_nvars, c_nlocs, c_hofx, c_key_obs_diags) &
                                             bind(c,name='ufo_gnssro_bendmetoffice_simobs_f90')

implicit none
integer(c_int), intent(in)     :: c_key_self                 ! Key giving pointer to self object
integer(c_int), intent(in)     :: c_key_geovals              ! Key giving pointer to geovals object
type(c_ptr), value, intent(in) :: c_obsspace                 ! Pointer to obs-space object
integer(c_int), intent(in)     :: c_nvars                    ! Number of variables (levels) in the data
integer(c_int), intent(in)     :: c_nlocs                    ! Number of observations (locations)
real(c_double), intent(inout)  :: c_hofx(c_nvars, c_nlocs)   ! Array of calculated H(x) object
integer(c_int), intent(in)     :: c_key_obs_diags            ! Key giving pointer to obs diagnostics object

type(ufo_gnssro_BendMetOffice), pointer :: self            ! Self object
type(ufo_geovals),  pointer             :: obs_diags       ! Observations diagnostics
type(ufo_geovals), pointer              :: geovals         ! Geovals object
character(len=*), parameter             :: myname_="ufo_gnssro_bendmetoffice_simobs_c"
character(len=200)                      :: output_message  ! Message to be output

write(output_message, *) 'TRACE: Beginning interface', c_key_obs_diags, c_key_geovals, c_key_self
call fckit_log % info(output_message)

call ufo_gnssro_BendMetOffice_registry % get(c_key_self, self)
call ufo_geovals_registry % get(c_key_obs_diags, obs_diags)
call ufo_geovals_registry % get(c_key_geovals, geovals)

call self%simobs(geovals, c_obsspace, c_nvars, c_nlocs, c_hofx, obs_diags)

write(output_message, *) 'TRACE: Finishing interface'
call fckit_log % info(output_message)

end subroutine ufo_gnssro_bendmetoffice_simobs_c

! ------------------------------------------------------------------------------

end module ufo_gnssro_bendmetoffice_mod_c
