! (C) Copyright 2017-2018 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

! Fortran module for identity observation operator
!---------------------------------------------------------------------------------------------------
module ufo_identity_mod

 use iso_c_binding
 use config_mod
 use kinds

 use ufo_geovals_mod, only: &
    ufo_geovals,            &
    ufo_geoval,             &
    ufo_geovals_get_var
 
 use ufo_geovals_mod_c, only: &
    ufo_geovals_registry
 
 use ufo_basis_mod, only: &
     ufo_basis
     
 use ufo_vars_mod
 use obsspace_mod

 implicit none
 private

 integer, parameter :: max_string=800

! Fortran derived type for the observation type
!---------------------------------------------------------------------------------------------------
 type, extends(ufo_basis) :: ufo_identity
 private
    integer, public :: nvars
    character(len=max_string), public, allocatable :: varin(:)
    character(len=max_string), public, allocatable :: varout(:)
 contains
   procedure :: setup  => ufo_identity_setup_
   procedure :: simobs => ufo_identity_simobs_
   final :: destructor
 end type ufo_identity

contains

! ------------------------------------------------------------------------------
subroutine ufo_identity_setup_(self, c_conf)
   use config_mod
   implicit none
   class(ufo_identity), intent(inout) :: self
   type(c_ptr),        intent(in)     :: c_conf

   integer :: ii

  !> Size of variables
  self%nvars = size(config_get_string_vector(c_conf, max_string, "variables"))
  
  !> Allocate varout: variables in the observation vector
  allocate(self%varout(self%nvars))
  
  !> Read variable list and store in varout
  self%varin = config_get_string_vector(c_conf, max_string, "variables")
  
  ! -----------------------------------------------------------------------------
  !> TODO: No need to set-up varout here
  !> -----------------------------------------------------------------------------
  !> Allocate varin: variables we need from the model
  !> need additional slot to hold vertical coord.
  allocate(self%varout(self%nvars+1))
  
  !> Set vars_in based on vars_out
  do ii = 1, self%nvars
     self%varout(ii) = self%varin(ii)
  enddo
  

end subroutine ufo_identity_setup_


! ------------------------------------------------------------------------------
subroutine ufo_identity_simobs_(self, geovals, hofx, obss)
  implicit none
  class(ufo_identity), intent(in)    :: self
  type(ufo_geovals),  intent(in)     :: geovals
  real(c_double),     intent(inout)  :: hofx(:)
  type(c_ptr), value, intent(in)     :: obss
  
  
  integer :: iobs, ivar, nlocs
  type(ufo_geoval), pointer :: point
  character(len=MAXVARLEN) :: geovar
  
  !> Get the observation vertical coordinates
  nlocs = obsspace_get_nlocs(obss)

  do ivar = 1, self%nvars
    !> Get the name of input variable in geovals
    geovar = self%varin(ivar)

    !> Get profile for this variable from geovals
    call ufo_geovals_get_var(geovals, geovar, point)
    
    !> Here we just apply a identity hofx
    do iobs = 1, nlocs
        hofx(ivar + (iobs-1)*self%nvars) = point%vals(1,iobs)
    enddo
  enddo

end subroutine ufo_identity_simobs_


! ------------------------------------------------------------------------------
subroutine  destructor(self)
  type(ufo_identity), intent(inout) :: self
  if (allocated(self%varout)) deallocate(self%varout)
  if (allocated(self%varin)) deallocate(self%varin)
end subroutine destructor



end module ufo_identity_mod
