! (C) Copyright 2017-2018 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

!> Fortran module for radiancerttov observation operator

module ufo_radiancerttov_mod

 use iso_c_binding
 use config_mod
 use kinds

 use ufo_geovals_mod, only: ufo_geovals, ufo_geoval, ufo_geovals_get_var
 use ufo_basis_mod, only: ufo_basis
 use ufo_vars_mod
 use obsspace_mod

 implicit none
 private

!> Fortran derived type for the observation type
! TODO: fill in if needed
 type, extends(ufo_basis), public :: ufo_radiancerttov
 private
 contains
   procedure :: setup  => ufo_radiancerttov_setup
   procedure :: delete => ufo_radiancerttov_delete
   procedure :: simobs => ufo_radiancerttov_simobs
 end type ufo_radiancerttov

contains

! ------------------------------------------------------------------------------
! TODO: add setup of your observation operator (optional)
subroutine ufo_radiancerttov_setup(self, c_conf)
implicit none
class(ufo_radiancerttov), intent(inout) :: self
type(c_ptr),        intent(in)    :: c_conf

end subroutine ufo_radiancerttov_setup

! ------------------------------------------------------------------------------
! TODO: add cleanup of your observation operator (optional)
subroutine ufo_radiancerttov_delete(self)
implicit none
class(ufo_radiancerttov), intent(inout) :: self

end subroutine ufo_radiancerttov_delete

! ------------------------------------------------------------------------------
! TODO: put code for your nonlinear observation operator in this routine
! Code in this routine is for radiancerttov only, please remove and replace
subroutine ufo_radiancerttov_simobs(self, geovals, hofx, obss)
implicit none
class(ufo_radiancerttov), intent(in)    :: self
type(ufo_geovals),  intent(in)    :: geovals
real(c_double),     intent(inout) :: hofx(:)
type(c_ptr), value, intent(in)    :: obss


end subroutine ufo_radiancerttov_simobs

! ------------------------------------------------------------------------------

end module ufo_radiancerttov_mod
