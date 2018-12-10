! (C) Copyright 2017-2018 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

!> Fortran module for identity observation operator

module ufo_identity_mod

 use iso_c_binding
 use config_mod
 use kinds

 use ufo_geovals_mod, only: ufo_geovals, ufo_geoval, ufo_geovals_get_var
 use ufo_basis_mod, only: ufo_basis
 use ufo_vars_mod
 use obsspace_mod

 implicit none
 private

 integer, parameter :: max_string=800

!> Fortran derived type for the observation type
 type, extends(ufo_basis), public :: ufo_identity
 private
 contains
   procedure :: setup  => ufo_identity_setup
   procedure :: delete => ufo_identity_delete
   procedure :: simobs => ufo_identity_simobs
 end type ufo_identity

contains

! ------------------------------------------------------------------------------
subroutine ufo_identity_setup(self, c_conf)
implicit none
class(ufo_identity), intent(inout) :: self
type(c_ptr),        intent(in)    :: c_conf

end subroutine ufo_identity_setup

! ------------------------------------------------------------------------------
subroutine ufo_identity_delete(self)
implicit none
class(ufo_identity), intent(inout) :: self

end subroutine ufo_identity_delete

! ------------------------------------------------------------------------------
subroutine ufo_identity_simobs(self, geovals, hofx, obss)
implicit none
class(ufo_identity), intent(in)    :: self
type(ufo_geovals),  intent(in)    :: geovals
real(c_double),     intent(inout) :: hofx(:)
type(c_ptr), value, intent(in)    :: obss

type(ufo_geoval), pointer :: geoval
integer :: iobs

! get the variable from geovals
call ufo_geovals_get_var(geovals, var_ocn_sst, geoval)

do iobs = 1, size(hofx,1)
  hofx(iobs) = geoval%vals(1,iobs)
enddo

end subroutine ufo_identity_simobs

end module ufo_identity_mod
