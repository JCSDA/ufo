! (C) Copyright 2017-2018 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

!> Fortran module for example observation operator

module ufo_example_mod

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
! TODO: fill in if needed
 type, extends(ufo_basis), public :: ufo_example
 private
   integer, public :: nvars_in, nvars_out
   character(len=max_string), public, allocatable :: varin(:)
   character(len=max_string), public, allocatable :: varout(:)
 contains
   procedure :: setup  => ufo_example_setup
   procedure :: delete => ufo_example_delete
   procedure :: simobs => ufo_example_simobs
 end type ufo_example

contains

! ------------------------------------------------------------------------------
! TODO: add setup of your observation operator (optional)
subroutine ufo_example_setup(self, c_conf)
implicit none
class(ufo_example), intent(inout) :: self
type(c_ptr),        intent(in)    :: c_conf
! TODO: add input variables (requested from the model) and
!           output variables (simulated by the obs operator)
  self%nvars_in  = 0
  self%nvars_out = 0

end subroutine ufo_example_setup

! ------------------------------------------------------------------------------
! TODO: add cleanup of your observation operator (optional)
subroutine ufo_example_delete(self)
implicit none
class(ufo_example), intent(inout) :: self

  if (allocated(self%varout)) deallocate(self%varout)
  if (allocated(self%varin))  deallocate(self%varin)

end subroutine ufo_example_delete

! ------------------------------------------------------------------------------
! TODO: put code for your nonlinear observation operator in this routine
! Code in this routine is for example only, please remove and replace
subroutine ufo_example_simobs(self, geovals, hofx, obss)
implicit none
class(ufo_example), intent(in)    :: self
type(ufo_geovals),  intent(in)    :: geovals
real(c_double),     intent(inout) :: hofx(:)
type(c_ptr), value, intent(in)    :: obss

! Local variables
type(ufo_geoval), pointer :: geoval
integer :: nlocs
real(kind_real), dimension(:), allocatable :: obss_metadata

! check if some variable is in geovals and get it (var_tv is defined in ufo_vars_mod)
call ufo_geovals_get_var(geovals, var_tv, geoval)

! get some metadata from obsspace
nlocs = obsspace_get_nlocs(obss)
allocate(obss_metadata(nlocs))
call obsspace_get_db(obss, "MetaData", "some_metadata", obss_metadata)

! put observation operator code here


end subroutine ufo_example_simobs


! ------------------------------------------------------------------------------

end module ufo_example_mod
