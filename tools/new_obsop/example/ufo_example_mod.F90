! (C) Copyright 2017-2018 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

!> Fortran module for example observation operator

module ufo_example_mod

 use oops_variables_mod
 use obs_variables_mod
 use ufo_vars_mod

 implicit none
 private

!> Fortran derived type for the observation type
! TODO: fill in if needed
 type, public :: ufo_example
 private
   type(obs_variables), public :: obsvars
   type(oops_variables), public :: geovars
 contains
   procedure :: setup  => ufo_example_setup
   procedure :: simobs => ufo_example_simobs
   final :: destructor
 end type ufo_example

contains

! ------------------------------------------------------------------------------
! TODO: add setup of your observation operator (optional)
subroutine ufo_example_setup(self, f_conf)
use fckit_configuration_module, only: fckit_configuration
implicit none
class(ufo_example), intent(inout)     :: self

! TODO: consider whether passing the Configuration object to this function
! is necessary. If only a small number of parameters are used,
! you could pass them in directly instead. In that case you can modify the
! interface appropriately.
type(fckit_configuration), intent(in) :: f_conf


! TODO: add input variables (requested from the model)
! self%geovars%push_back("variable name")

end subroutine ufo_example_setup

! ------------------------------------------------------------------------------
! TODO: add cleanup of your observation operator (optional)
subroutine destructor(self)
implicit none
type(ufo_example), intent(inout) :: self

end subroutine destructor

! ------------------------------------------------------------------------------
! TODO: put code for your nonlinear observation operator in this routine
! Code in this routine is for example only, please remove and replace
subroutine ufo_example_simobs(self, geovals, obss, nvars, nlocs, hofx)
use kinds
use ufo_geovals_mod, only: ufo_geovals, ufo_geoval, ufo_geovals_get_var
use iso_c_binding
use obsspace_mod
implicit none
class(ufo_example), intent(in)    :: self
integer, intent(in)               :: nvars, nlocs
type(ufo_geovals),  intent(in)    :: geovals
real(c_double),     intent(inout) :: hofx(nvars, nlocs)
type(c_ptr), value, intent(in)    :: obss

! Local variables
!type(ufo_geoval), pointer :: geoval
!real(kind_real), dimension(:), allocatable :: obss_metadata

! check if some variable is in geovals and get it (var_tv is defined in ufo_vars_mod)
!call ufo_geovals_get_var(geovals, var_tv, geoval)

! get some metadata from obsspace
!allocate(obss_metadata(nlocs))
!call obsspace_get_db(obss, "MetaData", "some_metadata", obss_metadata)

! put observation operator code here


end subroutine ufo_example_simobs


! ------------------------------------------------------------------------------

end module ufo_example_mod
