! (C) Copyright 2017-2018 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

!> Fortran module to implement example check

module ufo_example_mod

use fckit_configuration_module, only: fckit_configuration
use iso_c_binding
use kinds
use ufo_geovals_mod
use obsspace_mod
use ufo_vars_mod

implicit none
public :: ufo_example_create, ufo_example_delete, ufo_example_prior, ufo_example_post
private
integer, parameter :: max_string=800

! ------------------------------------------------------------------------------
!> TODO: fill in this type
type, public :: ufo_example
private
  character(len=max_string), public, allocatable :: geovars(:)
end type ufo_example

! ------------------------------------------------------------------------------
contains
! ------------------------------------------------------------------------------

subroutine ufo_example_create(self, f_conf)
implicit none
type(ufo_example), intent(inout)      :: self

! TODO: consider whether passing the Configuration object to this function
! is necessary. If only a small number of parameters are used,
! you could pass them in directly instead. In that case you can modify the
! interface appropriately.
type(fckit_configuration), intent(in) :: f_conf

! TODO: set self%geovars (list of variables to use from GeoVaLs) if needed

end subroutine ufo_example_create

! ------------------------------------------------------------------------------

subroutine ufo_example_delete(self)
implicit none
type(ufo_example), intent(inout) :: self

if (allocated(self%geovars))   deallocate(self%geovars)

end subroutine ufo_example_delete

! ------------------------------------------------------------------------------

subroutine ufo_example_prior(self, obspace, geovals)
implicit none
type(ufo_example),  intent(in) :: self
type(c_ptr), value, intent(in) :: obspace
type(ufo_geovals),  intent(in) :: geovals

end subroutine ufo_example_prior

! ------------------------------------------------------------------------------

subroutine ufo_example_post(self, obspace, nvars, nlocs, hofx, hofxbias, hofxdiags)
implicit none
type(ufo_example),  intent(in) :: self
type(c_ptr), value, intent(in) :: obspace
integer,            intent(in) :: nvars, nlocs
real(c_double),     intent(in) :: hofx(nvars, nlocs)
real(c_double),     intent(in) :: hofxbias(nvars, nlocs)
type(ufo_geovals),  intent(in) :: hofxdiags

end subroutine ufo_example_post

! ------------------------------------------------------------------------------

end module ufo_example_mod
