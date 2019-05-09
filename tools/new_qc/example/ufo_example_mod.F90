! (C) Copyright 2017-2018 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

!> Fortran module to implement example check

module ufo_example_mod

use iso_c_binding
use kinds
use ufo_geovals_mod
use obsspace_mod
use config_mod

implicit none
public :: ufo_example, ufo_example_create, ufo_example_delete, ufo_example_prior, ufo_example_post
private

! ------------------------------------------------------------------------------
!> TODO: fill in this type
type :: ufo_example
end type ufo_example

! ------------------------------------------------------------------------------
contains
! ------------------------------------------------------------------------------

subroutine ufo_example_create(self, conf)
implicit none
type(ufo_example), intent(inout) :: self
type(c_ptr), intent(in)          :: conf

end subroutine ufo_example_create

! ------------------------------------------------------------------------------

subroutine ufo_example_delete(self)
implicit none
type(ufo_example), intent(inout) :: self

end subroutine ufo_example_delete

! ------------------------------------------------------------------------------

subroutine ufo_example_prior(self, obspace, geovals)
implicit none
type(ufo_example),  intent(in) :: self
type(c_ptr), value, intent(in) :: obspace
type(ufo_geovals),  intent(in) :: geovals

end subroutine ufo_example_prior

! ------------------------------------------------------------------------------

subroutine ufo_example_post(self, obspace, hofx)
implicit none
type(ufo_example),  intent(in) :: self
type(c_ptr), value, intent(in) :: obspace
real(c_double),     intent(in) :: hofx(:)

end subroutine ufo_example_post

! ------------------------------------------------------------------------------

end module ufo_example_mod
