! (C) Copyright 2017-2018 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

!> Fortran module to implement example check

module ufo_examplecheck_mod

use iso_c_binding
use kinds
use ufo_geovals_mod
use obsspace_mod
use config_mod

implicit none
public :: ufo_examplecheck, ufo_examplecheck_create, ufo_examplecheck_delete, ufo_examplecheck_prior, ufo_examplecheck_post
private

! ------------------------------------------------------------------------------
!> TODO: fill in this type
type :: ufo_examplecheck
end type ufo_examplecheck

! ------------------------------------------------------------------------------
contains
! ------------------------------------------------------------------------------

subroutine ufo_examplecheck_create(self, obspace, conf)
implicit none
type(ufo_examplecheck), intent(inout) :: self
type(c_ptr), value, intent(in)   :: obspace
type(c_ptr), intent(in)          :: conf

end subroutine ufo_examplecheck_create

! ------------------------------------------------------------------------------

subroutine ufo_examplecheck_delete(self)
implicit none
type(ufo_examplecheck), intent(inout) :: self

end subroutine ufo_examplecheck_delete

! ------------------------------------------------------------------------------

subroutine ufo_examplecheck_prior(self, geovals)
implicit none
type(ufo_examplecheck), intent(in) :: self
type(ufo_geovals), intent(in) :: geovals

end subroutine ufo_examplecheck_prior

! ------------------------------------------------------------------------------

subroutine ufo_examplecheck_post(self, hofx)
implicit none
type(ufo_examplecheck), intent(in) :: self
real(c_double),  intent(in)        :: hofx(:)

end subroutine ufo_examplecheck_post

! ------------------------------------------------------------------------------

end module ufo_examplecheck_mod
