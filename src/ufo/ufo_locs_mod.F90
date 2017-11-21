! (C) Copyright 2009-2016 ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
! In applying this licence, ECMWF does not waive the privileges and immunities 
! granted to it by virtue of its status as an intergovernmental organisation nor
! does it submit to any jurisdiction.

!> Fortran module handling observation locations

module ufo_locs_mod

use ufo_obs_vectors
use kinds

implicit none
private
public :: ufo_locs, ufo_locs_setup, ufo_locs_delete

! ------------------------------------------------------------------------------

!> Fortran derived type to hold observation locations
type :: ufo_locs
  integer :: nloc
  real(kind=kind_real), allocatable :: xyz(:,:)
end type ufo_locs

! ------------------------------------------------------------------------------
contains

! ------------------------------------------------------------------------------

subroutine ufo_locs_setup(self, lvec)
implicit none
type(ufo_locs), intent(inout) :: self
type(obs_vector), intent(in) :: lvec
integer :: jc, jo

self%nloc=lvec%nobs
allocate(self%xyz(3,self%nloc))
do jo=1,self%nloc
  do jc=1,3
    self%xyz(jc,jo)=lvec%values(jo)
  enddo
enddo

end subroutine ufo_locs_setup

! ------------------------------------------------------------------------------

subroutine ufo_locs_delete(self)
implicit none
type(ufo_locs), intent(inout) :: self

deallocate(self%xyz)

end subroutine ufo_locs_delete

! ------------------------------------------------------------------------------

end module ufo_locs_mod
