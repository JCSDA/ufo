!
! (C) Copyright 2017 UCAR
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 

!> Fortran example module for observation space

! TODO: replace example with your_obsspace_name through the file

module ufo_obs_example_mod

use kinds
use fckit_log_module, only : fckit_log

implicit none
private
integer, parameter :: max_string=800

public ufo_obs_example
public ufo_obs_example_setup, ufo_obs_example_delete
public ufo_obs_example_read, ufo_obs_example_generate
public ufo_obs_example_getlocs

! ------------------------------------------------------------------------------

!> Fortran derived type to hold observation space info
! TODO: fill in, below is just an example
type :: ufo_obs_example
  integer :: nobs
end type ufo_obs_example

! ------------------------------------------------------------------------------

contains
! ------------------------------------------------------------------------------
! TODO: replace the below function with your constructor of obsspace
subroutine ufo_obs_example_setup(self, nobs)
implicit none
type(ufo_obs_example), intent(inout) :: self
integer, intent(in) :: nobs

call ufo_obs_example_delete(self)

self%nobs = nobs

end subroutine ufo_obs_example_setup

! ------------------------------------------------------------------------------
! TODO: replace the below function with your destructor of obsspace
subroutine ufo_obs_example_delete(self)
implicit none
type(ufo_obs_example), intent(inout) :: self

end subroutine ufo_obs_example_delete

! ------------------------------------------------------------------------------
! TODO: replace the below function with your random obs generator
subroutine ufo_obs_example_generate(self, nobs)
implicit none
type(ufo_obs_example), intent(inout) :: self
integer, intent(in) :: nobs

end subroutine ufo_obs_example_generate

! ------------------------------------------------------------------------------
! TODO: replace the below function with your obsspace read
subroutine ufo_obs_example_read(filename, self)
implicit none
character(max_string), intent(in)   :: filename
type(ufo_obs_example), intent(inout), target :: self

end subroutine ufo_obs_example_read

! ------------------------------------------------------------------------------
! TODO: replace the below function with your obsspace get locations function
subroutine ufo_obs_example_getlocs(self, locs)
use ufo_locs_mod
implicit none
type(ufo_obs_example), intent(in) :: self
type(ufo_locs), intent(inout) :: locs

end subroutine ufo_obs_example_getlocs

! ------------------------------------------------------------------------------

end module ufo_obs_example_mod
