!
!  (C) Copyright 2017 UCAR
!  
!  This software is licensed under the terms of the Apache Licence Version 2.0
!  which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
!

module ufo_vars_mod

implicit none
private
public :: ufo_vars, ufo_vars_setup, ufo_vars_clone, ufo_vars_delete
public :: ufo_vars_getindex, ufo_vars_nvars
public :: MAXVARLEN

integer, parameter :: MAXVARLEN=24
! ------------------------------------------------------------------------------

!> Fortran derived type to represent model variables
type :: ufo_vars
  integer :: nv
  character(len=MAXVARLEN), allocatable :: fldnames(:) !< Variable identifiers
end type ufo_vars

! ------------------------------------------------------------------------------
contains
! ------------------------------------------------------------------------------

subroutine ufo_vars_setup(self, cvars)
implicit none
type(ufo_vars), intent(inout) :: self
character(len=MAXVARLEN), intent(in) :: cvars(:)

self%nv = size(cvars)

! TODO: a check on whether this var is in the list of defined vars
allocate(self%fldnames(self%nv))
self%fldnames(:) = cvars(:)

end subroutine ufo_vars_setup

! ------------------------------------------------------------------------------

subroutine ufo_vars_clone(self, other)
implicit none
type(ufo_vars), intent(in)    :: self
type(ufo_vars), intent(inout) :: other

call ufo_vars_delete(other)
other%nv = self%nv
allocate(other%fldnames(other%nv))
other%fldnames(:) = self%fldnames(:)

end subroutine ufo_vars_clone

! ------------------------------------------------------------------------------

subroutine ufo_vars_delete(self)
implicit none
type(ufo_vars), intent(inout) :: self

if (allocated(self%fldnames)) deallocate(self%fldnames)
self%nv = 0

end subroutine ufo_vars_delete

! ------------------------------------------------------------------------------

integer function ufo_vars_getindex(self, varname)
implicit none
type(ufo_vars), intent(in)       :: self
character(MAXVARLEN), intent(in) :: varname

integer :: ivar

ufo_vars_getindex = -1

do ivar = 1, self%nv
  if (self%fldnames(ivar) == varname) then
    ufo_vars_getindex = ivar
    exit
  endif
enddo

end function ufo_vars_getindex

! ------------------------------------------------------------------------------

integer function ufo_vars_nvars(self) 
implicit none
type(ufo_vars), intent(in) :: self

ufo_vars_nvars = self%nv

end function ufo_vars_nvars

! ------------------------------------------------------------------------------

end module ufo_vars_mod
