! (C) Copyright 2017-2018 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

!> Fortran module for identity tl/ad observation operator

module ufo_identity_tlad_mod

 use iso_c_binding
 use config_mod
 use kinds

 use ufo_geovals_mod, only: &
   ufo_geovals,             &
   ufo_geoval,              &
   ufo_geovals_get_var
   
 use ufo_basis_tlad_mod, only: &
   ufo_basis_tlad
 
 use ufo_geovals_mod_c, only: &
   ufo_geovals_registry
   
 use ufo_vars_mod
 use obsspace_mod

 implicit none

 integer, parameter :: max_string=800
 
 ! ------------------------------------------------------------------------------

 !> Fortran derived type for the tl/ad observation operator
 type, extends(ufo_basis_tlad), public :: ufo_identity_tlad
 private
   integer          :: nval, nlocs
   integer, public  :: nvars
   character(len=max_string), public, allocatable :: varin(:)
   character(len=max_string), public, allocatable :: varout(:)
 contains
  procedure :: setup  => identity_tlad_setup_
  procedure :: delete  => identity_tlad_delete_
  procedure :: settraj => identity_tlad_settraj_
  procedure :: simobs_tl  => identity_simobs_tl_
  procedure :: simobs_ad  => identity_simobs_ad_
  final :: destructor
 end type ufo_identity_tlad

contains

! ------------------------------------------------------------------------------
subroutine identity_tlad_setup_(self, c_conf)

  use config_mod
  implicit none
  class(ufo_identity_tlad), intent(inout) :: self
  type(c_ptr), intent(in)    :: c_conf

  integer :: ii

  !> Size of variables
  self%nvars = size(config_get_string_vector(c_conf, max_string, "variables"))
  
  !> Allocate varout
  allocate(self%varin(self%nvars))
  
  !> Read variable list and store in varout
  self%varin = config_get_string_vector(c_conf, max_string, "variables")
  
  !> Allocate varin, need additional slot to hold vertical coord.
  allocate(self%varout(self%nvars))
  
  !> Set vars_in based on vars_out
  do ii = 1, self%nvars
       self%varout(ii) = self%varin(ii)
  enddo

end subroutine identity_tlad_setup_


! ------------------------------------------------------------------------------
subroutine identity_tlad_settraj_(self, geovals, obss)
implicit none
class(ufo_identity_tlad), intent(inout) :: self
type(ufo_geovals),       intent(in)    :: geovals
type(c_ptr), value,      intent(in)    :: obss

! since observation operator is linear, don't care about trajectory itself

end subroutine identity_tlad_settraj_

! ------------------------------------------------------------------------------
subroutine identity_simobs_tl_(self, geovals, hofx, obss)
  implicit none
  class(ufo_identity_tlad), intent(in)    :: self
  type(ufo_geovals),  intent(in)     :: geovals
  real(c_double),     intent(inout)  :: hofx(:)
  type(c_ptr), value, intent(in)     :: obss
  
  
  integer :: iobs, ivar, nlocs
  type(ufo_geoval), pointer :: point
  character(len=MAXVARLEN) :: geovar
  
  !> Get the observation vertical coordinates
  nlocs = obsspace_get_nlocs(obss)

  do ivar = 1, self%nvars
    !> Get the name of input variable in geovals
    geovar = self%varin(ivar)

    !> Get profile for this variable from geovals
    call ufo_geovals_get_var(geovals, geovar, point)
    
    !> Here we just apply a identity hofx
    do iobs = 1, nlocs
        hofx(ivar + (iobs-1)*self%nvars) = point%vals(1,iobs)
    enddo
  enddo

end subroutine identity_simobs_tl_

! ------------------------------------------------------------------------------
subroutine identity_simobs_ad_(self, geovals, hofx, obss)
  implicit none
  class(ufo_identity_tlad), intent(in)    :: self
  type(ufo_geovals),       intent(inout) :: geovals
  real(c_double),          intent(in)    :: hofx(:)
  type(c_ptr), value,      intent(in)    :: obss

  integer :: iobs, ivar, nlocs
  type(ufo_geoval), pointer :: point
  character(len=MAXVARLEN)  :: geovar

  if (.not. geovals%linit ) geovals%linit=.true.

  nlocs = obsspace_get_nlocs(obss)
 
  do ivar = 1, self%nvars
    !> Get the name of input variable in geovals
    geovar = self%varin(ivar)

    !> Get profile for this variable from geovals
    call ufo_geovals_get_var(geovals, geovar, point)
    
    
    if (.not.(allocated(point%vals))) then
       point%nval=1
       allocate(point%vals(1,size(hofx,1)))
       point%vals = 0.0
    end if

    ! backward obs operator
    do iobs = 1, nlocs
      point%vals(1,iobs) = hofx(ivar + (iobs-1)*self%nvars) 
    enddo
  enddo
  
end subroutine identity_simobs_ad_



! ------------------------------------------------------------------------------
subroutine identity_tlad_delete_(self)
  implicit none
  class(ufo_identity_tlad), intent(inout) :: self
  self%nval = 0
  self%ltraj = .false.
end subroutine identity_tlad_delete_

! ------------------------------------------------------------------------------
subroutine  destructor(self)
  type(ufo_identity_tlad), intent(inout)  :: self
  self%nval = 0
  self%nlocs = 0
  self%nvars = 0
  self%ltraj = .false.
  if (allocated(self%varin)) deallocate(self%varin)
  if (allocated(self%varout)) deallocate(self%varout)
end subroutine destructor


end module ufo_identity_tlad_mod
