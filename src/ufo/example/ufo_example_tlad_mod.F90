! (C) Copyright 2017-2018 UCAR
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 

!> Fortran example module for tl/ad observation operator

! TODO: replace example with your_observation_operator_name through the file

module ufo_example_tlad_mod
  
  use ioda_obsdb_mod
  use ioda_obs_vectors
  use ufo_vars_mod
  use ioda_locs_mod
  use ufo_geovals_mod
  use kinds
  use vert_interp_mod

  implicit none
  public :: ufo_example_tlad
  public :: ufo_example_tlad_settraj
  public :: ufo_example_tlad_tl
  public :: ufo_example_tlad_ad
  public :: ufo_example_tlad_delete
  private

!> Fortran derived type for the tl/ad observation operator
! TODO: replace below type with what you need for your tl/ad observation operator
!       this type can hold information on trajectory, for example
! TODO: the below code is only an example and should be removed/replaced/altered 
!       to your needs
type :: ufo_example_tlad
   logical :: ltraj = .false. !< flag if trajectory was set
   type(geoval) :: gv         !< geoval for some variable from the trajectory
                              !  that is needed in tl/ad
   integer :: nval
   integer, allocatable :: wi(:) !< some array that is needed in tl/ad
end type ufo_example_tlad

! ------------------------------------------------------------------------------

contains
    
! ------------------------------------------------------------------------------
! TODO: replace below function with destructing your ufo_example_tlad type
!       (deallocating arrays, calling delete for type members, etc.)
! Some sample code is provided and should be removed/replaced/altered to your needs
subroutine ufo_example_tlad_delete(self)
implicit none
type(ufo_example_tlad), intent(inout) :: self

if (allocated(self%wi)) deallocate(self%wi)
self%ltraj = .false.

end subroutine ufo_example_tlad_delete

! ------------------------------------------------------------------------------
! TODO: replace below function with your set trajectory for tl/ad code
! TODO: the below code is only an example and should be removed/replaced/altered 
!       to your needs
subroutine ufo_example_tlad_settraj(self, geovals, obss)
implicit none
type(ufo_example_tlad), intent(inout) :: self
type(ufo_geovals),      intent(in)    :: geovals
type(ioda_obsdb),       intent(in)    :: obss

character(len=MAXVARLEN)    :: varname
character(len=*), parameter :: myname_="ufo_example_tlad_settraj"

type(ufo_geoval), pointer :: geoval

!Check if some variable is in geovals and get it
varname = "some_variable_name"
call ufo_geovals_get_var(geovals, varname, geoval)

!Copy the variable to the ufo_example_tlad type (save for future tl/ad)
self%gv = geoval
!Flag that trajectory was set
self%ltraj = .true.

end subroutine ufo_example_tlad_settraj

! ------------------------------------------------------------------------------
! TODO: replace below function with your tl observation operator.
! Note: this can use information saved from trajectory in your ufo_example_tlad type
! Input geovals parameter represents dx for tangent linear model
subroutine ufo_example_tlad_eqv_tl(self, geovals, hofx, obss)
implicit none
type(ufo_example_tlad), intent(in)     :: self
type(ufo_geovals),      intent(in)     :: geovals
type(obs_vector),       intent(inout)  :: hofx
type(ioda_obsdb),       intent(in) :: obss


end subroutine ufo_example_tlad_eqv_tl

! ------------------------------------------------------------------------------
! TODO: replace below function with your ad observation operator.
! Note: this can use information saved from trajectory in your ufo_example_tlad type
subroutine ufo_example_tlad_eqv_ad(self, geovals, hofx, obss)
implicit none
type(ufo_example_tlad), intent(in)     :: self
type(ufo_geovals),      intent(in)     :: geovals
type(obs_vector),       intent(inout)  :: hofx
type(ioda_obsdb),       intent(in)     :: obss


end subroutine ufo_example_tlad_eqv_ad

! ------------------------------------------------------------------------------

end module ufo_example_tlad_mod
