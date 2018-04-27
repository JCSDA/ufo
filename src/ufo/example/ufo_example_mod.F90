!a(C) Copyright 2017 UCAR
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 

!> Fortran example module for observation operator

! TODO: replace example with your_observation_operator_name through the file

module ufo_example_mod
  
  use ufo_obs_example_mod
  use ufo_obs_vectors
  use ufo_vars_mod
  use ufo_locs_mod
  use ufo_geovals_mod
  use kinds
  use vert_interp_mod

  implicit none
  public :: ufo_example
  public :: ufo_example_eqv
  private

!> Fortran derived type for the observation type
! TODO: fill in if needed
type :: ufo_example
end type ufo_example

! ------------------------------------------------------------------------------

contains
    
! ------------------------------------------------------------------------------
! TODO: replace below function with your observation operator.
! Some sample code is provided and should be removed/replaced/altered to your needs
subroutine ufo_example_eqv(self, geovals, hofx, obss)
implicit none
type(ufo_example), intent(in)     :: self
type(ufo_geovals), intent(in)      :: geovals
type(obs_vector),  intent(inout)   :: hofx
type(ufo_obs_example), intent(in) :: obss

character(len=MAXVARLEN) :: varname

character(len=*), parameter :: myname_="ufo_example_eqv"

type(ufo_geoval), pointer :: geoval

! check if nobs is consistent in geovals & hofx
if (geovals%nobs /= hofx%nobs) then
  write(err_msg,*) myname_, ' error: nobs inconsistent!'
  call abor1_ftn(err_msg)
endif

varname = "some_variable_name"
! check if some variable is in geovals and get it
if (.not. ufo_geovals_get_var(geovals, varname, geoval)) then
  write(err_msg,*) myname_, trim(varname), ' doesnt exist'
  call abor1_ftn(err_msg)
endif

! put observation operator code here

end subroutine ufo_example_eqv

! ------------------------------------------------------------------------------

end module ufo_example_mod
