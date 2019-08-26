! (C) Copyright 2017-2018 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

!> Fortran module for radarreflectivity observation operator

module ufo_radarreflectivity_mod

 use fckit_configuration_module, only: fckit_configuration
 use iso_c_binding
 use kinds

 use ufo_geovals_mod, only: ufo_geovals, ufo_geoval, ufo_geovals_get_var
 use ufo_vars_mod
 use obsspace_mod

 implicit none
 private
 integer, parameter :: max_string=800

!> Fortran derived type for the observation type
! DONE 
 type, public :: ufo_radarreflectivity
 private
   integer, public :: nvars_in, nvars_out
   character(len=max_string), public, allocatable :: varin(:)
   character(len=max_string), public, allocatable :: varout(:)
 contains
   procedure :: setup  => ufo_radarreflectivity_setup
   procedure :: simobs => ufo_radarreflectivity_simobs
   final :: destructor
 end type ufo_radarreflectivity

 character(len=maxvarlen), dimension(1), parameter :: varin_default = (/var_refl/)

contains

! ------------------------------------------------------------------------------
! Done
subroutine ufo_radarreflectivity_setup(self, f_conf, vars)
implicit none
class(ufo_radarreflectivity), intent(inout)     :: self
type(fckit_configuration), intent(in) :: f_conf
character(len=MAXVARLEN), dimension(:), intent(inout) :: vars

  self%nvars_out = size(vars)
  allocate(self%varout(self%nvars_out))
  self%varout = vars

  self%nvars_in  = size(varin_default)
  allocate(self%varin(self%nvars_in))
  self%varin     = varin_default

end subroutine ufo_radarreflectivity_setup

! ------------------------------------------------------------------------------
! Done
subroutine destructor(self)
implicit none
type(ufo_radarreflectivity), intent(inout) :: self

  if (allocated(self%varout)) deallocate(self%varout)
  if (allocated(self%varin))  deallocate(self%varin)

end subroutine destructor

! ------------------------------------------------------------------------------
! TODO: put code for your nonlinear observation operator in this routine
! Code in this routine is for radarreflectivity only, please remove and replace
subroutine ufo_radarreflectivity_simobs(self, geovals, obss, nvars, nlocs, hofx)
implicit none
class(ufo_radarreflectivity), intent(in)    :: self
integer, intent(in)               :: nvars, nlocs
type(ufo_geovals),  intent(in)    :: geovals
real(c_double),     intent(inout) :: hofx(nvars, nlocs)
type(c_ptr), value, intent(in)    :: obss

! Local variables
type(ufo_geoval), pointer :: geoval
real(kind_real), dimension(:), allocatable :: obss_metadata

! check if some variable is in geovals and get it (var_tv is defined in ufo_vars_mod)
call ufo_geovals_get_var(geovals, var_tv, geoval)

! get some metadata from obsspace
allocate(obss_metadata(nlocs))
call obsspace_get_db(obss, "MetaData", "some_metadata", obss_metadata)

! put observation operator code here


end subroutine ufo_radarreflectivity_simobs


! ------------------------------------------------------------------------------

end module ufo_radarreflectivity_mod
