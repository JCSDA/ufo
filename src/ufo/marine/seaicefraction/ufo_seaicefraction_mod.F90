! (C) Copyright 2017-2018 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

!> Fortran module for seaicefraction observation operator

module ufo_seaicefraction_mod

 use fckit_configuration_module, only: fckit_configuration 
 use iso_c_binding
 use kinds

 use ufo_geovals_mod, only: ufo_geovals, ufo_geoval, ufo_geovals_get_var
 use ufo_basis_mod, only: ufo_basis
 use ufo_vars_mod
 use obsspace_mod

 implicit none
 private

 integer, parameter :: max_string=800

!> Fortran derived type for the observation type
 type, extends(ufo_basis), public :: ufo_seaicefraction
 private
 contains
   procedure :: setup  => ufo_seaicefraction_setup
   procedure :: delete => ufo_seaicefraction_delete
   procedure :: simobs => ufo_seaicefraction_simobs
 end type ufo_seaicefraction

contains

! ------------------------------------------------------------------------------
subroutine ufo_seaicefraction_setup(self, f_conf)
implicit none
class(ufo_seaicefraction), intent(inout) :: self
type(fckit_configuration), intent(in)    :: f_conf

end subroutine ufo_seaicefraction_setup

! ------------------------------------------------------------------------------
subroutine ufo_seaicefraction_delete(self)
implicit none
class(ufo_seaicefraction), intent(inout) :: self

end subroutine ufo_seaicefraction_delete

! ------------------------------------------------------------------------------
! Code in this routine is for seaicefraction only, please remove and replace
subroutine ufo_seaicefraction_simobs(self, geovals, hofx, obss)
use ufo_marine_ncutils
implicit none
class(ufo_seaicefraction), intent(in)    :: self
type(ufo_geovals),  intent(in)    :: geovals
real(c_double),     intent(inout) :: hofx(:)
type(c_ptr), value, intent(in)    :: obss

    character(len=*), parameter :: myname_="ufo_seaicefrac_simobs"
    character(max_string) :: err_msg

    integer :: iobs
    type(ufo_geoval), pointer :: geoval

    ! Netcdf stuff to write out geovals
    integer(kind=4) :: iNcid
    integer(kind=4) :: iDimStation_ID, iDimLev_ID
    integer(kind=4) :: iVarLev_ID, iVarGOM_ID
    integer :: ncat,nlocs

    ! check if nlocs is consistent in geovals & hofx
    if (geovals%nlocs /= size(hofx,1)) then
       write(err_msg,*) myname_, ' error: nlocs inconsistent!'
       call abor1_ftn(err_msg)
    endif

    ! check if sea ice fraction variables is in geovals and get it
    call ufo_geovals_get_var(geovals, var_seaicefrac, geoval)

    ! total sea ice fraction obs operator
    do iobs = 1, size(hofx,1)
       hofx(iobs) = sum(geoval%vals(:,iobs))
    enddo

end subroutine ufo_seaicefraction_simobs

end module ufo_seaicefraction_mod
