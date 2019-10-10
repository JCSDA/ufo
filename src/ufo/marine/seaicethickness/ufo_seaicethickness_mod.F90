! (C) Copyright 2017-2018 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

!> Fortran module for seaicethickness observation operator

module ufo_seaicethickness_mod

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
 type, extends(ufo_basis), public :: ufo_seaicethickness
 private
 contains
   procedure :: setup  => ufo_seaicethickness_setup
   procedure :: delete => ufo_seaicethickness_delete
   procedure :: simobs => ufo_seaicethickness_simobs
 end type ufo_seaicethickness

contains

! ------------------------------------------------------------------------------
subroutine ufo_seaicethickness_setup(self, f_conf)
implicit none
class(ufo_seaicethickness), intent(inout) :: self
type(fckit_configuration),  intent(in)    :: f_conf

end subroutine ufo_seaicethickness_setup

! ------------------------------------------------------------------------------
subroutine ufo_seaicethickness_delete(self)
implicit none
class(ufo_seaicethickness), intent(inout) :: self

end subroutine ufo_seaicethickness_delete

! ------------------------------------------------------------------------------
subroutine ufo_seaicethickness_simobs(self, geovals, hofx, obss)
implicit none
class(ufo_seaicethickness), intent(in)    :: self
type(ufo_geovals),  intent(in)    :: geovals
real(c_double),     intent(inout) :: hofx(:)
type(c_ptr), value, intent(in)    :: obss

    character(len=*), parameter :: myname_="ufo_seaicethick_simobs"
    character(max_string) :: err_msg

    integer :: iobs, icat, ncat
    type(ufo_geoval), pointer :: icethick, icefrac
    
    ! check if nlocs is consistent in geovals & hofx
    if (geovals%nlocs /= size(hofx,1)) then
       write(err_msg,*) myname_, ' error: nlocs inconsistent!'
       call abor1_ftn(err_msg)
    endif

    ! check if sea ice fraction variable is in geovals and get it
    call ufo_geovals_get_var(geovals, var_seaicefrac, icefrac)
    ! check if sea ice thickness variable is in geovals and get it
    call ufo_geovals_get_var(geovals, var_seaicethick, icethick)

    ncat = icefrac%nval
    hofx = 0.0
    ! total sea ice fraction obs operator
    do iobs = 1, size(hofx,1)
       do icat = 1, ncat
          hofx(iobs) = hofx(iobs) + icefrac%vals(icat,iobs) * icethick%vals(icat,iobs)
       enddo
    enddo

end subroutine ufo_seaicethickness_simobs

end module ufo_seaicethickness_mod
