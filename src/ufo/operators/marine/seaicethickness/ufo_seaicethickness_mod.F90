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
 use obs_variables_mod
 use missing_values_mod

 implicit none
 private

 integer, parameter :: max_string=800

!> Fortran derived type for the observation type
 type, extends(ufo_basis), public :: ufo_seaicethickness
    type(obs_variables), public :: obsvars
    real(kind=kind_real) :: rho_ice  = 905.0 !< [kg/m3]
    real(kind=kind_real) :: rho_snow = 330.0 !< [kg/m3]
    real(kind=kind_real) :: rho_water= 1000.0!< [kg/m3]   
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
integer :: ivar, nvars
character(max_string)  :: err_msg

nvars = self%obsvars%nvars()
if (nvars /= 1) then
  write(err_msg,*) 'ufo_seaicethickness_setup error: only variables size 1 supported!'
  call abor1_ftn(err_msg)
endif

end subroutine ufo_seaicethickness_setup

! ------------------------------------------------------------------------------
subroutine ufo_seaicethickness_delete(self)
implicit none
class(ufo_seaicethickness), intent(inout) :: self

end subroutine ufo_seaicethickness_delete

! ------------------------------------------------------------------------------
subroutine ufo_seaicethickness_simobs(self, geovals, hofx, obss)
use ufo_utils_mod, only: cmp_strings
implicit none
class(ufo_seaicethickness), intent(in)    :: self
type(ufo_geovals),  intent(in)    :: geovals
real(c_double),     intent(inout) :: hofx(:)
type(c_ptr), value, intent(in)    :: obss

    character(len=*), parameter :: myname_="ufo_seaicethick_simobs"
    character(max_string) :: err_msg
    real(c_double) :: missing

    integer :: iobs, icat, ncat
    type(ufo_geoval), pointer :: icethick, icefrac, snowthick
    real(kind=kind_real) :: rho_wiw, rho_wsw

    missing = missing_value(missing)

    ! check if nlocs is consistent in geovals & hofx
    if (geovals%nlocs /= size(hofx,1)) then
       write(err_msg,*) myname_, ' error: nlocs inconsistent!'
       call abor1_ftn(err_msg)
    endif 

    if (cmp_strings(self%obsvars%variable(1), "seaIceFreeboard")) then
       rho_wiw = (self%rho_water-self%rho_ice)/self%rho_water
       rho_wsw = (-self%rho_snow)/self%rho_water  
    endif

    ! check if sea ice fraction variable is in geovals and get it
    call ufo_geovals_get_var(geovals, var_seaicefrac, icefrac)
    ! check if snow thickness variable is in geovals and get it
    if (cmp_strings(self%obsvars%variable(1), "seaIceFreeboard")) &
       call ufo_geovals_get_var(geovals, var_seaicesnowthick, snowthick)
    ! check if sea ice thickness variable is in geovals and get it
    call ufo_geovals_get_var(geovals, var_seaicethick, icethick)

    ncat = icefrac%nval
    hofx = 0.0

    ! total sea ice fraction obs operator
    select case (trim(self%obsvars%variable(1)))
    case ("seaIceFreeboard")
       do iobs = 1, size(hofx,1)
          do icat = 1, ncat
             ! check for missing input values
             if (icefrac%vals(icat,iobs) == missing .or. &
                 icethick%vals(icat,iobs) == missing .or. &
                 snowthick%vals(icat,iobs) == missing ) then
                  hofx(iobs) = missing
                  exit
             end if

             hofx(iobs) = hofx(iobs)+ rho_wiw*icefrac%vals(icat,iobs) * icethick%vals(icat,iobs)&
                                    + rho_wsw*icefrac%vals(icat,iobs) * snowthick%vals(icat,iobs)
          enddo
       enddo
    case ("iceThickness")
       do iobs = 1, size(hofx,1)
          do icat = 1, ncat
             ! check for missing input values
             if (icefrac%vals(icat,iobs) == missing .or. &
                 icethick%vals(icat,iobs) == missing) then                  
               hofx(iobs) = missing
               exit
             end if

             hofx(iobs) = hofx(iobs) + icefrac%vals(icat,iobs) * icethick%vals(icat,iobs)
          enddo
       enddo
    case default
       write(err_msg,*) myname_, ' error: no match seaice thickness_option!'
       call abor1_ftn(err_msg)
    end select

end subroutine ufo_seaicethickness_simobs

end module ufo_seaicethickness_mod
