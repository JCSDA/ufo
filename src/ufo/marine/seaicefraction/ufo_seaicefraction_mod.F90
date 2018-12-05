! (C) Copyright 2017-2018 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

!> Fortran module for seaicefraction observation operator

module ufo_seaicefraction_mod

 use iso_c_binding
 use config_mod
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
subroutine ufo_seaicefraction_setup(self, c_conf)
implicit none
class(ufo_seaicefraction), intent(inout) :: self
type(c_ptr),        intent(in)    :: c_conf

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
    integer :: ncat,nobs

    ! netcdf stuff
    character(len=120) :: filename !< name of outpu file for omf, lon, lat, ...
    character(len=MAXVARLEN) :: dim_name    
    type(diag_marine_obs) :: sic_out    

    ! check if nobs is consistent in geovals & hofx
    if (geovals%nobs /= size(hofx,1)) then
       write(err_msg,*) myname_, ' error: nobs inconsistent!'
       call abor1_ftn(err_msg)
    endif

    ! check if sea ice fraction variables is in geovals and get it
    call ufo_geovals_get_var(geovals, var_seaicefrac, geoval)

    ! Information for temporary output file
    filename='sic-test.nc'    
    call sic_out%init(size(hofx,1),filename)
    
    ! total sea ice fraction obs operator
    do iobs = 1, size(hofx,1)
       hofx(iobs) = sum(geoval%vals(:,iobs))
    enddo

    dim_name="ncat"
    call sic_out%write_geoval(var_seaicefrac,geoval,arg_dim_name=dim_name)
    call sic_out%finalize()

end subroutine ufo_seaicefraction_simobs

end module ufo_seaicefraction_mod
