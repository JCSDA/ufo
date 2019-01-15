! (C) Copyright 2017-2018 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

!> Fortran module for seasurfacetemp observation operator

module ufo_seasurfacetemp_mod

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
 type, extends(ufo_basis), public :: ufo_seasurfacetemp
 private
 contains
   procedure :: setup  => ufo_seasurfacetemp_setup
   procedure :: delete => ufo_seasurfacetemp_delete
   procedure :: simobs => ufo_seasurfacetemp_simobs
 end type ufo_seasurfacetemp

contains

! ------------------------------------------------------------------------------
subroutine ufo_seasurfacetemp_setup(self, c_conf)
implicit none
class(ufo_seasurfacetemp), intent(inout) :: self
type(c_ptr),        intent(in)    :: c_conf

end subroutine ufo_seasurfacetemp_setup

! ------------------------------------------------------------------------------
subroutine ufo_seasurfacetemp_delete(self)
implicit none
class(ufo_seasurfacetemp), intent(inout) :: self

end subroutine ufo_seasurfacetemp_delete

! ------------------------------------------------------------------------------
subroutine ufo_seasurfacetemp_simobs(self, geovals, hofx, obss)
use ufo_marine_ncutils
implicit none
class(ufo_seasurfacetemp), intent(in)    :: self
type(ufo_geovals),  intent(in)    :: geovals
real(c_double),     intent(inout) :: hofx(:)
type(c_ptr), value, intent(in)    :: obss

    character(len=*), parameter :: myname_="ufo_seasurfacetemp_simobs"
    character(max_string) :: err_msg

    integer :: iobs
    type(ufo_geoval), pointer :: geoval_sst

    ! Netcdf stuff to write out geovals
    character(len=120) :: filename="sst_obs-2018-04-15_geovals.nc"
    integer(kind=4) :: iNcid
    integer(kind=4) :: iDimStation_ID, iDimLev_ID
    integer(kind=4) :: iVarLev_ID, iVarGOM_ID
    integer :: nlev,nobs

    type(diag_marine_obs) :: sst_out 
    
    ! check if nobs is consistent in geovals & hofx
    if (geovals%nobs /= size(hofx,1)) then
       write(err_msg,*) myname_, ' error: nobs inconsistent!'
       call abor1_ftn(err_msg)
    endif

    ! check if sst variables is in geovals and get it
    call ufo_geovals_get_var(geovals, var_ocn_sst, geoval_sst)

    ! Information for temporary output file ---------------------------------------!
    filename='sst-test.nc'    
    call sst_out%init(size(hofx,1),filename)
    
    ! sst obs operator
    do iobs = 1, size(hofx,1)
       hofx(iobs) = geoval_sst%vals(1,iobs)

       ! Output information:
       sst_out%diag(iobs)%Station_ID         = 9999
       sst_out%diag(iobs)%Observation_Type   = 999.9
       sst_out%diag(iobs)%Latitude           = 999.9
       sst_out%diag(iobs)%Longitude          = 999.9
       sst_out%diag(iobs)%Depth              = 999.9
       sst_out%diag(iobs)%Time               = 999.9
       sst_out%diag(iobs)%Observation        = 999.9
       sst_out%diag(iobs)%Observation_error  = 999.9
       sst_out%diag(iobs)%Obs_Minus_Forecast = 999.9
    enddo

    call sst_out%write_diag()
    call sst_out%write_geoval(var_ocn_sst,geoval_sst)
    call sst_out%finalize()

  end subroutine ufo_seasurfacetemp_simobs

end module ufo_seasurfacetemp_mod
