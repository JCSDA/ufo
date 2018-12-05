! (C) Copyright 2017-2018 UCAR
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 

!> Fortran module to handle ice concentration observations

module ufo_seaicefrac_mod

  use iso_c_binding
  use ufo_vars_mod
  use ufo_geovals_mod
  use kinds

  implicit none
  public :: ufo_seaicefrac
  public :: ufo_seaicefrac_simobs
  private
  integer, parameter :: max_string=800

  !> Fortran derived type for sea ice fraction observation operator
  type :: ufo_seaicefrac
  end type ufo_seaicefrac


  ! ------------------------------------------------------------------------------

contains

  ! ------------------------------------------------------------------------------

  subroutine ufo_seaicefrac_simobs(self, geovals, hofx)
    use ufo_marine_ncutils
    use ufo_vars_mod
        
    implicit none
    type(ufo_seaicefrac), intent(in) :: self
    type(ufo_geovals), intent(in)    :: geovals
    real(c_double),  intent(inout) :: hofx(:)

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

  end subroutine ufo_seaicefrac_simobs

end module ufo_seaicefrac_mod
