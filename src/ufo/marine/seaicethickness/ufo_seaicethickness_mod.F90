! (C) Copyright 2017-2018 UCAR
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 

!> Fortran module to handle ice concentration observations

module ufo_seaicethickness_mod

  use iso_c_binding
  use ufo_vars_mod
  use ufo_geovals_mod
  use kinds

  implicit none
  public :: ufo_seaicethickness
  public :: ufo_seaicethickness_simobs
  private
  integer, parameter :: max_string=800

  !> Fortran derived type for sea ice fraction observation operator
  type :: ufo_seaicethickness
  end type ufo_seaicethickness


  ! ------------------------------------------------------------------------------

contains

  ! ------------------------------------------------------------------------------

  subroutine ufo_seaicethickness_simobs(self, geovals, hofx)
    use ufo_marine_ncutils
    use ufo_vars_mod
    
    implicit none
    type(ufo_seaicethickness), intent(in) :: self
    type(ufo_geovals), intent(in)    :: geovals
    real(c_double),  intent(inout) :: hofx(:)

    character(len=*), parameter :: myname_="ufo_seaicethickness_simobs"
    character(max_string) :: err_msg

    integer :: iobs, icat, ncat
    type(ufo_geoval), pointer :: icethick, icefrac

    ! Netcdf stuff 
    character(len=120) :: filename !< name of outpu file for omf, lon, lat, ...
    character(len=MAXVARLEN) :: dim_name    
    type(diag_marine_obs) :: sit_out    
    
    ! check if nobs is consistent in geovals & hofx
    if (geovals%nobs /= size(hofx,1)) then
       write(err_msg,*) myname_, ' error: nobs inconsistent!'
       call abor1_ftn(err_msg)
    endif

    ! check if sea ice fraction variable is in geovals and get it
    call ufo_geovals_get_var(geovals, var_seaicefrac, icefrac)
    ! check if sea ice thickness variable is in geovals and get it
    call ufo_geovals_get_var(geovals, var_seaicethick, icethick)

    ! Information for temporary output file
    filename='sit-test.nc'    
    call sit_out%init(size(hofx,1),filename)
    
    ncat = icefrac%nval
    hofx = 0.0
    ! total sea ice fraction obs operator
    do iobs = 1, size(hofx,1)
       do icat = 1, ncat
          hofx(iobs) = hofx(iobs) + icefrac%vals(icat,iobs) * icethick%vals(icat,iobs) / 905.0
       enddo
    enddo

    dim_name="ncat"
    call sit_out%write_geoval(var_seaicefrac,icefrac,arg_dim_name=dim_name)
    call sit_out%write_geoval(var_seaicethick,icethick,arg_dim_name=dim_name)    
    call sit_out%finalize()

  end subroutine ufo_seaicethickness_simobs

end module ufo_seaicethickness_mod
