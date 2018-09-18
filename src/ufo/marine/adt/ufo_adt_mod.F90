! (C) Copyright 2017 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

!> Fortran module to handle adt observations

module ufo_adt_mod

  use ioda_obs_adt_mod
  use ioda_obs_vectors
  use ufo_vars_mod
  use ioda_locs_mod
  use ufo_geovals_mod
  use kinds
  use ncd_kinds, only:  i_kind,r_single,r_kind,r_double

  implicit none
  public :: ufo_adt
  public :: ufo_adt_simobs
  private
  integer, parameter :: max_string=800

  !> Fortran derived type for adt observation operator
  type :: ufo_adt
  end type ufo_adt

  type diag_adt_altimeter
     integer      :: Station_ID
     real(r_kind) :: Observation_Type
     real(r_kind) :: Latitude
     real(r_kind) :: Longitude
     real(r_kind) :: Time
     real(r_kind) :: Observation
     real(r_kind) :: Obs_Minus_Forecast
  end type diag_adt_altimeter

  ! ------------------------------------------------------------------------------

contains

  ! ------------------------------------------------------------------------------

  subroutine ufo_adt_simobs(self, geovals, hofx, obs_adt)

    use ufo_marine_ncutils
    
    implicit none
    
    type(ufo_adt), intent(in) :: self
    type(ufo_geovals), intent(in)    :: geovals
    type(ioda_obs_adt), intent(in) :: obs_adt     !< adt observations
    type(obs_vector),  intent(inout) :: hofx

    character(len=*), parameter :: myname_="ufo_adt_simobs"
    character(max_string) :: err_msg

    ! nc_diag stuff
    logical :: append
    character(len=120) :: filename !< name of outpu file for omf, lon, lat, ...
    !type(diag_adt_altimeter), allocatable :: adt_out(:)
    type(diag_marine_obs) :: adt_out    

    integer :: iobs
    real :: sum_obs,sum_hofx
    type(ufo_geoval), pointer :: geoval_adt

    ! check if nobs is consistent in geovals & hofx
    if (geovals%nobs /= hofx%nobs) then
       write(err_msg,*) myname_, ' error: nobs inconsistent!',geovals%nobs,hofx%nobs
       call abor1_ftn(err_msg)
    endif

    ! check if adt variable is in geovals and get it
    call ufo_geovals_get_var(geovals, var_abs_topo, geoval_adt)

    ! Information for temporary output file ---------------------------------------!
    filename='adt-test.nc'    
    call adt_out%init(hofx%nobs,filename)

    hofx%values = 0.0

    ! Compute offset
    sum_hofx=sum(geoval_adt%vals(1,:))
    sum_obs=sum(obs_adt%adt(:))
    print *,'ssh offset',(sum_obs-sum_hofx)/hofx%nobs

    ! adt obs operator
    do iobs = 1, hofx%nobs
       ! remove offset from hofx
       hofx%values(iobs) = geoval_adt%vals(1,iobs)+(sum_obs-sum_hofx)/hofx%nobs

       ! Output information:
       adt_out%diag(iobs)%Station_ID            = 1
       adt_out%diag(iobs)%Observation_Type      = 1.0
       adt_out%diag(iobs)%Latitude              = obs_adt%lat(iobs)
       adt_out%diag(iobs)%Longitude             = obs_adt%lon(iobs)
       adt_out%diag(iobs)%Time                  = 1.0
       adt_out%diag(iobs)%Observation           = obs_adt%adt(iobs)
       adt_out%diag(iobs)%Obs_Minus_Forecast    = obs_adt%adt(iobs) - hofx%values(iobs)
    enddo

    call adt_out%write_diag()
    call adt_out%write_geoval(var_abs_topo,geoval_adt)
    call adt_out%finalize()

  end subroutine ufo_adt_simobs

end module ufo_adt_mod
