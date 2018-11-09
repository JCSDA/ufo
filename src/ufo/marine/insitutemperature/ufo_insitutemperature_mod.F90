! (C) Copyright 2017-2018 UCAR
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 

!> Fortran module to handle temperature profile observations

module ufo_insitutemperature_mod

  use iso_c_binding
  use obsspace_mod
  use ufo_vars_mod
  use ioda_locs_mod
  use ufo_geovals_mod
  use kinds

  implicit none
  public :: ufo_insitutemperature
  public :: ufo_insitutemperature_simobs
  private
  integer, parameter :: max_string=800

  !> Fortran derived type for insitu temperature profile observation operator
  type :: ufo_insitutemperature
  end type ufo_insitutemperature

  ! ------------------------------------------------------------------------------

contains

  ! ------------------------------------------------------------------------------

  subroutine ufo_insitutemperature_simobs(self, geovals, hofx, obss)
    use gsw_pot_to_insitu
    use vert_interp_mod    
    use ufo_tpsp2ti_mod
    use ufo_marine_ncutils
    
    implicit none
    type(ufo_insitutemperature), intent(in)  :: self       !< Trajectory
    type(ufo_geovals), intent(in)            :: geovals    !< Model's Tp, Sp, h interpolated at obs location 
    type(c_ptr), value, intent(in)           :: obss       !< Insitu temperature observations
    real(c_double),  intent(inout)           :: hofx(:)    !< Ti(Tp,Sp,h)

    character(len=*), parameter :: myname_="ufo_insitutemperature_simobs"
    character(max_string)  :: err_msg

    integer :: iobs, ilev, nlev, nobs
    type(ufo_geoval), pointer :: temp, salt, h
    real (kind_real), allocatable :: depth(:,:)
    real(kind_real) :: lono, lato, deptho
    real(kind_real), allocatable :: obs_lon(:)
    real(kind_real), allocatable :: obs_lat(:)
    real(kind_real), allocatable :: obs_depth(:)
    real(kind_real), allocatable :: obs_val(:)
    integer :: obss_nobs
    
    ! Vertical interpolation
    real(kind_real) :: wf, tp, sp, prs
    integer :: wi
    
    ! netcdf stuff
    character(len=120) :: filename !< name of outpu file for omf, lon, lat, ...
    type(diag_marine_obs) :: insitu_out    

    ! check if nobs is consistent in geovals & hofx
    if (geovals%nobs /= size(hofx,1)) then
       write(err_msg,*) myname_, ' error: nobs inconsistent!'
       call abor1_ftn(err_msg)
    endif

    ! check if sea temperature profile variable is in geovals and get it
    call ufo_geovals_get_var(geovals, var_ocn_pot_temp, temp)
    
    ! check if sea salinity profile variable is in geovals and get it
    call ufo_geovals_get_var(geovals, var_ocn_salt, salt)

    ! check if ocean layer thickness variable is in geovals and get it
    call ufo_geovals_get_var(geovals, var_ocn_lay_thick, h)

    ! Read in obs data
    obss_nobs = obsspace_get_nobs(obss)
    allocate(obs_lon(obss_nobs))
    allocate(obs_lat(obss_nobs))
    allocate(obs_depth(obss_nobs))
    allocate(obs_val(obss_nobs))

    call obsspace_get_db(obss, "", "longitude", obs_lon)
    call obsspace_get_db(obss, "", "latitude", obs_lat)
    call obsspace_get_db(obss, "", "ocean_depth", obs_depth)
    call obsspace_get_db(obss, "", "insitu_temperature", obs_val)

    nlev = temp%nval
    nobs = temp%nobs        
    allocate(depth(nlev,nobs))
    do iobs = 1,size(hofx,1)
       !< Depth from layer thickness
       depth(1,iobs)=0.5*h%vals(1,iobs)
       do ilev = 2, nlev
          depth(ilev,iobs)=sum(h%vals(1:ilev-1,iobs))+0.5*h%vals(ilev,iobs)
       end do          
    end do

    ! Information for temporary output file
    
    hofx = 0.0
    ! insitu temperature profile obs operator
    do iobs = 1,size(hofx,1)

       lono = obs_lon(iobs)
       lato = obs_lat(iobs)
       deptho = obs_depth(iobs)
    
       !< Interpolation weight
       call vert_interp_weights(nlev, deptho, depth(:,iobs), wi, wf)
       if (deptho.ge.maxval(depth)) then
          wi=nlev-1
          wf=0.0
       end if

       ! Interpolate temp_p, salt_p to deptho
       call vert_interp_apply(nlev, temp%vals(:,iobs), tp, wi, wf)
       call vert_interp_apply(nlev, salt%vals(:,iobs), sp, wi, wf)

       ! Get insitu temp at model levels and obs location (lono, lato, zo)
       call insitu_t_nl(hofx(iobs), tp, sp, lono, lato, deptho)

    enddo

    deallocate(depth)
    deallocate(obs_lon)
    deallocate(obs_lat)
    deallocate(obs_depth)
    deallocate(obs_val)
    
  end subroutine ufo_insitutemperature_simobs

end module ufo_insitutemperature_mod
