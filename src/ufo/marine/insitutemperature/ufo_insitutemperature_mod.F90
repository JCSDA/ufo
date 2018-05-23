! (C) Copyright 2017 UCAR
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 

!> Fortran module to handle temperature profile observations

module ufo_insitutemperature_mod

  use ioda_obs_insitutemperature_mod
  use ioda_obs_vectors
  use ufo_vars_mod
  use ioda_locs_mod
  use ufo_geovals_mod
  use kinds

  implicit none
  public :: ufo_insitutemperature
  public :: ufo_insitutemperature_eqv
  private
  integer, parameter :: max_string=800

  !> Fortran derived type for temperature profile observation operator
  type :: ufo_insitutemperature
  end type ufo_insitutemperature

  ! ------------------------------------------------------------------------------

contains

  ! ------------------------------------------------------------------------------

  subroutine ufo_insitutemperature_eqv(self, geovals, hofx, obs_ti)
    use m_diag_marine_conv
    use nc_diag_write_mod, only: nc_diag_init, nc_diag_metadata, nc_diag_write
    use gsw_pot_to_insitu
    !use insitu_temperature_mod
    use vert_interp_mod    
    use ufo_tpsp2ti_mod
    
    implicit none
    type(ufo_insitutemperature), intent(in)     :: self       !< Trajectory !!!! REMOVE, NOT USED!!!! 
    type(ufo_geovals), intent(in)                :: geovals    !< Model's Tp, Sp, h interpolated at obs location 
    type(ioda_obs_insitutemperature), intent(in) :: obs_ti     !< Insitu temperature observations
    type(obs_vector),  intent(inout)             :: hofx       !< Ti(Tp,Sp,h)

    character(len=*), parameter :: myname_="ufo_insitutemperature_eqv"
    character(max_string)  :: err_msg

    integer :: iobs, ilev, nlev, nobs
    type(ufo_geoval), pointer :: temp, salt, h
    real(kind_real), allocatable :: tempi(:,:)
    real (kind_real), allocatable :: pressure(:,:), depth(:,:)
    real(kind_real) :: lono, lato, deptho
    
    ! Vertical interpolation
    real(kind_real) :: wf, tp, sp, prs
    integer :: wi
    
    ! nc_diag stuff
    logical :: append
    character(len=120) :: filename !< name of outpu file for omf, lon, lat, ...
    type(diag_marine_conv_tracer), allocatable :: Argo(:)

    ! check if nobs is consistent in geovals & hofx
    if (geovals%nobs /= hofx%nobs) then
       write(err_msg,*) myname_, ' error: nobs inconsistent!'
       call abor1_ftn(err_msg)
    endif

    ! check if sea temperature profile variable is in geovals and get it
    if (.not. ufo_geovals_get_var(geovals, var_ocn_pot_temp, temp)) then
       write(err_msg,*) myname_, trim(var_ocn_pot_temp), ' does not exist'
       call abor1_ftn(err_msg)
    endif

    ! check if sea salinity profile variable is in geovals and get it
    if (.not. ufo_geovals_get_var(geovals, var_ocn_salt, salt)) then
       write(err_msg,*) myname_, trim(var_ocn_salt), ' does not exist'
       call abor1_ftn(err_msg)
    endif

    ! check if ocean layer thickness variable is in geovals and get it
    if (.not. ufo_geovals_get_var(geovals, var_ocn_lay_thick, h)) then
       write(err_msg,*) myname_, trim(var_ocn_lay_thick), ' does not exist'
       call abor1_ftn(err_msg)
    endif

    nlev = temp%nval
    nobs = temp%nobs        
    allocate(tempi(nlev,nobs))
    allocate(pressure(nlev,nobs), depth(nlev,nobs))
    do iobs = 1,hofx%nobs
       !< Depth from layer thickness
       depth(1,iobs)=0.5*h%vals(1,iobs)
       do ilev = 2, nlev
          depth(ilev,iobs)=sum(h%vals(1:ilev-1,iobs))+0.5*h%vals(ilev,iobs)
       end do          
    end do

    do iobs = 1,hofx%nobs
       do ilev = 1, nlev
          lono = obs_ti%lon(iobs)
          lato = obs_ti%lat(iobs)          
          call insitu_t_nl(tempi(ilev,iobs),temp%vals(ilev,iobs),salt%vals(ilev,iobs),lono,lato,depth(ilev,iobs))
       end do
    end do

    allocate(Argo(hofx%nobs))    
    hofx%values = 0.0
    ! insitu temperature profile obs operator
    do iobs = 1,hofx%nobs

       lono = obs_ti%lon(iobs)
       lato = obs_ti%lat(iobs)
       deptho = obs_ti%depth(iobs)
    
       !< Interpolation weight
       call vert_interp_weights(nlev,deptho,depth(:,iobs),wi,wf)
       if (deptho.ge.maxval(depth)) then
          wi=nlev-1
          wf=0.0
       end if

       ! Interpolate temp_p, salt_p to deptho
       call vert_interp_apply(nlev, temp%vals(:,iobs), tp, wi, wf)
       call vert_interp_apply(nlev, salt%vals(:,iobs), sp, wi, wf)

       ! Get insitu temp at model levels and obs location (lono, lato, zo)
       call insitu_t_nl(hofx%values(iobs),tp,sp,lono,lato,deptho)

       if (isnan(hofx%values(iobs))) then !!!!!! HACK !!!!!!!!!!!!!!!!!!!!!
          hofx%values(iobs)=0.0 !!!! NEED TO QC OUT BAD OBS LOCATION !!!!!!
       end if
       
       Argo(iobs)%Station_ID                     = obs_ti%idx(iobs)
       Argo(iobs)%Observation_Type               = 1.0
       Argo(iobs)%Latitude                       = obs_ti%lat(iobs)
       Argo(iobs)%Longitude                      = obs_ti%lon(iobs)
       Argo(iobs)%Station_Depth                  = obs_ti%depth(iobs)
       Argo(iobs)%Time                           = 1.0
       Argo(iobs)%Observation                    = obs_ti%val(iobs)
       Argo(iobs)%Obs_Minus_Forecast_adjusted    = obs_ti%val(iobs) - hofx%values(iobs) 
       Argo(iobs)%Obs_Minus_Forecast_unadjusted  = obs_ti%val(iobs) - hofx%values(iobs) 

       write(402,*)hofx%values(iobs)
       
    enddo

    filename='test-obs-geovals.nc'
    call tonc(filename,Argo,tempi,temp,salt,h)
    deallocate(Argo)

  end subroutine ufo_insitutemperature_eqv

  subroutine tonc(filename,argo,tempi,temp,salt,h)
    use netcdf
    use m_diag_marine_conv
    implicit none

    character(len=120), intent(in)           :: filename !< name of outpu file for omf, lon, lat, ...
    type(diag_marine_conv_tracer), intent(in)       :: Argo(:)  !< Argo obs
    real(kind_real), allocatable, intent(in) :: tempi(:,:)
    type(ufo_geoval), pointer, intent(in)    :: temp, salt, h

    integer :: iobs

    !netcdf stuff
    integer(kind=4) :: iNcid,i,iStationNo
    integer(kind=4) :: iDimStation_ID,iDimTime_ID,iDimLenStringName_ID, iDimLev_ID
    integer(kind=4) :: iVarLON_ID, iVarLAT_ID, iVarLev_ID, iVarObs_ID
    integer(kind=4) :: iVarOMF_ID, iVarOMA_ID, iVarSIGO_ID, iVarOBSID_ID
    integer(kind=4) :: iVartemp_ID, iVarsalt_ID, iVarh_ID, iVartempi_ID    

    integer, allocatable :: obsid(:)
    integer :: nlev,nobs

    ! Create file.
    call check(nf90_create(filename, NF90_CLOBBER, iNcid))

    ! Define the dimensions. The Station-record dimension is defined to have
    call check(nf90_def_dim(iNcid, "nobs", NF90_UNLIMITED, iDimStation_ID))
    nlev = temp%nval
    nobs = temp%nobs    
    call check(nf90_def_dim(iNcid, "nlev", nlev, iDimLev_ID)) !< Number of model levels

    ! Define of variables.
    ! Obs space
    call check( nf90_def_var(iNcid, "OBSID", NF90_INT, (/ iDimStation_ID /), iVarOBSID_ID) )    
    call check( nf90_def_var(iNcid, "lon", NF90_REAL,  (/ iDimStation_ID /), iVarLON_ID) )
    call check( nf90_def_var(iNcid, "lat", NF90_REAL,  (/ iDimStation_ID /), iVarLAT_ID) )
    call check( nf90_def_var(iNcid, "lev", NF90_REAL,  (/ iDimStation_ID /), iVarLev_ID) )
    call check( nf90_def_var(iNcid, "obs", NF90_REAL,  (/ iDimStation_ID /), iVarObs_ID) )
    call check( nf90_def_var(iNcid, "sigo", NF90_REAL, (/ iDimStation_ID /), iVarSIGO_ID) )    
    call check( nf90_def_var(iNcid, "omf", NF90_REAL,  (/ iDimStation_ID /), iVarOMF_ID) )
    call check( nf90_def_var(iNcid, "oma", NF90_REAL,  (/ iDimStation_ID /), iVarOMA_ID) )
    ! At model levels
    call check( nf90_def_var(iNcid, "ocean_potential_temperature", NF90_REAL,  (/ iDimLev_ID, iDimStation_ID /), iVartemp_ID) )
    call check( nf90_def_var(iNcid, "ocean_salinity", NF90_REAL,  (/ iDimLev_ID, iDimStation_ID /), iVarsalt_ID) )
    call check( nf90_def_var(iNcid, "ocean_layer_thickness", NF90_REAL,  (/ iDimLev_ID, iDimStation_ID /), iVarh_ID) )
    call check( nf90_def_var(iNcid, "tempi", NF90_REAL,  (/ iDimLev_ID, iDimStation_ID /), iVartempi_ID) )    

    ! End define mode.
    call check(nf90_enddef(iNcid))

!!$    Argo(iobs)%Station_ID                     = "argo"
!!$    Argo(iobs)%Observation_Type               = 1.0
!!$    Argo(iobs)%Latitude                       = obs_ti%lon(iobs)
!!$    Argo(iobs)%Longitude                      = obs_ti%lat(iobs)
!!$    Argo(iobs)%Station_Depth                  = obs_ti%depth(iobs)
!!$    Argo(iobs)%Time                           = 1.0
!!$    Argo(iobs)%Observation                    = obs_ti%val(iobs)
!!$    Argo(iobs)%Obs_Minus_Forecast_adjusted    = obs_ti%val(iobs) - hofx%values(iobs) 
!!$    Argo(iobs)%Obs_Minus_Forecast_unadjusted  = obs_ti%val(iobs) - hofx%values(iobs) 

    ! Writing
    allocate(obsid(size(Argo(:),1)))
    obsid=3073
    call check(nf90_put_var(iNcid, iVarLON_ID , obsid(:)))!Argo(:)%Station_ID))
    call check(nf90_put_var(iNcid, iVarLON_ID , Argo(:)%Longitude))
    call check(nf90_put_var(iNcid, iVarLAT_ID , Argo(:)%Latitude))    
    call check(nf90_put_var(iNcid, iVarLev_ID , Argo(:)%Station_Depth))
    call check(nf90_put_var(iNcid, iVarObs_ID , Argo(:)%Observation))
    call check(nf90_put_var(iNcid, iVarOMF_ID , Argo(:)%Obs_Minus_Forecast_adjusted))
    call check(nf90_put_var(iNcid, iVarOMA_ID , Argo(:)%Obs_Minus_Forecast_adjusted))
    call check(nf90_put_var(iNcid, iVartemp_ID , temp%vals))
    call check(nf90_put_var(iNcid, iVarsalt_ID , salt%vals))
    call check(nf90_put_var(iNcid, iVarh_ID , h%vals))
    !call check(nf90_put_var(iNcid, iVarSIGO_ID , Argo%sigo))    
    call check(nf90_put_var(iNcid, iVartempi_ID , tempi))

    ! Close file.
    call check(nf90_close(iNcid))

  end subroutine tonc

  subroutine check(status)
    use netcdf
    IMPLICIT NONE
    integer(4), intent ( in) :: status
    if(status /= nf90_noerr) then
       print *, trim(nf90_strerror(status))
       stop "Stopped"
    end if
  end subroutine check
end module ufo_insitutemperature_mod
