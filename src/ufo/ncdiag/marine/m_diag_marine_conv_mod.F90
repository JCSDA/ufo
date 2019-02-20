module m_diag_marine_conv
  use fckit_log_module, only: fckit_log
  use ncd_kinds, only:  i_kind,r_single,r_kind,r_double
  use nc_diag_write_mod, only: nc_diag_init, nc_diag_metadata, nc_diag_write

  use nc_diag_read_mod, only: nc_diag_read_get_dim
  use nc_diag_read_mod, only: nc_diag_read_init, nc_diag_read_close
  use nc_diag_read_mod, only: nc_diag_read_get_var
  use nc_diag_read_mod, only: nc_diag_read_get_global_attr

  implicit none

  private
  save

  public :: diag_marine_conv_header
  public :: diag_marine_conv_tracer  !generic name for tao T,S obs - obs w/ single obs associated

  public :: write_split_marine_conv_diag_nc

  public :: read_marine_conv_diag_nc_header
  public :: read_marine_conv_diag_nc_tracer

  type diag_marine_conv_header
     character(3),    dimension(:), allocatable :: ObsType
     integer(i_kind)                            :: n_ObsType
     integer(i_kind)                            :: n_Observations_Tracer
     integer(i_kind)                            :: date
  end type diag_marine_conv_header

  type diag_marine_conv_tracer 
     integer      :: Station_ID
     real(r_kind) :: Observation_Type
     real(r_kind) :: Latitude
     real(r_kind) :: Longitude
     real(r_kind) :: Station_Depth
     real(r_kind) :: Time
     real(r_kind) :: Observation
     real(r_kind) :: Obs_Minus_Forecast_adjusted
     real(r_kind) :: Obs_Minus_Forecast_unadjusted
  end type diag_marine_conv_tracer

  integer,parameter :: maxobstype=30
  integer           :: nobstype
  integer           :: tao_tracer_type = 120 ! type for TAO T and S
  integer           :: t_qcmark       = 0   ! 0=tv; 1=tdry

contains

  subroutine write_split_marine_conv_diag_nc(infn,tao_header, tao_tracer, append_suffix)
    character(120),                                                 intent(in)    :: infn
    type(diag_marine_conv_header),                                         intent(in)    :: tao_header
    type(diag_marine_conv_tracer),dimension(tao_header%n_Observations_Tracer),intent(in)    :: tao_tracer
    logical,                                                        intent(in)    :: append_suffix

  end subroutine write_split_marine_conv_diag_nc

  subroutine read_marine_conv_diag_nc_header(infn,tao_header)
    character(len=*),       intent(in)    :: infn
    type(diag_marine_conv_header), intent(inout) :: tao_header

  end subroutine read_marine_conv_diag_nc_header

  subroutine read_marine_conv_diag_nc_tracer(infn, tao_tracer, ierr)
    use netcdf
    implicit none
    character(len=*),             intent(in)    :: infn
    type(diag_marine_conv_tracer),pointer, intent(inout) :: tao_tracer(:)
    integer,                      intent(out)   :: ierr
    character(20)   :: str, str2
    integer(i_kind) :: i, itype, fid, nobs
    real,    allocatable :: obs(:,:,:,:), lon(:), lat(:), depth(:), time(:)
    integer :: nlon, nlat, ndepth, ntime
    integer :: ilon, ilat, idepth, itime, cnt, varid

    call nc_diag_read_init(infn,fid)

    nlon = nc_diag_read_get_dim(fid,'lon')
    nlat = nc_diag_read_get_dim(fid,'lat')
    ndepth = nc_diag_read_get_dim(fid,'depth')
    ntime = nc_diag_read_get_dim(fid,'time')

    nobs=nlon*nlat*ndepth*ntime

    allocate(tao_tracer(nobs))

    !allocate(obs(ntime, ndepth, nlat, nlon))
    allocate(obs(nlon,nlat,ndepth,ntime))
    allocate(time(ntime),lon(nlon),lat(nlat),depth(ndepth))

    call nc_diag_read_get_var(fid,"lat",                          lat   )
    call nc_diag_read_get_var(fid,"lon",                          lon   )
    call nc_diag_read_get_var(fid,"depth",                        depth )
    call nc_diag_read_get_var(fid,"time",                         time  )
    !call nc_diag_read_get_var(fid,"T_20",                         obs   ) ! reading of 4d array not im[lemented yet
    call nc_diag_read_close(infn)

    call check( nf90_open(infn, NF90_NOWRITE, fid) )
    call check( nf90_inq_varid(fid, "T_20", varid) )
    call check( nf90_get_var(fid, varid, obs) )
    call check( nf90_close(fid) )
    cnt=1

    do ilon=1,nlon
       do ilat=1,nlat
          do idepth=1,ndepth
             tao_tracer(cnt)%Latitude      = lat(ilat)                 
             tao_tracer(cnt)%Longitude     = lon(ilon)
             tao_tracer(cnt)%Station_Depth = depth(idepth)
             tao_tracer(cnt)%Time          = time(1)
             !tao_tracer(cnt)%Observation   = obs(1,idepth,ilat,ilon)!obs(ilon,ilat,idepth,1)
             tao_tracer(cnt)%Observation   = obs(ilon,ilat,idepth,1)             
             cnt=cnt+1
          end do
       end do
    end do

    deallocate(lat,lon,depth,time,obs)

  end subroutine read_marine_conv_diag_nc_tracer

  subroutine check(status)
    use netcdf
    implicit none
    integer, intent ( in) :: status

    if(status /= nf90_noerr) then
       print *, trim(nf90_strerror(status))
       stop "Stopped"
    end if
  end subroutine check
end module m_diag_marine_conv
