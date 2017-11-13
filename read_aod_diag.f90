!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:    read_aoddiag                       read rad diag file
!   prgmmr: tahara           org: np20                date: 2003-01-01
!
! abstract:  This module contains code to process radiance
!            diagnostic files.  The module defines structures
!            to contain information from the radiance
!            diagnostic files and then provides two routines
!            to access contents of the file.
!
! program history log:
!   2005-07-22 treadon - add this doc block
!   2010-10-05 treadon - refactor code to GSI standard
!   2010-10-08 zhu     - use data_tmp to handle various npred values
!   2011-02-22 kleist  - changes related to memory allocate/deallocate
!   2011-04-08 li      - add tref, dtw, dtc to diag_data_fix_list, add tb_tz to diag_data_chan_list
!   2011-07-24 safford - make structure size for reading data_fix data version dependent 
!   2013-11-21 todling - revisit how versions are set (add set/get_aoddiag)
!   2014-01-27 todling - add ob sensitivity index
!   2017-07-13 mccarty - incorporate hooks for nc4/binary diag reading
!   2017-11-10 pagowski - converted radiance to aod 
!
! contains
!   read_aoddiag_header - read radiance diagnostic file header
!   read_aoddiag_data   - read radiance diagnostic file data
!   set_netcdf_read    - call set_netcdf_read(.true.) to use nc4 hooks, otherwise read file as 
!                           traditional binary format
!
! attributes:
!   language: f90
!   machine:  ibm RS/6000 SP
!
!$$$

MODULE read_aod_diag

  USE ncd_kinds, ONLY:  i_kind,r_single,r_kind
  USE nc_diag_read_mod, ONLY: nc_diag_read_get_var,  nc_diag_read_get_global_attr
  USE ncdr_dims, ONLY: nc_diag_read_get_dim
  IMPLICIT NONE

! Declare public and private
  PRIVATE

  PUBLIC :: diag_header_fix_list
  PUBLIC :: diag_header_chan_list
  PUBLIC :: diag_data_name_list
  PUBLIC :: diag_data_fix_list
  PUBLIC :: diag_data_chan_list
  PUBLIC :: read_aoddiag_header
  PUBLIC :: read_aoddiag_data
  PUBLIC :: set_netcdf_read
  PUBLIC :: ireal_aod
  PUBLIC :: ipchan_aod
  PUBLIC :: set_aoddiag
  PUBLIC :: get_aoddiag

  INTERFACE set_aoddiag
     MODULE PROCEDURE set_aoddiag_int_ ! internal procedure for integers
  END INTERFACE
  INTERFACE get_aoddiag
     MODULE PROCEDURE get_aoddiag_int_ ! internal procedure for integers
  END INTERFACE

  INTEGER(i_kind),PARAMETER :: ireal_aod  = 6   ! number of real entries per spot in aod diagnostic file
  INTEGER(i_kind),PARAMETER :: ipchan_aod = 4    ! number of entries per channel per spot in aod diagnostic file

!@for aod remove npred i/jextra some other terms
! Declare structures for radiance diagnostic file information
  TYPE diag_header_fix_list
     CHARACTER(len=20) :: isis           ! sat and sensor type
     CHARACTER(len=10) :: id             ! sat type
     CHARACTER(len=10) :: obstype        ! observation type
     INTEGER(i_kind) :: jiter            ! outer loop counter
     INTEGER(i_kind) :: nchan            ! number of channels in the sensor
     INTEGER(i_kind) :: idate            ! time (yyyymmddhh)
     INTEGER(i_kind) :: ireal            ! # of real elements in the fix part of a data record
     INTEGER(i_kind) :: ipchan           ! # of elements for each channel except for bias correction terms
     INTEGER(i_kind) :: nsig             ! # of sigma levels (layers?)
     INTEGER(i_kind) :: isens            ! sensitivity index
  END TYPE diag_header_fix_list

  TYPE diag_data_name_list
     CHARACTER(len=10),DIMENSION(ireal_aod) :: fix
     CHARACTER(len=10),DIMENSION(:),ALLOCATABLE :: chn
  END TYPE diag_data_name_list


!for aod diag_header_chan_list is same as for radiance  
  TYPE diag_header_chan_list
     REAL(r_single) :: freq              ! frequency (Hz)
     REAL(r_single) :: polar             ! polarization
     REAL(r_single) :: wave              ! wave number (cm^-1)
     REAL(r_single) :: varch             ! error variance (or SD error?)
     REAL(r_single) :: tlapmean          ! mean lapse rate
     INTEGER(i_kind):: iuse              ! use flag
     INTEGER(i_kind):: nuchan            ! sensor relative channel number
     INTEGER(i_kind):: iochan            ! satinfo relative channel number
  END TYPE diag_header_chan_list

!@some changes
  TYPE diag_data_fix_list
     REAL(r_single) :: lat               ! latitude (deg)
     REAL(r_single) :: lon               ! longitude (deg)
     REAL(r_single) :: psfc              ! psfc (hPa)
     REAL(r_single) :: obstime           ! observation time relative to analysis
     REAL(r_single) :: solzen_ang        ! solar zenith angle (deg)
     REAL(r_single) :: solazm_ang        ! solar azimumth angle (deg)
  END TYPE diag_data_fix_list

!@some changes to aod
  TYPE diag_data_chan_list
     REAL(r_single) :: aodobs             ! AOD (obs) 
     REAL(r_single) :: omgaod             ! AOD (obs) - AOD (guess)
     REAL(r_single) :: errinv             ! inverse error (K**(-1))
     REAL(r_single) :: qcmark             ! quality control mark
  END TYPE diag_data_chan_list

  REAL(r_single),PARAMETER::  rmiss_aoddiag    = -9.9e11_r_single

  LOGICAL,SAVE            ::  netcdf = .FALSE. 
  LOGICAL,SAVE            ::  nc_read = .FALSE.
  INTEGER,SAVE            ::  cur_ob_idx = -9999
  INTEGER,SAVE            ::  num_records = -9999

  TYPE(diag_data_fix_list)   ,ALLOCATABLE, SAVE :: all_data_fix(:)
  TYPE(diag_data_chan_list)  ,ALLOCATABLE, SAVE :: all_data_chan(:,:)

CONTAINS

  SUBROUTINE set_aoddiag_int_ (what,iv,ier)
    CHARACTER(len=*),INTENT(in) :: what
    INTEGER(i_kind),INTENT(in) :: iv
    INTEGER(i_kind),INTENT(out):: ier
    ier=0
  END SUBROUTINE set_aoddiag_int_

  SUBROUTINE get_aoddiag_int_ (what,iv,ier)
    CHARACTER(len=*),INTENT(in) :: what
    INTEGER(i_kind),INTENT(out):: iv
    INTEGER(i_kind),INTENT(out):: ier
    ier=0
  END SUBROUTINE get_aoddiag_int_

  SUBROUTINE set_netcdf_read(use_netcdf)
!                .      .    .                                       .
! subprogram:    read_diag_header_bin             read rad diag header
!   prgmmr: mccarty           org: gmao                date: 2015-08-06
!
! abstract:  This routine sets the routines to read from a netcdf file.
!            The default currently is to read binary files
!
! program history log:
!   2015-08-06 mccarty - created routine 
!
! input argument list:
!   use_netcdf    - logical .true. tells routine to read netcdf diag 
! attributes:
!   language: f90
!   machine:  ibm RS/6000 SP
!
!$$$
    LOGICAL,INTENT(in)                     :: use_netcdf

    netcdf = use_netcdf
  END SUBROUTINE set_netcdf_read


  SUBROUTINE read_aoddiag_header(ftin,header_fix,header_chan,data_name,iflag,lverbose)
!                .      .    .                                       .
! subprogram:    read_diag_header_bin             read rad diag header
!   prgmmr: mccarty           org: gmao                date: 2015-08-06
!
! abstract:  This routine reads the header record from a radiance
!            diagnostic file
!
! program history log:
!   2015-08-06 mccarty - created routine w/ fork for ncdiag or binary
!
! input argument list:
!   ftin          - unit number connected to diagnostic file 
!
! output argument list:
!   header_fix    - header information structure
!   header_chan   - channel information structure
!   data_name     - diag file data names
!   iflag         - error code
!   lverbose      - optional flag to turn off default output to standard out 
!
! attributes:
!   language: f90
!   machine:  ibm RS/6000 SP
!
!$$$

! Declare passed arguments
    INTEGER(i_kind),INTENT(in)             :: ftin
    TYPE(diag_header_fix_list ),INTENT(out):: header_fix
    TYPE(diag_header_chan_list),ALLOCATABLE :: header_chan(:)
    TYPE(diag_data_name_list)              :: data_name
    INTEGER(i_kind),INTENT(out)            :: iflag
    LOGICAL,OPTIONAL,INTENT(in)            :: lverbose

    iflag = 0
    IF (netcdf) THEN
       PRINT *,'netcdf slot'
       CALL read_aoddiag_header_nc(ftin,header_fix,header_chan,data_name,iflag,lverbose)
    ELSE
       CALL read_aoddiag_header_bin(ftin,header_fix,header_chan,data_name,iflag,lverbose)
    ENDIF

  END SUBROUTINE read_aoddiag_header

  SUBROUTINE read_aoddiag_header_nc(ftin,header_fix,header_chan,data_name,iflag,lverbose)

!nc not tested
!                .      .    .                                       .
! subprogram:    read_diag_header_nc              read rad diag header
!   prgmmr: mccarty           org: gmao                date: 2003-01-01
!
! abstract:  This routine reads the header record from a radiance
!            diagnostic file
!
! program history log:
!   2015-08-06 mccarty - Created routine for ncdiag header reading
!
! input argument list:
!   ftin          - unit number connected to diagnostic file 
!
! output argument list:
!   header_fix    - header information structure
!   header_chan   - channel information structure
!   data_name     - diag file data names
!   iflag         - error code
!   lverbose      - optional flag to turn off default output to standard out 
!
! attributes:
!   language: f90
!   machine:  ibm RS/6000 SP
!
!$$$
! Declare passed arguments
    INTEGER(i_kind),INTENT(in)             :: ftin
    TYPE(diag_header_fix_list ),INTENT(out):: header_fix
    TYPE(diag_header_chan_list),ALLOCATABLE :: header_chan(:)
    TYPE(diag_data_name_list)              :: data_name
    INTEGER(i_kind),INTENT(out)            :: iflag
    LOGICAL,OPTIONAL,INTENT(in)            :: lverbose

! local variables
    INTEGER(i_kind)                        :: nchan_dim
    REAL(r_kind),ALLOCATABLE,DIMENSION(:)               :: r_var_stor
    INTEGER(i_kind),ALLOCATABLE,DIMENSION(:)            :: i_var_stor
    CHARACTER(20)                          :: isis
    CHARACTER(10)                          :: id, obstype
!  integer(i_kind),dimension(:),allocatable           :: jiter, nchan_diag, idate, &
    INTEGER(i_kind)                        :: jiter, nchan_diag,  idate, &
         ireal, ipchan, isens

    iflag = 0
!  allocate(nchan_diag(1) )
    nchan_dim = nc_diag_read_get_dim(ftin,'nchans')
    header_fix%nchan = nchan_dim
    WRITE(*,*)'Number of channels=',nchan_dim

    CALL nc_diag_read_get_global_attr(ftin, "Number_of_channels", nchan_diag)

    IF (nchan_dim .NE. nchan_diag) THEN
       WRITE(*,*)'ERROR: Number of channels from dimension do not match those from header, aborting.'
       CALL abort
    ENDIF

    CALL nc_diag_read_get_global_attr(ftin, "Satellite_Sensor", isis)      ; header_fix%isis = isis
    CALL nc_diag_read_get_global_attr(ftin, "Satellite", id)               ; header_fix%id = id
    CALL nc_diag_read_get_global_attr(ftin, "Observation_type", obstype)   ; header_fix%obstype = obstype
    CALL nc_diag_read_get_global_attr(ftin, "Outer_Loop_Iteration", jiter) ; header_fix%jiter = jiter
    CALL nc_diag_read_get_global_attr(ftin, "date_time", idate)            ; header_fix%idate = idate
    CALL nc_diag_read_get_global_attr(ftin, "ireal_aod", ireal)         ; header_fix%ireal = ireal
    CALL nc_diag_read_get_global_attr(ftin, "ipchan_aod", ipchan)       ; header_fix%ipchan = ipchan
    CALL nc_diag_read_get_global_attr(ftin, "ioff0", isens)                ; header_fix%isens = isens


    ALLOCATE(header_chan(nchan_dim)   )

    ALLOCATE(r_var_stor(nchan_dim), &
         i_var_stor(nchan_dim)  )
    CALL nc_diag_read_get_var('frequency',r_var_stor)      ; header_chan%freq     = r_var_stor
    CALL nc_diag_read_get_var('polarization',i_var_stor)   ; header_chan%polar    = i_var_stor
    CALL nc_diag_read_get_var('wavenumber',r_var_stor)     ; header_chan%wave     = r_var_stor
    CALL nc_diag_read_get_var('error_variance',r_var_stor) ; header_chan%varch    = r_var_stor
    CALL nc_diag_read_get_var('use_flag',i_var_stor)       ; header_chan%iuse     = i_var_stor
    CALL nc_diag_read_get_var('sensor_chan',i_var_stor)    ; header_chan%nuchan   = i_var_stor
    CALL nc_diag_read_get_var('satinfo_chan',i_var_stor)   ; header_chan%iochan   = i_var_stor


  END SUBROUTINE read_aoddiag_header_nc
  
  SUBROUTINE read_aoddiag_header_bin(ftin,header_fix,header_chan,data_name,iflag,lverbose)
!                .      .    .                                       .
! subprogram:    read_diag_header_bin             read rad diag header
!   prgmmr: tahara           org: np20                date: 2003-01-01
!
! abstract:  This routine reads the header record from a radiance
!            diagnostic file
!
! program history log:
!   2010-10-05 treadon - add this doc block
!   2011-02-22 kleist  - changes related to memory allocation and standard output
!   2014-07-25 sienkiewicz - supress warning if npred_aoddiag == 0
!   2017-07-17 mccarty - renamed routine to _bin suffix for ncdiag
!
! input argument list:
!   ftin          - unit number connected to diagnostic file 
!
! output argument list:
!   header_fix    - header information structure
!   header_chan   - channel information structure
!   data_name     - diag file data names
!   iflag         - error code
!   lverbose      - optional flag to turn off default output to standard out 
!
! attributes:
!   language: f90
!   machine:  ibm RS/6000 SP
!
!$$$

! Declare passed arguments
    INTEGER(i_kind),INTENT(in)             :: ftin
    TYPE(diag_header_fix_list ),INTENT(out):: header_fix
    TYPE(diag_header_chan_list),ALLOCATABLE :: header_chan(:)
    TYPE(diag_data_name_list)              :: data_name
    INTEGER(i_kind),INTENT(out)            :: iflag
    LOGICAL,OPTIONAL,INTENT(in)            :: lverbose    

!  Declare local variables
    CHARACTER(len=2):: string
    CHARACTER(len=10):: satid,sentype
    CHARACTER(len=20):: sensat
    INTEGER(i_kind) :: i,ich
    INTEGER(i_kind):: jiter,nchanl,ianldate,ireal,ipchan,nsig,isens
    INTEGER(i_kind):: iuse_tmp,nuchan_tmp,iochan_tmp
    REAL(r_single) :: freq_tmp,polar_tmp,wave_tmp,varch_tmp,tlapmean_tmp
    LOGICAL loutall

    loutall=.TRUE.
    IF(PRESENT(lverbose)) loutall=lverbose

! Read header (fixed_part).
    READ(ftin,IOSTAT=iflag)  sensat,satid,sentype,jiter,nchanl,ianldate,&
         ireal,ipchan,nsig,isens

    IF (iflag/=0) THEN
       WRITE(6,*)'READ_AODDIAG_HEADER:  ***ERROR*** Unknown file format.  Cannot read'
       RETURN
    ENDIF

    header_fix%isis    = sensat
    header_fix%id      = satid
    header_fix%obstype = sentype
    header_fix%jiter   = jiter
    header_fix%nchan   = nchanl
    header_fix%idate   = ianldate
    header_fix%ireal   = ireal
    header_fix%ipchan  = ipchan
    header_fix%nsig   = nsig
    header_fix%isens   = isens

    IF (loutall) THEN
       WRITE(6,*)'READ_AODDIAG_HEADER:  isis=',header_fix%isis,&
            ' nchan=',header_fix%nchan,&
            ' isens=',header_fix%isens
    ENDIF

! Allocate and initialize as needed
    IF (ALLOCATED(header_chan)) DEALLOCATE(header_chan)
    IF (ALLOCATED(data_name%chn))  DEALLOCATE(data_name%chn)

    ALLOCATE(header_chan( header_fix%nchan))
    ALLOCATE(data_name%chn(header_fix%ipchan))

    data_name%fix(1) ='lat       '
    data_name%fix(2) ='lon       '
    data_name%fix(3) ='psfc      '
    data_name%fix(4) ='obstim    '
    data_name%fix(5) ='solzen    '
    data_name%fix(6) ='solazm    '
    data_name%chn(1)='obs       '
    data_name%chn(2)='omg       '
    data_name%chn(3)='errinv    '
    data_name%chn(4)='qcmark    '

! Read header (channel part)
    DO ich=1, header_fix%nchan
       READ(ftin,IOSTAT=iflag) freq_tmp,polar_tmp,wave_tmp,varch_tmp,iuse_tmp,nuchan_tmp,iochan_tmp

       header_chan(ich)%freq     = freq_tmp
       header_chan(ich)%polar    = polar_tmp
       header_chan(ich)%wave     = wave_tmp
       header_chan(ich)%varch    = varch_tmp
       header_chan(ich)%iuse     = iuse_tmp
       header_chan(ich)%nuchan   = nuchan_tmp
       header_chan(ich)%iochan   = iochan_tmp
       IF (iflag/=0) RETURN
    END DO

! Construct array containing menonics for data record entries

  END SUBROUTINE read_aoddiag_header_bin

  SUBROUTINE read_aoddiag_data(ftin,header_fix,data_fix,data_chan,iflag )
!                .      .    .                                       .
! subprogram:    read_aoddiag_dat                    read rad diag data
!   prgmmr: tahara           org: np20                date: 2003-01-01
!
! abstract:  This routine reads the data record from a radiance
!            diagnostic file
!
! program history log:
!   2010-10-05 treadon - add this doc block
!   2011-02-22 kleist  - changes related to memory allocation
!   2017-07-17 mccarty - change routine to be generalized for bin/nc4 files
!
! input argument list:
!   ftin - unit number connected to diagnostic file
!   header_fix - header information structure
!
! output argument list:
!   data_fix   - spot header information structure
!   data_chan  - spot channel information structure
!   iflag      - error code
!
! attributes:
!   language: f90
!   machine:  ibm RS/6000 SP
!
!$$$


! Declare passed arguments
    INTEGER(i_kind),INTENT(in)             :: ftin
    TYPE(diag_header_fix_list ),INTENT(in) :: header_fix
    TYPE(diag_data_fix_list)   ,INTENT(out):: data_fix
    TYPE(diag_data_chan_list)  ,ALLOCATABLE :: data_chan(:)
    INTEGER(i_kind),INTENT(out)            :: iflag

    IF (netcdf) THEN
       IF (.NOT. nc_read) CALL read_aoddiag_data_nc_init(ftin, header_fix)

       IF (cur_ob_idx .EQ. num_records ) THEN
          iflag = 0
       ELSE IF (cur_ob_idx .GT. num_records) THEN
          iflag = -1
       ELSE
          iflag = 1
       ENDIF

       IF (iflag .GE. 0) CALL read_aoddiag_data_nc(ftin,header_fix,data_fix,data_chan,iflag)

    ELSE
       CALL read_aoddiag_data_bin(ftin,header_fix,data_fix,data_chan,iflag )
    ENDIF

  END SUBROUTINE read_aoddiag_data

  SUBROUTINE read_aoddiag_data_nc_init(ftin, header_fix)
!                .      .    .                                       .
! subprogram:    read_aoddiag_data_nc_init           read rad diag data
!   prgmmr: mccarty          org: np20                date: 2015-08-10
!
! abstract:  This routine reads the data record from a netcdf radiance
!            diagnostic file
!
! program history log:
!   2015-06-10 mccarty  - create routine
!
! input argument list:
!   ftin - unit number connected to diagnostic file
!   header_fix - header information structure
!
! output argument list:
!   data_fix   - spot header information structure
!   data_chan  - spot channel information structure
!   iflag      - error code
!
! attributes:
!   language: f90
!   machine:  ibm RS/6000 SP
!
!$$$

! Declare passed arguments
    INTEGER(i_kind),INTENT(in)             :: ftin
    TYPE(diag_header_fix_list ),INTENT(in) :: header_fix

! Declare local variables
    INTEGER(i_kind)                          :: nrecord, ndatum, nangord
    INTEGER(i_kind)                          :: cch, ic, ir, cdatum
    REAL(r_kind), ALLOCATABLE, DIMENSION(:)  :: Latitude, Longitude, &
         &Obs_Time,  Psfc, Sol_Zenith_Angle, Sol_Azimuth_Angle,&
         &Observation, Obs_Minus_Forecast_unadjusted,  &
         &Inverse_Observation_Error, QC_Flag

    INTEGER(i_kind), ALLOCATABLE, DIMENSION(:)  :: Channel_Index

    REAL(r_kind)                                :: clat, clon

    ndatum = nc_diag_read_get_dim(ftin,'nobs')
    nrecord = ndatum / header_fix%nchan
    num_records = nrecord

    WRITE(*,*)'Reading ndatum, nrecord=',ndatum,nrecord

    ALLOCATE( &
         &Channel_Index(ndatum), Latitude(ndatum), Longitude(ndatum), &
         &Obs_Time(ndatum), Psfc(ndatum),&
         &Sol_Zenith_Angle(ndatum), Sol_Azimuth_Angle(ndatum),&
         &Observation(ndatum), Obs_Minus_Forecast_unadjusted(ndatum),&
         &Inverse_Observation_Error(ndatum),QC_Flag(ndatum))

    ALLOCATE(  all_data_fix(nrecord)        )
    ALLOCATE(  all_data_chan(nrecord, header_fix%nchan))

    CALL nc_diag_read_get_var('Channel_Index', Channel_Index)
    CALL nc_diag_read_get_var('Latitude', Latitude)
    CALL nc_diag_read_get_var('Longitude', Longitude)
    CALL nc_diag_read_get_var('Psfc', Psfc)
    CALL nc_diag_read_get_var('Obs_Time', Obs_Time)
    CALL nc_diag_read_get_var('Sol_Zenith_Angle', Sol_Zenith_Angle)
    CALL nc_diag_read_get_var('Sol_Azimuth_Angle', Sol_Azimuth_Angle)
    CALL nc_diag_read_get_var('Observation', Observation)
    CALL nc_diag_read_get_var('Obs_Minus_Forecast_unadjusted', Obs_Minus_Forecast_unadjusted)
    CALL nc_diag_read_get_var('Inverse_Observation_Error', Inverse_Observation_Error)
    CALL nc_diag_read_get_var('QC_Flag', QC_Flag)
    cdatum = 1

!  allocate(  all_data_fix(nrecord)        )
!  allocate(  all_data_chan(nrecord, nchan))


    DO ir=1,nrecord
       clat = Latitude(cdatum)
       clon = Longitude(cdatum)
       all_data_fix(ir)%lat               = Latitude(cdatum)
       all_data_fix(ir)%lon               = Longitude(cdatum)
       all_data_fix(ir)%psfc              = psfc(cdatum)
       all_data_fix(ir)%obstime           = Obs_Time(cdatum)
       all_data_fix(ir)%solzen_ang        = Sol_Zenith_Angle(cdatum)
       all_data_fix(ir)%solazm_ang        = Sol_Azimuth_Angle(cdatum)

       DO ic=1,header_fix%nchan
          IF (clat .NE. Latitude(cdatum) .OR. clon .NE. Longitude(cdatum)) THEN
             WRITE(*,*) 'ERROR: Lats & Lons are mismatched.  This is bad'
             PRINT *,'irecord=',ir
             PRINT *,'clat,clon=',clat,clon
             PRINT *,'lat/lon(datum)=',Latitude(cdatum), Longitude(cdatum)
             CALL abort
          ENDIF
          cch = Channel_Index(cdatum)
          all_data_chan(ir,cch)%aodobs = Observation(cdatum)
          all_data_chan(ir,cch)%omgaod = Obs_Minus_Forecast_unadjusted(cdatum)
          all_data_chan(ir,cch)%errinv= Inverse_Observation_Error(cdatum)
          all_data_chan(ir,cch)%qcmark= QC_Flag(cdatum)

          cdatum = cdatum + 1
       ENDDO
    ENDDO

    nc_read = .TRUE.
    cur_ob_idx = 1
  END SUBROUTINE read_aoddiag_data_nc_init

  SUBROUTINE read_aoddiag_data_nc(ftin,header_fix,data_fix,data_chan,iflag )
!                .      .    .                                       .
! subprogram:    read_aoddiag_dat                    read rad diag data
!   prgmmr: tahara           org: np20                date: 2015-08-10
!
! abstract:  This routine reads the data record from a netcdf radiance
!            diagnostic file
!
! program history log:
!   2015-08-10 mccarty  - create routine
!
! input argument list:
!   ftin - unit number connected to diagnostic file
!   header_fix - header information structure
!
! output argument list:
!   data_fix   - spot header information structure
!   data_chan  - spot channel information structure
!   iflag      - error code
!
! attributes:
!   language: f90
!   machine:  ibm RS/6000 SP
!
!$$$

! Declare passed arguments
    INTEGER(i_kind),INTENT(in)             :: ftin
    TYPE(diag_header_fix_list ),INTENT(in) :: header_fix
    TYPE(diag_data_fix_list)   ,INTENT(out):: data_fix
    TYPE(diag_data_chan_list)  ,ALLOCATABLE :: data_chan(:)
    INTEGER(i_kind),INTENT(out)            :: iflag

    iflag = 0
    IF (.NOT. ALLOCATED(data_chan)) ALLOCATE(data_chan(header_fix%nchan) )

    data_fix     = all_data_fix(cur_ob_idx)
    data_chan(:) = all_data_chan(cur_ob_idx,:) 

    cur_ob_idx = cur_ob_idx + 1

  END SUBROUTINE read_aoddiag_data_nc

  SUBROUTINE read_aoddiag_data_bin(ftin,header_fix,data_fix,data_chan,iflag )
!                .      .    .                                       .
! subprogram:    read_aoddiag_dat                    read rad diag data
!   prgmmr: tahara           org: np20                date: 2003-01-01
!
! abstract:  This routine reads the data record from a radiance
!            diagnostic file
!
! program history log:
!   2010-10-05 treadon - add this doc block
!   2011-02-22 kleist  - changes related to memory allocation
!   2017-07-17 mccarty - rename binary-specific procedure
!
! input argument list:
!   ftin - unit number connected to diagnostic file
!   header_fix - header information structure
!
! output argument list:
!   data_fix   - spot header information structure
!   data_chan  - spot channel information structure
!   iflag      - error code
!
! attributes:
!   language: f90
!   machine:  ibm RS/6000 SP
!
!$$$


! Declare passed arguments
    INTEGER(i_kind),INTENT(in)             :: ftin
    TYPE(diag_header_fix_list ),INTENT(in) :: header_fix
    TYPE(diag_data_fix_list)   ,INTENT(out):: data_fix
    TYPE(diag_data_chan_list)  ,ALLOCATABLE :: data_chan(:)
    INTEGER(i_kind),INTENT(out)            :: iflag

    INTEGER(i_kind) :: ich,iang,i,j
    REAL(r_single),DIMENSION(:,:),ALLOCATABLE :: data_tmp
    REAL(r_single),DIMENSION(:),ALLOCATABLE   :: fix_tmp

! Allocate arrays as needed
    IF (ALLOCATED(data_chan)) DEALLOCATE(data_chan)
    ALLOCATE(data_chan(header_fix%nchan))

! Allocate arrays to hold data record
    ALLOCATE(data_tmp(header_fix%ipchan,header_fix%nchan))

    ALLOCATE( fix_tmp( ireal_aod ) )

! Read data record

    READ(ftin,IOSTAT=iflag) fix_tmp, data_tmp

! Transfer fix_tmp record to output structure
    data_fix%lat = fix_tmp(1)
    data_fix%lon = fix_tmp(2)
    data_fix%psfc = fix_tmp(3)
    data_fix%obstime = fix_tmp(4) 
    data_fix%solzen_ang = fix_tmp(5)
    data_fix%solazm_ang = fix_tmp(6)

! Transfer data record to output structure
    DO ich=1,header_fix%nchan
       data_chan(ich)%aodobs =data_tmp(1,ich)
       data_chan(ich)%omgaod =data_tmp(2,ich)
       data_chan(ich)%errinv=data_tmp(3,ich)
       data_chan(ich)%qcmark=data_tmp(4,ich)
    ENDDO

    DEALLOCATE(data_tmp, fix_tmp)

  END SUBROUTINE read_aoddiag_data_bin

END MODULE read_aod_diag

