PROGRAM convert_aod_diag

  USE read_aod_diag,ONLY:read_aoddiag_header, diag_header_fix_list, diag_header_chan_list, diag_data_name_list
  USE read_aod_diag,ONLY:read_aoddiag_data, diag_data_fix_list, diag_data_chan_list
  USE ncd_kinds,ONLY: r_quad, r_single

  USE nc_diag_write_mod, ONLY: nc_diag_init, nc_diag_header, nc_diag_metadata, &
       nc_diag_write, nc_diag_data2d, nc_diag_chaninfo_dim_set, nc_diag_chaninfo

  IMPLICIT NONE 

  REAL,PARAMETER::     missing = -9.99e9
  INTEGER,PARAMETER:: imissing = -999999

  INTEGER nargs, iargc, n
  CHARACTER*256, ALLOCATABLE ::   arg(:)

  TYPE(diag_header_fix_list )          ::  header_fix
  TYPE(diag_header_chan_list),ALLOCATABLE  ::  header_chan(:)
  TYPE(diag_data_name_list)            ::  headname

  TYPE(diag_data_fix_list)             ::  data_fix
  TYPE(diag_data_chan_list)  ,ALLOCATABLE  ::  data_chan(:)


  INTEGER i
!!!  real(r_quad)                         ::  ret_var
!!!  real(r_quad)                         ::  ret_stddev

! commandline variables
  LOGICAL                              ::  debug 
  LOGICAL                              ::  append_suffix 

  CHARACTER*256 infn, outfn
  LOGICAL linfile, loutfile

  INTEGER,PARAMETER                    ::  inlun = 51
  INTEGER,PARAMETER                    ::  outlun= 52
  INTEGER,PARAMETER                    ::  nllun = 53

  INTEGER strlen, iflag
  INTEGER iuse, ich, nch, ipr

  LOGICAL,DIMENSION(:),ALLOCATABLE     :: luse
  LOGICAL lqcpass

  INTEGER,PARAMETER                    ::  nvar = 4 ! number of positions in array needed for inline variance calc
! variables for output, all to be allocated as nchan.  Variances will be calculated using inline algorithm
! accredited to Welford, according
!!!  real(r_quad),dimension(:),allocatable   :: nobstotal, nobsassim, tbtotal, tbassim, omf_nbc , omf_bc , sigo, jo
!!!  real(r_quad),dimension(:,:),allocatable ::                                         vomf_nbc, vomf_bc
! total bias and fixed bias terms.  When Yanqui's variational angle correction is brought in, this may need to be updated.
!!!  real(r_quad),dimension(:),allocatable   :: totbias , fixbias
!!!  real(r_quad),dimension(:,:),allocatable :: vtotbias, vfixbias
! variational bias correction variables, which will be allocated as nchan and npred_read
!!!  real(r_quad),dimension(:,:),allocatable   :: biasterms
!!!  real(r_quad),dimension(:,:,:),allocatable :: vbiasterms
! Definitions for above variables - 
!   nobstotal        - number of observations considered - total
!   nobsassim        - number of observations that are assimilated
!   tbtotal          - mean brightness temperature of all observations
!   tbassim          - mean brightness temperature of assimilated observations
!   omf_nbc/vomf_nbc - mean/variance of O-F without bias correction applied  
!   omf_bc/vomf_bc   - mean/variance of O-F with bias correction applied
!   sigo             - mean observation error of assimilated observations
!   jo               - mean cost (Jo) of assimilated observations
!   totbias/vtotbias - mean/variance of bias correction applied to assimilated observations
!   fixbias/vfixbias - mean/variance of fixed scan angle position bias correction (sac.x-derived)
!   biasterms
!        /vbiasterms - means/variances of the variational terms as defined by npred_read.  

! single variables used later for printing purposes
  INTEGER        :: inobstotal, inobsassim
  CHARACTER(len=13),DIMENSION(:),ALLOCATABLE :: chfrwn 

  REAL(r_quad) :: cvar, rch

  nargs = iargc()
  IF( nargs.EQ.0 ) THEN
     CALL usage
  ELSE
     debug = .FALSE.
     append_suffix = .FALSE.

     ALLOCATE(arg(nargs))
     DO n=1,nargs
        CALL getarg(n,arg(n))
     ENDDO
     DO n=1,nargs
        IF (TRIM(arg(n)).EQ.'-debug'     ) debug=.TRUE.
        IF (TRIM(arg(n)).EQ.'-append_nc4') append_suffix=.TRUE.
     ENDDO
  ENDIF

  IF (debug) WRITE(*,*)'Debugging on - Verbose Printing'

! get infn from command line
  CALL getarg(nargs, infn)

  strlen = LEN(TRIM(infn))

  WRITE(*,*)'Input bin diag:  ',TRIM(infn)
  INQUIRE(file=TRIM(infn), exist=linfile)
  IF (.NOT. linfile) THEN
     WRITE(*,*)TRIM(infn) // ' does not exist - exiting'
     CALL abort
  ENDIF

  IF (.NOT. append_suffix) THEN
     outfn = infn(1:strlen-3) // 'nc4'  ! assumes GMAO diag filename format ending with .bin, and replaces it
  ELSE
     outfn = infn(1:strlen) // '.nc4'    ! if not GMAO format, use append_suffix = .true. in namelist 
!   to simply append infile with .nc4 suffix
  ENDIF

  WRITE(*,*)'Output NC4 diag: ',TRIM(outfn)
  INQUIRE(file=TRIM(outfn), exist=loutfile)
  IF (loutfile) WRITE(*,*)'WARNING: ' // TRIM(infn) // ' exists - overwriting'

  iflag = 0

  OPEN(inlun,file=infn,form='unformatted',convert='big_endian')
  CALL nc_diag_init(outfn)

  CALL read_aoddiag_header( inlun, header_fix, header_chan, headname, iflag, debug )

  CALL nc_diag_chaninfo_dim_set(header_fix%nchan)

  CALL nc_diag_header("Satellite_Sensor",     header_fix%isis    )
  CALL nc_diag_header("Satellite",            header_fix%id      ) ! sat type
  CALL nc_diag_header("Observation_type",     header_fix%obstype ) ! observation type
  CALL nc_diag_header("Outer_Loop_Iteration", header_fix%jiter   )
  CALL nc_diag_header("Number_of_channels",   header_fix%nchan   ) ! number of channels in the sensor
  CALL nc_diag_header("date_time",            header_fix%idate   ) ! time (yyyymmddhh)
  CALL nc_diag_header("ireal_aoddiag",         header_fix%ireal   )
  CALL nc_diag_header("ipchan_aoddiag",        header_fix%ipchan  )
  CALL nc_diag_header("ioff0",                header_fix%isens   ) ! i think ioff0 = isens


  nch = header_fix%nchan

  ALLOCATE(luse(nch))

  IF (debug) THEN
     WRITE(*,*)'Number of Channels:                 ',nch
  ENDIF

  DO i=1,nch
     CALL nc_diag_chaninfo("chaninfoidx",     i                   )
     CALL nc_diag_chaninfo("frequency",       header_chan(i)%freq    )
     CALL nc_diag_chaninfo("polarization",    header_chan(i)%polar   )
     CALL nc_diag_chaninfo("wavenumber",      header_chan(i)%wave    )
     CALL nc_diag_chaninfo("error_variance",  header_chan(i)%varch   )
     CALL nc_diag_chaninfo("use_flag",        header_chan(i)%iuse    )
     CALL nc_diag_chaninfo("sensor_chan",     header_chan(i)%nuchan  )
     CALL nc_diag_chaninfo("satinfo_chan",    header_chan(i)%iochan  )
  END DO


  DO WHILE (iflag .GE. 0) ! iflag == 0 means the end of the file
     CALL read_aoddiag_data  ( inlun, header_fix, data_fix, data_chan, iflag )

     IF (iflag .LT. 0) CYCLE

     DO ich=1,nch
        lqcpass = luse(ich) .AND. data_chan(ich)%qcmark .EQ. 0 

        CALL nc_diag_metadata("Channel_Index",         i                                   )
        CALL nc_diag_metadata("Observation_Class",     '    aod'                           )
        CALL nc_diag_metadata("Latitude",              data_fix%lat                              ) ! observation latitude (degrees)
        CALL nc_diag_metadata("Longitude",             data_fix%lon                        ) ! observation longitude (degrees)

        CALL nc_diag_metadata("Psfc",             data_fix%psfc                        ) ! observation surface pressure (hPa)

        CALL nc_diag_metadata("Obs_Time",              data_fix%obstime                   ) ! observation time (hours relative to analysis time)

        CALL nc_diag_metadata("Sol_Zenith_Angle",      data_fix%solzen_ang                              ) ! solar zenith angle (degrees)
        CALL nc_diag_metadata("Sol_Azimuth_Angle",     data_fix%solazm_ang                 ) ! solar azimuth angle (degrees)
        CALL nc_diag_metadata("Observation",                           data_chan(ich)%aodobs  )     ! observed aod
        CALL nc_diag_metadata("Obs_Minus_Forecast_unadjusted",         data_chan(ich)%omgaod )     ! observed - simulated Tb with no bias correction (K)
!       errinv = sqrt(varinv(ich_diag(i)))
        CALL nc_diag_metadata("Inverse_Observation_Error",             data_chan(ich)%errinv            )

!       useflag=one
!       if (iuse_aod(ich(ich_diag(i))) < 1) useflag=-one

        CALL nc_diag_metadata("QC_Flag",                               data_chan(ich)%qcmark  )          ! quality control mark or event indicator

     ENDDO

  ENDDO

! finalize NCDIAG
  CALL nc_diag_write
END PROGRAM convert_aod_diag

SUBROUTINE usage
  WRITE(6,100)
100 FORMAT( "Usage:  ",/,/ &
       "  convert_aod_diag.x <options> <filename>",/,/ &
       "where options:",/ &
       "  -debug              :  Set debug verbosity",/ &
       "  -append_txt         :  Append .txt suffix, instead of replace last three",/ &
       "                             characters (default: replaced)",/ &
       "                             Note:  The GMAO diag files end with .bin or .nc4,",/ &
       "                               which is where fixed 3-char truncation originates",/,/,/ &
       "  Example:",/ &
       "     convert_aod_diag.x nc_4emily.diag_hirs4_n19_ges.20161202_00z.bin",/ &
       "  Output file:",/ &
       "     nc_4emily.diag_hirs4_n19_ges.20161202_00z.nc4",/ &
       )
  STOP

END SUBROUTINE usage

