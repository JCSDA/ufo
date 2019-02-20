program convert_rad_diag

  use read_diag,only:read_radiag_header, diag_header_fix_list, diag_header_chan_list, diag_data_name_list
  use read_diag,only:read_radiag_data, diag_data_fix_list, diag_data_extra_list, diag_data_chan_list
  use ncd_kinds,only: r_quad, r_single

  use nc_diag_write_mod, only: nc_diag_init, nc_diag_header, nc_diag_metadata, &
       nc_diag_write, nc_diag_data2d, nc_diag_chaninfo_dim_set, nc_diag_chaninfo

  implicit none 

  real,parameter::     missing = -9.99e9
  integer,parameter:: imissing = -999999

  integer nargs, iargc, n
  character*256, allocatable ::   arg(:)

  type(diag_header_fix_list )          ::  header_fix
  type(diag_header_chan_list),allocatable  ::  header_chan(:)
  type(diag_data_name_list)            ::  headname

  type(diag_data_fix_list)             ::  data_fix
  type(diag_data_chan_list)  ,allocatable  ::  data_chan(:)
  type(diag_data_extra_list) ,allocatable  ::  dataextra(:,:)


  integer i
!!!  real(r_quad)                         ::  ret_var
!!!  real(r_quad)                         ::  ret_stddev

! commandline variables
  logical                              ::  debug 
  integer                              ::  npred_read
  logical                              ::  sst_ret 
  integer                              ::  iversion 
  logical                              ::  append_suffix 

  character*256 infn, outfn
  logical linfile, loutfile

  integer,parameter                    ::  inlun = 51
  integer,parameter                    ::  outlun= 52
  integer,parameter                    ::  nllun = 53

  integer strlen, iflag
  integer iuse, ich, nch, ipr

  logical,dimension(:),allocatable     :: luse
  logical lqcpass

  integer,parameter                    ::  nvar = 4 ! number of positions in array needed for inline variance calc
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
  integer        :: inobstotal, inobsassim
  real(r_single) :: rvomf_nbc, rvomf_bc, rvtotbias, rvfixbias
  real(r_single),dimension(:),allocatable   :: rvbiasterms
  character(len=13),dimension(:),allocatable :: chfrwn 

  integer,parameter              :: max_npred = 9

  character(len=10),dimension(max_npred) :: biasnames  
  data biasnames   / 'bc_const  ', &
                     'bc_satang ', &
                     'bc_tlap   ', &
                     'bc_tlap2  ', &
                     'bc_clw    ', &
                     'bc_coslat ', &
                     'bc_sinlat ', &
                     'bc_emis   ', &
                     'bc_sst    ' /

  real(r_quad) :: cvar, rch


  nargs = iargc()
  if( nargs.eq.0 ) then
    call usage
  else
    debug = .false.
    npred_read = 7
    sst_ret = .false.
    iversion = -9999
    append_suffix = .false.

    allocate(arg(nargs))
    do n=1,nargs
      call getarg(n,arg(n))
    enddo
    do n=1,nargs
      if (trim(arg(n)).eq.'-debug'     ) debug=.true.
      if (trim(arg(n)).eq.'-sst_ret'   ) sst_ret=.true.
      if (trim(arg(n)).eq.'-append_nc4') append_suffix=.true.
      if (trim(arg(n)).eq.'-npred'     ) read ( arg(n+1),* ) npred_read
      if (trim(arg(n)).eq.'-iversion'  ) read ( arg(n+1),* ) iversion
    enddo
  endif


  if (debug) write(*,*)'Debugging on - Verbose Printing'

  ! get infn from command line
  call getarg(nargs, infn)

  strlen = len(trim(infn))

  write(*,*)'Input bin diag:  ',trim(infn)
  inquire(file=trim(infn), exist=linfile)
  if (.not. linfile) then
    write(*,*)trim(infn) // ' does not exist - exiting'
    call abort
  endif
  
  if (.not. append_suffix) then
    outfn = infn(1:strlen-3) // 'nc4'  ! assumes GMAO diag filename format ending with .bin, and replaces it
  else
    outfn = infn(1:strlen) // '.nc4'    ! if not GMAO format, use append_suffix = .true. in namelist 
                                       !   to simply append infile with .nc4 suffix
  endif

  write(*,*)'Output NC4 diag: ',trim(outfn)
  inquire(file=trim(outfn), exist=loutfile)
  if (loutfile) write(*,*)'WARNING: ' // trim(infn) // ' exists - overwriting'

  iflag = 0

  open(inlun,file=infn,form='unformatted',convert='big_endian')
  call nc_diag_init(outfn)

  call read_radiag_header( inlun, npred_read, sst_ret, header_fix, header_chan, headname, iflag, debug )

  call nc_diag_chaninfo_dim_set(header_fix%nchan)

  call nc_diag_header("Satellite_Sensor",     header_fix%isis    )
  call nc_diag_header("Satellite",            header_fix%id      ) ! sat type
  call nc_diag_header("Observation_type",     header_fix%obstype ) ! observation type
  call nc_diag_header("Outer_Loop_Iteration", header_fix%jiter   )
  call nc_diag_header("Number_of_channels",   header_fix%nchan   ) ! number of channels in the sensor
  call nc_diag_header("Number_of_Predictors", header_fix%npred   ) ! number of updating bias correction predictors
  call nc_diag_header("date_time",            header_fix%idate   ) ! time (yyyymmddhh)
  call nc_diag_header("ireal_radiag",         header_fix%ireal   )
  call nc_diag_header("ipchan_radiag",        header_fix%ipchan  )
  call nc_diag_header("iextra",               header_fix%iextra  )
  call nc_diag_header("jextra",               header_fix%jextra  )
  call nc_diag_header("idiag",                header_fix%idiag   )
  call nc_diag_header("angord",               header_fix%angord  )
  call nc_diag_header("iversion_radiag",      header_fix%iversion)
  call nc_diag_header("New_pc4pred",          header_fix%inewpc  ) ! indicator of newpc4pred (1 on, 0 off)
  call nc_diag_header("ioff0",                header_fix%isens   ) ! i think ioff0 = isens


  nch = header_fix%nchan

  allocate(luse(nch))

  if (debug) then
    write(*,*)'Number of Channels:                 ',nch
    write(*,*)'Number of variationalbc predictors: ',npred_read
    write(*,*)' predictors: ',biasnames(1:npred_read)
    write(*,*)' iversion=',header_fix%iversion
  endif

  if (iversion .gt. 0) then
    write(*,*)'BE AWARE THAT iversion IS BEING OVERRIDEN!'
    write(*,*)' iversion diag, override=',header_fix%iversion,iversion
    write(*,*)' (this was made necessary w/ emis bc...hopefully only temporary)'
    header_fix%iversion = iversion
  endif

  do i=1,nch
    call nc_diag_chaninfo("chaninfoidx",     i                   )
    call nc_diag_chaninfo("frequency",       header_chan(i)%freq    )
    call nc_diag_chaninfo("polarization",    header_chan(i)%polar   )
    call nc_diag_chaninfo("wavenumber",      header_chan(i)%wave    )
    call nc_diag_chaninfo("error_variance",  header_chan(i)%varch   )
    call nc_diag_chaninfo("mean_lapse_rate", header_chan(i)%tlapmean)
    call nc_diag_chaninfo("use_flag",        header_chan(i)%iuse    )
    call nc_diag_chaninfo("sensor_chan",     header_chan(i)%nuchan  )
    call nc_diag_chaninfo("satinfo_chan",    header_chan(i)%iochan  )
  end do


  do while (iflag .ge. 0) ! iflag == 0 means the end of the file
    call read_radiag_data  ( inlun, header_fix, .false., data_fix, data_chan, &
                             dataextra, iflag )

    if (iflag .lt. 0) cycle
    
    do ich=1,nch
       lqcpass = luse(ich) .and. data_chan(ich)%qcmark .eq. 0 

       call nc_diag_metadata("Channel_Index",         i                                   )
       call nc_diag_metadata("Observation_Class",     '    rad'                           )
       call nc_diag_metadata("Latitude",              data_fix%lat                              ) ! observation latitude (degrees)
       call nc_diag_metadata("Longitude",             data_fix%lon                        ) ! observation longitude (degrees)

       call nc_diag_metadata("Elevation",             data_fix%zsges                        ) ! model (guess) elevation at observation location

       call nc_diag_metadata("Obs_Time",              data_fix%obstime                   ) ! observation time (hours relative to analysis time)

       call nc_diag_metadata("Scan_Position",         data_fix%senscn_pos                 ) ! sensor scan position
       call nc_diag_metadata("Sat_Zenith_Angle",      data_fix%satzen_ang                       ) ! satellite zenith angle (degrees)
       call nc_diag_metadata("Sat_Azimuth_Angle",     data_fix%satazm_ang                 ) ! satellite azimuth angle (degrees)
       call nc_diag_metadata("Sol_Zenith_Angle",      data_fix%solzen_ang                              ) ! solar zenith angle (degrees)
       call nc_diag_metadata("Sol_Azimuth_Angle",     data_fix%solazm_ang                 ) ! solar azimuth angle (degrees)
       call nc_diag_metadata("Sun_Glint_Angle",       data_fix%sungln_ang                 ) ! sun glint angle (degrees) (sgagl)

       call nc_diag_metadata("Water_Fraction",        data_fix%water_frac           ) ! fractional coverage by water
       call nc_diag_metadata("Land_Fraction",         data_fix%land_frac           ) ! fractional coverage by land
       call nc_diag_metadata("Ice_Fraction",          data_fix%ice_frac             ) ! fractional coverage by ice
       call nc_diag_metadata("Snow_Fraction",         data_fix%snow_frac            ) ! fractional coverage by snow

       if (.not. sst_ret) then
         call nc_diag_metadata("Water_Temperature",     data_fix%water_temp     ) ! surface temperature over water (K)
         call nc_diag_metadata("Land_Temperature",      data_fix%land_temp      ) ! surface temperature over land (K)
         call nc_diag_metadata("Ice_Temperature",       data_fix%ice_temp      ) ! surface temperature over ice (K)
         call nc_diag_metadata("Snow_Temperature",      data_fix%snow_temp      ) ! surface temperature over snow (K)
         call nc_diag_metadata("Soil_Temperature",      data_fix%soil_temp      ) ! soil temperature (K)
         call nc_diag_metadata("Soil_Moisture",         data_fix%soil_mois ) ! soil moisture
         call nc_diag_metadata("Land_Type_Index",       data_fix%land_type          ) ! surface land type
  
         call nc_diag_metadata("tsavg5",                missing                          ) ! SST first guess used for SST retrieval
         call nc_diag_metadata("sstcu",                 missing                          ) ! NCEP SST analysis at t
         call nc_diag_metadata("sstph",                 missing                          ) ! Physical SST retrieval
         call nc_diag_metadata("sstnv",                 missing                          ) ! Navy SST retrieval
         call nc_diag_metadata("dta",                   missing                          ) ! d(ta) corresponding to sstph
         call nc_diag_metadata("dqa",                   missing                          ) ! d(qa) corresponding to sstph
         call nc_diag_metadata("dtp_avh",               missing                          ) ! data type
       else
         call nc_diag_metadata("Water_Temperature",     missing     ) ! surface temperature over water (K)
         call nc_diag_metadata("Land_Temperature",      missing      ) ! surface temperature over land (K)
         call nc_diag_metadata("Ice_Temperature",       missing      ) ! surface temperature over ice (K)
         call nc_diag_metadata("Snow_Temperature",      missing      ) ! surface temperature over snow (K)
         call nc_diag_metadata("Soil_Temperature",      missing      ) ! soil temperature (K)
         call nc_diag_metadata("Soil_Moisture",         missing ) ! soil moisture
         call nc_diag_metadata("Land_Type_Index",       missing          ) ! surface land type
  
         call nc_diag_metadata("tsavg5",                data_fix%water_temp                           ) ! SST first guess used for SST retrieval
         call nc_diag_metadata("sstcu",                 data_fix%land_temp               ) ! NCEP SST analysis at t
         call nc_diag_metadata("sstph",                 data_fix%ice_temp                           ) ! Physical SST retrieval
         call nc_diag_metadata("sstnv",                 data_fix%snow_temp               ) ! Navy SST retrieval
         call nc_diag_metadata("dta",                   data_fix%soil_temp               ) ! d(ta) corresponding to sstph
         call nc_diag_metadata("dqa",                   data_fix%soil_mois               ) ! d(qa) corresponding to sstph
         call nc_diag_metadata("dtp_avh",               data_fix%land_type               ) ! data type
       endif

       call nc_diag_metadata("Vegetation_Fraction",   data_fix%veg_frac      )
       call nc_diag_metadata("Snow_Depth",            data_fix%snow_depth                )
!      qcdiag1 = slot 25 ; qcdiag2 = slot 26 - simply mapping. not attempting to add logic for missing vals
       call nc_diag_metadata("tpwc_amsua",            missing                )
       call nc_diag_metadata("clw_guess_retrieval",   data_fix%qcdiag1                )

       call nc_diag_metadata("Sfc_Wind_Speed",        data_fix%sfc_wndspd               )
       call nc_diag_metadata("Cloud_Frac",            data_fix%qcdiag1                                 )
       call nc_diag_metadata("CTP",                   data_fix%qcdiag2                                )
       call nc_diag_metadata("CLW",                   data_fix%qcdiag1                             )
       call nc_diag_metadata("TPWC",                  data_fix%qcdiag2                              )
       call nc_diag_metadata("clw_obs",               data_fix%qcdiag1                             )
       call nc_diag_metadata("clw_guess",             data_fix%qcdiag2                           )

       call nc_diag_metadata("Foundation_Temperature",   data_fix%tref                   )       ! reference temperature (Tr) in NSST
       call nc_diag_metadata("SST_Warm_layer_dt",        data_fix%dtw                    )       ! dt_warm at zob
       call nc_diag_metadata("SST_Cool_layer_tdrop",     data_fix%dtc                    )       ! dt_cool at zob
       call nc_diag_metadata("SST_dTz_dTfound",          data_fix%tz_tr                  )       ! d(Tz)/d(Tr)

       call nc_diag_metadata("Observation",                           data_chan(ich)%tbobs  )     ! observed brightness temperature (K)
       call nc_diag_metadata("Obs_Minus_Forecast_adjusted",           data_chan(ich)%omgbc  )     ! observed - simulated Tb with bias corrrection (K)
       call nc_diag_metadata("Obs_Minus_Forecast_unadjusted",         data_chan(ich)%omgnbc )     ! observed - simulated Tb with no bias correction (K)
!       errinv = sqrt(varinv(ich_diag(i)))
       call nc_diag_metadata("Inverse_Observation_Error",             data_chan(ich)%errinv            )

!       useflag=one
!       if (iuse_rad(ich(ich_diag(i))) < 1) useflag=-one

       call nc_diag_metadata("QC_Flag",                               data_chan(ich)%qcmark  )          ! quality control mark or event indicator

       call nc_diag_metadata("Emissivity",                            data_chan(ich)%emiss     )           ! surface emissivity
       call nc_diag_metadata("Weighted_Lapse_Rate",                   data_chan(ich)%tlap        )           ! stability index
       call nc_diag_metadata("dTb_dTs",                               data_chan(ich)%tb_tz         )           ! d(Tb)/d(Ts)

       call nc_diag_metadata("BC_Constant",                           data_chan(ich)%bicons         )             ! constant bias correction term
       call nc_diag_metadata("BC_Scan_Angle",                         data_chan(ich)%biang          )             ! scan angle bias correction term
       call nc_diag_metadata("BC_Cloud_Liquid_Water",                 data_chan(ich)%biclw          )             ! CLW bias correction term
       call nc_diag_metadata("BC_Lapse_Rate_Squared",                 data_chan(ich)%bilap2         )             ! square lapse rate bias correction term
       call nc_diag_metadata("BC_Lapse_Rate",                         data_chan(ich)%bilap          )             ! lapse rate bias correction term
       call nc_diag_metadata("BC_Cosine_Latitude_times_Node",         data_chan(ich)%bicos          )             ! node*cos(lat) bias correction term
       call nc_diag_metadata("BC_Sine_Latitude",                      data_chan(ich)%bisin          )             ! sin(lat) bias correction term
       call nc_diag_metadata("BC_Emissivity",                         data_chan(ich)%biemis         )             ! emissivity sensitivity bias correction term
       if (header_fix%angord .eq. 1) then
          call nc_diag_metadata("BC_Fixed_Scan_Position",             data_chan(ich)%bifix(1)       )
       else if (header_fix%angord .ge. 2) then
          call nc_diag_data2d('BC_angord ',                           data_chan(ich)%bifix          )
       endif
! 
    enddo

  enddo

! finalize NCDIAG
  call nc_diag_write
end program convert_rad_diag

subroutine usage
     write(6,100)
100  format( "Usage:  ",/,/ &
             "  convert_rad_diag.x <options> <filename>",/,/ &
             "where options:",/ &
             "  -debug              :  Set debug verbosity",/ &
             "  -sst_ret            :  SST BC term is included (default: not included)",/ &
             "  -npred   INT        :  Number of preductors (default: 7)",/ &
             "  -iversion INT       :  Override iversion with INT (default: use internal iversion)",/   &
             "  -append_txt         :  Append .txt suffix, instead of replace last three",/ &
             "                             characters (default: replaced)",/ &
             "                             Note:  The GMAO diag files end with .bin or .nc4,",/ &
             "                               which is where fixed 3-char truncation originates",/,/,/ &
             "  Example:",/ &
             "     convert_rad_diag.x nc_4emily.diag_hirs4_n19_ges.20161202_00z.bin",/ &
             "  Output file:",/ &
             "     nc_4emily.diag_hirs4_n19_ges.20161202_00z.nc4",/ &
    )
    stop

end subroutine usage

