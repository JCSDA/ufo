program gsidiag_bin2txt

  use nc_diag_read_mod,only: nc_diag_read_init
  use read_diag,only:read_radiag_header, diag_header_fix_list, diag_header_chan_list, diag_data_name_list
  use read_diag,only:read_radiag_data, diag_data_fix_list, diag_data_extra_list, diag_data_chan_list, set_netcdf_read
  use kinds,only: r_quad, r_single

  implicit none 

  integer nargs, iargc, n
  character*256, allocatable ::   arg(:)

  type(diag_header_fix_list )          ::  headfix
  type(diag_header_chan_list),allocatable  ::  headchan(:)
  type(diag_data_name_list)            ::  headname

  type(diag_data_fix_list)             ::  datafix
  type(diag_data_chan_list)  ,allocatable  ::  datachan(:)
  type(diag_data_extra_list) ,allocatable  ::  dataextra(:,:)

  real(r_quad)                         ::  ret_var
  real(r_quad)                         ::  ret_stddev

! optional namelist inputs - can be overriden in a radmon_diag_bin2txt.nl
  logical                              ::  debug = .false.
  integer                              ::  npred_read = 7
  logical                              ::  sst_ret = .false.
  integer                              ::  iversion = -9999
  logical                              ::  append_txt_suffix = .false.

  logical                              ::  netcdf = .false.
  character*256 infn, outfn

!  integer,parameter                    ::  inlun = 51
  integer                              ::  inlun
  integer,parameter                    ::  outlun= 52
  integer,parameter                    ::  nllun = 53

  integer strlen, iflag
  integer iuse, ich, nch, ipr, counter

  logical,dimension(:),allocatable     :: luse
  logical lqcpass

  real(r_single),parameter             ::  missing = -9999.999
  integer,parameter                    ::  imissing = -9999
  integer,parameter                    ::  nvar = 4 ! number of positions in array needed for inline variance calc
! variables for output, all to be allocated as nchan.  Variances will be calculated using inline algorithm
! accredited to Welford, according
  real(r_quad),dimension(:),allocatable   :: nobstotal, nobsassim, tbtotal, tbassim, omf_nbc , omf_bc , sigo, jo
  real(r_quad),dimension(:,:),allocatable ::                                         vomf_nbc, vomf_bc
! total bias and fixed bias terms.  When Yanqui's variational angle correction is brought in, this may need to be updated.
  real(r_quad),dimension(:),allocatable   :: totbias , fixbias
  real(r_quad),dimension(:,:),allocatable :: vtotbias, vfixbias
! variational bias correction variables, which will be allocated as nchan and npred_read
  real(r_quad),dimension(:,:),allocatable   :: biasterms
  real(r_quad),dimension(:,:,:),allocatable :: vbiasterms
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

  logical linfile
  character*80                         ::  nlfn = './gsidiag_bin2txt.nl'


  nargs = iargc()
  if( nargs.eq.0 ) then
    call usage
  else
    netcdf = .false.
    debug = .false.
    npred_read = 7
    sst_ret = .false.
    iversion = -9999
    append_txt_suffix = .false.

    allocate(arg(nargs))
    do n=1,nargs
      call getarg(n,arg(n))
    enddo 
    do n=1,nargs
      if (trim(arg(n)).eq.'-nc4'       ) netcdf=.true.
      if (trim(arg(n)).eq.'-debug'     ) debug=.true.
      if (trim(arg(n)).eq.'-sst_ret'   ) sst_ret=.true.
      if (trim(arg(n)).eq.'-append_txt') append_txt_suffix=.true.
      if (trim(arg(n)).eq.'-npred'     ) read ( arg(n+1),* ) npred_read
      if (trim(arg(n)).eq.'-iversion'  ) read ( arg(n+1),* ) iversion
    enddo
  endif

  if (debug) write(*,*)'Debugging on - Verbose Printing'

  ! get infn from command line
  infn = arg(nargs)

  strlen = len(trim(infn))

  write(*,*)'Input diag file:     ',trim(infn)
  inquire(file=trim(infn), exist=linfile)
  if (.not. linfile) then
    write(*,*)trim(infn) // ' does not exist - exiting'
    call abort
  endif
  
  if (.not. append_txt_suffix) then
    outfn = infn(1:strlen-3) // 'txt'  ! assumes GMAO diag filename format ending with .bin, and replaces it
  else
    outfn = infn(1:strlen) // '.txt'    ! if not GMAO format, use append_txt_suffix = .true. in namelist 
                                       !   to simply append infile with .txt suffix
  endif

  write(*,*)'Output text summary: ',trim(outfn)

  iflag = 0

  if (netcdf) then
    call set_netcdf_read(.true.)
    call nc_diag_read_init(infn, inlun)
  else
    open(inlun,file=infn,form='unformatted',convert='big_endian')
  endif

  write(*,*)'File opened on lun=',inlun
!  open(inlun,file=infn,form='unformatted',convert='big_endian')

  call read_radiag_header( inlun, npred_read, sst_ret, headfix, headchan, headname, iflag, debug )

  nch = headfix%nchan
  allocate(luse(nch))

  if (debug) then
    write(*,*)'Number of Channels:                 ',nch
    write(*,*)'Number of variationalbc predictors: ',npred_read
    write(*,*)' predictors: ',biasnames(1:npred_read)
    write(*,*)' iversion=',headfix%iversion
  endif

  if (iversion .gt. 0) then
    write(*,*)'BE AWARE THAT iversion IS BEING OVERRIDEN!'
    write(*,*)' iversion diag, override=',headfix%iversion,iversion
    write(*,*)' (this was made necessary w/ emis bc...hopefully only temporary)'
    headfix%iversion = iversion
  endif

  allocate(nobstotal(nch),                  &
           nobsassim(nch),                  &
           tbtotal(nch),                    &
           tbassim(nch),                    &
           omf_nbc(nch),                    &
           omf_bc(nch),                     &
           sigo(nch),                       &
           jo(nch),                         &
           totbias(nch),                    &
           fixbias(nch)                     )
  allocate(vomf_nbc(nvar,nch),              &
           vomf_bc(nvar,nch),               &
           vtotbias(nvar,nch),              &
           vfixbias(nvar,nch)               )
  allocate(biasterms(nch,max_npred)        )
  allocate(vbiasterms(nvar,nch,max_npred)  )
  allocate(rvbiasterms(max_npred)  )
  allocate(chfrwn(nch)                     )

  nobstotal = 0.0
  nobsassim = 0.0
  tbtotal = 0.0
  tbassim = 0.0
  omf_nbc = 0.0
  omf_bc = 0.0
  sigo = 0.0
  jo = 0.0
  totbias = 0.0
  fixbias = 0.0
  vomf_nbc = 0.0
  vomf_bc = 0.0
  vtotbias = 0.0
  vfixbias = 0.0
  biasterms = 0.0
  vbiasterms = 0.0
  rvbiasterms = 0.0

  do ich=1,nch
    luse(ich) = (headchan(ich)%iuse .gt. 0) 
    if (headchan(ich)%wave .gt. 100.0) then
      write(chfrwn(ich),fmt='(F9.3,A4)')headchan(ich)%wave,'cm-1'
    else
      write(chfrwn(ich),fmt='(F9.3,A4)')headchan(ich)%freq,'GHz '
    endif
    if (debug) write(*,*)'ich,chfreq or wn=',ich,chfrwn(ich)
  enddo
  counter = 0

  do while (iflag .ge. 0) ! iflag == 0 means the end of the file
    call read_radiag_data  ( inlun, headfix, .false., datafix, datachan, &
                             dataextra, iflag )

    if (iflag .lt. 0) cycle
    counter = counter + 1
!    print *,counter,datafix%lon,datafix%lat

    do ich=1,nch
      lqcpass = luse(ich) .and. datachan(ich)%qcmark .eq. 0 

      ! check to make sure ob is realistic - SSMI seems to have the occasional bad ob sneak in
      if (datachan(ich)%tbobs .gt. 0.0 .and. datachan(ich)%tbobs .lt. 450) then 

        ! first, operations for all observations regardless of luse
        nobstotal(ich) = nobstotal(ich) + 1
        if (debug .and. nobstotal(ich) .lt. 15) print *,nobstotal(ich),ich,datachan(ich)%tbobs

        tbtotal(ich)   = tbtotal(ich) + datachan(ich)%tbobs
        if (luse(ich)) then
          if (lqcpass) then
            nobsassim(ich) = nobsassim(ich) + 1
            tbassim(ich)   = tbassim(ich) + datachan(ich)%tbobs
            omf_nbc(ich)   = omf_nbc(ich) + datachan(ich)%omgnbc
            call inc_var(datachan(ich)%omgnbc, vomf_nbc(:,ich))
            omf_bc(ich)    = omf_bc(ich) + datachan(ich)%omgbc
            call inc_var(datachan(ich)%omgbc, vomf_bc(:,ich))
            sigo(ich)      = sigo(ich) + 1.0 / datachan(ich)%errinv
            jo(ich)        = jo(ich) + (datachan(ich)%omgbc * datachan(ich)%errinv)**2
            totbias(ich)   = totbias(ich) + ( datachan(ich)%omgnbc - datachan(ich)%omgbc )
            call inc_var(datachan(ich)%omgnbc - datachan(ich)%omgbc, vtotbias(:,ich))
            fixbias(ich)   = fixbias(ich) + datachan(ich)%bifix(1)
            call inc_var(datachan(ich)%bifix(1), vfixbias(:,ich))
            biasterms(ich,1) = biasterms(ich,1) + datachan(ich)%bicons
            call inc_var(datachan(ich)%bicons, vbiasterms(:,ich,1))
            biasterms(ich,2) = biasterms(ich,2) + datachan(ich)%biang
            call inc_var(datachan(ich)%biang, vbiasterms(:,ich,2))
            biasterms(ich,3) = biasterms(ich,3) + datachan(ich)%bilap
            call inc_var(datachan(ich)%bilap, vbiasterms(:,ich,3))
            biasterms(ich,4) = biasterms(ich,4) + datachan(ich)%bilap2
            call inc_var(datachan(ich)%bilap2, vbiasterms(:,ich,4))
            biasterms(ich,5) = biasterms(ich,5) + datachan(ich)%biclw
            call inc_var(datachan(ich)%biclw, vbiasterms(:,ich,5))
            biasterms(ich,6) = biasterms(ich,6) + datachan(ich)%bicos
            call inc_var(datachan(ich)%bicos, vbiasterms(:,ich,6))
            biasterms(ich,7) = biasterms(ich,7) + datachan(ich)%bisin
            call inc_var(datachan(ich)%bisin, vbiasterms(:,ich,7))
            biasterms(ich,8) = biasterms(ich,8) + datachan(ich)%biemis
            call inc_var(datachan(ich)%biemis, vbiasterms(:,ich,8))
            biasterms(ich,9) = biasterms(ich,9) + datachan(ich)%bisst
            call inc_var(datachan(ich)%bisst, vbiasterms(:,ich,9))
          endif
        else
          omf_nbc(ich)   = omf_nbc(ich) + datachan(ich)%omgnbc
          call inc_var(datachan(ich)%omgnbc, vomf_nbc(:,ich))
        endif
 
      endif 
    enddo

  enddo

  open(unit=outlun,file=outfn)


  do ich=1,nch 
    inobstotal   = nobstotal(ich)
    if (nobstotal(ich) .gt. 1) then
      tbtotal(ich) = tbtotal(ich) / nobstotal(ich)
  
      if (luse(ich)) then
        inobsassim   = nobsassim(ich)
        if (nobsassim(ich) .gt. 0) then
          tbassim(ich) = tbassim(ich) / nobsassim(ich)
          omf_nbc(ich) = omf_nbc(ich) / nobsassim(ich)
          rvomf_nbc    = ret_stddev(vomf_nbc(:,ich)) 
          omf_bc(ich)  = omf_bc(ich)  / nobsassim(ich)
          rvomf_bc     = ret_stddev(vomf_bc(:,ich))
          sigo(ich)    = sigo(ich) / nobsassim(ich)
          jo(ich)      = jo(ich) / nobsassim(ich)
          totbias(ich) = totbias(ich) / nobsassim(ich)
          rvtotbias    = ret_stddev(vtotbias(:,ich))
          fixbias(ich) = fixbias(ich) / nobsassim(ich)
          rvfixbias    = ret_stddev(vfixbias(:,ich))
          do ipr=1,max_npred
            biasterms(ich,ipr) = biasterms(ich,ipr) / nobsassim(ich)
            rvbiasterms(ipr) = ret_stddev(vbiasterms(:,ich,ipr))
          enddo
        else  ! if zero obs assimilated, pass missings
          tbassim(ich) = missing
          omf_nbc(ich) = missing
          rvomf_nbc    = missing
          omf_bc(ich)  = missing
          rvomf_bc     = missing
          sigo(ich)    = missing
          jo(ich)      = missing
          totbias(ich) = missing
          rvtotbias    = missing
          fixbias(ich) = missing
          rvfixbias    = missing
          do ipr=1,max_npred
            biasterms(ich,ipr) = missing
            rvbiasterms(ipr) = missing
          enddo
        endif
      else
        inobsassim   = imissing
        tbassim(ich) = missing
        omf_nbc(ich) = omf_nbc(ich) / nobstotal(ich)
        rvomf_nbc    = ret_stddev(vomf_nbc(:,ich))
        omf_bc(ich)  = missing
        rvomf_bc     = missing
        sigo(ich)    = missing
        jo(ich)      = missing
        totbias(ich) = missing
        rvtotbias    = missing
        fixbias(ich) = missing
        rvfixbias    = missing
        do ipr=1,max_npred
          biasterms(ich,ipr) = missing
          rvbiasterms(ipr) = missing
        enddo
      endif
    else
      tbtotal(ich) = missing
      inobsassim   = imissing
      tbassim(ich) = missing
      omf_nbc(ich) = missing
      rvomf_nbc    = missing
      omf_bc(ich)  = missing
      rvomf_bc     = missing
      sigo(ich)    = missing
      jo(ich)      = missing
      totbias(ich) = missing
      rvtotbias    = missing
      fixbias(ich) = missing
      rvfixbias    = missing
      do ipr=1,max_npred
        biasterms(ich,ipr) = missing
        rvbiasterms(ipr) = missing
      enddo
    endif
    if (npred_read .lt. max_npred) then
      biasterms(ich,npred_read+1:max_npred) = missing
      rvbiasterms(npred_read+1:max_npred) = missing
    endif
    if (ich .eq. 1) then
      ! write header
      write(unit=outlun,fmt='(A1,A19,3x,A10,3x,A5)'),'!','Satellite/Sensor','YYYYMMDDHH','#chan'
      write(unit=outlun,fmt='(A20,3x,I10,3x,I5)'),trim(headfix%isis), headfix%idate, headfix%nchan

!_RT  write(unit=outlun,fmt='(A6,A1,A13,A1,A4,A1,A12,A1,A12,30(A1,A9))')'!ichan','|','freq/wavenum','|','iuse','|','#total obs','|', &
!_RT             '#assim obs','|','Tb-Total','|','Tb-Assim','|','O-F noBC','|','','|','O-F BC','|','','|','Obs Error','|','Cost (Jo)','|','bc_total','|','','|',   &
!_RT             'bc_fixang','|','',('|',biasnames(ipr),'|','',ipr=1,max_npred)
!_RT  write(unit=outlun,fmt='(A6,A1,A13,A1,A4,A1,A12,A1,A12,30(A1,A9))')'!     ','|',''            ,'|',''    ,'|',''          ,'|', &
!_RT             ''          ,'|','mean','|','mean','|','mean','|','stddev','|','mean','|','stddev','|','mean','|','mean','|','mean','|','stddev','|','mean','|','stddev',('|','mean','|','stddev',ipr=1,max_npred)
    endif
!_RTwrite(unit=outlun,fmt='(I6,1x,A13,1x,I4,1x,I12,1x,I12,1x,30(f9.3,1x))')headchan(ich)%nuchan,chfrwn(ich),headchan(ich)%iuse,inobstotal, &
!_RT       inobsassim,tbtotal(ich),tbassim(ich),omf_nbc(ich),rvomf_nbc,omf_bc(ich),rvomf_bc,sigo(ich),jo(ich),  &
!_RT       totbias(ich),rvtotbias,fixbias(ich),rvfixbias,(biasterms(ich,ipr),rvbiasterms(ipr),ipr=1,max_npred)
  enddo
end program gsidiag_bin2txt

subroutine inc_var(x,arr)
  use kinds,only: r_quad, r_single

  real(r_single)           ,intent(in)    :: x
  real(r_quad),dimension(4),intent(inout) :: arr

  arr(1) = arr(1) + 1
  arr(2) = x - arr(3)
  arr(3) = arr(3) + arr(2)/arr(1)
  arr(4) = arr(4) + arr(2)*(x-arr(3))

end subroutine inc_var

real(r_quad) function ret_var(arr)
  use kinds,only: r_quad

  real(r_quad),dimension(4),intent(in)  :: arr

  ret_var = arr(4) / (arr(1)- 1)
  return
end function ret_var

real(r_quad) function ret_stddev(arr)
  use kinds,only: r_quad

  real(r_quad),dimension(4),intent(in)  :: arr

  ret_stddev = (arr(4) / (arr(1)- 1))**(0.5)
  return
end function ret_stddev

subroutine usage
     write(6,100)
100  format( "Usage:  ",/,/ &
             "  gsidiag_bin2txt.x <options> <filename>",/,/ &
             "where options:",/ &
             "  -nc4                :  Read NC4 Diag (instead of binary)",/ &
             "  -debug              :  Set debug verbosity",/ &
             "  -sst_ret            :  SST BC term is included (default: not included)",/ &
             "  -npred   INT        :  Number of preductors (default: 7)",/ &
             "  -iversion INT       :  Override iversion with INT (default: use internal iversion)",/   &
             "  -append_txt         :  Append .txt suffix, instead of replace last three",/ &
             "                             characters (default: replaced)",/ &
             "                             Note:  The GMAO diag files end with .bin or .nc4,",/ &
             "                               which is where fixed 3-char truncation originates",/,/,/ &
             "  Example:",/ &
             "     gsidiag_bin2txt.x nc_4emily_nc4.diag_airs_aqua_ges.20161202_06z.nc4",/ &
             "  Output file:",/ &
             "     nc_4emily_nc4.diag_airs_aqua_ges.20161202_06z.txt",/ &
    )
    stop            
end subroutine usage
