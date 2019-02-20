module m_diag_conv
  use ncd_kinds, only:  i_kind,r_single,r_kind,r_double
  use nc_diag_write_mod, only: nc_diag_init, nc_diag_metadata, nc_diag_write

  use nc_diag_read_mod, only: nc_diag_read_get_dim
  use nc_diag_read_mod, only: nc_diag_read_init, nc_diag_read_close
  use nc_diag_read_mod, only: nc_diag_read_get_var
  use nc_diag_read_mod, only: nc_diag_read_get_global_attr

  implicit none

  private
  save

  public :: diag_conv_header
  public :: diag_conv_mass  !generic name for non-wind obs - obs w/ single obs associated
  public :: diag_conv_wind  !        name for wind obs - obs w/ two obs (u,v) associated

  public :: open_conv_diag
  public :: read_conv_diag

  public :: write_split_conv_diag_nc

  public :: read_conv_diag_nc_header
  public :: read_conv_diag_nc_mass

  type diag_conv_header
    character(3),    dimension(:), allocatable :: ObsType
    integer(i_kind)                            :: n_ObsType
    integer(i_kind), dimension(:), allocatable :: n_Observations
    integer(i_kind)                            :: n_Observations_Mass
    integer(i_kind)                            :: n_Observations_Wind
    integer(i_kind)                            :: n_Observations_Total
    integer(i_kind)                            :: date
  end type diag_conv_header

  type diag_conv_mass 
    character(8) :: Station_ID
    character(3) :: Observation_Class
    real(r_kind) :: Observation_Type
    real(r_kind) :: Observation_Subtype
    real(r_kind) :: Latitude
    real(r_kind) :: Longitude
    real(r_kind) :: Station_Elevation
    real(r_kind) :: Pressure
    real(r_kind) :: Height
    real(r_kind) :: Time
    real(r_kind) :: Prep_QC_Mark
    real(r_kind) :: Setup_QC_Mark
    real(r_kind) :: Prep_Use_Flag
    real(r_kind) :: Analysis_Use_Flag
    real(r_kind) :: Nonlinear_QC_Rel_Wgt
    real(r_kind) :: Errinv_Input
    real(r_kind) :: Errinv_Adjust
    real(r_kind) :: Errinv_Final
    real(r_kind) :: Observation
    real(r_kind) :: Obs_Minus_Forecast_adjusted
    real(r_kind) :: Obs_Minus_Forecast_unadjusted

  end type diag_conv_mass

  type diag_conv_wind
    character(8) :: Station_ID
    character(3) :: Observation_Class   
    real(r_kind) :: Observation_Type
    real(r_kind) :: Observation_Subtype
    real(r_kind) :: Latitude
    real(r_kind) :: Longitude
    real(r_kind) :: Station_Elevation
    real(r_kind) :: Pressure
    real(r_kind) :: Height
    real(r_kind) :: Time
    real(r_kind) :: Prep_QC_Mark
    real(r_kind) :: Setup_QC_Mark
    real(r_kind) :: Prep_Use_Flag
    real(r_kind) :: Analysis_Use_Flag
    real(r_kind) :: Nonlinear_QC_Rel_Wgt
    real(r_kind) :: Errinv_Input
    real(r_kind) :: Errinv_Adjust
    real(r_kind) :: Errinv_Final
    real(r_kind) :: u_Observation
    real(r_kind) :: u_Obs_Minus_Forecast_adjusted
    real(r_kind) :: u_Obs_Minus_Forecast_unadjusted
    real(r_kind) :: v_Observation
    real(r_kind) :: v_Obs_Minus_Forecast_adjusted
    real(r_kind) :: v_Obs_Minus_Forecast_unadjusted 
    real(r_kind) :: Wind_Reduction_Factor_at_10m
  end type diag_conv_wind

  integer,parameter :: maxobstype=30
  integer           :: nobstype

  integer,parameter :: lun=413
contains

  subroutine read_conv_diag(fn, conv_header, conv_mass, conv_wind, nobs_mass, nobs_wind, ncep)
     character(120), intent(in)              :: fn
     type(diag_conv_header), intent(inout)   :: conv_header
     type(diag_conv_mass), dimension(:), allocatable, intent(inout)     :: conv_mass
     type(diag_conv_wind), dimension(:), allocatable, intent(inout)     :: conv_wind
     integer(i_kind), intent(out)                                       :: nobs_mass, nobs_wind
     logical, intent(in)                                                :: ncep

     integer(i_kind) :: ios
     integer(i_kind) :: date
     character(3)    :: obstype
     integer(i_kind) :: nchar, ninfo, nobs, mype, ioff
     character(8),allocatable, dimension(:)   :: cdiagbuf
     real(r_single),     allocatable, dimension(:,:) :: diagbuf

     integer i, cobmass, cobwind

     nobs_wind = conv_header%n_Observations( get_obstype_index(' uv',conv_header%ObsType))
     nobs_mass = conv_header%n_Observations_Total - nobs_wind

     print *,'Mass, Wind obs count:',nobs_mass, nobs_wind

     open(unit=lun,file=trim(fn), iostat=ios, form='unformatted')

     read(lun) date
     print *,'Date=',date

     allocate( conv_mass(nobs_mass), conv_wind(nobs_wind) )

     cobmass = 1
     cobwind = 1

     do while (ios .eq. 0)
        if (ncep) then 
           read(lun,iostat=ios) obstype,nchar,ninfo,nobs,mype
        else
           read(lun,iostat=ios) obstype,nchar,ninfo,nobs,mype,ioff
        endif

        if (ios .eq. 0) then
           allocate( cdiagbuf(nobs), diagbuf(ninfo, nobs) )
           read(lun,iostat=ios) cdiagbuf, diagbuf

           do i=1,nobs
              if (obstype .eq. ' uv') then
                 conv_wind(cobwind)%Station_ID                      = cdiagbuf(i)
                 conv_wind(cobwind)%Observation_Class               = obstype
                 conv_wind(cobwind)%Observation_Type                = diagbuf( 1,i)
                 conv_wind(cobwind)%Observation_Subtype             = diagbuf( 2,i)
                 conv_wind(cobwind)%Latitude                        = diagbuf( 3,i)
                 conv_wind(cobwind)%Longitude                       = diagbuf( 4,i)
                 conv_wind(cobwind)%Station_Elevation               = diagbuf( 5,i)
                 conv_wind(cobwind)%Pressure                        = diagbuf( 6,i)
                 conv_wind(cobwind)%Height                          = diagbuf( 7,i)
                 conv_wind(cobwind)%Time                            = diagbuf( 8,i)
                 conv_wind(cobwind)%Prep_QC_Mark                    = diagbuf( 9,i)
                 conv_wind(cobwind)%Setup_QC_Mark                   = diagbuf(10,i)
                 conv_wind(cobwind)%Prep_Use_Flag                   = diagbuf(11,i)
                 conv_wind(cobwind)%Analysis_Use_Flag               = diagbuf(12,i)
                 conv_wind(cobwind)%Nonlinear_QC_Rel_Wgt            = diagbuf(13,i)
                 conv_wind(cobwind)%Errinv_Input                    = diagbuf(14,i)
                 conv_wind(cobwind)%Errinv_Adjust                   = diagbuf(15,i)
                 conv_wind(cobwind)%Errinv_Final                    = diagbuf(16,i)
                 conv_wind(cobwind)%u_Observation                   = diagbuf(17,i)
                 conv_wind(cobwind)%u_Obs_Minus_Forecast_adjusted   = diagbuf(18,i)
                 conv_wind(cobwind)%u_Obs_Minus_Forecast_unadjusted = diagbuf(19,i)
                 conv_wind(cobwind)%v_Observation                   = diagbuf(20,i)
                 conv_wind(cobwind)%v_Obs_Minus_Forecast_adjusted   = diagbuf(21,i)
                 conv_wind(cobwind)%v_Obs_Minus_Forecast_unadjusted = diagbuf(22,i)
                 conv_wind(cobwind)%Wind_Reduction_Factor_at_10m    = diagbuf(23,i)
                 cobwind = cobwind + 1
              else
                 conv_mass(cobmass)%Station_ID                      = cdiagbuf(i)
                 conv_mass(cobmass)%Observation_Class               = obstype
                 conv_mass(cobmass)%Observation_Type                = diagbuf( 1,i)
                 conv_mass(cobmass)%Observation_Subtype             = diagbuf( 2,i)
                 conv_mass(cobmass)%Latitude                        = diagbuf( 3,i)
                 conv_mass(cobmass)%Longitude                       = diagbuf( 4,i)
                 conv_mass(cobmass)%Station_Elevation               = diagbuf( 5,i)
                 conv_mass(cobmass)%Pressure                        = diagbuf( 6,i)
                 conv_mass(cobmass)%Height                          = diagbuf( 7,i)
                 conv_mass(cobmass)%Time                            = diagbuf( 8,i)
                 conv_mass(cobmass)%Prep_QC_Mark                    = diagbuf( 9,i)
                 conv_mass(cobmass)%Setup_QC_Mark                   = diagbuf(10,i)
                 conv_mass(cobmass)%Prep_Use_Flag                   = diagbuf(11,i)
                 conv_mass(cobmass)%Analysis_Use_Flag               = diagbuf(12,i)
                 conv_mass(cobmass)%Nonlinear_QC_Rel_Wgt            = diagbuf(13,i)
                 conv_mass(cobmass)%Errinv_Input                    = diagbuf(14,i)
                 conv_mass(cobmass)%Errinv_Adjust                   = diagbuf(15,i)
                 conv_mass(cobmass)%Errinv_Final                    = diagbuf(16,i)
                 conv_mass(cobmass)%Observation                   = diagbuf(17,i)
                 conv_mass(cobmass)%Obs_Minus_Forecast_adjusted   = diagbuf(18,i)
                 conv_mass(cobmass)%Obs_Minus_Forecast_unadjusted = diagbuf(19,i)
                 cobmass = cobmass + 1
              endif
           enddo





           deallocate( cdiagbuf, diagbuf)


        endif
     end do
     close(lun)

     

  end subroutine read_conv_diag

  subroutine open_conv_diag(fn, conv_header,ncep)
     character(120), intent(in)            :: fn
     type(diag_conv_header), intent(inout)   :: conv_header
     logical, intent(in)                     :: ncep

     character(3),   dimension(maxobstype) :: cobstype
     integer(i_kind),dimension(maxobstype) :: cnobs

     integer(i_kind) :: ios
     integer(i_kind) :: date
     character(3)    :: obstype
     integer(i_kind) :: nchar, ninfo, nobs, mype, ioff
     character(8),allocatable, dimension(:)   :: cdiagbuf
     real(4),     allocatable, dimension(:,:) :: diagbuf


     integer(i_kind) :: idx
     ios = 0

     nobstype=0

     open(unit=lun,file=trim(fn), iostat=ios, form='unformatted')

     read(lun) date
     print *,'Date=',date

     conv_header%date = date
     conv_header%n_Observations_Total = 0

     cnobs(:) = 0

     do while (ios .eq. 0)
        if (ncep) then
           read(lun,iostat=ios) obstype,nchar,ninfo,nobs,mype
        else
           read(lun,iostat=ios) obstype,nchar,ninfo,nobs,mype,ioff
        endif

        if (ios .eq. 0) then
           conv_header%n_Observations_Total = conv_header%n_Observations_Total + nobs
           idx = get_obstype_index(obstype,cobstype)
           cnobs(idx) = cnobs(idx) + nobs
        
           allocate( cdiagbuf(nobs), diagbuf(ninfo, nobs) )        
           read(lun,iostat=ios) cdiagbuf, diagbuf
           deallocate( cdiagbuf, diagbuf)

        endif
     end do
     close(lun)    

     conv_header%n_ObsType = nobstype 
     allocate( conv_header%ObsType(nobstype), conv_header%n_Observations(nobstype) )

     print *,'n_ObsType=',conv_header%n_ObsType
     print *,'obstype, count='
     do idx=1,nobstype
        conv_header%ObsType(idx)        = cobstype(idx)
        conv_header%n_Observations(idx) = cnobs(idx)
        print *,conv_header%ObsType(idx),conv_header%n_Observations(idx)
     enddo

     conv_header%n_Observations_Wind = conv_header%n_Observations( get_obstype_index(' uv',conv_header%ObsType))
     conv_header%n_Observations_Mass = conv_header%n_Observations_Total - conv_header%n_Observations_Wind

  end subroutine open_conv_diag

  integer(i_kind) function get_obstype_index(obstype, obstypearr)
!     integer(i_kind) :: get_obstype_index
     character(3),intent(in)                 :: obstype
     character(3),intent(inout),dimension(*) :: obstypearr

     integer :: i, idx
     logical :: matched

     matched = .false.

     if (nobstype .eq. 0) then
        nobstype = 1
        obstypearr(1) = obstype
        idx = nobstype
        print *,'obstype=',obstype,' set to index',idx
     else
        do i=1,nobstype
           if (obstype .eq. obstypearr(i)) then
              idx = i
              matched = .true.
           endif
        enddo
        if (.not. matched) then
           nobstype = nobstype + 1
           obstypearr(nobstype) = obstype
           idx = nobstype
           print *,'obstype=',obstype,' set to index',idx
        endif
     endif


     get_obstype_index = idx
  
  end function get_obstype_index    


  subroutine write_split_conv_diag_nc(infn,conv_header, conv_mass, conv_wind, append_suffix)
     character(120),                                                 intent(in)    :: infn
     type(diag_conv_header),                                         intent(in)    :: conv_header
     type(diag_conv_mass),dimension(conv_header%n_Observations_Mass),intent(in)    :: conv_mass
     type(diag_conv_wind),dimension(conv_header%n_Observations_Wind),intent(in)    :: conv_wind
     logical,                                                        intent(in)    :: append_suffix

     character(120)  :: outfn
     character(20)   :: str, str2
     integer         :: strlen
     integer         :: i, itype

     do itype=1, conv_header%n_ObsType
        str = conv_header%ObsType(itype)
        if (.not. append_suffix) then
           str2 = 'diag_conv_' // trim(adjustl(str))
           outfn = replace_text(trim(infn),'diag_conv',str2)
           strlen = len(trim(outfn))
           outfn = outfn(1:strlen-3) // 'nc4'  
        else
           outfn = trim(infn) // '.' // trim(adjustl(str)) // '.nc4'
        endif

        print *,outfn
 
        call nc_diag_init(outfn)

        if (conv_header%ObsType(itype) .eq. ' uv') then
           do i=1,conv_header%n_Observations_Wind
              call nc_diag_metadata("Station_ID",              conv_wind(i)%Station_ID)
              call nc_diag_metadata("Observation_Class",             conv_wind(i)%Observation_Class              )
              call nc_diag_metadata("Observation_Type",              conv_wind(i)%Observation_Type               )
              call nc_diag_metadata("Observation_Subtype",           conv_wind(i)%Observation_Subtype            )
              call nc_diag_metadata("Latitude",                      conv_wind(i)%Latitude                       )
              call nc_diag_metadata("Longitude",                     conv_wind(i)%Longitude                      )
              call nc_diag_metadata("Station_Elevation",             conv_wind(i)%Station_Elevation              )
              call nc_diag_metadata("Pressure",                      conv_wind(i)%Pressure                       )
              call nc_diag_metadata("Height",                        conv_wind(i)%Height                         )
              call nc_diag_metadata("Time",                          conv_wind(i)%Time                           )
              call nc_diag_metadata("Prep_QC_Mark",                  conv_wind(i)%Prep_QC_Mark                   )
              call nc_diag_metadata("Setup_QC_Mark",                 conv_wind(i)%Setup_QC_Mark                  )
              call nc_diag_metadata("Prep_Use_Flag",                 conv_wind(i)%Prep_Use_Flag                  )
              call nc_diag_metadata("Analysis_Use_Flag",             conv_wind(i)%Analysis_Use_Flag              )
              call nc_diag_metadata("Nonlinear_QC_Rel_Wgt",          conv_wind(i)%Nonlinear_QC_Rel_Wgt           )
              call nc_diag_metadata("Errinv_Input",                  conv_wind(i)%Errinv_Input                   )
              call nc_diag_metadata("Errinv_Adjust",                 conv_wind(i)%Errinv_Adjust                  )
              call nc_diag_metadata("Errinv_Final",                  conv_wind(i)%Errinv_Final                   )
              call nc_diag_metadata("u_Observation",                   conv_wind(i)%u_Observation                    )
              call nc_diag_metadata("u_Obs_Minus_Forecast_adjusted",   conv_wind(i)%u_Obs_Minus_Forecast_adjusted    )
              call nc_diag_metadata("u_Obs_Minus_Forecast_unadjusted", conv_wind(i)%u_Obs_Minus_Forecast_unadjusted  )
              call nc_diag_metadata("v_Observation",                   conv_wind(i)%v_Observation                    )
              call nc_diag_metadata("v_Obs_Minus_Forecast_adjusted",   conv_wind(i)%v_Obs_Minus_Forecast_adjusted    )
              call nc_diag_metadata("v_Obs_Minus_Forecast_unadjusted", conv_wind(i)%v_Obs_Minus_Forecast_unadjusted  )
              call nc_diag_metadata("Wind_Reduction_Factor_at_10m",    conv_wind(i)%Wind_Reduction_Factor_at_10m     )
           enddo
        else
           do i=1,conv_header%n_Observations_Mass
               if (conv_mass(i)%Observation_Class .eq. conv_header%ObsType(itype) ) then
                  call nc_diag_metadata("Station_ID",                    conv_mass(i)%Station_ID                           )
                  call nc_diag_metadata("Observation_Class",             conv_mass(i)%Observation_Class              ) 
                  call nc_diag_metadata("Observation_Type",              conv_mass(i)%Observation_Type               ) 
                  call nc_diag_metadata("Observation_Subtype",           conv_mass(i)%Observation_Subtype            ) 
                  call nc_diag_metadata("Latitude",                      conv_mass(i)%Latitude                       ) 
                  call nc_diag_metadata("Longitude",                     conv_mass(i)%Longitude                      ) 
                  call nc_diag_metadata("Station_Elevation",             conv_mass(i)%Station_Elevation              ) 
                  call nc_diag_metadata("Pressure",                      conv_mass(i)%Pressure                       ) 
                  call nc_diag_metadata("Height",                        conv_mass(i)%Height                         ) 
                  call nc_diag_metadata("Time",                          conv_mass(i)%Time                           ) 
                  call nc_diag_metadata("Prep_QC_Mark",                  conv_mass(i)%Prep_QC_Mark                   ) 
                  call nc_diag_metadata("Setup_QC_Mark",                 conv_mass(i)%Setup_QC_Mark                  ) 
                  call nc_diag_metadata("Prep_Use_Flag",                 conv_mass(i)%Prep_Use_Flag                  ) 
                  call nc_diag_metadata("Analysis_Use_Flag",             conv_mass(i)%Analysis_Use_Flag              ) 
                  call nc_diag_metadata("Nonlinear_QC_Rel_Wgt",          conv_mass(i)%Nonlinear_QC_Rel_Wgt           ) 
                  call nc_diag_metadata("Errinv_Input",                  conv_mass(i)%Errinv_Input                   ) 
                  call nc_diag_metadata("Errinv_Adjust",                 conv_mass(i)%Errinv_Adjust                  ) 
                  call nc_diag_metadata("Errinv_Final",                  conv_mass(i)%Errinv_Final                   ) 
                  call nc_diag_metadata("Observation",                   conv_mass(i)%Observation                    ) 
                  call nc_diag_metadata("Obs_Minus_Forecast_adjusted",   conv_mass(i)%Obs_Minus_Forecast_adjusted    ) 
                  call nc_diag_metadata("Obs_Minus_Forecast_unadjusted", conv_mass(i)%Obs_Minus_Forecast_unadjusted  ) 

               endif
           enddo
        endif

        call nc_diag_write

     enddo
  end subroutine write_split_conv_diag_nc

  subroutine read_conv_diag_nc_header(infn,conv_header,nobs)
     character(len=*),       intent(in)    :: infn
     type(diag_conv_header), intent(inout) :: conv_header
     integer(i_kind),        intent(inout) :: nobs

     integer(i_kind) :: fid

     nobs=0
     call nc_diag_read_init(infn,fid)
     nobs = nc_diag_read_get_dim(fid,'nobs')
     call nc_diag_read_get_global_attr(fid,"date_time", conv_header%date  )
     conv_header%n_Observations_Mass = nobs
     call nc_diag_read_close(infn)

  end subroutine read_conv_diag_nc_header

  subroutine read_conv_diag_nc_mass(infn, conv_header, conv_mass)
     character(len=*),                                               intent(in)    :: infn
     type(diag_conv_header),                                         intent(inout) :: conv_header
     type(diag_conv_mass),dimension(conv_header%n_Observations_Mass),intent(inout) :: conv_mass

     character(20)   :: str, str2
     integer(i_kind) :: i, itype, fid, nobs
     character(len=8),allocatable, dimension(:) :: c_var
     integer(i_kind), allocatable, dimension(:) :: i_var
     real(r_kind),    allocatable, dimension(:) :: r_var

     call nc_diag_read_init(infn,fid)
     
     nobs=conv_header%n_Observations_Mass
     allocate(c_var(nobs))
     allocate(i_var(nobs))
     allocate(r_var(nobs))

!    if (conv_mass(i)%Observation_Class .eq. conv_header%ObsType(itype) ) then
!       call nc_diag_read_get_var(fid,"Station_ID",                   c_var ); conv_mass(:)%Station_ID = c_var
!       call nc_diag_read_get_var(fid,"Observation_Class",            c_var ); conv_mass(:)%Observation_Class = c_var
        call nc_diag_read_get_var(fid,"Observation_Type",             i_var ); conv_mass(:)%Observation_Type    = i_var
        call nc_diag_read_get_var(fid,"Observation_Subtype",          i_var ); conv_mass(:)%Observation_Subtype = i_var
        call nc_diag_read_get_var(fid,"Latitude",                     r_var ); conv_mass(:)%Latitude   = r_var
        call nc_diag_read_get_var(fid,"Longitude",                    r_var ); conv_mass(:)%Longitude  = r_var
        call nc_diag_read_get_var(fid,"Pressure",                     r_var ); conv_mass(:)%Pressure  = r_var
        call nc_diag_read_get_var(fid,"Height",                       r_var ); conv_mass(:)%Height  = r_var
        call nc_diag_read_get_var(fid,"Time",                         r_var ); conv_mass(:)%Time  = r_var
        call nc_diag_read_get_var(fid,"Prep_QC_Mark",                 r_var ); conv_mass(:)%Prep_QC_Mark  = r_var
        call nc_diag_read_get_var(fid,"Setup_QC_Mark",                r_var ); conv_mass(:)%Setup_QC_Mark  = r_var
        call nc_diag_read_get_var(fid,"Prep_Use_Flag",                r_var ); conv_mass(:)%Prep_Use_Flag  = r_var
        call nc_diag_read_get_var(fid,"Analysis_Use_Flag",            r_var ); conv_mass(:)%Analysis_Use_Flag  = r_var
        call nc_diag_read_get_var(fid,"Nonlinear_QC_Rel_Wgt",         r_var ); conv_mass(:)%Nonlinear_QC_Rel_Wgt  = r_var
        call nc_diag_read_get_var(fid,"Errinv_Input",                 r_var ); conv_mass(:)%Errinv_Input         = r_var
        call nc_diag_read_get_var(fid,"Errinv_Adjust",                r_var ); conv_mass(:)%Errinv_Adjust       = r_var
        call nc_diag_read_get_var(fid,"Errinv_Final",                 r_var ); conv_mass(:)%Errinv_Final       = r_var
        call nc_diag_read_get_var(fid,"Observation",                  r_var ); conv_mass(:)%Observation      = r_var
        call nc_diag_read_get_var(fid,"Obs_Minus_Forecast_adjusted",  r_var ); conv_mass(:)%Obs_Minus_Forecast_adjusted   = r_var
        call nc_diag_read_get_var(fid,"Obs_Minus_Forecast_unadjusted",r_var ); conv_mass(:)%Obs_Minus_Forecast_unadjusted   = r_var
!    endif
     deallocate(r_var)
     deallocate(i_var)
     deallocate(c_var)

     call nc_diag_read_close(infn)
  end subroutine read_conv_diag_nc_mass

  function replace_text (s,text,rep)  result(outs)
     character(*)        :: s,text,rep
     character(len(s)+100) :: outs     ! provide outs with extra 100 char len
     integer             :: i, nt, nr
    
     outs = s ; nt = len_trim(text) ; nr = len_trim(rep)
     i = INDEX(outs,text(:nt)) !; IF (i == 0) EXIT
     outs = outs(:i-1) // rep(:nr) // outs(i+nt:)
  end function replace_text


end module m_diag_conv
