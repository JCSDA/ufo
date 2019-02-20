module m_diag_aircft
  use ncd_kinds, only:  i_kind,r_single,r_kind,r_double
  use nc_diag_write_mod, only: nc_diag_init, nc_diag_metadata, nc_diag_write

  use nc_diag_read_mod, only: nc_diag_read_get_dim
  use nc_diag_read_mod, only: nc_diag_read_init, nc_diag_read_close
  use nc_diag_read_mod, only: nc_diag_read_get_var
  use nc_diag_read_mod, only: nc_diag_read_get_global_attr

  implicit none

  private
  save

  public :: diag_aircft_header
  public :: diag_aircft_mass  !generic name for non-wind obs - obs w/ single obs associated
  public :: diag_aircft_wind  !        name for wind obs - obs w/ two obs (u,v) associated

  public :: write_split_aircft_diag_nc

  public :: read_aircft_diag_nc_header
  public :: read_aircft_diag_nc_mass

  type diag_aircft_header
    character(3),    dimension(:), allocatable :: ObsType
    integer(i_kind)                            :: n_ObsType
    integer(i_kind), dimension(:), allocatable :: n_Observations
    integer(i_kind)                            :: n_Observations_Mass
    integer(i_kind)                            :: n_Observations_Wind
    integer(i_kind)                            :: n_Observations_Total
    integer(i_kind)                            :: date
  end type diag_aircft_header

  type diag_aircft_mass 
    character(8) :: Station_ID
    character(3) :: Observation_Class
    real(r_single) :: Observation_Type
!   real(r_single) :: Observation_Subtype
    real(r_single) :: Latitude
    real(r_single) :: Longitude
    real(r_single) :: Station_Elevation
    real(r_single) :: Pressure
    real(r_single) :: Height
    real(r_single) :: Time
    real(r_single) :: Prep_QC_Mark
    real(r_single) :: Setup_QC_Mark
    real(r_single) :: Prep_Use_Flag
    real(r_single) :: Analysis_Use_Flag
    real(r_single) :: Nonlinear_QC_Rel_Wgt
    real(r_single) :: Errinv_Input
    real(r_single) :: Errinv_Adjust
    real(r_single) :: Errinv_Final
    real(r_single) :: Observation
    real(r_single) :: Obs_Minus_Forecast_adjusted
    real(r_single) :: Obs_Minus_Forecast_unadjusted

  end type diag_aircft_mass

  type diag_aircft_wind
    character(8) :: Station_ID
    character(3) :: Observation_Class   
    real(r_single) :: Observation_Type
!   real(r_single) :: Observation_Subtype
    real(r_single) :: Latitude
    real(r_single) :: Longitude
    real(r_single) :: Station_Elevation
    real(r_single) :: Pressure
    real(r_single) :: Height
    real(r_single) :: Time
    real(r_single) :: Prep_QC_Mark
    real(r_single) :: Setup_QC_Mark
    real(r_single) :: Prep_Use_Flag
    real(r_single) :: Analysis_Use_Flag
    real(r_single) :: Nonlinear_QC_Rel_Wgt
    real(r_single) :: Errinv_Input
    real(r_single) :: Errinv_Adjust
    real(r_single) :: Errinv_Final
    real(r_single) :: u_Observation
    real(r_single) :: u_Obs_Minus_Forecast_adjusted
    real(r_single) :: u_Obs_Minus_Forecast_unadjusted
    real(r_single) :: v_Observation
    real(r_single) :: v_Obs_Minus_Forecast_adjusted
    real(r_single) :: v_Obs_Minus_Forecast_unadjusted 
    real(r_single) :: Wind_Reduction_Factor_at_10m
  end type diag_aircft_wind

  integer,parameter :: maxobstype=30
  integer           :: nobstype
  integer           :: aircft_mass_type = 130 ! type for RAOB T and Q
  integer           :: aircft_wind_type = 230 ! type for RAOB U and V
  integer           :: t_qcmark       = 1   ! 0=tv; 1=tdry

contains

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


  subroutine write_split_aircft_diag_nc(infn,aircft_header, aircft_mass, aircft_wind, append_suffix)
     character(120),                                                 intent(in)    :: infn
     type(diag_aircft_header),                                         intent(in)  :: aircft_header
     type(diag_aircft_mass),dimension(aircft_header%n_Observations_Mass),intent(in):: aircft_mass
     type(diag_aircft_wind),dimension(aircft_header%n_Observations_Wind),intent(in):: aircft_wind
     logical,                                                        intent(in)    :: append_suffix

     character(120)  :: outfn
     character(20)   :: str, str2
     integer         :: strlen
     integer         :: i, itype

     do itype=1, aircft_header%n_ObsType
        str = aircft_header%ObsType(itype)
        if (.not. append_suffix) then
           str2 = 'diag_aircft_' // trim(adjustl(str))
           outfn = replace_text(trim(infn),'diag_aircft',str2)
           strlen = len(trim(outfn))
           outfn = outfn(1:strlen-3) // 'nc4'  
        else
           outfn = trim(infn) // '.' // trim(adjustl(str)) // '.nc4'
        endif

        print *,outfn
 
        call nc_diag_init(outfn)

        if (aircft_header%ObsType(itype) .eq. ' uv') then
           do i=1,aircft_header%n_Observations_Wind
              call nc_diag_metadata("Station_ID",                    aircft_wind(i)%Station_ID                     )
              call nc_diag_metadata("Observation_Class",             aircft_wind(i)%Observation_Class              )
              call nc_diag_metadata("Observation_Type",              aircft_wind(i)%Observation_Type               )
!             call nc_diag_metadata("Observation_Subtype",           aircft_wind(i)%Observation_Subtype            )
              call nc_diag_metadata("Latitude",                      aircft_wind(i)%Latitude                       )
              call nc_diag_metadata("Longitude",                     aircft_wind(i)%Longitude                      )
              call nc_diag_metadata("Station_Elevation",             aircft_wind(i)%Station_Elevation              )
              call nc_diag_metadata("Pressure",                      aircft_wind(i)%Pressure                       )
              call nc_diag_metadata("Height",                        aircft_wind(i)%Height                         )
              call nc_diag_metadata("Time",                          aircft_wind(i)%Time                           )
              call nc_diag_metadata("Prep_QC_Mark",                  aircft_wind(i)%Prep_QC_Mark                   )
              call nc_diag_metadata("Setup_QC_Mark",                 aircft_wind(i)%Setup_QC_Mark                  )
              call nc_diag_metadata("Prep_Use_Flag",                 aircft_wind(i)%Prep_Use_Flag                  )
              call nc_diag_metadata("Analysis_Use_Flag",             aircft_wind(i)%Analysis_Use_Flag              )
              call nc_diag_metadata("Nonlinear_QC_Rel_Wgt",          aircft_wind(i)%Nonlinear_QC_Rel_Wgt           )
              call nc_diag_metadata("Errinv_Input",                  aircft_wind(i)%Errinv_Input                   )
              call nc_diag_metadata("Errinv_Adjust",                 aircft_wind(i)%Errinv_Adjust                  )
              call nc_diag_metadata("Errinv_Final",                  aircft_wind(i)%Errinv_Final                   )
              call nc_diag_metadata("u_Observation",                   aircft_wind(i)%u_Observation                    )
              call nc_diag_metadata("u_Obs_Minus_Forecast_adjusted",   aircft_wind(i)%u_Obs_Minus_Forecast_adjusted    )
              call nc_diag_metadata("u_Obs_Minus_Forecast_unadjusted", aircft_wind(i)%u_Obs_Minus_Forecast_unadjusted  )
              call nc_diag_metadata("v_Observation",                   aircft_wind(i)%v_Observation                    )
              call nc_diag_metadata("v_Obs_Minus_Forecast_adjusted",   aircft_wind(i)%v_Obs_Minus_Forecast_adjusted    )
              call nc_diag_metadata("v_Obs_Minus_Forecast_unadjusted", aircft_wind(i)%v_Obs_Minus_Forecast_unadjusted  )
              call nc_diag_metadata("Wind_Reduction_Factor_at_10m",    aircft_wind(i)%Wind_Reduction_Factor_at_10m     )
           enddo
        else
           do i=1,aircft_header%n_Observations_Mass
               if (aircft_mass(i)%Observation_Class .eq. aircft_header%ObsType(itype) ) then
                  call nc_diag_metadata("Station_ID",                    aircft_mass(i)%Station_ID                           )
                  call nc_diag_metadata("Observation_Class",             aircft_mass(i)%Observation_Class              ) 
                  call nc_diag_metadata("Observation_Type",              aircft_mass(i)%Observation_Type               ) 
!                 call nc_diag_metadata("Observation_Subtype",           aircft_mass(i)%Observation_Subtype            ) 
                  call nc_diag_metadata("Latitude",                      aircft_mass(i)%Latitude                       ) 
                  call nc_diag_metadata("Longitude",                     aircft_mass(i)%Longitude                      ) 
                  call nc_diag_metadata("Station_Elevation",             aircft_mass(i)%Station_Elevation              ) 
                  call nc_diag_metadata("Pressure",                      aircft_mass(i)%Pressure                       ) 
                  call nc_diag_metadata("Height",                        aircft_mass(i)%Height                         ) 
                  call nc_diag_metadata("Time",                          aircft_mass(i)%Time                           ) 
                  call nc_diag_metadata("Prep_QC_Mark",                  aircft_mass(i)%Prep_QC_Mark                   ) 
                  call nc_diag_metadata("Setup_QC_Mark",                 aircft_mass(i)%Setup_QC_Mark                  ) 
                  call nc_diag_metadata("Prep_Use_Flag",                 aircft_mass(i)%Prep_Use_Flag                  ) 
                  call nc_diag_metadata("Analysis_Use_Flag",             aircft_mass(i)%Analysis_Use_Flag              ) 
                  call nc_diag_metadata("Nonlinear_QC_Rel_Wgt",          aircft_mass(i)%Nonlinear_QC_Rel_Wgt           ) 
                  call nc_diag_metadata("Errinv_Input",                  aircft_mass(i)%Errinv_Input                   ) 
                  call nc_diag_metadata("Errinv_Adjust",                 aircft_mass(i)%Errinv_Adjust                  ) 
                  call nc_diag_metadata("Errinv_Final",                  aircft_mass(i)%Errinv_Final                   ) 
                  call nc_diag_metadata("Observation",                   aircft_mass(i)%Observation                    ) 
                  call nc_diag_metadata("Obs_Minus_Forecast_adjusted",   aircft_mass(i)%Obs_Minus_Forecast_adjusted    ) 
                  call nc_diag_metadata("Obs_Minus_Forecast_unadjusted", aircft_mass(i)%Obs_Minus_Forecast_unadjusted  ) 

               endif
           enddo
        endif

        call nc_diag_write

     enddo
  end subroutine write_split_aircft_diag_nc

  subroutine read_aircft_diag_nc_header(infn,aircft_header)
     character(len=*),       intent(in)    :: infn
     type(diag_aircft_header), intent(inout) :: aircft_header

     integer(i_kind) :: fid,nobs

     nobs=0
     call nc_diag_read_init(infn,fid)
     nobs = nc_diag_read_get_dim(fid,'nobs')
     call nc_diag_read_get_global_attr(fid,"date_time", aircft_header%date  )
     aircft_header%n_Observations_Mass = nobs
     aircft_header%n_Observations_Wind = nobs
     call nc_diag_read_close(infn)

  end subroutine read_aircft_diag_nc_header

  subroutine read_aircft_diag_nc_mass(infn, aircft_header, aircft_mass, ierr)
     character(len=*),             intent(in)    :: infn
     type(diag_aircft_header),       intent(inout) :: aircft_header
     type(diag_aircft_mass),pointer, intent(inout) :: aircft_mass(:)
     integer,                      intent(out)   :: ierr

     character(20)   :: str, str2
     integer(i_kind) :: ii, ic, fid, nobs, naircft, ncount(1)
     integer(i_kind), allocatable :: indx(:)
     character(len=8),allocatable, dimension(:) :: c_var
     integer(i_kind), allocatable, dimension(:) :: i_var
     real(r_single),  allocatable, dimension(:) :: r_var

     type(diag_aircft_mass), pointer :: rtmp_mass(:)

     ierr=0
     call nc_diag_read_init(infn,fid)
     
     nobs=aircft_header%n_Observations_Mass
     allocate(rtmp_mass(nobs))
     allocate(c_var(nobs))
     allocate(i_var(nobs))
     allocate(r_var(nobs))

!    if (aircft_mass(i)%Observation_Class .eq. aircft_header%ObsType(itype) ) then
!       call nc_diag_read_get_var(fid,"Station_ID",                   c_var ); rtmp_mass(:)%Station_ID = c_var
!       call nc_diag_read_get_var(fid,"Observation_Class",            c_var ); rtmp_mass(:)%Observation_Class = c_var
        call nc_diag_read_get_var(fid,"Observation_Type",             i_var ); rtmp_mass(:)%Observation_Type    = i_var
!       call nc_diag_read_get_var(fid,"Observation_Subtype",          i_var ); rtmp_mass(:)%Observation_Subtype = i_var
        call nc_diag_read_get_var(fid,"Latitude",                     r_var ); rtmp_mass(:)%Latitude   = r_var
        call nc_diag_read_get_var(fid,"Longitude",                    r_var ); rtmp_mass(:)%Longitude  = r_var
        call nc_diag_read_get_var(fid,"Pressure",                     r_var ); rtmp_mass(:)%Pressure  = r_var
        call nc_diag_read_get_var(fid,"Height",                       r_var ); rtmp_mass(:)%Height  = r_var
        call nc_diag_read_get_var(fid,"Time",                         r_var ); rtmp_mass(:)%Time  = r_var
        call nc_diag_read_get_var(fid,"Prep_QC_Mark",                 r_var ); rtmp_mass(:)%Prep_QC_Mark  = r_var
        call nc_diag_read_get_var(fid,"Setup_QC_Mark",                r_var ); rtmp_mass(:)%Setup_QC_Mark  = r_var
        call nc_diag_read_get_var(fid,"Prep_Use_Flag",                r_var ); rtmp_mass(:)%Prep_Use_Flag  = r_var
        call nc_diag_read_get_var(fid,"Analysis_Use_Flag",            r_var ); rtmp_mass(:)%Analysis_Use_Flag  = r_var
        call nc_diag_read_get_var(fid,"Nonlinear_QC_Rel_Wgt",         r_var ); rtmp_mass(:)%Nonlinear_QC_Rel_Wgt  = r_var
        call nc_diag_read_get_var(fid,"Errinv_Input",                 r_var ); rtmp_mass(:)%Errinv_Input         = r_var
        call nc_diag_read_get_var(fid,"Errinv_Adjust",                r_var ); rtmp_mass(:)%Errinv_Adjust       = r_var
        call nc_diag_read_get_var(fid,"Errinv_Final",                 r_var ); rtmp_mass(:)%Errinv_Final       = r_var
        call nc_diag_read_get_var(fid,"Observation",                  r_var ); rtmp_mass(:)%Observation      = r_var
        call nc_diag_read_get_var(fid,"Obs_Minus_Forecast_adjusted",  r_var ); rtmp_mass(:)%Obs_Minus_Forecast_adjusted   = r_var
        call nc_diag_read_get_var(fid,"Obs_Minus_Forecast_unadjusted",r_var ); rtmp_mass(:)%Obs_Minus_Forecast_unadjusted   = r_var
!    endif
     deallocate(r_var)
     deallocate(i_var)
     deallocate(c_var)
      
     ic=0
     do ii=1,nobs
        if(rtmp_mass(ii)%Observation_Type==aircft_mass_type.and.&
           rtmp_mass(ii)%Setup_QC_Mark==t_qcmark) then
           ic=ic+1
        endif
     enddo
     naircft=ic 
     allocate(indx(naircft))

     ic=0
     do ii=1,nobs
        if(rtmp_mass(ii)%Observation_Type==aircft_mass_type.and.&
             rtmp_mass(ii)%Setup_QC_Mark==t_qcmark) then
             ic=ic+1
             indx(ic)=ii 
          endif
     enddo

     print *, ' found this many aircft ', naircft
     if(ic /= naircft) then
       print *, 'error determining Aircraft, inconsistent naircft, ic=', naircft, ic
       deallocate(indx)
       deallocate(aircft_mass)
       ierr = 99
       return
     endif

     if(associated(aircft_mass)) deallocate(aircft_mass)
     allocate  (aircft_mass(naircft))
     aircft_header%n_Observations_Mass = naircft

!    aircft_mass(:)%Station_ID = rtmp_mass(indx)
!    aircft_mass(:)%Observation_Class = c_var
     aircft_mass(:)%Observation_Type    = rtmp_mass(indx)%Observation_Type
!    aircft_mass(:)%Observation_Subtype = rtmp_mass(indx)%Observation_Subtype
     aircft_mass(:)%Latitude            = rtmp_mass(indx)%Latitude 
     aircft_mass(:)%Longitude           = rtmp_mass(indx)%Longitude
     aircft_mass(:)%Pressure            = rtmp_mass(indx)%Pressure
     aircft_mass(:)%Height              = rtmp_mass(indx)%Height
     aircft_mass(:)%Time                = rtmp_mass(indx)%Time
     aircft_mass(:)%Prep_QC_Mark        = rtmp_mass(indx)%Prep_QC_Mark
     aircft_mass(:)%Setup_QC_Mark       = rtmp_mass(indx)%Setup_QC_Mark
     aircft_mass(:)%Prep_Use_Flag       = rtmp_mass(indx)%Prep_Use_Flag 
     aircft_mass(:)%Analysis_Use_Flag   = rtmp_mass(indx)%Analysis_Use_Flag 
     aircft_mass(:)%Nonlinear_QC_Rel_Wgt= rtmp_mass(indx)%Nonlinear_QC_Rel_Wgt
     aircft_mass(:)%Errinv_Input        = rtmp_mass(indx)%Errinv_Input
     aircft_mass(:)%Errinv_Adjust       = rtmp_mass(indx)%Errinv_Adjust 
     aircft_mass(:)%Errinv_Final        = rtmp_mass(indx)%Errinv_Final
     aircft_mass(:)%Observation         = rtmp_mass(indx)%Observation 
     aircft_mass(:)%Obs_Minus_Forecast_adjusted   = rtmp_mass(indx)%Obs_Minus_Forecast_adjusted
     aircft_mass(:)%Obs_Minus_Forecast_unadjusted = rtmp_mass(indx)%Obs_Minus_Forecast_unadjusted

     deallocate(indx)
     deallocate(rtmp_mass)

     call nc_diag_read_close(infn)
  end subroutine read_aircft_diag_nc_mass

  function replace_text (s,text,rep)  result(outs)
     character(*)        :: s,text,rep
     character(len(s)+100) :: outs     ! provide outs with extra 100 char len
     integer             :: i, nt, nr
    
     outs = s ; nt = len_trim(text) ; nr = len_trim(rep)
     i = INDEX(outs,text(:nt)) !; IF (i == 0) EXIT
     outs = outs(:i-1) // rep(:nr) // outs(i+nt:)
  end function replace_text


end module m_diag_aircft
