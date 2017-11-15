module m_diag_raob
  use ncd_kinds, only:  i_kind,r_single,r_kind,r_double
  use nc_diag_write_mod, only: nc_diag_init, nc_diag_metadata, nc_diag_write

  use nc_diag_read_mod, only: nc_diag_read_get_dim
  use nc_diag_read_mod, only: nc_diag_read_init, nc_diag_read_close
  use nc_diag_read_mod, only: nc_diag_read_get_var
  use nc_diag_read_mod, only: nc_diag_read_get_global_attr

  implicit none

  private
  save

  public :: diag_raob_header
  public :: diag_raob_mass  !generic name for non-wind obs - obs w/ single obs associated
  public :: diag_raob_wind  !        name for wind obs - obs w/ two obs (u,v) associated

  public :: write_split_raob_diag_nc

  public :: read_raob_diag_nc_header
  public :: read_raob_diag_nc_mass

  type diag_raob_header
    character(3),    dimension(:), allocatable :: ObsType
    integer(i_kind)                            :: n_ObsType
    integer(i_kind), dimension(:), allocatable :: n_Observations
    integer(i_kind)                            :: n_Observations_Mass
    integer(i_kind)                            :: n_Observations_Wind
    integer(i_kind)                            :: n_Observations_Total
    integer(i_kind)                            :: date
  end type diag_raob_header

  type diag_raob_mass 
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

  end type diag_raob_mass

  type diag_raob_wind
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
  end type diag_raob_wind

  integer,parameter :: maxobstype=30
  integer           :: nobstype
  integer           :: raob_mass_type = 120 ! type for RAOB T and Q
  integer           :: raob_wind_type = 220 ! type for RAOB U and V

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


  subroutine write_split_raob_diag_nc(infn,raob_header, raob_mass, raob_wind, append_suffix)
     character(120),                                                 intent(in)    :: infn
     type(diag_raob_header),                                         intent(in)    :: raob_header
     type(diag_raob_mass),dimension(raob_header%n_Observations_Mass),intent(in)    :: raob_mass
     type(diag_raob_wind),dimension(raob_header%n_Observations_Wind),intent(in)    :: raob_wind
     logical,                                                        intent(in)    :: append_suffix

     character(120)  :: outfn
     character(20)   :: str, str2
     integer         :: strlen
     integer         :: i, itype

     do itype=1, raob_header%n_ObsType
        str = raob_header%ObsType(itype)
        if (.not. append_suffix) then
           str2 = 'diag_raob_' // trim(adjustl(str))
           outfn = replace_text(trim(infn),'diag_raob',str2)
           strlen = len(trim(outfn))
           outfn = outfn(1:strlen-3) // 'nc4'  
        else
           outfn = trim(infn) // '.' // trim(adjustl(str)) // '.nc4'
        endif

        print *,outfn
 
        call nc_diag_init(outfn)

        if (raob_header%ObsType(itype) .eq. ' uv') then
           do i=1,raob_header%n_Observations_Wind
              call nc_diag_metadata("Station_ID",                    raob_wind(i)%Station_ID                     )
              call nc_diag_metadata("Observation_Class",             raob_wind(i)%Observation_Class              )
              call nc_diag_metadata("Observation_Type",              raob_wind(i)%Observation_Type               )
              call nc_diag_metadata("Observation_Subtype",           raob_wind(i)%Observation_Subtype            )
              call nc_diag_metadata("Latitude",                      raob_wind(i)%Latitude                       )
              call nc_diag_metadata("Longitude",                     raob_wind(i)%Longitude                      )
              call nc_diag_metadata("Station_Elevation",             raob_wind(i)%Station_Elevation              )
              call nc_diag_metadata("Pressure",                      raob_wind(i)%Pressure                       )
              call nc_diag_metadata("Height",                        raob_wind(i)%Height                         )
              call nc_diag_metadata("Time",                          raob_wind(i)%Time                           )
              call nc_diag_metadata("Prep_QC_Mark",                  raob_wind(i)%Prep_QC_Mark                   )
              call nc_diag_metadata("Setup_QC_Mark",                 raob_wind(i)%Setup_QC_Mark                  )
              call nc_diag_metadata("Prep_Use_Flag",                 raob_wind(i)%Prep_Use_Flag                  )
              call nc_diag_metadata("Analysis_Use_Flag",             raob_wind(i)%Analysis_Use_Flag              )
              call nc_diag_metadata("Nonlinear_QC_Rel_Wgt",          raob_wind(i)%Nonlinear_QC_Rel_Wgt           )
              call nc_diag_metadata("Errinv_Input",                  raob_wind(i)%Errinv_Input                   )
              call nc_diag_metadata("Errinv_Adjust",                 raob_wind(i)%Errinv_Adjust                  )
              call nc_diag_metadata("Errinv_Final",                  raob_wind(i)%Errinv_Final                   )
              call nc_diag_metadata("u_Observation",                   raob_wind(i)%u_Observation                    )
              call nc_diag_metadata("u_Obs_Minus_Forecast_adjusted",   raob_wind(i)%u_Obs_Minus_Forecast_adjusted    )
              call nc_diag_metadata("u_Obs_Minus_Forecast_unadjusted", raob_wind(i)%u_Obs_Minus_Forecast_unadjusted  )
              call nc_diag_metadata("v_Observation",                   raob_wind(i)%v_Observation                    )
              call nc_diag_metadata("v_Obs_Minus_Forecast_adjusted",   raob_wind(i)%v_Obs_Minus_Forecast_adjusted    )
              call nc_diag_metadata("v_Obs_Minus_Forecast_unadjusted", raob_wind(i)%v_Obs_Minus_Forecast_unadjusted  )
              call nc_diag_metadata("Wind_Reduction_Factor_at_10m",    raob_wind(i)%Wind_Reduction_Factor_at_10m     )
           enddo
        else
           do i=1,raob_header%n_Observations_Mass
               if (raob_mass(i)%Observation_Class .eq. raob_header%ObsType(itype) ) then
                  call nc_diag_metadata("Station_ID",                    raob_mass(i)%Station_ID                           )
                  call nc_diag_metadata("Observation_Class",             raob_mass(i)%Observation_Class              ) 
                  call nc_diag_metadata("Observation_Type",              raob_mass(i)%Observation_Type               ) 
                  call nc_diag_metadata("Observation_Subtype",           raob_mass(i)%Observation_Subtype            ) 
                  call nc_diag_metadata("Latitude",                      raob_mass(i)%Latitude                       ) 
                  call nc_diag_metadata("Longitude",                     raob_mass(i)%Longitude                      ) 
                  call nc_diag_metadata("Station_Elevation",             raob_mass(i)%Station_Elevation              ) 
                  call nc_diag_metadata("Pressure",                      raob_mass(i)%Pressure                       ) 
                  call nc_diag_metadata("Height",                        raob_mass(i)%Height                         ) 
                  call nc_diag_metadata("Time",                          raob_mass(i)%Time                           ) 
                  call nc_diag_metadata("Prep_QC_Mark",                  raob_mass(i)%Prep_QC_Mark                   ) 
                  call nc_diag_metadata("Setup_QC_Mark",                 raob_mass(i)%Setup_QC_Mark                  ) 
                  call nc_diag_metadata("Prep_Use_Flag",                 raob_mass(i)%Prep_Use_Flag                  ) 
                  call nc_diag_metadata("Analysis_Use_Flag",             raob_mass(i)%Analysis_Use_Flag              ) 
                  call nc_diag_metadata("Nonlinear_QC_Rel_Wgt",          raob_mass(i)%Nonlinear_QC_Rel_Wgt           ) 
                  call nc_diag_metadata("Errinv_Input",                  raob_mass(i)%Errinv_Input                   ) 
                  call nc_diag_metadata("Errinv_Adjust",                 raob_mass(i)%Errinv_Adjust                  ) 
                  call nc_diag_metadata("Errinv_Final",                  raob_mass(i)%Errinv_Final                   ) 
                  call nc_diag_metadata("Observation",                   raob_mass(i)%Observation                    ) 
                  call nc_diag_metadata("Obs_Minus_Forecast_adjusted",   raob_mass(i)%Obs_Minus_Forecast_adjusted    ) 
                  call nc_diag_metadata("Obs_Minus_Forecast_unadjusted", raob_mass(i)%Obs_Minus_Forecast_unadjusted  ) 

               endif
           enddo
        endif

        call nc_diag_write

     enddo
  end subroutine write_split_raob_diag_nc

  subroutine read_raob_diag_nc_header(infn,raob_header)
     character(len=*),       intent(in)    :: infn
     type(diag_raob_header), intent(inout) :: raob_header

     integer(i_kind) :: fid,nobs

     nobs=0
     call nc_diag_read_init(infn,fid)
     nobs = nc_diag_read_get_dim(fid,'nobs')
     call nc_diag_read_get_global_attr(fid,"date_time", raob_header%date  )
     raob_header%n_Observations_Mass = nobs
     call nc_diag_read_close(infn)

  end subroutine read_raob_diag_nc_header

  subroutine read_raob_diag_nc_mass(infn, raob_header, raob_mass, ierr)
     character(len=*),             intent(in)    :: infn
     type(diag_raob_header),       intent(inout) :: raob_header
     type(diag_raob_mass),pointer, intent(inout) :: raob_mass(:)
     integer,                      intent(out)   :: ierr

     character(20)   :: str, str2
     integer(i_kind) :: ii, ic, fid, nobs, nraob, ncount(1)
     integer(i_kind), allocatable :: indx(:)
     character(len=8),allocatable, dimension(:) :: c_var
     integer(i_kind), allocatable, dimension(:) :: i_var
     real(r_kind),    allocatable, dimension(:) :: r_var

     type(diag_raob_mass), pointer :: rtmp_mass(:)

     ierr=0
     call nc_diag_read_init(infn,fid)
     
     nobs=raob_header%n_Observations_Mass
     allocate(rtmp_mass(nobs))
     allocate(c_var(nobs))
     allocate(i_var(nobs))
     allocate(r_var(nobs))

!    if (raob_mass(i)%Observation_Class .eq. raob_header%ObsType(itype) ) then
!       call nc_diag_read_get_var(fid,"Station_ID",                   c_var ); rtmp_mass(:)%Station_ID = c_var
!       call nc_diag_read_get_var(fid,"Observation_Class",            c_var ); rtmp_mass(:)%Observation_Class = c_var
        call nc_diag_read_get_var(fid,"Observation_Type",             i_var ); rtmp_mass(:)%Observation_Type    = i_var
        call nc_diag_read_get_var(fid,"Observation_Subtype",          i_var ); rtmp_mass(:)%Observation_Subtype = i_var
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
      
     ncount = count(rtmp_mass(:)%Observation_Type==raob_mass_type,1); nraob = ncount(1) 
     allocate(indx(nraob))
     ic=0
     do ii=1,nobs
        if(rtmp_mass(ii)%Observation_Type==raob_mass_type) then
          ic=ic+1
          indx(ic)=ii 
        endif
     enddo
     print *, 'found this many raob ', nraob
     if(ic /= nraob) then
       print *, 'error determining Raob, inconsistent nraob, ic=', nraob, ic
       deallocate(indx)
       deallocate(raob_mass)
       ierr = 99
       return
     endif

     if(associated(raob_mass)) deallocate(raob_mass)
     allocate  (raob_mass(nraob))
     raob_header%n_Observations_Mass = nraob

!    raob_mass(:)%Station_ID = rtmp_mass(indx)
!    raob_mass(:)%Observation_Class = c_var
     raob_mass(:)%Observation_Type    = rtmp_mass(indx)%Observation_Type
     raob_mass(:)%Observation_Subtype = rtmp_mass(indx)%Observation_Subtype
     raob_mass(:)%Latitude            = rtmp_mass(indx)%Latitude 
     raob_mass(:)%Longitude           = rtmp_mass(indx)%Longitude
     raob_mass(:)%Pressure            = rtmp_mass(indx)%Pressure
     raob_mass(:)%Height              = rtmp_mass(indx)%Height
     raob_mass(:)%Time                = rtmp_mass(indx)%Time
     raob_mass(:)%Prep_QC_Mark        = rtmp_mass(indx)%Prep_QC_Mark
     raob_mass(:)%Setup_QC_Mark       = rtmp_mass(indx)%Setup_QC_Mark
     raob_mass(:)%Prep_Use_Flag       = rtmp_mass(indx)%Prep_Use_Flag 
     raob_mass(:)%Analysis_Use_Flag   = rtmp_mass(indx)%Analysis_Use_Flag 
     raob_mass(:)%Nonlinear_QC_Rel_Wgt= rtmp_mass(indx)%Nonlinear_QC_Rel_Wgt
     raob_mass(:)%Errinv_Input        = rtmp_mass(indx)%Errinv_Input
     raob_mass(:)%Errinv_Adjust       = rtmp_mass(indx)%Errinv_Adjust 
     raob_mass(:)%Errinv_Final        = rtmp_mass(indx)%Errinv_Final
     raob_mass(:)%Observation         = rtmp_mass(indx)%Observation 
     raob_mass(:)%Obs_Minus_Forecast_adjusted   = rtmp_mass(indx)%Obs_Minus_Forecast_adjusted
     raob_mass(:)%Obs_Minus_Forecast_unadjusted = rtmp_mass(indx)%Obs_Minus_Forecast_unadjusted

     deallocate(indx)
     deallocate(rtmp_mass)

     call nc_diag_read_close(infn)
  end subroutine read_raob_diag_nc_mass

  function replace_text (s,text,rep)  result(outs)
     character(*)        :: s,text,rep
     character(len(s)+100) :: outs     ! provide outs with extra 100 char len
     integer             :: i, nt, nr
    
     outs = s ; nt = len_trim(text) ; nr = len_trim(rep)
     i = INDEX(outs,text(:nt)) !; IF (i == 0) EXIT
     outs = outs(:i-1) // rep(:nr) // outs(i+nt:)
  end function replace_text


end module m_diag_raob
