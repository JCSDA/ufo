program convert_and_split_conv_diag

   use ncd_kinds
   use m_diag_conv

   integer nargs, iargc, n
   character*256, allocatable ::   arg(:)

   logical ncep, append_suffix

   type(diag_conv_header) :: hdr
   type(diag_conv_mass),dimension(:),allocatable   :: mass
   type(diag_conv_wind),dimension(:),allocatable   :: wind
   integer(i_kind)   :: nobs_mass, nobs_wind
   character(120)         :: fn

   nargs = iargc()
   if( nargs.eq.0 ) then
     call usage
   else
     ncep = .false.
     append_suffix = .false.

     allocate(arg(nargs))
     do n=1,nargs
       call getarg(n,arg(n))
     enddo
     do n=1,nargs
       if (trim(arg(n)).eq.'-ncep'     ) ncep=.true.
       if (trim(arg(n)).eq.'-append_suffix') append_suffix=.true.
     enddo
   endif


   fn = arg(nargs)
   call open_conv_diag(fn, hdr, ncep)
   call read_conv_diag(fn, hdr, mass, wind, nobs_mass, nobs_wind, ncep)

   print *,hdr%ObsType
   call write_split_conv_diag_nc(fn, hdr, mass, wind, append_suffix)


end

subroutine usage
     write(6,100)
100  format( "Usage:  ",/,/ &
             "  convert_and_split_conv_diag.x <options> <filename>",/,/ &
             "where options:",/ &
             "  -ncep               :  Read NCEP (or MERRA2) diag file (default: read GMAO w/ ioff)",/ &
             "  -append_suffix      :  add '.type.nc4' suffix instead of conforming to GMAO filename standard",/ &
             "",/ &
             "  Example:",/ &
             "     convert_and_split_conv_diag.x nc_4emily_nc4.diag_conv_ges.20161202_06z.nc4",/ &
             "  Output files:",/ &
             "     nc_4emily_nc4.diag_conv_uv_ges.20161202_06z.nc4",/ &
             "     nc_4emily_nc4.diag_conv_t_ges.20161202_06z.nc4",/ &
             "     nc_4emily_nc4.diag_conv_q_ges.20161202_06z.nc4",/ &
             "     nc_4emily_nc4.diag_conv_ps_ges.20161202_06z.nc4",/ &
    )
    stop
end subroutine usage

