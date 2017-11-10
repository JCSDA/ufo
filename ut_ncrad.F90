program ut_NCrad

use ncd_kinds, only : i_kind
use read_diag, only: set_radiag,&
                     diag_header_fix_list,&
                     diag_header_chan_list,&
                     diag_data_name_list,&
                     read_radiag_header,&
                     set_netcdf_read
use read_diag, only: read_radiag_data,&
                     diag_data_fix_list,&
                     diag_data_extra_list,&
                     diag_data_chan_list
use nc_diag_read_mod, only: nc_diag_read_init, nc_diag_read_close

implicit none
character(len=*), parameter :: myname='ut_NCrad'

integer(i_kind) :: iversion = 30303
integer(i_kind) :: npred = 7   ! number of predictors in file
integer(i_kind) :: luin  = 10  ! file unit - needs to be unwired
logical :: retrieval = .false. ! true when dealing with SST retrievals
logical :: lverbose  = .true.  ! control verbose
logical :: ncftype   = .true.  ! file type NC4

integer(i_kind) iarg, argc, iargc
integer(i_kind) ier,icnt
character(len=255) :: ncfname

type(diag_header_fix_list )              ::  header_fix
type(diag_header_chan_list),allocatable  ::  header_chan(:)
type(diag_data_name_list)                ::  header_name

type(diag_data_fix_list)                 ::  datafix
type(diag_data_chan_list)  ,allocatable  ::  datachan(:)
type(diag_data_extra_list) ,allocatable  ::  dataextra(:,:)

argc = iargc()
if ( argc < 1 ) then
endif

iarg = 1
call GetArg ( iarg, ncfname )

print *, myname, ': Input file: ', trim(ncfname)

call set_netcdf_read(.true.)
call nc_diag_read_init(ncfname, luin)
call set_radiag("version",iversion,ier)

call read_radiag_header(luin,npred,retrieval,header_fix,header_chan,header_name,ier,lverbose)

icnt=0
do while (ier .ge. 0)
   call read_radiag_data ( luin, header_fix, .false., datafix, datachan, &
                           dataextra, ier )

   if (ier .lt. 0) cycle
   icnt = icnt + 1
   print *, 'icnt = ', icnt
enddo
print*, myname, ': Found this many channels: ', header_fix%nchan
print*, myname, ': Observation type in file: ', header_fix%obstype
print*, myname, ': Date of input file:       ', header_fix%idate
call nc_diag_read_close(filename=ncfname)

close(luin)
end program ut_NCrad
