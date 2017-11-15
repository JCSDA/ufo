module radDiag_mod

use read_diag, only: set_radiag,&
                     diag_header_fix_list,&
                     diag_header_chan_list,&
                     diag_data_name_list,&
                     read_radiag_header,&
                     set_netcdf_read
use read_diag, only: read_radiag_data,&
                     diag_data_fix_list,&
                     diag_data_extra_list,&
                     diag_data_chan_list,&
                     open_radiag, &
                     close_radiag

implicit none
private

character(len=*),parameter :: myname ="radNode_mod"
integer, parameter :: max_string=800

public :: radDiag
public :: radDiag_read

type :: radDiag
  type(diag_header_fix_list )              ::  header_fix
  type(diag_header_chan_list),allocatable  ::  header_chan(:)
  type(diag_data_name_list)                ::  header_name

  type(diag_data_fix_list)                 ::  datafix
  type(diag_data_chan_list)  ,allocatable  ::  datachan(:)
  type(diag_data_extra_list) ,allocatable  ::  dataextra(:,:)
end type radDiag

interface radDiag_read  ; module procedure this_read_  ; end interface

contains

subroutine this_read_(self,ncfname,nobs)
use ncd_kinds, only: i_kind
implicit none
character(len=*),parameter :: myname_ =myname//"*rad_read"
type(radDiag), intent(inout)  :: self
integer(i_kind), intent(inout) :: nobs
character(len=*),intent(in)    :: ncfname
integer(i_kind) :: ier
integer(i_kind) :: luin=0
integer(i_kind) :: npred = 7   
integer(i_kind) :: iversion=30303
logical :: lverbose  = .true.  ! control verbose
logical :: retrieval = .false. ! true when dealing with SST retrievals


call set_netcdf_read(.true.)
call open_radiag(ncfname, luin)
call set_radiag("version",iversion,ier)

call read_radiag_header(luin,npred,retrieval,self%header_fix,self%header_chan,self%header_name,ier,lverbose)

print*, myname_, ': Found this many channels: ', self%header_fix%nchan
print*, myname_, ': Observation type in file: ', self%header_fix%obstype
print*, myname_, ': Date of input file:       ', self%header_fix%idate

nobs=0
do while (ier .ge. 0)
   call read_radiag_data ( luin, self%header_fix, .false., self%datafix, self%datachan, &
                           self%dataextra, ier )

   if (ier .lt. 0) cycle
   nobs = nobs + 1
enddo
print *, myname_, ' Total number of observations in file: ', nobs
call close_radiag(ncfname,luin)
end subroutine this_read_

subroutine this_write_(self)
type(radDiag), intent(inout) :: self
end subroutine this_write_

end module radDiag_mod
