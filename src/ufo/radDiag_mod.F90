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
                     close_radiag, &
                     read_all_radiag
use ufo_locs_mod, only: ufo_locs, &
                        ufo_locs_setup
use fckit_log_module, only : fckit_log
use ufo_obs_data_basis_mod, only:  BasisObsData

implicit none
private

character(len=*),parameter :: myname ="radNode_mod"
integer, parameter :: max_string=800

public :: radDiag
public :: radDiag_read

type, extends(BasisObsData) :: radDiag
  type(diag_header_fix_list )              ::  header_fix
  type(diag_header_chan_list),allocatable  ::  header_chan(:)
  type(diag_data_name_list)                ::  header_name

  type(diag_data_fix_list)   ,allocatable  ::  datafix(:)
  type(diag_data_chan_list)  ,allocatable  ::  datachan(:,:)
  type(diag_data_extra_list) ,allocatable  ::  dataextra(:,:,:)
  contains
    procedure :: Setup  => read_
    procedure :: Delete => delete_
    procedure :: GetLocs => getlocs_
end type radDiag

interface radDiag_read  ; module procedure read_  ; end interface

contains

subroutine read_(self,filein,obstype,nobs,nlocs)
use ncd_kinds, only: i_kind
implicit none
class(radDiag), intent(inout)  :: self
integer(i_kind), intent(inout) :: nobs
integer(i_kind), intent(inout) :: nlocs
character(len=*),intent(in)    :: filein
character(len=*),intent(in)    :: obstype

character(len=*),parameter :: myname_ =myname//"*rad_read"
integer(i_kind) :: ier
integer(i_kind) :: luin=0
integer(i_kind) :: npred = 7   
integer(i_kind) :: iversion=30303
logical :: lverbose  = .true.  ! control verbose
logical :: retrieval = .false. ! true when dealing with SST retrievals


call set_netcdf_read(.true.)
call open_radiag(filein, luin)
call set_radiag("version",iversion,ier)

call read_radiag_header(luin,npred,retrieval,self%header_fix,self%header_chan,self%header_name,ier,lverbose)

print*, myname_, ': Found this many channels: ', self%header_fix%nchan
print*, myname_, ': Observation type in file: ', self%header_fix%obstype
print*, myname_, ': Date of input file:       ', self%header_fix%idate


call read_all_radiag(luin, self%header_fix, retrieval, self%datafix, &
                     self%datachan, self%dataextra, nobs, ier)

nlocs = nobs
nobs  = nobs * self%header_fix%nchan
call close_radiag(filein,luin)
print *, myname_, ' Total number of observations in file: (nobs,nlocs) ', nobs, nlocs
end subroutine read_


subroutine getlocs_(self, nlocs, locs)
class(RadDiag), intent(in) :: self
type(ufo_locs), intent(inout) :: locs
integer, intent(in) :: nlocs

character(len=*),parameter:: myname_=myname//"*rad_getlocs"
character(len=255) :: record
integer :: failed


call ufo_locs_setup(locs, nlocs)

failed=0
if(failed==0 .and. size(self%datafix(:)%Lat)==nlocs) then
  locs%lat(:) = self%datafix(:)%Lat
else
  failed=1
endif
if(failed==0 .and. size(self%datafix(:)%Lon)==nlocs) then
  locs%lon(:) = self%datafix(:)%Lon
else
  failed=2
endif
if(failed==0 .and. size(self%datafix(:)%obstime)==nlocs) then
  locs%time(:) = self%datafix(:)%obstime
else
  failed=3
endif
if(failed==0)then
  write(record,*)myname_,': allocated/assinged obs-data'
  call fckit_log%info(record)
else
  write(record,*)myname_,': failed allocation/assignment of obs-data, ier: ', failed
  call fckit_log%info(record)
  ! should exit in error here
endif

end subroutine getlocs_


subroutine delete_(self)
implicit none
class(radDiag), intent(inout) :: self
end subroutine delete_

subroutine write_(self)
type(radDiag), intent(inout) :: self
end subroutine write_

end module radDiag_mod
