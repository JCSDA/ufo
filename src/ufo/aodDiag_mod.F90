module aodDiag_mod

use read_aod_diag, only: set_aoddiag,&
                     diag_header_fix_list_aod,&
                     diag_header_chan_list_aod,&
                     diag_data_name_list_aod,&
                     diag_data_chan_list_aod,&
                     diag_data_fix_list_aod

use read_aod_diag, only: set_netcdf_read_aod,&
                     read_all_aoddiag,&
                     read_aoddiag_header
                     
use nc_diag_read_mod, only: nc_diag_read_init, nc_diag_read_close

use ufo_locs_mod, only: ufo_locs, &
                        ufo_locs_setup
use fckit_log_module, only : fckit_log
use ufo_obs_data_basis_mod, only:  BasisObsData


implicit none
private

character(len=*),parameter :: myname ="aodNode_mod"

public :: aodDiag
public :: aodDiag_read

TYPE, extends(BasisObsData) :: aodDiag
  type(diag_header_fix_list_aod )              ::  header_fix
  type(diag_header_chan_list_aod),allocatable  ::  header_chan(:)
  type(diag_data_name_list_aod)                ::  header_name
  TYPE(diag_data_fix_list_aod), allocatable    ::  datafix(:)
  TYPE(diag_data_chan_list_aod) ,ALLOCATABLE   ::  datachan(:,:)
  contains
    procedure :: Setup  => read_
    procedure :: Delete => delete_
    procedure :: GetLocs => getlocs_
end type aodDiag

interface aodDiag_read  ; module procedure read_  ; end interface

contains

SUBROUTINE read_(self,filein,obstype,nobs,nlocs)
use ncd_kinds, only: i_kind
implicit none
class(aodDiag), intent(inout)  :: self
integer(i_kind), intent(inout) :: nobs
integer(i_kind), intent(inout) :: nlocs
character(len=*),intent(in)    :: filein
character(len=*),intent(in)    :: obstype

character(len=*),parameter :: myname_ =myname//"*aod_read"

integer(i_kind) :: ier
integer(i_kind) :: luin=0
logical :: lverbose  = .true.  ! control verbose

call set_netcdf_read_aod(.true.)
CALL nc_diag_read_init(filein, luin)

call read_aoddiag_header(luin,self%header_fix,self%header_chan,self%header_name,ier,lverbose)

print*, myname_, ': Found this many channels: ', self%header_fix%nchan
print*, myname_, ': Observation type in file: ', self%header_fix%obstype
print*, myname_, ': Date of input file:       ', self%header_fix%idate


CALL read_all_aoddiag ( luin, self%header_fix,self%datafix, self%datachan, nlocs, ier)

nobs  = nlocs * self%header_fix%nchan

print *, myname_, ' Total number of observations in file: ', nobs
CALL nc_diag_read_close(filename=filein)

end subroutine read_

SUBROUTINE getlocs_(self, nlocs, locs)
class(AodDiag), intent(in) :: self
TYPE(ufo_locs), intent(inout) :: locs
INTEGER, intent(in) :: nlocs

CHARACTER(len=*),PARAMETER:: myname_=myname//"*aod_getlocs"
character(len=255) :: record
integer :: failed

CALL ufo_locs_setup(locs, nlocs)

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
  WRITE(record,*)myname_,': allocated/assinged obs-data'
  call fckit_log%info(record)
else
  WRITE(record,*)myname_,': failed allocation/assignment of obs-data, ier: ', failed
  call fckit_log%info(record)
  ! should exit in error here
endif

end subroutine getlocs_

subroutine delete_(self)
implicit none
class(aodDiag), intent(inout) :: self
end subroutine delete_

subroutine write_(self)
TYPE(aodDiag), intent(inout) :: self
end subroutine write_

end module aodDiag_mod
