module raobDiag_mod

use m_diag_raob, only: diag_raob_header, diag_raob_mass
use m_diag_raob, only: read_raob_diag_nc_header, read_raob_diag_nc_mass

implicit none
private

character(len=*),parameter :: myname ="raobNode_mod"
integer, parameter :: max_string=800

public :: raobDiag
public :: raobDiag_read

type :: raobDiag
  type(diag_raob_header)        ::  header
  type(diag_raob_mass),pointer  ::  mass(:)
end type raobDiag

interface raobDiag_read  ; module procedure this_read_  ; end interface

contains

subroutine this_read_(self,ncfname,nobs,nlocs)
use ncd_kinds, only: i_kind
implicit none
character(len=*),parameter :: myname_ =myname//"*raod_read"
type(raobDiag),  intent(inout) :: self
integer(i_kind), intent(inout) :: nobs, nlocs
character(len=*),intent(in)    :: ncfname
integer(i_kind) :: ier

call read_raob_diag_nc_header(ncfname,self%header)
nobs=self%header%n_Observations_Mass
allocate(self%mass(nobs))
call read_raob_diag_nc_mass(ncfname,self%header,self%mass,ier)
nobs=self%header%n_Observations_Mass
nlocs = nobs

print*, myname_, ': Found this many observations: ', nobs
print*, myname_, ': Size of type holding RAOB:    ', size(self%mass)
print*, myname_, ': Date of input file:           ', self%header%date
if(nobs>0)&
print*, myname_, ': Mean observations:            ', sum(self%mass%Observation)/nobs

end subroutine this_read_

subroutine this_write_(self)
type(raobDiag), intent(inout) :: self
end subroutine this_write_

end module raobDiag_mod
