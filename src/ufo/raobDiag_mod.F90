module raobDiag_mod

use m_diag_conv, only: diag_conv_header, diag_conv_mass
use m_diag_conv, only: read_conv_diag_nc_header, read_conv_diag_nc_mass

implicit none
private

character(len=*),parameter :: myname ="raobNode_mod"
integer, parameter :: max_string=800

public :: raobDiag
public :: raobDiag_read

type :: raobDiag
  type(diag_conv_header)            ::  header
  type(diag_conv_mass),allocatable  ::  mass(:)
end type raobDiag

interface raobDiag_read  ; module procedure this_read_  ; end interface

contains

subroutine this_read_(self,ncfname,nobs)
use ncd_kinds, only: i_kind
implicit none
character(len=*),parameter :: myname_ =myname//"*raod_read"
type(raobDiag),  intent(inout) :: self
integer(i_kind), intent(inout) :: nobs
character(len=*),intent(in)    :: ncfname
integer(i_kind) :: ier

call read_conv_diag_nc_header(ncfname,self%header,nobs)
allocate(self%mass(nobs))
call read_conv_diag_nc_mass(ncfname,self%header,self%mass)

print*, myname_, ': Found this many observations: ', nobs
print*, myname_, ': Date of input file:           ', self%header%date
if(nobs>0)&
print*, myname_, ': Mean observations ', sum(self%mass%Observation)/nobs

end subroutine this_read_

subroutine this_write_(self)
type(raobDiag), intent(inout) :: self
end subroutine this_write_

end module raobDiag_mod
