module raobDiag_mod

use m_diag_raob, only: diag_raob_header, diag_raob_mass
use m_diag_raob, only: read_raob_diag_nc_header, read_raob_diag_nc_mass

use ufo_obs_data_basis_mod, only:  BasisObsData

implicit none
private

character(len=*),parameter :: myname ="raobNode_mod"
integer, parameter :: max_string=800

public :: raobDiag
public :: raobDiag_read

type, extends(BasisObsData) :: raobDiag
  type(diag_raob_header)        ::  header
  type(diag_raob_mass),pointer  ::  mass(:)
  contains
    procedure :: Setup => this_read_
    procedure :: Delete => this_delete_
end type raobDiag

interface raobDiag_read  ; module procedure this_read_  ; end interface

contains

subroutine this_read_(self,filein,obstype,nobs,nlocs)
use ncd_kinds, only: i_kind
implicit none
character(len=*),parameter :: myname_ =myname//"*raod_read"
class(raobDiag),  intent(inout) :: self
integer(i_kind), intent(inout) :: nobs
integer(i_kind), intent(inout) :: nlocs
character(len=*),intent(in)    :: filein
character(len=*),intent(in)    :: obstype
integer(i_kind) :: ier

call read_raob_diag_nc_header(filein,self%header)
nobs=self%header%n_Observations_Mass
allocate(self%mass(nobs))
call read_raob_diag_nc_mass(filein,self%header,self%mass,ier)
nobs=self%header%n_Observations_Mass
nlocs = nobs

print*, myname_, ': Found this many observations: ', nobs
print*, myname_, ': Found this many instances:    ', nlocs
print*, myname_, ': Size of type holding RAOB:    ', size(self%mass)
print*, myname_, ': Date of input file:           ', self%header%date
if(nobs>0)&
print*, myname_, ': Mean observations:            ', sum(self%mass(:)%Observation)/nobs
print*, self%mass(1:20)%Observation

end subroutine this_read_

subroutine this_delete_(self)
implicit none
class(raobDiag), intent(inout) :: self
end subroutine this_delete_

subroutine this_write_(self)
class(raobDiag), intent(inout) :: self
end subroutine this_write_

end module raobDiag_mod
