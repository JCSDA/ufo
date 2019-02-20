program ut_NCrad

use ncd_kinds, only : i_kind
use m_diag_conv, only: diag_conv_header, diag_conv_mass
use m_diag_conv, only: read_conv_diag_nc_header, read_conv_diag_nc_mass

use m_diag_raob, only: diag_raob_header, diag_raob_mass
use m_diag_raob, only: read_raob_diag_nc_header, read_raob_diag_nc_mass

implicit none
character(len=*), parameter :: myname='ut_NCmass'

integer(i_kind) iarg, argc, iargc
integer(i_kind) ier,icnt,nobs
character(len=255) :: ncfname

type(diag_conv_header)               ::  conv_header
type(diag_conv_mass)       ,pointer  ::  conv_mass(:)
type(diag_raob_header)               ::  raob_header
type(diag_raob_mass)       ,pointer  ::  raob_mass(:)

argc = iargc()
if ( argc < 1 ) then
endif

iarg = 1
call GetArg ( iarg, ncfname )
print *, myname, ': Input file: ', trim(ncfname)

call read_conv_diag_nc_header(ncfname,conv_header,nobs)
allocate(conv_mass(nobs))
call read_conv_diag_nc_mass(ncfname,conv_header,conv_mass)
print*, myname, ': Found this many observations: ', nobs
print*, myname, ': Date of input file:           ', conv_header%date
deallocate(conv_mass)

call read_raob_diag_nc_header(ncfname,raob_header)
allocate(raob_mass(raob_header%n_Observations_Mass))
call read_raob_diag_nc_mass(ncfname,raob_header,raob_mass,ier)
print*, myname, ': Found this many RAOB observations: ', raob_header%n_Observations_Mass
print*, myname, ': Date of input RAOB file:           ', raob_header%date
print*, myname, ': Size of type holding RAOB:         ', size(raob_mass)
deallocate(raob_mass)

end program ut_NCrad
