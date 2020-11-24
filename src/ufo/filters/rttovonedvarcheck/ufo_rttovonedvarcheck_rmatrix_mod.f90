! (C) British Crown Copyright 2017-2018 Met Office
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 

!> Fortran derived type to hold configuration data for the observation covariance

module ufo_rttovonedvarcheck_rmatrix_mod

use kinds
use netcdf
use ufo_rttovonedvarcheck_constants_mod, only: max_string

implicit none
private

type, public :: ufo_rttovonedvarcheck_rmatrix
  integer, allocatable         :: channels(:) !< array with instruemnt channel numbers
  integer                      :: wmo_id      !< wmo id for satellite
  integer                      :: rtype       !< type of r-matrix (1=full; 2=diagonal)
  integer                      :: nchans      !< number of channels in rmatrix
  real(kind_real), allocatable :: errors(:)   !< vector of errors

contains
  procedure :: setup  => rmatrix_setup
  procedure :: delete => rmatrix_delete

end type ufo_rttovonedvarcheck_rmatrix

! ------------------------------------------------------------------------------
contains
! ------------------------------------------------------------------------------
!> Setup for the full r_matrix
!!
!! \author Met Office
!!
!! \date 09/06/2020: Created
!!
subroutine rmatrix_setup(self, filename)

implicit none
class(ufo_rttovonedvarcheck_rmatrix), intent(inout) :: self !< Full R matrix structure
character(len=*), intent(in)  :: filename !< Path to input filename

integer :: error  ! error code for read
integer :: lun    ! value for identifying file
integer :: dimid  ! value for dimension id
integer :: varid  ! value for variable id
integer :: nchans ! number of channels

! defaults
self % nchans = 0

! Open netcdf file
error = nf90_open(trim(filename), nf90_nowrite, lun)
if (error /= 0) call abor1_ftn("error in opening file")

! Read wmo_id
error = nf90_inq_varid(lun, "wmo_id", varid)
error = nf90_get_var(lun, varid, self % wmo_id)
if (error /= 0) call abor1_ftn("error in reading the WMO ID")

! Read rtype
error = nf90_inq_varid(lun, "r_type", varid)
error = nf90_get_var(lun, varid, self % rtype)
if (error /= 0) call abor1_ftn("error in reading the rmatrix type")

! Get dimensions of arrays
error = nf90_inq_dimid(lun, "nchans", dimid)
error = nf90_inquire_dimension(lun, dimid, len=self % nchans)
if (error /= 0) call abor1_ftn("error in reading the dimensions from numchans")

! Allocate arrays
allocate(self % channels(self % nchans))
allocate(self % errors(self % nchans))

! Read channels from files
error = nf90_inq_varid(lun, "channels", varid)
error = nf90_get_var(lun, varid, self % channels)
if (error /= 0) call abor1_ftn("error in reading the rmatrix channels")

! Read channels from files
error = nf90_inq_varid(lun, "obs_error", varid)
error = nf90_get_var(lun, varid, self % errors)
if (error /= 0) call abor1_ftn("error in reading the rmatrix errors")

! Close netcdf file
error = nf90_close(lun)
if (error /= 0) call abor1_ftn("error in closing the rmatrix")

write(*,*) "Successfully opened and close r matrix netcdf file"
call rmatrix_print(self)

end subroutine rmatrix_setup

! ------------------------------------------------------------------------------
!> Delete method for the full r_matrix
!!
!! \author Met Office
!!
!! \date 09/06/2020: Created
!!
subroutine rmatrix_delete(self)

implicit none
class(ufo_rttovonedvarcheck_rmatrix), intent(inout) :: self  !< Full R matrix structure

if (allocated(self % channels)) deallocate(self % channels)
if (allocated(self % errors))   deallocate(self % errors)

end subroutine rmatrix_delete

! ------------------------------------------------------------------------------
!> Print method for the full r_matrix
!!
!! \author Met Office
!!
!! \date 09/06/2020: Created
!!
subroutine rmatrix_print(self)

implicit none
class(ufo_rttovonedvarcheck_rmatrix), intent(inout) :: self  !< Full R matrix structure

write(*,*) "wmo_id = ", self % wmo_id
write(*,*) "rtype = ", self % rtype
write(*,*) "channels = ", self % channels(:)
write(*,*) "errors = ", self % errors(:)

end subroutine rmatrix_print

end module ufo_rttovonedvarcheck_rmatrix_mod
