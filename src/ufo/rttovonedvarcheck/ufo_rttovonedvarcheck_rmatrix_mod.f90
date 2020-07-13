! (C) British Crown Copyright 2017-2018 Met Office
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 

!> Fortran derived type to hold configuration data for the observation covariance

module ufo_rttovonedvarcheck_rmatrix_mod

use kinds
use ufo_rttovonedvarcheck_constants_mod, only: max_string

implicit none
private

type, public :: rmatrix_type

  integer :: nchans !< number of channels used in current r matrix
  real(kind_real), allocatable :: matrix(:,:) !< full matrix
  real(kind_real), allocatable :: inv_matrix(:,:) !< inverse full matrix
  real(kind_real), allocatable :: diagonal(:) !< diagonal matrix
  logical :: diagonal_flag !< flag to use diagonal r-matrix
  logical :: full_flag !< flag to use full r-matrix

contains
  procedure :: setup  => rmatrix_setup
  procedure :: delete => rmatrix_delete
  procedure :: info => rmatrix_print
  procedure :: multiply_vector => rmatrix_multiply
  procedure :: multiply_matrix => rmatrix_multiply_matrix
  procedure :: multiply_inverse_vector => rmatrix_inv_multiply
  procedure :: multiply_inverse_matrix => rmatrix_multiply_inv_matrix
  procedure :: add_to_matrix => rmatrix_add_to_u

end type rmatrix_type

! ------------------------------------------------------------------------------
contains
! ------------------------------------------------------------------------------
!> Setup for the r_matrix
!!
!! \author Met Office
!!
!! \date 09/06/2020: Created
!!
subroutine rmatrix_setup(self,mat_type,nchans,obs_error)

implicit none
class(rmatrix_type), intent(inout) :: self
character(len=*), intent(in)      :: mat_type
integer, intent(in)               :: nchans
real(kind_real), intent(in)       :: obs_error(nchans)

character(len=*), parameter :: &
      routinename = "rmatrix_setup"
character(len=max_string)   :: err_msg
integer                     :: ii

self % nchans = nchans
self % full_flag = .false.
self % diagonal_flag = .false.

select case (trim(mat_type))
   case ("full")
      allocate(self % matrix(nchans,nchans))
      allocate(self % inv_matrix(nchans,nchans))
      self % matrix(:,:) = 0.0_kind_real
      self % inv_matrix(:,:) = 0.0_kind_real
      self % full_flag = .true.
      do ii=1,nchans
        self % matrix(ii,ii) = obs_error(ii) * obs_error(ii)
        self % inv_matrix(ii,ii) = 1.0_kind_real / &
                           (obs_error(ii) * obs_error(ii))
      end do
   case ("diagonal")
      allocate(self % diagonal(nchans))
      self % diagonal(:) = 0.0_kind_real
      self % diagonal_flag = .true.
      self % diagonal(:) = obs_error(:) * obs_error(:)
   case default
      write(err_msg,*) 'Unknown r matrix type'
      call abor1_ftn(err_msg)
end select

end subroutine rmatrix_setup

! ------------------------------------------------------------------------------
!> Delete method for the r_matrix
!!
!! \author Met Office
!!
!! \date 09/06/2020: Created
!!
subroutine rmatrix_delete(self)

implicit none
class(rmatrix_type), intent(inout) :: self  !< R mtrix structure

if (allocated(self % matrix))       deallocate(self % matrix)
if (allocated(self % inv_matrix))   deallocate(self % inv_matrix)
if (allocated(self % diagonal))     deallocate(self % diagonal)

end subroutine rmatrix_delete

! ------------------------------------------------------------------------------
!> Multiply a vector by the r-matrix
!!
!! \author Met Office
!!
!! \date 09/06/2020: Created
!!
subroutine rmatrix_multiply(self,xin,xout)

implicit none
class(rmatrix_type), intent(in)  :: self
real(kind_real), intent(in)     :: xin(:)
real(kind_real), intent(inout)  :: xout(:)

if (size(xout) /= self % nchans) then
  call abor1_ftn("rmatrix_multiply: arrays incompatible sizes")
end if

! Full R matrix
if (self % full_flag) xout(:) = matmul(xin(:), self % matrix(:,:))

! Diagonal R matrix
if (self % diagonal_flag) xout(:) = xin(:) * self % diagonal(:)

end subroutine rmatrix_multiply

! ------------------------------------------------------------------------------
!> Multiply a matrix by the r-matrix
!!
!! \author Met Office
!!
!! \date 09/06/2020: Created
!!
subroutine rmatrix_multiply_matrix(self,xin,xout)

implicit none
class(rmatrix_type), intent(in)  :: self
real(kind_real), intent(in)     :: xin(:,:)
real(kind_real), intent(inout)  :: xout(:,:)

integer :: ii

if (size(xout, 2) /= self % nchans) then
  call abor1_ftn("rmatrix_multiply_matrix: arrays incompatible sizes")
end if

! Full R matrix
if (self % full_flag) xout = matmul(xin, self % matrix(:,:))

! Diagonal R matrix
if (self % diagonal_flag) then
  do ii=1, self % nchans
    xout(:,ii) = xin(:,ii) * self % diagonal(ii)
  end do
end if

end subroutine rmatrix_multiply_matrix

! ------------------------------------------------------------------------------
!> Multiply a vector by the inverse of the r-matrix
!!
!! \author Met Office
!!
!! \date 09/06/2020: Created
!!
subroutine rmatrix_inv_multiply(self,xin,xout)

implicit none
class(rmatrix_type), intent(in)  :: self
real(kind_real), intent(in)     :: xin(:)
real(kind_real), intent(inout)  :: xout(:)

if (size(xout) /= self % nchans) then
  call abor1_ftn("rmatrix_inv_multiply: arrays incompatible sizes")
end if

! Full R matrix
if (self % full_flag) xout(:) = matmul(xin(:), self % inv_matrix(:,:))

! Diagonal R matrix
if (self % diagonal_flag) xout(:) = xin(:) / self % diagonal(:)

end subroutine rmatrix_inv_multiply

! ------------------------------------------------------------------------------
!> Multiply a matrix by the inverse of the r-matrix
!!
!! \author Met Office
!!
!! \date 09/06/2020: Created
!!
subroutine rmatrix_multiply_inv_matrix(self,xin,xout)

implicit none
class(rmatrix_type), intent(in) :: self
real(kind_real), intent(in)    :: xin(:,:)
real(kind_real), intent(out)   :: xout(:,:)

integer :: ii

if (size(xout, 2) /= self % nchans) then
  call abor1_ftn("rmatrix_multiply_inv_matrix: arrays incompatible sizes")
end if

! Full R matrix
if (self % full_flag)  xout = matmul(xin, self % inv_matrix(:,:))

! Diagonal R matrix
if (self % diagonal_flag) then
  do ii=1, self % nchans
    xout(:,ii) = xin(:,ii) / self % diagonal(ii)
  end do
end if

end subroutine rmatrix_multiply_inv_matrix

! ------------------------------------------------------------------------------
!> Add a matrix to the r-matrix
!!
!! \author Met Office
!!
!! \date 09/06/2020: Created
!!
subroutine rmatrix_add_to_u(self,uin,uout)

implicit none
class(rmatrix_type), intent(in)  :: self
real(kind_real), intent(in)     :: uin(:,:)
real(kind_real), intent(inout)  :: uout(:,:)

integer :: ii

if (size(uout) /= self % nchans * self % nchans) then
  call abor1_ftn("rmatrix_add_to_u: arrays incompatible sizes")
end if

! Full R matrix
if (self % full_flag) uout = uin + self % matrix

! Diagonal R matrix
if (self % diagonal_flag) then
  do ii=1, self % nchans
    uout(ii,ii) = uin(ii,ii) + self % diagonal(ii)
  end do
end if

end subroutine rmatrix_add_to_u

! ------------------------------------------------------------------------------
!> Print the contents of the r-matrix
!!
!! \author Met Office
!!
!! \date 09/06/2020: Created
!!
subroutine rmatrix_print(self)

implicit none
class(rmatrix_type), intent(in)  :: self

integer :: ii

if (self % full_flag) then

  write(*,*) "Full R matrix used printing diagonal"
  write(*,*) "nchans = ",self % nchans
  write(*,*) "Matrix diagonal elements = "
  do ii = 1, self % nchans
    write(*,*) self % matrix(ii,ii)
  end do
  write(*,*) "Inverse Matrix diagonal elements = "
  do ii = 1, self % nchans
    write(*,*) self % inv_matrix(ii,ii)
  end do
  
end if

if (self % diagonal_flag) then

  write(*,*) "Diagonal R matrix used"
  write(*,*) "nchans = ",self % nchans
  write(*,*) "Diagonal = ",self % diagonal(:)

end if

end subroutine rmatrix_print

end module ufo_rttovonedvarcheck_rmatrix_mod
