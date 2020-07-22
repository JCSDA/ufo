! (C) British Crown Copyright 2017-2018 Met Office
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 

!> Fortran derived type to hold configuration data for the observation covariance

module ufo_rttovonedvarcheck_rsubmatrix_mod

use kinds
use ufo_rttovonedvarcheck_constants_mod, only: max_string
use ufo_rttovonedvarcheck_rmatrix_mod

implicit none
private

type, public :: ufo_rttovonedvarcheck_rsubmatrix

  integer :: nchans !< number of channels used in current r matrix
  real(kind_real), allocatable :: matrix(:,:) !< full matrix
  real(kind_real), allocatable :: inv_matrix(:,:) !< inverse full matrix
  real(kind_real), allocatable :: diagonal(:) !< diagonal matrix
  logical :: diagonal_flag !< flag to use diagonal r-matrix
  logical :: full_flag !< flag to use full r-matrix

contains
  procedure :: setup  => rsubmatrix_setup
  procedure :: delete => rsubmatrix_delete
  procedure :: info => rsubmatrix_print
  procedure :: multiply_vector => rsubmatrix_multiply
  procedure :: multiply_matrix => rsubmatrix_multiply_matrix
  procedure :: multiply_inverse_vector => rsubmatrix_inv_multiply
  procedure :: multiply_inverse_matrix => rsubmatrix_multiply_inv_matrix
  procedure :: add_to_matrix => rsubmatrix_add_to_u

end type ufo_rttovonedvarcheck_rsubmatrix

! ------------------------------------------------------------------------------
contains
! ------------------------------------------------------------------------------
!> Setup for the r sub-matrix
!!
!! \author Met Office
!!
!! \date 09/06/2020: Created
!!
subroutine rsubmatrix_setup(self, nchans, channels, full_rmatrix)

implicit none
class(ufo_rttovonedvarcheck_rsubmatrix), intent(inout) :: self
integer, intent(in)                             :: nchans
integer, intent(in)                             :: channels(:)
type(ufo_rttovonedvarcheck_rmatrix), intent(in) :: full_rmatrix

character(len=*), parameter :: &
      routinename = "rsubmatrix_setup"
character(len=max_string)   :: err_msg
integer                     :: ii, ff, ss
character(len=max_string)   :: mat_type

self % nchans = nchans
self % full_flag = .false.
self % diagonal_flag = .false.

! Get r matrix type from full matrix
if (full_rmatrix % rtype == 1) then
  mat_type = "full"
else if (full_rmatrix % rtype == 2) then
  mat_type = "diagonal"
else
  call abor1_ftn('Unknown r matrix type')
end if

! Setup correct r matrix
select case (trim(mat_type))
   case ("full")
      ! full rmatrix not setup yet but in principle this is
      ! a start as to how it would be initialised.
      !allocate(self % matrix(nchans,nchans))
      !allocate(self % inv_matrix(nchans,nchans))
      !self % matrix(:,:) = 0.0_kind_real
      !self % inv_matrix(:,:) = 0.0_kind_real
      !self % full_flag = .true.
      !do ii=1,nchans
      !  self % matrix(ii,ii) = obs_error(ii) * obs_error(ii)
      !  self % inv_matrix(ii,ii) = 1.0_kind_real / &
      !                     (obs_error(ii) * obs_error(ii))
      !end do
      call abor1_ftn('full r matrix under development - use a diagonal')
   case ("diagonal")
      allocate(self % diagonal(nchans))
      self % diagonal(:) = 0.0_kind_real
      self % diagonal_flag = .true.

      do ff=1,full_rmatrix % nchans
        do ss=1,self % nchans
          if (full_rmatrix % channels(ff) == channels(ss)) then
            self % diagonal(ss) = full_rmatrix % errors(ff) * full_rmatrix % errors(ff)
          end if
        end do
      end do
   case default
      call abor1_ftn('Unknown r matrix type')
end select

end subroutine rsubmatrix_setup

! ------------------------------------------------------------------------------
!> Delete method for the r_matrix
!!
!! \author Met Office
!!
!! \date 09/06/2020: Created
!!
subroutine rsubmatrix_delete(self)

implicit none
class(ufo_rttovonedvarcheck_rsubmatrix), intent(inout) :: self  !< R mtrix structure

if (allocated(self % matrix))       deallocate(self % matrix)
if (allocated(self % inv_matrix))   deallocate(self % inv_matrix)
if (allocated(self % diagonal))     deallocate(self % diagonal)

end subroutine rsubmatrix_delete

! ------------------------------------------------------------------------------
!> Multiply a vector by the r-matrix
!!
!! \author Met Office
!!
!! \date 09/06/2020: Created
!!
subroutine rsubmatrix_multiply(self,xin,xout)

implicit none
class(ufo_rttovonedvarcheck_rsubmatrix), intent(in) :: self
real(kind_real), intent(in)        :: xin(:)
real(kind_real), intent(inout)     :: xout(:)

if (size(xout) /= self % nchans) then
  call abor1_ftn("rsubmatrix_multiply: arrays incompatible sizes")
end if

! Full R matrix
if (self % full_flag) xout(:) = matmul(xin(:), self % matrix(:,:))

! Diagonal R matrix
if (self % diagonal_flag) xout(:) = xin(:) * self % diagonal(:)

end subroutine rsubmatrix_multiply

! ------------------------------------------------------------------------------
!> Multiply a matrix by the r-matrix
!!
!! \author Met Office
!!
!! \date 09/06/2020: Created
!!
subroutine rsubmatrix_multiply_matrix(self,xin,xout)

implicit none
class(ufo_rttovonedvarcheck_rsubmatrix), intent(in) :: self
real(kind_real), intent(in)        :: xin(:,:)
real(kind_real), intent(inout)     :: xout(:,:)

integer :: ii

if (size(xout, 2) /= self % nchans) then
  call abor1_ftn("rsubmatrix_multiply_matrix: arrays incompatible sizes")
end if

! Full R matrix
if (self % full_flag) xout = matmul(xin, self % matrix(:,:))

! Diagonal R matrix
if (self % diagonal_flag) then
  do ii=1, self % nchans
    xout(:,ii) = xin(:,ii) * self % diagonal(ii)
  end do
end if

end subroutine rsubmatrix_multiply_matrix

! ------------------------------------------------------------------------------
!> Multiply a vector by the inverse of the r-matrix
!!
!! \author Met Office
!!
!! \date 09/06/2020: Created
!!
subroutine rsubmatrix_inv_multiply(self,xin,xout)

implicit none
class(ufo_rttovonedvarcheck_rsubmatrix), intent(in) :: self
real(kind_real), intent(in)        :: xin(:)
real(kind_real), intent(inout)     :: xout(:)

if (size(xout) /= self % nchans) then
  call abor1_ftn("rsubmatrix_inv_multiply: arrays incompatible sizes")
end if

! Full R matrix
if (self % full_flag) xout(:) = matmul(xin(:), self % inv_matrix(:,:))

! Diagonal R matrix
if (self % diagonal_flag) xout(:) = xin(:) / self % diagonal(:)

end subroutine rsubmatrix_inv_multiply

! ------------------------------------------------------------------------------
!> Multiply a matrix by the inverse of the r-matrix
!!
!! \author Met Office
!!
!! \date 09/06/2020: Created
!!
subroutine rsubmatrix_multiply_inv_matrix(self,xin,xout)

implicit none
class(ufo_rttovonedvarcheck_rsubmatrix), intent(in) :: self
real(kind_real), intent(in)        :: xin(:,:)
real(kind_real), intent(out)       :: xout(:,:)

integer :: ii

if (size(xout, 2) /= self % nchans) then
  call abor1_ftn("rsubmatrix_multiply_inv_matrix: arrays incompatible sizes")
end if

! Full R matrix
if (self % full_flag)  xout = matmul(xin, self % inv_matrix(:,:))

! Diagonal R matrix
if (self % diagonal_flag) then
  do ii=1, self % nchans
    xout(:,ii) = xin(:,ii) / self % diagonal(ii)
  end do
end if

end subroutine rsubmatrix_multiply_inv_matrix

! ------------------------------------------------------------------------------
!> Add a matrix to the r-matrix
!!
!! \author Met Office
!!
!! \date 09/06/2020: Created
!!
subroutine rsubmatrix_add_to_u(self,uin,uout)

implicit none
class(ufo_rttovonedvarcheck_rsubmatrix), intent(in) :: self
real(kind_real), intent(in)        :: uin(:,:)
real(kind_real), intent(inout)     :: uout(:,:)

integer :: ii

if (size(uout) /= self % nchans * self % nchans) then
  call abor1_ftn("rsubmatrix_add_to_u: arrays incompatible sizes")
end if

! Full R matrix
if (self % full_flag) uout = uin + self % matrix

! Diagonal R matrix
if (self % diagonal_flag) then
  do ii=1, self % nchans
    uout(ii,ii) = uin(ii,ii) + self % diagonal(ii)
  end do
end if

end subroutine rsubmatrix_add_to_u

! ------------------------------------------------------------------------------
!> Print the contents of the r-matrix
!!
!! \author Met Office
!!
!! \date 09/06/2020: Created
!!
subroutine rsubmatrix_print(self)

implicit none
class(ufo_rttovonedvarcheck_rsubmatrix), intent(in)  :: self

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

end subroutine rsubmatrix_print

end module ufo_rttovonedvarcheck_rsubmatrix_mod
