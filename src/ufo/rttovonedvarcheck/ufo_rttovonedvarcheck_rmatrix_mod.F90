! (C) British Crown Copyright 2017-2018 Met Office
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 

module ufo_rttovonedvarcheck_rmatrix_mod

use kinds

implicit none

!> Fortran derived type to hold configuration data for the observation covariance
type :: rmatrix_type

  integer :: nchans
  real(kind_real), allocatable :: matrix(:,:)
  real(kind_real), allocatable :: inv_matrix(:,:)
  real(kind_real), allocatable :: diagonal(:)
  logical :: diagonal_flag
  logical :: full_flag

end type rmatrix_type

! ------------------------------------------------------------------------------
contains
! ------------------------------------------------------------------------------

!> Setup for the r_matrix

subroutine rmatrix_setup(self,mat_type,nchans,obs_error)

implicit none
type(rmatrix_type), intent(inout) :: self
character(len=*), intent(in)      :: mat_type
integer, intent(in)               :: nchans
real(kind_real), intent(in)       :: obs_error(nchans)

character(len=*), parameter :: &
      routinename = "ufo_rttovonedvarcheck_InitRmatrix"
integer :: ii

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
      write(*,*) "Unknown r matrix type"
      stop
end select

end subroutine rmatrix_setup

! ------------------------------------------------------------------------------

subroutine rmatrix_delete(self)

implicit none
type(rmatrix_type), intent(inout) :: self  !< R mtrix structure

if (allocated(self % matrix))       deallocate(self % matrix)
if (allocated(self % inv_matrix))   deallocate(self % inv_matrix)
if (allocated(self % diagonal))     deallocate(self % diagonal)

end subroutine rmatrix_delete

! ------------------------------------------------------------------------------

subroutine rmatrix_multiply(self,xin,xout)

implicit none
type(rmatrix_type), intent(in)  :: self
real(kind_real), intent(in)     :: xin(:)
real(kind_real), intent(inout)  :: xout(:)

if (size(xout) /= self % nchans) then
  write(*,*) "R-matrix and increment are not the same size stopping"
  write(*,*) "R-matrix and increment sizes are: ",self % nchans, &
                                                  size(xin), &
                                                  size(xout)
  stop
end if

! Full R matrix
if (self % full_flag) xout(:) = matmul(xin(:), self % matrix(:,:))

! Diagonal R matrix
if (self % diagonal_flag) xout(:) = xin(:) * self % diagonal(:)

end subroutine rmatrix_multiply

! ------------------------------------------------------------------------------

subroutine rmatrix_multiply_matrix(self,xin,xout)

implicit none
type(rmatrix_type), intent(in)  :: self
real(kind_real), intent(in)     :: xin(:,:)
real(kind_real), intent(inout)  :: xout(:,:)

integer :: ii

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

subroutine rmatrix_inv_multiply(self,xin,xout)

implicit none
type(rmatrix_type), intent(in)  :: self
real(kind_real), intent(in)     :: xin(:)
real(kind_real), intent(inout)  :: xout(:)

if (size(xout) /= self % nchans) then
  write(*,*) "R-matrix and increment are not the same size stopping"
  write(*,*) "R-matrix and increment sizes are: ",self % nchans, &
                                                  size(xin), &
                                                  size(xout)
  stop
end if

! Full R matrix
if (self % full_flag) xout(:) = matmul(xin(:), self % inv_matrix(:,:))

! Diagonal R matrix
if (self % diagonal_flag) xout(:) = xin(:) / self % diagonal(:)

end subroutine rmatrix_inv_multiply

! ------------------------------------------------------------------------------

subroutine rmatrix_multiply_inv_matrix(self,xin,xout)

implicit none
type(rmatrix_type), intent(in) :: self
real(kind_real), intent(in)    :: xin(:,:)
real(kind_real), intent(out)   :: xout(:,:)

integer :: ii

! Full R matrix
write(*,*) "xin = ",xin
write(*,*) "self % inv_matrix(:,:) = ",self % inv_matrix(:,:)
if (self % full_flag)  xout = matmul(xin, self % inv_matrix(:,:))

! Diagonal R matrix
if (self % diagonal_flag) then
  do ii=1, self % nchans
    xout(:,ii) = xin(:,ii) / self % diagonal(ii)
  end do
end if

end subroutine rmatrix_multiply_inv_matrix

! -----------------------------------------------------------------------------

subroutine rmatrix_add_to_u(self,uin,uout)

implicit none
type(rmatrix_type), intent(in)  :: self
real(kind_real), intent(in)     :: uin(:,:)
real(kind_real), intent(inout)  :: uout(:,:)

integer :: ii

write(*,*) "Using rmatrix_add_to_u ",self % full_flag,self % diagonal_flag

if (size(uout) /= self % nchans*self % nchans) then
  write(*,*) "R-matrix and U are not the same size stopping"
  write(*,*) "R-matrix and U sizes are: ",self % nchans, &
                                          size(uin), &
                                          size(uout)
  stop
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

subroutine rmatrix_print(self)

implicit none
type(rmatrix_type), intent(in)  :: self

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
