!-------------------------------------------------------------------------------
! (C) British Crown Copyright 2020 Met Office
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
!-------------------------------------------------------------------------------

module ufo_gnssroonedvarcheck_utils_mod

use missing_values_mod
use kinds

implicit none
private
public :: singleob_type, singlebg_type
public :: find_unique
public :: Ops_RealSortQuick
public :: allocate_singleob, deallocate_singleob
public :: allocate_singlebg, deallocate_singlebg

! Add the ability to hold various data and meta-data for a variable
type element_type
  real(kind_real) :: value      ! Observed value
  real(kind_real) :: OBErr      ! Observation error
  real(kind_real) :: PGEFinal   ! Probability of gross error
end type element_type

! Structure for the observation information
type :: singleob_type
  integer                          :: niter
  integer                          :: id
  real(kind_real)                  :: jcost
  real(kind_real)                  :: latitude
  real(kind_real)                  :: longitude
  type (element_type), allocatable :: p(:)
  type (element_type), allocatable :: q(:)
  real(kind_real), allocatable     :: solutbendingangle(:)
  type (element_type), allocatable :: bendingangle(:)
  type (element_type), allocatable :: impactparam(:)
  type (element_type)              :: ro_rad_curv
  type (element_type)              :: ro_geoid_und
  integer, allocatable             :: qc_flags(:)
end type

! Structure for the background (model) information
type :: singlebg_type
  real(kind_real), allocatable :: za(:)
  real(kind_real), allocatable :: zb(:)
  real(kind_real), allocatable :: p(:)
  real(kind_real), allocatable :: q(:)
end type

contains

!------------------------------------------------------------------------------
!> Find the unique entries in the input list
!!
!! \author Met Office
!!
!! \date 15/10/2020: Created
!!
subroutine find_unique(input, output)

integer, intent(in) :: input(:)
integer, allocatable, intent(out) :: output(:)

integer, allocatable :: unique_vals(:)    ! The list of unique elements
integer :: nfound = 0                     ! The number of unique elements found
integer :: cur_val                        ! The current value being considered
integer :: max_val                        ! The maximum value in the input list

allocate(unique_vals(1:size(input)))

cur_val = minval(input) - 1
max_val = maxval(input)
do while (cur_val < max_val)
    nfound = nfound + 1
    cur_val = minval(input, mask=input>cur_val)
    unique_vals(nfound) = cur_val
end do
allocate(output(nfound))
output = unique_vals(1:nfound)

end subroutine find_unique

!------------------------------------------------------------------------------
!> Generates a index array pointing to the elements of the array 'key'
!  in increasing order
!!
!! \author Met Office
!!
!! \date 15/10/2020: Created
!!
!
! Method:
!   The heap sort invented by J.W.J.Williams is used.
!   A description of the method can be found in 'Numerical Recipes'
!   The group information array is used to allow easy sorting on several
!   parameters of different types. For details see the Parent Module
!   OpsMod_Sort
!
! Inputs:
!   key : An array of character strings, to be sorted
!
! Input/Output:
!   index : An integer array pointing to the sorted items.
!-------------------------------------------------------------------------------
SUBROUTINE Ops_RealSortQuick(key,   &
                             index)

IMPLICIT NONE

! Subroutine arguments:
REAL(kind_real), INTENT(IN)       :: key(:)
INTEGER, ALLOCATABLE, INTENT(OUT) :: index(:)

! Local declarations:
INTEGER                           :: n     ! The number of items
INTEGER                           :: head  ! heaps are tree structures: head and child refer
INTEGER                           :: child ! to related items within the tree
INTEGER                           :: j
INTEGER                           :: dum   ! used to swap index items
CHARACTER(len=*), PARAMETER       :: RoutineName = 'Ops_RealSortQuick'

! Could put in an optional mask

n = SIZE (key)
ALLOCATE (Index(n))
DO j = 1, n
  Index(j) = j
END DO

! Do heapsort: Create the heap...

makeheap : DO j = n / 2, 1, -1
  head = j
  sift1 : DO

    ! find the largest out of the head and its two children...

    child = head * 2
    IF (child > n) EXIT sift1
    IF (child < n) THEN
      IF (key(Index(child + 1)) > key(Index(child))) child = child + 1
    END IF

    ! if the head is the largest, then sift is done...

    IF (key(Index(head)) >= key(Index(child))) EXIT sift1

    ! otherwise swap to put the largest child at the head,
    ! and prepare to repeat the procedure for the head in its new
    ! subordinate position.

    dum = Index(child)
    Index(child) = Index(head)
    Index(head) = dum
    head = child
  END DO sift1
END DO makeheap

! Retire heads of the heap, which are the largest, and
! stack them at the end of the array.

retire : DO j = n, 2, -1
  dum = Index(1)
  Index(1) = Index(j)
  Index(j) = dum
  head = 1

  ! second sift is similar to first...

  sift2: DO
    child = head * 2
    IF (child > (j - 1)) EXIT sift2
    IF (child < (j - 1)) THEN
      IF (key(Index(child + 1)) > key(Index(child))) child = child + 1
    END IF
    IF (key(Index(head)) >= key(Index(child))) EXIT sift2
    dum = Index(child)
    Index(child) = Index(head)
    Index(head) = dum
    head = child
  END DO sift2
END DO retire

END SUBROUTINE Ops_RealSortQuick


!------------------------------------------------------------------------------
!> Allocate the singleob_type structure, given a certain number of observations,
!  and model levels for pressure and specific humidity.
!!
!! \author Met Office
!!
!! \date 15/10/2020: Created
!!
subroutine allocate_singleob(singleob, nobs, nlevp, nlevq)

implicit none

type(singleob_type), intent(out) :: singleob
integer, intent(in) :: nobs
integer, intent(in) :: nlevp
integer, intent(in) :: nlevq

allocate(singleob % p(1:nlevp))
allocate(singleob % q(1:nlevq))
allocate(singleob % solutbendingangle(1:nobs))
allocate(singleob % bendingangle(1:nobs))
allocate(singleob % impactparam(1:nobs))
allocate(singleob % qc_flags(1:nobs))

singleob % solutbendingangle(:) = missing_value(singleob % solutbendingangle(1))
singleob % qc_flags(:) = 0  ! Set to unflagged

singleob % p(:) % value = missing_value(singleob % p(1) % value)
singleob % p(:) % oberr = missing_value(singleob % p(1) % oberr)
singleob % p(:) % pgefinal = missing_value(singleob % p(1) % pgefinal)

singleob % q(:) % value = missing_value(singleob % q(1) % value)
singleob % q(:) % oberr = missing_value(singleob % q(1) % oberr)
singleob % q(:) % pgefinal = missing_value(singleob % q(1) % pgefinal)

singleob % bendingangle(:) % value = missing_value(singleob % bendingangle(1) % value)
singleob % bendingangle(:) % oberr = missing_value(singleob % bendingangle(1) % oberr)
singleob % bendingangle(:) % pgefinal = missing_value(singleob % bendingangle(1) % pgefinal)

singleob % impactparam(:) % value = missing_value(singleob % impactparam(1) % value)
singleob % impactparam(:) % oberr = missing_value(singleob % impactparam(1) % oberr)
singleob % impactparam(:) % pgefinal = missing_value(singleob % impactparam(1) % pgefinal)

end subroutine allocate_singleob


!------------------------------------------------------------------------------
!> Deallocate the singleob_type structure
!!
!! \author Met Office
!!
!! \date 15/10/2020: Created
!!
subroutine deallocate_singleob(singleob)

implicit none

type(singleob_type), intent(inout) :: singleob

deallocate(singleob % p)
deallocate(singleob % q)
deallocate(singleob % solutbendingangle)
deallocate(singleob % bendingangle)
deallocate(singleob % impactparam)
deallocate(singleob % qc_flags)

end subroutine deallocate_singleob


!------------------------------------------------------------------------------
!> Allocate the structure to hold background information from a single profile.
!!
!! \author Met Office
!!
!! \date 15/10/2020: Created
!!
subroutine allocate_singlebg(singlebg, nlevp, nlevq)

implicit none

type(singlebg_type), intent(out) :: singlebg
integer, intent(in) :: nlevp
integer, intent(in) :: nlevq

allocate(singlebg % za(1:nlevp))
allocate(singlebg % zb(1:nlevq))
allocate(singlebg % p(1:nlevp))
allocate(singlebg % q(1:nlevq))

singlebg % za(:) = missing_value(singlebg % za(1))
singlebg % zb(:) = missing_value(singlebg % zb(1))
singlebg % p(:) = missing_value(singlebg % p(1))
singlebg % q(:) = missing_value(singlebg % q(1))

end subroutine allocate_singlebg


!------------------------------------------------------------------------------
!> Dealloate the singlebg_type structure
!!
!! \author Met Office
!!
!! \date 15/10/2020: Created
!!
subroutine deallocate_singlebg(singlebg)

implicit none

type(singlebg_type), intent(inout) :: singlebg

deallocate(singlebg % za)
deallocate(singlebg % zb)
deallocate(singlebg % p)
deallocate(singlebg % q)

end subroutine deallocate_singlebg

! ------------------------------------------------------------------------------
end module ufo_gnssroonedvarcheck_utils_mod
