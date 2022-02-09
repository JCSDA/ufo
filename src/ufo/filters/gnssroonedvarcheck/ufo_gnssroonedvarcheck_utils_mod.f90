!-------------------------------------------------------------------------------
! (C) British Crown Copyright 2020 Met Office
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
!-------------------------------------------------------------------------------

module ufo_gnssroonedvarcheck_utils_mod

use, intrinsic :: iso_c_binding
use missing_values_mod
use kinds

implicit none
private
public :: singleob_type, singlebg_type
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
