!-------------------------------------------------------------------------------
! (C) British Crown Copyright 2025 Met Office
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
!-------------------------------------------------------------------------------

module ufo_refractivityonedvarcheck_utils_mod

use fckit_exception_module, only: fckit_exception
use, intrinsic :: iso_c_binding
use missing_values_mod, only: missing_value
use kinds

implicit none
private
public :: singlerefob_type, singlerefbg_type
public :: allocate_singlerefob, deallocate_singlerefob
public :: allocate_singlerefbg, deallocate_singlerefbg

! Add the ability to hold various data and meta-data for a variable
type element_type
  real(kind_real) :: value      ! Observed value
  real(kind_real) :: OBErr      ! Observation error
  real(kind_real) :: PGEFinal   ! Probability of gross error
end type element_type

! Structure for the observation information
type :: singlerefob_type
  integer                          :: niter
  integer                          :: id
  real(kind_real)                  :: jcost
  real(kind_real)                  :: latitude
  real(kind_real)                  :: longitude
  type (element_type), allocatable :: p(:)
  type (element_type), allocatable :: q(:)
  real(kind_real), allocatable     :: solutrefractivity(:)
  type (element_type), allocatable :: refractivity(:)
  type (element_type), allocatable :: z(:)
  integer, allocatable             :: qc_flags(:)
end type

! Structure for the background (model) information
type :: singlerefbg_type
  real(kind_real), allocatable :: za(:)
  real(kind_real), allocatable :: zb(:)
  real(kind_real), allocatable :: p(:)
  real(kind_real), allocatable :: q(:)
end type

contains

!------------------------------------------------------------------------------
!> Allocate the singlerefob_type structure, given a certain number of observations,
!  and model levels for pressure and specific humidity.
!!
!! \author Met Office
!!
!! \date 13/05/2025: Created
!!
subroutine allocate_singlerefob(singleob, nobs, nlevp, nlevq)

implicit none

type(singlerefob_type), intent(out) :: singleob
integer, intent(in) :: nobs
integer, intent(in) :: nlevp
integer, intent(in) :: nlevq

allocate(singleob % p(1:nlevp))
allocate(singleob % q(1:nlevq))
allocate(singleob % solutrefractivity(1:nobs))
allocate(singleob % refractivity(1:nobs))
allocate(singleob % z(1:nobs))
allocate(singleob % qc_flags(1:nobs))

singleob % solutrefractivity(:) = missing_value(singleob % solutrefractivity(1))
singleob % qc_flags(:) = 0  ! Set to unflagged

singleob % p(:) % value = missing_value(singleob % p(1) % value)
singleob % p(:) % oberr = missing_value(singleob % p(1) % oberr)
singleob % p(:) % pgefinal = missing_value(singleob % p(1) % pgefinal)

singleob % q(:) % value = missing_value(singleob % q(1) % value)
singleob % q(:) % oberr = missing_value(singleob % q(1) % oberr)
singleob % q(:) % pgefinal = missing_value(singleob % q(1) % pgefinal)

singleob % refractivity(:) % value = missing_value(singleob % refractivity(1) % value)
singleob % refractivity(:) % oberr = missing_value(singleob % refractivity(1) % oberr)
singleob % refractivity(:) % pgefinal = missing_value(singleob % refractivity(1) % pgefinal)

singleob % z(:) % value = missing_value(singleob % z(1) % value)
singleob % z(:) % oberr = missing_value(singleob % z(1) % oberr)
singleob % z(:) % pgefinal = missing_value(singleob % z(1) % pgefinal)

end subroutine allocate_singlerefob


!------------------------------------------------------------------------------
!> Deallocate the singlerefob_type structure
!!
!! \author Met Office
!!
!! \date 13/05/2025: Created
!!
subroutine deallocate_singlerefob(singleob)

implicit none

type(singlerefob_type), intent(inout) :: singleob

deallocate(singleob % p)
deallocate(singleob % q)
deallocate(singleob % solutrefractivity)
deallocate(singleob % refractivity)
deallocate(singleob % z)
deallocate(singleob % qc_flags)

end subroutine deallocate_singlerefob


!------------------------------------------------------------------------------
!> Allocate the structure to hold background information from a single profile.
!!
!! \author Met Office
!!
!! \date 13/05/2025: Created
!!
subroutine allocate_singlerefbg(singlebg, nlevp, nlevq)

implicit none

type(singlerefbg_type), intent(out) :: singlebg
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

end subroutine allocate_singlerefbg


!------------------------------------------------------------------------------
!> Dealloate the singlerefbg_type structure
!!
!! \author Met Office
!!
!! \date 13/05/2025: Created
!!
subroutine deallocate_singlerefbg(singlebg)

implicit none

type(singlerefbg_type), intent(inout) :: singlebg

deallocate(singlebg % za)
deallocate(singlebg % zb)
deallocate(singlebg % p)
deallocate(singlebg % q)

end subroutine deallocate_singlerefbg

! ------------------------------------------------------------------------------
end module ufo_refractivityonedvarcheck_utils_mod
