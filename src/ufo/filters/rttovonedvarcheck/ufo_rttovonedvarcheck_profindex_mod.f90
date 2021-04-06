! (C) British Crown Copyright 2017-2018 Met Office
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 

!> Fortran module containing profile index

module ufo_rttovonedvarcheck_profindex_mod

use kinds
use fckit_log_module, only : fckit_log
use ufo_metoffice_bmatrixstatic_mod
use ufo_rttovonedvarcheck_constants_mod

implicit none
private

type, public :: ufo_rttovonedvarcheck_profindex

  ! general
  integer :: nprofelements !< number of profile elements being used
  integer :: nlevels       !< number of model levels

  ! atmosphere (locate start and end points of relevant fields)
  integer :: t(2)           !< temperature profile
  integer :: q(2)           !< water vapour profile (specific humidity)
  integer :: ql(2)          !< liquid water profile
  integer :: qt(2)          !< total water profile
  integer :: qi(2)          !< frozen ice  profile
  integer :: cf(2)          !< cloud fraction profile
  integer :: o3total        !< total column ozone
  integer :: o3profile(2)   !< ozone profile
  integer :: lwp            !< liquid water path

  ! surface
  integer :: t2             !< screen temperature
  integer :: q2             !< screen specific humidity
  integer :: rh2            !< screen relative humidity
  integer :: tstar          !< skin temperature
  integer :: pstar          !< surface pressure
  integer :: mwemiss(2)     !< microwave emissivity
  integer :: emisspc(2)     !< emissivity principal components
  integer :: windspeed      !< surface windspeed

  ! cloud (single-level grey cloud only)
  integer :: cloudtopp      !< single-level cloud top pressure
  integer :: cloudfrac      !< effective cloud fraction

  ! other
  integer :: t70hpa         !< temperature at 70hpa - used for ozone profile, not currently implemented
  integer :: t700hpa        !< temperature at 700hpa - used for ozone profile, not currently implemented
  integer :: t950hpa        !< temperature at 950hpa - used for ozone profile, not currently implemented
  integer :: t1000hpa       !< temperature at 1000hpa - used for ozone profile, not currently implemented
  integer :: qsurf          !< surface specific humidity

contains
  procedure :: setup  => ufo_rttovonedvarcheck_profindex_setup
  procedure :: delete => ufo_rttovonedvarcheck_profindex_delete

end type

contains

!-------------------------------------------------------------------------------
!> Profile index setup
!!
!! \details Heritage: Ops_SatRad_MapProfileToB
!!
!! Setup the profile index which requires the bmatrix object.
!!
!! \author Met Office
!!
!! \date 09/06/2020: Created
!!
subroutine ufo_rttovonedvarcheck_profindex_setup(self, bmatrix, nlevels)

implicit none

! subroutine arguments:
class(ufo_rttovonedvarcheck_profindex), intent(inout) :: self    !< profindex structure
type(ufo_metoffice_bmatrixstatic), intent(in)         :: bmatrix !< state error covariances
integer, intent(in)                                   :: nlevels !< number of model levels

! local constants:
character(len=*), parameter :: routinename = "ufo_rttovonedvarcheck_profindex_setup"

! local variables:
integer :: i,j
integer :: firstelement
integer :: lastelement
integer :: nelements

! Initialise to zeros
call self % delete()
nelements = 0
self % nlevels = nlevels

!loop through fields in bmatrix. the number of elements that the field is composed of
!is held in bmatrix % fields(i,2).
do j = 1, bmatrix % nfields
  firstelement = nelements + 1
  lastelement  = nelements + bmatrix % fields(j,2)

  !assign start and end points. if the field wasn't found then assign a value of
  !zero, which indicates absence.

  select case( bmatrix % fields(j,1) )

   !----------
   !atmosphere (set start and end points for multi-level fields)
   !----------

    case( ufo_metoffice_fieldtype_t )
      self % t(1)         = firstelement
      self % t(2)         = lastelement

    case( ufo_metoffice_fieldtype_q )
      self % q(1)         = firstelement
      self % q(2)         = lastelement

    case( ufo_metoffice_fieldtype_ql )
      call abor1_ftn("ufo_metoffice_fieldtype_ql: Not currently implemented aborting")

    case( ufo_metoffice_fieldtype_qi )
      call abor1_ftn("ufo_metoffice_fieldtype_qi: Not currently implemented aborting")

    case( ufo_metoffice_fieldtype_cf )
      call abor1_ftn("ufo_metoffice_fieldtype_cf: Not currently implemented aborting")

    case( ufo_metoffice_fieldtype_qt )
      self % qt(1)        = firstelement
      self % qt(2)        = lastelement

    case( ufo_metoffice_fieldtype_o3profile )
      call abor1_ftn("ufo_metoffice_fieldtype_o3profile: Not currently implemented aborting")

    case( ufo_metoffice_fieldtype_o3total )
      call abor1_ftn("ufo_metoffice_fieldtype_o3total: Not currently implemented aborting")

   !-------
   !surface
   !-------

    case( ufo_metoffice_fieldtype_t2 )
      self % t2         = firstelement

    case( ufo_metoffice_fieldtype_q2 )
      self % q2         = firstelement

    case( ufo_metoffice_fieldtype_tstar )
      self % tstar      = firstelement

    case( ufo_metoffice_fieldtype_pstar )
      self % pstar      = firstelement

    case( ufo_metoffice_fieldtype_windspeed )
      call abor1_ftn("ufo_metoffice_fieldtype_windspeed: Not currently implemented aborting")

    case( ufo_metoffice_fieldtype_mwemiss )
      call abor1_ftn("ufo_metoffice_fieldtype_mwemiss: Not currently implemented aborting")

    case( ufo_metoffice_fieldtype_emisspc )
      call abor1_ftn("ufo_metoffice_fieldtype_emisspc: Not currently implemented aborting")

   !------------------------------------
   !cloud (single-level grey cloud only)
   !------------------------------------

    case( ufo_metoffice_fieldtype_cloudtopp )
      call abor1_ftn("ufo_metoffice_fieldtype_cloudtopp: Not currently implemented aborting")

    case( ufo_metoffice_fieldtype_cloudfrac )
      call abor1_ftn("ufo_metoffice_fieldtype_cloudfrac: Not currently implemented aborting")

    case( ufo_metoffice_fieldtype_not_used ) ! currently unused
      continue

    case default
      write(*,*) 'invalid field type in b matrix file: ',j
      cycle

  end select

  if ( firstelement /= 0 ) nelements = nelements + bmatrix % fields(j,2)

end do

self % nprofelements = nelements

end subroutine ufo_rttovonedvarcheck_profindex_setup

!-------------------------------------------------------------------------------
!> Delete profile index
!!
!! \details Heritage: Ops_SatRad_InitProfInfo.f90
!!
!! Reset profile index
!!
!! \author Met Office
!!
!! \date 09/06/2020: Created
!!
subroutine ufo_rttovonedvarcheck_profindex_delete(self)

implicit none
class(ufo_rttovonedvarcheck_profindex), intent(inout) :: self !< profile index structure

! Zero all values
self % nprofelements = 0
self % t = 0
self % q = 0
self % ql = 0
self % qt = 0
self % qi = 0
self % cf = 0
self % o3total = 0
self % o3profile = 0
self % t2 = 0
self % q2 = 0
self % rh2 = 0
self % tstar = 0
self % pstar = 0
self % windspeed = 0
self % t70hpa = 0
self % t700hpa = 0
self % t950hpa = 0
self % t1000hpa = 0
self % qsurf = 0
self % lwp = 0
self % mwemiss = 0
self % cloudtopp = 0
self % cloudfrac = 0
self % emisspc = 0

end subroutine ufo_rttovonedvarcheck_profindex_delete

! ------------------------------------------------------------------------------

end module ufo_rttovonedvarcheck_profindex_mod
