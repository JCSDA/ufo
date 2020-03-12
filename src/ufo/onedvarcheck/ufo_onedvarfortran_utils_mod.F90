! (C) copyright 2018 UCAR
!
! this software is licensed under the terms of the apache licence version 2.0
! which can be obtained at http://www.apache.org/licenses/license-2.0.

!> fortran module to provide code shared between nonlinear and tlm/adm radiance calculations

module ufo_onedvarfortran_utils_mod

use iso_c_binding
use config_mod
use kinds
use ufo_geovals_mod

implicit none

! defaults are public

integer, parameter :: max_string_length=99 !

!===============================================================================
! type definitions
!===============================================================================

!-------------
! 1. B matrix
!-------------

!The B-matrix is read from the file into this structure, but see also the
!fieldtype definitions in section 13 which help define the format of the
!contents.

type bmatrix_type
  logical :: status                    ! status indicator
  integer :: nbands                    ! number of latitude bands
  integer :: nsurf                     ! number of surface type variations
  integer :: nfields                   ! number of fields
  integer, pointer :: fields(:,:)      ! fieldtypes and no. elements in each
  real, pointer :: store(:,:,:)        ! original b-matrices read from the file
  real, pointer :: inverse(:,:,:)      ! inverse of above
  real, pointer :: sigma(:,:)          ! diagonal elements
  real, pointer :: proxy(:,:)          ! copy of original for manipulation
  real, pointer :: inv_proxy(:,:)      ! copy of inverse
  real, pointer :: sigma_proxy(:,:)    ! copy of diagonal
  real, pointer :: south(:)            ! s limit of each latitude band
  real, pointer :: north(:)            ! n limit of each latitude band
end type bmatrix_type

!----------------------
! 2. Profile Information
!----------------------

!The profile variables locate a particular field, or element thereof, in the
!1d-var profile vector. absence of a field will be designated with zero. note
!that retrieval fields also require b matrix fieldtype definitions (see
!section 13 below).

!also, remember to update ops_satrad_initprofinfo and ops_satrad_mapprofiletob
!if adding fields to this structure.

type profileinfo_type

 !general

  integer :: nelements

 !atmosphere (locate start and end points of relevant fields)

  integer :: t(2)           ! temperature profile
  integer :: q(2)           ! water vapour profile (specific humidity)
  integer :: ql(2)          ! liquid water profile
  integer :: qt(2)          ! total water profile
  integer :: qi(2)          ! frozen ice  profile
  integer :: cf(2)          ! cloud fraction profile
  integer :: o3total        ! total column ozone
  integer :: o3profile(2)   ! ozone profile
  integer :: lwp            ! liquid water path

 !surface

  integer :: t2             ! screen temperature
  integer :: q2             ! screen specific humidity
  integer :: rh2            ! screen relative humidity
  integer :: tstar          ! skin temperature
  integer :: pstar          ! surface pressure
  integer :: mwemiss(2)     ! microwave emissivity
  integer :: emisspc(2)     ! emissivity principal components
  integer :: windspeed      ! surface windspeed

 !cloud (single-level grey cloud only)
  integer :: cloudtopp      ! single-level cloud top pressure
  integer :: cloudfrac      ! effective cloud fraction

 !other

  integer :: t70hpa         ! temperature at 70hpa
  integer :: t700hpa        ! temperature at 700hpa
  integer :: t950hpa        ! temperature at 950hpa
  integer :: t1000hpa       ! temperature at 1000hpa
  integer :: qsurf          ! surface specific humidity

end type

!---------------------------------------------------------
! 2. container for information about a single observation
!---------------------------------------------------------

type obinfo_type

  character(len=max_string_length) :: forward_mod_name
  integer         :: nlocs
  real(kind_real) :: latitude
  real(kind_real) :: longitude
  real(kind_real) :: elevation
  real(kind_real) :: sensor_zenith_angle
  real(kind_real) :: sensor_azimuth_angle
  real(kind_real) :: solar_zenith_angle
  real(kind_real) :: solar_azimuth_angle
  real(kind_real),allocatable :: yobs(:)

end type

!===============================================================================
! variables definitions
!===============================================================================

!---------------------------
! 1. 1d-var profile elements
!---------------------------

! define id codes for 1d-var retrieval fields.
! a list of these fieldtype codes is always present in the header of the bmatrix
! file and it's that list which decides the form of the retrieval vector.
!
! new definitions should be made in conjunction with the profileinfo_type
! structure found in section 3 above.

integer, parameter :: nfieldtypes = 19
integer, parameter :: &
  fieldtype_t          =  1, &   ! temperature
  fieldtype_q          =  2, &   ! specific humidity profile
  fieldtype_t2         =  3, &   ! surface air temperature
  fieldtype_q2         =  4, &   ! surface spec humidity
  fieldtype_tstar      =  5, &   ! surface skin temperature
  fieldtype_pstar      =  6, &   ! surface pressure
  fieldtype_o3total    =  7, &   ! total column ozone
  fieldtype_not_used   =  8, &   ! not currently in use
  fieldtype_ql         =  9, &   ! liquid water profile
  fieldtype_qt         = 10, &   ! total water profile
  fieldtype_windspeed  = 11, &   ! surface wind speed
  fieldtype_o3profile  = 12, &   ! ozone
  fieldtype_lwp        = 13, &   ! liquid water path
  fieldtype_mwemiss    = 14, &   ! microwave emissivity
  fieldtype_qi         = 15, &   ! ice profile
  fieldtype_cloudtopp  = 16, &   ! single-level cloud top pressure
  fieldtype_cloudfrac  = 17, &   ! effective cloud fraction
  fieldtype_emisspc    = 18, &   ! emissivity prinipal components
  fieldtype_cf         = 19      ! cloud fraction profile

character(len=*), parameter :: fieldtype_text(nfieldtypes) = &
  (/ 't                   ', &
     'q water vapour      ', &
     't2                  ', &
     'q2                  ', &
     'tstar               ', &
     'pstar               ', &
     'ozone (total column)', &
     '[unused field type] ', &
     'q liquid            ', &
     'q total             ', &
     'wind speed          ', &
     'ozone (profile)     ', &
     'liquid water path   ', &
     'microwave emissivity', &
     'q ice               ', &
     'cloud top pressure  ', &
     'cloud fraction      ', &
     'emissivity pcs      ', &
     'cloud fraction prof ' /)

!------------------------------
! 2. miscellaneous definitions
!------------------------------

!data source codes, used to assign background ozone information. the numbers
!are arbitrary.

integer, parameter :: &
  source_bg        = 1, &
  source_estimate  = 2, &
  source_product   = 3, &
  source_reference = 4

! for now default ozone profile to be provided by model
integer :: ozonesource = source_bg

contains

!===============================================================================
! subroutines
! order is inportant at the moment will look into interface blocks
!===============================================================================

subroutine ufo_onedvarfortran_iogetfreeunit(unit)

!-------------------------------------------------------------------------------
! initialize profile indices to zero. the profinfo structure is used to store
! the locations of fields in the retrieval vector. zero will indicate absence of
! a field.
!-------------------------------------------------------------------------------

implicit none

integer, intent(out) :: unit

integer, parameter :: unit_min=10
integer, parameter :: unit_max=1000
logical            :: opened
integer            :: lun
integer            :: newunit

newunit=-1
do lun=unit_min,unit_max
  inquire(unit=lun,opened=opened)
  if (.not. opened) then
      newunit=lun
    exit
  end if
end do
unit=newunit

end subroutine ufo_onedvarfortran_iogetfreeunit

!===============================================================================
! functions
!===============================================================================

! this function should be moved to a variables_mod to handle variables in fortran
function find_index(vars, var)
  
  character(len=max_string_length)  :: vars(:)
  character(len=max_string_length)  :: var
  integer             :: find_index
  integer             :: jv
  character(len=250)  :: buf

  find_index = -1
  do jv = 1, size(vars)
    if (trim(vars(jv)) == trim(var)) find_index = jv
  end do

  write(buf,*)'ufo background check: unknown variable: ',trim(var)
  if (find_index < 1) call abor1_ftn(buf)

end function find_index

end module ufo_onedvarfortran_utils_mod
