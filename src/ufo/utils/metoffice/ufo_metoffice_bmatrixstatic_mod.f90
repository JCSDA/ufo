! (C) British Crown Copyright 2020 Met Office
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

!> Fortran module containing the full b-matrix data type and methods for the 1D-Var.

module ufo_metoffice_bmatrixstatic_mod

use fckit_log_module, only : fckit_log
use kinds
use ufo_constants_mod, only: zero, one
use ufo_utils_mod, only : ufo_utils_iogetfreeunit, InvertMatrix
use ufo_vars_mod

implicit none
private

type, public :: ufo_metoffice_bmatrixstatic
  logical :: status                                 !< status indicator
  integer :: nbands                                 !< number of latitude bands
  integer :: nsurf                                  !< number of surface type variations
  integer :: nfields                                !< number of fields
  integer, allocatable :: fields(:,:)               !< fieldtypes and no. elements in each
  real(kind=kind_real), allocatable :: store(:,:,:)     !< original b-matrices read from the file
  real(kind=kind_real), allocatable :: inverse(:,:,:)   !< inverse of above
  real(kind=kind_real), allocatable :: sigma(:,:)       !< diagonal elements
  real(kind=kind_real), allocatable :: proxy(:,:)       !< copy of original for manipulation
  real(kind=kind_real), allocatable :: inv_proxy(:,:)   !< copy of inverse
  real(kind=kind_real), allocatable :: sigma_proxy(:,:) !< copy of diagonal
  real(kind=kind_real), allocatable :: south(:)         !< s limit of each latitude band
  real(kind=kind_real), allocatable :: north(:)         !< n limit of each latitude band
contains
  procedure :: setup  => ufo_metoffice_bmatrixstatic_setup
  procedure :: delete => ufo_metoffice_bmatrixstatic_delete
  procedure :: reset  => ufo_metoffice_bmatrixstatic_reset
end type ufo_metoffice_bmatrixstatic

character(len=200)     :: message

!-----------------------------------------------------------------------------
! 1. 1d-var profile elements
!-----------------------------------------------------------------------------

! define id codes for 1d-var retrieval fields.
! a list of these fieldtype codes is always present in the header of the bmatrix
! file and it's that list which decides the form of the retrieval vector.
!
! new definitions should be made in conjunction with the profileinfo_type
! structure found in ufo_rttovonedvarcheck_profindex_mod.F90.

integer, parameter, public :: nfieldtypes_ukmo = 19 !< number of fieldtypes
integer, parameter, public :: &
  ufo_metoffice_fieldtype_t          =  1, &   !< temperature
  ufo_metoffice_fieldtype_q          =  2, &   !< specific humidity profile
  ufo_metoffice_fieldtype_t2         =  3, &   !< surface air temperature
  ufo_metoffice_fieldtype_q2         =  4, &   !< surface spec humidity
  ufo_metoffice_fieldtype_tstar      =  5, &   !< surface skin temperature
  ufo_metoffice_fieldtype_pstar      =  6, &   !< surface pressure
  ufo_metoffice_fieldtype_o3total    =  7, &   !< total column ozone - not currently setup
  ufo_metoffice_fieldtype_not_used   =  8, &   !< not currently in use - not currently setup
  ufo_metoffice_fieldtype_ql         =  9, &   !< liquid water profile - not currently setup
  ufo_metoffice_fieldtype_qt         = 10, &   !< total water profile
  ufo_metoffice_fieldtype_windspeed  = 11, &   !< surface wind speed
  ufo_metoffice_fieldtype_o3profile  = 12, &   !< ozone - not currently setup
  ufo_metoffice_fieldtype_lwp        = 13, &   !< liquid water path - not currently setup
  ufo_metoffice_fieldtype_mwemiss    = 14, &   !< microwave emissivity - not currently setup
  ufo_metoffice_fieldtype_qi         = 15, &   !< ice profile - not currently setup
  ufo_metoffice_fieldtype_cloudtopp  = 16, &   !< single-level cloud top pressure
  ufo_metoffice_fieldtype_cloudfrac  = 17, &   !< effective cloud fraction
  ufo_metoffice_fieldtype_emisspc    = 18, &   !< emissivity prinipal components - not currently setup
  ufo_metoffice_fieldtype_cf         = 19      !< cloud fraction profile - not currently setup

character(len=*), parameter, public :: ufo_metoffice_fieldtype_text(nfieldtypes_ukmo) = &
  [character(len=MAXVARLEN) :: var_ts,                 &
     var_q,                  &
     var_sfc_t2m,            &
     var_sfc_q2m,            &
     var_sfc_tskin,          &
     var_ps,                 &
     'ozone (total column)', &
     '[unused field type] ', &
     var_clw,                &
     'q total             ', &
     var_sfc_wspeed,         &
     'ozone (profile)     ', &
     'liquid water path   ', &
     var_sfc_emiss,          &
     var_cli,                &
     'cloud top pressure  ', &
     'cloud fraction      ', &
     'emissivity pcs      ', &
     var_cldfrac_vol ]

contains

! ------------------------------------------------------------------------------------------------
!> Routine to read and setup the 1D-Var B-matrix
!!
!! \author Met Office
!!
!! \date 09/06/2020: Created
!!
subroutine ufo_metoffice_bmatrixstatic_setup(self, variables, filepath, qtotal_flag)

implicit none
class(ufo_metoffice_bmatrixstatic), intent(inout) :: self !< B-matrix Covariance
character(len=*), intent(in)       :: variables(:) !< Model variables in B matrix
character(len=*), intent(in)       :: filepath     !< Path to B matrix file
logical, intent(in)                :: qtotal_flag  !< Flag for qtotal

logical                       :: file_exists  ! Check if a file exists logical
integer                       :: fileunit     ! Unit number for reading in files
integer, allocatable          :: fields_in(:) ! Fields_in used to subset b-matrix for testing.
integer, allocatable          :: nonzero_fields(:) ! fields_in with nonzero entries removed
integer, allocatable          :: nonzero_fields_in(:) ! fields_in with nonzero entries removed
real(kind=kind_real)          :: t1,t2        ! Time values for logging
character(len=:), allocatable :: str
logical                       :: testing = .false.
integer                       :: ii, jj
logical                       :: match

call fckit_log % info("ufo_metoffice_bmatrixstatic_setup start")

call cpu_time(t1)

! Open file and read in b-matrix
inquire(file=trim(filepath), exist=file_exists)
if (file_exists) then
  fileunit = ufo_utils_iogetfreeunit()
  open(unit = fileunit, file = trim(filepath))
  call self % delete()
  call rttovonedvarcheck_create_fields_in(fields_in, variables, qtotal_flag)
  if (testing) then
    call rttovonedvarcheck_covariance_GetBmatrix(self, fileunit, fieldlist=fields_in)
  else
    call rttovonedvarcheck_covariance_GetBmatrix(self, fileunit)
  end if
  close(unit = fileunit)
  call fckit_log % info("rttovonedvarcheck bmatrix file exists and read in")
else
  call abor1_ftn("rttovonedvarcheck bmatrix file not found")
end if

call cpu_time(t2)

! Check the yaml input contains all required b-matrix elements
do ii = 1, size(self % fields(:,1)) ! loop over b-matrix elements
  if (self % fields(ii,1) == 0) cycle
  match = .false.
  do jj = 1, size(fields_in) ! loop over array generated from yaml
    if (self % fields(ii,1) == fields_in(jj)) match = .true.
  end do
  if (.not. match) then
    write(*,*) "input model variables do not have ",ufo_metoffice_fieldtype_text(self % fields(ii,1))
    call abor1_ftn("rttovonedvarcheck not all the model data is available for the b-matrix")
  end if
end do

! Check the yaml input is in the same order as b-matrix elements
! (this is to ensure the b-matrix rows & columns are in the anticipated order)
nonzero_fields_in = pack(fields_in, fields_in /= 0)
nonzero_fields = pack(self % fields(:,1), self % fields(:,1) /= 0)
if (.not. all(nonzero_fields(1:min(size(nonzero_fields_in),size(nonzero_fields))) == &
              nonzero_fields_in(1:min(size(nonzero_fields_in),size(nonzero_fields))))) then
  write(*,*) "Supplied field list [L] is not in the order contained in the B-matrix [R]:"
  do ii = 1, size(self % fields(:,1))
    if (self % fields(ii,1) == 0) exit
    write(*,*) ufo_metoffice_fieldtype_text(nonzero_fields_in(ii)), &
               ufo_metoffice_fieldtype_text(self % fields(ii,1))
  end do
  call abor1_ftn("Supplied field list is expected in the B-matrix order")
endif

write(message,*) "ufo_metoffice_bmatrixstatic_setup cpu time = ",(t2-t1)
call fckit_log % info(message)

end subroutine ufo_metoffice_bmatrixstatic_setup

! ------------------------------------------------------------------------------------------------
!> Routine to delete and squash the 1D-Var B-matrix
!!
!! \details Met Office OPS Heritage: Ops_SatRad_SquashBmatrix.f90
!!
!! \author Met Office
!!
!! \date 09/06/2020: Created
!!
subroutine ufo_metoffice_bmatrixstatic_delete(self)

implicit none
class(ufo_metoffice_bmatrixstatic), intent(inout) :: self  !< B-matrix Covariance

character(len=*), parameter :: RoutineName = "ufo_metoffice_bmatrixstatic_delete"

self % status = .false.
self % nbands = 0
self % nsurf = 0
if ( allocated(self % fields)      ) deallocate( self % fields      )
if ( allocated(self % store)       ) deallocate( self % store       )
if ( allocated(self % inverse)     ) deallocate( self % inverse     )
if ( allocated(self % sigma)       ) deallocate( self % sigma       )
if ( allocated(self % proxy)       ) deallocate( self % proxy       )
if ( allocated(self % inv_proxy)   ) deallocate( self % inv_proxy   )
if ( allocated(self % sigma_proxy) ) deallocate( self % sigma_proxy )
if ( allocated(self % south)       ) deallocate( self % south       )
if ( allocated(self % north)       ) deallocate( self % north       )

end subroutine ufo_metoffice_bmatrixstatic_delete

! ------------------------------------------------------------------------------------------------
!> Routine to initialize the 1D-Var B-matrix
!!
!! \details Met Office OPS Heritage: Ops_SatRad_GetBmatrix.f90
!!
!! read the input file and allocate and fill in all the components of the bmatrix
!! structure, with the exception of proxy variables.
!!
!! notes:
!!
!! it is assumed a valid bmatrix file is available and has been opened ready for
!! reading, accessed via the input unit number.
!!
!! two optional arguments are provided to allow a submatrix to be formed from the
!! original file, either depending on a list of fields to be used for retrieval,
!! or a list of specific element numbers to be retained.
!!
!! the file header may begin with any number of comment lines, which are defined
!! as those using either # or ! as the first non-blank character.
!!
!! immediately following any comments, the header should then contain information
!! on the size of the matrix to be read in, and a description of each matrix
!! element. memory will be allocated depending on these matrix specifications.
!!
!! each matrix should then have a one line header containing the latitude bounds
!! and a latitude band id. band id numbers should be sequential in latitude. i.e.
!! band 1 will begin at -90.0, band 2 from -60.0 ... (or however wide we choose
!! the bands).
!!
!! we also calculate an inverse of each b matrix, used in the current 1d-var for
!! cost function monitoring only but it may be required for other minimization
!! methods at some point.
!!
!! standard deviations are also stored in a separate vector.
!!
!! proxy variables are provided to allow manipulation of the chosen matrix during
!! processing without affecting the original. we might wish to do this, for
!! example, to take account of the different surface temperature errors for land
!! and sea. this routine makes no assumption about the use of these variables,
!! hence no space is allocated.
!!
!! \author Met Office
!!
!! \date 09/06/2020: Created
!!
subroutine rttovonedvarcheck_covariance_GetBmatrix (self,           &
                                                    fileunit,       &
                                                    b_elementsused, &
                                                    fieldlist)

implicit none

! subroutine arguments:
type (ufo_metoffice_bmatrixstatic), intent(inout) :: self !< B-matrix covariance
integer, intent(in)                :: fileunit           !< free file unit number
integer, optional, intent(in)      :: b_elementsused(:)  !< optional: list of elements used
integer, optional, intent(inout)   :: fieldlist(:)       !< optional: list of fields used

! local declarations:
character(len=*), parameter        :: routinename = "rttovonedvarcheck_covariance_GetBmatrix"
integer                            :: i
integer                            :: j
integer                            :: k
integer                            :: readstatus
integer                            :: status
integer                            :: nbands
integer                            :: matrixsize
integer                            :: band
integer                            :: nmatrix
integer                            :: nbfields
integer                            :: nelements
integer                            :: nelements_total
real(kind=kind_real)               :: southlimit
real(kind=kind_real)               :: northlimit
character(len=80)                  :: line
character(len=3)                   :: fieldtype
real(kind=kind_real), allocatable  :: bfromfile(:,:)
logical, allocatable               :: list(:)
integer, allocatable               :: bfields(:,:)
integer, allocatable               :: elementsused(:)

!--------------------------
! 1. read header information
!--------------------------

!----
! 1.1) header comments
!----

! read individual lines until one is found where the first non-blank character
! is not # or !

do
  read (fileunit, '(a)') line
  line = adjustl (line)
  if (verify (line, '#!') == 1) then
    backspace (fileunit)
    exit
  end if
end do

!----
! 1.2) matrix information
!----

nmatrix = 0
nbfields = 0
read (fileunit, *) matrixsize, nbands, nbfields  ! matrix size, no. of latitude
if (nbfields > 0) then                           ! bands, no. of fields
  allocate (bfields(nbfields,2))
  do i = 1, nbfields
    read (fileunit, *) bfields(i,1), bfields(i,2) ! field id, number of elements
  end do
end if

!----
! 1.3) output initial messages
!----

call fckit_log % debug('reading b matrix file:')
write (message, '(a,i0)') 'number of latitude bands = ', nbands
call fckit_log % debug(message)
write (message, '(a,i0)') 'matrix size = ', matrixsize
call fckit_log % debug(message)
if (nbfields > 0) then
  call fckit_log % debug('order of fields and number of elements in each:')
  do i = 1, nbfields
    write (message, '(i0,a)') bfields(i,2), ' x ' // ufo_metoffice_fieldtype_text(bfields(i,1))
    call fckit_log % debug(message)
  end do
end if

!-------------------------------------------------------
! 2. read matrix field mappings and generate element list
!-------------------------------------------------------

! each field in the bmatrix has a designated id number (defined in
! opsmod_satrad_info) that is included in the file header. under certain
! circumstances there may be a requirement to exclude some of the fields from
! the matrix and this can be achieved by including the fieldlist argument which
! should contain the required id values.
! this section checks the fields in fieldlist against those in the file and makes
! up a final list of those common to both. subsequently, this is translated to
! the exact element numbers that we wish to keep. note that an alternative to
! this method is to provide the element numbers explicitly in the argument
! b_elementsused.

allocate (self % fields(nfieldtypes_ukmo,2))
self % fields(:,:) = 0
self % status = .true.
nelements = 0
nelements_total = 0
self % nfields = 0

if (present (fieldlist)) then

  !----
  ! 2.1) use retrieval field list if optional argument fieldlist present
  !----

  allocate (elementsused(matrixsize))
  blist: do j = 1, nbfields
    do i = 1, size (fieldlist)
      if (fieldlist(i) == 0) cycle
      if (bfields(j,1) == fieldlist(i)) then
        do k = 1, bfields(j,2)
          elementsused(nelements + k) = nelements_total + k
        end do
        nelements_total = nelements_total + bfields(j,2)
        nelements = nelements + bfields(j,2)
        fieldlist(i) = 0  ! store only fields that are not used in here
        cycle blist
      end if
    end do
    bfields(j,1) = 0  ! store only fields that are used in here
    nelements_total = nelements_total + bfields(j,2)
  end do blist

  nbfields = count (bfields(:,1) /= 0)
  bfields(1:nbfields,2) = pack (bfields(:,2), bfields(:,1) /= 0)
  bfields(1:nbfields,1) = pack (bfields(:,1), bfields(:,1) /= 0)

  ! 2.1.1) write messages

  if (any (fieldlist /= 0)) then
    call fckit_log % debug('the following requested retrieval fields are not in the b matrix:')
    do i = 1, size (fieldlist)
      if (fieldlist(i) /= 0) then
        write (fieldtype, '(i0)') fieldlist(i)
        if (fieldlist(i) > 0 .and. fieldlist(i) <= nfieldtypes_ukmo) then
          write (message, '(a)') trim (ufo_metoffice_fieldtype_text(fieldlist(i))) // &
            ' (fieldtype ' // trim (adjustl (fieldtype)) // ')'
          call fckit_log % debug(message)
        else
          write (message, '(a)') 'fieldtype ' // trim (adjustl (fieldtype)) // &
            ' which is invalid'
          call fckit_log % debug(message)
        end if
      end if
    end do
  end if

  call fckit_log % debug('b matrix fields used to define the retrieval profile vector:')
  do i = 1, nbfields
    write (message, '(a)') ufo_metoffice_fieldtype_text(bfields(i,1))
    call fckit_log % debug(message)
  end do

  ! 2.1.2) reset fieldlist so that only used fields are now stored

  fieldlist(1:nbfields) = bfields(1:nbfields,1)
  if (nbfields < size (fieldlist)) fieldlist(nbfields + 1:) = 0

else if (present (b_elementsused)) then

  !----
  ! 2.2) or use exact element numbers if optional argument b_elementsused present
  !----

  allocate (elementsused(size (b_elementsused)))
  if (any (b_elementsused < 0 .and. b_elementsused > matrixsize)) then
    write(*,*) routinename // ' : invalid b matrix elements present in input list'
  end if
  nelements = count (b_elementsused > 0 .and. b_elementsused <= matrixsize)
  elementsused(1:nelements) = pack (b_elementsused, b_elementsused > 0 .and. b_elementsused <= matrixsize)

else

  !----
  ! 2.3) or use everything by default
  !----

  nelements = matrixsize

end if

!---------------------------------------------------------
! 3. read the file and store in bmatrix structure variables
!---------------------------------------------------------

! set field list in structure variable

if (nbfields > 0) self % fields(1:nbfields,:) = bfields(1:nbfields,:)
self % nfields = nbfields

! initialize the rest

self % nbands = nbands
allocate (self % store(nelements,nelements,nbands))
allocate (self % south(nbands))
allocate (self % north(nbands))
self % store(:,:,:) = zero

! allocate dummy variables for file input

allocate (list(nbands))
allocate (bfromfile(matrixsize,matrixsize))
list(:) = .false.

readallb : do

  !----
  ! 3.1) read matrix
  !----

  ! latitude band information

  read (fileunit, '(i3,2f8.2)', iostat = readstatus) band, southlimit, northlimit
  if (readstatus < 0) exit

  ! matrix data

  do j = 1, matrixsize
    read (fileunit, '(5e16.8)' ) (bfromfile(i,j), i = 1, matrixsize)
  end do

  write (message, '(a,i0,a,2f8.2)') &
    'band no.', band, ' has southern and northern latitude limits of', &
    southlimit, northlimit
  call fckit_log % debug(message)

  !----
  ! 3.2) transfer into 1d-var arrays
  !----

  if (band > 0 .and. band <= nbands) then
    list(band) = .true.
    nmatrix = nmatrix + 1
    if (nelements < matrixsize) then
      self % store(:,:,band) = &
        bfromfile(elementsused(1:nelements),elementsused(1:nelements))
    else
      self % store(:,:,band) = bfromfile(:,:)
    end if
    self % south(band) = southlimit
    self % north(band) = northlimit
  else
    write (message, '(a,i0)') 'skipped matrix with band number ', band
    call fckit_log % debug(message)
  end if

end do readallb

!---------------------
! 4. store error vector
!---------------------

allocate (self % sigma(nelements,nbands))
do i = 1, nelements
  self % sigma(i,:) = sqrt (self % store(i,i,:))
end do

!----------------
! 5. invert matrix
!----------------

allocate (self % inverse(nelements,nelements,nbands))
self % inverse(:,:,:) = self % store(:,:,:)

do k = 1, nbands
  call InvertMatrix (nelements,              & ! in
                     nelements,              & ! in
                     self % inverse(:,:,k),  & ! inout
                     status)                   ! out
  if (status /= 0) then
    call abor1_ftn("rttovonedvarcheck: bmatrix is not invertible")
  end if
end do

!----------
! 6. tidy up
!----------

if (nmatrix < nbands) then
  call abor1_ftn("rttovonedvarcheck: too few b matrices found in input file")
end if

end subroutine rttovonedvarcheck_covariance_GetBmatrix

! ------------------------------------------------------------------------------------------------
!> Routine to create error covariances for a single observation
!!
!! \details Met Office OPS Heritage: Ops_SatRad_ResetCovariances.f90
!!
!! \author Met Office
!!
!! \date 08/09/2020: Created
!!
subroutine ufo_metoffice_bmatrixstatic_reset(self, latitude, & ! in
                                      b_matrix, b_inverse, b_sigma ) ! out

implicit none

! Subroutine arguments
class(ufo_metoffice_bmatrixstatic), intent(in) :: self !< B-matrix covariance
real(kind_real), intent(in)  :: latitude
real(kind_real), intent(out) :: b_matrix(:,:)
real(kind_real), intent(out) :: b_inverse(:,:)
real(kind_real), intent(out) :: b_sigma(:)

! Local Variables
integer :: band, i

! select appropriate b matrix for latitude of observation
b_matrix(:,:) = zero
b_inverse(:,:) = zero
b_sigma(:) = zero
do band = 1, self % nbands
  if (latitude < self % north(band)) exit
end do
b_matrix(:,:) = self % store(:,:,band)
b_inverse(:,:) = self % inverse(:,:,band)
b_sigma(:) = self % sigma(:,band)

end subroutine ufo_metoffice_bmatrixstatic_reset

! ------------------------------------------------------------------------------------------------
!> Create a subset of the b-matrix.  Used for testing.
!!
!! \author Met Office
!!
!! \date 09/06/2020: Created
!!
subroutine rttovonedvarcheck_create_fields_in(fields_in, variables, qtotal_flag)

implicit none

integer, allocatable, intent(inout) :: fields_in(:) !< Array to specify fields used
character(len=*), intent(in)        :: variables(:) !< Model variables in B matrix
logical, intent(in)                 :: qtotal_flag  !< Flag for qtotal

character(len=MAXVARLEN)  :: varname
integer                   :: jvar
integer                   :: nmvars, counter
character(len=200)        :: message
logical                   :: clw_present = .false.
logical                   :: ciw_present = .false.

call fckit_log % info("rttovonedvarcheck_create_fields_in: starting")

! Which fields are being used from b matrix file - temporary until rttov can handle all fields
nmvars = size(variables)
allocate(fields_in(nmvars))
fields_in(:) = 0

counter = 0
do jvar = 1, nmvars

  counter = counter + 1
  varname = variables(jvar)

  select case (trim(varname))

    case (var_ts)
      fields_in(counter) = ufo_metoffice_fieldtype_t ! air_temperature

    case (var_q)
      if (qtotal_flag) then
        fields_in(counter) = ufo_metoffice_fieldtype_qt ! total water profile
      else
        fields_in(counter) = ufo_metoffice_fieldtype_q ! water profile
      end if

    case(var_sfc_t2m)
      fields_in(counter) = ufo_metoffice_fieldtype_t2 ! 2m air_temperature

    case(var_sfc_q2m)
      fields_in(counter) = ufo_metoffice_fieldtype_q2 ! 2m specific_humidity

    case(var_sfc_tskin)
      fields_in(counter) = ufo_metoffice_fieldtype_tstar ! surface skin temperature

    case(var_ps)
      fields_in(counter) = ufo_metoffice_fieldtype_pstar ! surface air pressure

    ! 7 - o3total is not implmented yet
    ! 8 - not used is not implmented yet

    case (var_clw)
      clw_present = .true.
      if (.NOT. qtotal_flag) then
        call abor1_ftn("rttovonedvarcheck not setup for independent clw yet")
      end if

    case (var_sfc_u10, var_sfc_v10)
      fields_in(counter) = ufo_metoffice_fieldtype_windspeed ! surface wind speed

    ! 12 - o3profile is not implmented yet
    ! 13 - lwp (liquid water path) is not implmented yet

    case (var_sfc_emiss) ! microwave emissivity
      fields_in(counter) = ufo_metoffice_fieldtype_mwemiss

    case (var_cli)
      ciw_present = .true.
      if (.NOT. qtotal_flag) then
        call abor1_ftn("rttovonedvarcheck not setup for independent ciw yet")
      end if

    case ("cloud_top_pressure")
      fields_in(counter) = ufo_metoffice_fieldtype_cloudtopp

    case ("cloud_fraction") ! cloud fraction
      fields_in(counter) = ufo_metoffice_fieldtype_cloudfrac

    case ("emissivity_pc") ! emissivity prinipal components
      fields_in(counter) = ufo_metoffice_fieldtype_emisspc

    ! 19 cloud fraction profile - not currently used

    case default
      write(message,*) trim(varname)," not implemented yet in rttovonedvarcheck Covariance"
      call abor1_ftn(message)

  end select

end do

! Check clw and ciw are both in the list of required variables
if (qtotal_flag .and. ((.not. clw_present) .or. (.not. ciw_present))) then
  call abor1_ftn("rttovonedvarcheck using qtotal but clw or ciw not in variables")
end if

end subroutine rttovonedvarcheck_create_fields_in

! ------------------------------------------------------------------------------------------------

end module ufo_metoffice_bmatrixstatic_mod
