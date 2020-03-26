! (C) Copyright 2018 UCAR
!
! This software is licensed under the terms of the apache licence version 2.0
! which can be obtained at http://www.apache.org/licenses/license-2.0.

!> Fortran module to provide code shared between nonlinear and tlm/adm radiance calculations

module ufo_rttovonedvarcheck_init_mod

use iso_c_binding
use config_mod
use kinds
use ufo_geovals_mod
use ufo_rttovonedvarcheck_utils_mod

implicit none

private

! subroutines
public ufo_rttovonedvarcheck_setup
public ufo_rttovonedvarcheck_InitBmatrix
public ufo_rttovonedvarcheck_InitProfInfo
public ufo_rttovonedvarcheck_InitObInfo
public ufo_rttovonedvarcheck_GetBmatrix
public ufo_rttovonedvarcheck_MapProfileToB
public ufo_rttovonedvarcheck_InvertMatrix

contains

!===============================================================================
! Subroutines
!===============================================================================

subroutine ufo_rttovonedvarcheck_setup(self, channels)

implicit none

! subroutine arguments
type(ufo_rttovonedvarcheck), intent(inout) :: self
integer(c_int), intent(in)                 :: channels(:)

! local variables
character(len=max_string_length) :: tmp

! Setup core paths and names
self % qcname = "rttovonedvarcheck"
self % b_matrix_path = config_get_string(self % conf, max_string_length, "BMatrix")
self % forward_mod_name = config_get_string(self % conf, max_string_length, "ModName")
self % nlevels = config_get_int(self % conf, "nlevels")

! Variables for profile (x,xb)
self % nmvars = size(config_get_string_vector(self % conf, max_string_length, "model_variables"))
allocate(self % model_variables(self % nmvars))
self % model_variables = config_get_string_vector(self % conf, max_string_length, "model_variables")

! Satellite channels
self % nchans = size(channels)
allocate(self % channels(self % nchans))
self % channels(:) = channels(:)
write(*,*) "nchans setup = ",self%nchans
write(*,*) "channels setup = ",self%channels

! Set defaults for 1D-var
self % qtotal = .false.
self % RTTOV_mwscattSwitch = .false.
self % use_totalice = .false.
self % UseMLMinimization = .false.
self % UseJforConvergence = .true.
self % Max1DVarIterations = 7
self % JConvergenceOption = 1
self % IterNumForLWPCheck = 2
self % ConvergenceFactor = 0.40
self % Cost_ConvergenceFactor = 0.01
self % Mqstart = 0.001 ! Marquardt starting parameter
self % Mqstep = 5.0    ! Marquardt step parameter

! Flag for total humidity
if (config_element_exists(self % conf, "qtotal")) then
  tmp = config_get_string(self % conf, max_string_length, "qtotal")
  if (trim(tmp) == 'true') self % qtotal = .true.
end if

! Flag for RTTOV MW scatt
if (config_element_exists(self % conf, "RTTOV_mwscattSwitch")) then
  tmp = config_get_string(self % conf, max_string_length, "RTTOV_mwscattSwitch")
  if (trim(tmp) == 'true') self % RTTOV_mwscattSwitch = .true.
end if

! Flag for use of total ice in RTTOV MW scatt
if (config_element_exists(self % conf, "use_totalice")) then
  tmp = config_get_string(self % conf, max_string_length, "use_totalice")
  if (trim(tmp) == 'true') self % use_totalice = .true.
end if

! Flag to turn on marquardt-levenberg minimiser
if (config_element_exists(self % conf, "UseMLMinimization")) then
  tmp = config_get_string(self % conf, max_string_length, "UseMLMinimization")
  if (trim(tmp) == 'true') self % UseMLMinimization = .true.
end if

! Flag to Use J for convergence
if (config_element_exists(self % conf, "UseJforConvergence")) then
  tmp = config_get_string(self % conf, max_string_length, "UseJforConvergence")
  if (trim(tmp) == 'true') self % UseJforConvergence = .true.
end if

! Flag to specify if delta_x has to be negative for converg. to be true
if (config_element_exists(self % conf, "Max1DVarIterations")) then
  self % Max1DVarIterations = config_get_int(self % conf, "Max1DVarIterations")
end if

! Flag to specify if delta_x has to be negative for converg. to be true
if (config_element_exists(self % conf, "JConvergenceOption")) then
  self % UseJforConvergence = config_get_int(self % conf, "UseJforConvergence")
end if

! Choose which iteration to start checking LWP
if (config_element_exists(self % conf, "IterNumForLWPCheck")) then
  self % IterNumForLWPCheck = config_get_int(self % conf, "IterNumForLWPCheck")
end if

! Flag to specify if delta_x has to be negative for converg. to be true
if (config_element_exists(self % conf, "ConvergenceFactor")) then
  self % ConvergenceFactor = config_get_real(self % conf, "ConvergenceFactor")
end if

! Flag to specify if delta_x has to be negative for converg. to be true
if (config_element_exists(self % conf, "Cost_ConvergenceFactor")) then
  self % ConvergenceFactor = config_get_real(self % conf, "Cost_ConvergenceFactor")
end if

! Flag to specify if delta_x has to be negative for converg. to be true
if (config_element_exists(self % conf, "Mqstart")) then
  self % Mqstart = config_get_real(self % conf, "Mqstart")
end if

! Flag to specify if delta_x has to be negative for converg. to be true
if (config_element_exists(self % conf, "Mqstep")) then
  self % Mqstep = config_get_real(self % conf, "Mqstep")
end if

end subroutine

!-------------------------------------------------------------------------------

subroutine ufo_rttovonedvarcheck_InitBmatrix(bmatrix)

!-------------------------------------------------------------------------------
! nullify B matrix pointers.
!-------------------------------------------------------------------------------

implicit none

! subroutine arguments:
type(bmatrix_type), intent(out) :: bmatrix

character(len=*), parameter :: routinename = "ufo_rttovonedvarcheck_InitBmatrix"

!-------------------------------------------------------------------------------

bmatrix % status = .false.
bmatrix % nbands = 0
bmatrix % nsurf = 0

nullify( bmatrix % fields      )
nullify( bmatrix % store       )
nullify( bmatrix % inverse     )
nullify( bmatrix % sigma       )
nullify( bmatrix % proxy       )
nullify( bmatrix % inv_proxy   )
nullify( bmatrix % sigma_proxy )
nullify( bmatrix % south       )
nullify( bmatrix % north       )

end subroutine ufo_rttovonedvarcheck_InitBmatrix

subroutine ufo_rttovonedvarcheck_InitProfInfo( profinfo ) ! out

!-------------------------------------------------------------------------------
! initialize profile indices to zero. the profinfo structure is used to store
! the locations of fields in the retrieval vector. zero will indicate absence of
! a field.
!-------------------------------------------------------------------------------

implicit none

! subroutine arguments:
type(profileinfo_type), intent(out) :: profinfo

character(len=*), parameter :: routinename = "ufo_rttovonedvarcheck_InitProfInfo"

!-------------------------------------------------------------------------------

profinfo % t = 0
profinfo % q = 0
profinfo % ql = 0
profinfo % qt = 0
profinfo % qi = 0
profinfo % cf = 0
profinfo % o3total = 0
profinfo % o3profile = 0
profinfo % t2 = 0
profinfo % q2 = 0
profinfo % rh2 = 0
profinfo % tstar = 0
profinfo % pstar = 0
profinfo % windspeed = 0
profinfo % t70hpa = 0
profinfo % t700hpa = 0
profinfo % t950hpa = 0
profinfo % t1000hpa = 0
profinfo % qsurf = 0
profinfo % lwp = 0
profinfo % mwemiss = 0
profinfo % cloudtopp = 0
profinfo % cloudfrac = 0
profinfo % emisspc = 0

end subroutine ufo_rttovonedvarcheck_InitProfInfo

!-------------------------------------------------------------------------------

subroutine ufo_rttovonedvarcheck_InitObInfo(ob_info, & ! out
                                       nchans)    ! in

implicit none

! subroutine arguments:
type(obinfo_type), intent(out) :: ob_info
integer, intent(in)  :: nchans

character(len=*), parameter :: routinename = "ufo_rttovonedvarcheck_InitObInfo"

ob_info % nlocs = 1
ob_info % latitude = 0.0
ob_info % longitude = 0.0
ob_info % elevation = 0.0
ob_info % sensor_zenith_angle = 0.0
ob_info % sensor_azimuth_angle = 0.0
ob_info % solar_zenith_angle = 0.0
ob_info % solar_azimuth_angle = 0.0

! Moved this to obs loop to allow channel selection
!allocate(ob_info%yobs(nchans))
!ob_info % yobs = 0.0

end subroutine ufo_rttovonedvarcheck_InitObInfo

!-------------------------------------------------------------------------------

subroutine ufo_rttovonedvarcheck_GetBmatrix (fileunit, &
                                  bmatrix,        &
                                  b_elementsused, &
                                  fieldlist)

!-------------------------------------------------------------------------------
! read the input file and allocate and fill in all the components of the bmatrix
! structure, with the exception of proxy variables.
!
! notes:
!
! it is assumed a valid bmatrix file is available and has been opened ready for
! reading, accessed via the input unit number.
!
! two optional arguments are provided to allow a submatrix to be formed from the
! original file, either depending on a list of fields to be used for retrieval,
! or a list of specific element numbers to be retained.
!
! the file header may begin with any number of comment lines, which are defined
! as those using either # or ! as the first non-blank character.
!
! immediately following any comments, the header should then contain information
! on the size of the matrix to be read in, and a description of each matrix
! element. memory will be allocated depending on these matrix specifications.
!
! each matrix should then have a one line header containing the latitude bounds
! and a latitude band id. band id numbers should be sequential in latitude. i.e.
! band 1 will begin at -90.0, band 2 from -60.0 ... (or however wide we choose
! the bands).
!
! we also calculate an inverse of each b matrix, used in the current 1d-var for
! cost function monitoring only but it may be required for other minimization
! methods at some point.
!
! standard deviations are also stored in a separate vector.
!
! proxy variables are provided to allow manipulation of the chosen matrix during
! processing without affecting the original. we might wish to do this, for
! example, to take account of the different surface temperature errors for land
! and sea. this routine makes no assumption about the use of these variables,
! hence no space is allocated. nullification of unused pointers should take
! place outside (use the ops_satrad_initbmatrix routine).
!-------------------------------------------------------------------------------


implicit none

! subroutine arguments:
integer, intent(in)                :: fileunit
type (bmatrix_type), intent(inout) :: bmatrix
integer, optional, intent(in)      :: b_elementsused(:)
integer, optional, intent(inout)   :: fieldlist(:)

! local declarations:
character(len=*), parameter        :: routinename = "ufo_rttovonedvarcheck_GetBmatrix"
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
real                               :: southlimit
real                               :: northlimit
character(len=80)                  :: line
character(len=3)                   :: fieldtype
real, allocatable                  :: bfromfile(:,:)
logical, allocatable               :: list(:)
integer, allocatable               :: bfields(:,:)
integer, allocatable               :: elementsused(:)
logical                            :: rttovonedvarcheck_debug = .true.

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

if (rttovonedvarcheck_debug) then
  write (*, '(a)') 'reading b matrix file:'
  write (*, '(a,i0)') 'number of latitude bands = ', nbands
  write (*, '(a,i0)') 'matrix size = ', matrixsize
  if (nbfields > 0) then
    write (*, '(a)') 'order of fields and number of elements in each:'
    do i = 1, nbfields
      write (*, '(i0,a)') bfields(i,2), ' x ' // fieldtype_text(bfields(i,1))
    end do
  end if
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

allocate (bmatrix % fields(nfieldtypes,2))
bmatrix % fields(:,:) = 0
bmatrix % status = .true.
nelements = 0
nelements_total = 0
bmatrix % nfields = 0

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

  if (rttovonedvarcheck_debug .and. any (fieldlist /= 0)) then
    write (*, '(a)') 'the following requested retrieval fields are not in the b matrix:'
    do i = 1, size (fieldlist)
      if (fieldlist(i) /= 0) then
        write (fieldtype, '(i0)') fieldlist(i)
        if (fieldlist(i) > 0 .and. fieldlist(i) <= nfieldtypes) then
          write (*, '(a)') trim (fieldtype_text(fieldlist(i))) // &
            ' (fieldtype ' // trim (adjustl (fieldtype)) // ')'
        else
          write (*, '(a)') 'fieldtype ' // trim (adjustl (fieldtype)) // &
            ' which is invalid'
        end if
      end if
    end do
  end if

  if (rttovonedvarcheck_debug) then
    write (*, '(a)') 'b matrix fields used to define the retrieval profile vector:'
    do i = 1, nbfields
      write (*, '(a)') fieldtype_text(bfields(i,1))
    end do
  end if

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

if (nbfields > 0) bmatrix % fields(1:nbfields,:) = bfields(1:nbfields,:)
bmatrix % nfields = nbfields

! initialize the rest

bmatrix % nbands = nbands
allocate (bmatrix % store(nelements,nelements,nbands))
allocate (bmatrix % south(nbands))
allocate (bmatrix % north(nbands))
bmatrix % store(:,:,:) = 0.0

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

  if (rttovonedvarcheck_debug) then
    write (*, '(a,i0,a,2f8.2)') &
      'band no.', band, ' has southern and northern latitude limits of', &
      southlimit, northlimit
  end if

  !----
  ! 3.2) transfer into 1d-var arrays
  !----

  if (band > 0 .and. band <= nbands) then
    list(band) = .true.
    nmatrix = nmatrix + 1
    if (nelements < matrixsize) then
      bmatrix % store(:,:,band) = &
        bfromfile(elementsused(1:nelements),elementsused(1:nelements))
    else
      bmatrix % store(:,:,band) = bfromfile(:,:)
    end if
    bmatrix % south(band) = southlimit
    bmatrix % north(band) = northlimit
  else if (rttovonedvarcheck_debug) then
    write (*, '(a,i0)') 'skipped matrix with band number ', band
  end if

end do readallb

!---------------------
! 4. store error vector
!---------------------

allocate (bmatrix % sigma(nelements,nbands))
do i = 1, nelements
  bmatrix % sigma(i,:) = sqrt (bmatrix % store(i,i,:))
end do

!----------------
! 5. invert matrix
!----------------

allocate (bmatrix % inverse(nelements,nelements,nbands))
bmatrix % inverse(:,:,:) = bmatrix % store(:,:,:)

do k = 1, nbands
  call ufo_rttovonedvarcheck_InvertMatrix (nelements,                & ! in
                                      nelements,                & ! in
                                      bmatrix % inverse(:,:,k), & ! inout
                                      status)                     ! out
  if (status /= 0) then
    write(*,*) routinename // ' : b matrix is not invertible'
  end if
end do

!----------
! 6. tidy up
!----------

if (nmatrix < nbands) then
  write(*,*) routinename // ' : too few b matrices found in input file'
end if

end subroutine ufo_rttovonedvarcheck_GetBmatrix

subroutine ufo_rttovonedvarcheck_MapProfileToB( &
  bmatrix,  & ! in
  profinfo, & ! inout
  nelements ) ! out

!-------------------------------------------------------------------------------
! use b matrix elements to define a 1d-var profile vector. the profinfo data
! structure contains variables specifying the locations of all possible fields.
!
! in order to allow as much flexibility as possible, the profile vector is
! undefined at the beginning of processing so that different combinations of
! fields can be used. the essential criterion for deciding what should be
! present is that we have suitable background error covariances and therefore
! the form of the bmatrix can be used to construct the 1d-var profile. when the
! matrix is read in, we keep a record of the fields that are present, and we use
! that record in here.
!-------------------------------------------------------------------------------

implicit none

! subroutine arguments:
type(bmatrix_type),     intent(in)    :: bmatrix   !background error covariances
type(profileinfo_type), intent(inout) :: profinfo  !profile information
integer,                intent(out)   :: nelements !size of profile vector

! local constants:
character(len=*), parameter :: routinename = "ufo_rttovonedvarcheck_MapProfileToB"

! local variables:
integer :: i,j
integer :: firstelement
integer :: lastelement

!loop through fields in bmatrix. the number of elements that the field is composed of
!is held in bmatrix % fields(i,2).
nelements = 0

do j = 1, bmatrix % nfields
  firstelement = nelements + 1
  lastelement  = nelements + bmatrix % fields(j,2)

  !assign start and end points. if the field wasn't found then assign a value of
  !zero, which indicates absence.

  select case( bmatrix % fields(j,1) )

   !----------
   !atmosphere (set start and end points for multi-level fields)
   !----------

    case( fieldtype_t )
      profinfo % t(1)         = firstelement
      profinfo % t(2)         = lastelement

    case( fieldtype_q )
      profinfo % q(1)         = firstelement
      profinfo % q(2)         = lastelement

    case( fieldtype_ql )
      profinfo % ql(1)        = firstelement
      profinfo % ql(2)        = lastelement

    case( fieldtype_qi )
      profinfo % qi(1)        = firstelement
      profinfo % qi(2)        = lastelement

    case( fieldtype_cf )
      profinfo % cf(1)        = firstelement
      profinfo % cf(2)        = lastelement

    case( fieldtype_qt )
      profinfo % qt(1)        = firstelement
      profinfo % qt(2)        = lastelement

    case( fieldtype_o3profile )
      profinfo % o3profile(1) = firstelement
      profinfo % o3profile(2) = lastelement
      if ( ozonesource /= source_bg .and. profinfo % o3profile(1) > 0 ) then
        write(*,*) 'o3 profile retrieval requested but no cx % ozone available'
        write(*,*) 'o3 profile retrieval may result in high rejection rate'
      end if

    case( fieldtype_o3total )
      profinfo % o3total      = firstelement

    case( fieldtype_lwp )
      profinfo % lwp          = firstelement

   !-------
   !surface
   !-------

    case( fieldtype_t2 )
      profinfo % t2         = firstelement

    case( fieldtype_q2 )
      profinfo % q2         = firstelement

    case( fieldtype_tstar )
      profinfo % tstar      = firstelement

    case( fieldtype_pstar )
      profinfo % pstar      = firstelement

    case( fieldtype_windspeed )
      profinfo % windspeed  = firstelement

    case( fieldtype_mwemiss )
      profinfo % mwemiss(1) = firstelement
      profinfo % mwemiss(2) = lastelement

    case( fieldtype_emisspc )
      profinfo % emisspc(1) = firstelement
      profinfo % emisspc(2) = lastelement

   !------------------------------------
   !cloud (single-level grey cloud only)
   !------------------------------------

    case( fieldtype_cloudtopp )
      profinfo % cloudtopp         = firstelement

    case( fieldtype_cloudfrac )
      profinfo % cloudfrac         = firstelement

    case( fieldtype_not_used ) ! currently unused
      continue

    case default
      write(*,*) 'invalid field type in b matrix file: ',j
      cycle

  end select

  if ( firstelement /= 0 ) nelements = nelements + bmatrix % fields(j,2)

end do

end subroutine ufo_rttovonedvarcheck_MapProfileToB

!===============================================================================
! subroutines
! order is inportant at the moment will look into interface blocks
!===============================================================================

subroutine ufo_rttovonedvarcheck_InvertMatrix (n,      &
                                          m,      &
                                          a,      &
                                          status, &
                                          matrix)

!-------------------------------------------------------------------------------
! invert a matrix and optionally premultiply by another matrix.
!
! variables with intent in:
!
!     n:  size of the matrix being inverted
!     m:  if matrix is not present this is the same as n, else this is
!         the other dimension of matrix.
!
! variables with intent inout:
!
!     a:               real matrix (assumed square and symmetrical)
!                      overwritten by its inverse if matrix is not
!                      present.
!
! variables with intent out:
!
!     status:          0: ok, 1: a is not positive definite.
!
! optional variables with intent inout:
!
!     matrix:          if present this input matrix is replaced by
!                      (matrix).a^-1 on exit (leaving a unchanged).
!
! uses cholesky decomposition - a method particularly suitable for real
! symmetric matrices.
!
! cholesky decomposition solves the linear equation uq=v for q where u is a
! symmetric positive definite matrix and u and q are vectors of length n.
!
! the method follows that in golub and van loan although this is pretty
! standard.
!
! if u is not positive definite this will be detected by the program and flagged
! as an error.  u is assumed to be symmetric as only the upper triangle is in
! fact used.
!
! if the the optional parameter matrix is present, it is replaced by
! (matrix).a^-1 on exit and a is left unchanged.
!-------------------------------------------------------------------------------

implicit none

! subroutine arguments:
integer, intent(in)           :: n           ! order of a
integer, intent(in)           :: m           ! order of optional matrix, if required
real, intent(inout)           :: a(n,n)      ! square mx, overwritten by its inverse
integer, intent(out)          :: status      ! 0 if all ok 1 if matrix is not positive definite
real, optional, intent(inout) :: matrix(n,m) ! replaced by (matrix).a^-1 on exit

! local declarations:
character(len=*), parameter   :: routinename = "ufo_rttovonedvarcheck_InvertMatrix"
character(len=80)             :: errormessage(2)
real, parameter               :: tolerance = tiny (0.0) * 100.0
integer                       :: i
integer                       :: j
integer                       :: k
integer                       :: mm
real                          :: g(n,n)      ! the cholesky triangle matrix
real                          :: q(n)
real                          :: tmp(n,m)
real                          :: v(n)
real                          :: x(n)        ! temporary array

status = 0

! determine the cholesky triangle matrix.

do j = 1, n
  x(j:n) = a(j:n,j)
  if (j /= 1) then
    do k = 1, j - 1
      x(j:n) = x(j:n) - g(j,k) * g(j:n,k)
    end do
  end if
  if (x(j) <= tolerance) then
    write(*,*) routinename // ' : matrix is not positive definite'
    status = 1
    goto 9999
  end if
  g(j:n,j) = x(j:n) / sqrt (x(j))
end do

! now solve the equation g.g^t.q=v for the set of
! vectors, v, with one element = 1 and the rest zero.
! the solutions q are brought together at the end to form
! the inverse of g.g^t (i.e., the inverse of a).

if (present (matrix)) then
   mm = m
else

  ! make sure that the dimensions of tmp were correctly specified

  if (m /= n) then
    errormessage(1) = '2nd and 3rd arguments of routine must be'
    errormessage(2) = 'identical if the matrix option is not present'
    write(*,*) routinename // errormessage(1:2)
    status = 2
    goto 9999
  end if
  mm = n
end if

main_loop : do i = 1, mm
  if (.not. (present (matrix))) then
    v(:) = 0.0
    v(i) = 1.0
  else
    v(:) = matrix(i,:)
  end if

  ! solve gx=v for x by forward substitution

  x = v
  x(1) = x(1) / g(1,1)
  do j = 2, n
    x(j) = (x(j) - dot_product (g(j,1:j - 1), x(1:j - 1))) / g(j,j)
  end do

  ! solve g^t.q=x for q by backward substitution

  q = x
  q(n) = q(n) / g(n,n)
  do j = n - 1, 1, -1
    q(j) = (q(j) - dot_product (g(j + 1:n,j), q(j + 1:n))) / g(j,j)
  end do

  tmp(:,i) = q(:)

end do main_loop

if (.not. (present (matrix))) then
  a(:,:) = tmp(:,:)
else
  matrix(:,:) = tmp(:,:)
end if

9999 continue

end subroutine ufo_rttovonedvarcheck_InvertMatrix

!------------------------------------------------------------------------------

end module ufo_rttovonedvarcheck_init_mod
