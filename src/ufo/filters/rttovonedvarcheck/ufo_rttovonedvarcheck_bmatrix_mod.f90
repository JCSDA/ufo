! (C) British Crown Copyright 2020 Met Office
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 

!> Fortran module containing the full b-matrix data type and methods for the 1D-Var.

module ufo_rttovonedvarcheck_bmatrix_mod

use fckit_log_module, only : fckit_log
use kinds
use ufo_rttovonedvarcheck_constants_mod
use ufo_rttovonedvarcheck_utils_mod
use ufo_vars_mod

implicit none
private

type, public :: ufo_rttovonedvarcheck_bmatrix
  logical :: status                                 !< status indicator
  logical :: debug                                  !< flag for printing verbose output
  integer :: nbands                                 !< number of latitude bands
  integer :: nsurf                                  !< number of surface type variations
  integer :: nfields                                !< number of fields
  integer, pointer      :: fields(:,:)              !< fieldtypes and no. elements in each
  real(kind=kind_real), pointer :: store(:,:,:)     !< original b-matrices read from the file
  real(kind=kind_real), pointer :: inverse(:,:,:)   !< inverse of above
  real(kind=kind_real), pointer :: sigma(:,:)       !< diagonal elements
  real(kind=kind_real), pointer :: proxy(:,:)       !< copy of original for manipulation
  real(kind=kind_real), pointer :: inv_proxy(:,:)   !< copy of inverse
  real(kind=kind_real), pointer :: sigma_proxy(:,:) !< copy of diagonal
  real(kind=kind_real), pointer :: south(:)         !< s limit of each latitude band
  real(kind=kind_real), pointer :: north(:)         !< n limit of each latitude band
contains
  procedure :: setup  => ufo_rttovonedvarcheck_bmatrix_setup
  procedure :: delete => ufo_rttovonedvarcheck_bmatrix_delete
  procedure :: reset  => ufo_rttovonedvarcheck_reset_covariances
end type ufo_rttovonedvarcheck_bmatrix

contains

! ------------------------------------------------------------------------------------------------
!> Routine to read and setup the 1D-Var B-matrix
!!
!! \author Met Office
!!
!! \date 09/06/2020: Created
!!
subroutine ufo_rttovonedvarcheck_bmatrix_setup(self, variables, filepath, qtotal_flag)

implicit none
class(ufo_rttovonedvarcheck_bmatrix), intent(inout) :: self !< B-matrix Covariance
character(len=*), intent(in)       :: variables(:) !< Model variables in B matrix
character(len=*), intent(in)       :: filepath     !< Path to B matrix file
logical, intent(in)                :: qtotal_flag  !< Flag for qtotal

logical                       :: file_exists  ! Check if a file exists logical
integer                       :: fileunit     ! Unit number for reading in files
integer, allocatable          :: fields_in(:) ! Fields_in used to subset b-matrix for testing.
real(kind=kind_real)          :: t1,t2        ! Time values for logging
character(len=max_string)     :: message
character(len=:), allocatable :: str
logical                       :: testing = .true.

call fckit_log % info("ufo_rttovonedvarcheck_bmatrix_setup start")

call CPU_TIME(t1)

! Open file and read in b-matrix
inquire(file=trim(filepath), exist=file_exists)
if (file_exists) then
  call ufo_rttovonedvarcheck_iogetfreeunit(fileunit)
  open(unit = fileunit, file = trim(filepath))
  call rttovonedvarcheck_covariance_InitBmatrix(self)
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

self % debug = .true.

call cpu_time(t2)

write(message,*) "ufo_rttovonedvarcheck_bmatrix_setup cpu time = ",(t2-t1)
call fckit_log % info(message)

end subroutine ufo_rttovonedvarcheck_bmatrix_setup

! ------------------------------------------------------------------------------------------------
!> Routine to delete and squash the 1D-Var B-matrix
!!
!! \details Met Office OPS Heritage: Ops_SatRad_SquashBmatrix.f90
!!
!! \author Met Office
!!
!! \date 09/06/2020: Created
!!
subroutine ufo_rttovonedvarcheck_bmatrix_delete(self)

implicit none
class(ufo_rttovonedvarcheck_bmatrix), intent(inout) :: self  !< B-matrix Covariance

character(len=*), parameter :: RoutineName = "ufo_rttovonedvarcheck_bmatrix_delete"

self % status = .false.
self % debug = .false.
self % nbands = 0
self % nsurf = 0
if ( associated(self % fields)      ) deallocate( self % fields      )
if ( associated(self % store)       ) deallocate( self % store       )
if ( associated(self % inverse)     ) deallocate( self % inverse     )
if ( associated(self % sigma)       ) deallocate( self % sigma       )
if ( associated(self % proxy)       ) deallocate( self % proxy       )
if ( associated(self % inv_proxy)   ) deallocate( self % inv_proxy   )
if ( associated(self % sigma_proxy) ) deallocate( self % sigma_proxy )
if ( associated(self % south)       ) deallocate( self % south       )
if ( associated(self % north)       ) deallocate( self % north       )

end subroutine ufo_rttovonedvarcheck_bmatrix_delete

! ------------------------------------------------------------------------------------------------
!> \brief Routine to initialize the 1D-Var B-matrix
!!
!! \details Met Office OPS Heritage: Ops_SatRad_InitBmatrix.f90
!!
!! \author Met Office
!!
!! \date 09/06/2020: Created
!!
subroutine rttovonedvarcheck_covariance_InitBmatrix(self)

implicit none

! subroutine arguments:
type(ufo_rttovonedvarcheck_bmatrix), intent(out) :: self !< B-matrix Covariance 

character(len=*), parameter     :: routinename = "rttovonedvarcheck_covariance_InitBmatrix"

self % status = .false.
self % debug = .false.
self % nbands = 0
self % nsurf = 0
nullify( self % fields      )
nullify( self % store       )
nullify( self % inverse     )
nullify( self % sigma       )
nullify( self % proxy       )
nullify( self % inv_proxy   )
nullify( self % sigma_proxy )
nullify( self % south       )
nullify( self % north       )

end subroutine rttovonedvarcheck_covariance_InitBmatrix

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
!! hence no space is allocated. nullification of unused pointers should take
!! place outside (use the ops_satrad_initbmatrix routine).
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
type (ufo_rttovonedvarcheck_bmatrix), intent(inout) :: self !< B-matrix covariance
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

if (self % debug) then
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

allocate (self % fields(nfieldtypes,2))
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

  if (self % debug .and. any (fieldlist /= 0)) then
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

  if (self % debug) then
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

if (nbfields > 0) self % fields(1:nbfields,:) = bfields(1:nbfields,:)
self % nfields = nbfields

! initialize the rest

self % nbands = nbands
allocate (self % store(nelements,nelements,nbands))
allocate (self % south(nbands))
allocate (self % north(nbands))
self % store(:,:,:) = 0.0_kind_real

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

  if (self % debug) then
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
      self % store(:,:,band) = &
        bfromfile(elementsused(1:nelements),elementsused(1:nelements))
    else
      self % store(:,:,band) = bfromfile(:,:)
    end if
    self % south(band) = southlimit
    self % north(band) = northlimit
  else if (self % debug) then
    write (*, '(a,i0)') 'skipped matrix with band number ', band
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
  call rttovonedvarcheck_covariance_InvertMatrix (nelements,       & ! in
                                           nelements,              & ! in
                                           self % inverse(:,:,k),  & ! inout
                                           status)                   ! out
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
subroutine ufo_rttovonedvarcheck_reset_covariances(self, latitude, & ! in
                                      b_matrix, b_inverse, b_sigma ) ! out

implicit none

! Subroutine arguments
class(ufo_rttovonedvarcheck_bmatrix), intent(in) :: self !< B-matrix covariance
real(kind_real), intent(in)  :: latitude
real(kind_real), intent(out) :: b_matrix(:,:)
real(kind_real), intent(out) :: b_inverse(:,:)
real(kind_real), intent(out) :: b_sigma(:)

! Local Variables
integer :: band, i

! select appropriate b matrix for latitude of observation
b_matrix(:,:) = 0.0
b_inverse(:,:) = 0.0
b_sigma(:) = 0.0
do band = 1, self % nbands
  if (latitude < self % north(band)) exit
end do
b_matrix(:,:) = self % store(:,:,band)
b_inverse(:,:) = self % inverse(:,:,band)
b_sigma(:) = self % sigma(:,band)

! This has been left in for future development
!! Use errors associated with microwave emissivity atlas
!if (profindex % mwemiss(1) > 0) then
!
!  ! Atlas uncertainty stored in Ob % MwEmErrAtlas, use this to scale each
!  ! row/column of the B matrix block.
!  ! The default B matrix, see e.g. file ATOVS_Bmatrix_43, contains error
!  ! covariances representing a global average. Here, those elements are scaled
!  ! by a factor MwEmissError/SQRT(diag(B_matrix)) for each channel.
!
!  do i = profindex % mwemiss(1), profindex % mwemiss(2)
!    MwEmissError = Ob % MwEmErrAtlas(i - profindex % mwemiss(1) + 1)
!    if (MwEmissError > 1.0E-4 .and. MwEmissError < 1.0) then
!      bscale = MwEmissError / sqrt (B_matrix(i,i))
!      B_matrix(:,i) = B_matrix(:,i) * bscale
!      B_matrix(i,:) = B_matrix(i,:) * bscale
!      B_inverse(:,i) = B_inverse(:,i) / bscale
!      B_inverse(i,:) = B_inverse(i,:) / bscale
!      B_sigma(i) = B_sigma(i) * bscale
!    end if
!  end do
!
!end if
!
!! Scale the background skin temperature error covariances over land
!if (Ob % Surface == RTland .and. SkinTempErrorLand >= 0.0) then
!
!  bscale = SkinTempErrorLand / sqrt (B_matrix(profindex % tstar,profindex % tstar))
!  B_matrix(:,profindex % tstar) = B_matrix(:,profindex % tstar) * bscale
!  B_matrix(profindex % tstar,:) = B_matrix(profindex % tstar,:) * bscale
!  B_inverse(:,profindex % tstar) = B_inverse(:,profindex % tstar) / bscale
!  B_inverse(profindex % tstar,:) = B_inverse(profindex % tstar,:) / bscale
!  B_sigma(profindex % tstar) = B_sigma(profindex % tstar) * bscale
!
!end if

end subroutine ufo_rttovonedvarcheck_reset_covariances

! ------------------------------------------------------------------------------------------------
!> Routine to invert the 1D-Var B-matrix
!!
!! \details Met Office OPS Heritage: Ops_SatRad_InvertMatrix.f90
!!
!! invert a matrix and optionally premultiply by another matrix.
!!
!! variables with intent in:
!!
!!     n:  size of the matrix being inverted <br>
!!     m:  if matrix is not present this is the same as n, else this is
!!         the other dimension of matrix.
!!
!! variables with intent inout:
!!
!!     a:               real matrix (assumed square and symmetrical)
!!                      overwritten by its inverse if matrix is not
!!                      present.
!!
!! variables with intent out:
!!
!!     status:          0: ok, 1: a is not positive definite.
!!
!! optional variables with intent inout:
!!
!!     matrix:          if present this input matrix is replaced by
!!                      (matrix).a^-1 on exit (leaving a unchanged).
!!
!! uses cholesky decomposition - a method particularly suitable for real
!! symmetric matrices.
!!
!! cholesky decomposition solves the linear equation uq=v for q where u is a
!! symmetric positive definite matrix and u and q are vectors of length n.
!!
!! the method follows that in golub and van loan although this is pretty
!! standard.
!!
!! if u is not positive definite this will be detected by the program and flagged
!! as an error.  u is assumed to be symmetric as only the upper triangle is in
!! fact used.
!!
!! if the the optional parameter matrix is present, it is replaced by
!! (matrix).a^-1 on exit and a is left unchanged.
!!
!! \author Met Office
!!
!! \date 09/06/2020: Created
!!
subroutine rttovonedvarcheck_covariance_InvertMatrix (n,      &
                                                 m,      &
                                                 a,      &
                                                 status, &
                                                 matrix)

implicit none

! subroutine arguments:
integer, intent(in)                           :: n           !< order of a
integer, intent(in)                           :: m           !< order of optional matrix, if required
real(kind=kind_real), intent(inout)           :: a(n,n)      !< square mx, overwritten by its inverse
integer, intent(out)                          :: status      !< 0 if all ok 1 if matrix is not positive definite
real(kind=kind_real), optional, intent(inout) :: matrix(n,m) !< replaced by (matrix).a^-1 on exit

! local declarations:
character(len=*), parameter   :: routinename = "rttovonedvarcheck_covariance_InvertMatrix"
character(len=80)             :: errormessage(2)
real(kind=kind_real), parameter :: tolerance = tiny (0.0_kind_real) * 100.0_kind_real
integer                       :: i
integer                       :: j
integer                       :: k
integer                       :: mm
real(kind=kind_real)          :: g(n,n)      ! the cholesky triangle matrix
real(kind=kind_real)          :: q(n)
real(kind=kind_real)          :: tmp(n,m)
real(kind=kind_real)          :: v(n)
real(kind=kind_real)          :: x(n)        ! temporary array

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
    v(:) = 0.0_kind_real
    v(i) = 1.0_kind_real
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

end subroutine rttovonedvarcheck_covariance_InvertMatrix

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

character(len=max_string) :: varname
integer                   :: jvar
integer                   :: nmvars
character(len=max_string) :: message

call fckit_log % info("rttovonedvarcheck_create_fields_in: starting")

! Which fields are being used from b matrix file - temporary until rttov can handle all fields
allocate(fields_in(19))
fields_in(:) = 0

nmvars = size(variables)
do jvar = 1, nmvars

  varname = variables(jvar)

  select case (trim(varname))

    case (var_ts)
      fields_in(1) = 1 ! air_temperature

    case (var_q)
      if (qtotal_flag) then
        fields_in(10) = 10 ! total water profile
      else
        fields_in(2) = 2 ! water profile
      end if

    case(var_sfc_t2m)
      fields_in(3) = 3 ! 2m air_temperature

    case(var_sfc_q2m)
      fields_in(4) = 4 ! 2m specific_humidity

    case(var_sfc_tskin)
      fields_in(5) = 5 ! surface skin temperature

    case(var_sfc_p2m)
      fields_in(6) = 6 ! surface air pressure

    ! 7 - o3total is not implmented yet
    ! 8 - not used is not implmented yet

    case (var_clw)
      if (.NOT. qtotal_flag) then 
        fields_in(9) = 9  ! liquid water profile
        call abor1_ftn("rttovonedvarcheck not setup for independent clw yet")
      end if

    case ("surface_wind_speed") ! surface wind speed
      fields_in(11) = 11
      call abor1_ftn("rttovonedvarcheck not setup for surface windspeed yet")

    ! 12 - o3profile is not implmented yet
    ! 13 - lwp (liquid water path) is not implmented yet

    case ("surface_emissivity") ! microwave emissivity
      fields_in(14) = 14

    case (var_cli)
      if (.NOT. qtotal_flag) then
        fields_in(15) = 15 ! ice profile
        call abor1_ftn("rttovonedvarcheck not setup for independent ciw yet")
      end if

    case ("cloud_top_pressure")
      fields_in(16) = 16
      call abor1_ftn("rttovonedvarcheck not setup for cloud retrievals yet")

    case ("effective_cloud_fraction") ! effective cloud fraction
      fields_in(17) = 17
      call abor1_ftn("rttovonedvarcheck not setup for cloud retrievals yet")

    case ("emissivity_pc") ! emissivity prinipal components
      fields_in(18) = 18
      call abor1_ftn("rttovonedvarcheck not setup for pc emissivity yet")

    ! 19 cloud fraction profile - not currently used

    case default
      write(message,*) trim(varname)," not implemented yet in rttovonedvarcheck Covariance"
      call abor1_ftn(message)

  end select

end do

end subroutine rttovonedvarcheck_create_fields_in

! ------------------------------------------------------------------------------------------------

end module ufo_rttovonedvarcheck_bmatrix_mod
