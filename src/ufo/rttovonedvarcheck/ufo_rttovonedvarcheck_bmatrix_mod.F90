! (C) British Crown Copyright 2017-2018 Met Office
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 

module ufo_rttovonedvarcheck_bmatrix_mod

use kinds
use iso_c_binding
use fckit_configuration_module, only: fckit_configuration
use fckit_log_module, only : fckit_log
use ufo_rttovonedvarcheck_utils_mod, only: &
        max_string, &
        fieldtype_text, &
        nfieldtypes

implicit none

private

! public types and subroutines

! private list for easy naviagtion
! onedvarcheck_covariance_InitBmatrix
! onedvarcheck_covariance_GetBmatrix
! onedvarcheck_covariance_InvertMatrix
! onedvarcheck_covariance_iogetfreeunit
! public ufo_rttovonedvarcheck_bmatrix_setup
! public ufo_rttovonedvarcheck_bmatrix_delete

type, public :: bmatrix_type
  !private
  logical :: status                                 ! status indicator
  integer :: nbands                                 ! number of latitude bands
  integer :: nsurf                                  ! number of surface type variations
  integer :: nfields                                ! number of fields
  integer, pointer      :: fields(:,:)              ! fieldtypes and no. elements in each
  real(kind=kind_real), pointer :: store(:,:,:)     ! original b-matrices read from the file
  real(kind=kind_real), pointer :: inverse(:,:,:)   ! inverse of above
  real(kind=kind_real), pointer :: sigma(:,:)       ! diagonal elements
  real(kind=kind_real), pointer :: proxy(:,:)       ! copy of original for manipulation
  real(kind=kind_real), pointer :: inv_proxy(:,:)   ! copy of inverse
  real(kind=kind_real), pointer :: sigma_proxy(:,:) ! copy of diagonal
  real(kind=kind_real), pointer :: south(:)         ! s limit of each latitude band
  real(kind=kind_real), pointer :: north(:)         ! n limit of each latitude band
contains
  procedure :: setup  => ufo_rttovonedvarcheck_bmatrix_setup
  procedure :: delete => ufo_rttovonedvarcheck_bmatrix_delete
end type bmatrix_type

contains

! ------------------------------------------------------------------------------

subroutine ufo_rttovonedvarcheck_bmatrix_setup(self, variables, filepath, qtotal_flag)

implicit none
class(bmatrix_type), intent(inout) :: self         !< Covariance structure
character(len=*), intent(in)      :: variables(:)  !< Model variables in B matrix
character(len=*), intent(in)      :: filepath      !< Path to B matrix file
logical                           :: qtotal_flag   !< Flag for qtotal

logical                       :: file_exists  !< Check if a file exists logical
integer                       :: fileunit     !< Unit number for reading in files
integer, allocatable          :: fields_in(:)
integer                       :: nx, ny, jvar
integer                       :: nmvars
character(len=max_string)     :: varname
real(kind=kind_real)          :: t1,t2
logical                       :: flag
character(len=max_string)     :: message
character(len=:), allocatable :: str

call fckit_log % info("ufo_rttovonedvarcheck_bmatrix_setup start")

call CPU_TIME(t1)

! Which fields are being used from b matrix file - temporary until rttov can handle all fields
allocate(fields_in(8))
fields_in(:) = 0

nmvars = size(variables)
do jvar = 1, nmvars

  varname = variables(jvar)

  select case (trim(varname))

    case ("air_temperature")
      fields_in(1) = 1 ! air_temperature

    case ("specific_humidity")
      if (qtotal_flag) then
        fields_in(2) = 10 ! total water profile in bmatrix - specific humidity in geovals!?!
      else
        fields_in(2) = 2 ! water profile in bmatrix - specific humidity in geovals!?!
      end if

    case("air_temperature_at_two_meters_above_surface")
      fields_in(3) = 3 ! 2m air_temperature

    case("specific_humidity_at_two_meters_above_surface")
      fields_in(4) = 4 ! 2m specific_humidity

    case("skin_temperature")
      fields_in(5) = 5 ! surface skin temperature

    case("surface_air_pressure")
      fields_in(6) = 6 ! surface air pressure

    case ("mass_content_of_cloud_ice_in_atmosphere_layer")
      if (.NOT. qtotal_flag) fields_in(7) = 15 ! solid water profile

    case ("mass_content_of_cloud_liquid_water_in_atmosphere_layer")
      if (.NOT. qtotal_flag) fields_in(8) = 9  ! liquid water profile

    case default
      call abor1_ftn("Variable not implemented yet in OneDVarCheck Covariance")

  end select

end do

inquire(file=trim(filepath), exist=file_exists)
if (file_exists) then
  call onedvarcheck_covariance_IOGetFreeUnit(fileunit)
  open(unit = fileunit, file = trim(filepath))
  call onedvarcheck_covariance_InitBmatrix(self)
  !call onedvarcheck_covariance_GetBmatrix(fileunit, self, fieldlist=fields_in) ! used for development
  call onedvarcheck_covariance_GetBmatrix(fileunit, self)
  close(unit = fileunit)
  call fckit_log % info("onedvarcheck bmatrix file exists and read in")
else
  call abor1_ftn("onedvarcheck bmatrix file not found")
end if

nx = size(self % store,1)
ny = size(self % store,2)

call cpu_time(t2)

write(message,*) "ufo_rttovonedvarcheck_bmatrix_setup cpu time = ",(t2-t1)
call fckit_log % info(message)

end subroutine ufo_rttovonedvarcheck_bmatrix_setup

! ------------------------------------------------------------------------------

subroutine ufo_rttovonedvarcheck_bmatrix_delete(self)
implicit none
class(bmatrix_type), intent(inout) :: self  !< Covariance structure

end subroutine ufo_rttovonedvarcheck_bmatrix_delete

! ------------------------------------------------------------------------------

subroutine onedvarcheck_covariance_InitBmatrix(bmatrix)

! nullify B matrix pointers.

implicit none

! subroutine arguments:
type(bmatrix_type), intent(out) :: bmatrix
character(len=*), parameter     :: routinename = "onedvarcheck_covariance_InitBmatrix"

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

end subroutine onedvarcheck_covariance_InitBmatrix

!-------------------------------------------------------------------------------

subroutine onedvarcheck_covariance_GetBmatrix (fileunit, &
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
character(len=*), parameter        :: routinename = "onedvarcheck_covariance_GetBmatrix"
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
logical                            :: onedvarcheck_debug = .true.

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

if (onedvarcheck_debug) then
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

  if (onedvarcheck_debug .and. any (fieldlist /= 0)) then
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

  if (onedvarcheck_debug) then
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
bmatrix % store(:,:,:) = 0.0_kind_real

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

  if (onedvarcheck_debug) then
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
  else if (onedvarcheck_debug) then
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
  call onedvarcheck_covariance_InvertMatrix (nelements,          & ! in
                                             nelements,                 & ! in
                                             bmatrix % inverse(:,:,k),  & ! inout
                                             status)                      ! out
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

end subroutine onedvarcheck_covariance_GetBmatrix

!-------------------------------------------------------------------------------

subroutine onedvarcheck_covariance_InvertMatrix (n,      &
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
integer, intent(in)                           :: n           ! order of a
integer, intent(in)                           :: m           ! order of optional matrix, if required
real(kind=kind_real), intent(inout)           :: a(n,n)      ! square mx, overwritten by its inverse
integer, intent(out)                          :: status      ! 0 if all ok 1 if matrix is not positive definite
real(kind=kind_real), optional, intent(inout) :: matrix(n,m) ! replaced by (matrix).a^-1 on exit

! local declarations:
character(len=*), parameter   :: routinename = "onedvarcheck_covariance_InvertMatrix"
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

end subroutine onedvarcheck_covariance_InvertMatrix

!-------------------------------------------------------------------------------

subroutine onedvarcheck_covariance_iogetfreeunit(unit)

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

end subroutine onedvarcheck_covariance_iogetfreeunit

! ------------------------------------------------------------------------------

end module ufo_rttovonedvarcheck_bmatrix_mod
