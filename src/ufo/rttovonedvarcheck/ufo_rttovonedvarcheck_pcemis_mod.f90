! (C) copyright 2020 Met Office
!
! this software is licensed under the terms of the apache licence version 2.0
! which can be obtained at http://www.apache.org/licenses/license-2.0.

!> Fortran module which contains the methods for Infrared
!> Principal Component Emissivity
!> The science can be found at:
!> Pavelin, E., Candy, B., 2014. Assimilation of surface-sensitive infrared radiances over land: 
!> Estimation of land surface temperature and emissivity. Q. J. R. Metorol. Soc., 140, 1198-1208.

module ufo_rttovonedvarcheck_pcemis_mod

use fckit_log_module, only : fckit_log
use kinds
use ufo_rttovonedvarcheck_constants_mod
use ufo_rttovonedvarcheck_utils_mod

implicit none
private

!< Emissivity eigen vector type definition
type, public :: ufo_rttovonedvarcheck_EmisEigenvec
   integer          :: NChans
   integer          :: NumEV
   integer, pointer :: Channels(:) => null()
   real, pointer    :: Mean(:) => null()
   real, pointer    :: PCmin(:) => null()
   real, pointer    :: PCmax(:) => null()
   real, pointer    :: PCguess(:) => null()
   real, pointer    :: EV(:,:) => null()
   real, pointer    :: EV_Inverse(:,:) => null()
end type ufo_rttovonedvarcheck_EmisEigenvec

!< Emissivity eigen vector atlas type definition
TYPE ufo_rttovonedvarcheck_EmisAtlas
   INTEGER          :: Nlat
   INTEGER          :: Nlon
   INTEGER          :: Npc
   REAL             :: gridstep
   REAL, POINTER    :: EmisPC(:,:,:)
END TYPE ufo_rttovonedvarcheck_EmisAtlas

!< Principal component emissivity type definition
type, public :: ufo_rttovonedvarcheck_pcemis

  type(ufo_rttovonedvarcheck_EmisEigenvec) :: emis_eigen
  type(ufo_rttovonedvarcheck_EmisAtlas)    :: emis_atlas
  logical :: initialised = .false.

contains
  procedure :: setup  => ufo_rttovonedvarcheck_InitPCemis
  procedure :: delete => ufo_rttovonedvarcheck_DeletePCemis
  procedure :: info   => ufo_rttovonedvarcheck_PrintPCemis
  procedure :: emistopc => ufo_rttovonedvarcheck_EmisToPC
  procedure :: pctoemis => ufo_rttovonedvarcheck_PCToEmis
  procedure :: emisktopc => ufo_rttovonedvarcheck_EmisKToPC

end type ufo_rttovonedvarcheck_pcemis

contains

!-------------------------------------------------------------------------------
!> Initialize PC emissivity object
!!
!! \author Met Office
!!
!! \date 04/08/2020: Created
!!
subroutine ufo_rttovonedvarcheck_InitPCEmis(self, filepath, atlaspath)

implicit none

! subroutine arguments:
class(ufo_rttovonedvarcheck_pcemis), intent(out) :: self !< PC emissivity type
character(len=*), intent(in) :: filepath
character(len=*), intent(in), optional :: atlaspath

character(len=*), parameter :: routinename = "ufo_rttovonedvarcheck_InitPCEmis"
logical                     :: file_exists  ! Check if a file exists logical
integer                     :: fileunit     ! Unit number for reading in files

! Read eigenvector file
inquire(file=trim(filepath), exist=file_exists)
if (file_exists) then
  call ufo_rttovonedvarcheck_iogetfreeunit(fileunit)
  open(unit = fileunit, file = trim(filepath))
  call ufo_rttovonedvarcheck_GetEmisEigenVec(self, fileunit)
  close(unit = fileunit)
  call fckit_log % info("rttovonedvarcheck EmisEigenVec file exists and read in")
else
  call abor1_ftn("rttovonedvarcheck EmisEigenVec file not found: aborting")
end if

self % initialised = .true.

! Read in emissivity atlas if file path present - 
! if not a first guess will be used from eigenvector file
if (present(atlaspath)) then
  inquire(file=trim(atlaspath), exist=file_exists)
  if (file_exists) then
    call ufo_rttovonedvarcheck_iogetfreeunit(fileunit)
    open(unit = fileunit, file = trim(filepath))
    call ufo_rttovonedvarcheck_GetEmisAtlas(self, fileunit)
    close(unit = fileunit)
    call fckit_log % info("rttovonedvarcheck Emis Atlas file exists and read in")
  else
    call abor1_ftn("rttovonedvarcheck Emis Atlas file not found but requested: aborting")
  end if
end if

call self % info()

end subroutine ufo_rttovonedvarcheck_InitPCEmis

!-------------------------------------------------------------------------------
!> Read the emissivity eigen vector from file
!!
!! Heritage: Ops_SatRad_GetEmisEigenVec
!!
!! \author Met Office
!!
!! \date 04/08/2020: Created
!!
subroutine ufo_rttovonedvarcheck_GetEmisEigenVec (self,     &
                                                  fileunit  )

implicit none

! Subroutine arguments:
type(ufo_rttovonedvarcheck_pcemis), intent(out) :: self !< PC emissivity type
integer, intent(in)                     :: fileunit

! Local declarations:
character(len=*), parameter :: RoutineName = "ufo_rttovonedvarcheck_GetEmisEigenVec"
integer :: eigversion
integer :: i
integer :: readstatus
character(len=max_string)   :: message

!----------------------------------------------
! 1. Read header information and allocate arrays
!----------------------------------------------

read(fileunit, *, iostat = readstatus) eigversion, self % emis_eigen % Nchans, &
                                        self % emis_eigen % NumEV

allocate (self % emis_eigen % Channels( self % emis_eigen % Nchans))
allocate (self % emis_eigen % Mean( self % emis_eigen % Nchans))
allocate (self % emis_eigen % PCmin( self % emis_eigen % NumEV))
allocate (self % emis_eigen % PCmax( self % emis_eigen % NumEV))
allocate (self % emis_eigen % PCguess( self % emis_eigen % NumEV))
allocate (self % emis_eigen % EV( self % emis_eigen % NumEV, self % emis_eigen % Nchans))
if (eigversion >= 2) then
  allocate (self % emis_eigen % EV_Inverse( self % emis_eigen % Nchans, self % emis_eigen % NumEV))
end if

!--------------------------------------------------------
! 2. Read the channels, mean emissivities and eigenvectors
!--------------------------------------------------------

read (fileunit, *, iostat = readstatus) self % emis_eigen % Channels(:)
read (fileunit, *, iostat = readstatus) self % emis_eigen % Mean(:)
read (fileunit, *, iostat = readstatus) self % emis_eigen % PCmin(:)
read (fileunit, *, iostat = readstatus) self % emis_eigen % PCmax(:)
read (fileunit, *, iostat = readstatus) self % emis_eigen % PCguess(:)
do i = 1, self % emis_eigen % NumEV
  read (fileunit, *, iostat = readstatus) self % emis_eigen % EV(i,:)
end do

! Has there been an error in the read?
if (readstatus /= 0) then
  write(message,*) RoutineName,  &
       'Problem reading in emis eigenvectors - please check the file '
  call abor1_ftn(message)
end if

if (eigversion >= 2) then
  do i = 1, self % emis_eigen % Nchans
    read (fileunit, *, iostat = readstatus) self % emis_eigen % EV_Inverse(i,:)
  end do
  ! Has there been an error in the read?
  if (readstatus /= 0) then
    write(message,*) RoutineName,  &
         'Problem reading in inverse emis eigenvectors - please check the file '
    call abor1_ftn(message)
  end if
end if

write(*, '(A,I0,A,I0,A)') 'Finished reading ',self % emis_eigen % NumEV, &
                          ' emissivity eigenvectors on ', &
                           self % emis_eigen % Nchans,' channels.'

end subroutine ufo_rttovonedvarcheck_GetEmisEigenVec

!-------------------------------------------------------------------------------
!> Read the emissivity eigen atlas from file
!!
!! Heritage: Ops_SatRad_GetEmisAtlas
!!
!! \author Met Office
!!
!! \date 16/09/2020: Created
!!
subroutine ufo_rttovonedvarcheck_GetEmisAtlas (self, fileunit)

implicit none

! Subroutine arguments:
type(ufo_rttovonedvarcheck_pcemis), intent(out) :: self !< PC emissivity type
integer, intent(in) :: fileunit

! Local declarations:
character(len=*), parameter          :: RoutineName = "ufo_rttovonedvarcheck_GetEmisAtlas"
character(len=max_string)            :: message
integer                              :: readstatus
integer                              :: i
integer                              :: j

!----------------------------------------------
! 1. Read header information and allocate arrays
!----------------------------------------------

read (fileunit, *, iostat = readstatus) self % emis_atlas % Nlat, &
                                        self % emis_atlas % Nlon, &
                                        self % emis_atlas % Npc, &
                                        self % emis_atlas % gridstep

allocate (self % emis_atlas % EmisPC(self % emis_atlas % Nlon, &
                                     self % emis_atlas % Nlat, &
                                     self % emis_atlas % Npc))

!--------------------------------------------------------
! 2. Read the emissivity PCs
!--------------------------------------------------------

do i = 1, self % emis_atlas % nlon
  do j = 1, self % emis_atlas % nlat
    read (fileunit, '(12F10.6)', iostat = readstatus) self % emis_atlas % EmisPC(i,j,:)
  end do
end do

! Has there been an error in the read?
if (readstatus /= 0) then
  write(message,*) RoutineName,  &
       'Problem reading in EmisAtlas - please check the file'
  call abor1_ftn(message)
else
  write (*, '(A,I0,A)') 'Finished reading IR emissivity atlas with ', &
                         self % emis_atlas % Npc, ' principal components.'
end if

end subroutine ufo_rttovonedvarcheck_GetEmisAtlas

!------------------------------------------------------------------------------
!> Delete the PC emissivity object
!!
!! \details Heritage: Ops_SatRad_KillEmisEigenVec
!!
!! \author Met Office
!!
!! \date 04/08/2020: Created
!!
subroutine ufo_rttovonedvarcheck_DeletePCEmis(self)    ! inout

implicit none

! subroutine arguments:
class(ufo_rttovonedvarcheck_pcemis), intent(inout) :: self !< PC emissivity type

character(len=*), parameter :: routinename = "ufo_rttovonedvarcheck_DeletePCEmis"

self % emis_eigen % NChans = 0
self % emis_eigen % NumEV = 0
if (associated (self % emis_eigen % Channels)) deallocate (self % emis_eigen % Channels)
if (associated (self % emis_eigen % Mean)) deallocate (self % emis_eigen % Mean)
if (associated (self % emis_eigen % PCmin)) deallocate (self % emis_eigen % PCmin)
if (associated (self % emis_eigen % PCmax)) deallocate (self % emis_eigen % PCmax)
if (associated (self % emis_eigen % PCguess)) deallocate (self % emis_eigen % PCguess)
if (associated (self % emis_eigen % EV)) deallocate (self % emis_eigen % EV)
if (associated (self % emis_eigen % EV_Inverse)) deallocate (self % emis_eigen % EV_Inverse)

end subroutine ufo_rttovonedvarcheck_DeletePCEmis

!------------------------------------------------------------------------------
!> Print information about the PC Emissivity object
!!
!! \author Met Office
!!
!! \date 04/08/2020: Created
!!
subroutine ufo_rttovonedvarcheck_PrintPCEmis(self)

implicit none

! subroutine arguments:
class(ufo_rttovonedvarcheck_pcemis), intent(inout) :: self !< PC emissivity type

write(*,*) "Printing contents of PC emiss"
write(*,*) "emis_eigen % Channels = ",self % emis_eigen % Channels

end subroutine ufo_rttovonedvarcheck_PrintPCEmis

!------------------------------------------------------------------------------
!> Convert from spectral emissivity to principal component weights.
!! This is used to convert CAMEL emissivities to PC weights for the 1D-Var.
!!
!! \details Heritage: Ops_SatRad_EmissToPC.f90
!!
!! \author Met Office
!!
!! \date 04/08/2020: Created
!!
subroutine ufo_rttovonedvarcheck_EmisToPC (self,       &
                                           Channels,   &
                                           Emissivity, &
                                           PC)

implicit none

! Subroutine arguments
class(ufo_rttovonedvarcheck_pcemis), intent(inout) :: self !< PC emissivity type
integer, intent(in)          :: Channels(:)
real(kind_real), intent(in)  :: Emissivity(:)
real(kind_real), intent(out) :: PC(:)

! Local declarations:
character(len=*), parameter :: RoutineName = "ufo_rttovonedvarcheck_EmisToPC"
real(kind_real)             :: Temp_Emissivity(size(Emissivity))
integer                     :: npc
character(len=max_string)   :: message

npc = size(PC)

if (associated(self % emis_eigen % EV_Inverse)) then

  ! Convert from emissivity to sine transform
  Temp_Emissivity(:) = asin( 2.0_kind_real * Emissivity(:) - 1.0_kind_real )

  ! Subtract means from transformed emissivities
  Temp_Emissivity(:) = Temp_Emissivity(:) - self % emis_eigen % Mean(Channels(:))

  ! Calculate PC weights from emissivity spectrum
  PC(1:npc) = matmul(Temp_Emissivity(:), self % emis_eigen % EV_Inverse(Channels(:),1:npc))

else

  write(message, *) RoutineName,                             &
                 'Missing inverse eigenvector matrix - cannot convert emissivities to PCs'
  call abor1_ftn(message)

end if

end subroutine ufo_rttovonedvarcheck_EmisToPC

!-------------------------------------------------------------------------------
!> Transform from principal components to emissivity spectrum.
!!
!! \details Heritage: Ops_SatRad_PCToEmiss.f90
!!
!! \author Met Office
!!
!! \date 04/08/2020: Created
!!
subroutine ufo_rttovonedvarcheck_PCToEmis (self, &
                                NumChans,   &
                                Channels,   &
                                NumPC,      &
                                PC,         &
                                Emissivity)

implicit none

! Subroutine arguments:
class(ufo_rttovonedvarcheck_pcemis), intent(inout) :: self !< PC emissivity type
integer, intent(in)          :: NumChans
integer, intent(in)          :: Channels(NumChans)
integer, intent(in)          :: NumPC
real(kind_real), intent(in)  :: PC(NumPC)
real(kind_real), intent(out) :: Emissivity(NumChans)

! Local declarations:
character(len=*), parameter :: RoutineName = "ufo_rttovonedvarcheck_PCToEmis"
real(kind_real)             :: BigEOF(NumChans,NumChans)
real(kind_real)             :: BigPC(NumChans)

! Populate PC array with nchans elements
BigPC(:) = 0.0_kind_real
BigPC(1:NumPC) = PC(:)

! Populate EOF array with nchans elements
BigEOF(:,:) = 0.0_kind_real
BigEOF(1:NumPC,:) = self % emis_eigen % EV(1:NumPC,Channels)

! Calculate reconstructed emissivity spectrum
Emissivity = matmul(BigPC, BigEOF)

! Add means (these may have been subtracted off, otherwise they are zero)
Emissivity(:) = Emissivity(:) + self % emis_eigen % Mean(Channels)

! Convert from sine transform to physical emissivity
Emissivity(1:NumChans) = 0.5_kind_real * (sin(Emissivity(1:NumChans)) + 1.0_kind_real)

end subroutine ufo_rttovonedvarcheck_PCToEmis

!-------------------------------------------------------------------------------
!> Transform from principal components to emissivity spectrum.
!!
!! \details Heritage: Ops_SatRad_EmissKToPC.f90
!!
!! \author Met Office
!!
!! \date 04/08/2020: Created
!!
subroutine ufo_rttovonedvarcheck_EmisKToPC (self, &
                                 NumChans,     &
                                 Channels,     &
                                 NumPC,        &
                                 Emissivity,   &
                                 Emissivity_K, &
                                 PC_K)

implicit none

! Subroutine arguments
class(ufo_rttovonedvarcheck_pcemis), intent(inout) :: self !< PC emissivity type
integer, intent(in)          :: NumChans
integer, intent(in)          :: Channels(NumChans)
integer, intent(in)          :: NumPC
real(kind_real), intent(in)  :: Emissivity(NumChans)
real(kind_real), intent(in)  :: Emissivity_K(NumChans)
real(kind_real), intent(out) :: PC_K(NumChans,NumPC)

! Local declarations:
character(len=*), parameter :: RoutineName = "ufo_rttovonedvarcheck_EmisKToPC"
real                        :: JEMatrix_element
integer                     :: ichan

do ichan = 1, NumChans
  ! Calculate diagonal matrix elements of emissivity Jacobians
  ! Convert Emissivity_K to sine transform
  ! cos(asin(2x-1)) === sqrt(1-(1-2x)^2)
  JEMatrix_element = Emissivity_K(ichan) * 0.5_kind_real* &
                     sqrt (1.0_kind_real - (1.0_kind_real - 2.0_kind_real * Emissivity(ichan)) ** 2)

! EOF array === EmisEigenvec % EV for the appropriate channel selection
! N.B. Note 'manual' transposition of matrix here
  PC_K(ichan,1:NumPC) = self % emis_eigen % EV(1:NumPC,Channels(ichan)) * JEMatrix_element
end do

end subroutine ufo_rttovonedvarcheck_EmisKToPC

!-------------------------------------------------------------------------------

end module ufo_rttovonedvarcheck_pcemis_mod
