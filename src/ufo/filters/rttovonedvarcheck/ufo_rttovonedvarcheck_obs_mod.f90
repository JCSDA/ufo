! (C) copyright 2020 Met Office
!
! this software is licensed under the terms of the apache licence version 2.0
! which can be obtained at http://www.apache.org/licenses/license-2.0.

!> Fortran module which contains the observation data for a all the obs space

module ufo_rttovonedvarcheck_obs_mod

use datetime_mod, only: datetime
use fckit_exception_module, only: fckit_exception
use kinds
use iso_c_binding
use missing_values_mod
use obsspace_mod
use oops_variables_mod
use ufo_constants_mod, only: zero, one, Pa_to_hPa
use ufo_geovals_mod
use ufo_rttovonedvarcheck_constants_mod
use ufo_rttovonedvarcheck_pcemis_mod
use ufo_rttovonedvarcheck_profindex_mod
use ufo_rttovonedvarcheck_setup_mod
use ufo_vars_mod

implicit none
private

! Ob info type definition
type, public :: ufo_rttovonedvarcheck_obs

integer                      :: iloc
integer, allocatable         :: QCflags(:,:)    ! current qc flags needed for channel selection
real(kind_real), allocatable :: yobs(:,:)       ! observation value from obs files
real(kind_real), allocatable :: ybias(:,:)      ! observation bias from obs files
real(kind_real), allocatable :: lat(:)          ! observation latitude
real(kind_real), allocatable :: lon(:)          ! observation longitude
type(datetime), allocatable  :: date(:)         ! read in the date and time which is needed for ozone
real(kind_real), allocatable :: elevation(:)    ! observation elevation
real(kind_real), allocatable :: sat_zen(:)      ! observation satellite zenith angle
real(kind_real), allocatable :: sat_azi(:)      ! observation satellite azimuth angle
real(kind_real), allocatable :: sol_zen(:)      ! observation solar zenith angle
real(kind_real), allocatable :: sol_azi(:)      ! observation solar azimuth angle
real(kind_real), allocatable :: cloudtopp(:)    !< cloud top pressure (used in if cloudy retrieval used)
real(kind_real), allocatable :: cloudfrac(:)    !< cloud fraction (used in if cloudy retrieval used)
real(kind_real), allocatable :: cloudtopp_error(:)  !< cloud top pressure error calculated from retrieved profile
integer, allocatable         :: surface_type(:) ! surface type
integer, allocatable         :: satellite_identifier(:) ! WMO_ID
integer, allocatable         :: niter(:)        ! number of iterations
integer, allocatable         :: channels(:)     ! channel numbers
real(kind_real), allocatable :: final_cost(:)   ! final cost at solution
real(kind_real), allocatable :: LWP(:)          ! liquid water path from final iteration
real(kind_real), allocatable :: IWP(:)          ! liquid water path from final iteration
real(kind_real), allocatable :: clw(:,:)        ! cloud liquid water profile from final iteration
real(kind_real), allocatable :: skinTemperature(:) ! skin temperature for first iteration
real(kind_real), allocatable :: transmittance(:,:) ! surface to space transmittance for each channel
real(kind_real), allocatable :: output_profile(:,:) ! output profile
real(kind_real), allocatable :: output_BT(:,:)   ! output brightness temperature
real(kind_real), allocatable :: recalc_BT(:,:)   ! recalculate BT using retrieved variables for surface
real(kind_real), allocatable :: background_BT(:,:)   ! 1st iteration brightness temperature
logical                      :: Store1DVarLWP   ! flag to output the LWP if the profile converges
logical                      :: Store1DVarIWP   ! flag to output the IWP if the profile converges
logical                      :: Store1DVarCLW   ! flag to output the clw if the profile converges
logical                      :: Store1DVarTransmittance !  flag to output the surfacespace transmittance if profile convergences
real(kind_real), allocatable :: emiss(:,:)      ! initial surface emissivity
logical, allocatable         :: calc_emiss(:)   ! flag to request RTTOV calculate first guess emissivity
real(kind_real), allocatable :: mwemisserr(:,:) ! surface emissivity error from atlas
type(ufo_rttovonedvarcheck_pcemis), pointer :: pcemiss_object ! Infrared principal components object
real(kind_real), allocatable :: pcemiss(:,:)    ! principal component emissivity array
logical, allocatable         :: output_to_db(:)   ! flag to output data for this profile

contains
  procedure :: setup  => ufo_rttovonedvarcheck_obs_setup
  procedure :: delete => ufo_rttovonedvarcheck_obs_delete
  procedure :: output => ufo_rttovonedvarcheck_obs_output

end type

interface put_1d_indb
  module procedure put_1dint_indb
  module procedure put_1dfloat_indb
end interface

contains

!-------------------------------------------------------------------------------
!> Initialize observation object
!!
!! \author Met Office
!!
!! \date 09/06/2020: Created
!!
subroutine ufo_rttovonedvarcheck_obs_setup(self,       & ! out
                                           config,     & ! in
                                           prof_index, & ! in
                                           geovals,    & ! in
                                           vars)         ! in

implicit none

! subroutine arguments:
class(ufo_rttovonedvarcheck_obs), intent(out) :: self !< observation metadata type
type(ufo_rttovonedvarcheck), intent(in)       :: config !< observation metadata type
type(ufo_rttovonedvarcheck_profindex), intent(in) :: prof_index !< index to elements in the profile
type(ufo_geovals), intent(in)                 :: geovals  !< model data at obs location
type(oops_variables), intent(in)              :: vars     !< channels for 1D-Var

character(len=*), parameter :: routinename = "ufo_rttovonedvarcheck_obs_init"
real(kind_real)             :: missing_real
integer                     :: missing_int
integer                     :: jvar, iloc , jobs   !< counters
character(len=max_string)   :: var
character(len=max_string)   :: varname
character(len=200)          :: message
logical                     :: variable_present = .false.
logical                     :: model_surface_present = .false.
type(ufo_geoval), pointer   :: geoval
integer                     :: numpc

missing_real = missing_value(missing_real)
missing_int = missing_value(missing_int)
self % iloc = obsspace_get_nlocs(config % obsdb)

! allocate arrays
allocate(self % yobs(config % nchans, self % iloc))
allocate(self % ybias(config % nchans, self % iloc))
allocate(self % QCflags(config % nchans, self % iloc))
allocate(self % lat(self % iloc))
allocate(self % lon(self % iloc))
allocate(self % date(self % iloc))
allocate(self % elevation(self % iloc))
allocate(self % sat_zen(self % iloc))
allocate(self % sat_azi(self % iloc))
allocate(self % sol_zen(self % iloc))
allocate(self % sol_azi(self % iloc))
allocate(self % surface_type(self % iloc))
allocate(self % satellite_identifier(self % iloc))
allocate(self % niter(self % iloc))
allocate(self % channels(config % nchans))
allocate(self % final_cost(self % iloc))
allocate(self % LWP(self % iloc))
allocate(self % IWP(self % iloc))
allocate(self % skinTemperature(self % iloc))
if (config % Store1DVarCLW) allocate(self % CLW(config % nlevels, self % iloc))
if (config % Store1DVarTransmittance) allocate(self % transmittance(config % nchans, self % iloc))
allocate(self % output_profile(prof_index % nprofelements, self % iloc))
allocate(self % output_BT(config % nchans, self % iloc))
if (config % RecalculateBT) allocate(self % recalc_BT(config % nchans, self % iloc))
allocate(self % background_BT(config % nchans, self % iloc))
allocate(self % emiss(config % nchans, self % iloc))
allocate(self % calc_emiss(self % iloc))
allocate(self % mwemisserr(config % nchans, self % iloc))
allocate(self % output_to_db(self % iloc))

! initialize arrays
self % yobs(:,:) = missing_real
self % ybias(:,:) = zero
self % QCflags(:,:) = 0
self % lat(:) = missing_real
self % lon(:) = missing_real
self % elevation(:) = missing_real
self % sat_zen(:) = missing_real
self % sat_azi(:) = missing_real
self % sol_zen(:) = missing_real
self % sol_azi(:) = missing_real
self % surface_type(:) = missing_int
self % satellite_identifier(:) = missing_int
self % niter(:) = 0
self % channels(:) = 0
self % final_cost(:) = missing_real
self % LWP(:) = missing_real
self % IWP(:) = missing_real
self % skinTemperature(:) = missing_real
if (allocated(self % CLW)) self % CLW(:,:) = missing_real
if (allocated(self % transmittance)) self % transmittance(:,:) = missing_real
self % emiss(:,:) = missing_real
self % mwemisserr(:,:) = missing_real
self % output_profile(:,:) = missing_real
self % output_BT(:,:) = missing_real
if (allocated(self % recalc_BT)) self % recalc_BT(:,:) = missing_real
self % background_BT(:,:) = missing_real
self % calc_emiss(:) = .true.
self % Store1DVarLWP = config % Store1DVarLWP
self % Store1DVarIWP = config % Store1DVarIWP
self % Store1DVarCLW = config % Store1DVarCLW
self % Store1DVarTransmittance = config % Store1DVarTransmittance
self % pcemiss_object => null()
self % output_to_db(:) = .false.

! read in observations and associated errors / biases for full ObsSpace
do jvar = 1, config % nchans

  var = vars % variable(jvar)
  call obsspace_get_db(config % obsdb, "FortranQC", trim(var), self % QCflags(jvar,:))
  call obsspace_get_db(config % obsdb, "ObsValue",  trim(var), self % yobs(jvar,:))
  call obsspace_get_db(config % obsdb, "ObsBias",   trim(var), self % ybias(jvar,:))

end do

! The obs bias may contain missing values, and the intention is for those
! entries to be zero. Go through the bias data and replace the missing values
! with zero.
do jobs = 1, self % iloc
  do jvar = 1, config % nchans
    if (self % ybias(jvar, jobs) == missing_real) then
      self % ybias(jvar, jobs) = zero
    endif
  enddo
enddo

! Subtract bias from the observations (apply bias correction)
self % yobs = self % yobs - self % ybias

! Read in prerequisites
call obsspace_get_db(config % obsdb, "MetaData", "latitude", self % lat(:))
call obsspace_get_db(config % obsdb, "MetaData", "longitude", self % lon(:))
call obsspace_get_db(config % obsdb, "MetaData", "dateTime", self % date(:))
call obsspace_get_db(config % obsdb, "MetaData", "sensorZenithAngle", self % sat_zen(:))

! Read in optional angles
variable_present = obsspace_has(config % obsdb, "MetaData", "sensorAzimuthAngle")
if (variable_present) then
  call obsspace_get_db(config % obsdb, "MetaData", "sensorAzimuthAngle", self % sat_azi(:))
end if

variable_present = obsspace_has(config % obsdb, "MetaData", "solarZenithAngle")
if (variable_present) then
  call obsspace_get_db(config % obsdb, "MetaData", "solarZenithAngle", self % sol_zen(:))
end if

variable_present = obsspace_has(config % obsdb, "MetaData", "solarAzimuthAngle")
if (variable_present) then
  call obsspace_get_db(config % obsdb, "MetaData", "solarAzimuthAngle", self % sol_azi(:))
end if

! Read in initial cloud top pressure and cloud fraction if doing cloud retrieval
if (config % cloud_retrieval) then
  allocate(self % cloudtopp(self % iloc))
  allocate(self % cloudfrac(self % iloc))
  allocate(self % cloudtopp_error(self % iloc))
  self % cloudtopp(:) = 850.0_kind_real
  self % cloudfrac(:) = zero
  self % cloudtopp_error(:) = missing_real

  variable_present = obsspace_has(config % obsdb, "MetaData", "pressureAtTopOfCloud")
  if (variable_present) then
    call obsspace_get_db(config % obsdb, "MetaData", "pressureAtTopOfCloud", self % cloudtopp(:))
    where (self % cloudtopp /= missing_real)
      self % cloudtopp = self % cloudtopp * Pa_to_hPa
    end where
  end if

  variable_present = obsspace_has(config % obsdb, "MetaData", "cloudAmount")
  if (variable_present) then
    call obsspace_get_db(config % obsdb, "MetaData", "cloudAmount", self % cloudfrac(:))
  end if

  where(self % cloudfrac < zero .or. self % cloudfrac > one)
    self % cloudtopp = 850.0_kind_real
    self % cloudfrac = zero
  end where

end if

! Read in elevation for all obs
if (obsspace_has(config % obsdb, "MetaData", "heightOfSurface")) then
  call obsspace_get_db(config % obsdb, "MetaData", "heightOfSurface", self % elevation(:))
else if (ufo_vars_getindex(geovals % variables, 'surface_altitude') > 0) then
  call ufo_geovals_get_var(geovals, 'surface_altitude', geoval)
  self % elevation(:) = geoval%vals(1, :)
else
  self % elevation(:) = zero
endif

! Copy channels from config
self % channels(:) = config % channels(:)

! Read in surface type from ObsSpace or model data (deprecated)
if (obsspace_has(config % obsdb, "MetaData", "surfaceQualifier")) then
  call obsspace_get_db(config % obsdb, "MetaData", "surfaceQualifier", self % surface_type(:))
else
  call ufo_geovals_get_var(geovals, "surface_type", geoval)
  self % surface_type(:) = geoval%vals(1, :)
endif

! Read in satellite identifier
if (obsspace_has(config % obsdb, "MetaData", "satelliteIdentifier")) then
  call obsspace_get_db(config % obsdb, "MetaData", "satelliteIdentifier", self % satellite_identifier(:))
endif

! Read in skin temperature from the obs space
if (config % skinTemperatureFromObsSpace) then
  if (obsspace_has(config % obsdb, "MetaData", "skinTemperature")) then
    call obsspace_get_db(config % obsdb, "MetaData", "skinTemperature", self % skinTemperature(:))
  else
     write(message,*) "You have requested the initial skin temperature from the ObsSpace ", &
                      "but the array MetaData/skinTemperature is not present ", &
                      "=> aborting."
    call fckit_exception % throw(message)
  endif
end if

! Setup surface emissivity
! default self % emiss = zero
! default self % calc_emiss = true
select case (trim(config % EmissivityType))

  case("rttovtocalculate")
    write(*,*) "RTTOV will calculate emissivity first guess for sea, seaice and land"

  case("fixed")
    do iloc = 1, self % iloc
      self % calc_emiss(iloc) = .false.
      if(self % surface_type(iloc) == RTSea) then
        self % emiss(:,iloc) = config % EmissSeaDefault
      else if(self % surface_type(iloc) == RTLand) then
        self % emiss(:,iloc) = config % EmissLandDefault
      else if(self % surface_type(iloc) == RTIce) then
        self % emiss(:,iloc) = config % EmissSeaIceDefault
      end if
      if (self % emiss(1,iloc) == zero) self % calc_emiss(iloc) = .true.
    end do

  case("readfromdb")
    call ufo_rttovonedvarcheck_obs_ReadFromDB(self, config, vars, .false.)

  case("readfromdbwitherror")
    call ufo_rttovonedvarcheck_obs_ReadFromDB(self, config, vars, .true.)

  case("principalcomponent")
    numpc = prof_index % emisspc(2) - prof_index % emisspc(1) + 1
    call ufo_rttovonedvarcheck_obs_InitIRPCEmiss(self, config, numpc)

  case default
    write(message,*) "Land emissivity type not correctly defined it should be either ", &
                     "rttovtocalculate, fixed, readfromdb, readfromdbwitherror, or ", &
                     "princialcomponent"
    call abor1_ftn(message)

end select

end subroutine ufo_rttovonedvarcheck_obs_setup

!------------------------------------------------------------------------------
!> Delete the observation object
!!
!! \author Met Office
!!
!! \date 09/06/2020: Created
!!
subroutine ufo_rttovonedvarcheck_obs_delete(self)    ! inout

implicit none

! subroutine arguments:
class(ufo_rttovonedvarcheck_obs), intent(inout) :: self !< observation metadata type

character(len=*), parameter :: routinename = "ufo_rttovonedvarcheck_obs_delete"

! deallocate arrays
if (allocated(self % QCflags))        deallocate(self % QCflags)
if (allocated(self % yobs))           deallocate(self % yobs)
if (allocated(self % ybias))          deallocate(self % ybias)
if (allocated(self % lat))            deallocate(self % lat)
if (allocated(self % lon))            deallocate(self % lon)
if (allocated(self % date))           deallocate(self % date)
if (allocated(self % elevation))      deallocate(self % elevation)
if (allocated(self % sat_zen))        deallocate(self % sat_zen)
if (allocated(self % sat_azi))        deallocate(self % sat_azi)
if (allocated(self % sol_zen))        deallocate(self % sol_zen)
if (allocated(self % sol_azi))        deallocate(self % sol_azi)
if (allocated(self % cloudtopp))      deallocate(self % cloudtopp)
if (allocated(self % cloudfrac))      deallocate(self % cloudfrac)
if (allocated(self % cloudtopp_error)) deallocate(self % cloudtopp_error)
if (allocated(self % surface_type))   deallocate(self % surface_type)
if (allocated(self % satellite_identifier)) deallocate(self % satellite_identifier)
if (allocated(self % niter))          deallocate(self % niter)
if (allocated(self % final_cost))     deallocate(self % final_cost)
if (allocated(self % LWP))            deallocate(self % LWP)
if (allocated(self % IWP))            deallocate(self % IWP)
if (allocated(self % skinTemperature)) deallocate(self % skinTemperature)
if (allocated(self % CLW))            deallocate(self % CLW)
if (allocated(self % transmittance))  deallocate(self % transmittance)
if (allocated(self % emiss))          deallocate(self % emiss)
if (allocated(self % mwemisserr))     deallocate(self % mwemisserr)
if (allocated(self % output_profile)) deallocate(self % output_profile)
if (allocated(self % output_BT))      deallocate(self % output_BT)
if (allocated(self % recalc_BT))      deallocate(self % recalc_BT)
if (allocated(self % background_BT))  deallocate(self % background_BT)
if (allocated(self % calc_emiss))     deallocate(self % calc_emiss)
if (associated(self % pcemiss_object)) then
  call self % pcemiss_object % delete()
  self % pcemiss_object => null()
end if
if (allocated(self % pcemiss))        deallocate(self % pcemiss)

end subroutine ufo_rttovonedvarcheck_obs_delete

!------------------------------------------------------------------------------
!> Initialize the microwave emissivity array
!!
!! \details Heritage: Ops_SatRad_InitEmissivity.f90 - the MW part only
!!
!! \author Met Office
!!
!! \date 06/08/2020: Created
!!
subroutine ufo_rttovonedvarcheck_obs_ReadFromDB(self, config, vars, readerror)

implicit none

! subroutine arguments:
type(ufo_rttovonedvarcheck_obs), intent(inout) :: self !< observation metadata type
type(ufo_rttovonedvarcheck), intent(in) :: config !< main rttovonedvarcheck type
type(oops_variables), intent(in) :: vars     !< channels for 1D-Var
logical, intent(in) :: readerror

integer :: jvar, iloc
character(len=max_string)   :: var

! Read in the initial values from the db
do jvar = 1, size(self % channels)
  ! Read in from the db
  write(var,"(A11,I0)") "emissivity_", self % channels(jvar)
  call obsspace_get_db(config % obsdb, trim(config % EmissGroupInObsSpace), trim(var), self % emiss(jvar,:))
  if (readerror) then
    write(var,"(A16,I0)") "emissivityError_", self % channels(jvar)
    call obsspace_get_db(config % obsdb, trim(config % EmissGroupInObsSpace), trim(var), self % mwemisserr(jvar,:))
  end if

  ! Set calc emiss to true if emissivity is zero for any channel
  do iloc = 1, self % iloc
    self % calc_emiss(iloc) = .false.
    if(self % emiss(jvar,iloc) == zero) self % calc_emiss(iloc) = .true.
  end do
end do

end subroutine ufo_rttovonedvarcheck_obs_ReadFromDB

!------------------------------------------------------------------------------
!> Initialize the infrared emissivity array
!!
!! \details Heritage: Ops_SatRad_InitEmissivity.f90 - the IR part only
!!
!! \author Met Office
!!
!! \date 06/08/2020: Created
!!
subroutine ufo_rttovonedvarcheck_obs_InitIRPCEmiss(self, config, nemisspc)

implicit none

! subroutine arguments:
class(ufo_rttovonedvarcheck_obs), intent(inout) :: self !< observation metadata type
type(ufo_rttovonedvarcheck) :: config  !< main object containing configuration
integer, intent(in) :: nemisspc

! local variables
integer :: i
integer :: emis_x
integer :: emis_y

! Allocate and setup defaults - get RTTOV to calculate
self % emiss(:,:) = zero
self % calc_emiss(:) = .true.
allocate(self % pcemiss(nemisspc, self % iloc))
allocate(self % pcemiss_object)

! Setup pc emissivity object and initialise atlas (if needed)
if (len(trim(config % EmisAtlas)) > 0) then
  call self % pcemiss_object % setup(config % EmisEigVecPath, config % EmisAtlas)
else
  call self % pcemiss_object % setup(config % EmisEigVecPath)
end if

!-------------------------
! 1.2 Principal components
!-------------------------

! Initialise emissivity using principal components
if (self % pcemiss_object % initialised) then

  ! Don't do this if the RTTOV CAMEL atlas is being used.
  ! Defaults above will let the RTTOV camel atlas be used.
  !if (allocated (pcemiss_object % emis_eigen % PCGuess) .and. (.not. RTTOV_UseAtlas)) then
  if (allocated (self % pcemiss_object % emis_eigen % PCGuess)) then

    do i = 1, self % iloc

      ! Skip obs that may have invalid geographical coordinates
      !IF (BTEST (Obs % QCflags(i), QC_ModelDomain)) CYCLE

      ! If there's an atlas available, try to use it.
      if (allocated (self % pcemiss_object % emis_atlas % EmisPC)) then

        ! Find the nearest lat/lon
        emis_y = nint((self % lat(i) + 90.0_kind_real) / &
                   self % pcemiss_object % emis_atlas % gridstep + 1)
        emis_x = nint((self % lon(i) + 180.0_kind_real) / &
                   self % pcemiss_object % emis_atlas % gridstep + 1)
        if (emis_x > self % pcemiss_object % emis_atlas % nlon) then
          emis_x = emis_x - self % pcemiss_object % emis_atlas % nlon
        end if

        ! If the atlas is valid at this point, then use it,
        ! otherwise use PCGuess. NB: missing or sea points are
        ! flagged as -9.99 in the atlas.
        if (any (self % pcemiss_object % emis_atlas % EmisPC(emis_x,emis_y,:) > -9.99_kind_real)) then
          self % pcemiss(:,i) = self % pcemiss_object % emis_atlas % EmisPC(emis_x,emis_y,1:nemisspc)
        else
          self % pcemiss(:,i) = self % pcemiss_object % emis_eigen % PCGuess(1:nemisspc)
          !! Flag invalid atlas points over land as bad surface
          !IF (Obs % surface(i) == RTland) THEN
          !  Obs % QCflags(i) = IBSET (Obs % QCflags(i), QC_BadSurface)
          !END IF
        end if

      else ! If no atlas present, use PCGuess.
        self % pcemiss(:,i) = self % pcemiss_object % emis_eigen % PCGuess(1:nemisspc)
      end if

      ! If over sea make sure that emissivity is done by rttov
      if (self % surface_type(i) == RTSea) then
        self % calc_emiss(i) = .true.
        self % emiss(:, i) = zero
        self % pcemiss(:, i) = zero
      else
        self % calc_emiss(i) = .false.
        call self % pcemiss_object % pctoemis(size(self % channels), self % channels, nemisspc, &
                                         self % pcemiss(:, i), self % emiss(:, i))
      end if

    end do

  else
    call abor1_ftn("If emisspc is requested then EmisEigenVec file must be provided: aborting")
  end if

end if

end subroutine ufo_rttovonedvarcheck_obs_InitIRPCEmiss

!------------------------------------------------------------------------------
!> Store the 1D-Var analysis variables in obsspace for future assessment
!!
!! \details Heritage: Ops_SatRad_SetOutput_RTTOV12
!!
!! \author Met Office
!!
!! \date 02/09/2020: Created
!!
subroutine ufo_rttovonedvarcheck_obs_output(self, obsdb, prof_index, vars, nchans)

implicit none

! subroutine arguments:
class(ufo_rttovonedvarcheck_obs), intent(inout) :: self !< observation metadata type
type(c_ptr), value, intent(in)        :: obsdb          !< pointer to the observation space
type(ufo_rttovonedvarcheck_profindex), intent(in) :: prof_index !< index to elements in the profile
type(oops_variables), intent(in)      :: vars           !< channels for 1D-Var
integer, intent(in)                   :: nchans         ! number of channels in the obsspace

! local variables
integer :: jvar ! counter
integer :: nobs ! number of observations to be written to database
character(len=max_string)    :: var
real(kind_real), allocatable :: surface_pressure(:), ctp(:)
real(kind_real) :: missing_real

missing_real = missing_value(missing_real)

! Put QC flags and retrieved BT's back in database
do jvar = 1, nchans
  var = vars % variable(jvar)
  call put_1d_indb(self % output_to_db(:), obsdb, trim(var), "FortranQC", self % QCflags(jvar,:))
  call put_1d_indb(self % output_to_db(:), obsdb, trim(var), "OneDVar", self % output_BT(jvar,:))
  call put_1d_indb(self % output_to_db(:), obsdb, trim(var), "OneDVarBack", self % background_BT(jvar,:))
  if (allocated(self % recalc_BT)) then
    call put_1d_indb(self % output_to_db(:), obsdb, trim(var), "OneDVarRecalc", self % recalc_BT(jvar,:))
  end if
  write(var,"(A11,I0)") "emissivity_", self % channels(jvar)
  call put_1d_indb(self % output_to_db(:), obsdb, trim(var), "OneDVar", self % emiss(jvar,:))
  if (self % Store1DVarTransmittance) then
    write(var,"(A14,I0)") "transmittance_",self % channels(jvar)
    call put_1d_indb(self % output_to_db(:), obsdb, trim(var), "OneDVar", self % transmittance(jvar,:))
  endif
end do

! Output Diagnostics
call put_1d_indb(self % output_to_db(:), obsdb, "finalCost", "OneDVar", self % final_cost(:))
nobs = size(self % final_cost(:))
call put_1d_indb(self % output_to_db(:), obsdb, "numberOfIterations", "OneDVar", self % niter(:))

if (self % Store1DVarLWP) then
  call put_1d_indb(self % output_to_db(:), obsdb, "liquidWaterPath", "OneDVar", self % LWP(:))
end if

if (self % Store1DVarIWP) then
  call put_1d_indb(self % output_to_db(:), obsdb, "iceWaterPath", "OneDVar", self % IWP(:))
end if

!--
! Output Retrieved profiles into ObsSpace
!--
! 1) Temperature levels
! 2) Humidity levels: convert q to rh (%)
! 3) Ozone levels - no planned implementation
! 4b) Cloud Ice Water if rttovscat is activated from qtotal

!--
! 4) clw from qtotal
!--

if (self % Store1DVarCLW) then
  do jvar = 1, size(self % CLW,1)
    write(var,"(A,I0)") "lev",jvar
    call put_1d_indb(self % output_to_db(:), obsdb, var, "OneDVar/liquidWaterContent", self % CLW(jvar,:))
  end do
end if

!--
! 5) Surface pressure
!--
if (prof_index % pstar > 0) THEN
  allocate(surface_pressure(nobs))
  surface_pressure(:) = self % output_profile(prof_index % pstar, :)
  where (surface_pressure /= missing_real)
    surface_pressure = surface_pressure / Pa_to_hPa ! hPa to Pa
  end where
  call put_1d_indb(self % output_to_db(:), obsdb, "surfacePressure", "OneDVar", &
                   surface_pressure(:))

  deallocate(surface_pressure)
end if

!--
! 6) Surface temperature
!--
if (prof_index % t2 > 0) THEN
  call put_1d_indb(self % output_to_db(:), obsdb, "surfaceTemperature", "OneDVar", &
                   self % output_profile(prof_index % t2, :))
end if

!--
! 7) Surface humidity
!--

!--
! 8) Surface Windspeed
!--
! Windspeed retrieval is directionless, i.e., there are no separate u and v
! components.
if (prof_index % windspeed > 0) then
  call put_1d_indb(self % output_to_db(:), obsdb, "windSpeed", "OneDVar", &
                   self % output_profile(prof_index % windspeed, :))
end if

!--
! 9) Skin temperature
!--
if (prof_index % tstar > 0) then
  call put_1d_indb(self % output_to_db(:), obsdb, "skinTemperature", "OneDVar", &
                   self % output_profile(prof_index % tstar, :))
end if

!--
! 10) Total ozone - no planned implementation

!--
! 11) Cloud top pressure
! In the 1D-Var its hPa but needs to be Pa for the ObsSpace
!--
if (prof_index % cloudtopp > 0) then
  allocate(ctp(nobs))
  ctp(:) = self % output_profile(prof_index % cloudtopp, :)
  where (ctp /= missing_real)
    ctp = ctp / Pa_to_hPa ! hPa to Pa
  end where
  call put_1d_indb(self % output_to_db(:), obsdb, "pressureAtTopOfCloud", "OneDVar", ctp)
  ! Add the error estimate
  ctp(:) = self % cloudtopp_error(:)
  where (ctp /= missing_real)
    ctp = ctp / Pa_to_hPa ! hPa to Pa
  end where
  call put_1d_indb(self % output_to_db(:), obsdb, "pressureAtTopOfCloudError", "OneDVar", ctp)
  deallocate(ctp)
end if

! 12) Cloud fraction
if (prof_index % cloudfrac > 0) then
  call put_1d_indb(self % output_to_db(:), obsdb, "cloudAmount", "OneDVar", &
                   self % output_profile(prof_index % cloudfrac, :))
end if

! 14) IWP only meaningful is mwscattswitch is activated which means model levels too....
! 16/17) QC related not being ported.
! 18) Cloud type
! 19) Channel information
! 20) RTTOVSCATT cloud profiles
! 21) IR cloud profiles

end subroutine ufo_rttovonedvarcheck_obs_output

!-------------------------------------------------------------------------------

subroutine put_1dfloat_indb(apply, obsdb, variable, group, outputdata)
implicit none
logical, intent(in)             :: apply(:)  !< apply
type(c_ptr), value, intent(in)  :: obsdb !< pointer to the observation space
character(len=*), intent(in)    :: variable
character(len=*), intent(in)    :: group
real(kind_real), intent(in)     :: outputdata(:)

real(kind_real), allocatable :: tmp(:)
real(kind_real) :: missing_real
logical :: array_present
integer :: iprof, nobs

missing_real = missing_value(missing_real)
nobs = size(outputdata)

allocate(tmp(nobs))
tmp(:) = missing_real

! Get array from db if present
array_present = obsspace_has(obsdb, group, variable)
if (array_present) then
  call obsspace_get_db(obsdb, group, variable, tmp(:))
end if

! Update the tmp array with new data
do iprof = 1, nobs
  if (apply(iprof)) then
    tmp(iprof) = outputdata(iprof)
  end if
end do

! Write data to db
call obsspace_put_db(obsdb, group, variable, tmp)

deallocate(tmp)

end subroutine put_1dfloat_indb

!-------------------------------------------------------------------------------

subroutine put_1dint_indb(apply, obsdb, variable, group, outputdata)
implicit none
logical, intent(in)             :: apply(:)  !< apply
type(c_ptr), value, intent(in)  :: obsdb     !< pointer to the observation space
character(len=*), intent(in)    :: variable
character(len=*), intent(in)    :: group
integer, intent(in)             :: outputdata(:)

integer, allocatable :: tmp(:)
integer :: missing_int
logical :: array_present
integer :: iprof, nobs

missing_int = missing_value(missing_int)
nobs = size(outputdata)

allocate(tmp(nobs))
tmp(:) = missing_int

! Get array from db if present
array_present = obsspace_has(obsdb, group, variable)
if (array_present) then
  call obsspace_get_db(obsdb, group, variable, tmp(:))
end if

! Update the tmp array with new data
do iprof = 1, nobs
  if (apply(iprof)) then
    tmp(iprof) = outputdata(iprof)
  end if
end do

! Write data to db
call obsspace_put_db(obsdb, group, variable, tmp)

deallocate(tmp)

end subroutine put_1dint_indb

!-------------------------------------------------------------------------------

end module ufo_rttovonedvarcheck_obs_mod
