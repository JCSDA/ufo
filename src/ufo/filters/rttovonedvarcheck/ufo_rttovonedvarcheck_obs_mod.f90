! (C) copyright 2020 Met Office
!
! this software is licensed under the terms of the apache licence version 2.0
! which can be obtained at http://www.apache.org/licenses/license-2.0.

!> Fortran module which contains the observation data for a all the obs space

module ufo_rttovonedvarcheck_obs_mod

use kinds
use iso_c_binding
use missing_values_mod
use obsspace_mod
use oops_variables_mod
use ufo_geovals_mod
use ufo_rttovonedvarcheck_constants_mod
use ufo_rttovonedvarcheck_pcemis_mod
use ufo_rttovonedvarcheck_profindex_mod
use ufo_rttovonedvarcheck_utils_mod
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
real(kind_real), allocatable :: elevation(:)    ! observation elevation
real(kind_real), allocatable :: sat_zen(:)      ! observation satellite zenith angle
real(kind_real), allocatable :: sat_azi(:)      ! observation satellite azimuth angle
real(kind_real), allocatable :: sol_zen(:)      ! observation solar zenith angle
real(kind_real), allocatable :: sol_azi(:)      ! observation solar azimuth angle
integer, allocatable         :: surface_type(:) ! surface type
real(kind_real), allocatable :: final_cost(:)   ! final cost at solution
real(kind_real), allocatable :: emiss(:,:)      ! initial surface emissivity
real(kind_real), allocatable :: output_profile(:,:) ! output profile
real(kind_real), allocatable :: output_BT(:,:)  ! output brightness temperature
logical, allocatable         :: calc_emiss(:)   ! flag to request RTTOV calculate first guess emissivity

contains
  procedure :: setup  => ufo_rttovonedvarcheck_obs_setup
  procedure :: delete => ufo_rttovonedvarcheck_obs_delete
  procedure :: output => ufo_rttovonedvarcheck_obs_output

end type

contains

!-------------------------------------------------------------------------------
!> Initialize observation object
!!
!! \author Met Office
!!
!! \date 09/06/2020: Created
!!
subroutine ufo_rttovonedvarcheck_obs_setup(self,     & ! out
                                           config,   & ! in
                                           nprofelements, & ! in
                                           geovals,  & ! in
                                           vars,     & ! in
                                           ir_pcemis )

implicit none

! subroutine arguments:
class(ufo_rttovonedvarcheck_obs), intent(out) :: self !< observation metadata type
type(ufo_rttovonedvarcheck), intent(in)       :: config !< observation metadata type
integer, intent(in)                           :: nprofelements !< number of profile elements
type(ufo_geovals), intent(in)                 :: geovals  !< model data at obs location
type(oops_variables), intent(in)              :: vars     !< channels for 1D-Var
type(ufo_rttovonedvarcheck_pcemis)            :: ir_pcemis  !< Infrared principal components object

character(len=*), parameter :: routinename = "ufo_rttovonedvarcheck_obs_init"
real(kind_real)             :: missing
integer                     :: jvar    !< counters
character(len=max_string)   :: var
character(len=max_string)   :: varname
logical                     :: variable_present = .false.
logical                     :: model_surface_present = .false.
type(ufo_geoval), pointer   :: geoval

missing = missing_value(missing)
self % iloc = obsspace_get_nlocs(config % obsdb)

! allocate arrays
allocate(self % yobs(config % nchans, self % iloc))
allocate(self % ybias(config % nchans, self % iloc))
allocate(self % QCflags(config % nchans, self % iloc))
allocate(self % lat(self % iloc))
allocate(self % lon(self % iloc))
allocate(self % elevation(self % iloc))
allocate(self % sat_zen(self % iloc))
allocate(self % sat_azi(self % iloc))
allocate(self % sol_zen(self % iloc))
allocate(self % sol_azi(self % iloc))
allocate(self % surface_type(self % iloc))
allocate(self % final_cost(self % iloc))
allocate(self % emiss(config % nchans, self % iloc))
allocate(self % output_profile(nprofelements, self % iloc))
allocate(self % output_BT(config % nchans, self % iloc))
allocate(self % calc_emiss(self % iloc))

! initialize arrays
self % yobs(:,:) = missing
self % ybias(:,:) = 0.0
self % QCflags(:,:) = 0
self % lat(:) = missing
self % lon(:) = missing
self % elevation(:) = missing
self % sat_zen(:) = missing
self % sat_azi(:) = missing
self % sol_zen(:) = missing
self % sol_azi(:) = missing
self % surface_type = RTSea
self % final_cost = missing
self % emiss(:,:) = 0.0
self % output_profile(:,:) = missing
self % output_BT(:,:) = missing
self % calc_emiss(:) = .true.

! read in observations and associated errors / biases for full ObsSpace
do jvar = 1, config % nchans

  var = vars % variable(jvar)
  call obsspace_get_db(config % obsdb, "FortranQC", trim(var), self % QCflags(jvar,:))
  call obsspace_get_db(config % obsdb, "ObsValue",  trim(var), self % yobs(jvar,:))

  ! Optionally get the observation bias
  variable_present = obsspace_has(config % obsdb, "ObsBias", trim(var))
  if (variable_present) then
    call obsspace_get_db(config % obsdb, "ObsBias",   trim(var), self % ybias(jvar,:))
  end if

end do

if (.not. variable_present) write(*,*) "Using uncorrected brightness temperature"

! Add bias correction to the observations
self % yobs = self % yobs + self % ybias

! Read in prerequisites
call obsspace_get_db(config % obsdb, "MetaData", "latitude", self % lat(:))
call obsspace_get_db(config % obsdb, "MetaData", "longitude", self % lon(:))
call obsspace_get_db(config % obsdb, "MetaData", "sensor_zenith_angle", self % sat_zen(:))

! Read in optional angles
variable_present = obsspace_has(config % obsdb, "MetaData", "sensor_azimuth_angle")
if (variable_present) then
  call obsspace_get_db(config % obsdb, "MetaData", "sensor_azimuth_angle", self % sat_azi(:))
end if

variable_present = obsspace_has(config % obsdb, "MetaData", "solar_zenith_angle")
if (variable_present) then
  call obsspace_get_db(config % obsdb, "MetaData", "solar_zenith_angle", self % sol_zen(:))
end if

variable_present = obsspace_has(config % obsdb, "MetaData", "solar_azimuth_angle")
if (variable_present) then
  call obsspace_get_db(config % obsdb, "MetaData", "solar_azimuth_angle", self % sol_azi(:))
end if

! Read in elevation for all obs
variable_present = obsspace_has(config % obsdb, "MetaData", "elevation")
if (variable_present) then
  call obsspace_get_db(config % obsdb, "MetaData", "elevation", self % elevation(:))
else
  model_surface_present = obsspace_has(config % obsdb, "MetaData", "model_surface")
  if (model_surface_present) then
    call obsspace_get_db(config % obsdb, "MetaData", "model_surface", self % elevation(:))
  else
    self % elevation(:) = 0.0
  end if
end if

! Read in surface type from model data
call ufo_geovals_get_var(geovals, "surface_type", geoval)
self % surface_type(:) = geoval%vals(1,:)

! Setup emissivity
if (config % pcemiss) then
  write(*,*) "PC emissivity being used"
  call ufo_rttovonedvarcheck_obs_InitIREmiss(self, config % nchans, ir_pcemis)
else
  write(*,*) "Conventional emissivity being used"
  call ufo_rttovonedvarcheck_obs_InitMWEmiss(self, config)
end if

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
if (allocated(self % QCflags))    deallocate(self % QCflags)
if (allocated(self % yobs))       deallocate(self % yobs)
if (allocated(self % ybias))      deallocate(self % ybias)
if (allocated(self % lat))        deallocate(self % lat)
if (allocated(self % lon))        deallocate(self % lon)
if (allocated(self % elevation))  deallocate(self % elevation)
if (allocated(self % sat_zen))    deallocate(self % sat_zen)
if (allocated(self % sat_azi))    deallocate(self % sat_azi)
if (allocated(self % sol_zen))    deallocate(self % sol_zen)
if (allocated(self % sol_azi))    deallocate(self % sol_azi)
if (allocated(self % surface_type)) deallocate(self % surface_type)
if (allocated(self % final_cost))   deallocate(self % final_cost)
if (allocated(self % emiss))        deallocate(self % emiss)
if (allocated(self % output_profile)) deallocate(self % output_profile)
if (allocated(self % output_BT))  deallocate(self % output_BT)
if (allocated(self % calc_emiss)) deallocate(self % calc_emiss)

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
subroutine ufo_rttovonedvarcheck_obs_InitMWEmiss(self, config)

implicit none

! subroutine arguments:
type(ufo_rttovonedvarcheck_obs), intent(inout) :: self !< observation metadata type
type(ufo_rttovonedvarcheck), intent(in) :: config !< main rttovonedvarcheck type

integer :: i

!-------------
! 2.1 Defaults
!-------------

do i = 1, self % iloc

  ! Only calculate in RTTOV over sea
  if (self % surface_type(i) == RTSea) then
    self % calc_emiss(i) = .true.
  else
    self % calc_emiss(i) = .false.
  end if

  ! The default emissivity for land is a very crude estimate - the same
  ! for all surface types and all frequencies. However, we do not use
  ! channels which see the surface over land where we rely on this default.
  self % emiss(:,i) = 0.0
  if (self % surface_type(i) == RTLand) then
    self % emiss(:,i) = config % EmissLandDefault
  else if (self % surface_type(i) == RTIce) then
    self % emiss(:,i) = config % EmissSeaIceDefault
  end if

end do

end subroutine

!------------------------------------------------------------------------------
!> Initialize the infrared emissivity array
!!
!! \details Heritage: Ops_SatRad_InitEmissivity.f90 - the IR part only
!!
!! \author Met Office
!!
!! \date 06/08/2020: Created
!!
subroutine ufo_rttovonedvarcheck_obs_InitIREmiss(self, nchans, ir_pcemis)

implicit none

! subroutine arguments:
class(ufo_rttovonedvarcheck_obs), intent(inout) :: self !< observation metadata type
integer, intent(in) :: nchans   !< total number of channels
type(ufo_rttovonedvarcheck_pcemis) :: ir_pcemis  !< Infrared principal components object

! local variables
real(kind_real), allocatable :: EmissPC(:,:)
integer :: i
integer :: nemisspc = 5
integer :: emis_x
integer :: emis_y

! Allocate and setup defaults - get RTTOV to calculate
self % emiss(:,:) = 0.0
self % calc_emiss(:) = .true.

allocate(EmissPC(self % iloc,nemisspc))

!-------------------------
! 1.2 Principal components
!-------------------------

! Initialise emissivity using principal components
if (ir_pcemis % initialised) then

  ! Don't do this if the RTTOV CAMEL atlas is being used.
  ! Defaults above will let the RTTOV camel atlas be used.
  !if (associated (ir_pcemiss % emis_eigen % PCGuess) .and. (.not. RTTOV_UseAtlas)) then
  if (associated (ir_pcemis % emis_eigen % PCGuess)) then

    do i = 1, self % iloc

      ! Skip obs that may have invalid geographical coordinates
      !IF (BTEST (Obs % QCflags(i), QC_ModelDomain)) CYCLE

      ! If there's an atlas available, try to use it.
      if (associated (ir_pcemis % emis_atlas % EmisPC)) then

        ! Find the nearest lat/lon
        emis_y = nint((self % lat(i) + 90.0) / ir_pcemis % emis_atlas % gridstep + 1)
        emis_x = nint((self % lon(i) + 180.0) / ir_pcemis % emis_atlas % gridstep + 1)
        if (emis_x > ir_pcemis % emis_atlas % nlon) then
          emis_x = emis_x - ir_pcemis % emis_atlas % nlon
        end if

        ! If the atlas is valid at this point, then use it,
        ! otherwise use PCGuess. NB: missing or sea points are
        ! flagged as -9.99 in the atlas.
        if (any (ir_pcemis % emis_atlas % EmisPC(emis_x,emis_y,:) > -9.99)) then
          EmissPC(i,:) = ir_pcemis % emis_atlas % EmisPC(emis_x,emis_y,1:nemisspc)
        else
          EmissPC(i,:) = ir_pcemis % emis_eigen % PCGuess(1:nemisspc)
          !! Flag invalid atlas points over land as bad surface
          !IF (Obs % surface(i) == RTland) THEN
          !  Obs % QCflags(i) = IBSET (Obs % QCflags(i), QC_BadSurface)
          !END IF
        END IF

      ELSE ! If no atlas present, use PCGuess.
        EmissPC(i,:) = ir_pcemis % emis_eigen % PCGuess(1:nemisspc)
      END IF

    end do

  end if

end if

if (allocated(EmissPC)) deallocate(EmissPC)

end subroutine

!------------------------------------------------------------------------------
!> Initialize the infrared emissivity array
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
real(kind_real), allocatable :: surface_pressure(:)
real(kind_real) :: missing

missing = missing_value(missing)

! Put QC flags and retrieved BT's back in database
do jvar = 1, nchans
  var = vars % variable(jvar)
  call obsspace_put_db(obsdb, "FortranQC", trim(var), self % QCflags(jvar,:))
  call obsspace_put_db(obsdb, "OneDVar", trim(var), self % output_BT(jvar,:))
end do

! Output final cost at solution
call obsspace_put_db(obsdb, "OneDVar", "FinalCost", self % final_cost(:))
nobs = size(self % final_cost(:))

!--
! Output Retrieved profiles into ObsSpace
!--
! 1) Temperature levels
! 2) Humidity levels: convert q to rh (%)
! 3) Ozone levels - no planned implementation
! 4) Cloud Liquid Water from qtotal
! 4b) Cloud Ice Water if rttovscat is activated from qtotal

!--
! 5) Surface pressure
!--
if (prof_index % pstar > 0) THEN
  allocate(surface_pressure(nobs))
  surface_pressure(:) = self % output_profile(prof_index % pstar, :)
  where (surface_pressure /= missing)
    surface_pressure = surface_pressure * 100 ! Pa to hPa
  end where
  call obsspace_put_db(obsdb, "OneDVar", trim(var_ps), surface_pressure)
  deallocate(surface_pressure)
end if

!--
! 6) Surface temperature
!--
if (prof_index % t2 > 0) THEN
  call obsspace_put_db(obsdb, "OneDVar", trim(var_sfc_t2m), &
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
  call obsspace_put_db(obsdb, "OneDVar", "surface_wind_speed", &
                       self % output_profile(prof_index % windspeed, :))
end if

!--
! 9) Skin temperature
!--
if (prof_index % tstar > 0) then
  call obsspace_put_db(obsdb, "OneDVar", trim(var_sfc_tskin), &
                       self % output_profile(prof_index % tstar, :))
end if

!--
! 10) Total ozone - no planned implementation
! 11) Cloud top pressure
! 12) Cloud fraction
! 13) LWP
! 14) IWP only meaningful is mwscattswitch is activated which means model levels too....
! 15) Microwave emissivity
! 16/17) QC related not being ported.
! 18) Cloud type
! 19) Channel information
! 20) RTTOVSCATT cloud profiles
! 21) IR cloud profiles

end subroutine

!-------------------------------------------------------------------------------

end module ufo_rttovonedvarcheck_obs_mod
