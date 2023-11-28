!-------------------------------------------------------------------------------
! (C) British Crown Copyright 2020 Met Office
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
!-------------------------------------------------------------------------------

!> Fortran module for gnssro bending angle Met Office forward operator

module ufo_gnssro_bendmetoffice_mod

use iso_c_binding
use kinds
use ufo_vars_mod
use ufo_geovals_mod
use ufo_geovals_mod_c, only: ufo_geovals_registry
use vert_interp_mod
use lag_interp_mod,    only: lag_interp_const, lag_interp_smthWeights
use obsspace_mod
use missing_values_mod
use ufo_gnssro_ukmo1d_utils_mod
use ufo_utils_refractivity_calculator, only: ufo_calculate_refractivity
use fckit_log_module,  only : fckit_log
use fckit_exception_module, only: fckit_exception

implicit none
public             :: ufo_gnssro_bendmetoffice
private

  !> Fortran derived type for gnssro trajectory
type :: ufo_gnssro_BendMetOffice
  logical :: vert_interp_ops
  logical :: pseudo_ops
  logical :: noSuperCheck
  real(kind_real) :: min_temp_grad
  integer, allocatable :: chanList(:)
  contains
    procedure :: setup     => ufo_gnssro_bendmetoffice_setup
    procedure :: simobs    => ufo_gnssro_bendmetoffice_simobs
end type ufo_gnssro_BendMetOffice

contains

! ------------------------------------------------------------------------------
! Get the optional settings for the forward model, and save them in the object
! so that they can be used in the code.
! ------------------------------------------------------------------------------
subroutine ufo_gnssro_bendmetoffice_setup(self, vert_interp_ops, pseudo_ops, min_temp_grad, chanList, noSuperCheck)

implicit none

class(ufo_gnssro_BendMetOffice), intent(inout) :: self
logical(c_bool), intent(in) :: vert_interp_ops
logical(c_bool), intent(in) :: pseudo_ops
real(c_float), intent(in) :: min_temp_grad
integer(c_int), intent(in) :: chanList(:)
logical(c_bool), intent(in) :: noSuperCheck

character(len=*), parameter  :: myname_ = "ufo_gnssro_bendmetoffice_setup"
integer, parameter           :: max_string = 800
character(max_string)        :: message                       ! General message for output
integer                      :: i                             ! Loop variable

self % vert_interp_ops = vert_interp_ops
self % pseudo_ops = pseudo_ops
self % min_temp_grad = min_temp_grad
allocate(self % chanList(1:size(chanList)))
self % chanList = chanList
self % noSuperCheck = noSuperCheck

write(message, *) myname_, ' Setting up Met Office GNSS-RO forward operator with'
call fckit_log%info(message)
write(message, *) 'vert_interp_ops =', self % vert_interp_ops
call fckit_log%info(message)
write(message, *) 'pseudo_ops =', self % pseudo_ops
call fckit_log%info(message)
write(message, *) 'min_temp_grad =', self % min_temp_grad
call fckit_log%info(message)
write(message, *) 'no super check =', self % noSuperCheck
call fckit_log%info(message)
write(message, '(A)') 'chanList = '
call fckit_log % debug(message)
do i = 1, SIZE(chanList), 100
  write(message, '(100I5)') chanList(i:min(i+99, size(chanList)))
  call fckit_log % debug(message)
end do

end subroutine ufo_gnssro_bendmetoffice_setup

! ------------------------------------------------------------------------------
! 1-dimensional GNSS-RO forward operator for the Met Office system
! ------------------------------------------------------------------------------
subroutine ufo_gnssro_bendmetoffice_simobs(self, geovals, obss, nlevels, nlocs, hofx, obs_diags)

  implicit none

  ! Arguments to this routine
  class(ufo_gnssro_BendMetOffice), intent(in)    :: self                   ! The object in which this operator is contained
  type(ufo_geovals),               intent(in)    :: geovals                ! The model values, interpolated to the observation locations
  integer(c_int),                  intent(in)    :: nlevels                ! Number of variables (levels) in the data
  integer(c_int),                  intent(in)    :: nlocs                  ! Number of observations (locations)
  real(kind_real),                 intent(inout) :: hofx(nlevels, nlocs)   ! The model forecast of the observations
  type(c_ptr), value,              intent(in)    :: obss                   ! The observations, and meta-data for those observations
  type(ufo_geovals),               intent(inout) :: obs_diags              ! Observations diagnostics

  character(len=*), parameter  :: myname_ = "ufo_gnssro_bendmetoffice_simobs"
  integer, parameter           :: max_string = 800

  character(max_string)        :: err_msg                       ! Error message for output
  character(max_string)        :: message                       ! General message for output
  integer                      :: nlocs_check                   ! Number of observations (profiles)
  integer                      :: nlevels_check                 ! Number of observations in a profile
  integer                      :: ilev                          ! Loop variable, level number
  integer                      :: iloc                          ! Loop variable, observation number
  type(ufo_geoval), pointer    :: q                             ! Model background values of specific humidity
  type(ufo_geoval), pointer    :: prs                           ! Model background values of air pressure
  type(ufo_geoval), pointer    :: theta_heights                 ! Model heights of levels containing specific humidity
  type(ufo_geoval), pointer    :: rho_heights                   ! Model heights of levels containing air pressure
  real(kind_real), allocatable :: obsLat(:)                     ! Latitude of the observation
  real(kind_real), allocatable :: obsLon(:)                     ! Longitude of the observation
  real(kind_real), allocatable :: impact_param(:)               ! Impact parameter of the observation
                                                                ! Note: impact_param can be represented as 2D data, but the array
                                                                ! is 1D as the interface to get_db is 1D.
  real(kind_real), allocatable :: radius_curv(:)                ! Earth's radius of curvature at the observation tangent point
  real(kind_real), allocatable :: undulation(:)                 ! Undulation - height of the geoid above the ellipsoid
  logical                      :: BAErr                         ! Was there an error in the calculation?
  integer                      :: iVar                          ! Loop variable, obs diagnostics variable number
  real(kind_real), allocatable :: refractivity(:)               ! Refractivity on various model levels
  real(kind_real), allocatable :: model_heights(:)              ! Geopotential heights that refractivity is calculated on
  real(kind_real)              :: calculated_hofx(nlevels)      ! Array to receive the calculated h(x) on levels
  real(kind_real), allocatable :: tobs(:)                       ! Virtual temperature at observation locations

  write(err_msg,*) "TRACE: ufo_gnssro_bendmetoffice_simobs: begin"
  call fckit_log%info(err_msg)

  ! If output to refractivity (and heights of the refractivity levels) is needed,
  ! then use nval as a way to check whether the array has been initialised (since
  ! it is called in a loop).
  DO iVar = 1, obs_diags % nvar
    IF (obs_diags % variables(ivar) == "atmosphericRefractivity_model" .OR. &
        obs_diags % variables(ivar) == "geopotentialHeight_model" .OR. &
        obs_diags % variables(ivar) == "virtualTemperature") THEN
      write(err_msg,*) "TRACE: ufo_gnssro_bendmetoffice_simobs: initialising obs_diags for " // &
        obs_diags % variables(ivar)
      call fckit_log%info(err_msg)
      obs_diags % geovals(iVar) % nval = 0
    END IF
  END DO

! check if nlocs is consistent in geovals and what was passed in
  if (geovals%nlocs /= nlocs) then
      write(err_msg,*) myname_, ' error: nlocs inconsistent between geovals and what was passed in!'
      call abor1_ftn(err_msg)
  endif

! get variables from geovals
  call ufo_geovals_get_var(geovals, var_q, q)               ! specific humidity
  call ufo_geovals_get_var(geovals, var_prsi, prs)          ! pressure
  call ufo_geovals_get_var(geovals, var_z, theta_heights)   ! Geopotential height of the normal model levels
  call ufo_geovals_get_var(geovals, var_zi, rho_heights)    ! Geopotential height of the pressure levels

  write(message, '(A,10I6)') 'Q: ', q%nval, q%nprofiles, shape(q%vals)
  call fckit_log%info(message)
  write(message, '(A,10I6)') 'Pressure: ', prs%nval, prs%nprofiles, shape(prs%vals)
  call fckit_log%info(message)

  nlocs_check = obsspace_get_nlocs(obss)
  ! nchans is used to define the number of vertical levels available in the data
  ! (i.e. the length of an observed profile)
  nlevels_check = obsspace_get_nchans(obss)
  ! If the data has been read in a 1D, then nchans may not have been set, so will
  ! default to zero, rather than one.  So reset this:
  if (nlevels_check == 0) nlevels_check = 1
  if (nlocs /= nlocs_check .OR. nlevels_check /= nlevels) then
    write(err_msg,'(2A,4I8)') myname_, ' error: nlocs or nlevels doesnt match', nlocs, nlocs_check, nlevels, nlevels_check
    call fckit_exception%throw(err_msg)
  end if

  if (nlevels > 1 .AND. nlevels /= size(self % chanList)) then
    write(err_msg,'(2A,4I8)') myname_, ' error: channel list must match nlevels', nlevels, size(self%chanList)
    call fckit_exception%throw(err_msg)
  end if

  if (prs%vals(1,1) .gt. prs%vals(prs%nval,1) ) then
    write(err_msg,'(a)') 'Geovals should be ordered top to bottom'
    call fckit_exception%throw(err_msg)
  end if

  allocate(obsLat(nlocs))
  allocate(obsLon(nlocs))
  allocate(impact_param(nlevels * nlocs))
  allocate(radius_curv(nlocs))
  allocate(undulation(nlocs))
  allocate(tobs(nlevels))

  call obsspace_get_db(obss, "MetaData", "longitude", obsLon)
  call obsspace_get_db(obss, "MetaData", "latitude", obsLat)
  if (nlevels > 1) then
    call obsspace_get_db(obss, "MetaData", "impactParameterRO", impact_param, self % chanList)
  else
    call obsspace_get_db(obss, "MetaData", "impactParameterRO", impact_param)
  end if
  call obsspace_get_db(obss, "MetaData", "earthRadiusCurvature", radius_curv)
  call obsspace_get_db(obss, "MetaData", "geoidUndulation", undulation)

  obs_loop: do iloc = 1, nlocs

    ! Note: The forward operator is coded for bottom-to-top geovals.  Therefore
    ! the arguments are flipped as they are passed.  This will not create a memory
    ! overhead since the geovals are processed one at a time.
    call Ops_GPSRO_ForwardModel(prs % nval, &
                                q % nval, &
                                rho_heights % vals(rho_heights%nval:1:-1, iloc), &
                                theta_heights % vals(theta_heights%nval:1:-1, iloc), &
                                prs % vals(prs%nval:1:-1, iloc), &
                                q % vals(q%nval:1:-1, iloc), &
                                self % pseudo_ops, &
                                self % vert_interp_ops, &
                                self % min_temp_grad, &
                                nlevels, &
                                impact_param(1+(iloc-1)*nlevels:iloc*nlevels), &
                                radius_curv(iloc), &
                                obsLat(iloc), &
                                undulation(iloc), &
                                calculated_hofx, &
                                BAErr, &
                                refractivity, &
                                model_heights, &
                                self % noSuperCheck, &
                                tobs)
    hofx(:, iloc) = calculated_hofx

    if (BAErr) then
      write(err_msg,*) "Error with observation processing ", iloc
      call fckit_log % info(err_msg)
    end if

    ! If output to refractivity is needed, then initialise things
    DO iVar = 1, obs_diags % nvar
        IF (obs_diags % variables(ivar) == "atmosphericRefractivity_model") THEN
            IF (iloc == 1) THEN
                obs_diags % geovals(iVar) % nval = SIZE(refractivity)
                ALLOCATE(obs_diags % geovals(iVar) % vals(SIZE(refractivity), obs_diags % nlocs))
            END IF

            IF (BAerr) THEN
                obs_diags % geovals(iVar) % vals(:,iloc) = missing_value(obs_diags % geovals(iVar) % vals(1,1))
            ELSE
                ! Flip the order of the calculated refractivity, since the geovals are oriented
                ! top-to-bottom, but the code works bottom-to-top
                obs_diags % geovals(iVar) % vals(:,iloc) = refractivity(SIZE(refractivity):1:-1)
            END IF
        END IF

        IF (obs_diags % variables(ivar) == "geopotentialHeight_model") THEN
            IF (iloc == 1) THEN
                obs_diags % geovals(iVar) % nval = SIZE(model_heights)
                ALLOCATE(obs_diags % geovals(iVar) % vals(SIZE(model_heights), obs_diags % nlocs))
            END IF

            IF (BAerr) THEN
                obs_diags % geovals(iVar) % vals(:,iloc) = missing_value(obs_diags % geovals(iVar) % vals(1,1))
            ELSE
                ! Flip the order of the calculated refractivity, since the geovals are oriented
                ! top-to-bottom, but the code works bottom-to-top
                obs_diags % geovals(iVar) % vals(:,iloc) = model_heights(SIZE(model_heights):1:-1)
            END IF
        END IF

        IF (obs_diags % variables(ivar) == "virtualTemperature") THEN
            IF (iloc == 1) THEN
                obs_diags % geovals(iVar) % nval = nlevels
                ALLOCATE(obs_diags % geovals(iVar) % vals(nlevels, obs_diags % nlocs))
            END IF

            IF (BAerr) THEN
                obs_diags % geovals(iVar) % vals(:,iloc) = missing_value(obs_diags % geovals(iVar) % vals(1,1))
            ELSE
                obs_diags % geovals(iVar) % vals(:,iloc) = tobs(:)
            END IF
        END IF
    END DO
  end do obs_loop

  deallocate(obsLat)
  deallocate(obsLon)
  deallocate(impact_param)
  deallocate(radius_curv)
  deallocate(undulation)

  write(err_msg,*) "TRACE: ufo_gnssro_bendmetoffice_simobs: completed"
  call fckit_log%info(err_msg)

end subroutine ufo_gnssro_bendmetoffice_simobs
! ------------------------------------------------------------------------------


SUBROUTINE Ops_GPSRO_ForwardModel(nlevp, &
                                  nlevq, &
                                  za, &
                                  zb, &
                                  pressure, &
                                  humidity, &
                                  GPSRO_pseudo_ops, &
                                  GPSRO_vert_interp_ops, &
                                  GPSRO_min_temp_grad, &
                                  nobs, &
                                  zobs, &
                                  RO_Rad_Curv, &
                                  Latitude, &
                                  RO_geoid_und, &
                                  ycalc, &
                                  BAErr, &
                                  refractivity, &
                                  model_heights, &
                                  noSuperCheck, &
                                  tobs)

INTEGER, INTENT(IN)            :: nlevp                  ! no. of p levels in state vec.
INTEGER, INTENT(IN)            :: nlevq                  ! no. of theta levels
REAL(kind_real), INTENT(IN)    :: za(1:nlevp)            ! heights of rho levs
REAL(kind_real), INTENT(IN)    :: zb(1:nlevq)            ! heights of theta levs
REAL(kind_real), INTENT(IN)    :: pressure(1:nlevp)      ! Model background pressure
REAL(kind_real), INTENT(IN)    :: humidity(1:nlevq)      ! Model background specific humidity
LOGICAL, INTENT(IN)            :: GPSRO_pseudo_ops       ! Option: Use pseudo-levels in vertical interpolation?
LOGICAL, INTENT(IN)            :: GPSRO_vert_interp_ops  ! Option: Use ln(p) for vertical interpolation? (rather than exner)
REAL(kind_real), INTENT(IN)    :: GPSRO_min_temp_grad    ! The minimum temperature gradient which is used
INTEGER, INTENT(IN)            :: nobs                   ! Number of observations in the profile
REAL(kind_real), INTENT(IN)    :: zobs(1:nobs)           ! Impact parameter for the obs
REAL(kind_real), INTENT(IN)    :: RO_Rad_Curv            ! Earth's radius of curvature for these observations
REAL(kind_real), INTENT(IN)    :: Latitude               ! Latitude of this profile
REAL(kind_real), INTENT(IN)    :: RO_geoid_und           ! Undulation - difference between the geoid and the ellipsoid
REAL(kind_real), INTENT(INOUT) :: ycalc(1:nobs)          ! Model forecast of the observations
LOGICAL, INTENT(OUT)           :: BAErr                  ! Was an error encountered during the calculation?
REAL(kind_real), INTENT(INOUT), ALLOCATABLE :: refractivity(:)  ! Refractivity as calculated
REAL(kind_real), INTENT(INOUT), ALLOCATABLE :: model_heights(:) ! Height of the levels for refractivity
LOGICAL, INTENT(IN)            :: noSuperCheck           ! Do we skip a super-refraction check in the operator?
REAL(kind_real), INTENT(OUT)   :: tobs(1:nobs)           ! Virtual temperature on model levels
!
! Things that may need to be output, as they are used by the TL/AD calculation
! 
INTEGER                      :: nRefLevels          ! Number of levels in refractivity calculation
REAL(kind_real), ALLOCATABLE :: nr(:)               ! Model calculation of impact parameters
REAL(kind_real), ALLOCATABLE :: temperature(:)      ! Calculated virtual temperature on pseudo levels
! 
! Local parameters
! 
integer, parameter           :: max_string = 800  ! Length of strings
character(len=*), parameter  :: myname_ = "Ops_GPSRO_ForwardModel"
!
! Local variables
! 
INTEGER                      :: num_pseudo        ! Number of levels, including pseudo levels
REAL(kind_real)              :: x(1:nlevp+nlevq)  ! state vector
character(max_string)        :: err_msg           ! Error message to be output
integer                      :: ilevel            ! Loop variable, model level number
integer                      :: iobs              ! Loop variable, observation number

! The model data must be on a staggered grid, with nlevp = nlevq+1
IF (nlevp /= nlevq + 1) THEN
    write(err_msg,*) myname_ // ':' // ' Data must be on a staggered grid nlevp, nlevq = ', nlevp, nlevq
    call fckit_log % warning(err_msg)
    write(err_msg,*) myname_ // ':' // ' error: number of levels inconsistent!'
    call abor1_ftn(err_msg)
END IF

BAErr = .FALSE.

CALL ufo_calculate_refractivity (nlevp,                 &
                                 nlevq,                 &
                                 za,                    &
                                 zb,                    &
                                 pressure,              &
                                 humidity,              &
                                 GPSRO_pseudo_ops,      &
                                 GPSRO_vert_interp_ops, &
                                 GPSRO_min_temp_grad,   &
                                 BAerr,                 &
                                 nRefLevels,            &
                                 refractivity,          &
                                 model_heights,         &
                                 tpseudo=temperature)

ALLOCATE(nr(1:nRefLevels))

! no point proceeding further if ...
IF (.NOT. BAerr) THEN
    !  2.  Calculate the refractive index * radius on theta model levels (or model impact parameter)
    CALL Ops_GPSROcalc_nr (nRefLevels,    &           ! number of levels for refractivity calculation
                           model_heights, &           ! geopotential heights of refractivity levels
                           refractivity,  &           ! refractivity of model
                           RO_Rad_Curv,   &           ! radius of curvature of earth at observation
                           Latitude,      &           ! latitude at observation
                           RO_geoid_und,  &           ! geoid undulation above WGS-84
                           nr)                        ! Calculated model impact parameters

    !  3.  Calculate model bending angle on observation impact parameters
    CALL Ops_GPSROcalc_alpha (nobs,         &      ! size of ob. vector
                              nRefLevels,   &      ! no. of refractivity levels
                              zobs,         &      ! obs impact parameters
                              refractivity, &      ! refractivity values
                              nr,           &      ! index * radius product
                              ycalc,        &      ! forward modelled bending angle
                              noSuperCheck)        ! Don't use super-refraction check in operator?

    ! 4. Linearly interpolate the virtual temperature to the observation levels
    DO iobs = 1, nobs
        DO iLevel = 1, nRefLevels-1
            IF (nr(iLevel) < zobs(iobs) .AND. nr(iLevel+1) > zobs(iobs)) EXIT
        END DO
        IF (iLevel == nRefLevels) THEN
            tobs(iobs) = missing_value(tobs(iobs))
        ELSE
            tobs(iobs) = temperature(iLevel) + (temperature(iLevel+1) - temperature(iLevel)) * (zobs(iobs) - nr(iLevel)) / (nr(iLevel+1) - nr(iLevel))
        END IF
    END DO
END IF

END SUBROUTINE Ops_GPSRO_ForwardModel

end module ufo_gnssro_bendmetoffice_mod
