!-------------------------------------------------------------------------------
! (C) British Crown Copyright 2020 Met Office
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
!-------------------------------------------------------------------------------
!> Fortran module for gnssro bending angle Met Office's tangent linear and adjoint

module ufo_gnssro_bendmetoffice_tlad_mod

use iso_c_binding
use kinds
use ufo_vars_mod
use ufo_geovals_mod
use ufo_geovals_mod_c,   only: ufo_geovals_registry
use vert_interp_mod
use ufo_basis_tlad_mod,  only: ufo_basis_tlad
use obsspace_mod
use gnssro_mod_conf
use missing_values_mod
use fckit_log_module, only : fckit_log
use ufo_gnssro_bendmetoffice_tlad_utils_mod, only: &
    Ops_GPSROcalc_alphaK, Ops_GPSROcalc_nrK
use ufo_gnssro_ukmo1d_utils_mod, only: Ops_GPSROcalc_nr
use ufo_utils_refractivity_calculator, only: &
    ufo_calculate_refractivity, ufo_refractivity_kmat

private
public :: ufo_gnssro_bendmetoffice_tlad
public :: ufo_gnssro_bendmetoffice_setup
public :: ufo_gnssro_bendmetoffice_tlad_settraj
public :: ufo_gnssro_bendmetoffice_simobs_tl
public :: ufo_gnssro_bendmetoffice_simobs_ad
public :: ufo_gnssro_bendmetoffice_tlad_delete

integer, parameter         :: max_string=800

!> Fortran derived type for gnssro trajectory
type, extends(ufo_basis_tlad)   ::  ufo_gnssro_bendmetoffice_tlad
  private
  logical                       :: vert_interp_ops
  logical                       :: pseudo_ops
  logical                       :: noSuperCheck
  real(kind_real)               :: min_temp_grad
  integer, allocatable          :: chanList(:)
  integer                       :: nlevp
  integer                       :: nlevq
  integer                       :: nlocs
  integer                       :: nlevels            ! Number of vertical levels in the observations
  real(kind_real), allocatable  :: K(:,:)
  contains
    procedure :: setup      => ufo_gnssro_bendmetoffice_setup
    procedure :: delete     => ufo_gnssro_bendmetoffice_tlad_delete
    procedure :: settraj    => ufo_gnssro_bendmetoffice_tlad_settraj
    procedure :: simobs_tl  => ufo_gnssro_bendmetoffice_simobs_tl
    procedure :: simobs_ad  => ufo_gnssro_bendmetoffice_simobs_ad
end type ufo_gnssro_bendmetoffice_tlad

contains

!-------------------------------------------------------------------------------
!> \brief Get the optional settings for the forward model, and save them in the
!         object so that they can be used in the code.
!!
!! \author Neill Bowler (Met Office)
!!
!! \date 26 Aug 2021
!!
!-------------------------------------------------------------------------------
subroutine ufo_gnssro_bendmetoffice_setup(self, vert_interp_ops, pseudo_ops, &
                                          min_temp_grad, chanList, noSuperCheck)

implicit none

class(ufo_gnssro_bendmetoffice_tlad), intent(inout) :: self             !< The object in which to save the variables
logical(c_bool),                         intent(in) :: vert_interp_ops  !< Whether to vertically interpolate using ln(p)
logical(c_bool),                         intent(in) :: pseudo_ops       !< Whether to use pseudo-levels in the calculation
real(c_float),                           intent(in) :: min_temp_grad    !< The minimum temperature gradient in the vertical
integer(c_int),                          intent(in) :: chanList(:)      !< List of channels (vertical levels) to use
logical(c_bool),                         intent(in) :: noSuperCheck     !< If true the don't perform super-refraction check

self % vert_interp_ops = vert_interp_ops
self % pseudo_ops = pseudo_ops
self % min_temp_grad = min_temp_grad
allocate(self % chanList(1:size(chanList)))
self % chanList = chanList
self % noSuperCheck = noSuperCheck

end subroutine ufo_gnssro_bendmetoffice_setup


!-------------------------------------------------------------------------------
!> \brief Calculate the K-matrix (Jacobian) for the observation.
!!
!! \details **ufo_gnssro_bendmetoffice_tlad_settraj**
!! * It is necessary to run this routine before calling the TL or AD routines.
!! * Get the geovals specifying the state around which to linearise, flipping
!!   the vertical order if required.
!! * Call the helper function to calculate the K-matrix for each observation.
!!
!! \author Neill Bowler (Met Office)
!!
!! \date 26 Aug 2021
!!
!-------------------------------------------------------------------------------
subroutine ufo_gnssro_bendmetoffice_tlad_settraj(self, geovals, obss)
       
  use fckit_exception_module, only: fckit_exception

  implicit none
! Subroutine arguments
  class(ufo_gnssro_bendmetoffice_tlad), intent(inout) :: self  !< The object that we use to save data in
  type(ufo_geovals),                intent(in)    :: geovals   !< The input geovals
  type(c_ptr), value,               intent(in)    :: obss      !< The input observations

! Local parameters
  character(len=*), parameter :: myname_="ufo_gnssro_bendmetoffice_tlad_settraj"

! Local variables
  character(max_string)       :: err_msg                       ! Messages to be output to the user
  type(ufo_geoval), pointer   :: q                             ! The model geovals - specific humidity
  type(ufo_geoval), pointer   :: prs                           ! The model geovals - atmospheric pressure
  type(ufo_geoval), pointer   :: rho_heights                   ! The model geovals - heights of the pressure-levels
  type(ufo_geoval), pointer   :: theta_heights                 ! The model geovals - heights of the theta-levels (stores q)
  integer                     :: iobs                          ! Loop variable, observation number
  integer                     :: ilevel                        ! Loop variable, level number
  integer                     :: min_ob                        ! Minimum observation number when passing
  integer                     :: max_ob                        ! Maximum observation number when passing
  integer                     :: this_ob                       ! Ob number, used in flipping geovals

  real(kind_real), allocatable       :: obsLat(:)              ! Latitude of the observation
  real(kind_real), allocatable       :: impact_param(:)        ! Impact parameter of the observation
  real(kind_real), allocatable       :: obsLocR(:)             ! Earth's radius of curvature at the observation tangent point
  real(kind_real), allocatable       :: obsGeoid(:)            ! Undulation - height of the geoid above the ellipsoid

  write(err_msg,*) "TRACE: ufo_gnssro_bendmetoffice_tlad_settraj: begin"
  call fckit_log%info(err_msg)

! Make sure that any previous values of geovals don't get carried over
  call self%delete()

! get model state variables from geovals
  call ufo_geovals_get_var(geovals, var_q,    q)             ! specific humidity
  call ufo_geovals_get_var(geovals, var_prsi, prs)           ! pressure
  call ufo_geovals_get_var(geovals, var_z,    theta_heights) ! Geopotential height of the normal model levels
  call ufo_geovals_get_var(geovals, var_zi,   rho_heights)   ! Geopotential height of the pressure levels

! make sure that the geovals are in the correct vertical order (top-to-bottom)
  if( prs%vals(1,1) > prs%vals(prs%nval,1) ) then 
    write(err_msg,'(a)') 'Geovals should be ordered top to bottom'
    call fckit_exception%throw(err_msg)
  endif

! Keep copy of dimensions
  self % nlevp = prs % nval
  self % nlevq = q % nval
  self % nlocs = obsspace_get_nlocs(obss)
  self % nlevels = max(1, obsspace_get_nchans(obss))
  
! Check that the number of vertical levels in the observations is consistent
! with what we were given in setup
  if (self % nlevels > 1 .AND. self % nlevels /= size(self % chanList)) then
    write(err_msg,'(2A,4I8)') myname_, ' error: channel list must match nlevels', self % nlevels, size(self%chanList)
    call fckit_exception%throw(err_msg)
  end if

! Get the meta-data from the observations
  allocate(obsLat(self%nlocs))
  allocate(impact_param(self%nlevels * self%nlocs))
  allocate(obsLocR(self%nlocs))
  allocate(obsGeoid(self%nlocs))

  call obsspace_get_db(obss, "MetaData", "latitude",         obsLat)
  if (self % nlevels > 1) then
    call obsspace_get_db(obss, "MetaData", "impactParameterRO", impact_param, self%chanList)
  else
    call obsspace_get_db(obss, "MetaData", "impactParameterRO", impact_param)
  end if
  call obsspace_get_db(obss, "MetaData", "earthRadiusCurvature", obsLocR)
  call obsspace_get_db(obss, "MetaData", "geoidUndulation", obsGeoid)
  ALLOCATE(self % K(1:self%nlevels * self%nlocs, 1:prs%nval + q%nval))

! For each observation, calculate the K-matrix
  obs_loop: do iobs = 1, self % nlocs
    min_ob = 1 + (iobs - 1) * self % nlevels
    max_ob = iobs * self % nlevels
    ! Note: The Geovals are passed bottom-to-top as this is the way the bending angle code works
    CALL jacobian_interface(prs % nval, &                              ! Number of pressure levels
                            q % nval, &                                ! Number of specific humidity levels
                            rho_heights % vals(prs%nval:1:-1,iobs), &  ! Heights of the pressure levels
                            theta_heights % vals(q%nval:1:-1,iobs), &  ! Heights of the specific humidity levels
                            q % vals(q%nval:1:-1,iobs), &              ! Values of the specific humidity
                            prs % vals(prs%nval:1:-1,iobs), &          ! Values of the pressure
                            self % pseudo_ops, &                       ! Whether to use pseudo-levels in the calculation
                            self % vert_interp_ops, &                  ! Whether to interpolate using log(pressure)
                            self % min_temp_grad, &                    ! Minimum allowed vertical temperature gradient
                            obsLocR(iobs), &                           ! Local radius of curvature of the earth
                            obsLat(iobs), &                            ! Latitude of the observation
                            obsGeoid(iobs), &                          ! Geoid undulation at the tangent point
                            self % nlevels, &                          ! Number of observations in the profile
                            impact_param(min_ob:max_ob), &             ! Impact parameter for these observations
                            self % K(min_ob:max_ob, &
                                     1:prs%nval+q%nval), &             ! K-matrix (Jacobian of the observation with respect to the inputs)
                            self % noSuperCheck)                       ! If true then don't use super-refraction check
    do ilevel = 1, self % nlevels
      this_ob = ilevel + (iobs - 1) * self % nlevels
      ! Flip the K-matrix back the right way around
      self % K(this_ob, 1:prs%nval) = self % K(this_ob, prs%nval:1:-1)
      self % K(this_ob, prs%nval+1:prs%nval+q%nval) = self % K(this_ob, prs%nval+q%nval:prs%nval+1:-1)
    end do
  end do obs_loop

! Note that this routine has been run.
  self%ltraj = .true.

  deallocate(obsLat)
  deallocate(impact_param)
  deallocate(obsLocR)
  deallocate(obsGeoid)

end subroutine ufo_gnssro_bendmetoffice_tlad_settraj


!-------------------------------------------------------------------------------
!> \brief Given an increment to the model state, calculate an increment to the
!         observation
!!
!! \details **ufo_gnssro_bendmetoffice_simobs_tl**
!! * Get the values from the geovals
!! * Multiply increments by K matrix
!!
!! \author Neill Bowler (Met Office)
!!
!! \date 26 Aug 2021
!!
!-------------------------------------------------------------------------------
subroutine ufo_gnssro_bendmetoffice_simobs_tl(self, geovals, hofx, obss)

  implicit none

! Subroutine arguments
  class(ufo_gnssro_bendmetoffice_tlad), intent(in) :: self      !< Object which is being used to transfer information
  type(ufo_geovals),                intent(in)     :: geovals   !< Model perturbations
  real(kind_real),                  intent(inout)  :: hofx(:)   !< Increment to the observations
  type(c_ptr),   value,             intent(in)     :: obss      !< Input - the observations

! Local parameters
  character(len=*), parameter  :: myname_="ufo_gnssro_bendmetoffice_simobs_tl"

! Local variables
  integer                      :: iobs      ! Loop variable, observation number
  integer                      :: nlocs     ! Number of observations
  character(max_string)        :: err_msg   ! Message to be output
  type(ufo_geoval), pointer    :: q_d       ! Increment to the specific humidity
  type(ufo_geoval), pointer    :: prs_d     ! Increment to the air pressure
  real(kind_real), allocatable :: x_d(:)    ! Increment to the complete state

  write(err_msg,*) "TRACE: ufo_gnssro_bendmetoffice_simobs_tl: begin"
  call fckit_log%info(err_msg)

! Check if trajectory was set
  if (.not. self%ltraj) then
     write(err_msg,*) myname_, ' trajectory wasnt set!'
     call abor1_ftn(err_msg)
  endif
      
! Check if nlocs is consistent in geovals & hofx
  if (geovals%nlocs /= size(hofx)) then
     write(err_msg,*) myname_, ' error: nlocs inconsistent!'
     call abor1_ftn(err_msg)
  endif

! Get variables from geovals
  call ufo_geovals_get_var(geovals, var_q,     q_d)         ! specific humidity
  call ufo_geovals_get_var(geovals, var_prsi,  prs_d)       ! pressure on rho levels

  nlocs = self % nlocs ! number of observations

  allocate(x_d(1:prs_d%nval+q_d%nval))
! Loop through the obs, calculating the increment to the observation
  obs_loop: do iobs = 1, nlocs   ! order of loop doesn't matter

    x_d(1:prs_d%nval) = prs_d % vals(1:prs_d%nval,iobs)
    x_d(prs_d%nval+1:prs_d%nval+q_d%nval) = q_d % vals(1:q_d%nval,iobs)
    hofx(iobs) = SUM(self % K(iobs,:) * x_d)

  end do obs_loop

  deallocate(x_d)

  write(err_msg,*) "TRACE: ufo_gnssro_bendmetoffice_simobs_tl: complete"
  call fckit_log%info(err_msg)

  return
    
end subroutine ufo_gnssro_bendmetoffice_simobs_tl
 

!-------------------------------------------------------------------------------
!> \brief Given an increment to the observation, find the equivalent increment
!         to the model state
!!
!! \details **ufo_gnssro_bendmetoffice_simobs_ad**
!! * Get the geovals to be returning the information in
!! * Multiply observation increments by K
!!
!! \author Neill Bowler (Met Office)
!!
!! \date 26 Aug 2021
!!
!-------------------------------------------------------------------------------
subroutine ufo_gnssro_bendmetoffice_simobs_ad(self, geovals, hofx, obss)

  use typesizes,     only: wp => EightByteReal

  implicit none

! Subroutine arguments
  class(ufo_gnssro_bendmetoffice_tlad), intent(in) :: self      !< Object which is being used to transfer information
  type(ufo_geovals),                intent(inout)  :: geovals   !< Calculated perturbations to model state
  real(kind_real),                  intent(in)     :: hofx(:)   !< Increment to the observations
  type(c_ptr),  value,              intent(in)     :: obss      !< Input - the observations

! Local parameters
  character(len=*), parameter     :: myname_="ufo_gnssro_bendmetoffice_simobs_ad"

! Local variables
  real(c_double)               :: missing  ! Missing data values
  type(ufo_geoval), pointer    :: q_d      ! Pointer to the specific humidity perturbations
  type(ufo_geoval), pointer    :: prs_d    ! Pointer to the pressure perturbations
  integer                      :: iobs     ! Loop variable, observation number
  real(kind_real), allocatable :: x_d(:)   ! Perturbation to the full model state
  character(max_string)        :: err_msg  ! Message to be output

  write(err_msg,*) "TRACE: ufo_gnssro_bendmetoffice_simobs_ad: begin"
  call fckit_log%info(err_msg)

! Check if trajectory was set
  if (.not. self%ltraj) then
     write(err_msg,*) myname_, ' trajectory wasnt set!'
     call abor1_ftn(err_msg)
  endif

! Check if nlocs is consistent in geovals & hofx
  if (geovals%nlocs /= size(hofx)) then
     write(err_msg,*) myname_, ' error: nlocs inconsistent!'
     call abor1_ftn(err_msg)
  endif
     
! Get variables from geovals
  call ufo_geovals_get_var(geovals, var_q,     q_d)         ! specific humidity
  call ufo_geovals_get_var(geovals, var_prsi,  prs_d)       ! pressure

  missing = missing_value(missing)
  allocate(x_d(1:prs_d%nval + q_d%nval))

! Loop through the obs, calculating the increment to the model state
  obs_loop: do iobs = 1, self % nlocs

    if (hofx(iobs) /= missing) then
        x_d = self % K(iobs,:) * hofx(iobs)
        prs_d % vals(:,iobs) = prs_d % vals(:,iobs) + x_d(1:prs_d%nval)
        q_d % vals(:,iobs) = q_d % vals(:,iobs) + x_d(prs_d%nval+1:prs_d%nval+q_d%nval)
    end if

  end do obs_loop

  deallocate(x_d)

  write(err_msg,*) "TRACE: ufo_gnssro_bendmetoffice_simobs_ad: complete"
  call fckit_log%info(err_msg)

  return

end subroutine ufo_gnssro_bendmetoffice_simobs_ad


!-------------------------------------------------------------------------------
!> \brief Tidy up the variables that are used for passing information
!!
!! \author Neill Bowler (Met Office)
!!
!! \date 26 Aug 2021
!!
!-------------------------------------------------------------------------------
subroutine ufo_gnssro_bendmetoffice_tlad_delete(self)

  implicit none
  class(ufo_gnssro_bendmetoffice_tlad), intent(inout) :: self  !< The object being tidied up
  character(len=*), parameter :: myname_="ufo_gnssro_bendmetoffice_tlad_delete"
      
  self%nlocs = 0
  self%nlevp = 0
  self%nlevq = 0
  if (allocated(self%K)) deallocate(self%K)
  if (allocated(self%chanlist)) deallocate(self%chanlist)
  self%ltraj = .false. 

end subroutine ufo_gnssro_bendmetoffice_tlad_delete


!-------------------------------------------------------------------------------
!> \brief Interface for calculating the K-matrix for calculating TL/AD
!!
!! \details **jacobian_interface**
!! * Calculate refractivity and impact parameter (nr)
!! * Call routine to calculate K-matrix
!!
!! \author Neill Bowler (Met Office)
!!
!! \date 26 Aug 2021
!!
!-------------------------------------------------------------------------------
SUBROUTINE jacobian_interface(nlevp, &
                              nlevq, &
                              za, &
                              zb, &
                              q, &
                              prs, &
                              pseudo_ops, &
                              vert_interp_ops, &
                              min_temp_grad, &
                              ro_rad_curv, &
                              latitude, &
                              ro_geoid_und, &
                              nobs, &
                              zobs, &
                              K, &
                              noSuperCheck)

IMPLICIT NONE

INTEGER, INTENT(IN)            :: nlevp            !< The number of model pressure levels
INTEGER, INTENT(IN)            :: nlevq            !< The number of model theta levels
REAL(kind_real), INTENT(IN)    :: za(:)            !< The geometric height of the model pressure levels
REAL(kind_real), INTENT(IN)    :: zb(:)            !< The geometric height of the model theta levels
REAL(kind_real), INTENT(IN)    :: q(1:nlevq)       !< The model values that are being perturbed
REAL(kind_real), INTENT(IN)    :: prs(1:nlevp)     !< The model values that are being perturbed
LOGICAL, INTENT(IN)            :: pseudo_ops       !< Whether to use pseudo levels in the calculation
LOGICAL, INTENT(IN)            :: vert_interp_ops  !< Whether to use exner for the vertical interpolation
REAL(kind_real), INTENT(IN)    :: min_temp_grad    !< The minimum allowed vertical temperature gradient
REAL(kind_real), INTENT(IN)    :: ro_rad_curv      !< The earth's radius of curvature at the ob location
REAL(kind_real), INTENT(IN)    :: latitude         !< The latitude of the ob location
REAL(kind_real), INTENT(IN)    :: ro_geoid_und     !< The geoid undulation at the ob location
INTEGER, INTENT(IN)            :: nobs             !< The number of observations in this column
REAL(kind_real), INTENT(IN)    :: zobs(:)          !< The impact parameters of the column of observations
REAL(kind_real), INTENT(INOUT) :: K(:,:)           !< The calculated K matrix
LOGICAL, INTENT(IN)            :: noSuperCheck     !< If true then don't perform super-refraction check
!
! Things that may need to be output, as they are used by the TL/AD calculation
!
REAL(kind_real), ALLOCATABLE :: model_heights(:)   ! Heights of the pseudo levels
REAL(kind_real), ALLOCATABLE :: refractivity(:)    ! Refractivity on the pseudo levels
INTEGER                      :: nRefLevels         ! Number of pseudo levels
REAL(kind_real)              :: T(1:nlevq)         ! Temperature on model levels
REAL(kind_real), ALLOCATABLE :: nr(:)              ! Model calculation of impact parameters
!
! Local variables
!
INTEGER                      :: num_pseudo        ! Number of levels, including pseudo levels
REAL(kind_real)              :: x(1:nlevp+nlevq)  ! state vector
LOGICAL                      :: BAErr             ! Whether we encountered an error in calculating the refractivity
CHARACTER(LEN=200)           :: err_msg           ! Output message

! Set up the size of the state
x(1:nlevp) = prs
x(nlevp+1:nlevp+nlevq) = q

BAErr = .FALSE.

CALL ufo_calculate_refractivity (nlevp,            &
                                 nlevq,            &
                                 za,               &
                                 zb,               &
                                 prs,              &
                                 q,                &
                                 pseudo_ops,       &
                                 vert_interp_ops,  &
                                 min_temp_grad,    &
                                 BAerr,            &
                                 nRefLevels,       &
                                 refractivity,     &
                                 model_heights)

ALLOCATE(nr(1:nRefLevels))

IF (.NOT. BAErr) THEN
    !  2.  Calculate the impact parameter (refractive index * radius) on refractivity levels
    CALL Ops_GPSROcalc_nr (nRefLevels,    &           ! number of model+pseudo-levels
                           model_heights, &           ! geopotential heights of pseudo levels
                           refractivity,  &           ! refractivity of model on model+pseudo levels
                           RO_Rad_Curv,   &           ! radius of curvature of earth at observation
                           Latitude,      &           ! latitude at observation
                           RO_geoid_und,  &           ! geoid undulation above WGS-84
                           nr)                        ! Calculated model impact parameters

    ! Calculate the K-matrix (Jacobian)
    CALL Ops_GPSRO_GetK(nlevp, &
                        nRefLevels, &
                        nlevq, &
                        za, &
                        zb, &
                        model_heights, &
                        pseudo_ops, &
                        vert_interp_ops, &
                        min_temp_grad, &
                        prs, &
                        q, &
                        ro_rad_curv, &
                        latitude, &
                        ro_geoid_und, &
                        refractivity, &
                        nobs, &
                        zobs, &
                        nr, &
                        K, &
                        noSuperCheck)
ELSE
    K = 0
    write(err_msg,*) "Error in refractivity calculation"
    CALL fckit_log % warning(err_msg)
END IF

DEALLOCATE(nr)
DEALLOCATE(refractivity)
DEALLOCATE(model_heights)

END SUBROUTINE jacobian_interface


!-------------------------------------------------------------------------------
!> \brief Calculate the K-matrix (Jacobian)
!!
!! \details **Ops_GPSRO_GetK**
!! * Calculate the gradient of ref wrt p (on rho levels) and q (on theta levels)
!! * Calculate the gradient of nr wrt ref
!! * Calculate the gradient of bending angle wrt ref and nr
!! * Calculate overall gradient of bending angle wrt p and q
!!
!! \author Neill Bowler (Met Office)
!!
!! \date 26 Aug 2021
!!
!-------------------------------------------------------------------------------
SUBROUTINE Ops_GPSRO_GetK(nlevp, &
                          nRefLevels, &
                          nlevq, &
                          za, &
                          zb, &
                          model_heights, &
                          pseudo_ops, &
                          vert_interp_ops, &
                          min_temp_grad, &
                          pressure, &
                          humidity, &
                          ro_rad_curv, &
                          latitude, &
                          ro_geoid_und, &
                          ref_model, &
                          nobs, &
                          zobs, &
                          nr, &
                          K, &
                          noSuperCheck)
!
! Return the K-matrix for calculating TL/AD
!
    IMPLICIT NONE

    INTEGER, INTENT(IN)          :: nlevp                 !< The number of model pressure levels
    INTEGER, INTENT(IN)          :: nRefLevels            !< Number of refractivity levels
    INTEGER, INTENT(IN)          :: nlevq                 !< The number of model theta levels
    REAL(kind_real), INTENT(IN)  :: za(:)                 !< The geometric height of the model pressure levels
    REAL(kind_real), INTENT(IN)  :: zb(:)                 !< The geometric height of the model theta levels
    REAL(kind_real), INTENT(IN)  :: model_heights(:)      !< The geometric height of the refractivity levels
    LOGICAL, INTENT(IN)          :: pseudo_ops            !< Whether to use pseudo levels in the calculation
    LOGICAL, INTENT(IN)          :: vert_interp_ops       !< Whether to use exner for the vertical interpolation
    REAL(kind_real), INTENT(IN)  :: min_temp_grad         !< Minimum allowed vertical temperature gradient
    REAL(kind_real), INTENT(IN)  :: pressure(nlevp)       !< Model pressure
    REAL(kind_real), INTENT(IN)  :: humidity(nlevq)       !< Model specific humidity
    REAL(kind_real), INTENT(IN)  :: ro_rad_curv           !< The earth's radius of curvature at the ob location
    REAL(kind_real), INTENT(IN)  :: latitude              !< The latitude of the ob location
    REAL(kind_real), INTENT(IN)  :: ro_geoid_und          !< The geoid undulation at the ob location
    REAL(kind_real), INTENT(IN)  :: ref_model(nRefLevels) !< Model refractivity on theta levels - returned from forward model
    INTEGER, INTENT(IN)          :: nobs                  !< The number of observations in this column
    REAL(kind_real), INTENT(IN)  :: zobs(:)               !< The impact parameters of the column of observations
    REAL(kind_real), INTENT(IN)  :: nr(nRefLevels)        !< The impact parameters of the model data
    REAL(kind_real), INTENT(OUT) :: K(nobs,nlevp+nlevq)   !< The calculated K matrix
    LOGICAL, INTENT(IN)          :: noSuperCheck          !< If true, then don't apply super-refraction check

    REAL(kind_real)              :: m1(nobs, nRefLevels)             ! Intermediate term in the K-matrix calculation
    REAL(kind_real), ALLOCATABLE :: dref_dp(:, :)                    ! Partial derivative of refractivity wrt. pressure
    REAL(kind_real), ALLOCATABLE :: dref_dq(:, :)                    ! Partial derivative of refractivity wrt. specific humidity
    REAL(kind_real)              :: dnr_dref(nRefLevels, nRefLevels) ! Partial derivative of impact parameter wrt. refractivity
    REAL(kind_real)              :: dalpha_dref(nobs, nRefLevels)    ! Partial derivative of bending angle wrt. refractivity
    REAL(kind_real)              :: dalpha_dnr(nobs, nRefLevels)     ! Partial derivative of bending angle wrt. impact parameter

    !  1.  Calculate the gradient of ref wrt p (on rho levels) and q (on theta levels)
    CALL ufo_refractivity_kmat(nlevp,      &
                               nlevq,      &
                               nRefLevels, &
                               za,         &
                               zb,         &
                               pressure,   &
                               humidity,   &
                               pseudo_ops, &
                               vert_interp_ops, &
                               min_temp_grad, &
                               dref_dp,    &       !out
                               dref_dq)            !out

    !  2.  Calculate the gradient of nr wrt ref
    CALL Ops_GPSROcalc_nrK (model_heights, &       ! geopotential heights of pseudo levels
                            nRefLevels,    &       ! number of refractivity levels
                            RO_Rad_Curv,   &       ! radius of curvature of earth at observation
                            Latitude,      &       ! latitude at observation
                            RO_geoid_und,  &       ! geoid undulation above WGS-84
                            ref_model,     &       ! refractivity of model on model levels
                            dnr_dref)              ! out

    !  3.  Calculate the gradient of bending angle wrt ref and nr
    CALL Ops_GPSROcalc_alphaK (nobs,        &      ! size of ob. vector
                               nRefLevels,  &      ! no. of refractivity levels
                               zobs,        &      ! obs impact parameters
                               ref_model,   &      ! refractivity values on model levels
                               nr,          &      ! index * radius product
                               dalpha_dref, &      ! out
                               dalpha_dnr,  &      ! out
                               noSuperCheck)       ! Don't use super-refraction check in operator?

    ! Calculate overall gradient of bending angle wrt p and q
    m1 = MATMUL (dalpha_dnr,dnr_dref)
    K(1:nobs, 1:nlevp) = MATMUL (dalpha_dref,dref_dp) + MATMUL (m1,dref_dp)    !P part
    K(1:nobs, nlevp+1:nlevp+nlevq) = MATMUL (dalpha_dref,dref_dq) + MATMUL (m1,dref_dq) !q part

    IF (ALLOCATED(dref_dp)) DEALLOCATE(dref_dp)
    IF (ALLOCATED(dref_dq)) DEALLOCATE(dref_dq)

END SUBROUTINE Ops_GPSRO_GetK

end module ufo_gnssro_bendmetoffice_tlad_mod
