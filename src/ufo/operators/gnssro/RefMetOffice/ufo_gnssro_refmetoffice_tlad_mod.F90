!-------------------------------------------------------------------------------
! (C) British Crown Copyright 2021 Met Office
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
!-------------------------------------------------------------------------------
!> Fortran module for gnssro refractivity Met Office's tangent linear and adjoint

module ufo_gnssro_refmetoffice_tlad_mod

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
use fckit_exception_module, only: fckit_exception
use ufo_utils_refractivity_calculator, only: &
    ufo_calculate_refractivity, ufo_refractivity_partial_derivatives
use ufo_constants_mod, only: &
    rd,                      &    ! Gas constant for dry air
    cp,                      &    ! Heat capacity at constant pressure for air
    rd_over_cp,              &    ! Ratio of gas constant to heat capacity
    grav,                    &    ! Gravitational field strength
    pref,                    &    ! Reference pressure for calculating exner
    mw_ratio,                &    ! Ratio of molecular weights of water and dry air
    c_virtual,               &    ! Related to mw_ratio
    n_alpha,                 &    ! Refractivity constant a
    n_beta                        ! Refractivity constant b

private
public :: ufo_gnssro_refmetoffice_tlad
public :: ufo_gnssro_refmetoffice_tlad_setup
public :: ufo_gnssro_refmetoffice_tlad_settraj
public :: ufo_gnssro_refmetoffice_simobs_tl
public :: ufo_gnssro_refmetoffice_simobs_ad
public :: ufo_gnssro_refmetoffice_tlad_delete

integer, parameter :: max_string=800

!> Fortran derived type for gnssro trajectory
type, extends(ufo_basis_tlad)   ::  ufo_gnssro_refmetoffice_tlad
  private
  logical :: vert_interp_ops              !< Do vertical interpolation using ln(p) or exner?
  logical :: pseudo_ops                   !< Use pseudo-levels in the refractivity calculation
  real(kind_real) :: min_temp_grad        !< The minimum temperature gradient before a profile is considered isothermal
                                          !  Used in pseudo-levels calculation
  integer                       :: nlevp  !< The number of pressure levels
  integer                       :: nlevq  !< The number of specific humidity (or temperature) levels
  integer                       :: nlocs  !< The number of observations
  real(kind_real), allocatable  :: K(:,:) !< The K-matrix (Jacobian)
  contains
    procedure :: setup      => ufo_gnssro_refmetoffice_tlad_setup
    procedure :: delete     => ufo_gnssro_refmetoffice_tlad_delete
    procedure :: settraj    => ufo_gnssro_refmetoffice_tlad_settraj
    procedure :: simobs_tl  => ufo_gnssro_refmetoffice_simobs_tl
    procedure :: simobs_ad  => ufo_gnssro_refmetoffice_simobs_ad
end type ufo_gnssro_refmetoffice_tlad

contains

!-------------------------------------------------------------------------------
!> \brief Set up the Met Office GNSS-RO refractivity TL/AD
!!
!! \details **ufo_gnssro_refmetoffice_tlad_setup**
!! * Get the optional settings for the forward model and its linear, and save
!!   them in the object so that they can be used in the code.
!!
!! \author Neill Bowler (Met Office)
!!
!! \date 26 May 2021
!!
!-------------------------------------------------------------------------------
subroutine ufo_gnssro_refmetoffice_tlad_setup(self, vert_interp_ops, pseudo_ops, min_temp_grad)

implicit none

class(ufo_gnssro_refmetoffice_tlad), intent(inout) :: self
logical(c_bool), intent(in) :: vert_interp_ops
logical(c_bool), intent(in) :: pseudo_ops
real(c_float), intent(in) :: min_temp_grad

self % vert_interp_ops = vert_interp_ops
self % pseudo_ops = pseudo_ops
self % min_temp_grad = min_temp_grad

end subroutine ufo_gnssro_refmetoffice_tlad_setup

!-------------------------------------------------------------------------------
!> \brief Set up the K-matrix for (Jacobian) for the Met Office's GNSS-RO
!!        refractivity operator
!!
!! \details **ufo_gnssro_refmetoffice_tlad_settraj**
!! * It is necessary to run this routine before calling the TL or AD routines.
!! * Get the geovals specifying the state around which to linearise, flipping
!!   the vertical order as refractivity is calculated starting from the surface.
!! * Call the helper function to calculate the K-matrix for each observation.
!!
!! \author Neill Bowler (Met Office)
!!
!! \date 26 May 2021
!!
!-------------------------------------------------------------------------------
subroutine ufo_gnssro_refmetoffice_tlad_settraj(self, geovals, obss)
       
  implicit none
! Subroutine arguments
  class(ufo_gnssro_refmetoffice_tlad), intent(inout) :: self     !< The object that we use to save data in
  type(ufo_geovals),                intent(in)       :: geovals  !< The input geovals
  type(c_ptr), value,               intent(in)       :: obss     !< The input observations

! Local parameters
  character(len=*), parameter :: myname_="ufo_gnssro_refmetoffice_tlad_settraj"

! Local variables
  character(max_string)       :: err_msg                       ! Messages to be output to the user
  type(ufo_geoval), pointer   :: q                             ! The model geovals - specific humidity
  type(ufo_geoval), pointer   :: prs                           ! The model geovals - atmospheric pressure
  type(ufo_geoval), pointer   :: rho_heights                   ! The model geovals - heights of the pressure-levels
  type(ufo_geoval), pointer   :: theta_heights                 ! The model geovals - heights of the theta-levels (stores q)
  integer                     :: iobs                          ! Loop variable, observation number

  real(kind_real), allocatable       :: obs_height(:)          ! Geopotential height of the observation

  write(err_msg,*) "TRACE: ufo_gnssro_refmetoffice_tlad_settraj: begin"
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
  
! Get the meta-data from the observations
  allocate(obs_height(self%nlocs))
  call obsspace_get_db(obss, "MetaData", "height", obs_height)
  allocate(self % K(1:self%nlocs, 1:prs%nval + q%nval))

! For each observation, calculate the K-matrix
  obs_loop: do iobs = 1, self % nlocs
    CALL jacobian_interface(prs % nval, &                              ! Number of pressure levels
                            q % nval, &                                ! Number of specific humidity levels
                            rho_heights % vals(prs%nval:1:-1,iobs), &  ! Heights of the pressure levels
                            theta_heights % vals(q%nval:1:-1,iobs), &  ! Heights of the specific humidity levels
                            prs % vals(prs%nval:1:-1,iobs), &          ! Values of the pressure
                            q % vals(q%nval:1:-1,iobs), &              ! Values of the specific humidity
                            self % pseudo_ops, &                       ! Whether to use pseudo-levels in the calculation
                            self % vert_interp_ops, &                  ! Whether to interpolate using log(pressure)
                            self % min_temp_grad, &                    ! Minimum allowed vertical temperature gradient
                            1, &                                       ! Number of observations in the profile
                            obs_height(iobs:iobs), &                   ! Impact parameter for this observation
                            self % K(iobs:iobs,1:prs%nval+q%nval))     ! K-matrix (Jacobian of the observation with respect to the inputs)
    ! Flip the K-matrix back the right way around
    self % K(iobs,1:prs%nval) = self % K(iobs, prs%nval:1:-1)
    self % K(iobs,prs%nval+1:prs%nval+q%nval) = self % K(iobs, prs%nval+q%nval:prs%nval+1:-1)
  end do obs_loop

! Note that this routine has been run.
  self%ltraj = .true.

  deallocate(obs_height)

end subroutine ufo_gnssro_refmetoffice_tlad_settraj

!-------------------------------------------------------------------------------
!> \brief Given an increment to the model state, calculate an increment to the
!!        observation
!!
!! \details **ufo_gnssro_refmetoffice_simobs_tl**
!! * Check that set trajectory has been previously called.
!! * Get the geovals for the increment.
!! * For each observation apply the K-matrix to calculate the increment to the
!!   observation.
!!
!! \author Neill Bowler (Met Office)
!!
!! \date 26 May 2021
!!
!-------------------------------------------------------------------------------

subroutine ufo_gnssro_refmetoffice_simobs_tl(self, geovals, hofx, obss)

  implicit none

! Subroutine arguments
  class(ufo_gnssro_refmetoffice_tlad), intent(in) :: self       !< Object which is being used to transfer information
  type(ufo_geovals),                intent(in)    :: geovals    !< Model perturbations
  real(kind_real),                  intent(inout) :: hofx(:)    !< Increment to the observations
  type(c_ptr),   value,             intent(in)    :: obss       !< Input - the observations

! Local parameters
  character(len=*), parameter  :: myname_="ufo_gnssro_refmetoffice_simobs_tl"

! Local variables
  integer                      :: iobs      ! Loop variable, observation number
  integer                      :: nlocs     ! Number of observations
  character(max_string)        :: err_msg   ! Message to be output
  type(ufo_geoval), pointer    :: q_d       ! Increment to the specific humidity
  type(ufo_geoval), pointer    :: prs_d     ! Increment to the air pressure
  real(kind_real), allocatable :: x_d(:)    ! Increment to the complete state

  write(err_msg,*) "TRACE: ufo_gnssro_refmetoffice_simobs_tl: begin"
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

    x_d(1:prs_d%nval) = prs_d % vals(:,iobs)
    x_d(prs_d%nval+1:prs_d%nval+q_d%nval) = q_d % vals(:,iobs)
    hofx(iobs) = SUM(self % K(iobs,:) * x_d)

  end do obs_loop

  deallocate(x_d)

  write(err_msg,*) "TRACE: ufo_gnssro_refmetoffice_simobs_tl: complete"
  call fckit_log%info(err_msg)

  return
    
end subroutine ufo_gnssro_refmetoffice_simobs_tl

!-------------------------------------------------------------------------------
!> \brief Given an increment to the observation, find the equivalent increment
!!        to the model state
!!
!! \details **ufo_gnssro_refmetoffice_simobs_ad**
!! * Check that set trajectory has previously been called.
!! * Get the geovals for the model increment, and allocate these if they have
!!   not already been allocated.
!! * For each observation calculate the observation increment using the K-matrix
!!
!! \author Neill Bowler (Met Office)
!!
!! \date 26 May 2021
!!
!-------------------------------------------------------------------------------

subroutine ufo_gnssro_refmetoffice_simobs_ad(self, geovals, hofx, obss)

  use typesizes,     only: wp => EightByteReal

  implicit none

! Subroutine arguments
  class(ufo_gnssro_refmetoffice_tlad), intent(in) :: self       !< Object which is being used to transfer information
  type(ufo_geovals),                intent(inout) :: geovals    !< Calculated perturbations to model state
  real(kind_real),                  intent(in)    :: hofx(:)    !< Increment to the observations
  type(c_ptr),  value,              intent(in)    :: obss       !< Input - the observations

! Local parameters
  character(len=*), parameter     :: myname_="ufo_gnssro_refmetoffice_simobs_ad"

! Local variables
  real(c_double)               :: missing  ! Missing data values
  type(ufo_geoval), pointer    :: q_d      ! Pointer to the specific humidity perturbations
  type(ufo_geoval), pointer    :: prs_d    ! Pointer to the pressure perturbations
  integer                      :: iobs     ! Loop variable, observation number
  real(kind_real), allocatable :: x_d(:)   ! Perturbation to the full model state
  character(max_string)        :: err_msg  ! Message to be output

  write(err_msg,*) "TRACE: ufo_gnssro_refmetoffice_simobs_ad: begin"
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

  write(err_msg,*) "TRACE: ufo_gnssro_refmetoffice_simobs_ad: complete"
  call fckit_log%info(err_msg)

  return

end subroutine ufo_gnssro_refmetoffice_simobs_ad

!-------------------------------------------------------------------------------
!> \brief Tidy up the variables that are used for passing information
!!
!! \details **ufo_gnssro_refmetoffice_tlad_delete**
!! * Set lengths to zero, and deallocate K-matrix
!!
!! \author Neill Bowler (Met Office)
!!
!! \date 26 May 2021
!!
!-------------------------------------------------------------------------------

subroutine ufo_gnssro_refmetoffice_tlad_delete(self)

  implicit none
  class(ufo_gnssro_refmetoffice_tlad), intent(inout) :: self
  character(len=*), parameter :: myname_="ufo_gnssro_refmetoffice_tlad_delete"
      
  self%nlocs = 0
  self%nlevp = 0
  self%nlevq = 0
  if (allocated(self%K)) deallocate(self%K)
  self%ltraj = .false. 

end subroutine ufo_gnssro_refmetoffice_tlad_delete

!-------------------------------------------------------------------------------
!> \brief Calculate the K-matrix used in the TL/AD
!!
!! \details **jacobian_interface**
!! * Allocate temporary arrays
!! * Call partial-derivatives code which calculates various quantities on model
!!   levels
!! * Loop over each observation, doing the following
!! * If using pseudo-levels, calculate specific humidity, pressure and
!!   temperature on intermediate levels, and interpolate these appropriately.
!!   Then calculate the refractivity gradients for each observation,
!!   interpolated from the pseudo-levels.
!! * If not using pseudo-levels, assume that refractivity varies exponentially
!!   with height and calculate gradients.
!! * Calculate K-matrix from the component gradients.
!!
!! \author Neill Bowler (Met Office)
!!
!! \date 26 May 2021
!!
!-------------------------------------------------------------------------------

SUBROUTINE jacobian_interface(nlevp,  &
                              nlevq,  &
                              za,     &
                              zb,     &
                              p,      &
                              q,      &
                              pseudo_ops, &
                              vert_interp_ops, &
                              min_temp_grad, &
                              nobs,   &
                              zobs,   &
                              Kmat)

IMPLICIT NONE

! Subroutine arguments:
integer, intent(in)          :: nlevp            !< Number of pressure levels
integer, intent(in)          :: nlevq            !< Number of specific humidity levels
real(kind_real), intent(in)  :: za(:)            !< Height of the pressure levels
real(kind_real), intent(in)  :: zb(:)            !< Height of the specific humidity levels
real(kind_real), intent(in)  :: p(:)             !< Input pressure
real(kind_real), intent(in)  :: q(:)             !< Input specific humidity
logical, intent(in)          :: vert_interp_ops  !< Use log(p) for vertical interpolation?
logical, intent(in)          :: pseudo_ops       !< Use pseudo-levels to reduce errors?
real(kind_real), intent(in)  :: min_temp_grad    !< Minimum value for the vertical temperature gradient
integer, intent(in)          :: nobs             !< Number of observations
real(kind_real), intent(in)  :: zobs(:)          !< Height of the observations
real(kind_real), intent(out) :: Kmat(:,:)        !< K-matrix (Jacobian)

! Local declarations:
integer                      :: n                          ! Loop variable, observation number
integer                      :: i                          ! Loop variable
real(kind_real)              :: Exner(nlevp)               ! Exner pressure, calculated on model pressure levels
real(kind_real)              :: Pb(nlevq)                  ! Pressure on model theta levels
real(kind_real)              :: Tv(nlevq)                  ! Virtual temperature on model theta levels
real(kind_real)              :: T(nlevq)                   ! Temperature on model theta levels
real(kind_real)              :: dExtheta_dPb(nlevq,nlevq)  ! Partial derivative of theta wrt pressure (on model theta levels)
real(kind_real), ALLOCATABLE :: drefob_dPob(:,:)           ! Partial derivative of refractivity wrt pressure (at ob location)
real(kind_real), ALLOCATABLE :: drefob_dTob(:,:)           ! Partial derivative of refractivity wrt temperature (at ob location)
real(kind_real), ALLOCATABLE :: drefob_dqob(:,:)           ! Partial derivative of refractivity wrt specific humidity (at ob location)
real(kind_real), ALLOCATABLE :: dqob_dqb(:,:)              ! Partial derivative of specific humidity at ob location wrt specific humidity on model theta levels
real(kind_real), ALLOCATABLE :: dTob_dTb(:,:)              ! Partial derivative of temperature at ob location wrt temperature on model theta levels
real(kind_real), ALLOCATABLE :: dPob_dT(:,:)               ! Partial derivative of pressure at ob location wrt temperature on model theta levels
real(kind_real), ALLOCATABLE :: dPob_dPb(:,:)              ! Partial derivative of pressure at ob location wrt pressure on model theta levels
real(kind_real), ALLOCATABLE :: dPob_dP(:,:)               ! Partial derivative of pressure at ob location wrt pressure on model pressure levels
real(kind_real), ALLOCATABLE :: dTob_dP(:,:)               ! Partial derivative of temperature at ob location wrt pressure on model pressure levels
real(kind_real), ALLOCATABLE :: dTb_dp(:,:)                ! Partial derivative of temperature on model theta levels wrt pressure on model pressure levels
real(kind_real)              :: c                          ! Continuity constant for hydrostatic pressure
real(kind_real)              :: P_ob                       ! Model pressure at the ob location
real(kind_real)              :: q_ob                       ! Model specific humidity at the ob location
real(kind_real)              :: T_ob                       ! Model temperature at the ob location
real(kind_real)              :: gamma                      ! Vertical gradient of log(q)
real(kind_real)              :: beta                       ! Vertical temperature gradient
real(kind_real)              :: Ndry                       ! Component of refractivity due to dry terms
real(kind_real)              :: Nwet                       ! Component of refractivity due to wet terms
real(kind_real)              :: refrac(nlevq)              ! Refractivity on model theta levels
real(kind_real)              :: dEx_dP(nlevp,nlevp)        ! Partial derivative of exner wrt pressure
real(kind_real)              :: dPb_dp(nlevq,nlevp)        ! Partial derivative of pressure on theta levels wrt pressure on pressure levels
real(kind_real)              :: dTv_dExtheta(nlevq,nlevq)  ! Virtual temperature divided by exner on theta levels
real(kind_real)              :: dTv_dEx(nlevq,nlevp)       ! Partial derivative of virtual temperature wrt exner
real(kind_real)              :: dT_dTv(nlevq,nlevq)        ! Partial derivative of temperature wrt virtual temperature
real(kind_real)              :: dT_dq(nlevq,nlevq)         ! Partial derivative of temperature wrt specific humidity
real(kind_real)              :: dref_dPb(nlevq,nlevq)      ! Partial derivative of refractivity wrt pressure on theta levels
real(kind_real)              :: dref_dT(nlevq,nlevq)       ! Partial derivative of refractivity wrt temperature
real(kind_real)              :: dref_dq(nlevq,nlevq)       ! Partial derivative of refractivity wrt specific humidity
real(kind_real)              :: dNref_dref(nobs,nlevq)     ! Partial derivative of refractivity at the ob loation wrt model refractivity
real(kind_real)              :: kmatP(nobs,nlevp)          ! K-matrix contribution for pressure terms
real(kind_real)              :: kmatq(nobs,nlevq)          ! K-matrix contribution for specific humidity terms
real(kind_real)              :: Nref(nobs)                 ! Model refractivity at observation locations
real(kind_real)              :: m1(nobs,nlevq)             ! Temporary matrix product used in calculation
real(kind_real)              :: m2(nobs,nlevp)             ! Temporary matrix product used in calculation
real(kind_real)              :: m3(nobs,nlevq)             ! Temporary matrix product used in calculation
real(kind_real)              :: m4(nobs,nlevq)             ! Temporary matrix product used in calculation
real(kind_real)              :: g_RB                       ! Frequently used term
real(kind_real)              :: c_ZZ                       ! Frequently used term
real(kind_real)              :: dPo_dT1                    ! dP_ob / dT_below
real(kind_real)              :: dPo_dTo                    ! dP_ob / dT_ob
real(kind_real)              :: dTo_dT1                    ! dT_ob / dT_below
real(kind_real)              :: dPo_dbeta                  ! dP_ob / dbeta
real(kind_real)              :: dbeta_dT1                  ! dbeta / dT_below
real(kind_real)              :: dTo_dbeta                  ! dT_ob / dbeta
real(kind_real)              :: dbeta_dT2                  ! dbeta / dT_above
real(kind_real)              :: dPo_dc                     ! dP_ob / dc
real(kind_real)              :: dc_dT1                     ! dc / dT_below
real(kind_real)              :: dc_dbeta                   ! dc / dbeta
real(kind_real)              :: dc_dt2                     ! dc / dT_above
real(kind_real)              :: dPo_dP1                    ! dP_ob / dP_below
real(kind_real)              :: dc_dP1                     ! dc / dP_below
real(kind_real)              :: dc_dP2                     ! dc / dP_above
real(kind_real)              :: model_height_diff          ! Difference in height between two model levels
real(kind_real)              :: obs_height_diff            ! Difference in height between the observation and the model level below it
real(kind_real)              :: height_diff_ratio          ! obs_height_diff / model_height_diff

!-----------------------
! 1. Allocate/initialise matrices
!-----------------------

IF (pseudo_ops) THEN
  ALLOCATE (drefob_dPob(nobs,nobs))
  ALLOCATE (drefob_dTob(nobs,nobs))
  ALLOCATE (drefob_dqob(nobs,nobs))
  ALLOCATE (dqob_dqb(nobs,nlevq))
  ALLOCATE (dTob_dTb(nobs,nlevq))
  ALLOCATE (dPob_dT(nobs,nlevq))
  ALLOCATE (dPob_dPb(nobs,nlevq))
  ALLOCATE (dPob_dP(nobs,nlevp))
  ALLOCATE (dTob_dP(nobs,nlevp))
  ALLOCATE (dTb_dP(nlevq,nlevp))

  drefob_dPob(:,:) = 0.0
  drefob_dTob(:,:) = 0.0
  drefob_dqob(:,:) = 0.0
  dqob_dqb(:,:) = 0.0
  dTob_dTb(:,:) = 0.0
  dPob_dT(:,:) = 0.0
  dPob_dPb(:,:) = 0.0
  dPob_dP(:,:) = 0.0
  dTob_dP(:,:) = 0.0
  dTb_dP(:,:) = 0.0
END IF

call ufo_refractivity_partial_derivatives(nlevP,           &
                                          nlevq,           &
                                          za,              &
                                          zb,              &
                                          P,               &
                                          q,               &
                                          vert_interp_ops, &
                                          dT_dTv,          &
                                          dT_dq,           &
                                          dref_dPb,        &
                                          dref_dT,         &
                                          dref_dq,         &
                                          refrac,          &
                                          T,               &
                                          Pb,              &
                                          dEx_dP,          &
                                          dExtheta_dPb,    &
                                          dTv_dExtheta,    &
                                          dPb_dP,          &
                                          dTv_dEx)


dNref_dref(:,:) = 0.0
kmatP(:,:) = 0.0
kmatq(:,:) = 0.0
kmat(:,:) = 0.0

!-------------------------------------------------
! 2. Calculate the partial derivatives at the observation locations
!-------------------------------------------------

DO n = 1, nobs
  i = 1

  if (zobs(n) >= zb(nlevq) .or. zobs(n) == missing_value(zobs(n))) cycle
  DO
    ! Find model layer containing ob
    IF (zobs(n) < zb(i + 1)) EXIT

    i = i + 1

  END DO

  IF (pseudo_ops) THEN
    ! Calculate some quantities that will be re-used
    model_height_diff = zb(i + 1) - zb(i)
    obs_height_diff = zobs(n) - zb(i)
    height_diff_ratio = obs_height_diff / model_height_diff
  
    ! Interpolate P,T,q separately
    IF (MIN (q(i), q(i + 1)) > 0.0) THEN
      ! q varies exponentially with height
      gamma = LOG (q(i) / q(i + 1)) / model_height_diff
      q_ob = q(i) * EXP (-gamma * obs_height_diff)

      dqob_dqb(n,i) = (q_ob / q(i)) * (1.0 - height_diff_ratio)
      dqob_dqb(n,i + 1) = (q_ob / q(i + 1)) * height_diff_ratio

    ! Assume linear variation if humidities are -ve to avoid singularity
    ELSE
      q_ob = q(i) + (q(i + 1) - q(i)) * height_diff_ratio

      dqob_dqb(n,i) = (zb(i + 1) - zobs(n)) / model_height_diff
      dqob_dqb(n,i + 1) = height_diff_ratio
    END IF

    ! T varies linearly with height
    beta = (T(i + 1) - T(i)) / model_height_diff
    T_ob = T(i) + beta * obs_height_diff

    dTob_dTb(n,i) = (zb(i + 1) - zobs(n)) / model_height_diff
    dTob_dTb(n,i + 1) = height_diff_ratio

    ! P varies to maintain hydrostatic balance
    IF (ABS (T(i + 1) - T(i)) > min_temp_grad) THEN
      c = ((Pb(i + 1) / Pb(i)) * (T(i + 1)/T(i)) ** (grav / (rd * beta)) - 1.0) / model_height_diff
      P_ob = (Pb(i) * (T_ob / T(i)) ** &
                  (-grav / (rd * beta))) * (1.0 + c * obs_height_diff)
    ELSE
      ! If layer is isothermal, assume exponential variation to
      ! avoid singularity
      P_ob = Pb(i) * EXP(LOG(Pb(i + 1) / pb(i)) * height_diff_ratio)
    END IF

    ! Temporary variables to keep equations neater
    g_RB = grav / (rd * beta)
    c_ZZ = c + 1.0 / model_height_diff

    dPo_dT1 = (g_RB / t(i)) * P_ob
    dPo_dTo = -(g_RB / t_ob) * p_ob
    dTo_dT1 = dTob_dTb(n,i)
    dPo_dbeta = (g_RB / beta) * LOG (T_ob / T(i)) * P_ob
    dbeta_dT1 = -1.0 / model_height_diff
    dTo_dbeta = obs_height_diff
    dbeta_dT2 = 1.0 / model_height_diff

    ! Pressure Jacobians for hydrostatic/exponential cases
    IF (ABS (T(i + 1) - T(i)) > min_temp_grad) THEN
      dPo_dc = obs_height_diff * Pb(i) * (T_ob / T(i)) ** (-g_RB)
      dc_dT1 = -(g_RB / (T(i))) * (c_ZZ)
      dc_dbeta = -(g_RB / beta) * LOG (T(i + 1) / T(i)) * c_ZZ
      dc_dT2 = (g_RB / T(i + 1)) * c_ZZ
      dPo_dP1 = P_ob / Pb(i)
      dPob_dT(n,i) = dPo_dT1 + dPo_dTo * dTo_dT1 + dPo_dbeta * dbeta_dT1 + dPo_dc * &
                     (dc_dT1 + dc_dbeta * dbeta_dT1)
      dPob_dT(n,i+1) = (dPo_dbeta + dPo_dTo * dTo_dbeta + dPo_dc * dc_dbeta) * &
                     dbeta_dT2 + dPo_dc * dc_dT2

      dc_dP1 = -(1.0 / Pb(i)) * c_ZZ
      dc_dP2 = (1.0 / Pb(i + 1)) * c_ZZ

      dPob_dPb(n,i) = dPo_dP1 + dPo_dc * dc_dP1
      dPob_dPb(n,i + 1) = dPo_dc * dc_dP2
    ELSE
      dPob_dPb(n,i) = EXP (LOG (Pb(i + 1) / Pb(i)) * height_diff_ratio) * &
        (1.0 - height_diff_ratio)
      dPob_dPb(n,i + 1) = (Pb(i) / Pb(i + 1)) * height_diff_ratio * &
        EXP(LOG(Pb(i + 1) / Pb(i)) * height_diff_ratio)
    END IF

    ! Calculate refractivity
    Ndry = n_alpha * P_ob / T_ob
    Nwet = n_beta * P_ob * q_ob / (T_ob ** 2 * (mw_ratio + (1.0 - mw_ratio) * q_ob))
    Nref(n) = Ndry + Nwet

    drefob_dPob(n,n) = Nref(n) / P_ob
    drefob_dTob(n,n) = -(Ndry + 2.0 * Nwet) / T_ob
    drefob_dqob(n,n) = n_beta * p_ob * mw_ratio / (T_ob * (mw_ratio + (1.0 - mw_ratio) * q_ob)) ** 2

  ELSE

    ! Use simple assumption of exponentially varying refractivity

    gamma = (zb(i + 1) - zobs(n)) / (zb(i + 1) - zb(i))

    beta = 1.0 - gamma

    Nref(n) = EXP (gamma * LOG (refrac(i)) + beta * LOG (refrac(i + 1)))

    dNref_dref(n,i) = Nref(n) * gamma / refrac(i)

    dNref_dref(n,i + 1) = Nref(n) * beta / refrac(i + 1)
  END IF

END DO

!-------------------------------------------------
! 3. Evaluate the Kmatrix by matrix multiplication
!-------------------------------------------------

IF (pseudo_ops) THEN

  ! Derivatives:
  ! dPob/dP = (dPob/dPb * dPb/dP) + (dPob/dT * dT/dTv * dTv/dEx * dEx/dP)
  dPob_dP = MATMUL (dPob_dPb, dPb_dP) + MATMUL (dPob_dT, MATMUL (dT_dTv, MATMUL (dTv_dEx, dEx_dP)))
  ! dTob/dP = (dTob/dT * dT/dTv * dTv/dEx * dEx/dP)
  dTob_dP = MATMUL (dTob_dTb, MATMUL (dT_dTv, MATMUL (dTv_dEx, dEx_dP)))

  ! calc K matrix for p on rho levels
  ! dNob/dP = (dNob/dPob * dPob/dP) + (dNob/dTob *dTob/dP)
  KmatP(:,:) = MATMUL (drefob_dPob, dPob_dP) + MATMUL (drefob_dTob, dTob_dP)

  ! calc Kmatrix for q on theta levels
  ! dNob/dq = (dNob/dqob * dqob/dq) + (dNob/dTob * dTob/dT * dT/dq) + (dNob/dPob * dPob/dT * dT/dq)
  Kmatq(:,:) = MATMUL (drefob_dqob, dqob_dqb) + MATMUL (MATMUL (drefob_dTob, dTob_dTb), dT_dq) + &
               MATMUL (MATMUL (drefob_dPob, dPob_dT), dT_dq)

ELSE

  ! calc K matrix for p on rho levels
  ! dNob/dP = (dNob/dN * dN/dPb * dPb/dP)   + ....
  KmatP(:,:) = MATMUL (dNref_dref, MATMUL (dref_dPb, dPb_dP))

  !  .... (dNob/dN * dN/dT * dT/dTv * dTv/dEx * dEx/dP) + .....
  m1(:,:) = MATMUL (dNref_dref, MATMUL (dref_dT, dT_dTv))
  m2(:,:) = MATMUL (m1, dTv_dEx)
  KmatP(:,:) = KmatP(:,:) + MATMUL (m2, dEx_dP)

  !  .... (dNob/dN * dN/dT * dT/dTv * dTv/dExtheta * dExtheta/dPb * dPb/dP)
  m3(:,:) = MATMUL (m1, dTv_dExtheta)
  m4(:,:) = MATMUL (m3, dExtheta_dPb)
  KmatP(:,:) = KmatP(:,:) + MATMUL (m4,dPb_dP)

  ! calc Kmatrix for q on theta levels
  !  dNob/dq = (dNob/dN * dN/dq) + (dNob/dN * dN/dT * dT/dq)
  Kmatq(:,:) = Kmatq(:,:) + MATMUL (dNref_dref, MATMUL (dref_dT, dT_dq))

END IF

! Calculate the full kmatrix in correct units

Kmat(1:nobs, 1:nlevp) = 1.0E2 * Kmatp(1:nobs, 1:nlevp)   ! h/Pa

Kmat(1:nobs, nlevp+1:nlevp+nlevq) = 1.0E-3 * Kmatq(1:nobs, 1:nlevq)   ! g/kg

!-------------------------------------------------------------------------------
! 4. Deallocate the temporary arrays
!-------------------------------------------------------------------------------

IF (ALLOCATED (drefob_dPob)) DEALLOCATE (drefob_dPob)
IF (ALLOCATED (drefob_dTob)) DEALLOCATE (drefob_dTob)
IF (ALLOCATED (drefob_dqob)) DEALLOCATE (drefob_dqob)
IF (ALLOCATED (dqob_dqb)) DEALLOCATE (dqob_dqb)
IF (ALLOCATED (dTob_dTb)) DEALLOCATE (dTob_dTb)
IF (ALLOCATED (dPob_dT)) DEALLOCATE (dPob_dT)
IF (ALLOCATED (dPob_dPb)) DEALLOCATE (dPob_dPb)
IF (ALLOCATED (dPob_dP)) DEALLOCATE (dPob_dP)
IF (ALLOCATED (dTob_dP)) DEALLOCATE (dTob_dP)
IF (ALLOCATED (dTb_dp)) DEALLOCATE (dTb_dp)

END SUBROUTINE jacobian_interface

end module ufo_gnssro_refmetoffice_tlad_mod
