!-------------------------------------------------------------------------------
! (C) British Crown Copyright 2020 Met Office
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
!-------------------------------------------------------------------------------
!> Fortran module for ground based GNSS Met Office's tangent linear and adjoint

module ufo_groundgnss_metoffice_tlad_mod
use iso_c_binding

use kinds
use ufo_vars_mod
use ufo_geovals_mod
use ufo_geovals_mod_c,   only: ufo_geovals_registry
use ufo_basis_tlad_mod,  only: ufo_basis_tlad
use obsspace_mod
use missing_values_mod
use fckit_log_module, only : fckit_log
use ufo_utils_refractivity_calculator, only: &
    ufo_calculate_refractivity, ufo_refractivity_kmat



use ufo_constants_mod, only: &
    rd,                      &    ! Gas constant for dry air
    grav,                    &    ! Gravitational field strength
    n_alpha                       ! Refractivity constant a


integer, parameter             :: max_string=800

!> Fortran derived type for groundgnss trajectory
type, extends(ufo_basis_tlad)  ::  ufo_groundgnss_metoffice_tlad
  private

  integer                      :: nlevp, nlevq, nlocs, iflip
  real(kind_real), allocatable :: K(:,:)
  real(kind_real), allocatable :: dztd_dp(:,:)
  real(kind_real), allocatable :: dztd_dq(:,:)
  logical :: vert_interp_ops
  logical :: pseudo_ops
  real(kind_real) :: min_temp_grad
  contains
    procedure :: setup      => ufo_groundgnss_metoffice_setup
    procedure :: delete     => ufo_groundgnss_metoffice_tlad_delete
    procedure :: settraj    => ufo_groundgnss_metoffice_tlad_settraj
    procedure :: simobs_tl  => ufo_groundgnss_metoffice_simobs_tl
    procedure :: simobs_ad  => ufo_groundgnss_metoffice_simobs_ad
end type ufo_groundgnss_metoffice_tlad

contains

! ------------------------------------------------------------------------------
! Get the optional settings for the forward model, and save them in the object
! so that they can be used in the code.
! ------------------------------------------------------------------------------
subroutine ufo_groundgnss_metoffice_setup(self, f_conf)

use fckit_configuration_module, only: fckit_configuration
implicit none
class(ufo_groundgnss_metoffice_tlad), intent(inout) :: self
type(fckit_configuration), intent(in)               :: f_conf

call f_conf%get_or_die("vert_interp_ops", self % vert_interp_ops)
call f_conf%get_or_die("pseudo_ops", self % pseudo_ops)
call f_conf%get_or_die("min_temp_grad", self % min_temp_grad)

end subroutine ufo_groundgnss_metoffice_setup


! ------------------------------------------------------------------------------
! Calculate the K-matrix (Jacobian) for the observation.  It is necessary to run
! this routine before calling the TL or AD routines.
! ------------------------------------------------------------------------------
subroutine ufo_groundgnss_metoffice_tlad_settraj(self, geovals, obss)

  use fckit_exception_module, only: fckit_exception

  implicit none
! Subroutine arguments
  class(ufo_groundgnss_metoffice_tlad), intent(inout) :: self      ! The object that we use to save data in
  type(ufo_geovals),                    intent(in)    :: geovals   ! The input geovals
  type(c_ptr), value,                   intent(in)    :: obss      ! The input observations

! Local parameters
  character(len=*), parameter :: myname_="ufo_groundgnss_metoffice_tlad_settraj"

! Local variables
  type(ufo_geoval), pointer    :: q                  ! The model geovals - specific humidity
  type(ufo_geoval), pointer    :: prs                ! The model geovals - atmospheric pressure
  type(ufo_geoval), pointer    :: rho_heights        ! The model geovals - heights of the pressure-levels
  type(ufo_geoval), pointer    :: theta_heights      ! The model geovals - heights of the theta-levels (q on theta)
  integer                      :: nstate             ! The size of the state vector
  integer                      :: iobs               ! Loop variable, observation number
  integer                      :: nobs               ! Number of observations
  integer                      :: ilev               ! Loop variable, vertical level number
  real(kind_real), allocatable :: zStation(:)        ! The station height
  real(kind_real), allocatable :: pressure(:)        ! Model background values of air pressure (monotonic order)
  real(kind_real), allocatable :: humidity(:)        ! Model background specific humidity  (in pressure monotonic order)
  real(kind_real), allocatable :: za(:)              ! Model heights of rho levs (in pressure monotonic order)
  real(kind_real), allocatable :: zb(:)              ! Model heights of theta levs (in pressure monotonic order)

  call fckit_log%info("TRACE: ufo_groundgnss_metoffice_tlad_settraj: begin")

! Make sure that any previous values of geovals don't get carried over
  call self%delete()

! get model state variables from geovals
  call ufo_geovals_get_var(geovals, var_q, q)               ! specific humidity
  call ufo_geovals_get_var(geovals, var_prsi, prs)          ! pressure
  call ufo_geovals_get_var(geovals, var_z, theta_heights)   ! Geopotential height of the normal model levels
  call ufo_geovals_get_var(geovals, var_zi, rho_heights)    ! Geopotential height of the pressure levels

! Check that the geovals are ordered top to bottom
  if( prs%vals(1,1) > prs%vals(prs%nval,1) ) then
    call fckit_exception%throw('Geovals should be ordered top to bottom')
  endif

! Keep copy of dimensions
  self % nlevp = prs % nval
  self % nlevq = q % nval
  self % nlocs = obsspace_get_nlocs(obss)

! Get the meta-data from the observations
  nobs  = obsspace_get_nlocs(obss)
  allocate(zStation(nobs))

  call obsspace_get_db(obss, "MetaData", "stationElevation", zStation)

  nstate = prs % nval + q % nval
  ALLOCATE(self % K(1:self%nlocs, 1:nstate))
  ALLOCATE(pressure(1:self%nlevp))
  ALLOCATE(humidity(1:self%nlevq))
  ALLOCATE(za(1:self%nlevp))
  ALLOCATE(zb(1:self%nlevq))

! For each observation, calculate the K-matrix
  obs_loop: do iobs = 1, self % nlocs

      pressure = prs % vals(:,iobs)
      humidity = q % vals(:,iobs)
      za = rho_heights % vals(:,iobs)
      zb = theta_heights % vals(:,iobs)

      CALL groundgnss_jacobian_interface(self % nlevp,                 &   ! Number of pressure levels
                                         self % nlevq,                 &   ! Number of specific humidity levels
                                         za(1:self%nlevp),             &   ! Heights of the pressure levels
                                         zb(1:self%nlevq),             &   ! Heights of the specific humidity levels
                                         humidity(1:self%nlevq),       &   ! Values of the specific humidity
                                         pressure(1:self%nlevp),       &   ! Values of the pressure
                                         zStation(iobs),               &   ! Station height
                                         iobs,                         &   ! Ob number
                                         self % vert_interp_ops,       &   ! Pressure varies exponentially with height?
                                         self % pseudo_ops,            &   ! Use pseudo-levels in calculation?
                                         self % min_temp_grad,         &   ! Minimum temperature gradient allowed
                                         self % K(:, 1:nstate))            ! K-matrix (Jacobian of the observation with respect to the inputs

  end do obs_loop

! Note that this routine has been run.
  self%ltraj = .true.

  deallocate(zStation)

end subroutine ufo_groundgnss_metoffice_tlad_settraj


! ------------------------------------------------------------------------------
! Given an increment to the model state, calculate an increment to the
! observation
! ------------------------------------------------------------------------------
subroutine ufo_groundgnss_metoffice_simobs_tl(self, geovals, hofx, obss)

  implicit none

! Subroutine arguments
  class(ufo_groundgnss_metoffice_tlad), intent(in) :: self      ! Object which is being used to transfer information
  type(ufo_geovals),                intent(in)     :: geovals   ! Model perturbations
  real(kind_real),                  intent(inout)  :: hofx(:)   ! Increment to the observations
  type(c_ptr),   value,             intent(in)     :: obss      ! Input - the observations

! Local parameters
  character(len=*), parameter  :: myname_="ufo_groundgnss_metoffice_simobs_tl"

! Local variables
  integer                      :: iobs          ! Loop variable, observation number
  integer                      :: nlocs         ! Number of observations
  integer                      :: ilev          ! Loop variable, pressure level number
  integer                      :: iflip         ! Index for vertical flip
  type(ufo_geoval), pointer    :: q_d           ! Increment to the specific humidity
  type(ufo_geoval), pointer    :: prs_d         ! Increment to the air pressure
  real(kind_real), allocatable :: pressure_d(:) ! Increment to the air pressure (monotonic order)
  real(kind_real), allocatable :: humidity_d(:) ! Increment to the specific humidity  (in pressure monotonic order)
  real(kind_real), allocatable :: x_d(:)        ! Increment to the complete state

  call fckit_log%info("TRACE: ufo_groundgnss_metoffice_simobs_tl: begin")

! Check if trajectory was set
  if (.not. self%ltraj) then
     call abor1_ftn(trim(myname_) // ' trajectory wasnt set!')
  endif

! Check if nlocs is consistent in geovals & hofx
  if (geovals%nlocs /= size(hofx)) then
     call abor1_ftn(trim(myname_) // ' error: nlocs inconsistent!')
  endif

! Get variables from geovals
  call ufo_geovals_get_var(geovals, var_q,     q_d)         ! specific humidity
  call ufo_geovals_get_var(geovals, var_prsi,  prs_d)       ! pressure on rho levels

  nlocs = self % nlocs ! number of observations

  allocate(x_d(1:prs_d%nval+q_d%nval))
  allocate(pressure_d(1:self % nlevp))
  allocate(humidity_d(1:self % nlevq))


! Loop through the obs, calculating the increment to the observation
  obs_loop: do iobs = 1, nlocs   ! order of loop doesn't matter

    pressure_d(1:self % nlevp) = prs_d % vals(:,iobs)
    humidity_d(1:self % nlevq) = q_d % vals(:,iobs)

    x_d(1:prs_d%nval) = pressure_d
    x_d(prs_d%nval+1:prs_d%nval+q_d%nval) = humidity_d
    hofx(iobs) = SUM(self % K(iobs,:) * x_d)

  end do obs_loop

  deallocate(x_d)
  deallocate(pressure_d)
  deallocate(humidity_d)

  call fckit_log%info("TRACE: ufo_groundgnss_metoffice_simobs_tl: complete")

  return

end subroutine ufo_groundgnss_metoffice_simobs_tl



! ------------------------------------------------------------------------------
! Given an increment to the observation, find the equivalent increment to the
! model state
! ------------------------------------------------------------------------------
subroutine ufo_groundgnss_metoffice_simobs_ad(self, geovals, hofx, obss)

  use typesizes,     only: wp => EightByteReal

  implicit none

! Subroutine arguments
  class(ufo_groundgnss_metoffice_tlad), intent(in) :: self      ! Object which is being used to transfer information
  type(ufo_geovals),                intent(inout)  :: geovals   ! Calculated perturbations to model state
  real(kind_real),                  intent(in)     :: hofx(:)   ! Increment to the observations
  type(c_ptr),  value,              intent(in)     :: obss      ! Input - the observations

! Local parameters
  character(len=*), parameter     :: myname_="ufo_groundgnss_metoffice_simobs_ad"

! Local variables
  real(c_double)               :: missing       ! Missing data values
  type(ufo_geoval), pointer    :: q_d           ! Pointer to the specific humidity perturbations
  type(ufo_geoval), pointer    :: prs_d         ! Pointer to the pressure perturbations
  integer                      :: iobs          ! Loop variable, observation number
  integer                      :: ilev          ! Loop variable, pressure level number
  integer                      :: iflip         ! Index for vertical flip
  real(kind_real), allocatable :: x_d(:)        ! Perturbation to the full model state
  real(kind_real), allocatable :: pressure_d(:) ! Perturbation to pressure (monotonic order)
  real(kind_real), allocatable :: humidity_d(:) ! Perturbation to specific humidity  (in pressure monotonic order)

  call fckit_log%info("TRACE: ufo_groundgnss_metoffice_simobs_ad: begin")

! Check if trajectory was set
  if (.not. self%ltraj) then
     call abor1_ftn(trim(myname_) // ' trajectory wasnt set!')
  endif

! Check if nlocs is consistent in geovals & hofx
  if (geovals%nlocs /= size(hofx)) then
     call abor1_ftn(trim(myname_) // ' error: nlocs inconsistent!')
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

  call fckit_log%info("TRACE: ufo_groundgnss_metoffice_simobs_ad: complete")

  return

end subroutine ufo_groundgnss_metoffice_simobs_ad



!-------------------------------------------------------------------------
! Tidy up the variables that are used for passing information
!-------------------------------------------------------------------------
subroutine ufo_groundgnss_metoffice_tlad_delete(self)

  implicit none
  class(ufo_groundgnss_metoffice_tlad), intent(inout) :: self
  character(len=*), parameter :: myname_="ufo_groundgnss_metoffice_tlad_delete"

  self%nlocs = 0
  self%nlevp = 0
  self%nlevq = 0
  if (allocated(self%K)) deallocate(self%K)
  self%ltraj = .false.

end subroutine ufo_groundgnss_metoffice_tlad_delete


!-------------------------------------------------------------------------
! Interface for calculating the K-matrix for calculating TL/AD
!-------------------------------------------------------------------------
SUBROUTINE groundgnss_jacobian_interface(nlevp,     &
                              nlevq,                &
                              za,                   &
                              zb,                   &
                              q,                    &
                              prs,                  &
                              zStation,             &
                              iobs,                 &
                              vert_interp_ops,      &
                              pseudo_ops,           &
                              gbgnss_min_temp_grad, &
                              K)

IMPLICIT NONE

INTEGER, INTENT(IN)            :: nlevP                 ! The number of model pressure levels
INTEGER, INTENT(IN)            :: nlevq                 ! The number of model theta levels
REAL(kind_real), INTENT(IN)    :: za(:)                 ! The geometric height of the model pressure levels
REAL(kind_real), INTENT(IN)    :: zb(:)                 ! The geometric height of the model theta levels
REAL(kind_real), INTENT(IN)    :: q(1:nlevq)            ! The model values that are being perturbed
REAL(kind_real), INTENT(IN)    :: prs(1:nlevP)          ! The model values that are being perturbed
REAL(kind_real), INTENT(IN)    :: zStation              ! Station height
LOGICAL, INTENT(IN)            :: vert_interp_ops       ! Pressure varies exponentially with height?
LOGICAL, INTENT(IN)            :: pseudo_ops            ! Use pseudo-levels in calculation?
REAL(kind_real), INTENT(IN)    :: gbgnss_min_temp_grad  ! The minimum temperature gradient which is used
INTEGER, INTENT(IN)            :: iobs                  ! Ob number

REAL(kind_real), INTENT(INOUT) :: K(:,:)                ! The calculated K matrix
!
! Things that may need to be output, as they are used by the TL/AD calculation
!

REAL(kind_real)                :: T(1:nlevq)            ! Temperature on model levels
REAL(kind_real), ALLOCATABLE   :: refrac(:)             ! model refractivity on theta levels
!
! Local variables
!
INTEGER                      :: nstate             ! Number of levels in state vector
REAL(kind_real)              :: x(1:nlevp+nlevq)   ! state vector
LOGICAL                      :: refracerr          ! Whether we encountered an error in calculating the refractivity

REAL(kind_real)              :: pN(nlevq)          ! Presure on theta levels

REAL(kind_real), ALLOCATABLE :: model_heights(:)   ! Geopotential heights of the refractivity levels (not needed for this oper)
INTEGER                      :: nRefLevels         ! Number of levels in refractivity calculation

! Set up the size of the state
nstate = nlevP + nlevq

refracerr = .FALSE.

! Calculate the refractivity
! The ufo_utils_refractivity_calculator currently relies on the variables being provided
! bottom to top. The geovals are provided top to bottom and so the vertical variables need
! to be reveresed when they are passed into this routine and the outputs then need to be reversed
CALL ufo_calculate_refractivity(nlevp,                &
                                nlevq,                &
                                za(nlevp:1:-1),       &
                                zb(nlevq:1:-1),       &
                                prs(nlevp:1:-1),      &
                                q(nlevq:1:-1),        &
                                vert_interp_ops,      &
                                pseudo_ops,           &
                                gbgnss_min_temp_grad, &
                                refracerr,            &
                                nRefLevels,           &
                                refrac,               &
                                model_heights)

! Flip the vertical direction of refrac and model_height back to top to bottom.
refrac = refrac(nRefLevels:1:-1)
model_heights = model_heights(nRefLevels:1:-1)

IF (.NOT. refracerr) THEN
    ! Calculate the K-matrix (Jacobian)
    CALL Groundgnss_GetK(nstate,               &
                         nlevP,                &
                         nlevq,                &
                         za,                   &
                         zb,                   &
                         prs,                  &
                         q,                    &
                         zStation,             &
                         iobs,                 &
                         vert_interp_ops,      &
                         pseudo_ops,           &
                         gbgnss_min_temp_grad, &
                         refracerr,            &
                         refrac,               &
                         K)
ELSE
    K = 0
    CALL fckit_log % warning("Error in refractivity calculation")
END IF


END SUBROUTINE groundgnss_jacobian_interface



!-------------------------------------------------------------------------
! Calculate the K-matrix (Jacobian)
!-------------------------------------------------------------------------
SUBROUTINE Groundgnss_GetK(nstate,              &
                          nlevP,                &
                          nlevq,                &
                          za,                   &
                          zb,                   &
                          P,                    &
                          q,                    &
                          zStation,             &
                          iobs,                 &
                          vert_interp_ops,      &
                          pseudo_ops,           &
                          gbgnss_min_temp_grad, &
                          refracerr,            &
                          refrac,               &
                          K)

!
! Return the K-matrix for calculating TL/AD
!
IMPLICIT NONE

INTEGER, INTENT(IN)             :: nstate
INTEGER, INTENT(IN)             :: nlevP                  ! The number of model pressure levels
INTEGER, INTENT(IN)             :: nlevq                  ! The number of model theta levels
REAL(kind_real), INTENT(IN)     :: za(:)                  ! The geometric height of the model pressure levels
REAL(kind_real), INTENT(IN)     :: zb(:)                  ! The geometric height of the model theta levels
REAL(kind_real), INTENT(IN)     :: P(:)                   ! The model pressure values
REAL(kind_real), INTENT(IN)     :: q(:)                   ! The model humidity values
REAL(kind_real), INTENT(IN)     :: zStation               ! Station height
INTEGER, INTENT(IN)             :: iobs                   ! Ob number
LOGICAL, INTENT(IN)             :: vert_interp_ops        ! Pressure varies exponentially with height?
LOGICAL, INTENT(IN)             :: pseudo_ops             ! Use pseudo-levels in calculation?
REAL(kind_real), INTENT(IN)     :: gbgnss_min_temp_grad   ! The minimum temperature gradient which is used
LOGICAL, INTENT(INOUT)          :: refracerr              ! Whether we encountered an error in calculating the refractivity
REAL(kind_real), INTENT(IN)     :: refrac(:)              ! Model refractivity on theta levels - returned from forward model
REAL(kind_real), INTENT(INOUT)  :: K(:,:)                 ! The calculated K matrix

REAL(kind_real), ALLOCATABLE    :: dref_dP(:, :)          ! Partial derivative of refractivity wrt. pressure
REAL(kind_real), ALLOCATABLE    :: dref_dq(:, :)          ! Partial derivative of refractivity wrt. specific humidity


! Local constants

CHARACTER (LEN=*), PARAMETER :: RoutineName = "Groundgnss_GetK"
REAL, PARAMETER :: refrac_scale = 1.0E-6      ! Conversion factor between refractivity and refractive index

! Local variables

INTEGER           :: Level             ! Used for iteration over levels
INTEGER           :: Lowest_Level      ! Lowest height level
INTEGER           :: FirstNeg          ! First negative

REAL(kind_real)              :: p_local(nlevp)         ! pressure on rho levels (with no negative pressures)
REAL(kind_real)              :: pN(nlevq)              ! pressure on theta levels
REAL(kind_real)              :: LocalZenithDelay       ! Zenith Total Delay
REAL(kind_real)              :: h_diff, station_diff   ! Height diff, station height diff
REAL(kind_real)              :: c_rep                  ! 1/scale height
REAL(kind_real)              :: const, c, term1, term2 ! constant, scale height,
REAL(kind_real)              :: StationRefrac          ! Refraction at station height
REAL(kind_real)              :: z_weight1              ! For linear interpolation
REAL(kind_real)              :: z_weight2              ! For linear interpolation
REAL(kind_real)              :: dc_dref                ! derivative of scale height wrt refrac term1
REAL(kind_real)              :: dc_dref2               ! derivative of scale height wrt refrac term2
REAL(kind_real)              :: dztd_dc                ! derivative of ZTD wrt scale height
REAL(kind_real)              :: drefsta_dref           ! derivative of station refrac wrt refrac
REAL(kind_real)              :: drefsta_dc             ! derivative of station refrac wrt scale height
REAL(kind_real)              :: dztd_drefsta           ! derivative of ZTD wrt station refrac
REAL(kind_real), ALLOCATABLE :: dp_local_dPin(:,:)     ! derivative of pressure rho levels wrt pressure
REAL(kind_real), ALLOCATABLE :: dztd_dpN(:)            ! array for derivative w.r.t top theta level
REAL(kind_real)              :: dztd_dq(nlevq)         ! The calculated K matrix
REAL(kind_real)              :: dztd_dp(nlevP)         ! The calculated K matrix
REAL(kind_real), ALLOCATABLE :: dPb_dP(:,:)            ! derivative of pressure theta wrt pressure
REAL(kind_real), ALLOCATABLE :: dztd_dref(:)           ! derivative of ZTD wrt refrac
REAL(kind_real), ALLOCATABLE :: x1(:,:)                ! Matrix placeholder
REAL(kind_real), ALLOCATABLE :: x2(:)                  ! Matrix placeholder

REAL(kind_real), PARAMETER   :: PressScale = 6000.0    ! Pressure scale height

!----------------------------------------------------------------------

!-------------------------------------------------------
! 0. Initialise variables
!-------------------------------------------------------

ALLOCATE (dref_dP(nlevq,nlevp))
ALLOCATE (dref_dq(nlevq,nlevq))
ALLOCATE (dPb_dP(nlevq,nlevp))
ALLOCATE (dztd_dref(nlevq))
ALLOCATE (dztd_dpN(nlevq))
ALLOCATE (x2(nlevp))
ALLOCATE (x1(nlevq,nlevp))
ALLOCATE (dp_local_dPin(nlevp,nlevp))

! Initialise matrices

dref_dq(:,:)   = 0.0
dref_dP(:,:)   = 0.0
dztd_dref(:)   = 0.0
dztd_dpN(:)    = 0.0
dPb_dP(:,:)    = 0.0
dztd_dq(:)     = 0.0
dztd_dp(:)     = 0.0
x1(:,:)        = 0.0
x2(:)          = 0.0
dp_local_dpin(:,:) = 0.0
p_local(:)     = 0.0
pN(:)          = 0.0

LocalZenithDelay = 0.0
StationRefrac    = 0.0

! If negative pressures exist or pressure is greater at a level then
! the level below, replace these
! and any above with values calculated as an exponential decay (scale height
! = 6km) from the highest positive pressure.
! Start at the bottom level and increment the levels towards the top
FirstNeg = nlevp+1
DO Level=nlevp, 1, -1
  IF (Level==nlevp) THEN
    p_local(Level) = P(Level)
    dp_local_dPin(Level,Level) = 1.0
  ELSE
    IF ((P(Level) <= 0.0 .OR. (FirstNeg /= nlevp+1 .AND. Level < FirstNeg)) .OR. &
        (P(Level) > P(Level+1))) THEN

      IF (FirstNeg == nlevp+1) FirstNeg = Level

      p_local(Level) = P(FirstNeg+1) * EXP(-(za(Level) - za(FirstNeg+1)) / PressScale)
      dp_local_dPin(Level,FirstNeg+1) = p_local(Level) / P(FirstNeg+1)

    ELSE
      p_local(Level)             = P(Level)
      dp_local_dPin(Level,Level) = 1.0
    END IF
  END IF
END DO

! Calculate Pressure on theta
! Assume ln(p) linear with height

DO Level = nlevp-1, 1, -1

  z_weight1 = (za(Level) - zb(Level)) / (za(Level) - za(Level+1))
  z_weight2 = 1.0 - z_weight1

  pN(Level) = EXP(z_weight1 * LOG(p_local(Level+1)) + z_weight2 * LOG(p_local(Level)))

  dPb_dP(Level,Level+1) = pN(Level) * z_weight1 / p_local(Level+1)
  dPb_dP(Level,Level) = pN(Level) * z_weight2 / p_local(Level)

END DO

! Calculate the gradient of ref wrt p (on rho levels) and q (on theta levels)

! The ufo_utils_refractivity_calculator currently relies on the variables being provided
! bottom to top. The geovals are provided top to bottom and so the vertical variables need
! to be reveresed when they are passed into this routine and the outputs then need to be reversed
! The arrays need to flipped
call reverse_levels_in_matrix(nlevq, nlevp, dPb_dP)
CALL ufo_refractivity_kmat (nlevP,                &
                            nlevq,                &
                            nlevq,                &
                            za(nlevp:1:-1),       &
                            zb(nlevq:1:-1),       &
                            P(nlevp:1:-1),        &
                            q(nlevq:1:-1),        &
                            pseudo_ops,           &
                            vert_interp_ops,      &
                            gbgnss_min_temp_grad, &
                            dref_dP,              &
                            dref_dq,              &
                            dPb_dP)

call reverse_levels_in_matrix(nlevq, nlevp, dref_dP)
call reverse_levels_in_matrix(nlevq, nlevq, dref_dq)
call reverse_levels_in_matrix(nlevq, nlevp, dPb_dP)

! In Layer where station height lies, define lowest level required for
! iteration and integration

DO Level = nlevp-1, 1, -1
  IF (zb(Level) > zStation) THEN
    Lowest_Level = Level
    EXIT
  END IF
END DO

DO Level = 1, Lowest_Level
  IF (refrac(Level) <= 0.0) THEN

    CALL fckit_log % warning("Refractivity error. Refractivity < 0.0")

    RETURN

  END IF
END DO

!---------------------------
! 3. Calculate Zenith delays
!---------------------------

! Start at bottom of the model

DO Level = Lowest_Level, 1,  -1

  LocalZenithDelay = 0.0

  IF (Level == Lowest_Level .AND. Level /= nlevq) THEN

    ! If station lies above the lowest model level, interpolate refractivity
    ! to station height
    ! h_diff is expected to be a negative number
    h_diff             = zb(Level + 1) - zb(Level)
    station_diff       = zStation - zb(Level)
    c                  = (LOG (refrac(Level) / refrac(Level + 1))) / h_diff
    StationRefrac      = refrac(Level + 1) * EXP(-c * (zStation-zb(Level + 1)))
    const              = (-StationRefrac / c) * EXP(c * zStation)
    term1              = EXP(-c * (zb(Level)))
    term2              = EXP(-c * zStation)
    LocalZenithDelay   = refrac_scale * const * (term1 - term2)

    c_rep              = 1.0 / c
    dztd_dc            = (-refrac_scale * StationRefrac / c) * &
                         (c_rep + EXP(c * station_diff) * (station_diff - c_rep))
    dztd_dc            = dztd_dc - (station_diff * LocalZenithDelay)
    dc_dref            = -1.0 / (refrac(Level + 1) * h_diff)
    dc_dref2           = 1.0 / (refrac(Level) * h_diff)
    drefsta_dref       = StationRefrac / refrac(Level + 1)
    drefsta_dc         = -(zStation - zb(Level + 1)) * StationRefrac
    drefsta_dref       = drefsta_dref + drefsta_dc * dc_dref
    dztd_drefsta       = LocalZenithDelay / StationRefrac

    dztd_dref(Level+1) = dztd_drefsta * drefsta_dref + dztd_dc * dc_dref
    dztd_dref(Level)   = dztd_dc * dc_dref2

  ELSE IF (Level == nlevq) THEN

    ! If station lies below the bottom level (ie. the lowest level for which refractivity is
    ! calculated), then use c from the first full layer, but integrate down to height of
    ! station

    h_diff            = zb(Level) - zb(Level - 1)
    c                 = (LOG (refrac(Level - 1) / refrac(Level))) / h_diff
    const             = (-refrac(Level) / c) * EXP(c * (zb(Level)))
    term1             = EXP(-c * (zb(Level - 1)))
    term2             = EXP(-c * zStation)
    LocalZenithDelay  = refrac_scale * const * (term1 - term2)

    c_rep             = 1.0 / c
    dztd_dc           = (-refrac_scale * refrac(Level) / c) * &
                        (c_rep + EXP(c * h_diff) * (h_diff - c_rep))
    dc_dref           = -1.0 / (refrac(Level) * h_diff)
    dc_dref2          = 1.0 / (refrac(Level - 1) * h_diff)

    dztd_dref(Level)  = LocalZenithDelay / refrac(Level)
    dztd_dref(Level)  = dztd_dref(Level) + dztd_dc * dc_dref
    dztd_dref(Level-1)= dztd_dc * dc_dref2

  ELSE IF (Level >= 1 .AND. Level < nlevq-1 .AND. Level < Lowest_Level) THEN

    ! For the other Levels in the column above the station.

    h_diff            = zb(Level + 1) - zb(Level)
    c                 = (LOG (refrac(Level) / refrac(Level+1))) / h_diff
    const             = (-refrac(Level + 1) / c) * EXP(c * (zb(Level + 1)))
    term1             = EXP(-c * (zb(Level)))
    term2             = EXP(-c * (zb(Level + 1)))
    LocalZenithDelay  = refrac_scale * const * (term1 - term2)

    c_rep             = 1.0 / c
    dztd_dc           = (-refrac_scale * refrac(Level + 1) / c) * (c_rep + EXP(c * h_diff) * (h_diff - c_rep))
    dc_dref           = -1.0 / (refrac(Level + 1) * h_diff)
    dc_dref2          = 1.0 / (refrac(Level) * h_diff)

    dztd_dref(Level+1) = dztd_dref(Level + 1) + LocalZenithDelay / refrac(Level + 1) + dztd_dc * dc_dref
    dztd_dref(Level)   = dztd_dc * dc_dref2

  END IF

END DO

!-------------------------------------------
! Construct K Matrices
!-------------------------------------------

! Account for negative pressure
dref_dP = MATMUL(dref_dP,dp_local_dPin)

dztd_dq(:) = MATMUL(dztd_dref,dref_dq)
dztd_dp(:) = MATMUL(dztd_dref,dref_dP)

! Account for negative pressure
dztd_dp(:) = MATMUL(dztd_dp,dp_local_dPin)

! First add in dZTD/dp for the top correction, which only depends on top level theta pressure

dztd_dpN(1) = refrac_scale * n_alpha * rd / grav
x1 = MATMUL(dPb_dP, dp_local_dPin)
x2 = MATMUL(dztd_dpN, x1)
dztd_dp = x2 + dztd_dp

K(iobs, 1:nlevp)  = dztd_dp
K(iobs, nlevp+1:nstate) = dztd_dq

DEALLOCATE (dp_local_dPin)
DEALLOCATE (x1)
DEALLOCATE (x2)
DEALLOCATE (dztd_dpN)
DEALLOCATE (dztd_dref)
DEALLOCATE (dPb_dP)
DEALLOCATE (dref_dq)
DEALLOCATE (dref_dP)

END SUBROUTINE Groundgnss_GetK


subroutine reverse_levels_in_matrix(sizea, sizeb, array)

implicit none

integer, intent(in) :: sizea
integer, intent(in) :: sizeb
real(kind_real), intent(inout) :: array(1:sizea, 1:sizeb)

real(kind_real) :: temp(1:sizea, 1:sizeb)
integer :: i

temp = array
do i = 1, sizea
  array(i,:) = temp(sizea-i+1, sizeb:1:-1)
end do

end subroutine reverse_levels_in_matrix

end module ufo_groundgnss_metoffice_tlad_mod

