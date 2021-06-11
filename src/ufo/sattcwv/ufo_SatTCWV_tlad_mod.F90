!-------------------------------------------------------------------------------
! (C) Crown Copyright 2021 Met Office
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
!-------------------------------------------------------------------------------
!> Fortran module for precipitable water tangent linear and adjoint

module ufo_sattcwv_tlad_mod
use iso_c_binding

use kinds
use oops_variables_mod
use ufo_vars_mod
use ufo_geovals_mod
use ufo_geovals_mod_c,   only: ufo_geovals_registry
use ufo_basis_tlad_mod,  only: ufo_basis_tlad
use obsspace_mod
use fckit_log_module, only : fckit_log

! Generic routines from elsewhere in jedi
use missing_values_mod
use ufo_constants_mod, only: one, zero, half, grav ! Gravitational field strength 

integer, parameter         :: max_string=800

!> Fortran derived type for precipitable water trajectory
type, extends(ufo_basis_tlad)   ::  ufo_sattcwv_tlad
  private
  integer                       :: nlevp, nlevq, nlocs
  real(kind_real), allocatable  :: K(:,:)

  contains
    procedure :: setup      => ufo_sattcwv_setup
    procedure :: delete     => ufo_sattcwv_tlad_delete
    procedure :: settraj    => ufo_sattcwv_tlad_settraj
    procedure :: simobs_tl  => ufo_sattcwv_simobs_tl
    procedure :: simobs_ad  => ufo_sattcwv_simobs_ad
end type ufo_sattcwv_tlad

contains

! ------------------------------------------------------------------------------
! Get the optional settings for the forward model, and save them in the object
! so that they can be used in the code.
! ------------------------------------------------------------------------------
subroutine ufo_sattcwv_setup(self, f_conf)

use fckit_configuration_module, only: fckit_configuration
implicit none
class(ufo_sattcwv_tlad), intent(inout) :: self
type(fckit_configuration), intent(in)    :: f_conf

end subroutine ufo_sattcwv_setup

!-------------------------------------------------------------------------
! Tidy up the variables that are used for passing information
!-------------------------------------------------------------------------
subroutine ufo_sattcwv_tlad_delete(self)

  implicit none
  class(ufo_sattcwv_tlad), intent(inout) :: self
  character(len=*), parameter :: myname_="ufo_sattcwv_tlad_delete"

  self%nlocs = 0
  self%nlevP = 0
  self%nlevq = 0
  if (allocated(self%K)) deallocate(self%K)
  self%ltraj = .false.

end subroutine ufo_sattcwv_tlad_delete

! ------------------------------------------------------------------------------
! Calculate the K-matrix (Jacobian) for the observation. It is necessary to run
! this routine before calling the TL or AD routines.
! ------------------------------------------------------------------------------
subroutine ufo_sattcwv_tlad_settraj(self, geovals, obss)
!
  implicit none
! Subroutine arguments
  class(ufo_sattcwv_tlad), intent(inout) :: self      ! The object that we use to save data in
  type(ufo_geovals),       intent(in)    :: geovals   ! The input geovals
  type(c_ptr), value,      intent(in)    :: obss      ! The input observations

! Local parameters
  character(len=*), parameter :: myname_="ufo_sattcwv_tlad_settraj"

! Local variables
  type(ufo_geoval), pointer    :: prs          ! The model geovals - atmospheric pressure (Pa)
  type(ufo_geoval), pointer    :: psfc         ! The model geovals - surface pressure (Pa)
  type(ufo_geoval), pointer    :: q            ! Model background values of specific humidity (kg/kg)
  real(kind_real), allocatable :: dtcwv_dq(:)  ! The deriv of precip water vs layer q conc (kg/sq.m/kg/kg)
  logical :: ascend                         ! Flag on direction of model levels
  integer :: iobs                           ! Counter
  integer :: nlocs                          ! number of observations
  integer :: ilev1                          ! starting level for loop
  integer :: ilev2                          ! ending level for loop
  integer :: inc                            ! increment to level loop 
  integer :: ibot                           ! index of second lowest level
  integer :: isfc                           ! index of lowest level
  integer, parameter      :: max_string = 800
  character(max_string)   :: message      ! General message for output
  character(max_string)   :: err_msg      ! Error Messages to be output to the user

  write(err_msg,*) "TRACE: ufo_sattcwv_tlad_settraj: begin"
  call fckit_log%info(err_msg)

! Make sure that any previous values of geovals don't get carried over
  call self%delete()

! get model state variables from geovals
  call ufo_geovals_get_var(geovals, var_prsi, prs)  ! pressure (Pa)
  call ufo_geovals_get_var(geovals, var_ps, psfc)   ! Surface pressure (Pa)
  call ufo_geovals_get_var(geovals, var_q, q)       ! Specific humidity (kg/kg)

! Keep copy of dimensions for TL & AD
  self % nlevP = prs % nval
  self % nlevq = q % nval
  self % nlocs = obsspace_get_nlocs(obss)
!
! The model data must be on a staggered grid, with nlevp = nlevq+1
  IF (prs % nval /= q % nval + 1) THEN
    write(err_msg,*) myname_ // ':' // &
    ' Data must be on a staggered grid nlevp, nlevq = ',prs%nval,q%nval
    call fckit_log % warning(err_msg)
    write(err_msg,*) myname_ // ':' // ' error: number of levels inconsistent!'
    call abor1_ftn(err_msg)
  END IF
!
! Determine if models levels are ascending or descending and set up indices for level loop
! Note assumes model levels stay same for all obs
!
  ascend = .false.
  if((prs%vals(1,1)-prs%vals(2,1)) > zero )ascend = .true.
  if (ascend)then       ! Model level starts above surface
    ilev1 = 1
    ilev2 = q%nval
    inc = 1 
    ibot = 2
    isfc = 1
  else                  ! Model level starts at top of atmosphere
    ilev1 = q%nval
    ilev2 = 1
    inc = -1
    ibot = q%nval
    isfc = q%nval
  endif
!
  ALLOCATE(self % K(1:self%nlocs, 1:self%nlevq))
  ALLOCATE(dtcwv_dq(1:self%nlevq))
!
! For each observation, calculate the K-matrix
  obs_loop: do iobs = 1, self % nlocs
      CALL sattcwv_GetK(  self % nlevP,        &  ! Number of pressure levels
                          ilev1,               &  ! Starting pressure layer
                          ilev2,               &  ! Ending pressure layer
                          inc,                 &  ! Layer increment
                          ibot,                &  ! Lowest layer inc surface
                          isfc,                &  ! Lowest level next to surface
                          psfc % vals(1,iobs), &  ! Surface pressure (Pa)
                          prs % vals(:,iobs),  &  ! Pressure levels (Pa)
                          dtcwv_dq)               ! deriv of precip water vs layer q concentration
!
! Build K-matrix (Jacobian of precip water with respect to layer q)
      self%K(iobs,1:self%nlevq) = dtcwv_dq(1:self%nlevq) 
 !
  end do obs_loop

! Note that this routine has been run.
  self%ltraj = .true.
  deallocate(dtcwv_dq)

end subroutine ufo_sattcwv_tlad_settraj

! ------------------------------------------------------------------------------
! Given an increment to the model state, calculate an increment to the
! observation
! ------------------------------------------------------------------------------
subroutine ufo_sattcwv_simobs_tl(self, geovals, hofx, obss)
!
  implicit none
!
! Subroutine arguments
  class(ufo_sattcwv_tlad), intent(in)     :: self      ! Object which is being used to transfer information
  type(ufo_geovals),       intent(in)     :: geovals   ! Model perturbations
  real(kind_real),         intent(inout)  :: hofx(:)   ! Increment to the observations
  type(c_ptr),   value,    intent(in)     :: obss      ! Input - the observations

! Local parameters
  character(len=*), parameter  :: myname_="ufo_sattcwv_simobs_tl"

! Local variables
  integer                      :: iobs      ! Loop variable, observation number
  integer                      :: ilev      ! loop variable, level number
  integer                      :: nlev      ! number of levels
  integer                      :: nlocs     ! Number of observations
  character(max_string)        :: err_msg   ! Message to be output
  character(max_string)        :: message   ! General message for output
  type(ufo_geoval), pointer    :: q_d       ! Increment to the specific humidity
  real(kind_real), allocatable :: x_d(:)    ! Increment to the state vector

  write(err_msg,*) "TRACE: ufo_sattcwv_simobs_tl: begin"
  call fckit_log%info(err_msg)

! Check if trajectory was set
  if (.not. self%ltraj) then
     write(err_msg,*) myname_, ' trajectory wasnt set!'
     call abor1_ftn(err_msg)
  endif

! Check if nlocs is consistent in geovals & hofx
  if (geovals%nlocs /= size(hofx)) then
     write(err_msg,*) myname_, ' error: nlocs ',&
       geovals%nlocs,' inconsistent with HofX size',size(hofx)
     call abor1_ftn(err_msg)
  endif

! Get perturbed variables from geovals q (kg/kg)
!
  call ufo_geovals_get_var(geovals, var_q, q_d)  ! specific humidity

  nlev = self % nlevq  ! number of layers
  nlocs = self % nlocs ! number of observations
  hofx = zero
! Allocate state vector increment
  allocate(x_d(1:nlev))
  x_d = zero
!
! Loop through the obs, calculating the increment to the observation hofx
!  
  obs_loop: do iobs = 1, nlocs
    lev_loop: do ilev = 1 , nlev
      x_d(ilev) = q_d % vals(ilev,iobs)
      hofx(iobs) = hofx(iobs) + (self % K(iobs,ilev) * x_d(ilev) )
    end do lev_loop
  end do obs_loop
!
  deallocate(x_d)
!
  write(err_msg,*) "TRACE: ufo_sattcwv_simobs_tl: complete"
  call fckit_log%info(err_msg)
!
  return
end subroutine ufo_sattcwv_simobs_tl

! ------------------------------------------------------------------------------
! Given an increment to the observation, find the equivalent increment to the
! model state
! ------------------------------------------------------------------------------
subroutine ufo_sattcwv_simobs_ad(self, geovals, hofx, obss)
!
  use typesizes,     only: wp => EightByteReal
!
  implicit none
!
! Subroutine arguments
  class(ufo_sattcwv_tlad), intent(in)     :: self     ! Object which is being used to transfer information
  type(ufo_geovals),       intent(inout)  :: geovals  ! Calculated perturbations to model state
  real(kind_real),         intent(in)     :: hofx(:)  ! Increment to the observations
  type(c_ptr),  value,     intent(in)     :: obss     ! Input - the observations

! Local parameters
  character(len=*), parameter     :: myname_="ufo_sattcwv_simobs_ad"

! Local variables
  real(c_double)               :: missing  ! Missing data values
  type(ufo_geoval), pointer    :: q_d      ! Pointer to the specific humidity perturbations
  integer                      :: iobs     ! Loop variable, observation number
  integer                      :: nlocs    ! Number of observations
  integer                      :: ilev     ! Loop variable, level number
  integer                      :: nlev     ! level number
  real(kind_real), allocatable :: x_d(:)   ! Perturbation to the state vector
  character(max_string)        :: err_msg  ! Message to be output

  write(err_msg,*) "TRACE: ufo_sattcwv_simobs_ad: begin"
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

! setting up array for perturbed variables in geovals 
  call ufo_geovals_get_var(geovals, var_q,  q_d)         ! layer specific humidity
!
  nlev =  self % nlevq  ! Number of layers
  nlocs = self % nlocs  ! Number of observations

! Allocate the output for the specific humidity perturbations
  if (.not. allocated(q_d%vals)) then
      q_d % nlocs = self % nlocs
      q_d % nval = self % nlevq
      allocate(q_d % vals(q_d % nval, q_d % nlocs))
      q_d % vals = zero
  endif
!
  missing = missing_value(missing)
! Allocate state vector x_d
  allocate(x_d(1:nlev))
!
! Loop through the obs, calculating the increment to the model state layer q
!  
  obs_loop: do iobs = 1, nlocs
    if (hofx(iobs) /= missing) then
      level_loop: do ilev = 1 , nlev
          x_d(ilev) =  hofx(iobs) * self % K(iobs,ilev)
          q_d% vals(ilev,iobs) =  x_d(ilev)
      end do level_loop
    end if
  end do obs_loop
!
  deallocate(x_d)
!
  write(err_msg,*) "TRACE: ufo_sattcwv_simobs_ad: complete"
  call fckit_log%info(err_msg)
!
  return
end subroutine ufo_sattcwv_simobs_ad

!------------------------------------------------------------
! Calculate the partial derivatives dprecipwater/dq 
!------------------------------------------------------------
SUBROUTINE sattcwv_GetK(  nlevP,         &
                          ilev1,         &
                          ilev2,         &
                          inc,           &
                          ibot,          &
                          isfc,          &
                          psfc,          &
                          prs,           &
                          dtcwv_dq)
!
IMPLICIT NONE

INTEGER, INTENT(IN)            :: nlevP         ! The number of model pressure levels
INTEGER, INTENT(IN)            :: ilev1         ! Starting layer 
INTEGER, INTENT(IN)            :: ilev2         ! Ending layer
INTEGER, INTENT(IN)            :: inc           ! loop increment
INTEGER, INTENT(IN)            :: ibot          ! lowest layer including surface
INTEGER, INTENT(IN)            :: isfc          ! lowest level next to surface
REAL(kind_real), INTENT(IN)    :: psfc          ! surface pressure (Pa)
REAL(kind_real), INTENT(IN)    :: prs(1:nlevP)  ! Model background pressure (Pa) at levels
REAL(kind_real), INTENT(INOUT) :: dtcwv_dq(1:nlevP-1)! Partial deriv of HofX for precip water wrt layer water vapour
!
! Local constants
! 
character(max_string) :: err_msg           ! Error message to be output
character(max_string) :: message           ! General message for output
character(len=*), parameter :: myname_ = "ufo_sattcwv_tlad"
REAL(kind_real)       :: PDiff             ! Pressure diff across layer
REAL(kind_real)       :: GK                ! 1/gravity
INTEGER               :: i                 ! level counter
!-------------------------------------------------------------------------------
!1. Initialise variables 
!-------------------------------------------------------------------------------
PDiff = zero
GK = one / grav
dtcwv_dq(:) = zero
!
!-------------------------------------------------------------------------------
! Calculate derivative dtcwv wrt dq at all levels for K matrix
!-------------------------------------------------------------------------------
!
DO i = ilev1, ilev2, inc
   PDiff = prs(i) - prs(i+1)     ! prs is pressure on level i
! Include surface layer assuming surface q is same as q 10m but could use q2m in future
   IF (i == isfc)PDiff = psfc - prs(ibot)
   dtcwv_dq(i) =  GK * PDiff     ! compute dtcwv/dq on each layer 
END DO
!
END SUBROUTINE sattcwv_GetK

end module ufo_sattcwv_tlad_mod

