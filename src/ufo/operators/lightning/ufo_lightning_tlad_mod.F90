! (C) Copyright 2017-2018 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

!> Fortran module for lightning tl/ad observation operator
!

module ufo_lightning_tlad_mod

use obs_variables_mod
use oops_variables_mod
use ufo_vars_mod
use ufo_geovals_mod
use kinds
use missing_values_mod
use ufo_constants_mod, only: one, zero, half, grav ! Gravitational field strength

 !> Fortran derived type for the tl/ad observation operator
 type, public :: ufo_lightning_tlad
 private
  type(obs_variables), public  :: obsvars
  type(oops_variables), public :: geovars

  integer                       :: nlevdp, nlevq, nlocs
  real(kind_real), allocatable  :: K(:,:)
!  logical, public     ::  l_fed_nonlinear
  integer, public     ::  n_horiz, i_horiz

 contains
  procedure :: cleanup     => ufo_lightning_tlad_cleanup_
  procedure :: settraj => ufo_lightning_tlad_settraj_
  procedure :: simobs_tl  => ufo_lightning_simobs_tl_
  procedure :: simobs_ad  => ufo_lightning_simobs_ad_
  final :: destructor
 end type ufo_lightning_tlad

 integer, parameter         :: max_string=800

contains

! ------------------------------------------------------------------------------
! This code is for the TL/AD lightning operator.
subroutine ufo_lightning_tlad_settraj_(self, geovals, obss)
use iso_c_binding
use obsspace_mod
implicit none
class(ufo_lightning_tlad), intent(inout) :: self
type(ufo_geovals),         intent(in)    :: geovals
type(c_ptr), value,        intent(in)    :: obss
!type(ufo_geovals),        intent(inout) :: hofxdiags    !non-h(x) diagnostics

! Local parameters
character(len=*), parameter :: myname_="ufo_lightning_tlad_settraj"

! Local variables
  type(ufo_geoval), pointer      :: delp                ! The model geovals - atmospheric pressure (Pa)
  type(ufo_geoval), pointer      :: qg                  ! The model geovals - graupel mixing ratios (kg/kg)
  real(kind_real), allocatable   :: dfed_dqg(:,:)       ! The deriv of FED vs layer qg conc (kg/sq.m/kg/kg) 
  logical :: ascend                                     ! Flag on direction of model levels
  integer :: iobs                                       ! Counter
!  integer :: ivar                                      ! Counter
  integer :: nlocs                                      ! number of observations
  integer :: ilev1                                      ! starting level for loop
  integer :: ilev2                                      ! ending level for loop
  integer :: ilev                                       ! level for loop
  integer :: inc                                        ! increment to level loop
  integer :: ibot                                       ! index of second lowest level
  integer :: isfc                                       ! index of lowest level
  integer :: n_horiz                                    ! number of square grid boxes for graupel
  integer, parameter      :: max_string = 800
  character(max_string)   :: message                    ! General message for output
  character(max_string)   :: err_msg                    ! Error Messages to be output to the user

  n_horiz = self % n_horiz

! get model state variables from geovals
  call ufo_geovals_get_var(geovals, var_delp, delp)     ! pressure (Pa)
  call ufo_geovals_get_var(geovals, var_qg, qg)         ! graupel mixing ratio (kg/kg)
! Make sure nothing already allocated
  call self%cleanup()

! Keep copy of dimensions for TL & AD
  self % nlevdP = delp % nval
  self % nlevq = qg % nval
  self % nlocs = obsspace_get_nlocs(obss)
!
! The model data must be on a staggered grid, with nlevp = nlevq+1
  IF (delp % nval /= qg % nval ) THEN
    write(err_msg,*) myname_ // ':' // &
    ' Data must be on a staggered grid nlevp, nlevq = ',delp%nval,qg%nval
    write(err_msg,*) myname_ // ':' // ' error: number of levels inconsistent!'
    call abor1_ftn(err_msg)
  END IF
!

! Determine if models levels are ascending or descending and set up indices for
! level loop
! Note assumes model levels stay same for all obs
!
  ascend = .false.
  if(delp%vals(1,1) > zero )ascend = .true.
  if (ascend)then       ! Model level starts above surface
    ilev1 = 1
    ilev2 = qg%nval
    inc = 1
    ibot = 2
    isfc = 1
  else                  ! Model level starts at top of atmosphere
    ilev1 = qg%nval
    ilev2 = 1
    inc = -1
    ibot = qg%nval
    isfc = qg%nval
  endif
!
  ALLOCATE(self % K(1:self%nlocs*self%n_horiz, 1:self%nlevq))
  ALLOCATE(dfed_dqg(1:self%nlevq,1:self%n_horiz))

! For each observation, calculate the K-matrix
  obs_loop: do iobs = 1, self % nlocs
      CALL lightning_GetK(  self % nlevdP,                                 &  ! Number of pressure levels
                          self % nlevq,                                    &  ! Number of hydrometeor levels
                          ilev1,                                           &  ! Starting pressure layer
                          ilev2,                                           &  ! Ending pressure layer
                          inc,                                             &  ! Layer increment
                          ibot,                                            &  ! Lowest layer inc surface
                          isfc,                                            &  ! Lowest level next to surface
                          delp % vals(:,(iobs-1)*n_horiz+1:iobs*n_horiz),  &  ! Pressure levels (Pa)
                          qg %vals (:,(iobs-1)*n_horiz+1:iobs*n_horiz),    &  ! horizontal integration of graupel mixing ratios (kg/kg)
                          dfed_dqg(:,1:n_horiz),                           &  ! deriv of FED vs layer qg concentration
                          !self % l_fed_nonlinear,                          &  ! option to use nonlinear observation operator
                          self % n_horiz)                                     ! No. of parterner grids in 15x15km^2 area

    ! Build K-matrix (Jacobian of FED with respect to layer qg)
     do ilev=1,self%nlevq
       self%K((iobs-1)*n_horiz+1 : iobs*n_horiz, ilev) = dfed_dqg(ilev,1:n_horiz)  
     enddo
  end do obs_loop
  deallocate(dfed_dqg)

end subroutine ufo_lightning_tlad_settraj_

! ------------------------------------------------------------------------------
! TODO: replace below function with your tl observation operator.
! Note: this can use information saved from trajectory in your ufo_lightning_tlad type
! Input geovals parameter represents dx for tangent linear model
subroutine ufo_lightning_simobs_tl_(self, geovals, obss, nvars, nlocs, hofx)
use iso_c_binding
implicit none
class(ufo_lightning_tlad), intent(in)  :: self
type(ufo_geovals),       intent(in)    :: geovals
integer,                 intent(in)    :: nvars, nlocs
real(c_double),          intent(inout) :: hofx(nvars,nlocs)
type(c_ptr), value,      intent(in)    :: obss

! Local parameters
  character(len=*), parameter  :: myname_="ufo_lightning_simobs_tl"

! Local variables
  integer                      :: iobs                 ! Loop variable, observation number
  integer                      :: ivar                 ! obervation variable
  integer                      :: ilev                 ! loop variable, level number
  integer                      :: nlev                 ! number of levels
  real(kind_real), allocatable :: x_d(:)               ! Increment to the state vector

  type(ufo_geoval), pointer    :: qg_d                 ! Pointer to the specific humidity perturbations
  type(ufo_geoval), pointer    :: delp_d               ! Pointer to the specific humidity perturbations
  integer                      :: n_horiz, i_horiz

! Get perturbed variables from geovals qghi (kg/kg)
!
  call ufo_geovals_get_var(geovals, var_qg, qg_d)      ! horizontal integration of graupel mixing ratio
  call ufo_geovals_get_var(geovals, var_delp, delp_d)  ! horizontal integration of graupel mixing ratio

  nlev = self % nlevq  ! number of layers
  hofx = zero
! Allocate state vector increment
!
 n_horiz = self%n_horiz
! Loop through the obs, calculating the increment to the observation hofx
  var_loop : do ivar =1,nvars
  obs_loop: do iobs = 1, nlocs
    lev_loop: do ilev = 1 , nlev
      nhoriz_loop: do i_horiz = 1, n_horiz
         hofx(ivar,iobs) =  hofx(ivar,iobs) + self % K((iobs-1)*n_horiz+i_horiz,ilev) * qg_d % vals(ilev, (iobs-1)*n_horiz+i_horiz)
      end do nhoriz_loop
    end do lev_loop
  end do obs_loop
  end do var_loop
!
!
  return
end subroutine ufo_lightning_simobs_tl_

! ------------------------------------------------------------------------------
! TODO: replace below function with your ad observation operator.
! Note: this can use information saved from trajectory in your ufo_lightning_tlad type
subroutine ufo_lightning_simobs_ad_(self, geovals, obss, nvars, nlocs, hofx)
use iso_c_binding
implicit none
class(ufo_lightning_tlad), intent(in)    :: self
type(ufo_geovals),         intent(inout) :: geovals
integer,                   intent(in)    :: nvars, nlocs
real(c_double),            intent(in)    :: hofx(nvars,nlocs)
type(c_ptr), value,        intent(in)    :: obss

! Local parameters
  character(len=*), parameter  :: myname_="ufo_lightning_simobs_ad"

! Local variables
  real(c_double)               :: missing  ! Missing data values

  type(ufo_geoval), pointer    :: qg_d      ! Pointer to the specific humidity perturbations
  type(ufo_geoval), pointer    :: delp_d    ! Pointer to the pressure interval perturbations

  integer                      :: iobs     ! Loop on obs, observation number
  integer                      :: ivar     ! Loop on observed variables
  integer                      :: ilev     ! Loop variable, level number
  integer                      :: nlev     ! level number
  real(kind_real), allocatable :: x_d(:,:) ! Perturbation to the state vector
  integer, parameter           :: max_string = 800
  character(max_string)        :: err_msg  ! Message to be output
  integer                      :: n_horiz, i_horiz


! setting up array for perturbed variables in geovals
  call ufo_geovals_get_var(geovals, var_qg,  qg_d)         ! layer graupel mixing ratio
  call ufo_geovals_get_var(geovals, var_delp,  delp_d)     ! layer pressure change
!
  n_horiz = self%n_horiz                                   ! Number of partner grids within a 15x15km^2 area
  nlev =  self % nlevq                                     ! Number of layers

! Allocate the output for the specific humidity perturbations
  if (.not. allocated(qg_d%vals)) then
      qg_d % nval = self % nlevq
      allocate(qg_d % vals(qg_d % nval, nlocs*n_horiz))
      qg_d % vals = zero
  endif

  if (.not. allocated(delp_d%vals)) then
      !delp_d % nlocs = self % nlocs
      delp_d % nval = self % nlevdP
      allocate(delp_d % vals(delp_d % nval, nlocs*n_horiz))
      delp_d % vals = zero
  endif


!
  missing = missing_value(missing)
! Allocate state vector x_d
  allocate(x_d(nlocs*n_horiz,1:nlev))
  
!
! Loop through the obs, calculating the increment to the model state layer q
!
  var_loop: do ivar = 1, nvars
  obs_loop: do iobs = 1, nlocs
    if (hofx(ivar,iobs) /= missing) then
      level_loop: do ilev = 1 , nlev
          nhoriz_loop: do i_horiz = 1, n_horiz
            x_d((iobs-1)*n_horiz+i_horiz,ilev) =  hofx(ivar,iobs) * self % K((iobs-1)*n_horiz+i_horiz,ilev)
            qg_d% vals(ilev,(iobs-1)*n_horiz+i_horiz) =  qg_d% vals(ilev,(iobs-1)*n_horiz+i_horiz) + x_d((iobs-1)*n_horiz+i_horiz,ilev)
          end do nhoriz_loop
      end do level_loop
    end if
  end do obs_loop
  end do var_loop
  deallocate(x_d)
!
!  write(err_msg,*) "TRACE: ufo_lightning_simobs_ad: complete"
!
  return


end subroutine ufo_lightning_simobs_ad_


! ------------------------------------------------------------------------------
subroutine ufo_lightning_tlad_cleanup_(self)
  implicit none
  class(ufo_lightning_tlad), intent(inout) :: self

  self%nlocs = 0
  if (allocated(self%K)) deallocate(self%K)

end subroutine ufo_lightning_tlad_cleanup_

subroutine  destructor(self)
  type(ufo_lightning_tlad), intent(inout)  :: self

  call self%cleanup()

end subroutine destructor

! ------------------------------------------------------------------------------



!------------------------------------------------------------
! Calculate the partial derivatives dFED/dqg
!------------------------------------------------------------
SUBROUTINE lightning_GetK(  nlevdP,         &
                          nlevq,            &
                          ilev1,            &
                          ilev2,            &
                          inc,              &
                          ibot,             &
                          isfc,             &
                          delp,             &
                          qg,               &
                          dfed_dqg,         &
                          !l_fed_nonlinear,  &
                          n_horiz)
!
IMPLICIT NONE

INTEGER, INTENT(IN)            :: nlevdP                    ! The number of model pressure levels
INTEGER, INTENT(IN)            :: nlevq                     ! The number of model hydrometeor levels
INTEGER, INTENT(IN)            :: ilev1                     ! Starting layer
INTEGER, INTENT(IN)            :: ilev2                     ! Ending layer
INTEGER, INTENT(IN)            :: inc                       ! loop increment
INTEGER, INTENT(IN)            :: ibot                      ! lowest layer including surface
INTEGER, INTENT(IN)            :: isfc                      ! lowest level next to surface
INTEGER, INTENT(IN)            :: n_horiz
REAL(kind_real), INTENT(IN)    :: qg(1:nlevq,n_horiz)       ! graupel mixing ratio
REAL(kind_real), INTENT(IN)    :: delp(1:nlevdP,n_horiz)    ! Model background pressure (Pa) at levels
REAL(kind_real), INTENT(INOUT) :: dfed_dqg(1:nlevq,n_horiz) ! Partial deriv of HofX for precip water wrt layer water vapour
REAL(kind_real)                :: dag_dqg
DOUBLE PRECISION               :: ag, answer
!
! Local constants
!
integer, parameter    :: max_string = 800
character(max_string) :: err_msg           ! Error message to be output
character(max_string) :: message           ! General message for output
character(len=*), parameter :: myname_ = "ufo_lightning_tlad"
REAL(kind_real)       :: PDiff             ! Pressure diff across layer
REAL(kind_real)       :: GK                ! 1/gravity
INTEGER               :: ilev              ! level counter

REAL(kind_real)                :: coeff3rdorder(4)
!LOGICAL, INTENT(IN)            :: l_fed_nonlinear
REAL(kind_real)                :: gridarea 
REAL(kind_real), parameter     :: glmcoeff = 2.088e-8
INTEGER                        :: i_horiz

gridarea = 15000*15000/n_horiz !the representaive area for each partner point
!-------------------------------------------------------------------------------
! The operator is calibrated using GSL's RFFS 1-3 hour forecasts, provided hourly from 00 to 23 Z,
! specifically for the period June 11 to June 15, 2023.
!-------------------------------------------------------------------------------
coeff3rdorder(1) = 0.030104143076945
coeff3rdorder(2) =  -0.441067352788872
coeff3rdorder(3) =  6.250421068636704
coeff3rdorder(4) =  0.0 !0.003091491729412
!-------------------------------------------------------------------------------
!1. Initialise variables
!-------------------------------------------------------------------------------
PDiff = zero
GK = one / grav
dfed_dqg = zero
dag_dqg = zero
ag = zero

!
!-------------------------------------------------------------------------------
! Calculate derivative dfed wrt dqg at all levels for K matrix
!-------------------------------------------------------------------------------
!
DO ilev = ilev1, ilev2, inc  !Used for vertical integration
   DO i_horiz = 1, n_horiz   !Used for horizontal integration
     PDiff = -1.0 * delp(ilev,i_horiz)     ! prs is pressure on level i
     if(PDiff>0)then
       write(0,*) 'PDiff is larger than zero, STOP HERE2'
     endif
     ag = ag - GK * qg(ilev,i_horiz) * PDiff * gridarea ! Accumulate layer graupel mass
   END DO
END DO
ag = ag * 1.0E-9  !The coefficients are fitted based on the relscaled data to avoid excessively small values 
!-------------------------------------------------------------------------------
! Then linearize around the nonlonlinear trajectory, here ag is used in calculation of nonlinear trajectory
!-------------------------------------------------------------------------------
DO ilev = ilev1, ilev2, inc
   DO i_horiz = 1, n_horiz
     PDiff = -1.0 * delp(ilev,i_horiz)     
     dag_dqg = -1.0*  GK * PDiff * gridarea
     answer = (3*coeff3rdorder(1)*ag**2 + 2*coeff3rdorder(2)*ag + coeff3rdorder(3))*dag_dqg
     dfed_dqg(ilev,i_horiz) = SNGL(answer)
  END DO
END DO

dfed_dqg = dfed_dqg * 1.0E-9 !don't forget this part since the observation operator is a composite function 
!
END SUBROUTINE lightning_GetK


! ------------------------------------------------------------------------------

end module ufo_lightning_tlad_mod
