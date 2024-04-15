! (C) Copyright 2020 Met Office
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

!> Fortran module containing the routines to perform a Marquardt-Levenberg 
!! minimization.

module ufo_rttovonedvarcheck_minimize_ml_mod

use kinds
use ufo_constants_mod, only: zero, ten
use fckit_log_module, only : fckit_log
use ufo_geovals_mod
use ufo_radiancerttov_mod
use ufo_rttovonedvarcheck_constants_mod
use ufo_rttovonedvarcheck_minimize_jacobian_mod
use ufo_rttovonedvarcheck_minimize_utils_mod
use ufo_rttovonedvarcheck_profindex_mod
use ufo_rttovonedvarcheck_rsubmatrix_mod
use ufo_rttovonedvarcheck_ob_mod
use ufo_rttovonedvarcheck_setup_mod
use ufo_rttovonedvarcheck_utils_mod
use ufo_utils_mod, only: Ops_Cholesky
use ufo_vars_mod

implicit none
private

public ufo_rttovonedvarcheck_minimize_ml

contains

!-------------------------------------------------------------------------------
!> Perform a minimization using the Marquardt-Levenberg method.
!!
!! \details Heritage: Ops_SatRad_MinimizeML_RTTOV12.f90
!!
!! Find the most probable atmospheric state vector by minimizing a cost function
!! through a series of iterations. if a solution exists, the iterations will
!! converge when the iterative increments are acceptably small. A limit on the
!! total number of iterations allowed is imposed.
!!
!! Using the formulation given by Rodgers (1976) :
!!
!!   Delta_x = xn + (xb-xn).I' + Wn.(ym-y(xn) - H.(xb-xn)) <br>
!!   where: <br>
!!          x is an atmospheric state vector, subscripted b=background,n=nth
!!           iteration <br>
!!          I' is a diagonal matrix with I'(J,J) = B_damped(J,J)/B_undamped(J,J)
!!           although, however, damping will no longer be used <br>
!!          Wn = B.Hn'.(Hn.B.Hn'+R)^-1 <br>
!!          B is the background error covariance matrix <br>
!!          R is the combined forward model and ob error covariance matrix <br>
!!
!!   Delta_x is checked for convergence after each iteration
!!
!! The loop is exited with convergence if either of the following conditions are
!! true, depending on whether UseJforConvergence is true or false
!!   -# if UseJforConvergence is true then the change in total cost function from
!!      one iteration to the next is sufficiently small to imply convergence
!!   -# if UseJforConvergence is false then the increments to the atmospheric
!!      state vector are sufficiently small to imply convergence at an acceptable
!!      solution
!!
!! Either of the following two conditions will cause the 1dvar to stop and exit
!! with an error.
!!   -# The increments are sufficiently large to suppose a solution will not
!!      be found.
!!   -# The maximum number of allowed iterations has been reached. In most
!!      cases, one of the above criteria will have occurred.
!!
!! References:
!!
!!   Rodgers, Retrieval of atmospheric temperature and composition from
!!   remote measurements of thermal radiation, Rev. Geophys.Sp.Phys. 14, 1976.
!!
!!   Eyre, Inversion of cloudy satellite sounding radiances by nonlinear
!!   optimal estimation. I: Theory and simulation for TOVS,QJ,July 89.
!!
!! \author Met Office
!!
!! \date 09/06/2020: Created
!!
subroutine ufo_rttovonedvarcheck_minimize_ml(config,    &
                                         ob,            &
                                         r_matrix,      &
                                         b_matrix,      &
                                         b_inv,         &
                                         b_sigma,       &
                                         firstguess_geovals, &
                                         hofxdiags,     &
                                         rttov_simobs,  &
                                         profile_index, &
                                         onedvar_success)

implicit none

type(ufo_rttovonedvarcheck), intent(inout) :: config !< structure containing settings
type(ufo_rttovonedvarcheck_ob), intent(inout) :: ob  !< satellite metadata
type(ufo_rttovonedvarcheck_rsubmatrix), intent(inout) :: r_matrix !< observation error covariance
real(kind_real), intent(in)       :: b_matrix(:,:)   !< state error covariance
real(kind_real), intent(in)       :: b_inv(:,:)      !< inverse state error covariance
real(kind_real), intent(in)       :: b_sigma(:)      !< standard deviations of the state error covariance diagonal
type(ufo_geovals), intent(inout)  :: firstguess_geovals !< background model data at obs location
type(ufo_geovals), intent(inout)  :: hofxdiags       !< model data containing the jacobian
type(ufo_radiancerttov), intent(inout) :: rttov_simobs !< rttov simulated obs object
type(ufo_rttovonedvarcheck_profindex), intent(in) :: profile_index !< index array for x vector
logical, intent(out)              :: onedvar_success !< convergence flag

! Local declarations:
character(len=*), parameter     :: RoutineName = "ufo_rttovonedvarcheck_minimize_ml"
integer                         :: inversionstatus ! status of minimisation
logical                         :: outOfRange      ! out of range flag for check iteration
logical                         :: Converged       ! converged flag
logical                         :: Error           ! error flag
integer                         :: iter            ! iteration counter
integer                         :: nchans          ! number of satellite channels used
integer                         :: ii, jj          ! counters
integer                         :: usedchan, allchans ! counters
integer                         :: nprofelements   ! number of profile elements
real(kind_real)                 :: Gamma           ! steepness of descent parameter
real(kind_real)                 :: Jcost           ! current cost value
real(kind_real)                 :: JcostOld        ! previous iteration cost value
real(kind_real)                 :: JcostOrig       ! initial cost value
real(kind_real)                 :: DeltaJ          ! change in cost during iteration
real(kind_real)                 :: DeltaJo         ! change in cost from original

real(kind_real), allocatable    :: OldProfile(:)
real(kind_real), allocatable    :: GuessProfile(:)
real(kind_real), allocatable    :: GuessProfileBefore(:)
real(kind_real), allocatable    :: BackProfile(:)
real(kind_real), allocatable    :: H_matrix(:,:)
real(kind_real), allocatable    :: Diffprofile(:)
real(kind_real), allocatable    :: AbsDiffProfile(:)
real(kind_real), allocatable    :: Ydiff(:)
real(kind_real), allocatable    :: Y(:)
real(kind_real), allocatable    :: Y0(:)
real(kind_real), allocatable    :: output_BT_usedchans(:)
real(kind_real)                 :: Jout(3)
real(kind_real)                 :: tskin
type(ufo_geovals)               :: geovals
type(ufo_geoval), pointer       :: geoval

!------------
! Setup
!------------

Converged = .false.
onedvar_success = .false.
Error = .false.
nchans = size(ob % channels_used)
Gamma = 1.0e-4_kind_real
Jcost = 1.0e4_kind_real
JcostOld = 1.0e4_kind_real
nprofelements = profile_index % nprofelements
allocate(OldProfile(nprofelements))
allocate(GuessProfile(nprofelements))
allocate(GuessProfileBefore(nprofelements))
allocate(BackProfile(nprofelements))
allocate(H_matrix(nchans,nprofelements))
allocate(Diffprofile(nprofelements))
allocate(AbsDiffprofile(nprofelements))
allocate(Ydiff(nchans))
allocate(Y(nchans))
allocate(Y0(nchans))
call ufo_geovals_copy(firstguess_geovals, geovals)

call fckit_log % debug("Using ML solver")

! Map GeovaLs to 1D-var profile using B matrix profile structure
call ufo_rttovonedvarcheck_GeoVaLs2ProfVec(geovals, config, profile_index, ob, GuessProfile(:))

Iterations: do iter = 1, config % max1DVarIterations

  !-----------------------------------------------------------
  ! 0.5. Reset obs errors of specified channels in Rmatrix
  !      if non-converging. Prevent offending channels passed
  !      to Var if SatInfo%DontAssimSlowConvergedObsVar set.
  !      This behaviour is disabled when retrieving MW emissivity
  !-----------------------------------------------------------
  if (iter > config % ConvergeCheckChansAfterIteration .and. &
      allocated(config % ConvergeCheckChans) .and. &
      profile_index % mwemiss(1) == 0) then
    if (size(config % ConvergeCheckChans) > 0) then
      call r_matrix % reset_errors(config % ConvergeCheckChans, 100000.0_kind_real)
      ob % QC_SlowConvChans = .true.
    end if
  endif

  !-------------------------
  ! 1. Generate new profile
  !-------------------------

  ! Save current profile
  OldProfile(:) = GuessProfile(:)

  ! Get jacobian and new hofx
  call ufo_rttovonedvarcheck_get_jacobian(config, geovals, ob, ob % channels_used, &
                                       profile_index, GuessProfile(:), &
                                       hofxdiags, rttov_simobs, Y(:), H_matrix)

  if (ob % rterror) then
    call fckit_log % warning("Radiative transfer error exit and reject profile")
    exit Iterations
  end if

  if (iter == 1) then    
    Diffprofile(:) = zero
    BackProfile(:) = GuessProfile(:)
    Y0(:) = Y(:)
    call ufo_rttovonedvarcheck_subset_to_all_by_channels(ob % channels_used, Y0, &
                                     ob % channels_all, ob % background_BT)
  end if

  ! Calculate Obs diff and initial cost
  Ydiff(:) = ob % yobs(:) - Y(:)
  if (iter == 1) then
    Diffprofile(:) = GuessProfile(:) - BackProfile(:)
    call ufo_rttovonedvarcheck_CostFunction(Diffprofile, b_inv, Ydiff, r_matrix, Jout)
    Jcost = Jout(1)
    JCostOrig = Jcost
  end if

  !-----------------------------------------------------
  ! 1a. if UseJForConvergence is true
  !     then compute costfunction for latest profile
  !     and determine convergence using change in cost fn
  !-----------------------------------------------------

  if (config % UseJForConvergence) then

    ! exit on error
    !if (inversionStatus /= 0) exit Iterations

    ! check for convergence
    if (iter > 1) then

      if (config % JConvergenceOption == 1) then

        ! percentage change tested between iterations
        DeltaJ = abs ((Jcost - JcostOld) / max (Jcost, tiny (zero)))

        ! default test for checking that overall cost is getting smaller
        DeltaJo = -1.0_kind_real

      else

        ! absolute change tested between iterations
        DeltaJ = abs (Jcost - JcostOld)

        ! change between current cost and initial
        DeltaJo = Jcost - JCostorig

      end if

      if (DeltaJ < config % cost_convergencefactor .and. &
          DeltaJo < zero)  then ! overall is cost getting smaller?
        converged = .true.
        exit iterations
      end if

    end if

  end if ! end of specific code for cost test convergence

  ! store cost function from previous cycle
  JcostOld = Jcost

  ! Useful diagnostics to check the minimization step by step
  if (config % FullDiagnostics) then
    write(*,*) "Iter info outer loop"
    call ufo_rttovonedvarcheck_PrintIterInfo(ob % yobs(:), Y(:), ob % channels_used, &
                                             guessprofile, backprofile, &
                                             diffprofile, b_inv, H_matrix, &
                                             r_matrix % diagonal(:))
  end if

  ! Linearly iterate (Guess) profile vector
  call ufo_rttovonedvarcheck_ML ( config,              &
                                  Ydiff,               &
                                  nChans,              &
                                  ob,                  &
                                  H_Matrix,            &
                                  transpose(H_Matrix), &
                                  nprofelements,       &
                                  profile_index,       &
                                  DiffProfile,         &
                                  GuessProfile,        &
                                  BackProfile,         &
                                  b_inv,               &
                                  r_matrix,            &
                                  geovals,             &
                                  hofxdiags,           &
                                  rttov_simobs,        &
                                  gamma,               &
                                  Jcost,               &
                                  inversionStatus)

  if (inversionStatus /= 0) then
    inversionStatus = 1
    call fckit_log % warning("inversion failed exiting iterations")
    exit Iterations
  end if

  !---------------------------------------------------------
  ! 2. Check new profile and transfer to forward model format
  !---------------------------------------------------------

  ! Check profile and constrain humidity variables
  call ufo_rttovonedvarcheck_CheckIteration (config, & ! in
                                  geovals,           & ! in
                                  profile_index,     & ! in
                                  ob,                & ! in
                                  GuessProfile(:),   & ! inout
                                  outOfRange)          ! out

  ! Update geovals to be the same as guess profile
  call ufo_rttovonedvarcheck_ProfVec2GeoVaLs(geovals, config, profile_index, &
                                             ob, GuessProfile)
  
  ! if qtotal in retrieval vector check cloud
  ! variables for current iteration

  if ((.not. outofRange) .and. profile_index % qt(1) > 0) then

    if (iter >= config % IterNumForLWPCheck) then

        call ufo_rttovonedvarcheck_CheckCloudyIteration( config,  & ! in
                                              geovals,            & ! in
                                              profile_index,      & ! in
                                              config % nlevels,   & ! in
                                              OutOfRange )          ! out
    end if

  end if

  !-------------------------------------------------
  ! 3. Check for convergence using change in profile
  !    This is performed if UseJforConvergence is false
  !-------------------------------------------------

  absDiffprofile(:) = zero

  if ((.not. outOfRange) .and. (.not. config % UseJForConvergence)) then
    absDiffProfile(:) = abs(GuessProfile(:) - OldProfile(:))
    if (ALL (absDiffProfile(:) <= B_sigma(:) * config % ConvergenceFactor)) then
      call fckit_log % debug("Profile used for convergence")
      Converged = .true.
    end if
  end if

  !---------------------
  ! 4. output diagnostics
  !---------------------

  if (config % FullDiagnostics) then
    write (*, '(A,I0)') 'Iteration', iter
    write (*, '(A)') '------------'
    write (*, '(A,L1)') 'Status: converged = ', Converged
    if (outOfRange) write (*, '(A)') 'exiting with bad increments'
    write (*, '(A)') 'New profile:'
    call ufo_geovals_print(geovals, 1)
    call ob % info()
    write (*, '(A)')
  end if

  ! exit conditions

  if (Converged .OR. outOfRange) exit iterations

end do Iterations

! Pass convergence flag out
onedvar_success = converged

! Put final values in the ob structure
ob % final_cost = Jcost
ob % niter = iter

! Pass output profile, final BTs and final cost out
if (converged) then
  ob % output_profile(:) = GuessProfile(:)

  ! If lwp output required then recalculate
  if (config % Store1DVarLWP .or. config % Store1DVarIWP) then
    call ufo_rttovonedvarcheck_CheckCloudyIteration( config,    & ! in
                                            geovals,            & ! in
                                            profile_index,      & ! in
                                            config % nlevels,   & ! in
                                            OutOfRange,         & ! out
                                            OutLWP = ob % LWP,  & ! out
                                            OutIWP = ob % IWP)    ! out
  end if

  ! store final clw if required
  if (allocated(ob % clw)) then
    call ufo_geovals_get_var(geovals, var_clw, geoval)
    ob % clw = geoval%vals(:, 1)
  end if
  
  ! Recalculate final BTs for all channels
  call ufo_rttovonedvarcheck_get_bts(config, geovals, ob, ob % channels_all, &
                                     rttov_simobs, ob % output_BT)

  ! Recalculate BTs for all channels using the background profile and
  ! update the Tskin over land, surface emissivity, ctp, eca if retrieved.
  ! surface emissivity, ctp and eca are in the ob structure
  if (config % RecalculateBT) then
    ! Update Tskin
    call ufo_geovals_get_var(geovals, var_sfc_tskin, geoval)
    tskin = geoval%vals(1, 1)
    call ufo_geovals_get_var(firstguess_geovals, var_sfc_tskin, geoval)
    geoval%vals(1, 1) = tskin
    ! Calculate BT
    call ufo_rttovonedvarcheck_get_bts(config, firstguess_geovals, ob, ob % channels_all, &
                                       rttov_simobs, ob % recalc_BT)
  end if

  ! Fill the final BT diff
  allocate(output_BT_usedchans(size(ob % channels_used)))
  call ufo_rttovonedvarcheck_all_to_subset_by_channels(ob % channels_all, &
                             ob % output_BT, ob % channels_used, output_BT_usedchans)
  ob % final_bt_diff(:) = ob % yobs(:) - output_BT_usedchans(:)
  deallocate(output_BT_usedchans)

  ! If ctp retrieved then check 1 % of the integrated Jacobian for a channel is not below the
  ! cloud top.  If it is then reject that channel.
  if(config % cloud_retrieval .and. ob % cloudtopp > zero .and. &
    ob % cloudfrac > 0.05_kind_real) then
    call ufo_rttovonedvarcheck_ctp_error(geovals, r_matrix, profile_index, H_matrix, ob)
    call ufo_rttovonedvarcheck_cloudy_channel_rejection( &
         config, profile_index, geovals, H_matrix, ob)
  end if
end if

!----------------------
! 4. output diagnostics
!----------------------

if (config % UseJForConvergence .and. config % FullDiagnostics) then
  write(*,'(A45,3F10.3,I5,L5)') "ML J initial, final, lowest, iter, converged = ", &
                                 JCostorig, Jcost,  Jcost, iter, onedvar_success
end if

! ----------
! Tidy up
! ----------
if (allocated(OldProfile))         deallocate(OldProfile)
if (allocated(GuessProfile))       deallocate(GuessProfile)
if (allocated(GuessProfileBefore)) deallocate(GuessProfileBefore)
if (allocated(BackProfile))        deallocate(BackProfile)
if (allocated(H_matrix))           deallocate(H_matrix)
if (allocated(Diffprofile))        deallocate(Diffprofile)
if (allocated(AbsDiffprofile))     deallocate(AbsDiffprofile)
if (allocated(Ydiff))              deallocate(Ydiff)
if (allocated(Y))                  deallocate(Y)
if (allocated(Y0))                 deallocate(Y0)

end subroutine ufo_rttovonedvarcheck_minimize_ml

!---------------------------------------------------------------------
!> Updates the profile vector during the Marquardt-Levenberg 
!! minimization.
!!
!! \details Heritage: Ops_SatRad_MarquardtLevenberg_RTTOV12.f90
!!
!! Updates the profile vector DeltaProfile according to Rodgers (1976), Eqn. 100,
!! extended to allow for additional cost function terms and Marquardt-Levenberg
!! descent.
!!
!!   x_(n+1) = xb + U^-1.V <br>
!!   where <br>
!!         U=(B^-1 + H^T R^-1 H + J2) <br>
!!         V=H^T R^-1 [(ym-y(x_n))+H(x_n-xb)] - J1 <br>
!!   and   <br>
!!         J_extra=J0+J1.(x-xb)+(x-xb)^T.J2.(x-xb) is the additional cost
!!         function <br>
!!         x is an atmospheric state vector, subscripted b=background,n=nth
!!         iteration <br>
!!         ym is the measurement vector (i.e. observed brightness temperatures) <br>
!!         y(xn) is the observation vector calculated for xn <br>
!!         ym and y(xn) are not used individually at all, hence these are input <br>
!!         as a difference vector DeltaBT. <br>
!!         B is the background error covariance matrix <br>
!!         R is the combined forward model and ob error covariance matrix <br>
!!         H is the forward model gradient (w.r.t. xn) matrix <br>
!!         H' is the transpose of H <br>
!!
!!   When J_extra is zero this is simply Rogers (1976), Eqn. 100.
!!
!!   U^-1.V is solved using Cholesky decomposition.
!!
!! This routine should be used when:
!!     -# The length of the observation vector is greater than the length
!!        of the state vector,
!!     -# Additional cost function terms are provided and
!!     -# where Newtonian minimisation is desired.
!!
!! References:
!!
!!   Rodgers, Retrieval of atmospheric temperature and composition from
!!   remote measurements of thermal radiation, Rev. Geophys.Sp.Phys. 14, 1976.
!!
!!   Rodgers, Inverse Methods for Atmospheres: Theory and Practice.  World
!!            Scientific Publishing, 2000.
!!
!! \author Met Office
!!
!! \date 09/06/2020: Created
!!
subroutine ufo_rttovonedvarcheck_ML ( config,        &
                                      DeltaBT,       &
                                      nChans,        &
                                      ob,            &
                                      H_Matrix,      &
                                      H_Matrix_T,    &
                                      nprofelements, &
                                      profile_index, &
                                      DeltaProfile,  &
                                      GuessProfile,  &
                                      BackProfile,   &
                                      b_inv,         &
                                      r_matrix,      &
                                      geovals,       &
                                      hofxdiags,     &
                                      rttov_simobs,  &
                                      gamma,         &
                                      Jold,          &
                                      Status)

implicit none

! subroutine arguments:
type(ufo_rttovonedvarcheck), intent(in) :: config          !< structure containing settings
real(kind_real), intent(in)             :: DeltaBT(:)      !< y-y(x)
integer, intent(in)                     :: nChans          !< number of channels
type(ufo_rttovonedvarcheck_ob), intent(inout) :: ob        !< satellite metadata
real(kind_real), intent(in)             :: H_Matrix(:,:)   !< Jacobian
real(kind_real), intent(in)             :: H_Matrix_T(:,:) !< transpose  of the Jacobian
integer, intent(in)                     :: nprofelements   !< number of profile elements
type(ufo_rttovonedvarcheck_profindex), intent(in) :: profile_index !< index array for x vector
real(kind_real), intent(inout)          :: DeltaProfile(:) !< x_(n+1) - xb
real(kind_real), intent(inout)          :: GuessProfile(:) !< x_(n+1)
real(kind_real), intent(in)             :: BackProfile(:)  !< xb
real(kind_real), intent(in)             :: b_inv(:,:)      !< inverse of the state error covariance
type(ufo_rttovonedvarcheck_rsubmatrix), intent(in) :: r_matrix !< observation error covariance
type(ufo_geovals), intent(inout)        :: geovals         !< model data at obs location
type(ufo_geovals), intent(inout)        :: hofxdiags       !< model data containing the jacobian
type(ufo_radiancerttov), intent(inout)  :: rttov_simobs
real(kind_real), intent(inout)          :: gamma           !< steepness of descent parameter
real(kind_real), intent(inout)          :: Jold            !< previous steps cost which gets update by good step
integer, intent(out)                    :: Status          !< code to capture failed Cholesky decomposition

! Local declarations:
character(len=*), parameter         :: RoutineName = 'ufo_rttovonedvarcheck_ML'
integer                             :: i
integer                             :: Iter_ML             ! M-L iteration count
real(kind_real)                     :: JOut(3)             ! J, Jb, Jo
real(kind_real)                     :: JCost               ! Cost function,
logical                             :: OutOfRange
real(kind_real)                     :: HTR(nprofelements, nChans)              ! Scratch vector
real(kind_real)                     :: U(nprofelements, nprofelements)         ! U = H.B.H^T + R
real(kind_real)                     :: V(nprofelements)                        ! V = (y-y(x_n))-H^T(xb-x_n)
real(kind_real)                     :: USave(nprofelements, nprofelements)     ! U = H.B.H^T + R
real(kind_real)                     :: VSave(nprofelements)                    ! V = (y-y(x_n))-H^T(xb-x_n)
real(kind_real)                     :: New_DeltaProfile(nprofelements)
real(kind_real)                     :: BriTemp(nchans)                 ! Forward modelled brightness temperatures
real(kind_real)                     :: Ydiff(nchans)
real(kind_real)                     :: Xdiff(nprofelements)
real(kind_real)                     :: Emiss(nchans)
real(kind_real)                     :: H_matrix_tmp(nchans,nprofelements)
logical                             :: CalcEmiss(nchans)
integer, parameter                  :: StatusOK = 0
integer                             :: RTStatus

Status = StatusOK
RTStatus = 0

!---------------------------------------------------------------------------
! 1. Calculate HTR-1 for the three allowed forms of R matrix. if required, the
!    matrix is tested to determine whether it is stored as an inverse and
!    inverted if not.
!---------------------------------------------------------------------------

!HTR(-1) = matmul(H_matrix_T, R_inverse)
call r_matrix % multiply_inverse_matrix(H_matrix_T, HTR)

!---------------------------------------------------------------------------
! 2. Calculate U and V
!---------------------------------------------------------------------------

U = matmul (HTR, H_matrix)
V = matmul (HTR, DeltaBT)
V = V + matmul (U, DeltaProfile)

!---------------------------------------------------------------------------
! 3. Add on inverse of B-matrix to U.
!---------------------------------------------------------------------------

U = U + b_inv

!---------------------------------------------------------------------------
! 5. Start Marquardt-Levenberg loop.
!    Search for a value of gamma that reduces the cost function.
!---------------------------------------------------------------------------

JCost = 1.0e10_kind_real
Gamma = Gamma / ten
USave = U
VSave = V
Iter_ML = 0

DescentLoop : do while (JCost > JOld .and.                       &
                        Iter_ML < config % MaxMLIterations .and. &
                        Gamma <= 1.0e25_kind_real)

  Iter_ML = Iter_ML + 1

  !------------------------------------------------------------------------
  ! 5.1 Add on Marquardt-Levenberg terms to U and V.
  !------------------------------------------------------------------------

  Gamma = Gamma * ten

  do I = 1, nprofelements
    U(I,I) = USave(I,I) + Gamma
  end do
  V = VSave + Gamma * DeltaProfile

  !------------------------------------------------------------------------
  ! 5.2 Calculate new profile increment.
  !------------------------------------------------------------------------

  call Ops_Cholesky (U,                &
                     V,                &
                     nprofelements,    &
                     New_DeltaProfile, &
                     Status)
  if (Status /= 0) then
     call fckit_log % warning("Error in Cholesky decomposition")
  end if

  !------------------------------------------------------------------------
  ! 5.3. Calculate the radiances for the new profile.
  !------------------------------------------------------------------------

  GuessProfile(:) = BackProfile(:) + New_DeltaProfile(:)

  call ufo_rttovonedvarcheck_CheckIteration (config, & ! in
                                  geovals,           & ! in
                                  profile_index,     & ! in
                                  ob,                & ! in
                                  GuessProfile(:),   & ! inout
                                  outOfRange)          ! out

  if (.not. OutOfRange) then
    call ufo_rttovonedvarcheck_ProfVec2GeoVaLs(geovals, config, profile_index, &
                                               ob, GuessProfile)

    ! Get new hofx
    call ufo_rttovonedvarcheck_get_bts(config, geovals, ob, ob % channels_used, &
                                       rttov_simobs, BriTemp(:))
    if (ob % rterror) RTStatus = 1

    !------------------------------------------------------------------------
    ! 5.6. Calculate the new cost function.
    !------------------------------------------------------------------------

    if (RTStatus > 0) then

      ! In case of RT error, set cost to large value and
      ! continue with next iteration

      Jcost = 1.0e10_kind_real
    else
      Xdiff(:) = GuessProfile(:) - BackProfile(:)
      Ydiff(:) = ob % yobs(:) - BriTemp(:)
      call ufo_rttovonedvarcheck_CostFunction(Xdiff, b_inv, Ydiff, r_matrix, Jout)
      Jcost = Jout(1)
      if (Status /= StatusOK) return
    end if

  else

    ! if profile out of range, set cost to large value and
    ! continue with next iteration

    Jcost = 1.0e10_kind_real

  end if

  ! Useful diagnostics to check the minimization step by step
  if (config % FullDiagnostics) then
    write(*,*) "Iter info inner loop"
    call ufo_rttovonedvarcheck_PrintIterInfo(ob % yobs(:), BriTemp, ob % channels_used, &
                                             guessprofile, backprofile, &
                                             Xdiff, b_inv, H_matrix, r_matrix % diagonal(:))
  end if

end do DescentLoop

DeltaProfile = New_DeltaProfile
Gamma = Gamma / ten
JOld = JCost

end subroutine ufo_rttovonedvarcheck_ML

end module ufo_rttovonedvarcheck_minimize_ml_mod
