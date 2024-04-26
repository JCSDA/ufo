! (C) Copyright 2020 Met Office
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

!> Fortran module to perform newton steepest descent minimization

module ufo_rttovonedvarcheck_minimize_newton_mod

use kinds
use ufo_constants_mod, only: zero
use fckit_log_module, only : fckit_log
use ufo_geovals_mod
use ufo_radiancerttov_mod
use ufo_rttovonedvarcheck_constants_mod
use ufo_rttovonedvarcheck_minimize_jacobian_mod
use ufo_rttovonedvarcheck_minimize_utils_mod
use ufo_rttovonedvarcheck_ob_mod
use ufo_rttovonedvarcheck_profindex_mod
use ufo_rttovonedvarcheck_rsubmatrix_mod
use ufo_rttovonedvarcheck_setup_mod
use ufo_rttovonedvarcheck_utils_mod
use ufo_utils_mod, only: Ops_Cholesky
use ufo_vars_mod

implicit none
private

! public subroutines
public ufo_rttovonedvarcheck_minimize_newton

contains

!------------------------------------------------------------------------------
!> Get the jacobian used in the 1D-Var.
!!
!! \details Heritage: Ops_SatRad_MinimizeNewton_RTTOV12.f90
!!
!! Find the most probable atmospheric state vector by minimizing a cost function
!! through a series of iterations. if a solution exists, the iterations will
!! converge when the iterative increments are acceptably small. A limit on the
!! total number of iterations allowed, is imposed.
!!
!! Using the formulation given by Rodgers (1976) :
!!
!!   Delta_x = xn + (xb-xn).I' + Wn.(ym-y(xn) - H.(xb-xn)) <br>
!!   where: <br> x is an atmospheric state vector, subscripted b=background,n=nth
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
!! true, depending on whether UseJforConvergence is true or false <br>
!!   -# if UseJforConvergence is true then the change in total cost function from
!!      one iteration to the next is sufficiently small to imply convergence
!!   -# if UseJforConvergence is false then the increments to the atmospheric
!!      state vector are sufficiently small to imply convergence at an acceptable
!!      solution <br>
!!   Either of the following two conditions will cause the 1dvar to stop and exit
!!   with an error. <br>
!!   -# The increments are sufficiently large to suppose a solution will not
!!      be found.
!!   -# The maximum number of allowed iterations has been reached. in most
!!      cases, one of the above criteria will have occurred.
!!
!! References:
!!
!!   Rodgers, Retrieval of atmospheric temperature and composition from
!!   remote measurements of thermal radiation, Rev. Geophys.Sp.Phys. 14, 1976.
!!
!!   Eyre, inversion of cloudy satellite sounding radiances by nonlinear
!!   optimal estimation. I: Theory and simulation for TOVS,QJ,July 89.
!!
!! \author Met Office
!!
!! \date 09/06/2020: Created
!!
subroutine ufo_rttovonedvarcheck_minimize_newton(config,  &
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

type(ufo_rttovonedvarcheck), intent(inout) :: config !< Main 1D-Var object
type(ufo_rttovonedvarcheck_ob), intent(inout) :: ob  !< satellite metadata
type(ufo_rttovonedvarcheck_rsubmatrix), intent(inout) :: r_matrix !< observation error covariance
real(kind_real), intent(in)       :: b_matrix(:,:)   !< state error covariance
real(kind_real), intent(in)       :: b_inv(:,:)      !< inverse of the state error covariance
real(kind_real), intent(in)       :: b_sigma(:)      !< standard deviations of the state error covariance diagonal
type(ufo_geovals), intent(inout)  :: firstguess_geovals !< background model data at obs location
type(ufo_geovals), intent(inout)  :: hofxdiags       !< model data containing the jacobian
type(ufo_radiancerttov), intent(inout) :: rttov_simobs
type(ufo_rttovonedvarcheck_profindex), intent(in) :: profile_index !< index array for x vector
logical, intent(out)              :: onedvar_success !< convergence flag

! Local declarations:
character(len=*), parameter     :: RoutineName = "ufo_rttovonedvarcheck_minimize_newton"
integer                         :: inversionstatus
logical                         :: outOfRange
logical                         :: Converged
logical                         :: Error
integer                         :: iter
integer                         :: nchans
integer                         :: nprofelements
integer                         :: usedchan, allchans ! counters
real(kind_real)                 :: Jcost     ! current value
real(kind_real)                 :: JcostOld  ! previous iteration value
real(kind_real)                 :: JcostOrig ! initial value
real(kind_real)                 :: DeltaJ
real(kind_real)                 :: DeltaJo

real(kind_real), allocatable    :: OldProfile(:)
real(kind_real), allocatable    :: GuessProfile(:)
real(kind_real), allocatable    :: BackProfile(:)
real(kind_real), allocatable    :: H_matrix(:,:)
real(kind_real), allocatable    :: Diffprofile(:)
real(kind_real), allocatable    :: AbsDiffProfile(:)
real(kind_real), allocatable    :: Xdiff(:)
real(kind_real), allocatable    :: Ydiff(:)
real(kind_real), allocatable    :: Y(:)
real(kind_real), allocatable    :: Y0(:)
real(kind_real), allocatable    :: output_BT_usedchans(:)
real(kind_real)                 :: Jout(3)
real(kind_real)                 :: tskin
type(ufo_geovals)               :: geovals
type(ufo_geoval), pointer       :: geoval

integer                         :: ii, jj

! ---------
! Setup
! ---------
Converged = .false.
onedvar_success = .false.
Error = .false.
nchans = size(ob % channels_used)
inversionstatus = 0
nprofelements = profile_index % nprofelements
allocate(OldProfile(nprofelements))
allocate(GuessProfile(nprofelements))
allocate(BackProfile(nprofelements))
allocate(H_matrix(nchans,nprofelements))
allocate(Diffprofile(nprofelements))
allocate(AbsDiffprofile(nprofelements))
allocate(Xdiff(nprofelements))
allocate(Ydiff(nchans))
allocate(Y(nchans))
allocate(Y0(nchans))
call ufo_geovals_copy(firstguess_geovals, geovals)

if (config % FullDiagnostics) call ufo_geovals_print(geovals,1)
call fckit_log % debug("Using Newton solver")

JCost = 1.0e4_kind_real

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
  ! Save cost from previous iteration
  if (config % UseJForConvergence) then
    JcostOld = Jcost
  end if

  ! initialise profile increments on first iteration
  if (iter == 1) then

    Diffprofile(:) = zero
    JcostOld = 1.0e4_kind_real

    ! Map GeovaLs to 1D-var profile using B matrix profile structure
    call ufo_rttovonedvarcheck_GeoVaLs2ProfVec(geovals, config, profile_index, &
                                               ob, GuessProfile(:))

    if (config % FullDiagnostics) &
      write(*,*) "Humidity GuessProfile 1st iteration = ",GuessProfile(profile_index % qt(1):profile_index % qt(2))

  end if

  ! Save current profile
  OldProfile(:) = GuessProfile(:)

  ! Get jacobian and hofx
  call ufo_rttovonedvarcheck_get_jacobian(config, geovals, ob, ob % channels_used, &
                                          profile_index, GuessProfile(:), &
                                          hofxdiags, rttov_simobs, Y(:), H_matrix)

  if (ob % rterror) then
    call fckit_log % warning("Radiative transfer error exit and reject profile")
    exit Iterations
  end if

  if (iter == 1) then
    BackProfile(:) = GuessProfile(:)
    Y0(:) = Y(:)
    call ufo_rttovonedvarcheck_subset_to_all_by_channels(ob % channels_used, Y0, &
                                     ob % channels_all, ob % background_BT)
  end if

  !-----------------------------------------------------
  ! 1a. if UseJForConvergence is true
  !     then compute costfunction for latest profile
  !     and determine convergence using change in cost fn
  !-----------------------------------------------------

  ! Profile differences
  Xdiff(:) = GuessProfile(:) - BackProfile(:)
  Ydiff(:) = ob % yobs(:) - Y(:)

  if (config % FullDiagnostics) then
    write(*,*) "Ob BT = "
    write(*,'(10F10.3)') ob % yobs(:)
    write(*,*) "HofX BT = "
    write(*,'(10F10.3)') Y(:)
    call ufo_rttovonedvarcheck_PrintIterInfo(ob % yobs(:), Y(:), ob % channels_used, &
                                             guessprofile, backprofile, &
                                             Xdiff, b_inv, H_matrix, r_matrix % diagonal(:))
  end if

  if (config % UseJForConvergence) then

    call ufo_rttovonedvarcheck_CostFunction(Xdiff, b_inv, Ydiff, r_matrix, Jout)
    Jcost = Jout(1)

    ! exit on error
    if (inversionStatus /= 0) exit Iterations

    ! store initial cost value
    if (iter == 1) JCostOrig = jcost

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

      if (config % FullDiagnostics) THEN
        write (*, '(A,F12.5)') 'Cost Function = ', Jcost
        write (*, '(A,F12.5)') 'Cost Function old = ', JcostOld
        write (*, '(A,F12.5)') 'Cost Function Increment = ', deltaj
      end if

      if (DeltaJ < config % cost_convergencefactor .and. &
          DeltaJo < zero)  then ! overall is cost getting smaller?
        converged = .true.
        if (config % FullDiagnostics) then
          write (*, '(A,I0)') 'Iteration', iter
          write (*, '(A)') '------------'
          write (*, '(A,L1)') 'Status: converged = ', Converged
          write (*, '(A)') 'New profile:'
          call ufo_geovals_print(geovals, 1)
          call ob % info()
          write (*, '(A)')
          write (*, '(A,3F12.5)') 'Cost Function, increment, cost_convergencefactor = ', &
                                   Jcost, deltaj, config % cost_convergencefactor
        end if 
        exit iterations
      end if

    end if

  end if ! end of specific code for cost test convergence

  ! Iterate (Guess) profile vector
  if (nchans > nprofelements) then
    call fckit_log % debug("Many Chans")
    call ufo_rttovonedvarcheck_NewtonManyChans (Ydiff,          &
                                     nchans,                    &
                                     H_matrix(:,:),             & ! in
                                     transpose (H_matrix(:,:)), & ! in
                                     nprofelements,             &
                                     Diffprofile,               &
                                     b_inv,                     &
                                     r_matrix,                  &
                                     inversionStatus)
  else ! nchans <= nprofelements
    call fckit_log % debug("Few Chans")
    call ufo_rttovonedvarcheck_NewtonFewChans (Ydiff,          &
                                    nchans,                    &
                                    H_matrix(:,:),             & ! in
                                    transpose (H_matrix(:,:)), & ! in
                                    nprofelements,             &
                                    Diffprofile,               &
                                    b_matrix,                  &
                                    r_matrix,                  &
                                    inversionStatus)
  end if

  if (inversionStatus /= 0) then
    inversionStatus = 1
    call fckit_log % warning("inversion failed exiting iterations")
    exit Iterations
  end if

  GuessProfile(:) = BackProfile(:) + Diffprofile(:)

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

  ! Update geovals with guess profile
  call ufo_rttovonedvarcheck_ProfVec2GeoVaLs(geovals, config, profile_index, &
                                             ob, GuessProfile)

  ! if qtotal in retrieval vector check cloud
  ! variables for current iteration

  if ((.NOT. outofRange) .and. profile_index % qt(1) > 0) then

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

  if ((.NOT. outOfRange) .and. (.NOT. config % UseJForConvergence))then
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

! Recalculate final cost - to make sure output when profile has not converged
call ufo_rttovonedvarcheck_CostFunction(Xdiff, b_inv, Ydiff, r_matrix, Jout)
ob % final_cost = Jout(1)
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

!---------------------
! 4. output diagnostics
!---------------------

if (config % UseJForConvergence .and. config % FullDiagnostics) then
  write(*,'(A70,3F10.3,I5,2L5)') "Newton J initial, final, lowest, iter, converged, outofrange = ", &
                                 JCostorig, Jcost,  Jcost, iter, onedvar_success, outOfRange
end if

! ----------
! Tidy up
! ----------
if (allocated(OldProfile))         deallocate(OldProfile)
if (allocated(GuessProfile))       deallocate(GuessProfile)
if (allocated(BackProfile))        deallocate(BackProfile)
if (allocated(H_matrix))           deallocate(H_matrix)
if (allocated(Diffprofile))        deallocate(Diffprofile)
if (allocated(AbsDiffprofile))     deallocate(AbsDiffprofile)
if (allocated(Xdiff))              deallocate(Xdiff)
if (allocated(Ydiff))              deallocate(Ydiff)
if (allocated(Y))                  deallocate(Y)
if (allocated(Y0))                 deallocate(Y0)

call fckit_log % debug("finished with ufo_rttovonedvarcheck_minimize_newton")

end subroutine ufo_rttovonedvarcheck_minimize_newton

!------------------------------------------------------------------------------
!> Update the profile if newber of channels is less than number of elements in 
!! the profile
!!
!! \details Heritage: Ops_SatRad_NewtonFewChans.f90
!!
!! Updates the profile vector DeltaProfile according to Rodgers (1976), Eqn. 101:
!!
!!   x_(n+1) = xb + B.Hn'.Q <br>
!!   Q = U^-1.V <br>
!!   where: <br> x is an atmospheric state vector, subscripted b=background,n=nth
!!           iteration <br>
!!          U = (Hn.B.Hn'+R) <br>
!!          V = (ym-y(xn) - H.(xb-xn)) <br>
!!          ym is the measurement vector (i.e. observed brightness temperatures) <br>
!!          y(xn) is the observation vector calculated for xn <br>
!!          ym and y(xn) are not used individually at all, hence these are input
!!          as a difference vector DeltaBT. <br>
!!          B is the background error covariance matrix <br>
!!          R is the combined forward model and ob error covariance matrix <br>
!!          H is the forward model gradient (w.r.t. xn) matrix <br>
!!          H' is the transpose of H <br>
!!
!!   Q = U^-1.V is solved by Cholesky decomposition.
!!
!! This routine should be used when: <br>
!!     -# The length of the observation vector is less than the length
!!        of the state vector.
!!
!! Note on input/output variable DeltaProfile:
!!
!!   On input, DeltaProfile is x(n-1)-xb. <br>
!!   in construction of variable v, the sign is reversed: <br>
!!   V = (ym-y(xn) + H.(xn-xb)) -- see equation in description above. <br>
!!   On output, DeltaProfile is xn-xb and should be ADDED to the background
!!   profile
!!
!! References:
!!
!!   Rodgers, Retrieval of atmospheric temperature and composition from
!!   remote measurements of thermal radiation, Rev. Geophys.Sp.Phys. 14, 1976.
!!
!!   Rodgers, inverse Methods for Atmospheres: Theory and Practice.  World
!!            Scientific Publishing, 2000.
!!
!! \author Met Office
!!
!! \date 09/06/2020: Created
!!
subroutine ufo_rttovonedvarcheck_NewtonFewChans (DeltaBT,       &
                                      nChans,        &
                                      H_Matrix,      &
                                      H_Matrix_T,    &
                                      nprofelements, &
                                      DeltaProfile,  &
                                      B_matrix,      &
                                      r_matrix,      &
                                      Status)

implicit none

! subroutine arguments:
real(kind_real), intent(in)      :: DeltaBT(:)      !< y-y(x)
integer, intent(in)              :: nChans          !< number of channels
real(kind_real), intent(in)      :: H_Matrix(:,:)   !< Jacobian
real(kind_real), intent(in)      :: H_Matrix_T(:,:) !< (Jacobian)^T
integer, intent(in)              :: nprofelements   !< number of elements in x profile
real(kind_real), intent(inout)   :: DeltaProfile(:) !< x-xb
real(kind_real), intent(in)      :: B_matrix(:,:)   !< state error covariance
type(ufo_rttovonedvarcheck_rsubmatrix), intent(in) :: r_matrix !< observation error covariance
integer, intent(out)              :: Status          !< check if Cholesky decomposition fails

! Local declarations:
character(len=*), parameter :: RoutineName = "ufo_rttovonedvarcheck_NewtonFewChans"
integer                     :: Element
integer                     :: i
real(kind_real)             :: HB(nChans, nprofelements)  ! Scratch matrix
real(kind_real)             :: HBT(nprofelements, nChans) ! Transpose of scratch matrix
real(kind_real)             :: Q(nChans)                  ! Q = U^-1.V
real(kind_real)             :: U(nChans, nChans)          ! U = H.B.H^T + R
real(kind_real)             :: V(nChans)                  ! V = (y-y(x_n))-H^T(xb-x_n)

Status = 0

!---------------------------------------------------------------------------
! 1. Calculate the U and V vectors for the three forms of R matrix allowed
!    for.
!---------------------------------------------------------------------------
HB = matmul(H_matrix, B_matrix)
U = matmul(HB, H_matrix_T)

! Plus sign is not in error - it reverses sign of DeltaProfile to change
! from xn-xb to xb-xn

V = DeltaBT + matmul(H_matrix, DeltaProfile)

!---------------------------------------------------------------------------
! 1.1. Add the R matrix into the U matrix. U = U + R
!---------------------------------------------------------------------------
!U = U + R_Matrix
call r_matrix % add_to_matrix(U, U)

! Calculate Q=(U^-1).V
!------

call Ops_Cholesky (U,      &
                   V,      &
                   nChans, &
                   Q,      &
                   Status)

! Delta profile is (HB)^T.Q
!------
if (Status == 0) then
  HBT = transpose(HB)
  DeltaProfile = matmul(HBT, Q)
end if

end subroutine ufo_rttovonedvarcheck_NewtonFewChans

!------------------------------------------------------------------------------
!> Update the profile if number of channels is more than number of elements in 
!! the profile
!!
!! \details Heritage: Ops_SatRad_NewtonManyChans.f90
!!
!! Updates the profile vector Delta_Profile according to Rodgers (1976), Eqn.
!! 100, extended to allow for additional cost function terms.
!!
!!   x_(n+1) = xb + U^-1.V <br>
!!   where <br> U=(B^-1 + H^T R^-1 H + J2) <br>
!!         V=H^T R^-1 [(ym-y(x_n))+H(x_n-xb)] - J1 <br>
!!   and   <br> J_extra=J0+J1.(x-xb)+(x-xb)^T.J2.(x-xb) is the additional cost
!!         function <br>
!!         x is an atmospheric state vector, subscripted b=background,n=nth
!!         iteration <br>
!!         ym is the measurement vector (i.e. observed brightness temperatures) <br>
!!         y(xn) is the observation vector calculated for xn <br>
!!         ym and y(xn) are not used individually at all, hence these are input
!!         as a difference vector DeltaBT. <br>
!!         B is the background error covariance matrix <br>
!!         R is the combined forward model and ob error covariance matrix <br>
!!         H is the forward model gradient (w.r.t. xn) matrix <br>
!!         H' is the transpose of H <br>
!!
!! When J_extra is zero this is simply Rogers (1976), Eqn. 100.
!!
!!   U^-1.V is solved using Cholesky decomposition.
!!
!! This routine should be used when:
!!     -# The length of the observation vector is greater than the length
!!        of the state vector
!!
!! References:
!!
!!   Rodgers, Retrieval of atmospheric temperature and composition from
!!   remote measurements of thermal radiation, Rev. Geophys.Sp.Phys. 14, 1976.
!!
!!   Rodgers, inverse Methods for Atmospheres: Theory and Practice.  World
!!            Scientific Publishing, 2000.
!!
!! \author Met Office
!!
!! \date 09/06/2020: Created
!!
subroutine ufo_rttovonedvarcheck_NewtonManyChans (DeltaBT, &
                                       nChans,        &
                                       H_Matrix,      &
                                       H_Matrix_T,    &
                                       nprofelements, &
                                       DeltaProfile,  &
                                       B_inverse,     &
                                       r_matrix,      &
                                       Status)

implicit none

! subroutine arguments:
real(kind_real), intent(in)       :: DeltaBT(:)      !< y-y(x)
integer, intent(in)               :: nChans          !< number of channels
real(kind_real), intent(in)       :: H_Matrix(:,:)   !< Jacobian
real(kind_real), intent(in)       :: H_Matrix_T(:,:) !< (Jacobian)^T
integer, intent(in)               :: nprofelements   !< number of elements in profile vector
real(kind_real), intent(inout)    :: DeltaProfile(:) !< x-xb
real(kind_real), intent(in)       :: B_inverse(:,:)  !< inverse state error covariance
type(ufo_rttovonedvarcheck_rsubmatrix), intent(in) :: r_matrix !< observation error covariance
integer, intent(out)              :: Status          !< check if Cholesky decomposition fails

! Local declarations:
character(len=*), parameter :: RoutineName = 'ufo_rttovonedvarcheck_NewtonManyChans'
real(kind_real)             :: HTR(nprofelements, nChans)      ! Scratch vector
real(kind_real)             :: U(nprofelements, nprofelements) ! U = H.B.H^T + R
real(kind_real)             :: V(nprofelements)                ! V = (y-y(x_n))-H^T(xb-x_n)

Status = 0

!---------------------------------------------------------------------------
! 1. Calculate HTR for the three allowed forms of R matrix. if required, the
!    matrix is tested to determine whether it is stored as an inverse and
!    inverted if not.
!---------------------------------------------------------------------------
!HTR = matmul(H_matrix_T, R_inverse)
call r_matrix % multiply_inverse_matrix(H_matrix_T,HTR)

!---------------------------------------------------------------------------
! 2. Calculate U and V
!---------------------------------------------------------------------------

U = matmul(HTR, H_matrix)
V = matmul(HTR, DeltaBT)
V = V + matmul(U, DeltaProfile)

!---------------------------------------------------------------------------
! 3. Add on inverse of B-matrix to U.
!---------------------------------------------------------------------------

U = U + B_inverse

!---------------------------------------------------------------------------
! 5. Calculate new profile increment.
!---------------------------------------------------------------------------

call Ops_Cholesky (U,             &
                   V,             &
                   nprofelements, &
                   DeltaProfile, &
                   Status)

end subroutine ufo_rttovonedvarcheck_NewtonManyChans

!---------------------------------------------------

end module ufo_rttovonedvarcheck_minimize_newton_mod
