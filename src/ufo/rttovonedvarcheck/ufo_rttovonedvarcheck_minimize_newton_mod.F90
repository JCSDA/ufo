! (C) Copyright 2018 UCAR
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

!> Fortran module to provide code shared between nonlinear and tlm/adm radiance calculations

module ufo_rttovonedvarcheck_minimize_newton_mod

use iso_c_binding
use config_mod
use kinds
use ufo_geovals_mod
use ufo_rttovonedvarcheck_utils_mod
use ufo_radiancerttov_tlad_mod

implicit none

private

! public subroutines
public Ops_SatRad_MinimizeNewton_RTTOV12

contains

!-------------------------------------------------------------------------------
! (C) Crown copyright Met Office. All rights reserved.
!     Refer to COPYRIGHT.txt of this distribution for details.
!-------------------------------------------------------------------------------
! Find the most probable atmospheric state vector by minimizing a cost function
! through a series of iterations. if a solution exists, the iterations will
! converge when the iterative increments are acceptably small. A limit on the
! total number of iterations allowed, is imposed.
!
! Using the formulation given by Rodgers (1976) :
!
!   Delta_x = xn + (xb-xn).I' + Wn.(ym-y(xn) - H.(xb-xn))
!   where: x is an atmospheric state vector, subscripted b=background,n=nth
!           iteration
!          I' is a diagonal matrix with I'(J,J) = B_damped(J,J)/B_undamped(J,J)
!           although, however, damping will no longer be used
!          Wn = B.Hn'.(Hn.B.Hn'+R)^-1
!          B is the background error covariance matrix
!          R is the combined forward model and ob error covariance matrix
!
!   Delta_x is checked for convergence after each iteration
!
! The loop is exited with convergence if either of the following conditions are
! true, depending on whether UseJforConvergence is true or false
!   1) if UseJforConvergence is true then the change in total cost function from
!      one iteration to the next is sufficiently small to imply convergence
!   2) if UseJforConvergence is false then the increments to the atmospheric
!      state vector are sufficiently small to imply convergence at an acceptable
!      solution
!   Either of the following two conditions will cause the 1dvar to stop and exit
!   with an error.
!   3) The increments are sufficiently large to suppose a solution will not
!      be found.
!   4) The maximum number of allowed iterations has been reached. in most
!      cases, one of the above criteria will have occurred.
!
! References:
!
!   Rodgers, Retrieval of atmospheric temperature and composition from
!   remote measurements of thermal radiation, Rev. Geophys.Sp.Phys. 14, 1976.
!
!   Eyre, inversion of cloudy satellite sounding radiances by nonlinear
!   optimal estimation. I: Theory and simulation for TOVS,QJ,July 89.
!-------------------------------------------------------------------------------
subroutine Ops_SatRad_MinimizeNewton_RTTOV12(ob_info,       &
                                         r_matrix,      &
                                         r_inv,         &
                                         b_matrix,      &
                                         b_inv,         &
                                         local_geovals, &
                                         profile_index, &
                                         nprofelements, &
                                         conf,          &
                                         obsdb,         &
                                         channels,      &
                                         onedvar_success)

use ufo_rttovonedvarcheck_process_mod, only: &
                    ufo_rttovonedvarcheck_GeoVaLs2ProfVec, &
                    ufo_rttovonedvarcheck_ProfVec2GeoVaLs, &
                    ufo_rttovonedvarcheck_CostFunction, &
                    Ops_SatRad_Qsplit

use ufo_rttovonedvarcheck_forward_model_mod, only: &
                    ufo_rttovonedvarcheck_ForwardModel

implicit none

type(Obinfo_type), intent(in)        :: ob_info
real(kind_real), intent(in)          :: r_matrix(:,:)
real(kind_real), intent(in)          :: r_inv(:,:)
real(kind_real), intent(in)          :: b_matrix(:,:)
real(kind_real), intent(in)          :: b_inv(:,:)
type(ufo_geovals), intent(inout)     :: local_geovals
type(Profileinfo_type), intent(in)   :: profile_index
integer, intent(in)                  :: nprofelements
type(c_ptr), value, intent(in)       :: conf
type(c_ptr), value, intent(in)       :: obsdb
integer(c_int), intent(in)           :: channels(:)
logical, intent(out)                 :: onedvar_success

! Local declarations:
character(len=*), parameter     :: RoutineName = "Ops_SatRad_MinimizeNewton_RTTOV12"
integer                         :: inversionstatus
LOGICAL                         :: outOfRange
LOGICAL                         :: Converged
LOGICAL                         :: Error
LOGICAL                         :: UseJForConvergence
integer                         :: Max1DVarIterations
integer                         :: JConvergenceOption
real                            :: cost_convergencefactor
integer                         :: iter
integer                         :: RTerrorcode
integer                         :: nchans
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
real(kind_real), allocatable    :: Ydiff(:)
real(kind_real), allocatable    :: Y(:)
real(kind_real), allocatable    :: Y0(:)
type(ufo_geovals)               :: geovals
real(kind_real)                 :: Jout(3)

! interface blocks:
!inCLUDE 'Ops_SatRad_CheckIteration.interface'
!inCLUDE 'Ops_SatRad_CheckCloudyIteration.interface'
!inCLUDE 'Ops_SatRad_PrintRTprofile_RTTOV12.interface'

Converged = .FALSE.
onedvar_success = .FALSE.
Error = .FALSE.
UseJForConvergence = .TRUE.
Max1DVarIterations = 10
JConvergenceOption = -1
cost_convergencefactor = 0.0001
nchans = size(channels)

! allocate arrays
allocate(OldProfile(nprofelements))
allocate(GuessProfile(nprofelements))
allocate(BackProfile(nprofelements))
allocate(H_matrix(nchans,nprofelements))
allocate(Diffprofile(nprofelements))
allocate(Ydiff(nchans))
allocate(Y(nchans))
allocate(Y0(nchans))

geovals = local_geovals

Iterations: do iter = 1, Max1DVarIterations

  !-------------------------
  ! 1. Generate new profile
  !-------------------------
  ! Save cost from previous iteration
  if (UseJForConvergence) then
    JcostOld = Jcost
  end if

  ! initialise RTerrorcode and profile increments on first iteration
  if (iter == 1) then

    RTerrorcode = 0
    Diffprofile(:) = 0.0
    JcostOld = 1.0e4

    ! Map GeovaLs to 1D-var profile using B matrix profile structure
    call ufo_rttovonedvarcheck_GeoVaLs2ProfVec(geovals, profile_index, nprofelements, GuessProfile(:))

  end if

  ! Save current profile
  OldProfile(:) = GuessProfile(:)

  ! call forward model to generate jacobian
  call ufo_rttovonedvarcheck_ForwardModel(geovals, ob_info, obsdb, &
                                       channels(:), conf, &
                                       profile_index, GuessProfile(:), &
                                       Y(:), H_matrix)
  write(*,*) "Observed BTs After bias correction: ",ob_info%yobs(:)
  write(*,*) "RTTOV BTs: = ",Y(:)

  if (iter == 1) then
    BackProfile(:) = GuessProfile(:)
    Y0(:) = Y(:)
  end if

  ! exit on error
  if (RTerrorcode /= 0) then
    write(*,*) "Radiative transfer error"
    exit Iterations
  end if

  !-----------------------------------------------------
  ! 1a. if UseJForConvergence is true
  !     then compute costfunction for latest profile
  !     and determine convergence using change in cost fn
  !-----------------------------------------------------

  if (UseJForConvergence) then

    Diffprofile(:) = GuessProfile(:) - BackProfile(:)
    Ydiff(:) = ob_info%yobs(:) - Y(:)
    call ufo_rttovonedvarcheck_CostFunction(Diffprofile, b_inv, Ydiff, r_inv, Jout)
    Jcost = Jout(1)

    ! exit on error
    !if (inversionStatus /= 0) exit Iterations

    ! store initial cost value
    if (iter == 1) JCostOrig = jcost

    write(*,*) "iter,Jcost,JcostOld,JCostorig = ",iter,Jcost,JcostOld,JCostorig

    ! check for convergence
    if (iter > 1) then

      if (JConvergenceOption == 1) then

        ! percentage change tested between iterations
        DeltaJ = ABS ((Jcost - JcostOld) / MAX (Jcost, TinY (0.0)))

        ! default test for checking that overall cost is getting smaller
        DeltaJo = -1.0

      ELSE

        ! absolute change tested between iterations
        DeltaJ = ABS (Jcost - JcostOld)

        ! change between current cost and initial
        DeltaJo = Jcost - JCostorig

      end if

!      if (SatRad_FullDiagnostics) then
        write (*, '(A,F12.5)') 'Cost Function=', Jcost
        write (*, '(A,F12.5)') 'Cost Function increment=', deltaj
        write (*, '(A,F12.5)') 'cost_convergencefactor=', cost_convergencefactor
!      end if

      if (DeltaJ < cost_convergencefactor .AND. &
          DeltaJo < 0.0)  then ! overall is cost getting smaller?
        converged = .TRUE.
        onedvar_success = .TRUE.

!        if (SatRad_FullDiagnostics) then
!          write (*, '(A,I0)') 'Iteration', iter
!          write (*, '(A)') '------------'
!          write (*, '(A,L1)') 'Status: converged = ', Converged
!          write (*, '(A)') 'New profile:'
!          CALL Ops_SatRad_PrintRTprofile_RTTOV12 (RTprof_Guess)
!          write (*, '(A)')
!        end if

        exit iterations
      end if

    end if

  end if ! end of specific code for cost test convergence

  ! Iterate (Guess) profile vector
  if (nchans > nprofelements) then
    write(*,*) "Many Chans"
    CALL Ops_SatRad_NewtonManyChans (Ydiff,                     &
                                     nchans,                    &
                                     H_matrix(:,:),             & ! in
                                     transpose (H_matrix(:,:)), & ! in
                                     nprofelements,             &
                                     Diffprofile,               &
                                     b_inv,                     &
                                     r_matrix,                  &
                                     inversionStatus)
  ELSE ! nchans <= nprofelements
    write(*,*) "Few Chans"
    CALL Ops_SatRad_NewtonFewChans (Ydiff,                     &
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
    write(*,*) "inversion failed"
    exit Iterations
  end if

  GuessProfile(:) = BackProfile(:) + Diffprofile(:)

!  !---------------------------------------------------------
!  ! 2. Check new profile and transfer to forward model format
!  !---------------------------------------------------------
!
!  ! Check profile and constrain humidity variables
!
!  CALL Ops_SatRad_CheckIteration (RTprof_Guess % rttov12_profile % s2m % p,            & ! in
!                                  RTprof_Guess % rttov12_profile % s2m % t,            & ! in
!                                  RTprof_Guess % rttov12_profile % t(:), & ! in
!                                  GuessProfile(:),                                     & ! inout
!                                  outOfRange)                                            ! out
!
  ! Update RT-format guess profile
  call ufo_rttovonedvarcheck_ProfVec2GeoVaLs(geovals, profile_index, nprofelements, GuessProfile)
  
!  ! if qtotal in retrieval vector check cloud
!  ! variables for current iteration
!
!  if ((.NOT. outofRange) .AND. &
!      profindex % qt(1) > 0) then
!
!    if (iter >= IterNumForLWPCheck) then
!    
!      if (RTTOV_mwscattSwitch) then
!      
!        !cloud information is from the scatt profile
!        if (RTprof_Guess % rttov12_profile_scatt % use_totalice) then
!        
!          CALL Ops_SatRad_CheckCloudyIteration (RTprof_Guess % rttov12_profile_scatt % clw(:),      & ! in
!                                                RTprof_Guess % rttov12_profile_scatt % totalice(:), & ! in
!                                                outOfRange)                                                         ! out
!        ELSE
!        
!          CALL Ops_SatRad_CheckCloudyIteration (RTprof_Guess % rttov12_profile_scatt % clw(:), & ! in
!                                                RTprof_Guess % rttov12_profile_scatt % ciw(:), & ! in
!                                                outOfRange)
!        
!        end if
!                                                                                               
!      ELSE
!      
!        !clear air profile is used for clw and cloudice diagnostic 
!        CALL Ops_SatRad_CheckCloudyIteration (RTprof_Guess % rttov12_profile % clw(:), & ! in
!                                              CloudIce(:),                                           & ! in
!                                              outOfRange)                                              ! out
!      
!      end if                                       
!    
!    end if                                                                                  
!
!  end if

  !-------------------------------------------------
  ! 3. Check for convergence using change in profile
  !    This is performed if UseJforConvergence is false
  !-------------------------------------------------

!  AbsDiffprofile(:) = 0.0
!
!  if ((.NOT. outOfRange) .AND. &
!      (.NOT. UseJForConvergence))then
!    ABSDiffProfile(:) = ABS (GuessProfile(:) - OldProfile(:))
!    if (ALL (AbsDiffProfile(:) <= B_sigma(:) * ConvergenceFactor)) then
!      Converged = .TRUE.
!    end if
!  end if

  !---------------------
  ! 4. output diagnostics
  !---------------------

!  if (SatRad_FullDiagnostics) then
!    write (*, '(A,I0)') 'Iteration', iter
!    write (*, '(A)') '------------'
!    write (*, '(A,L1)') 'Status: converged = ', Converged
!    if (outOfRange) write (*, '(A)') 'exiting with bad increments'
!    write (*, '(A)') 'New profile:'
!    CALL Ops_SatRad_PrintRTprofile_RTTOV12 (RTprof_Guess)
!    write (*, '(A)')
!  end if

  ! exit conditions

  if (Converged .OR. outOfRange) exit iterations

end do Iterations

write(*,*) "----------------------------"
write(*,*) "Starting cost = ",JCostorig
write(*,*) "Final cost = ",Jcost
write(*,*) "Converged? ", Converged
write(*,*) "out Profile = ",GuessProfile(:)
write(*,*) "----------------------------"

!Ob % Niter = iter
!if (RTerrorcode /= 0 .OR. .NOT. Converged .OR. outOfRange .OR. &
!    inversionStatus /= 0) then
!  Error = .TRUE.
!end if
!
!if (.NOT. Error .AND. UseJForConvergence) then
!  ! store final cost function and retrieved bts
!  Ob % Jcost = Jcost
!  Ob % Britemp(Channels_1dvar(:)) = Britemp(:)
!end if

end subroutine Ops_SatRad_MinimizeNewton_RTTOV12

!-------------------------------------------------------------------------------
! (C) Crown copyright Met Office. All rights reserved.
!     Refer to COPYRIGHT.txt of this distribution for details.
!-------------------------------------------------------------------------------
! Updates the profile vector DeltaProfile according to Rodgers (1976), Eqn. 101:
!
!   x_(n+1) = xb + B.Hn'.Q
!   Q = U^-1.V
!   where: x is an atmospheric state vector, subscripted b=background,n=nth
!           iteration
!          U = (Hn.B.Hn'+R)
!          V = (ym-y(xn) - H.(xb-xn))
!          ym is the measurement vector (i.e. observed brightness temperatures)
!          y(xn) is the observation vector calculated for xn
!          ym and y(xn) are not used individually at all, hence these are input
!          as a difference vector DeltaBT.
!          B is the background error covariance matrix
!          R is the combined forward model and ob error covariance matrix
!          H is the forward model gradient (w.r.t. xn) matrix
!          H' is the transpose of H
!
!   Q = U^-1.V is solved by Cholesky decomposition.
!
! This routine should be used when:
!     1) The length of the observation vector is less than the length
!        of the state vector,
!     2) Where no additional cost function terms are provided and
!     3) where Newtonian minimisation is desired.
!
! Note on input/output variable DeltaProfile:
!
!   On input, DeltaProfile is x(n-1)-xb.
!   in construction of variable v, the sign is reversed:
!   V = (ym-y(xn) + H.(xn-xb)) -- see equation in description above.
!   On output, DeltaProfile is xn-xb and should be ADDED to the background
!   profile
!
!
! References:
!
!   Rodgers, Retrieval of atmospheric temperature and composition from
!   remote measurements of thermal radiation, Rev. Geophys.Sp.Phys. 14, 1976.
!
!   Rodgers, inverse Methods for Atmospheres: Theory and Practice.  World
!            Scientific Publishing, 2000.
!-------------------------------------------------------------------------------

subroutine Ops_SatRad_NewtonFewChans (DeltaBT,       &
                                      nChans,        &
                                      H_Matrix,      &
                                      H_Matrix_T,    &
                                      nprofelements, &
                                      DeltaProfile,  &
                                      B_matrix,      &
                                      R_matrix,      &
                                      Status)

implicit none

! subroutine arguments:
real(kind_real), intent(in)     :: DeltaBT(:)        ! y-y(x)
integer, intent(in)             :: nChans
real(kind_real), intent(in)     :: H_Matrix(:,:)     ! Jacobian
real(kind_real), intent(in)     :: H_Matrix_T(:,:)   ! (Jacobian)^T
integer, intent(in)             :: nprofelements
real(kind_real), intent(inout)  :: DeltaProfile(:)   ! see note in header
real(kind_real), intent(in)     :: B_matrix(:,:)
real(kind_real), intent(in)     :: R_matrix(:,:)
integer, intent(out)            :: Status

! Local declarations:
character(len=*), parameter :: RoutineName = "Ops_SatRad_Minimize_101"
integer                     :: Element
integer                     :: i
real(kind_real)             :: HB(nChans,nprofelements) ! Scratch vector
real(kind_real)             :: Q(nChans)                ! Q = U^-1.V
real(kind_real)             :: U(nChans,nChans)         ! U = H.B.H^T + R
real(kind_real)             :: V(nChans)                ! V = (y-y(x_n))-H^T(xb-x_n)

Status = 0

!---------------------------------------------------------------------------
! 1. Calculate the U and V vectors for the three forms of R matrix allowed
!    for.
!---------------------------------------------------------------------------

write(*,*) "B_matrix shape = ",shape(B_matrix)
write(*,*) "H_matrix shape = ",shape(H_matrix)
write(*,*) "H_matrix_T shape = ",shape(H_matrix_T)
write(*,*) "HB shape = ",shape(HB)

HB = matmul(H_matrix, B_matrix)
U = matmul(HB, H_matrix_T)

! Plus sign is not in error - it reverses sign of DeltaProfile to change
! from xn-xb to xb-xn

V = DeltaBT + matmul(H_matrix, DeltaProfile)

!---------------------------------------------------------------------------
! 1.1. Add the R matrix into the U matrix.
!---------------------------------------------------------------------------
U = U + R_Matrix

! Calculate Q=(U^-1).V
!------

CALL Ops_Cholesky (U,      &
                   V,      &
                   nChans, &
                   Q,      &
                   Status)
if (Status /= 0) goto 9999

! Delta profile is (HB)^T.Q
!------

DeltaProfile = matmul(transpose(HB), Q)

9999 continue

end subroutine Ops_SatRad_NewtonFewChans

!-------------------------------------------------------------------------------
! (C) Crown copyright Met Office. All rights reserved.
!     Refer to COPYRIGHT.txt of this distribution for details.
!-------------------------------------------------------------------------------
! Updates the profile vector Delta_Profile according to Rodgers (1976), Eqn.
! 100, extended to allow for additional cost function terms.
!
!   x_(n+1) = xb + U^-1.V
!   where U=(B^-1 + H^T R^-1 H + J2)
!         V=H^T R^-1 [(ym-y(x_n))+H(x_n-xb)] - J1
!   and   J_extra=J0+J1.(x-xb)+(x-xb)^T.J2.(x-xb) is the additional cost
!         function
!         x is an atmospheric state vector, subscripted b=background,n=nth
!         iteration
!         ym is the measurement vector (i.e. observed brightness temperatures)
!         y(xn) is the observation vector calculated for xn
!         ym and y(xn) are not used individually at all, hence these are input
!         as a difference vector DeltaBT.
!         B is the background error covariance matrix
!         R is the combined forward model and ob error covariance matrix
!         H is the forward model gradient (w.r.t. xn) matrix
!         H' is the transpose of H
!
! When J_extra is zero this is simply Rogers (1976), Eqn. 100.
!
!   U^-1.V is solved using Cholesky decomposition.
!
! This routine should be used when:
!     1) The length of the observation vector is greater than the length
!        of the state vector,
!     2) Additional cost function terms are provided and
!     3) where Newtonian minimisation is desired.
!
! References:
!
!   Rodgers, Retrieval of atmospheric temperature and composition from
!   remote measurements of thermal radiation, Rev. Geophys.Sp.Phys. 14, 1976.
!
!   Rodgers, inverse Methods for Atmospheres: Theory and Practice.  World
!            Scientific Publishing, 2000.
!-------------------------------------------------------------------------------

subroutine Ops_SatRad_NewtonManyChans (DeltaBT,       &
                                       nChans,        &
                                       H_Matrix,      &
                                       H_Matrix_T,    &
                                       nprofelements, &
                                       DeltaProfile,  &
                                       B_inverse,     &
                                       R_matrix,      &
                                       Status)

implicit none

! subroutine arguments:
real(kind_real), intent(in)     :: DeltaBT(:)        ! y-y(x)
integer, intent(in)             :: nChans
real(kind_real), intent(in)     :: H_Matrix(:,:)     ! Jacobian
real(kind_real), intent(in)     :: H_Matrix_T(:,:)   ! (Jacobian)^T
integer, intent(in)             :: nprofelements
real(kind_real), intent(inout)  :: DeltaProfile(:)   ! see note in header
real(kind_real), intent(in)     :: B_inverse(:,:)
real(kind_real), intent(in)     :: R_matrix(:,:)
integer, intent(out)            :: Status

! Local declarations:
character(len=*), parameter :: RoutineName = 'Ops_SatRad_NewtonManyChans'
real(kind_real)             :: HTR(nprofelements, nChans)              ! Scratch vector
real(kind_real)             :: U(nprofelements, nprofelements)         ! U = H.B.H^T + R
real(kind_real)             :: V(nprofelements)                        ! V = (y-y(x_n))-H^T(xb-x_n)

Status = 0

!---------------------------------------------------------------------------
! 1. Calculate HTR for the three allowed forms of R matrix. if required, the
!    matrix is tested to determine whether it is stored as an inverse and
!    inverted if not.
!---------------------------------------------------------------------------
HTR = matmul(H_matrix_T, R_matrix(:,:))

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

CALL Ops_Cholesky (U,             &
                   V,             &
                   nprofelements, &
                   DeltaProfile, &
                   Status)

end subroutine Ops_SatRad_NewtonManyChans

!-------------------------------------------------------------------------------
! (C) Crown copyright Met Office. All rights reserved.
!     Refer to COPYRIGHT.txt of this distribution for details.
!-------------------------------------------------------------------------------
! Solves the Linear equation UQ=V for Q where U is a symmetric positive definite
! matrix and U and Q are vectors of length N.  The method follows that in Golub
! and Van Loan although this is pretty standard.
!
! if U is not positive definite this will be detected by the program and flagged
! as an error.  U is assumed to be symmetric as only the upper triangle is in
! fact used.
!-------------------------------------------------------------------------------

subroutine Ops_Cholesky (U,         &
                         V,         &
                         N,         &
                         Q,         &
                         ErrorCode)

implicit none

! subroutine arguments:
integer, intent(in)          :: n
real(kind_real), intent(in)  :: U(n,n)
real(kind_real), intent(in)  :: V(n)
real(kind_real), intent(out) :: Q(n)
integer, intent(out)         :: ErrorCode

! Local declarations:
character(len=*), parameter  :: RoutineName = "Ops_Cholesky"
real(kind_real), parameter   :: Tolerance = TinY (0.0) * 100.0
character(len=80)            :: ErrorMessage
integer                      :: j
integer                      :: k
real(kind_real)              :: G(n,n)   ! The Cholesky Triangle Matrix
real(kind_real)              :: X(n)     ! Temporary array used in calculating G

ErrorCode = 0

! Determine the Cholesky triangle matrix.

do j = 1, n
  X(j:n) = U(j:n,j)
  if (j /= 1) then
    do k = 1, j - 1
      X(j:n) = X(j:n) - G(j,k) * G(j:n,k)
    end do
  end if
  if (X(j) <= Tolerance) then
    ErrorCode = 1
    Errormessage = ' :U matrix is not positive definite'
    write(*,*) RoutineName,ErrorMessage
    goto 9999
  end if
  G(J:N,J) = X(J:N) / sqrt (X(J))
end do

! Solve Gx=v for x by forward substitution

X = V
X(1) = X(1) / G(1,1)
do j = 2, n
  X(j) = (X(j) - dot_product(G(j,1:j - 1), X(1:j - 1))) / G(j,j)
end do

! Solve G^T.q=x for q by backward substitution

Q = x
Q(n) = Q(n) / G(n,n)
do j = n - 1, 1, -1
  Q(j) = (Q(j) - dot_product(G(j + 1:n,j), Q(j + 1:n))) / G(j,j)
end do

9999 continue

end subroutine Ops_Cholesky

end module ufo_rttovonedvarcheck_minimize_newton_mod
