! (C) Copyright 2018 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

!> Fortran module to provide code shared between nonlinear and tlm/adm radiance calculations

module ufo_onedvarfortran_minimize_newton_mod

use iso_c_binding
use config_mod
use kinds
use ufo_geovals_mod
use ufo_onedvarfortran_utils_mod
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
! through a series of iterations. If a solution exists, the iterations will
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
!   1) If UseJforConvergence is true then the change in total cost function from
!      one iteration to the next is sufficiently small to imply convergence
!   2) If UseJforConvergence is false then the increments to the atmospheric
!      state vector are sufficiently small to imply convergence at an acceptable
!      solution
!   Either of the following two conditions will cause the 1dvar to stop and exit
!   with an error.
!   3) The increments are sufficiently large to suppose a solution will not
!      be found.
!   4) The maximum number of allowed iterations has been reached. In most
!      cases, one of the above criteria will have occurred.
!
! References:
!
!   Rodgers, Retrieval of atmospheric temperature and composition from
!   remote measurements of thermal radiation, Rev. Geophys.Sp.Phys. 14, 1976.
!
!   Eyre, Inversion of cloudy satellite sounding radiances by nonlinear
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

use ufo_onedvarfortran_process_mod, only: &
                    ufo_onedvarfortran_GeoVaLs2ProfVec, &
                    ufo_onedvarfortran_ProfVec2GeoVaLs, &
                    ufo_onedvarfortran_CostFunction, &
                    Ops_SatRad_Qsplit

use ufo_onedvarfortran_forward_model_mod, only: &
                    ufo_onedvarfortran_ForwardModel

implicit none

type(Obinfo_type), intent(in)        :: ob_info
real(kind_real), intent(in)          :: r_matrix(:,:)
real(kind_real), intent(in)          :: r_inv(:,:)
real(kind_real), intent(in)          :: b_matrix(:,:)
real(kind_real), intent(in)          :: b_inv(:,:)
type(ufo_geovals), intent(inout)     :: local_geovals
type(Profileinfo_type), intent(in)   :: profile_index
integer, intent(in)                  :: nprofelements
type(c_ptr), VALUE, intent(in)       :: conf
type(c_ptr), VALUE, intent(in)       :: obsdb
integer(c_int), intent(in)           :: channels(:)
logical, intent(out)                 :: onedvar_success

! Local declarations:
CHARACTER(len=*), PARAMETER     :: RoutineName = "Ops_SatRad_MinimizeNewton_RTTOV12"
INTEGER                         :: Inversionstatus
LOGICAL                         :: OutOfRange
LOGICAL                         :: Converged
LOGICAL                         :: Error
LOGICAL                         :: UseJForConvergence
INTEGER                         :: Max1DVarIterations
INTEGER                         :: JConvergenceOption
REAL                            :: cost_convergencefactor
INTEGER                         :: iter
INTEGER                         :: RTerrorcode
INTEGER                         :: nchans
REAL(kind_real)                 :: Jcost     ! current value
REAL(kind_real)                 :: JcostOld  ! previous iteration value
REAL(kind_real)                 :: JcostOrig ! initial value
REAL(kind_real)                 :: DeltaJ
REAL(kind_real)                 :: DeltaJo

REAL(kind_real), allocatable    :: OldProfile(:)
REAL(kind_real), allocatable    :: GuessProfile(:)
REAL(kind_real), allocatable    :: BackProfile(:)
REAL(kind_real), allocatable    :: H_matrix(:,:)
REAL(kind_real), allocatable    :: Diffprofile(:)
REAL(kind_real), allocatable    :: Ydiff(:)
REAL(kind_real), allocatable    :: Y(:)
REAL(kind_real), allocatable    :: Y0(:)
type(ufo_geovals)               :: geovals
real(kind_real)                 :: Jout(3)

! Interface blocks:
!INCLUDE 'Ops_SatRad_CheckIteration.interface'
!INCLUDE 'Ops_SatRad_CheckCloudyIteration.interface'
!INCLUDE 'Ops_SatRad_PrintRTprofile_RTTOV12.interface'

Converged = .FALSE.
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

Iterations: DO iter = 1, Max1DVarIterations

  !-------------------------
  ! 1. Generate new profile
  !-------------------------

  ! Initialise RTerrorcode and profile increments on first iteration
  IF (iter == 1) THEN

    RTerrorcode = 0
    Diffprofile(:) = 0.0
    JcostOld = 1.0e4

    ! Map GeovaLs to 1D-var profile using B matrix profile structure
    call ufo_onedvarfortran_GeoVaLs2ProfVec(geovals, profile_index, nprofelements, GuessProfile(:))
    OldProfile(:) = GuessProfile(:)

  END IF

  ! call forward model to generate jacobian
  call ufo_onedvarfortran_ForwardModel(geovals, ob_info, obsdb, &
                                       channels(:), conf, &
                                       profile_index, GuessProfile(:), &
                                       Y(:), H_matrix)
  write(*,*) "From forward model Y(:) = ",Y(:)

  IF (iter == 1) THEN
    BackProfile(:) = GuessProfile(:)
    Y0(:) = Y(:)
  END IF

  ! Exit on error
  IF (RTerrorcode /= 0) THEN
    write(*,*) "Radiative transfer error"
    EXIT Iterations
  END IF

  !-----------------------------------------------------
  ! 1a. If UseJForConvergence is true
  !     then compute costfunction for latest profile
  !     and determine convergence using change in cost fn
  !-----------------------------------------------------

  IF (UseJForConvergence) THEN

    ! store cost function from previous cycle
    IF (UseJForConvergence) THEN
      JcostOld = Jcost
    END IF

    Diffprofile(:) = GuessProfile(:) - BackProfile(:)
    write(*,*) "Diffprofile = ",Diffprofile(:)
    write(*,*) "GuessProfile = ",GuessProfile(:)
    write(*,*) "BackProfile(:) = ",BackProfile(:)
    Ydiff(:) = ob_info%yobs(:) - Y(:)
    call ufo_onedvarfortran_CostFunction(Diffprofile, b_inv, Ydiff, r_inv, Jout)
    Jcost = Jout(1)

    ! Exit on error
    !IF (InversionStatus /= 0) EXIT Iterations

    ! store initial cost value
    IF (iter == 1) JCostOrig = jcost

    ! check for convergence
    IF (iter > 1) THEN

      IF (JConvergenceOption == 1) THEN

        ! percentage change tested between iterations
        DeltaJ = ABS ((Jcost - JcostOld) / MAX (Jcost, TINY (0.0)))

        ! default test for checking that overall cost is getting smaller
        DeltaJo = -1.0

      ELSE

        ! absolute change tested between iterations
        DeltaJ = ABS (Jcost - JcostOld)

        ! change between current cost and initial
        DeltaJo = Jcost - JCostorig

      END IF

!      IF (SatRad_FullDiagnostics) THEN
        WRITE (*, '(A,F12.5)') 'Cost Function=', Jcost
        WRITE (*, '(A,F12.5)') 'Cost Function Increment=', deltaj
        WRITE (*, '(A,F12.5)') 'cost_convergencefactor=', cost_convergencefactor
!      END IF

      IF (DeltaJ < cost_convergencefactor .AND. &
          DeltaJo < 0.0)  THEN ! overall is cost getting smaller?
        converged = .TRUE.

!        IF (SatRad_FullDiagnostics) THEN
!          WRITE (*, '(A,I0)') 'Iteration', iter
!          WRITE (*, '(A)') '------------'
!          WRITE (*, '(A,L1)') 'Status: converged = ', Converged
!          WRITE (*, '(A)') 'New profile:'
!          CALL Ops_SatRad_PrintRTprofile_RTTOV12 (RTprof_Guess)
!          WRITE (*, '(A)')
!        END IF

        EXIT iterations
      END IF

    END IF

  END IF ! end of specific code for cost test convergence

  ! Iterate (Guess) profile vector
  IF (nchans > nprofelements) THEN
    write(*,*) "Many Chans"
    CALL Ops_SatRad_NewtonManyChans (Ydiff,                     &
                                     nchans,                    &
                                     H_matrix(:,:),             & ! in
                                     TRANSPOSE (H_matrix(:,:)), & ! in
                                     nprofelements,             &
                                     Diffprofile,               &
                                     b_inv,                     &
                                     r_matrix,                  &
                                     InversionStatus)
  ELSE ! nchans <= nprofelements
    write(*,*) "Few Chans"
    CALL Ops_SatRad_NewtonFewChans (Ydiff,                     &
                                    nchans,                    &
                                    H_matrix(:,:),             & ! in
                                    TRANSPOSE (H_matrix(:,:)), & ! in
                                    nprofelements,             &
                                    Diffprofile,               &
                                    b_matrix,                  &
                                    r_matrix,                  &
                                    InversionStatus)
  END IF

  IF (InversionStatus /= 0) THEN
    InversionStatus = 1
    write(*,*) "Inversion failed"
    EXIT Iterations
  END IF

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
!                                  OutOfRange)                                            ! out
!
  ! Update RT-format guess profile
  call ufo_onedvarfortran_ProfVec2GeoVaLs(geovals, profile_index, nprofelements, GuessProfile)
  
!  ! If qtotal in retrieval vector check cloud
!  ! variables for current iteration
!
!  IF ((.NOT. OutofRange) .AND. &
!      profindex % qt(1) > 0) THEN
!
!    IF (iter >= IterNumForLWPCheck) THEN
!    
!      IF (RTTOV_mwscattSwitch) THEN
!      
!        !cloud information is from the scatt profile
!        IF (RTprof_Guess % rttov12_profile_scatt % use_totalice) THEN
!        
!          CALL Ops_SatRad_CheckCloudyIteration (RTprof_Guess % rttov12_profile_scatt % clw(:),      & ! in
!                                                RTprof_Guess % rttov12_profile_scatt % totalice(:), & ! in
!                                                OutOfRange)                                                         ! out
!        ELSE
!        
!          CALL Ops_SatRad_CheckCloudyIteration (RTprof_Guess % rttov12_profile_scatt % clw(:), & ! in
!                                                RTprof_Guess % rttov12_profile_scatt % ciw(:), & ! in
!                                                OutOfRange)
!        
!        END IF
!                                                                                               
!      ELSE
!      
!        !clear air profile is used for clw and cloudice diagnostic 
!        CALL Ops_SatRad_CheckCloudyIteration (RTprof_Guess % rttov12_profile % clw(:), & ! in
!                                              CloudIce(:),                                           & ! in
!                                              OutOfRange)                                              ! out
!      
!      END IF                                       
!    
!    END IF                                                                                  
!
!  END IF

  !-------------------------------------------------
  ! 3. Check for convergence using change in profile
  !    This is performed if UseJforConvergence is false
  !-------------------------------------------------

!  AbsDiffprofile(:) = 0.0
!
!  IF ((.NOT. OutOfRange) .AND. &
!      (.NOT. UseJForConvergence))THEN
!    ABSDiffProfile(:) = ABS (GuessProfile(:) - OldProfile(:))
!    IF (ALL (AbsDiffProfile(:) <= B_sigma(:) * ConvergenceFactor)) THEN
!      Converged = .TRUE.
!    END IF
!  END IF

  !---------------------
  ! 4. Output diagnostics
  !---------------------

!  IF (SatRad_FullDiagnostics) THEN
!    WRITE (*, '(A,I0)') 'Iteration', iter
!    WRITE (*, '(A)') '------------'
!    WRITE (*, '(A,L1)') 'Status: converged = ', Converged
!    IF (OutOfRange) WRITE (*, '(A)') 'Exiting with bad increments'
!    WRITE (*, '(A)') 'New profile:'
!    CALL Ops_SatRad_PrintRTprofile_RTTOV12 (RTprof_Guess)
!    WRITE (*, '(A)')
!  END IF

  ! Exit conditions

  IF (Converged .OR. OutOfRange) EXIT iterations

END DO Iterations

write(*,*) "----------------------------"
write(*,*) "Starting cost = ",JCostorig
write(*,*) "Final cost = ",Jcost
write(*,*) "Converged? ", Converged
write(*,*) "Out Profile = ",GuessProfile(:)
write(*,*) "----------------------------"

!Ob % Niter = iter
!IF (RTerrorcode /= 0 .OR. .NOT. Converged .OR. OutOfRange .OR. &
!    InversionStatus /= 0) THEN
!  Error = .TRUE.
!END IF
!
!IF (.NOT. Error .AND. UseJForConvergence) THEN
!  ! store final cost function and retrieved bts
!  Ob % Jcost = Jcost
!  Ob % Britemp(Channels_1dvar(:)) = Britemp(:)
!END IF

END SUBROUTINE Ops_SatRad_MinimizeNewton_RTTOV12

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
!          as a dIfference vector DeltaBT.
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
!   In construction of variable v, the sign is reversed:
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
!   Rodgers, Inverse Methods for Atmospheres: Theory and Practice.  World
!            ScientIfic Publishing, 2000.
!-------------------------------------------------------------------------------

SUBROUTINE Ops_SatRad_NewtonFewChans (DeltaBT,       &
                                      nChans,        &
                                      H_Matrix,      &
                                      H_Matrix_T,    &
                                      nprofelements, &
                                      DeltaProfile,  &
                                      B_matrix,      &
                                      R_matrix,      &
                                      Status)

IMPLICIT NONE

! Subroutine arguments:
REAL(kind_real), INTENT(IN)     :: DeltaBT(:)        ! y-y(x)
INTEGER, INTENT(IN)             :: nChans
REAL(kind_real), INTENT(IN)     :: H_Matrix(:,:)     ! Jacobian
REAL(kind_real), INTENT(IN)     :: H_Matrix_T(:,:)   ! (Jacobian)^T
INTEGER, INTENT(IN)             :: nprofelements
REAL(kind_real), INTENT(INOUT)  :: DeltaProfile(:)   ! see note in header
REAL(kind_real), INTENT(IN)     :: B_matrix(:,:)
REAL(kind_real), INTENT(IN)     :: R_matrix(:,:)
INTEGER, INTENT(OUT)            :: Status

! Local declarations:
CHARACTER(len=*), PARAMETER :: RoutineName = "Ops_SatRad_Minimize_101"
INTEGER                     :: Element
INTEGER                     :: i
REAL(kind_real)             :: HB(nChans,nprofelements) ! Scratch vector
REAL(kind_real)             :: Q(nChans)                ! Q = U^-1.V
REAL(kind_real)             :: U(nChans,nChans)         ! U = H.B.H^T + R
REAL(kind_real)             :: V(nChans)                ! V = (y-y(x_n))-H^T(xb-x_n)

Status = 0

!---------------------------------------------------------------------------
! 1. Calculate the U and V vectors for the three forms of R matrix allowed
!    for.
!---------------------------------------------------------------------------

write(*,*) "B_matrix shape = ",shape(B_matrix)
write(*,*) "H_matrix shape = ",shape(H_matrix)
write(*,*) "H_matrix_T shape = ",shape(H_matrix_T)
write(*,*) "HB shape = ",shape(HB)

HB = MATMUL (H_matrix, B_matrix)
U = MATMUL (HB, H_matrix_T)

! Plus sign is not in error - it reverses sign of DeltaProfile to change
! from xn-xb to xb-xn

V = DeltaBT + MATMUL (H_matrix, DeltaProfile)

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
IF (Status /= 0) GOTO 9999

! Delta profile is (HB)^T.Q
!------

DeltaProfile = MATMUL (TRANSPOSE (HB), Q)

9999 CONTINUE

END SUBROUTINE Ops_SatRad_NewtonFewChans

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
!   Rodgers, Inverse Methods for Atmospheres: Theory and Practice.  World
!            Scientific Publishing, 2000.
!-------------------------------------------------------------------------------

SUBROUTINE Ops_SatRad_NewtonManyChans (DeltaBT,       &
                                       nChans,        &
                                       H_Matrix,      &
                                       H_Matrix_T,    &
                                       nprofelements, &
                                       DeltaProfile,  &
                                       B_Inverse,     &
                                       R_matrix,      &
                                       Status)

IMPLICIT NONE

! Subroutine arguments:
REAL(kind_real), INTENT(IN)     :: DeltaBT(:)        ! y-y(x)
INTEGER, INTENT(IN)             :: nChans
REAL(kind_real), INTENT(IN)     :: H_Matrix(:,:)     ! Jacobian
REAL(kind_real), INTENT(IN)     :: H_Matrix_T(:,:)   ! (Jacobian)^T
INTEGER, INTENT(IN)             :: nprofelements
REAL(kind_real), INTENT(INOUT)  :: DeltaProfile(:)   ! see note in header
REAL(kind_real), INTENT(IN)     :: B_Inverse(:,:)
REAL(kind_real), INTENT(IN)     :: R_matrix(:,:)
INTEGER, INTENT(OUT)            :: Status

! Local declarations:
CHARACTER(len=*), PARAMETER :: RoutineName = 'Ops_SatRad_NewtonManyChans'
REAL(kind_real)             :: HTR(nprofelements, nChans)              ! Scratch vector
REAL(kind_real)             :: U(nprofelements, nprofelements)         ! U = H.B.H^T + R
REAL(kind_real)             :: V(nprofelements)                        ! V = (y-y(x_n))-H^T(xb-x_n)

Status = 0

!---------------------------------------------------------------------------
! 1. Calculate HTR for the three allowed forms of R matrix. If required, the
!    matrix is tested to determine whether it is stored as an inverse and
!    inverted if not.
!---------------------------------------------------------------------------
HTR = MATMUL (H_matrix_T, R_matrix(:,:))

!---------------------------------------------------------------------------
! 2. Calculate U and V
!---------------------------------------------------------------------------

U = MATMUL (HTR, H_matrix)
V = MATMUL (HTR, DeltaBT)
V = V + MATMUL (U, DeltaProfile)

!---------------------------------------------------------------------------
! 3. Add on inverse of B-matrix to U.
!---------------------------------------------------------------------------

U = U + B_Inverse

!---------------------------------------------------------------------------
! 5. Calculate new profile increment.
!---------------------------------------------------------------------------

CALL Ops_Cholesky (U,             &
                   V,             &
                   nprofelements, &
                   DeltaProfile, &
                   Status)

END SUBROUTINE Ops_SatRad_NewtonManyChans

!-------------------------------------------------------------------------------
! (C) Crown copyright Met Office. All rights reserved.
!     Refer to COPYRIGHT.txt of this distribution for details.
!-------------------------------------------------------------------------------
! Solves the Linear equation UQ=V for Q where U is a symmetric positive definite
! matrix and U and Q are vectors of length N.  The method follows that in Golub
! and Van Loan although this is pretty standard.
!
! If U is not positive definite this will be detected by the program and flagged
! as an error.  U is assumed to be symmetric as only the upper triangle is in
! fact used.
!-------------------------------------------------------------------------------

SUBROUTINE Ops_Cholesky (U,         &
                         V,         &
                         N,         &
                         Q,         &
                         ErrorCode)

IMPLICIT NONE

! Subroutine arguments:
INTEGER, INTENT(IN)         :: n
REAL(kind_real), INTENT(IN)            :: U(n,n)
REAL(kind_real), INTENT(IN)            :: V(n)
REAL(kind_real), INTENT(OUT)           :: Q(n)
INTEGER, INTENT(OUT)        :: ErrorCode

! Local declarations:
CHARACTER(len=*), PARAMETER :: RoutineName = "Ops_Cholesky"
REAL(kind_real), PARAMETER             :: Tolerance = TINY (0.0) * 100.0
CHARACTER(len=80)           :: ErrorMessage
INTEGER                     :: j
INTEGER                     :: k
REAL(kind_real)                        :: G(n,n)   ! The Cholesky Triangle Matrix
REAL(kind_real)                        :: X(n)     ! Temporary array used in calculating G

ErrorCode = 0

! Determine the Cholesky triangle matrix.

DO j = 1, n
  X(j:n) = U(j:n,j)
  IF (j /= 1) THEN
    DO k = 1, j - 1
      X(j:n) = X(j:n) - G(j,k) * G(j:n,k)
    END DO
  END IF
  IF (X(j) <= Tolerance) THEN
    ErrorCode = 1
    Errormessage = ' :U matrix is not positive definite'
    write(*,*) RoutineName,ErrorMessage
    GOTO 9999
  END IF
  G(J:N,J) = X(J:N) / SQRT (X(J))
END DO

! Solve Gx=v for x by forward substitution

X = V
X(1) = X(1) / G(1,1)
DO j = 2, n
  X(j) = (X(j) - DOT_PRODUCT (G(j,1:j - 1), X(1:j - 1))) / G(j,j)
END DO

! Solve G^T.q=x for q by backward substitution

Q = x
Q(n) = Q(n) / G(n,n)
DO j = n - 1, 1, -1
  Q(j) = (Q(j) - DOT_PRODUCT (G(j + 1:n,j), Q(j + 1:n))) / G(j,j)
END DO

9999 CONTINUE

END SUBROUTINE Ops_Cholesky

end module ufo_onedvarfortran_minimize_newton_mod
