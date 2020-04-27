! (C) Copyright 2020 Met Office
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

!> Fortran module to provide code shared between nonlinear and tlm/adm radiance calculations

module ufo_rttovonedvarcheck_minimize_newton_mod

use iso_c_binding
use kinds
use ufo_geovals_mod
use ufo_rttovonedvarcheck_utils_mod
use ufo_radiancerttov_tlad_mod
use ufo_rttovonedvarcheck_rmatrix_mod, only: rmatrix_type
use ufo_rttovonedvarcheck_minimize_utils_mod, only: &
        ufo_rttovonedvarcheck_GeoVaLs2ProfVec, &
        ufo_rttovonedvarcheck_ProfVec2GeoVaLs, &
        ufo_rttovonedvarcheck_CostFunction, &
        ufo_rttovonedvarcheck_CheckIteration, &
        ufo_rttovonedvarcheck_CheckCloudyIteration, &
        ufo_rttovonedvarcheck_Cholesky
use ufo_rttovonedvarcheck_get_jacobian_mod, only: &
        ufo_rttovonedvarcheck_get_jacobian
use ufo_rttovonedvarcheck_profindex_mod, only: profindex_type


implicit none
private

! public subroutines
public ufo_rttovonedvarcheck_minimize_newton

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
subroutine ufo_rttovonedvarcheck_minimize_newton(self, &
                                         ob_info,       &
                                         r_matrix,      &
                                         b_matrix,      &
                                         b_inv,         &
                                         b_sigma,       &
                                         local_geovals, &
                                         profile_index, &
                                         channels,      &
                                         onedvar_success)

implicit none

type(ufo_rttovonedvarcheck), intent(inout) :: self
type(Obinfo_type), intent(in)    :: ob_info
type(rmatrix_type), intent(in)   :: r_matrix
real(kind_real), intent(in)      :: b_matrix(:,:)
real(kind_real), intent(in)      :: b_inv(:,:)
real(kind_real), intent(in)      :: b_sigma(:)
type(ufo_geovals), intent(inout) :: local_geovals
type(profindex_type), intent(in) :: profile_index
integer(c_int), intent(in)       :: channels(:)
logical, intent(out)             :: onedvar_success

! Local declarations:
character(len=*), parameter     :: RoutineName = "ufo_rttovonedvarcheck_minimize_newton"
integer                         :: inversionstatus
logical                         :: outOfRange
logical                         :: Converged
logical                         :: Error
integer                         :: iter
integer                         :: RTerrorcode
integer                         :: nchans
integer                         :: nprofelements
real(kind_real)                 :: Jcost     ! current value
real(kind_real)                 :: JcostOld  ! previous iteration value
real(kind_real)                 :: JcostOrig ! initial value
real(kind_real)                 :: DeltaJ
real(kind_real)                 :: DeltaJo

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
type(ufo_geovals)               :: geovals
real(kind_real)                 :: Jout(3)

integer                         :: ii

! ---------
! Setup
! ---------
Converged = .false.
onedvar_success = .false.
Error = .false.
nchans = size(channels)
inversionstatus = 0
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
geovals = local_geovals

call ufo_geovals_print(geovals,1)
write(*,*) "Using Newton solver"

Iterations: do iter = 1, self % max1DVarIterations

  !-------------------------
  ! 1. Generate new profile
  !-------------------------
  ! Save cost from previous iteration
  if (self % UseJForConvergence) then
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

  ! Get jacobian and hofx
  call ufo_rttovonedvarcheck_get_jacobian(geovals, ob_info, self % obsdb, &
                                       channels(:), self % conf, &
                                       profile_index, GuessProfile(:), &
                                       Y(:), H_matrix)

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

  ! Profile differences
  Ydiff(:) = ob_info%yobs(:) - Y(:)
  Diffprofile(:) = GuessProfile(:) - BackProfile(:)

  if (self % UseJForConvergence) then

    call ufo_rttovonedvarcheck_CostFunction(Diffprofile, b_inv, Ydiff, r_matrix, Jout)
    Jcost = Jout(1)

    ! exit on error
    if (inversionStatus /= 0) exit Iterations

    ! store initial cost value
    if (iter == 1) JCostOrig = jcost

    ! check for convergence
    if (iter > 1) then

      if (self % JConvergenceOption == 1) then

        ! percentage change tested between iterations
        DeltaJ = abs ((Jcost - JcostOld) / max (Jcost, tiny (0.0)))

        ! default test for checking that overall cost is getting smaller
        DeltaJo = -1.0

      else

        ! absolute change tested between iterations
        DeltaJ = abs (Jcost - JcostOld)

        ! change between current cost and initial
        DeltaJo = Jcost - JCostorig

      end if

      if (self % FullDiagnostics) then
        write (*, '(A,3F12.5)') 'Cost Function, increment, cost_convergencefactor = ', &
                                Jcost, deltaj, self % cost_convergencefactor
      end if

      if (DeltaJ < self % cost_convergencefactor .and. &
          DeltaJo < 0.0)  then ! overall is cost getting smaller?
        converged = .true.
        exit iterations
      end if

    end if

  end if ! end of specific code for cost test convergence

  ! Iterate (Guess) profile vector
  if (nchans > nprofelements) then
    write(*,*) "Many Chans"
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
    write(*,*) "Few Chans"
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
    write(*,*) "inversion failed"
    exit Iterations
  end if

  GuessProfile(:) = BackProfile(:) + Diffprofile(:)

  !---------------------------------------------------------
  ! 2. Check new profile and transfer to forward model format
  !---------------------------------------------------------

  ! Check profile and constrain humidity variables

  GuessProfileBefore(:) = GuessProfile(:)
  call ufo_rttovonedvarcheck_CheckIteration (geovals,         & ! in
                                  profile_index,   & ! in
                                  self % nlevels,  & ! in
                                  GuessProfile(:), & ! inout
                                  outOfRange)        ! out

  ! Update RT-format guess profile
  call ufo_rttovonedvarcheck_ProfVec2GeoVaLs(geovals, profile_index, nprofelements, GuessProfile)
  
  ! if qtotal in retrieval vector check cloud
  ! variables for current iteration

  if ((.NOT. outofRange) .and. profile_index % qt(1) > 0) then

    if (iter >= self % IterNumForLWPCheck) then

      if (self % RTTOV_mwscattSwitch) then
      
        !cloud information is from the scatt profile
        if (self % use_totalice) then

          call ufo_rttovonedvarcheck_CheckCloudyIteration( geovals,         & ! in
                                                profile_index,   & ! in
                                                self % nlevels,  & ! in
                                                OutOfRange )

!          call ufo_rttovonedvarcheck_CheckCloudyIteration (RTprof_Guess % rttov12_profile_scatt % clw(:),      & ! in
!                                                RTprof_Guess % rttov12_profile_scatt % totalice(:), & ! in
!                                                outOfRange)                                                         ! out
        else

          call ufo_rttovonedvarcheck_CheckCloudyIteration( geovals,         & ! in
                                                profile_index,   & ! in
                                                self % nlevels,  & ! in
                                                OutOfRange )

!          call ufo_rttovonedvarcheck_CheckCloudyIteration (RTprof_Guess % rttov12_profile_scatt % clw(:), & ! in
!                                                RTprof_Guess % rttov12_profile_scatt % ciw(:), & ! in
!                                                outOfRange)

        end if

      else

        call ufo_rttovonedvarcheck_CheckCloudyIteration( geovals,         & ! in
                                              profile_index,   & ! in
                                              self % nlevels,  & ! in
                                              OutOfRange )

      end if

    end if

  end if

  !-------------------------------------------------
  ! 3. Check for convergence using change in profile
  !    This is performed if UseJforConvergence is false
  !-------------------------------------------------

  absDiffprofile(:) = 0.0

  if ((.NOT. outOfRange) .and. (.NOT. self % UseJForConvergence))then
    absDiffProfile(:) = abs(GuessProfile(:) - OldProfile(:))
    if (ALL (absDiffProfile(:) <= B_sigma(:) * self % ConvergenceFactor)) then
      write(*,*) "Profile used for convergence"
      Converged = .true.
    end if
  end if

  !---------------------
  ! 4. output diagnostics
  !---------------------

  if (self % FullDiagnostics) then
    write (*, '(A,I0)') 'Iteration', iter
    write (*, '(A)') '------------'
    write (*, '(A,L1)') 'Status: converged = ', Converged
    if (outOfRange) write (*, '(A)') 'exiting with bad increments'
    write (*, '(A)') 'New profile:'
    call ufo_geovals_print(geovals, 1)
    write (*, '(A)')
  end if

  ! exit conditions

  if (Converged .OR. outOfRange) exit iterations

end do Iterations

! Pass convergence flag out
onedvar_success = converged

!---------------------
! 4. output diagnostics
!---------------------

if (self % UseJForConvergence) then
  write(*,'(A45,3F10.3,I5,L5)') "J initial, final, lowest, iter, converged = ", &
                                 JCostorig, Jcost,  Jcost, iter, onedvar_success
end if

!Ob % Niter = iter
!if (RTerrorcode /= 0 .OR. .NOT. Converged .OR. outOfRange .OR. &
!    inversionStatus /= 0) then
!  Error = .true.
!end if

!if (.NOT. Error .and. self % UseJForConvergence) then
!  ! store final cost function and retrieved bts
!  Ob % Jcost = Jcost
!  Ob % Britemp(Channels_1dvar(:)) = Britemp(:)
!end if

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

end subroutine ufo_rttovonedvarcheck_minimize_newton

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
real(kind_real), intent(in)     :: DeltaBT(:)        ! y-y(x)
integer, intent(in)             :: nChans
real(kind_real), intent(in)     :: H_Matrix(:,:)     ! Jacobian
real(kind_real), intent(in)     :: H_Matrix_T(:,:)   ! (Jacobian)^T
integer, intent(in)             :: nprofelements
real(kind_real), intent(inout)  :: DeltaProfile(:)   ! see note in header
real(kind_real), intent(in)     :: B_matrix(:,:)
type(rmatrix_type), intent(in)  :: r_matrix
integer, intent(out)            :: Status

! Local declarations:
character(len=*), parameter :: RoutineName = "ufo_rttovonedvarcheck_Minimize_101"
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
! 1.1. Add the R matrix into the U matrix. U = U + R
!---------------------------------------------------------------------------
!U = U + R_Matrix
call r_matrix % add_to_matrix(U, U)

! Calculate Q=(U^-1).V
!------

call ufo_rttovonedvarcheck_Cholesky (U,      &
                                     V,      &
                                     nChans, &
                                     Q,      &
                                     Status)
if (Status /= 0) goto 9999

! Delta profile is (HB)^T.Q
!------

DeltaProfile = matmul(transpose(HB), Q)

9999 continue

end subroutine ufo_rttovonedvarcheck_NewtonFewChans

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

subroutine ufo_rttovonedvarcheck_NewtonManyChans (DeltaBT,       &
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
real(kind_real), intent(in)     :: DeltaBT(:)        ! y-y(x)
integer, intent(in)             :: nChans
real(kind_real), intent(in)     :: H_Matrix(:,:)     ! Jacobian
real(kind_real), intent(in)     :: H_Matrix_T(:,:)   ! (Jacobian)^T
integer, intent(in)             :: nprofelements
real(kind_real), intent(inout)  :: DeltaProfile(:)   ! see note in header
real(kind_real), intent(in)     :: B_inverse(:,:)
type(rmatrix_type), intent(in)  :: r_matrix
integer, intent(out)            :: Status

! Local declarations:
character(len=*), parameter :: RoutineName = 'ufo_rttovonedvarcheck_NewtonManyChans'
real(kind_real)             :: HTR(nprofelements, nChans)              ! Scratch vector
real(kind_real)             :: U(nprofelements, nprofelements)         ! U = H.B.H^T + R
real(kind_real)             :: V(nprofelements)                        ! V = (y-y(x_n))-H^T(xb-x_n)

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

call ufo_rttovonedvarcheck_Cholesky (U,             &
                                     V,             &
                                     nprofelements, &
                                     DeltaProfile, &
                                     Status)

end subroutine ufo_rttovonedvarcheck_NewtonManyChans

end module ufo_rttovonedvarcheck_minimize_newton_mod
