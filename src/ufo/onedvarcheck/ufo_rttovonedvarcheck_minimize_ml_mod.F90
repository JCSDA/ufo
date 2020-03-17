! (C) Copyright 2018 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

!> Fortran module to provide code shared between nonlinear and tlm/adm radiance calculations

module ufo_rttovonedvarcheck_minimize_ml_mod

use iso_c_binding
use config_mod
use kinds
use ufo_geovals_mod
use ufo_rttovonedvarcheck_utils_mod
use ufo_radiancerttov_tlad_mod

implicit none

private

! subroutines - all listed for complete
public ufo_rttovonedvarcheck_MinimizeML

contains

!------------------------------------------------------------------------------

subroutine ufo_rttovonedvarcheck_MinimizeML(ob_info, r_matrix, Syinv, b_matrix, Sxinv, &
                                     local_geovals, profile_index, nprofelements, &
                                     conf, obsdb, channels, onedvar_success)

use ufo_rttovonedvarcheck_process_mod, only: &
                    ufo_rttovonedvarcheck_GeoVaLs2ProfVec, &
                    ufo_rttovonedvarcheck_ProfVec2GeoVaLs, &
                    ufo_rttovonedvarcheck_CostFunction, &
                    Ops_SatRad_Qsplit

use ufo_rttovonedvarcheck_forward_model_mod, only: &
                    ufo_rttovonedvarcheck_ForwardModel

implicit none

type(Obinfo_type), intent(in)               :: ob_info
real(kind_real), intent(in)                 :: r_matrix(:,:)
real(kind_real), intent(in)                 :: Syinv(:,:)
real(kind_real), intent(in)                 :: b_matrix(:,:)
real(kind_real), intent(in)                 :: Sxinv(:,:)
type(ufo_geovals), intent(inout)            :: local_geovals
type(Profileinfo_type), intent(in)          :: profile_index
integer, intent(in)                         :: nprofelements
type(c_ptr), VALUE, intent(in)              :: conf
type(c_ptr), VALUE, intent(in)              :: obsdb
integer(c_int), intent(in)                  :: channels(:)
logical, intent(out)                        :: onedvar_success

! Local Variables
real(kind_real), allocatable       :: X(:)    ! x vector for 1D-var
real(kind_real), allocatable       :: Xb(:)   ! Xb background vector for 1D-var
real(kind_real), allocatable       :: X0(:)   ! Xb background vector for 1D-var
real(kind_real)                    :: J_old        ! initial cost before
real(kind_real)                    :: J_new        ! initial cost before
real(kind_real), allocatable       :: H_matrix(:,:)
real(kind_real), allocatable       :: Kx(:,:)
integer                            :: nchans
real(kind_real), allocatable       :: Xdiff(:)      ! profile increments
real(kind_real), allocatable       :: Ydiff(:)       ! BT differences
real(kind_real), allocatable       :: Y(:)       ! update hofx
real(kind_real), allocatable       :: Y0(:)       ! initial hofx
integer                            :: iter              ! iterative loop counter
real(kind_real), allocatable       :: unit_matrix(:,:)
integer                            :: i
real(kind_real)                    :: Jzero
real(kind_real)                    :: Jcurrent
real(kind_real)                    :: Jinitial
real(kind_real)                    :: delta_J         ! Cost difference
real(kind_real), allocatable       :: Sx(:,:)         ! Error covariance matrix of the a priori (background) state vector
real(kind_real), allocatable       :: invVarRad(:)    ! inverse of variances of radiances
real(kind_real), allocatable       :: Xplus_dX(:)     ! State vector being checked in a given iteration (Xn + delta_X)
real(kind_real), allocatable       :: dJ_dX(:)        ! 1st derivative of J wrt state variables (Tau and Re)
real(kind_real), allocatable       :: d2J_dX2(:,:)    ! 2nd derivative of J wrt state variables
real(kind_real), allocatable       :: delta_X(:)      ! Step to be applied to state variables
real(kind_real), allocatable       :: J2plus_A(:,:)   ! Temporary array to hold sum of d2J_dX2 and alpha * unit_matrix.
real(kind_real), allocatable       :: KxT_Syinv(:,:)  ! intermediate matrix product

!real(kind_real), parameter         :: Ccj = 0.01      ! Cost convergence criterion
real(kind_real)                    :: Ccj_Ny          ! Cost convergence criterion adjusted for number of channels
logICAL                            :: convergence
integer                            :: CholeskyStatus
real(kind_real)                    :: Av_Hess         ! Average of Hessian (d2J_dX2) diagonal
real(kind_real)                    :: alpha           ! Marquardt control variable
integer                            :: m, ii
real(kind_real)                    :: sum_errors
real(kind_real)                    :: DeltaJ
type(ufo_geovals)                  :: geovals
real(kind_real)                    :: Jlowest
real(kind_real)                    :: Jout(3)

! Values to move to yaml file
real(kind_real), parameter         :: Ccj = 0.01      ! Cost convergence criterion
real(kind_real), parameter         :: Mqstart = 0.001 ! Marquardt starting parameter
real(kind_real), parameter         :: Mqstep = 5.0    ! Marquardt step parameter
real(kind_real), parameter         :: MaxIter = 10

onedvar_success = .false. ! Assume failed

write(*,*) "ufo_rttovonedvarcheck_MinimizeML start"

nchans = size(channels)
geovals = local_geovals

write(*,*) "ufo_rttovonedvarcheck_MinimizeML nchans = ",nchans

! allocate arrays
allocate(Xb(nprofelements))
allocate(X(nprofelements))
allocate(X0(nprofelements))
allocate(H_matrix(nchans,nprofelements))
allocate(Kx(nchans,nprofelements))
allocate(Xdiff(nprofelements))
allocate(Ydiff(nchans))
allocate(Y(nchans))
allocate(Y0(nchans))
allocate(unit_matrix(nprofelements,nprofelements))
allocate(Sx(nprofelements,nprofelements))
allocate(invVarRad(nchans))
allocate(Xplus_dX(nprofelements))
allocate(dJ_dX(nprofelements))
allocate(d2J_dX2(nprofelements,nprofelements))
allocate(delta_X(nprofelements))
allocate(J2plus_A(nprofelements,nprofelements))
allocate(KxT_Syinv(nprofelements,nchans))

! Set everything to zero
Xb(:) = 0.0
X(:) = 0.0
X0(:) = 0.0
H_matrix(:,:) = 0.0
Kx(:,:) = 0.0
Xdiff(:) = 0.0
Ydiff(:) = 0.0
Y(:) = 0.0
Y0(:) = 0.0
unit_matrix(:,:) = 0.0
Sx(:,:) = 0.0
invVarRad(:) = 0.0
Xplus_dX(:) = 0.0
dJ_dX(:) = 0.0
d2J_dX2(:,:) = 0.0
delta_X(:) = 0.0
J2plus_A(:,:) = 0.0
KxT_Syinv(:,:) = 0.0

! Map GeovaLs to 1D-var profile using B matrix profile structure
call ufo_rttovonedvarcheck_GeoVaLs2ProfVec(geovals, profile_index, nprofelements, X0(:))
Xb(:) = X0(:)

! call forward model
call ufo_rttovonedvarcheck_ForwardModel(geovals, ob_info, obsdb, &
                                   channels(:), conf, &
                                   profile_index, X(:), &
                                   Y0(:), H_matrix)

Y(:) = Y0(:)

Kx(:, :) = H_matrix(:, :)
write(*,*) "Kx(:,:) = ",Kx

!  Calculate Ydiff, the difference between the measurements and calculated values for X.
Ydiff(:) = ob_info%yobs(:) - Y(:) ! amended from previous.

! Update state vector for first iteration

X(:) = X0(:)

!  Calculate Xdiff, the difference between the current state and the a priori (background)

Xdiff(:) = X(:) - Xb(:) ! Will be zero initially in this case, but not generally true

!  Calculate cost at X0.
write(*,*) "Initial X(:) = ",X(:)
write(*,*) "Initial H(X)(:) = ",Y(:)
write(*,*) "Obs = ",ob_info%yobs(:)

call ufo_rttovonedvarcheck_CostFunction(Xdiff, Sxinv, Ydiff, Syinv, Jout)
Jzero = Jout(1)
Jinitial = Jzero

!  Adjust cost convergence criteria for the number of active channels,
!  since final costs will be divided by this for output.

Ccj_Ny = Ccj * nchans

!  Set up unit_matrix matrix

unit_matrix(:,:) = 0.0
do i = 1, nprofelements
  unit_matrix(i,i) = 1.0
end do

convergence = .false.
iter = 1
CholeskyStatus = 0
delta_X(:) = 0.0

! Calculate derivatives of starting cost (required for setting Marquardt
! parameters). Needs Kx, the scaled version of dY_dX.

KxT_Syinv = matmul(transpose(Kx), Syinv) ! intermediate product, used twice below.
dJ_dX   = 2*(matmul(Sxinv, Xdiff)) - 2*(matmul(KxT_Syinv, Ydiff)) ! amended from previous.
d2J_dX2 = 2*(matmul(KxT_Syinv, Kx) + Sxinv) ! amended from previous.

! Set Marquardt parameters (initial weighting favours steepest descent)
! unit_matrix matrix and inverse function required.

! Average leading diagonal of Hessian d2J_dX2 in order to set alpha.,
! then use alpha and d2J_dX2 to set delta_X.
! Solve_Cholesky is used to find delta_X. The equation to be solved is:
! dX = -(J'' + alphaI)^-1 * J' or delta_X = - (inv(J2plus_A) * dJ_dX)
! Multiplying through by J2_plusA we get:
!   J2plus_A * delta_X = -dJ_dX
! which we can solve for delta_X using Solve_Cholesky.

Av_Hess = 0
do m = 1, nprofelements
  Av_Hess = Av_Hess + d2J_dX2(m,m)
end do
Av_Hess  = Av_Hess / FLOAT(nprofelements)
alpha    = MqStart * Av_Hess
J2plus_A = d2J_dX2 + (alpha * unit_matrix)

write(*,*) "d2J_dX2 = ",d2J_dX2
write(*,*) "J2plus_A = ",J2plus_A
write(*,*) "dJ_dX = ",dJ_dX

call ufo_rttovonedvarcheck_SolveCholesky(J2plus_A, -dJ_dX, delta_X, nprofelements, CholeskyStatus)
if (CholeskyStatus /= 0) write(*,*) 'Minimize: Error in Solve_Cholesky'

! Main iteration loop
!Jzero = 1.0e4
Jlowest = Jzero

Main_Loop_index: do

  ! N.B. Jumps out of the loop if convergence occurs,
  ! or if Max no. of iterations is reached, or if CholeskyStatus not zero

  if (CholeskyStatus /= 0) exit Main_Loop_index

  ! Apply step delta_x to the active state variables
  ! Assumes Xn and delta_X are both unscaled.

  write(*,*) "X in loop = ",X(:)
  write(*,*) "delta_X in loop = ",delta_X(:)
  Xplus_dX(:)  = X(:) + delta_X(:)

! ------------------------------------------------------------------------------------------------------------
! do the next forward model calculations. i.e. Calculate Y for Xn + delta_X. Xplus_dX is currently un-scaled.
! ------------------------------------------------------------------------------------------------------------

  ! Profile to GeoVaLs
  call ufo_rttovonedvarcheck_ProfVec2GeoVaLs(geovals, profile_index, nprofelements, Xplus_dX)

  ! call forward model
  call ufo_rttovonedvarcheck_ForwardModel(geovals, ob_info, obsdb, &
                                       channels(:), conf, &
                                       profile_index, Xplus_dx(:), &
                                       Y(:), H_matrix)

  Kx(:, :) = H_matrix(:, :)

  ! Calculate new cost, J.

  Xdiff(:) = Xplus_dX(:) - Xb(:)
  Ydiff(:) = ob_info%yobs(:) - Y(:)  ! amended from previous.

  call ufo_rttovonedvarcheck_CostFunction(Xdiff, Sxinv, Ydiff, Syinv, Jout)
  Jcurrent = Jout(1)

  !        Check J vs previous value. Lower value means progress:
  !        "accept" the step delta_X and reset the Marquardt values.

  delta_J = Jcurrent-Jzero

  if (Jcurrent < Jlowest) Jlowest = Jcurrent

  if (Jcurrent <= Jzero) then
    write(*,*) "Good step"
    ! Calculate the Hessian (J'') as this is required for error analysis
    ! if convergence is reached or for setting new delta_X.

    X(:) = Xplus_dX(:)
    KxT_Syinv = matmul(transpose(Kx), Syinv)
    d2J_dX2 = 2*(matmul(KxT_Syinv, Kx) + Sxinv) ! amended from previous.

    ! Check for convergence.
    DeltaJ = ABS (Jcurrent - Jzero)
    !sum_errors = 0.0
    !do ii=1,nprofelements
    !  sum_errors = sum_errors + (dJ_dX(ii)*dJ_dX(ii))
    !end do
    !if (SQRT(sum_errors) <= Ccj_Ny) then ! Note: dJ_dX is scaled
    !  convergence = .true.
    !if (DeltaJ < Ccj_Ny .and. (Jcurrent - Jzero) < 0.0_kind_real) then
    !  convergence = .true.
    !else
    if (DeltaJ < Ccj .and. (Jcurrent - Jzero) < 0.0_kind_real) then
      convergence = .true.
    else
    ! Decrease steepest descent part for next iteration
    ! Save new J as J0 for checking next time round

      alpha = alpha / MqStep
      Jzero = Jcurrent

      ! Calculate values required to set the new step delta_X based on
      ! derivatives of J (delta_X itself is set after this "if")

      dJ_dX   = 2*(matmul(Sxinv, Xdiff)) - 2*(matmul(KxT_Syinv, Ydiff)) ! amended from previous.
      J2plus_A = d2J_dX2 + (alpha * unit_matrix)
      call ufo_rttovonedvarcheck_SolveCholesky(J2plus_A, -dJ_dX, delta_X, nprofelements, CholeskyStatus)
      if (CholeskyStatus /= 0) print *, ' Error in Solve_Cholesky'

    end if

  else ! No improvement in cost

    write(*,*) "Bad step"

  ! increase steepest descent part for next iteration and set values
  ! required for setting new deltaX using old cost derivatives.

    alpha = alpha * MqStep
    J2plus_A = d2J_dX2 + (alpha * unit_matrix)
    call ufo_rttovonedvarcheck_SolveCholesky(J2plus_A, -dJ_dX, delta_X, nprofelements, CholeskyStatus)
    !           Check CholeskyStatus returned by ufo_rttovonedvarcheck_SolveCholesky later
  end if

  if (convergence) then
    onedvar_success = .true.
    exit Main_Loop_index ! Drops out of the iteration do loop
  end if

  ! increment iteration counter and check whether the maximum has been
  ! passed. (increment after breakpoint output so that counter value
  ! is correct in the output).

  if (iter == MaxIter) then
    onedvar_success = .false.
    exit Main_Loop_index ! Drops out of the iteration do loop
  end if
  iter = iter + 1

  write(*,'(A40,4F10.3,I5)') "J initial, J current, Jb, Jo, iter = ",Jinitial,Jcurrent,Jout(2),Jout(3),iter

end do Main_Loop_index

write(*,'(A34,3F10.3,I5)') "J initial, final, lowest, iter = ",Jinitial,Jcurrent,Jlowest,iter
write(*,*) "Obs, initial hofx, final hofx = ",ob_info%yobs(1),Y0(1),Y(1)

! deallocate arrays
call ufo_geovals_delete(geovals)
if (allocated(X))           deallocate(X)
if (allocated(X0))          deallocate(X0)
if (allocated(Xb))          deallocate(Xb)
if (allocated(H_matrix))    deallocate(H_matrix)
if (allocated(Kx))          deallocate(Kx)
if (allocated(Xdiff))       deallocate(Xdiff)
if (allocated(Ydiff))       deallocate(Ydiff)
if (allocated(Y))           deallocate(Y)
if (allocated(Y0))          deallocate(Y0)
if (allocated(unit_matrix)) deallocate(unit_matrix)
if (allocated(Sx))          deallocate(Sx)
if (allocated(invVarRad))   deallocate(invVarRad)
if (allocated(Xplus_dX))    deallocate(Xplus_dX)
if (allocated(dJ_dX))       deallocate(dJ_dX)
if (allocated(d2J_dX2))     deallocate(d2J_dX2)
if (allocated(delta_X))     deallocate(delta_X)
if (allocated(J2plus_A))    deallocate(J2plus_A)
if (allocated(KxT_Syinv))   deallocate(KxT_Syinv)

write(*,*) "ufo_rttovonedvarcheck_MinimizeML end"

end subroutine ufo_rttovonedvarcheck_MinimizeML

!---------------------------------------------------------------------------

subroutine ufo_rttovonedvarcheck_SolveCholeskyDC(A, b, x, n)

   implicit none

   integer, intent(in) :: n
   real(kind_real), intent(in) :: A(n,n)
   real(kind_real), intent(in) :: b(n)
   real(kind_real), intent(out) :: x(n)
   integer :: i

  write(*,*) "ufo_rttovonedvarcheck_SolveCholeskyDC start"
!
! Solve L.y = b, storing y in x
!
   do i = 1, n
      x(i) = (b(i) - dot_product(A(i,1:i-1), x(1:i-1))) / A(i,i)
   end do
!
! Solve L^T·x = y
!
   do i = n, 1, -1
      x(i) = (x(i) - dot_product(A(i+1:n,i), x(i+1:n))) / A(i,i)
   end do

  write(*,*) "ufo_rttovonedvarcheck_SolveCholeskyDC end"

end subroutine ufo_rttovonedvarcheck_SolveCholeskyDC

!---------------------------------------------------------------------------

subroutine ufo_rttovonedvarcheck_DecomposeCholesky(A, n, Status)

   implicit none

   integer, intent(in) :: n
   real(kind_real), intent(inout) :: A(n,n)
   integer, intent(out) :: Status
   integer :: i

   write(*,*) "ufo_rttovonedvarcheck_DecomposeCholesky start"

   Status = 0

   do i = 1, n
      A(i,i) = A(i,i) - dot_product(A(i,1:i-1), A(i,1:i-1))
      if (A(i,i) <= 0.0_kind_real) then
         Status = 9999
         print *,"Error in decomposing cholesky"
         return
      end if
      A(i,i) = sqrt(A(i,i))
      A(i+1:n,i) = (A(i,i+1:n) - matmul(A(i+1:n,1:i-1), A(i,1:i-1))) / A(i,i)
   end do

   write(*,*) "ufo_rttovonedvarcheck_DecomposeCholesky end"

end subroutine ufo_rttovonedvarcheck_DecomposeCholesky

!---------------------------------------------------------------------------

subroutine ufo_rttovonedvarcheck_SolveCholesky(A, b, x, n, Status)

   implicit none

   integer, intent(in) :: n
   real(kind_real), intent(in) :: A(n,n)
   real(kind_real), intent(in) :: b(n)
   real(kind_real), intent(out) :: x(n)
   integer, intent(out) :: Status
   
   real(kind_real) :: D(n,n)

   write(*,*) "ufo_rttovonedvarcheck_SolveCholesky start"

   D = A
    
   call ufo_rttovonedvarcheck_DecomposeCholesky(D, n, Status)
   if (Status /= 0) return
   call ufo_rttovonedvarcheck_SolveCholeskyDC(D, b, x, n)

   write(*,*) "ufo_rttovonedvarcheck_SolveCholesky end"

end subroutine ufo_rttovonedvarcheck_SolveCholesky

end module ufo_rttovonedvarcheck_minimize_ml_mod
