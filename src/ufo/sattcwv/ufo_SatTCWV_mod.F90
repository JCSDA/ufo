!
! (C) Crown Copyright 2021 Met Office
!
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

!> Fortran module for satellite total column water vapour observation operator

module ufo_sattcwv_mod

 use iso_c_binding
 use kinds
 use oops_variables_mod
 use ufo_vars_mod
 use ufo_geovals_mod
 use ufo_geovals_mod_c, only: ufo_geovals_registry
 use ufo_basis_mod,     only: ufo_basis
 use obsspace_mod
 use fckit_log_module, only: fckit_log

! Generic routines from elsewhere in jedi
 use missing_values_mod
 use ufo_constants_mod, only: one, zero, half, grav     ! Gravitational field strength 

 implicit none
 public        :: ufo_sattcwv
 private

  !> Fortran derived type for sattcwv trajectory
 type, extends(ufo_basis) :: ufo_sattcwv
  contains
    procedure :: setup     => ufo_sattcwv_setup
    procedure :: simobs    => ufo_sattcwv_simobs
 end type ufo_sattcwv

contains

! ------------------------------------------------------------------------------
subroutine ufo_sattcwv_setup(self, f_conf)

  use fckit_configuration_module, only: fckit_configuration
  implicit none
  class(ufo_sattcwv), intent(inout)     :: self
  type(fckit_configuration), intent(in) :: f_conf

end subroutine ufo_sattcwv_setup

! ------------------------------------------------------------------------------
subroutine destructor(self)
  implicit none
  type(ufo_sattcwv), intent(inout) :: self
end subroutine destructor
! ------------------------------------------------------------------------------
! SatTCWV observation operator
subroutine ufo_sattcwv_simobs(self, geovals, hofx, obss)
  use kinds
  use ufo_geovals_mod, only: ufo_geovals, ufo_geoval, &
                             ufo_geovals_get_var, &
                             ufo_geovals_copy
  use satcolumn_mod, only: simulate_column_ob
  use iso_c_binding
  use obsspace_mod

  implicit none

  class(ufo_sattcwv), intent(in)    :: self ! The object in which this operator is 
  type(ufo_geovals),  intent(in)    :: geovals ! Model vals at Obs locs
  real(kind_real),    intent(inout) :: hofx(:)  ! Simulated Obs
  type(c_ptr), value, intent(in)    :: obss     ! The obs and metadata 

  ! Local variables
  !
  type(ufo_geoval), pointer :: prs    ! Model background values of air pressure
  type(ufo_geoval), pointer :: psfc   ! Model background values of surface pressure
  type(ufo_geoval), pointer :: q      ! Model background values of specific humidity

  integer :: iobs                           ! Counter
  integer :: nlocs                          ! number of observations
  integer, parameter    :: max_string = 800
  character(max_string) :: err_msg         ! Error message for output
  character(max_string) :: message         ! General message for output
  character(len=*), parameter :: myname_ = "ufo_sattcwv_simobs"

  ! check if nlocs is consistent in geovals & hofx
  if (geovals%nlocs /= size(hofx)) then
      write(err_msg,*) myname_, ' error: nlocs inconsistent!'
      call abor1_ftn(err_msg)
  endif

  ! get geovals of atmospheric pressure in (Pa) and q (kg/kg)
  ! 
  call ufo_geovals_get_var(geovals, var_ps, psfc)
  call ufo_geovals_get_var(geovals, var_prsi, prs)
  call ufo_geovals_get_var(geovals, var_q, q)

  ! get observation SatTCWV 

  nlocs = obsspace_get_nlocs(obss)

  write(err_msg,*) "TRACE:ufo_sattcwv_simobs: beg obs loop, nlocs=",nlocs
  call fckit_log%info(err_msg)

  obs_loop: do iobs = 1, nlocs

    call SatTCWV_ForwardModel(prs % nval, &
                              q % nval, &
                              psfc % vals(1,iobs), &
                              prs % vals(:,iobs), &
                              q % vals(:,iobs), &
                              hofx(iobs))
!
  end do obs_loop
!
  write(err_msg,*) "TRACE: ufo_sattcwv_simobs: completed"
  call fckit_log%info(err_msg)

end subroutine ufo_sattcwv_simobs

! ------------------------------------------------------------------------------
SUBROUTINE SatTCWV_ForwardModel(nlevP, &
                                nlevq, &
                                psfc, &
                                prs, &
                                q, &
                                Model_SatTCWV)
!-------------------------------------------------------------------------------
!> Compute Total column water vapour amount from geovals water vapour profile.
!!
!! \details Heritage: Var_SatTCWVOperator.f90
!!
!! \author Met Office
!!
!! \date 05/01/2021: Created
!!
INTEGER, INTENT(IN)            :: nlevP                  ! no. of p levels in state vec.
INTEGER, INTENT(IN)            :: nlevq                  ! no. of theta levels
REAL(kind_real), INTENT(IN)    :: psfc                   ! surface pressure (Pa)
REAL(kind_real), INTENT(IN)    :: prs(1:nlevP)           ! Model background pressure (Pa) at levels
REAL(kind_real), INTENT(IN)    :: q(1:nlevq)             ! Model background specific humidity (kg/kg)
REAL(kind_real), INTENT(INOUT) :: Model_SatTCWV          ! Model forecast of the observations
! 
! Local parameters
! 
integer, parameter           :: max_string = 800  ! Length of strings
character(len=*), parameter  :: myname_ = "SatTCWV_ForwardModel"
!
! Local variables
! 
character(max_string) :: err_msg           ! Error message to be output
character(max_string) :: message           ! General message for output
REAL(kind_real)       :: PDiff             ! Pressure diff across layer
REAL(kind_real)       :: Surf_Layer        ! Layer of humidity at surface
REAL(kind_real)       :: GK                ! 1/gravity
INTEGER               :: i                 ! level counter
!-------------------------------------------------------------------------------
!1. Initialise variables and check model levels
!-------------------------------------------------------------------------------
PDiff = zero
Surf_Layer = zero
GK = one / grav
Model_SatTCWV = zero

! The model data must be on a staggered grid, with nlevp = nlevq+1
IF (nlevP /= nlevQ + 1) THEN
    write(err_msg,*) myname_ // ':' // &
    ' Data must be on a staggered grid nlevp, nlevq = ',nlevp,nlevq
    call fckit_log % warning(err_msg)
    write(err_msg,*) myname_ // ':' // ' error: number of levels inconsistent!'
    call abor1_ftn(err_msg)
END IF
!-------------------------------------------------------------------------------
! 2. Calculate model equivalent of SatTCWV.
!-------------------------------------------------------------------------------
DO i = 1 , nlevP-1
   PDiff = prs(i) - prs(i+1) ! prs is pressure on level i
   Model_SatTCWV = Model_SatTCWV + GK * q(i) * PDiff ! Accumulate layer water vapour conc
END DO
! Add surface layer assuming surface q is same as q 10m but could use q2m in future
Surf_Layer = GK * q(1) * (prs(1) - psfc)
Model_SatTCWV = Model_SatTCWV + Surf_Layer

write(message,'(A,F10.4)') "Model_SatTCWV = ", Model_SatTCWV
call fckit_log%info(message)

END SUBROUTINE SatTCWV_forwardmodel

end module ufo_sattcwv_mod
