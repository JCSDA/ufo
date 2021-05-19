!
! (C) Crown Copyright 2021 Met Office
!
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

!> Fortran module for satellite precipitable water observation operator

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
    final :: destructor
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
! percipitable water observation operator
subroutine ufo_sattcwv_simobs(self, geovals, hofx, obss)
  use kinds
  use ufo_geovals_mod, only: ufo_geovals, ufo_geoval, ufo_geovals_get_var
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
  logical :: ascend                         ! Flag on direction of model levels
  integer :: iobs                           ! Counter
  integer :: nlocs                          ! number of observations
  integer :: ilev1                          ! starting level for loop
  integer :: ilev2                          ! ending level for loop
  integer :: inc                            ! increment to level loop 
  integer :: ibot                           ! index of second lowest level
  integer :: isfc                           ! index of lowest level
  integer, parameter    :: max_string = 800
  character(max_string) :: err_msg         ! Error message for output
  character(max_string) :: message         ! General message for output
  character(len=*), parameter :: myname_ = "ufo_sattcwv_simobs"

  ! check if nlocs is consistent in geovals & hofx
  if (geovals%nlocs /= size(hofx)) then
      write(err_msg,*) myname_, ' error: nlocs ',geovals%nlocs,' inconsistent with HofX size',size(hofx)
      call abor1_ftn(err_msg)
  endif

  ! get geovals of atmospheric pressure in (Pa) and q (kg/kg)
  ! 
  call ufo_geovals_get_var(geovals, var_ps, psfc)
  call ufo_geovals_get_var(geovals, var_prsi, prs)
  call ufo_geovals_get_var(geovals, var_q, q)
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

  ! get number of observations
  nlocs = obsspace_get_nlocs(obss)
  !
  obs_loop: do iobs = 1, nlocs
  !
    call SatTCWV_ForwardModel(prs%nval,             &  ! Number of pressure levels
                              ilev1,                &  ! Starting pressure layer
                              ilev2,                &  ! Ending pressure layer
                              inc,                  &  ! Layer increment
                              ibot,                 &  ! Lowest layer inc surface
                              isfc,                 &  ! Level next to surface
                              psfc % vals(1,iobs), &   ! Surface pressure (Pa)
                              prs % vals(:,iobs), &    ! Pressure levels (Pa)
                              q % vals(:,iobs), &      ! Mean layer specific humidity (kg/kg)
                              hofx(iobs))              ! Simulated precipitable water (kg/sq.m)
!
  end do obs_loop
!
  write(err_msg,*) "TRACE: ufo_sattcwv_simobs: completed"
  call fckit_log%info(err_msg)
!
end subroutine ufo_sattcwv_simobs
! ------------------------------------------------------------------------------
SUBROUTINE SatTCWV_ForwardModel(nlevP, &
                                ilev1, &
                                ilev2, &
                                inc,   &
                                ibot,  &
                                isfc,  &
                                psfc,  &
                                prs,   &
                                q,     &
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
INTEGER, INTENT(IN)            :: nlevp                  ! Number of pressure levels 
INTEGER, INTENT(IN)            :: ilev1                  ! Starting layer 
INTEGER, INTENT(IN)            :: ilev2                  ! Ending layer
INTEGER, INTENT(IN)            :: inc                    ! loop increment
INTEGER, INTENT(IN)            :: ibot                   ! lowest layer including surface
INTEGER, INTENT(IN)            :: isfc                   ! lowest level next to surface
REAL(kind_real), INTENT(IN)    :: psfc                   ! surface pressure (Pa)
REAL(kind_real), INTENT(IN)    :: prs(1:nlevP)           ! Model background pressure (Pa) at levels
REAL(kind_real), INTENT(IN)    :: q(1:nlevP-1)           ! Model background specific humidity (kg/kg)
REAL(kind_real), INTENT(INOUT) :: Model_SatTCWV          ! Model forecast of the observations (kg/sq.m)
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
REAL(kind_real)       :: GK                ! 1/gravity
INTEGER               :: i                 ! level counter
!-------------------------------------------------------------------------------
!1. Initialise variables and check model levels
!-------------------------------------------------------------------------------
PDiff = zero
GK = one / grav
Model_SatTCWV = zero
!
!-------------------------------------------------------------------------------
! Calculate model equivalent of satellite total column water vapour.
!-------------------------------------------------------------------------------
!
! Now integrate through atmosphere layers
DO i = ilev1, ilev2, inc
   PDiff = prs(i) - prs(i+1) ! prs is pressure on level i
! Include surface layer assuming surface q is same as q 10m but could use q2m in future
   IF (i == isfc)PDiff = psfc - prs(ibot)
!
   Model_SatTCWV = Model_SatTCWV + GK * q(i) * PDiff ! Accumulate layer water vapour conc
END DO
!
END SUBROUTINE SatTCWV_forwardmodel

end module ufo_sattcwv_mod
