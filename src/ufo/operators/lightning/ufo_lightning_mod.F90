! (C) Copyright 2017-2018 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

!> Fortran module for lightning observation operator

module ufo_lightning_mod

 use kinds
 use ufo_vars_mod
 use obs_variables_mod
 use oops_variables_mod
 use ufo_constants_mod, only: one, zero, half, grav     ! Gravitational field strength
 use ufo_basis_mod,     only: ufo_basis
 implicit none
 private
!> Fortran derived type for the observation type
! TODO: fill in if needed
 type, public :: ufo_lightning
   private
   type(obs_variables), public  :: obsvars
   type(oops_variables), public :: geovars
   !logical, public     ::  l_fed_nonlinear
   integer, public     ::  n_horiz
   !real(kind_real),allocatable,public  :: obsLon2d(:), obsLat2d(:)

   contains
   procedure :: simobs => ufo_lightning_simobs
 end type ufo_lightning

contains

! ------------------------------------------------------------------------------
! This code is for non-linear lightning operator.

subroutine ufo_lightning_simobs(self, geovals, obss, nvars, nlocs, hofx)
use ufo_geovals_mod, only: ufo_geovals, ufo_geoval, ufo_geovals_get_var
use iso_c_binding
use obsspace_mod

implicit none
class(ufo_lightning), intent(in)    :: self
integer, intent(in)               :: nvars, nlocs
type(ufo_geovals),  intent(in)    :: geovals
real(c_double),     intent(inout) :: hofx(nvars, nlocs)
type(c_ptr), value, intent(in)    :: obss

! Local variables
!type(ufo_geoval), pointer :: geoval
!real(kind_real), dimension(:), allocatable :: obss_metadata
  type(ufo_geoval), pointer :: delp    ! Model background values of air pressure
  type(ufo_geoval), pointer :: qg      ! Model background values of graupel mixing ratio
  logical :: ascend                         ! Flag on direction of model levels
  integer :: iobs,ivar                      ! Counter
  integer :: ilev1                          ! starting level for loop
  integer :: ilev2                          ! ending level for loop
  integer :: inc                            ! increment to level loop
  integer :: ibot                           ! index of second lowest level
  integer :: isfc                           ! index of lowest level
  integer, parameter    :: max_string = 800
  character(max_string) :: err_msg         ! Error message for output
  character(max_string) :: message         ! General message for output
  character(len=*), parameter :: myname_ = "ufo_lightning_simobs"
  !logical  ::  l_fed_nonlinear 
  integer  ::   n_horiz 

  ! check if nlocs is consistent in geovals & hofx
  if (geovals%geovals(1)%nprofiles /= size(hofx,2)*self%n_horiz) then
      write(err_msg,*) myname_, ' error: nlocs ',geovals%geovals(1)%nprofiles, &
              ' inconsistent with HofX size * n_horiz',size(hofx,2)*self%n_horiz
      call abor1_ftn(err_msg)
  endif

  ! check if some variable is in geovals and get it (var_tv is defined in ufo_vars_mod)
  call ufo_geovals_get_var(geovals, var_delp, delp)
  call ufo_geovals_get_var(geovals, var_qg, qg)


! get some metadata from obsspace
!allocate(obss_metadata(nlocs))
!call obsspace_get_db(obss, "MetaData", "some_metadata", obss_metadata)

  !l_fed_nonlinear = self%l_fed_nonlinear
  n_horiz = self%n_horiz

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
  var_loop: do ivar = 1, nvars
  obs_loop: do iobs = 1, nlocs
  !

    call Lightning_ForwardModel(delp%nval,                                     &  ! Number of pressure levels
                              qg%nval,                                         &  ! Number of mass levels
                              ilev1,                                           &  ! Starting pressure layer
                              ilev2,                                           &  ! Ending pressure layer
                              inc,                                             &  ! Layer increment
                              ibot,                                            &  ! Lowest layer inc surface
                              isfc,                                            &  ! obs-1)*Level next to surface
                              delp % vals(:,(iobs-1)*n_horiz+1:iobs*n_horiz),  &  ! Surface pressure (Pa)
                              qg %vals(:,(iobs-1)*n_horiz+1:iobs*n_horiz),     &  ! Mean layer graupel mass  (kg/kg)
                              hofx(ivar,iobs),                                 &  ! Simulated FED(flash min^-1 pixel^-1)
                              n_horiz)!
                              !l_fed_nonlinear,n_horiz)!
  end do obs_loop
  end do var_loop
  !write(err_msg,*) "TRACE: ufo_lightning_simobs: completed"

end subroutine ufo_lightning_simobs


! ------------------------------------------------------------------------------
SUBROUTINE Lightning_ForwardModel(nlevdP, &
                                nlevq,    &
                                ilev1,    &
                                ilev2,    &
                                inc,      &
                                ibot,     &
                                isfc,     &
                                delp,     &
                                qg,       &
                                Model_Lightning, &
                                !l_fed_nonlinear, &
                                n_horiz )
!-------------------------------------------------------------------------------
!> Compute Total column graupel mass amount from geovals graupel mixing ratio
!profile.
!!
!! \details Heritage: ufo_SatTCWV_mod.F90
!!
!! \author Rong Kong (CAPS/OU)
!!
!! \date 11/08/2021: Created
!!
INTEGER, INTENT(IN)            :: nlevdP                        ! Number of pressure levels
INTEGER, INTENT(IN)            :: nlevq                         ! Number of mass levels
INTEGER, INTENT(IN)            :: ilev1                         ! Starting layer
INTEGER, INTENT(IN)            :: ilev2                         ! Ending layer
INTEGER, INTENT(IN)            :: inc                           ! loop increment
INTEGER, INTENT(IN)            :: ibot                          ! lowest layer including surface
INTEGER, INTENT(IN)            :: isfc                          ! lowest level next to surface
INTEGER, INTENT(IN)            :: n_horiz
REAL(kind_real), INTENT(IN)    :: delp(1:nlevdP,n_horiz)        ! Model background pressure (Pa) at levels
REAL(kind_real), INTENT(IN)    :: qg(1:nlevq,n_horiz)           ! Model background graupel mixing ratio (kg/kg)
REAL(kind_real), INTENT(INOUT) :: Model_Lightning               ! Model forecast of the observations (flash min^-1 pixel^-1)
REAL(kind_real)                :: gridarea
REAL(kind_real), parameter     :: glmcoeff = 2.088e-8
!REAL(kind_real)                :: ag
DOUBLE PRECISION               :: ag, answer
REAL(kind_real)                :: coeff3rdorder(4)
!LOGICAL, INTENT(IN)            :: l_fed_nonlinear
!
! Local parameters
!
integer, parameter           :: max_string = 800  ! Length of strings
character(len=*), parameter  :: myname_ = "Lightning_ForwardModel"
!
! Local variables
!
character(max_string) :: err_msg           ! Error message to be output
character(max_string) :: message           ! General message for output
REAL(kind_real)       :: PDiff             ! Pressure diff across layer
REAL(kind_real)       :: GK                ! 1/gravity
INTEGER               :: ilev              ! level counter
INTEGER               :: i_horiz             

gridarea = 15000*15000/n_horiz  !the representaive area for each partner point 
!-------------------------------------------------------------------------------
! The operator is calibrated using GSL's RFFS 1-3 hour forecasts, provided hourly from 00 to 23 Z,
! specifically for the period June 11 to June 15, 2023.
!-------------------------------------------------------------------------------
coeff3rdorder(1) = 0.030104143076945   !a
coeff3rdorder(2) =  -0.441067352788872 !b
coeff3rdorder(3) =  6.250421068636704  !c
coeff3rdorder(4) =  0.0 !0.003091491729412, d
!-------------------------------------------------------------------------------
! Observation operator formula: FED = a*(ag*1.0E-9)**3 + b*(ag*1.0E-9)**2 + c*(ag*1.0E-9) + d
!-------------------------------------------------------------------------------
!1. Initialise variables and check model levels
!-------------------------------------------------------------------------------
PDiff = zero
GK = one / grav
Model_Lightning = zero
ag = zero
!
! If there is no graupel in any column, just return rapidly.
if (all(qg .eq. zero)) return
!-------------------------------------------------------------------------------
! Calculate model equivalent of total column graupel mass.
!-------------------------------------------------------------------------------
!
! Now integrate through atmosphere layers
DO ilev = ilev1, ilev2, inc
   !Pa = kg⋅m−1⋅s−2
   !the units of dp/g, Pa/(m s-2)=kg⋅m−1⋅s−2/(m s-2)=kg m-2,
   !reduce to kg m-2
   !ag = ∑(qg)dp/ g  (unit: kg m-2)
   !Horizontal integration of the graupel mass profile over a 15x15 km² area,
   !effectively reducing the extended geovals back to the observation space (i.e, from n_horiz*nlocs to nlocs).
   DO i_horiz = 1, n_horiz !Horizontal integration over n_horiz grids
     PDiff = -1.0* delp(ilev,i_horiz) 
     ag = ag - GK * qg(ilev,i_horiz) * PDiff * gridarea ! Accumulate layer graupel mass
   ENDDO
END DO
ag = ag * 1.0E-9  !The coefficients are fitted based on the relscaled data to avoid excessively small values
answer = coeff3rdorder(1)*ag**3  + coeff3rdorder(2)*ag**2+ &
                  coeff3rdorder(3)*ag     + coeff3rdorder(4)
Model_Lightning = SNGL(answer)
!print*,'Model_Lightning=',Model_Lightning
END SUBROUTINE Lightning_forwardmodel

! ------------------------------------------------------------------------------

end module ufo_lightning_mod
