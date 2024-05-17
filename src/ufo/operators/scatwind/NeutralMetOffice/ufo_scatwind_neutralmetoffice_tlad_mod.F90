!-------------------------------------------------------------------------------
! (C) British Crown Copyright 2022 Met Office
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
!-------------------------------------------------------------------------------

!> \brief Fortran module for Met Office scatwind neutral wind forward operator
!!
!! \details TLAD of Met Office observation operator for scatterometer
!! wind data.
!!
!! \author J.Cotton (Met Office)
!!
!! \date 27/09/2022: Created
!!
module ufo_scatwind_neutralmetoffice_tlad_mod

use iso_c_binding
use kinds
use ufo_vars_mod
use ufo_geovals_mod
use ufo_geovals_mod_c, only: ufo_geovals_registry
use ufo_basis_mod,     only: ufo_basis
use obsspace_mod
use oops_variables_mod
use obs_variables_mod
use missing_values_mod
use fckit_log_module,  only : fckit_log
use fckit_exception_module,  only : fckit_exception
use vert_interp_mod
use ufo_scatwind_neutralmetoffice_mod, only: &
    ops_scatwind_phi_m_sea

implicit none

  !> Fortran derived type for neutral wind
type, public :: ufo_scatwind_neutralmetoffice_tlad
private
    type(oops_variables), public     :: geovars
    type(obs_variables), public     :: obsvars
    integer, allocatable, public     :: channels(:)
    logical                          :: surface_type_check
    integer                          :: surface_type_sea
    character(len=maxvarlen), public :: v_coord      ! GeoVaL to use to interpolate in vertical
    integer                          :: nval, nlocs
    real(kind_real), allocatable     :: CDR10(:)     ! 10m neutral winds interpolation coefficient
    real(kind_real), allocatable     :: wf(:)        ! Vertical interpolation weight
    integer, allocatable             :: wi(:)        ! Vertical interpolation index
  contains
    procedure :: setup     => ufo_scatwind_neutralmetoffice_tlad_setup
    procedure :: cleanup   => ufo_scatwind_neutralmetoffice_tlad_cleanup
    procedure :: settraj   => ufo_scatwind_neutralmetoffice_tlad_settraj
    procedure :: simobs_tl => ufo_scatwind_neutralmetoffice_simobs_tl
    procedure :: simobs_ad => ufo_scatwind_neutralmetoffice_simobs_ad
    final :: destructor
end type ufo_scatwind_neutralmetoffice_tlad

character(len=maxvarlen), dimension(2), parameter :: geovars_default = (/ &
                                                             var_u,       &
                                                             var_v /)

! ------------------------------------------------------------------------------
contains
! ------------------------------------------------------------------------------
subroutine ufo_scatwind_neutralmetoffice_tlad_setup(self,               &
                                                    channels,           &
                                                    surface_type_check, &
                                                    surface_type_sea)
  implicit none
  class(ufo_scatwind_neutralmetoffice_tlad), intent(inout) :: self
  integer(c_int), intent(in)                               :: channels(:)  !List of channels to use
  logical(c_bool), intent(in)                              :: surface_type_check
  integer(c_int), intent(in)                               :: surface_type_sea

  call self%geovars%push_back(geovars_default)
  self%v_coord = var_zimo
  self % surface_type_check = surface_type_check
  self % surface_type_sea = surface_type_sea

  ! save channels
  allocate(self%channels(size(channels)))

  if (size(channels) > 0) self%channels(:) = channels(:)

end subroutine ufo_scatwind_neutralmetoffice_tlad_setup

! ------------------------------------------------------------------------------
subroutine ufo_scatwind_neutralmetoffice_tlad_cleanup(self)

  implicit none
  class(ufo_scatwind_NeutralMetOffice_tlad), intent(inout) :: self
  self%nval = 0
  self%nlocs = 0
  if (allocated(self%channels)) deallocate(self%channels)
  if (allocated(self%wi)) deallocate(self%wi)
  if (allocated(self%wf)) deallocate(self%wf)
  if (allocated(self%CDR10)) deallocate(self%CDR10)

end subroutine ufo_scatwind_neutralmetoffice_tlad_cleanup

! ------------------------------------------------------------------------------
!> Set trajectory. Must run this routine before calling the TL or AD routines.
!!
!! \author Met Office
!!
!! \date 27/09/2022: Created
!!
! ------------------------------------------------------------------------------

subroutine ufo_scatwind_neutralmetoffice_tlad_settraj(self, geovals, obss)
  implicit none

  ! Arguments to this routine
  class(ufo_scatwind_neutralmetoffice_tlad), intent(inout) :: self
  type(ufo_geovals),                         intent(in)    :: geovals
  type(c_ptr), value,                        intent(in)    :: obss

  integer, parameter           :: max_string = 800
  character(max_string)        :: err_msg             ! Error message for output
  type(ufo_geoval), pointer    :: cx_za               ! Model heights of wind levels
  type(ufo_geoval), pointer    :: cx_friction_vel     ! Model friction velocity
  type(ufo_geoval), pointer    :: cx_obukhov_length   ! Model obukhov length
  real(kind_real), parameter   :: scatt_height = 10.0 ! height of observation in m
  real(kind_real), allocatable :: za(:)
  real(kind_real)              :: missing
  integer                      :: iobs
  integer, allocatable         :: surface_type(:)     ! Surface type qualifier

  ! Make sure nothing already allocated
  call self%cleanup()

  ! Get variables from geovals
  call ufo_geovals_get_var(geovals, self%v_coord, cx_za)                ! heights
  call ufo_geovals_get_var(geovals, var_sea_fric_vel, cx_friction_vel)  ! Friction velocity
  call ufo_geovals_get_var(geovals, var_obk_length, cx_obukhov_length)  ! Obukhov length
  self%nval = cx_za%nval
  self%nlocs = obsspace_get_nlocs(obss)

  ! check GeoVaLs are in correct vertical order (top to bottom)
  if (cx_za % vals(1,1) .lt. cx_za % vals(cx_za % nval,1) ) then
    write(err_msg,*) "GeoVaLs are not ordered from model top to bottom", cx_za % vals(1,1), cx_za % vals(cx_za % nval,1)
    call fckit_exception%throw(err_msg)
  end if

  ! Allocate arrays for interpolation index and weight
  allocate(self%wi(self%nlocs))
  allocate(self%wf(self%nlocs))
  ! Allocate 10m neutral wind interpolation coefficients and initialise to missing
  allocate(self%CDR10(self%nlocs))
  self%CDR10(:) = missing_value(self%CDR10(1))

  ! Allocate array for surface type qualifier
  allocate(surface_type(self%nlocs))
  if (self % surface_type_check) then
    call obsspace_get_db(obss, "MetaData", "surfaceQualifier", surface_type)
  end if

  allocate(za(cx_za%nval))
  do iobs = 1, self%nlocs
    ! Calculate the vertical interpolation index wi and weight wf
    za = cx_za%vals(:,iobs)
    call vert_interp_weights(cx_za%nval, scatt_height, za, self%wi(iobs), self%wf(iobs))
    ! Calculate the 10m neutral wind interpolation coefficients CDR10
    call ops_scatwind_cdr10(cx_za % nval,                     &
                            cx_za % vals(:, iobs),            &
                            cx_friction_vel % vals(1,iobs),   &
                            cx_obukhov_length % vals(1,iobs), &
                            surface_type(iobs),               &
                            self % surface_type_check,        &
                            self % surface_type_sea,          &
                            self%CDR10(iobs))
  enddo

  ! Cleanup memory
  deallocate(za)
  deallocate(surface_type)

end subroutine ufo_scatwind_neutralmetoffice_tlad_settraj

! ------------------------------------------------------------------------------
!> \brief Given an increment to the model state, calculate an increment to the
!  observation
!!
!! \author Met Office
!!
!! \date 27/09/2022: Created
!!
! ------------------------------------------------------------------------------
subroutine ufo_scatwind_neutralmetoffice_simobs_tl(self, geovals, obss, nvars, &
                                                   nlocs, hofx)
  implicit none

  ! Arguments to this routine
  class(ufo_scatwind_NeutralMetOffice_tlad), intent(in) :: self           !< The object in which this operator is contained
  integer, intent(in)                              :: nvars, nlocs        !< The number of variables and locations
  type(ufo_geovals), intent(in)                    :: geovals             !< Model perturbations
  real(c_double), intent(inout)                    :: hofx(nvars, nlocs)  !< The output increment to the observations
  type(c_ptr), value, intent(in)                   :: obss                !< The observations, and meta-data for those observations

  character(len=*), parameter        :: myname_ = "ufo_scatwind_neutralmetoffice_simobs_tl"
  integer, parameter                 :: max_string = 800

  character(max_string)              :: err_msg           ! Error message for output
  character(max_string)              :: message           ! General message for output
  integer                            :: nobs              ! Number of observations
  integer                            :: iobs              ! Loop variable, observation number
  type(ufo_geoval), pointer          :: u_d               ! Increment to eastward wind
  type(ufo_geoval), pointer          :: v_d               ! Increment to northward wind
  real(c_double)                     :: hofx_u(nlocs)
  real(c_double)                     :: hofx_v(nlocs)
  integer                            :: nchans
  integer                            :: ichan
  real(kind_real)                    :: u10               ! eastward wind at 10m
  real(kind_real)                    :: v10               ! northward wind at 10m

  write(err_msg,*) "TRACE: ufo_scatwind_neutralmetoffice_simobs_tl: begin"
  call fckit_log%info(err_msg)

  ! check if nlocs is consistent in geovals & hofx
  if (geovals%nlocs /= size(hofx(1,:))) then
    write(err_msg,*) myname_, ' error: nlocs inconsistent!'
    call abor1_ftn(err_msg)
  endif

  ! number of channels
  nchans = size(self%channels)

  ! check that hofx is the correct size for simulated variables
  ! if we have channels as second dimension then we should have 2*nchans variables
  ! if we have a single dimension then we should have 2 variables
  if (nchans /= 0) then
    if (size(hofx(:,1)) /= 2*nchans) then
      write(err_msg, '(A,I5,A,I5)') "HofX should have nchans variables for both windEastward and windNorthward. Was given ", size(hofx(:,1)), " but expected ", 2*nchans
      call fckit_exception%throw(err_msg)
    endif
  else
    if (size(hofx(:,1)) /= 2) then
      call fckit_exception%throw("HofX should have 2 variables windEastward and windNorthward")
    endif
  end if

  write(message, *) myname_, ' Running TL of Met Office neutral wind operator'
  call fckit_log%info(message)

  ! get variables from geovals
  call ufo_geovals_get_var(geovals, var_u, u_d)                        ! Eastward wind
  call ufo_geovals_get_var(geovals, var_v, v_d)                        ! Northward wind

  write(err_msg,*) "TRACE: ufo_scatwind_neutralmetoffice_simobs_tl: begin observation loop, nobs =  ", nlocs
  call fckit_log%info(err_msg)

  obs_loop: do iobs = 1, nlocs
    ! Get u,v wind components at 10m using Tl of interpolate from geovals to observational location
    call vert_interp_apply_tl(u_d % nval, u_d % vals(:, iobs), u10, self%wi(iobs), self%wf(iobs))
    call vert_interp_apply_tl(v_d % nval, v_d % vals(:, iobs), v10, self%wi(iobs), self%wf(iobs))

    ! TL forward model
    if (self%CDR10(iobs) /= missing_value(self%CDR10(iobs))) then
      hofx(1,iobs) = self%CDR10(iobs) * u10
      hofx(2,iobs) = self%CDR10(iobs) * v10
    else
      ! ob over land/seaice
      hofx(1,iobs) = missing_value(hofx(1,iobs))
      hofx(2,iobs) = missing_value(hofx(2,iobs))
    end if
  end do obs_loop

  ! if we have channels then need to spread these values across the channels correctly
  if (nchans /= 0) then
    ! windEastward hofx is stored in slot 1
    hofx_u = hofx(1,:)
    ! windNorthward hofx is stored in slot 2
    hofx_v = hofx(2,:)
    chan_loop: do ichan = 1, nchans
      hofx(ichan,:) = hofx_u
      hofx(ichan+nchans,:) = hofx_v
    end do chan_loop
  end if

  write(err_msg,*) "TRACE: ufo_scatwind_neutralmetoffice_simobs_tl: completed"
  call fckit_log%info(err_msg)

end subroutine ufo_scatwind_neutralmetoffice_simobs_tl

! ------------------------------------------------------------------------------
!> \brief Given an increment to the observation, find the equivalent increment
!         to the model state
!!
!! \author Met Office
!!
!! \date 27/09/2022: Created
!!
! ------------------------------------------------------------------------------
subroutine ufo_scatwind_neutralmetoffice_simobs_ad(self, geovals, obss, nvars, &
                                                   nlocs, hofx)
  implicit none

  ! Arguments to this routine
  class(ufo_scatwind_NeutralMetOffice_tlad), intent(in) :: self           !< The object in which this operator is contained
  integer, intent(in)                              :: nvars, nlocs        !< The number of variables and locations
  type(ufo_geovals), intent(in)                    :: geovals             !< Calculated perturbations to model state
  real(c_double), intent(inout)                    :: hofx(nvars, nlocs)  !< Increment to the observations
  type(c_ptr), value, intent(in)                   :: obss                !< Input - the observations

  character(len=*), parameter        :: myname_ = "ufo_scatwind_neutralmetoffice_simobs_ad"
  integer, parameter                 :: max_string = 800

  character(max_string)              :: err_msg           ! Error message for output
  character(max_string)              :: message           ! General message for output
  integer                            :: nobs              ! Number of observations
  integer                            :: iobs              ! Loop variable, observation number
  type(ufo_geoval), pointer          :: u_d               ! Pointer to eastward wind perturbations
  type(ufo_geoval), pointer          :: v_d               ! Pointer to northward wind perturbations
  real(c_double)                     :: hofx_u(nlocs)
  real(c_double)                     :: hofx_v(nlocs)
  integer                            :: nchans
  integer                            :: ichan
  real(kind_real)                    :: u10               ! eastward wind at 10m
  real(kind_real)                    :: v10               ! northward wind at 10m

  write(err_msg,*) "TRACE: ufo_scatwind_neutralmetoffice_simobs_ad: begin"
  call fckit_log%info(err_msg)

  ! check if nlocs is consistent in geovals & hofx
  if (geovals%nlocs /= size(hofx(1,:))) then
    write(err_msg,*) myname_, ' error: nlocs inconsistent!'
    call abor1_ftn(err_msg)
  endif

  ! number of channels
  nchans = size(self%channels)

  ! check that hofx is the correct size for simulated variables
  ! if we have channels as second dimension then we should have 2*nchans variables
  ! if we have a single dimension then we should have 2 variables
  if (nchans /= 0) then
    if (size(hofx(:,1)) /= 2*nchans) then
      write(err_msg, '(A,I5,A,I5)') "HofX should have nchans variables for both windEastward and windNorthward. Was given ", size(hofx(:,1)), " but expected ", 2*nchans
      call fckit_exception%throw(err_msg)
    endif
  else
    if (size(hofx(:,1)) /= 2) then
      call fckit_exception%throw("HofX should have 2 variables windEastward and windNorthward")
    endif
  end if

  write(message, *) myname_, ' Running adjoint of Met Office neutral wind operator'
  call fckit_log%info(message)

  ! get variables from geovals
  call ufo_geovals_get_var(geovals, var_u, u_d)                        ! Eastward wind
  call ufo_geovals_get_var(geovals, var_v, v_d)                        ! Northward wind

  write(err_msg,*) "TRACE: ufo_scatwind_neutralmetoffice_simobs_ad: begin observation loop, nobs =  ", nlocs
  call fckit_log%info(err_msg)

  obs_loop: do iobs = 1, nlocs
    if (self%CDR10(iobs) /= missing_value(self%CDR10(iobs)) .and. self%CDR10(iobs) > 0.0) then
      if (nchans /= 0) then
        ! u and v are stored according to the number of channels
        chan_loop: do ichan = 1, nchans
          u10 = hofx(ichan,iobs) * self%CDR10(iobs)
          v10 = hofx(ichan+nchans,iobs) * self%CDR10(iobs)
          ! Adjoint of interpolate, from hofx into geovals
          call vert_interp_apply_ad(u_d % nval, u_d % vals(:, iobs), u10, self%wi(iobs), self%wf(iobs))
          call vert_interp_apply_ad(v_d % nval, v_d % vals(:, iobs), v10, self%wi(iobs), self%wf(iobs))
        end do chan_loop
      else
        ! u is stored in slot 1, v is stored in slot 2
        u10 = hofx(1,iobs) * self%CDR10(iobs)
        v10 = hofx(2,iobs) * self%CDR10(iobs)
        ! Adjoint of interpolate, from hofx into geovals
        call vert_interp_apply_ad(u_d % nval, u_d % vals(:, iobs), u10, self%wi(iobs), self%wf(iobs))
        call vert_interp_apply_ad(v_d % nval, v_d % vals(:, iobs), v10, self%wi(iobs), self%wf(iobs))
      end if
    end if

  end do obs_loop

  write(err_msg,*) "TRACE: ufo_scatwind_neutralmetoffice_simobs_tl: completed"
  call fckit_log%info(err_msg)

end subroutine ufo_scatwind_neutralmetoffice_simobs_ad


! ------------------------------------------------------------------------------
!> \brief Scatterometer neutral wind coefficent CDR10
!!
!! \details The main steps are as follows:
!! * For each observation we calculate the momentum roughness length using the
!! Charnock relation (Charnock, 1955) modified to include low-wind conditions
!! (Smith, 1988):
!! \f[
!! z_{0m(sea)} = \frac{0.11\nu}{1.0 \times 10^{-5} + u_*} +
!!               \alpha_{ch} \frac{u_*^2}{g}
!! \f]
!! where
!! \f$\nu = 14 \times 10^{-6}\f$ ms-1 is the dynamic viscosity of air,
!! \f$u_{*}\f$ is the friction velocity,
!! \f$\alpha_{ch} =0.018\f$ is the charnock parameter, and
!! \f$g\f$ is the acceleration due to gravity.
!! * Call a subroutine to calculate the Monin-Obukhov stability
!! function for momentum integrated to 10m,
!! \f$\Phi_{m}\f$.
!! The calculations are dependant upon on whether we have stable or unstable
!! conditions according to the Obukhov Length, \f$L\f$.
!! * The neutral 10m wind interpolation coefficients, \f$CDR10\f$ are then
!! calculated as
!! \f[
!! CDR10=\frac{\ln\left(\left(10+z_{0m}\right)/z_{0m}\right)}
!!                       {\Phi_m\left(L,10+z_{0m},z_{0m}\right)}
!! \f]
!!
!! References:
!!
!! Charnock, H. (1955). Wind stress on a water surface. Quart. J. Royal
!!   Meteorol. Soc., 81, 639-640.
!! Smith, R. N. B. (1988). Coefficient for sea surface wind stress, heat flux
!!   and wind profiles as a function of wind speed and temperature.
!!   J. Geophys. Res., 93, 15467-15472.
!!
!! \author Met Office
!!
!! \date 27/09/2022: Created
!!
! ------------------------------------------------------------------------------
subroutine ops_scatwind_cdr10(nlevz,                   &
                              za,                      &
                              ustr,                    &
                              oblen,                   &
                              scat_surface_type,       &
                              scat_surface_type_check, &
                              scat_surface_type_sea,   &
                              cdr10)

use ufo_constants_mod, only: &
    grav                                       ! Gravitational field strength

integer, intent(in)            :: nlevz        !< no. of height levels in state vec.
real(kind_real), intent(in)    :: za(:)        !< heights of rho levs
real(kind_real), intent(in)    :: ustr         !< Model friction velocity
real(kind_real), intent(in)    :: oblen        !< Model obukhov length
integer, intent(in)            :: scat_surface_type       !< Surface type
logical, intent(in)            :: scat_surface_type_check !< Option: check for surface type being sea?
integer, intent(in)            :: scat_surface_type_sea   !< Surface type value for sea
real(kind_real), intent(inout) :: cdr10        !< 10m interpolation coefficients
!
! Local parameters
!
integer, parameter           :: max_string = 800    ! Length of strings
real(kind_real), parameter   :: scatt_height = 10.0 ! height of observation in m
real, parameter              :: charnock = 0.018    ! Charnock parameter
character(len=*), parameter  :: myname_ = "ops_scatwind_cdr"
character(max_string)        :: message             ! General message for output
!
! Local variables
!
real                         :: oblen_1      ! ObLen checked for very small numbers
real                         :: recip_l_mo   ! Reciprocal of ObLen
real                         :: z1_uv        ! Height of lowest wind (rho) level
real                         :: z0m          ! Roughness length for momentum
real                         :: phi_m_10     ! Monin-Obukhov stability function
                                             ! integrated to 10m
real                         :: phi_mn_10    ! Neutral form of stability
                                             ! function integrated to 10m
logical                      :: over_sea

if (za(nlevz) == missing_value(za(nlevz))) then  ! height missing
  write(message, *) myname_, "Missing value z1_uv"
  call abor1_ftn(message)
end if

if (oblen == missing_value(oblen)) then  ! obukhov length missing
  write(message, *) myname_, "Missing value obukhov length"
  call abor1_ftn(message)
end if

if (ustr == missing_value(ustr)) then  ! friction vel missing
  write(message, *) myname_, "Missing value friction velocity"
  call abor1_ftn(message)
end if

! Height (m) of lowest wind (rho) level
z1_uv = za(nlevz)

! Obukhov length (m) and its reciprocal
oblen_1 = sign( max(1.0E-6, abs(oblen)),oblen)
recip_l_mo = 1.0 / oblen_1

! Optionally check the surface type
over_sea = .true.
if (scat_surface_type_check) then
  if (scat_surface_type /= scat_surface_type_sea) then
    over_sea = .false.
  end if
end if

if (over_sea .and. z1_uv > 0.0) then ! over sea only

  ! Calculate roughness height for momentum
  z0m = 1.54E-6 / (1.0E-5 + ustr) + (charnock / grav) * ustr * ustr

  ! check z0m > 0 before proceeding
  if (z0m <= 0) then
    write(message, *) myname_, "Invalid roughness height"
    call abor1_ftn(message)
  end if

  ! Calculate Monin-Obukhov stability function for momentum
  ! integrated to 10m (in case this is different to level 1)
  call Ops_Scatwind_phi_m_sea (recip_l_mo,   & ! in
                               scatt_height, & ! in
                               z0m,          & ! in
                               phi_m_10)       ! out

  ! Calculate model 10m neutral wind components
  phi_mn_10 = log ((scatt_height + z0m) / z0m)

  ! Store 10m interpolation coefficients to go from 10m real
  ! wind to 10m neutral wind (consistent with VAR using fixed 10m)
  ! rather than lowest model level
  if (phi_m_10 > 0.0) then
    ! avoid zero divide
    cdr10 = (phi_mn_10 / phi_m_10)
  end if

else
  ! ob over land/seaice
  cdr10 = missing_value(cdr10)

end if

end subroutine ops_scatwind_cdr10

! ------------------------------------------------------------------------------

subroutine  destructor(self)
  type(ufo_scatwind_NeutralMetOffice_tlad), intent(inout)  :: self

  call self%cleanup()

end subroutine destructor

end module ufo_scatwind_neutralmetoffice_tlad_mod
