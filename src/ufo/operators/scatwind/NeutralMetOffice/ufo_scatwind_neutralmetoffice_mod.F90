!-------------------------------------------------------------------------------
! (C) British Crown Copyright 2020 Met Office
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
!-------------------------------------------------------------------------------

!> \brief Fortran module for Met Office scatwind neutral wind forward operator
!!
!! \details This code introduces the Met Office observation operator for scatterometer
!! wind data. We assimilate the data as a "neutral" 10m wind, i.e. where the effects of
!! atmospheric stability are neglected. For each observation we calculate the momentum
!! roughness length using the Charnock relation. We then calculate the Monin-Obukhov
!! stability function for momentum, integrated to 10m.
!! The calculations are dependant upon on whether we have stable or unstable conditions
!! according to the Obukhov Length. The neutral 10m wind components are then calculated
!! from the 10m model winds.
!!
!! \author J.Cotton (Met Office)
!!
!! \date 22/12/2020: Created
!!
module ufo_scatwind_neutralmetoffice_mod

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

implicit none
public :: ops_scatwind_phi_m_sea
private

  !> Fortran derived type for neutral wind
type, public :: ufo_scatwind_neutralmetoffice
    type(oops_variables), public :: geovars
    type(obs_variables), public :: obsvars
    integer, allocatable, public :: channels(:)
    logical                      :: surface_type_check
    integer                      :: surface_type_sea
  contains
    procedure :: setup     => ufo_scatwind_neutralmetoffice_setup
    procedure :: delete    => ufo_scatwind_neutralmetoffice_delete
    procedure :: simobs    => ufo_scatwind_neutralmetoffice_simobs
end type ufo_scatwind_neutralmetoffice

character(len=maxvarlen), dimension(7), parameter :: geovars_default = (/ &
                                                             var_u,            &
                                                             var_v,            &
                                                             var_zimo,         &
                                                             var_sfc_ifrac,    &
                                                             var_sfc_geomz,    &
                                                             var_sea_fric_vel, &
                                                             var_obk_length /)

! ------------------------------------------------------------------------------
contains
! ------------------------------------------------------------------------------
subroutine ufo_scatwind_neutralmetoffice_setup(self,               &
                                               channels,           &
                                               surface_type_check, &
                                               surface_type_sea)
  implicit none
  class(ufo_scatwind_neutralmetoffice), intent(inout) :: self
  integer(c_int), intent(in)                          :: channels(:)  !List of channels to use
  logical(c_bool), intent(in)                         :: surface_type_check
  integer(c_int), intent(in)                          :: surface_type_sea

  call self%geovars%push_back(geovars_default)
  self % surface_type_check = surface_type_check
  self % surface_type_sea = surface_type_sea

  ! save channels
  allocate(self%channels(size(channels)))

  if (size(channels) > 0) self%channels(:) = channels(:)

end subroutine ufo_scatwind_neutralmetoffice_setup

! ------------------------------------------------------------------------------
subroutine ufo_scatwind_neutralmetoffice_delete(self)

  implicit none
  class(ufo_scatwind_NeutralMetOffice), intent(inout) :: self

  if (allocated(self%channels)) deallocate(self%channels)

end subroutine ufo_scatwind_neutralmetoffice_delete

! ------------------------------------------------------------------------------
!> Neutral wind forward operator for the Met Office system
!!
!! \author Met Office
!!
!! \date 22/12/2020: Created
!!
! ------------------------------------------------------------------------------
subroutine ufo_scatwind_neutralmetoffice_simobs(self, geovals, obss, nvars, &
                                                nlocs, hofx)
  implicit none

  ! Arguments to this routine
  class(ufo_scatwind_NeutralMetOffice), intent(in) :: self                !< The object in which this operator is contained
  integer, intent(in)                              :: nvars, nlocs        !< The number of variables and locations
  type(ufo_geovals), intent(in)                    :: geovals             !< The model values, interpolated to the obsevation locations
  real(c_double), intent(inout)                    :: hofx(nvars, nlocs)  !< The output model equivalent of the observations
  type(c_ptr), value, intent(in)                   :: obss                !< The observations, and meta-data for those observations

  character(len=*), parameter     :: myname_ = "ufo_scatwind_neutralmetoffice_simobs"
  integer, parameter              :: max_string = 800

  character(max_string)              :: err_msg           ! Error message for output
  character(max_string)              :: message           ! General message for output
  integer                            :: nobs              ! Number of observations
  integer                            :: iobs              ! Loop variable, observation number
  type(ufo_geoval), pointer          :: cx_u              ! Model column of eastward wind
  type(ufo_geoval), pointer          :: cx_v              ! Model column of northward wind
  type(ufo_geoval), pointer          :: cx_za             ! Model heights of wind levels
  type(ufo_geoval), pointer          :: cx_friction_vel   ! Model friction velocity
  type(ufo_geoval), pointer          :: cx_obukhov_length ! Model obukhov length
  real(kind_real), allocatable       :: CDR10(:)          ! 10m interpolation coefficients
  real(c_double)                     :: hofx_u(nlocs)     ! The model equivalent of windEastward
  real(c_double)                     :: hofx_v(nlocs)     ! The model equivalent of windNorthward
  integer, allocatable               :: surface_type(:)   ! Surface type qualifier
  integer                            :: nchans
  integer                            :: ichan

  write(err_msg,*) "TRACE: ufo_scatwind_neutralmetoffice_simobs: begin"
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

  write(message, *) myname_, ' Running Met Office neutral wind operator with'
  call fckit_log%info(message)

  write(message, *) 'surface_type_check =', self % surface_type_check, &
    'surface_type_sea =', self % surface_type_sea
  call fckit_log%info(message)

  ! get variables from geovals
  call ufo_geovals_get_var(geovals, var_u, cx_u)                        ! Eastward wind
  call ufo_geovals_get_var(geovals, var_v, cx_v)                        ! Northward wind
  call ufo_geovals_get_var(geovals, var_zimo, cx_za)                    ! Geopotential height of wind levels
  call ufo_geovals_get_var(geovals, var_sea_fric_vel, cx_friction_vel)  ! Friction velocity
  call ufo_geovals_get_var(geovals, var_obk_length, cx_obukhov_length)  ! Obukhov length

  ! check GeoVaLs are in correct vertical order (top to bottom)
  if (cx_za % vals(1,1) .lt. cx_za % vals(cx_za % nval,1) ) then
    call fckit_exception%throw("GeoVaLs are not ordered from model top to bottom")
  end if

  ! Allocate arrays for interpolation weights and initialise to missing data
  allocate(CDR10(nlocs))
  CDR10(:) = missing_value(CDR10(1))

  ! Allocate array for surface type qualifier
  allocate(surface_type(nlocs))
  if (self % surface_type_check) then
    call obsspace_get_db(obss, "MetaData", "surfaceQualifier", surface_type)
  end if

  write(err_msg,*) "TRACE: ufo_scatwind_neutralmetoffice_simobs: begin observation loop, nobs =  ", nlocs
  call fckit_log%info(err_msg)

  obs_loop: do iobs = 1, nlocs
    call ops_scatwind_forwardmodel(cx_za % nval,                     &
                                   cx_za % vals(:, iobs),            &
                                   cx_u % vals(:, iobs),             &
                                   cx_v % vals(:, iobs),             &
                                   cx_friction_vel % vals(1,iobs),   &
                                   cx_obukhov_length % vals(1,iobs), &
                                   surface_type(iobs),               &
                                   self % surface_type_check,        &
                                   self % surface_type_sea,          &
                                   hofx(:,iobs),                     &
                                   CDR10(iobs))
  end do obs_loop

  deallocate(surface_type)
  deallocate(CDR10)

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

  write(err_msg,*) "TRACE: ufo_scatwind_neutralmetoffice_simobs: completed"
  call fckit_log%info(err_msg)

end subroutine ufo_scatwind_neutralmetoffice_simobs

! ------------------------------------------------------------------------------
!> \brief Scatterometer forward model
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
!! * The neutral 10m wind components, \f$v_{10n}\f$ are then calculated from
!! the 10m model winds, \f$v_{10}\f$, as
!! \f[
!! \textbf{v}_{10n}=\frac{\ln\left(\left(10+z_{0m}\right)/z_{0m}\right)}
!!                       {\Phi_m\left(L,10+z_{0m},z_{0m}\right)}
!!                  \textbf{v}_{10}
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
!! \date 22/12/2020: Created
!!
! ------------------------------------------------------------------------------
subroutine ops_scatwind_forwardmodel(nlevz,                   &
                                     za,                      &
                                     u,                       &
                                     v,                       &
                                     ustr,                    &
                                     oblen,                   &
                                     scat_surface_type,       &
                                     scat_surface_type_check, &
                                     scat_surface_type_sea,   &
                                     ycalc,                   &
                                     cdr10)

use ufo_constants_mod, only: &
    grav                                       ! Gravitational field strength
use vert_interp_mod

integer, intent(in)            :: nlevz        !< no. of height levels in state vec.
real(kind_real), intent(in)    :: za(:)        !< heights of rho levs
real(kind_real), intent(in)    :: u(:)         !< Model eastward wind profile
real(kind_real), intent(in)    :: v(:)         !< Model northward wind profile
real(kind_real), intent(in)    :: ustr         !< Model friction velocity
real(kind_real), intent(in)    :: oblen        !< Model obukhov length
integer, intent(in)            :: scat_surface_type       !< Surface type
logical, intent(in)            :: scat_surface_type_check !< Option: check for surface type being sea?
integer, intent(in)            :: scat_surface_type_sea   !< Surface type value for sea
real(kind_real), intent(inout) :: ycalc(:)     !< Model equivalent of the obs
real(kind_real), intent(inout) :: cdr10        !< 10m interpolation coefficients
!
! Local parameters
!
integer, parameter           :: max_string = 800    ! Length of strings
real(kind_real), parameter   :: scatt_height = 10.0 ! height of observation in m
real, parameter              :: charnock = 0.018    ! Charnock parameter
character(len=*), parameter  :: myname_ = "Ops_Scatwind_ForwardModel"
character(max_string)        :: message             ! General message for output
!
! Local variables
!
integer                      :: wi           ! vertical interpolation index
real(kind_real)              :: wf           ! vertical interpolation weight
real(kind_real)              :: u10          ! eastward wind at 10m
real(kind_real)              :: v10          ! northward wind at 10m
real                         :: oblen_1      ! ObLen checked for very small numbers
real                         :: recip_l_mo   ! Reciprocal of ObLen
real                         :: z1_uv        ! Height of lowest wind (rho) level
real                         :: z0m          ! Roughness length for momentum
real                         :: phi_m_10     ! Monin-Obukhov stability function
                                             ! integrated to 10m
real                         :: phi_mn_10    ! Neutral form of stability
                                             ! function integrated to 10m
character(max_string)        :: err_msg      ! Error message to be output
logical                      :: over_sea

if (u(nlevz) == missing_value(u(nlevz))) then  ! u wind missing
  write(message, *) myname_, "Missing value u1"
  call abor1_ftn(message)
end if

if (v(nlevz) == missing_value(v(nlevz))) then  ! v wind missing
  write(message, *) myname_, "Missing value v1"
  call abor1_ftn(message)
end if

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

! Get u,v wind components at 10m
call vert_interp_weights(nlevz, scatt_height , za, wi, wf)
call vert_interp_apply(nlevz, u, u10, wi, wf)
call vert_interp_apply(nlevz, v, v10, wi, wf)

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
    ycalc(1) = cdr10 * u10
    ycalc(2) = cdr10 * v10
  end if

else
  ! ob over land/seaice
  ycalc(1) = missing_value(ycalc(1))
  ycalc(2) = missing_value(ycalc(2))

end if

end subroutine ops_scatwind_forwardmodel

! ------------------------------------------------------------------------------
!> \brief Calculate the integrated froms of the Monin-Obukhov stability functions
!! for surface exchanges.
!!
!! \details The main steps are as follows:
!! * In neutral conditions we have the logarithmic profile:
!! \f[
!! \Phi_{mn} = \ln\left(\frac{10 + z_{0m}}{z_{0m}}\right)
!! \f]
!! * In stable conditions the stability functions of Beljaars and Holtslag (1991)
!!   are used:
!! \f[
!! \Phi_{m} = \Phi_{mn} +
!!            a\left(\zeta_1 - \zeta_{0m}\right) +
!!            b\left(
!!                   \left(\zeta_1 - \frac{c}{d}\right)
!!                   \exp\left(-d\zeta_1\right) -
!!                   \left(\zeta_{0m} - \frac{c}{d}\right)
!!                   \exp\left(-d\zeta_{0m}\right)
!!             \right)
!! \f]
!! with
!! \f$\zeta_1=(10+z_{0m})/L\f$,
!! \f$\zeta_{0m}=z_{0m}/L\f$, and
!! a=1, b=2/3, c=5, d=0.35.
!!
!! * In unstable conditions the Dyer and Hicks forms (Dyer, 1974) are used:
!! \f[
!! \Phi_{m} = \Phi_{mn} -
!!            2\ln\left(\frac{1+X_1}{1+X_0}\right) -
!!             \ln\left(\frac{1+X_1^2}{1+X_0^2}\right) +
!!            2\left(\tan^{-1}(X_1) - \tan^{-1}(X_0)\right)
!! \f]
!! with
!! \f$X_1=(1-16\zeta_1)^{1/4}\f$
!! \f$X_0=(1-16\zeta_{0m})^{1/4}\f$
!!
!! References:
!!
!! Beljaars, A. C. M. and Holtslag, A. A. M. (1991). Flux parameterisation
!!   over land surfaces for atmospheric models. J. Appl. Meteor., 30,
!!   327-341.
!!
!! Dyer, A. J., 1974: A review of flux-profile relationships. Bound. Layer
!!   Meteor., 7, 363-372.
!!
!! \author Met Office
!!
!! \date 22/12/2020: Created
!!
! ------------------------------------------------------------------------------

subroutine ops_scatwind_phi_m_sea (recip_l_mo, &
                                   z_uv,       &
                                   z0m,        &
                                   phi_m)

implicit none
! Subroutine arguments:
real, intent(in)            :: recip_l_mo !< Reciprocal of Monin-Obukhov length (m^-1).
real(kind_real), intent(in) :: z_uv       !< Height of wind level above roughness height(m).
real, intent(in)            :: z0m        !< Roughness length for momentum (m).
real, intent(out)           :: phi_m      !< Stability function for momentum.
! Local declarations:
character(len=*), parameter :: RoutineName = 'ops_scatwind_phi_m_sea'
real, parameter             :: a = 1.0
real, parameter             :: b = 2.0 / 3.0
real, parameter             :: c = 5.0
real, parameter             :: d = 0.35
real, parameter             :: c_over_d = c / d
real                        :: phi_mn      ! Neutral value of stability function for momentum
real                        :: zeta_uv
real                        :: zeta_0m
real                        :: x_uv_sq
real                        :: x_0m_sq
real                        :: x_uv
real                        :: x_0m

! Calculate neutral value of PHI_M
phi_mn = log ((z_uv + z0m)/z0m)

! Calculate stability parameters
zeta_uv = (z_uv + z0m) * recip_l_mo
zeta_0m = z0m * recip_l_mo

if (recip_l_mo  >=  0.0) then
  ! Calculate PHI_M for neutral and stable conditions.
  ! Formulation of Beljaars and Holtslag (1991).
  phi_m = phi_mn  +                                         &
          a * (zeta_uv - zeta_0m) +                         &
          b * ((zeta_uv - c_over_d) * exp (-d * zeta_uv) -  &
               (zeta_0m - c_over_d) * exp (-d * zeta_0m))

else
  ! Calculate PHI_M for unstable conditions.
  x_uv_sq = sqrt (1.0 - 16.0 * zeta_uv)
  x_0m_sq = sqrt (1.0 - 16.0 * zeta_0m)
  x_uv = sqrt (x_uv_sq)
  x_0m = sqrt (x_0m_sq)

  phi_m = phi_mn - 2.0 * log ((1.0 + x_uv) / (1.0 + x_0m)) - &
                   log ((1.0 + x_uv_sq) / (1.0 + x_0m_sq)) + &
                   2.0 * (atan (x_uv) - atan (x_0m))

end if

end subroutine ops_scatwind_phi_m_sea

end module ufo_scatwind_neutralmetoffice_mod
