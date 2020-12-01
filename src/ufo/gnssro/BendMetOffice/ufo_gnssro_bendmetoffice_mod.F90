!-------------------------------------------------------------------------------
! (C) British Crown Copyright 2020 Met Office
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
!-------------------------------------------------------------------------------

!> Fortran module for gnssro bending angle Met Office forward operator

module ufo_gnssro_bendmetoffice_mod

use iso_c_binding
use kinds
use ufo_vars_mod
use ufo_geovals_mod
use ufo_geovals_mod_c, only: ufo_geovals_registry
use ufo_basis_mod,     only: ufo_basis
use vert_interp_mod
use lag_interp_mod,    only: lag_interp_const, lag_interp_smthWeights
use obsspace_mod  
use missing_values_mod
use ufo_gnssro_ukmo1d_utils_mod
use fckit_log_module,  only : fckit_log

implicit none
public             :: ufo_gnssro_bendmetoffice
private

  !> Fortran derived type for gnssro trajectory
type, extends(ufo_basis) :: ufo_gnssro_BendMetOffice
  logical :: vert_interp_ops
  logical :: pseudo_ops
  contains
    procedure :: setup     => ufo_gnssro_bendmetoffice_setup
    procedure :: simobs    => ufo_gnssro_bendmetoffice_simobs
end type ufo_gnssro_BendMetOffice

contains

! ------------------------------------------------------------------------------
! Get the optional settings for the forward model, and save them in the object
! so that they can be used in the code.
! ------------------------------------------------------------------------------
subroutine ufo_gnssro_bendmetoffice_setup(self, f_conf)

use fckit_configuration_module, only: fckit_configuration
implicit none
class(ufo_gnssro_BendMetOffice), intent(inout) :: self
type(fckit_configuration), intent(in)   :: f_conf

call f_conf%get_or_die("vert_interp_ops", self % vert_interp_ops)
call f_conf%get_or_die("pseudo_ops", self % pseudo_ops)

end subroutine ufo_gnssro_bendmetoffice_setup

! ------------------------------------------------------------------------------
! 1-dimensional GNSS-RO forward operator for the Met Office system
! ------------------------------------------------------------------------------
subroutine ufo_gnssro_bendmetoffice_simobs(self, geovals, hofx, obss)

  implicit none

  ! Arguments to this routine
  class(ufo_gnssro_BendMetOffice), intent(in)    :: self     ! The object in which this operator is contained
  type(ufo_geovals),               intent(in)    :: geovals  ! The model values, interpolated to the obsevation locations
  real(kind_real),                 intent(inout) :: hofx(:)  ! The model forecast of the observations
  type(c_ptr), value,              intent(in)    :: obss     ! The observations, and meta-data for those observations

  character(len=*), parameter     :: myname_ = "ufo_gnssro_bendmetoffice_simobs"
  integer, parameter              :: max_string = 800

  character(max_string)              :: err_msg         ! Error message for output
  character(max_string)              :: message         ! General message for output
  integer                            :: nobs            ! Number of observations
  integer                            :: ilev            ! Loop variable, level number
  integer                            :: iobs            ! Loop variable, observation number
  type(ufo_geoval), pointer          :: q               ! Model background values of specific humidity
  type(ufo_geoval), pointer          :: prs             ! Model background values of air pressure
  type(ufo_geoval), pointer          :: theta_heights   ! Model heights of levels containing specific humidity
  type(ufo_geoval), pointer          :: rho_heights     ! Model heights of levels containing air pressure
  real(kind_real), allocatable       :: obsLat(:)             ! Latitude of the observation
  real(kind_real), allocatable       :: obsLon(:)             ! Longitude of the observation
  real(kind_real), allocatable       :: impact_param(:)       ! Impact parameter of the observation
  real(kind_real), allocatable       :: radius_curv(:)        ! Earth's radius of curvature at the observation tangent point
  real(kind_real), allocatable       :: undulation(:)         ! Undulation - height of the geoid above the ellipsoid
  logical                            :: flip_data             ! Whether to reverse the order of the model data
  logical                            :: BAErr                 ! Was there an error in the calculation?

  write(err_msg,*) "TRACE: ufo_gnssro_bendmetoffice_simobs: begin"
  call fckit_log%info(err_msg)

! check if nlocs is consistent in geovals & hofx
  if (geovals%nlocs /= size(hofx)) then
      write(err_msg,*) myname_, ' error: nlocs inconsistent!'
      call abor1_ftn(err_msg)
  endif

  write(message, *) myname_, ' Running Met Office GNSS-RO forward operator with'
  call fckit_log%info(message)
  write(message, *) 'vert_interp_ops =', self % vert_interp_ops, &
    'pseudo_ops =', self % pseudo_ops
  call fckit_log%info(message)

! get variables from geovals
  call ufo_geovals_get_var(geovals, var_q, q)               ! specific humidity
  call ufo_geovals_get_var(geovals, var_prsi, prs)          ! pressure
  call ufo_geovals_get_var(geovals, var_z, theta_heights)   ! Geopotential height of the normal model levels
  call ufo_geovals_get_var(geovals, var_zi, rho_heights)    ! Geopotential height of the pressure levels

  write(message, '(A,10I6)') 'Q: ', q%nval, q%nlocs, shape(q%vals)
  call fckit_log%info(message)
  write(message, '(A,10I6)') 'Pressure: ', prs%nval, prs%nlocs, shape(prs%vals)
  call fckit_log%info(message)

  nobs  = obsspace_get_nlocs(obss)

  flip_data = .false.
  if (prs%vals(1,1) .lt. prs%vals(prs%nval,1) ) then
    write(err_msg,'(a)') '  ufo_gnssro_bendmetoffice_simobs:'//new_line('a')//                         &
                         '  Model vertical height profile is in descending order,'//new_line('a')// &
                         '  The data will be flipped for processing'
    call fckit_log%info(err_msg)
    flip_data = .true.
  end if

! set obs space struture
  allocate(obsLon(nobs))
  allocate(obsLat(nobs))
  allocate(impact_param(nobs))
  allocate(radius_curv(nobs))
  allocate(undulation(nobs))

  call obsspace_get_db(obss, "MetaData", "longitude", obsLon)
  call obsspace_get_db(obss, "MetaData", "latitude", obsLat)
  call obsspace_get_db(obss, "MetaData", "impact_parameter", impact_param)
  call obsspace_get_db(obss, "MetaData", "earth_radius_of_curvature", radius_curv)
  call obsspace_get_db(obss, "MetaData", "geoid_height_above_reference_ellipsoid", undulation)

  call fckit_log%info(err_msg)

  obs_loop: do iobs = 1, nobs 

    call fckit_log%info(err_msg)

    if (flip_data) then
        call Ops_GPSRO_ForwardModel(prs % nval, &
                                    q % nval, &
                                    rho_heights % vals(rho_heights%nval:1:-1, iobs), &
                                    theta_heights % vals(theta_heights%nval:1:-1, iobs), &
                                    prs % vals(prs%nval:1:-1, iobs), &
                                    q % vals(q%nval:1:-1, iobs), &
                                    self % pseudo_ops, &
                                    self % vert_interp_ops, &
                                    1, &
                                    impact_param(iobs:iobs), &
                                    radius_curv(iobs), &
                                    obsLat(iobs), &
                                    undulation(iobs), &
                                    hofx(iobs:iobs), &
                                    BAErr)
    else
        call Ops_GPSRO_ForwardModel(prs % nval, &
                                    q % nval, &
                                    rho_heights % vals(:,iobs), &
                                    theta_heights % vals(:,iobs), &
                                    prs % vals(:,iobs), &
                                    q % vals(:,iobs), &
                                    self % pseudo_ops, &
                                    self % vert_interp_ops, &
                                    1, &
                                    impact_param(iobs:iobs), &
                                    radius_curv(iobs), &
                                    obsLat(iobs), &
                                    undulation(iobs), &
                                    hofx(iobs:iobs), &
                                    BAErr)
    end if

    if (BAErr) then
      write(err_msg,*) "Error with observation processing ", iobs
      call fckit_log % info(err_msg)
    end if

  end do obs_loop

  deallocate(obsLat)
  deallocate(obsLon)
  deallocate(impact_param)
  deallocate(radius_curv)
  deallocate(undulation)

  write(err_msg,*) "TRACE: ufo_gnssro_bendmetoffice_simobs: completed"
  call fckit_log%info(err_msg)

end subroutine ufo_gnssro_bendmetoffice_simobs
! ------------------------------------------------------------------------------


SUBROUTINE Ops_GPSRO_ForwardModel(nlevP, &
                                  nlevq, &
                                  za, &
                                  zb, &
                                  pressure, &
                                  humidity, &
                                  GPSRO_pseudo_ops, &
                                  GPSRO_vert_interp_ops, &
                                  nobs, &
                                  zobs, &
                                  RO_Rad_Curv, &
                                  Latitude, &
                                  RO_geoid_und, &
                                  ycalc, &
                                  BAErr)

INTEGER, INTENT(IN)            :: nlevP                  ! no. of p levels in state vec.
INTEGER, INTENT(IN)            :: nlevq                  ! no. of theta levels
REAL(kind_real), INTENT(IN)    :: za(1:nlevP)            ! heights of rho levs
REAL(kind_real), INTENT(IN)    :: zb(1:nlevQ)            ! heights of theta levs
REAL(kind_real), INTENT(IN)    :: pressure(1:nlevP)      ! Model background pressure
REAL(kind_real), INTENT(IN)    :: humidity(1:nlevQ)      ! Model background specific humidity
LOGICAL, INTENT(IN)            :: GPSRO_pseudo_ops       ! Option: Use pseudo-levels in vertical interpolation?
LOGICAL, INTENT(IN)            :: GPSRO_vert_interp_ops  ! Option: Use ln(p) for vertical interpolation? (rather than exner)
INTEGER, INTENT(IN)            :: nobs                   ! Number of observations in the profile
REAL(kind_real), INTENT(IN)    :: zobs(1:nobs)           ! Impact parameter for the obs
REAL(kind_real), INTENT(IN)    :: RO_Rad_Curv            ! Earth's radius of curvature for these observations
REAL(kind_real), INTENT(IN)    :: Latitude               ! Latitude of this profile
REAL(kind_real), INTENT(IN)    :: RO_geoid_und           ! Undulation - difference between the geoid and the ellipsoid
REAL(kind_real), INTENT(INOUT) :: ycalc(1:nobs)          ! Model forecast of the observations
LOGICAL, INTENT(OUT)           :: BAErr                  ! Was an error encountered during the calculation?
! 
! Things that may need to be output, as they are used by the TL/AD calculation
! 
REAL(kind_real), ALLOCATABLE :: z_pseudo(:)        ! Heights of the pseudo levels       | Allocated by
REAL(kind_real), ALLOCATABLE :: N_pseudo(:)        ! Refractivity on the pseudo levels  | Ops_GPSRO_refrac
INTEGER                      :: nb_pseudo          ! Number of pseudo levels
REAL(kind_real)              :: T(1:nlevq)         ! Temperature on model levels
REAL(kind_real), ALLOCATABLE :: nr(:)              ! Model calculation of impact parameters
REAL(kind_real)              :: Refmodel(1:nlevq)  ! model refractivity on theta levels
! 
! Local parameters
! 
integer, parameter           :: max_string = 800  ! Length of strings
character(len=*), parameter  :: myname_ = "Ops_GPSRO_ForwardModel"
!
! Local variables
! 
INTEGER                      :: nstate            ! no. of levels in state vec.
INTEGER                      :: num_pseudo        ! Number of levels, including pseudo levels
INTEGER                      :: nb                ! no. of non-pseudo levs
REAL(kind_real)              :: x(1:nlevP+nlevQ)  ! state vector
character(max_string)        :: err_msg           ! Error message to be output

! The model data must be on a staggered grid, with nlevp = nlevq+1
IF (nlevP /= nlevQ + 1) THEN
    write(err_msg,*) myname_ // ':' // ' Data must be on a staggered grid nlevp, nlevq = ', nlevp, nlevq
    call fckit_log % warning(err_msg)
    write(err_msg,*) myname_ // ':' // ' error: number of levels inconsistent!'
    call abor1_ftn(err_msg)
END IF

nstate = nlevP + nlevq
nb = nlevq
x(1:nlevP) = pressure
x(nlevP+1:nstate) = humidity

IF (GPSRO_pseudo_ops) THEN
    num_pseudo = 2 * nlevq - 1
ELSE
    num_pseudo = nlevq
END IF
ALLOCATE(nr(1:num_pseudo))

BAErr = .FALSE.

CALL Ops_GPSRO_refrac (nstate,   &
                       nlevP,    &
                       nb,       &
                       nlevq,    &
                       za,       &
                       zb,       &
                       x,        &
                       GPSRO_pseudo_ops, &
                       GPSRO_vert_interp_ops, &
                       BAerr,    &
                       Refmodel, &
                       T,        &
                       z_pseudo, &
                       N_pseudo, &
                       nb_pseudo)

! no point proceeding further if ...
IF (.NOT. BAerr) THEN
    ! Pseudo levels
    IF (GPSRO_pseudo_ops) THEN
        !  2.  Calculate the refractive index * radius on theta model levels (or model impact parameter)
        CALL Ops_GPSROcalc_nr (z_pseudo,     &           ! geopotential heights of pseudo levels
                               nb_pseudo,    &           ! number of model+pseudo-levels
                               RO_Rad_Curv,  &           ! radius of curvature of earth at observation
                               Latitude,     &           ! latitude at observation
                               RO_geoid_und, &           ! geoid undulation above WGS-84
                               n_pseudo,     &           ! refractivity of model on model+pseudo levels
                               nr)                       ! Calculated model impact parameters

        !  3.  Calculate model bending angle on observation impact parameters
        CALL Ops_GPSROcalc_alpha (nobs,      &      ! size of ob. vector
                                  nb_pseudo, &      ! no. of refractivity levels
                                  zobs,      &      ! obs impact parameters
                                  n_pseudo,  &      ! refractivity values on model+pseudo levels
                                  nr,        &      ! index * radius product
                                  ycalc)            ! forward modelled bending angle
    ! Model levels only
    ELSE
        !  2.  Calculate the refractive index * radius on theta model levels (or model impact parameter)
        CALL Ops_GPSROcalc_nr (zb,           &           ! geopotential heights of model levels
                               nb,           &           ! number of levels in zb
                               RO_Rad_Curv,  &           ! radius of curvature of earth at observation
                               Latitude,     &           ! latitude at observation
                               RO_geoid_und, &           ! geoid undulation above WGS-84
                               Refmodel,     &           ! refractivity of model on model levels
                               nr)                       ! Calculated model impact parameters

        !  3.  Calculate model bending angle on observation impact parameters
        CALL Ops_GPSROcalc_alpha (nobs,     &      ! size of ob. vector
                                  nb,       &      ! no. of refractivity levels
                                  zobs,     &      ! obs impact parameters
                                  Refmodel, &      ! refractivity values on model levels
                                  nr,       &      ! index * radius product
                                  ycalc)           ! forward modelled bending angle
    END IF
END IF

END SUBROUTINE Ops_GPSRO_ForwardModel

end module ufo_gnssro_bendmetoffice_mod
