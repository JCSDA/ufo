! (C) Copyright 2025 Met Office
! 
! this software is licensed under the terms of the apache licence version 2.0
! which can be obtained at http://www.apache.org/licenses/license-2.0.

!> The main Fortran module for implementing the GNSS-RO onedvar check

module ufo_refractivityonedvarcheck_mod

use, intrinsic :: iso_c_binding
use fckit_configuration_module, only: fckit_configuration
use logger_mod, only : oops_log
use fckit_exception_module, only: fckit_exception
use iso_c_binding
use kinds
use missing_values_mod
use obsspace_mod
use oops_variables_mod
use ufo_geovals_mod
use ufo_vars_mod
use ufo_gnssro_refmetoffice_mod
use ufo_refractivityonedvarcheck_utils_mod, only: &
    deallocate_singlerefob, allocate_singlerefob, allocate_singlerefbg, &
    deallocate_singlerefbg, singlerefbg_type, singlerefob_type
use ufo_gnssroonedvarcheck_get_bmatrix_mod, only: bmatrix_type
use ufo_refractivityonedvarcheck_do1dvar_mod, only: Refractivity_Do1DVar
use ufo_roobserror_utils_mod, only: &
    rmatrix_type, ufo_roobserror_findnearest_rmatrix, ufo_roobserror_getrmatrix
use ufo_utils_mod, only: Ops_RealSortQuick, find_unique
use gnssro_mod_transform, only: geometric2geop

implicit none
private
public :: ufo_refractivityonedvarcheck_create
public :: ufo_refractivityonedvarcheck_delete
public :: ufo_refractivityonedvarcheck_apply

type, public :: ufo_refractivityonedvarcheck
  character(len=800)        :: qcname            !< name of the filter
  type(c_ptr)               :: obsdb             !< pointer to the observation space
  type(fckit_configuration) :: conf              !< contents of the yaml file
  integer(c_int)            :: onedvarflag       !< flag uased by the qc manager for a 1D-var check
  logical                   :: capsupersat       !< Whether to remove super-saturation (wrt ice?)
  real(kind_real)           :: cost_funct_test   !< Threshold value for the cost function convergence test
  real(kind_real)           :: Delta_ct2         !< Threshold used in calculating convergence
  real(kind_real)           :: Delta_factor      !< Threshold used in calculating convergence
  integer                   :: n_iteration_test  !< Maximum number of iterations in the 1DVar
  real(kind_real)           :: OB_test           !< Threshold for the O-B throughout the profile
  real(kind_real)           :: y_test            !< Threshold on distance between observed and solution refractivities
  real(kind_real)           :: minval_ytest      !< Minimum value for RHS of y-test
  real(kind_real)           :: maxval_ytest      !< Maximum value for RHS of y-test
  character(len=800)        :: bmatrix_filename  !< Location of the B-matrix file
  logical                   :: pseudo_ops        !< Whether to use pseudo levels in forward operator
  logical                   :: vert_interp_ops   !< Whether to use ln(p) or exner in vertical interpolation
  real(kind_real)           :: min_temp_grad     !< The minimum vertical temperature gradient allowed
  type(rmatrix_type), allocatable :: Rmatrix_list(:) !< All the R matrices read in from the file
  integer                   :: r_num_sats        !< The number of R matrices in the list
end type ufo_refractivityonedvarcheck

! ------------------------------------------------------------------------------
contains
! ------------------------------------------------------------------------------
!> Setup the main GNSS-RO onedvar object in Fortran
!!
!! \details Makes a call to the main setup routine.
!!
!! \author Met Office
!!
!! \date 09/06/2020: Created
!!
subroutine ufo_refractivityonedvarcheck_create( &
  self, &
  obsspace, &
  bmatrix_filename, &
  rmatrix_filename, &
  capsupersat, &
  cost_funct_test, &
  Delta_ct2, &
  Delta_factor, &
  min_temp_grad, &
  n_iteration_test, &
  OB_test, &
  pseudo_ops, &
  vert_interp_ops, &
  y_test, &
  minval_ytest, &
  maxval_ytest, &
  onedvarflag &
)

  implicit none
  type(ufo_refractivityonedvarcheck), intent(inout) :: self        !< refractivityonedvarcheck main object
  type(c_ptr), value, intent(in)              :: obsspace          !< observation database pointer
  character(len=*), intent(in)                :: bmatrix_filename  !< Location of the B-matrix file
  character(len=*), intent(in)                :: rmatrix_filename  !< Location of the R-matrix file
  logical(c_bool), intent(in)                 :: capsupersat       !< Whether to remove super-saturation (wrt ice?)
  real(c_float), intent(in)                   :: cost_funct_test   !< Threshold value for the cost function convergence test
  real(c_float), intent(in)                   :: Delta_ct2         !< Threshold used in calculating convergence
  real(c_float), intent(in)                   :: Delta_factor      !< Threshold used in calculating convergence
  real(c_float), intent(in)                   :: min_temp_grad     !< The minimum vertical temperature gradient allowed
  integer(c_int), intent(in)                  :: n_iteration_test  !< Maximum number of iterations in the 1DVar
  real(c_float), intent(in)                   :: OB_test           !< Threshold for the O-B throughout the profile
  logical(c_bool), intent(in)                 :: pseudo_ops        !< Whether to use pseudo levels in forward operator
  logical(c_bool), intent(in)                 :: vert_interp_ops   !< Whether to use ln(p) or exner in vertical interpolation
  real(c_float), intent(in)                   :: y_test            !< Threshold on distance between observed and solution refractivities
  real(c_float), intent(in)                   :: minval_ytest      !< Minimum value for RHS of y-test
  real(c_float), intent(in)                   :: maxval_ytest      !< Maximum value for RHS of y-test
  integer(c_int), intent(in)                  :: onedvarflag       !< flag for qc manager

  character(len=800) :: message
  integer :: i
  integer, parameter :: Rmax_num = 2000   !< Max number of R matrices to read in.  Hard-wire to a large number

  self % obsdb = obsspace
  self % onedvarflag = onedvarflag

  self % bmatrix_filename = bmatrix_filename
  self % capsupersat = capsupersat
  self % cost_funct_test = cost_funct_test
  self % Delta_ct2 = Delta_ct2
  self % Delta_factor = Delta_factor
  self % min_temp_grad = min_temp_grad
  self % n_iteration_test = n_iteration_test
  self % OB_test = OB_test
  self % pseudo_ops = pseudo_ops
  self % vert_interp_ops = vert_interp_ops
  self % y_test = y_test
  self % minval_ytest = minval_ytest
  self % maxval_ytest = maxval_ytest

  write(message, '(A)') 'GNSS-RO 1D-Var check: input parameters are:'
  call oops_log % debug(message)
  write(message, '(2A)') 'bmatrix_filename = ', bmatrix_filename
  call oops_log % debug(message)
  write(message, '(A,L1)') 'capsupersat = ', capsupersat
  call oops_log % debug(message)
  write(message, '(A,F16.8)') 'cost_funct_test = ', cost_funct_test
  call oops_log % debug(message)
  write(message, '(A,F16.8)') 'Delta_ct2 = ', Delta_ct2
  call oops_log % debug(message)
  write(message, '(A,F16.8)') 'Delta_factor = ', Delta_factor
  call oops_log % debug(message)
  write(message, '(A,F16.8)') 'min_temp_grad = ', min_temp_grad
  call oops_log % debug(message)
  write(message, '(A,I7)') 'n_iteration_test = ', n_iteration_test
  call oops_log % debug(message)
  write(message, '(A,F16.8)') 'OB_test = ', OB_test
  call oops_log % debug(message)
  write(message, '(A,L1)') 'pseudo_ops = ', pseudo_ops
  call oops_log % debug(message)
  write(message, '(A,L1)') 'vert_interp_ops = ', vert_interp_ops
  call oops_log % debug(message)
  write(message, '(A,F16.8)') 'y_test = ', y_test
  call oops_log % debug(message)
  write(message, '(A,F16.8)') 'minval_ytest = ', minval_ytest
  call oops_log % debug(message)
  write(message, '(A,F16.8)') 'maxval_ytest = ', maxval_ytest
  call oops_log % debug(message)

  write(message, '(2A)') 'Attempting to read rmatrix file: ', TRIM(rmatrix_filename)
  call oops_log % debug(message)

! Read in R matrix data
  call ufo_roobserror_getrmatrix(Rmax_num,            &  ! Max number of R matrices to read in
                                 rmatrix_filename,    &  ! The name of the file to be read in
                                 self % rmatrix_list, &  ! List of all R matrices to use
                                 self % R_num_sats)      ! Number of R matrices read in

end subroutine ufo_refractivityonedvarcheck_create

! ------------------------------------------------------------------------------
!> Delete the main GNSS-RO onedvar object in Fortran
!!
!! \author Met Office
!!
!! \date 09/06/2020: Created
!!
subroutine ufo_refractivityonedvarcheck_delete(self)

  implicit none
  type(ufo_refractivityonedvarcheck), intent(inout) :: self !< refractivityonedvarcheck main object

end subroutine ufo_refractivityonedvarcheck_delete

! ------------------------------------------------------------------------------
!> The main routine that applys the GNSS-RO onedvar filter
!!
!! \details Heritage : 
!!
!! This routine is called from the c++ apply method.  The filter performs 
!! a 1D-Var minimization
!!
!! \author Met Office
!!
!! \date 09/06/2020: Created
!!
subroutine ufo_refractivityonedvarcheck_apply(self, geovals, apply)

  implicit none

  ! Subroutine arguments
  type(ufo_refractivityonedvarcheck), intent(inout) :: self   !< refractivityonedvarcheck main object
  type(ufo_geovals), intent(in)              :: geovals       !< model values at observation space
  logical, intent(in)                        :: apply(:)      !< qc manager flags

  ! Local parameters
  logical, parameter :: verboseOutput = .FALSE.               ! Whether to output extra debugging information

  ! Local variables
  integer :: nobs                                             ! Number of observations to be processed
  type(ufo_geoval), pointer          :: q                     ! Model background values of specific humidity
  type(ufo_geoval), pointer          :: prs                   ! Model background values of air pressure
  type(ufo_geoval), pointer          :: theta_heights         ! Model heights of levels containing specific humidity
  type(ufo_geoval), pointer          :: rho_heights           ! Model heights of levels containing air pressure
  type(singlerefbg_type)             :: Back                  ! Model background fields
  type(singlerefob_type)             :: Ob                    ! The profile of observations
  type(bmatrix_type)                 :: b_matrix              ! Background-error covariance matrix
  real(kind_real), allocatable       :: obsLat(:)             ! Latitude of the observation
  real(kind_real), allocatable       :: obsLon(:)             ! Longitude of the observation
  real(kind_real), allocatable       :: obs_height(:)         ! Impact parameter of the observation
  real(kind_real), allocatable       :: obs_refractivity(:)   ! Observed refractivity
  real(kind_real), allocatable       :: obs_err(:)            ! Observation error, taken from a previous filter
  integer, allocatable               :: qc_flags(:)           ! QC flags to be updated
  integer, allocatable               :: obsSatid(:)           ! Satellite identifier for each observation
  integer, allocatable               :: obsOrigC(:)           ! Originating centre for each observation
  integer(c_size_t), allocatable     :: record_number(:)      ! Number used to identify unique profiles in the data
  real(kind_real), allocatable       :: sort_key(:)           ! Key for the sorting (based on record number and impact parameter)
  integer, allocatable               :: index_vals(:)         ! Indices of sorted observation
  integer, allocatable               :: unique(:)             ! Set of unique profile numbers
  integer                            :: start_point           ! Starting index of the current profile
  integer                            :: current_point         ! Ending index of the current profile
  integer                            :: iprofile              ! Loop variable, profile number
  integer                            :: iobs                  ! Loop variable, observation number
  integer                            :: nobs_profile          ! Number of observations in the profile
  character(len=800)                 :: Message               ! Message to be output
  logical                            :: BAerr                 ! Has there been an error in the refractivity calculation?
  integer                            :: iband                 ! Selected latitude band of the B-matrix
  integer                            :: iseason               ! Selected season of the B-matrix
  integer                            :: ipoint                ! Loop variable, observation point
  real(kind_real)                    :: dfs                   ! Degrees of freedom for signal in profile
  real(kind_real)                    :: O_Bdiff               ! Average RMS(O-B) for profile
  real(kind_real), allocatable       :: Tb(:)                 ! Calculated background temperature (derived from p,q)
  real(kind_real), allocatable       :: Ts(:)                 ! 1DVar solution temperature
  integer                            :: ilevel                ! Loop variable, level number
  real(kind_real)                    :: missing               ! Missing data indicator (for reals)
  type(rmatrix_type)                 :: rmatrix               ! The current R-matrix information
!
! Diagnostics to push back to the obs-space
!
  integer, allocatable               :: indices(:)            ! The indices of the diagnostic elements to be updated
  integer, allocatable               :: niter(:)              ! Number of iterations required to converge
  real(kind_real), allocatable       :: initial_cost(:)       ! Initial cost-function value
  real(kind_real), allocatable       :: final_cost(:)         ! Final cost-function value
  real(kind_real), allocatable       :: dfs_list(:)           ! Degrees of freedom for signal

  write(Message,*) "TRACE: ufo_refractivityonedvarcheck_refractivityonedvarcheck_apply: begin"
  call oops_log % debug(Message)

  ! Get the obs-space information
  nobs = obsspace_get_nlocs(self % obsdb)

  write(Message,*) "nobs = ", nobs
  call oops_log % debug(Message)

  allocate(obsLon(nobs))
  allocate(obsLat(nobs))
  allocate(obs_height(nobs))
  allocate(obs_refractivity(nobs))
  allocate(record_number(nobs))
  allocate(obsSatid(nobs))
  allocate(obsOrigC(nobs))
  allocate(qc_flags(nobs))
  allocate(obs_err(nobs))
! Allocate room for diagnostics
  allocate(niter(nobs))
  allocate(initial_cost(nobs))
  allocate(final_cost(nobs))
  allocate(dfs_list(nobs))

  call obsspace_get_db(self % obsdb, "MetaData", "longitude", obsLon)
  call obsspace_get_db(self % obsdb, "MetaData", "latitude", obsLat)
  call obsspace_get_recnum(self % obsdb, record_number)

  call obsspace_get_db(self % obsdb, "FortranQC", "atmosphericRefractivity", qc_flags)
  call obsspace_get_db(self % obsdb, "MetaData", "height", obs_height)
  call obsspace_get_db(self % obsdb, "ObsValue", "atmosphericRefractivity", obs_refractivity)
  call obsspace_get_db(self % obsdb, "FortranERR", "atmosphericRefractivity", obs_err)
  call obsspace_get_db(self % obsdb, "MetaData", "satelliteIdentifier", obsSatid)
  call obsspace_get_db(self % obsdb, "MetaData", "dataProviderOrigin", obsOrigC)

! Convert geometric height to geopotential height
  do iobs = 1, nobs
    if (obs_height(iobs) /= missing_value(obs_height(iobs))) then
      call geometric2geop(obsLat(iobs), obs_height(iobs), obs_height(iobs))
    end if
  end do

  ! get variables from geovals
  call ufo_geovals_get_var(geovals, var_q, q)               ! specific humidity
  call ufo_geovals_get_var(geovals, var_prsi, prs)          ! pressure
  call ufo_geovals_get_var(geovals, var_z, theta_heights)   ! Geopotential height of the normal model levels
  call ufo_geovals_get_var(geovals, var_zi, rho_heights)    ! Geopotential height of the pressure levels

  ! Read in the B-matrix and background profiles
  call b_matrix % get(self % bmatrix_filename, prs % nval, q % nval)
  call allocate_singlerefbg(Back, prs % nval, q % nval)
  allocate(Tb(q % nval), Ts(q % nval))

  ! Read through the record numbers in order to find a profile of observations
  ! Each profile shares the same record number
  missing = missing_value(missing)
  allocate(sort_key(nobs))
  do iobs = 1, nobs
    if (obs_height(iobs) == missing) then
      sort_key(iobs) = record_number(iobs) + 0.9
    else
      sort_key(iobs) = record_number(iobs) + 0.5 * (obs_height(iobs) / MAXVAL(obs_height))
    end if
  end do

  call Ops_RealSortQuick(sort_key, index_vals)
  call find_unique(record_number, unique)
  WRITE (Message, '(A,I0)') 'Number of unique profiles ', size(unique)
  call oops_log % debug(Message)

  ! For every profile that we have found, perform a 1DVar minimisation
  current_point = 1
  do iprofile = 1, size(unique)
    start_point = current_point
    iobs = index_vals(start_point)
    WRITE (Message, '(A,I0)') 'ObNumber ', iprofile
    call oops_log % info(Message)
    WRITE (Message, '(A,F12.2)') 'Latitude ', obsLat(iobs)
    call oops_log % info(Message)
    WRITE (Message, '(A,F12.2)') 'Longitude ', obsLon(iobs)
    call oops_log % info(Message)
    WRITE (Message, '(A,I0)') 'Processing centre ', obsOrigC(iobs)
    call oops_log % info(Message)
    WRITE (Message, '(A,I0)') 'Sat ID ', obsSatid(iobs)
    call oops_log % info(Message)

    ! Work out which observations belong to the current profile
    do current_point = start_point, nobs
      if (unique(iprofile) /= record_number(index_vals(current_point))) exit
    end do

    ! Load the geovals into the background structure
    ! Reverse the order of the geovals, since this routine (and the forward
    ! operator) works bottom-to-top
    Back % za(:) = rho_heights % vals(prs%nval:1:-1, iobs)
    Back % zb(:) = theta_heights % vals(q%nval:1:-1, iobs)
    Back % p(:) = prs % vals(prs % nval:1:-1, iobs)
    Back % q(:) = q % vals(q%nval:1:-1, iobs)
    
    ! Allocate the observations structure
    nobs_profile = current_point - start_point
    call allocate_singlerefob(Ob, nobs_profile, prs % nval, q % nval)

    ! Load the observations information into the obsevations structure
    Ob % id = record_number(index_vals(iobs))
    Ob % latitude = obsLat(index_vals(iobs))
    Ob % longitude = obsLon(index_vals(iobs))
    Ob % niter = 0  ! We haven't yet run 1DVar
    Ob % jcost = missing_value(Ob % jcost)

    Ob % refractivity(:) % value = obs_refractivity(index_vals(start_point:current_point-1))
    Ob % refractivity(:) % oberr = obs_err(index_vals(start_point:current_point-1))
    Ob % z(:) % value            = obs_height(index_vals(start_point:current_point-1))
    Ob % qc_flags(:)             = qc_flags(index_vals(start_point:current_point-1))

    !--------------------------------------------------------
    ! 6. Choose R matrix depending on satid, origctr and latitude
    ! No interpolation between matrices to match old code
    !--------------------------------------------------------
    CALL ufo_roobserror_findnearest_rmatrix(obsSatid(iobs), &
                                            obsOrigC(iobs), &
                                            obsLat(iobs), &
                                            self % R_num_sats, &
                                            self % RMatrix_list, &
                                            RMatrix)

    ! Choose the latitude band and season of the B-matrix information
    iseason = 1    ! Temporary -only one season at present!
    iband = 1
    DO
      IF (b_matrix % band_up_lim(iband) > Ob % latitude .OR. &
          iband == b_matrix % nband) EXIT
      iband = iband + 1
    END DO

!   Call the code to set up the 1D-Var calculation
    call Refractivity_Do1DVar( &
      prs % nval,              &   ! Number of pressure levels
      q % nval,                &   ! Number of specific humidity levels
      b_matrix % inverse(iseason,iband,:,:), &  ! Inverse of the b-matrix
      b_matrix % sigma(iseason,iband,:), &      ! Standard deviations of the b-matrix
      Back,                    &   ! Structure containing the model background information
      Ob,                      &   ! Structure containing the observation information
      RMatrix,                 &   ! Structure containing the observation uncertainties
      self % pseudo_ops,       &   ! Whether to use pseudo-levels in calculation
      self % vert_interp_ops,  &   ! Whether to interpolate using ln(p) or exner
      self % min_temp_grad,    &   ! Minimum vertical temperature gradient allowed
      self % cost_funct_test,  &   ! Threshold value for the cost function convergence test
      self % y_test,           &   ! Threshold value for the yobs-ysol tes
      self % n_iteration_test, &   ! Maximum number of iterations
      self % Delta_factor,     &   ! Delta
      self % Delta_ct2,        &   ! Delta observations
      self % OB_test,          &   ! Threshold value for the O-B test
      self % minval_ytest,     &   ! Minimum value for RHS of y-test
      self % maxval_ytest,     &   ! Maximum value for RHS of y-test
      self % capsupersat,      &   ! Whether to remove super-saturation
      BAerr,                   &   ! Whether there are errors in the refractivity calculation
      Tb,                      &   ! Calculated background temperature
      Ts,                      &   ! 1DVar solution temperature
      O_Bdiff,                 &   ! Difference between observations and background for profile
      DFS                      &   ! Estimated degrees of freedom for signal
    )

    ! Flag bad profiles
    do ipoint = 0, nobs_profile-1
      if (qc_flags(start_point + ipoint) > 0) then
        ! Do nothing, since the data are already flagged
      else if (Ob % refractivity(ipoint+1) % PGEFinal > 0.5) then
        qc_flags(start_point + ipoint) = self % onedvarflag
        Ob % qc_flags(ipoint+1) = self % onedvarflag
      end if
    end do

    if (verboseOutput) then
      do ipoint = 0, nobs_profile-1, 20
          write(Message,'(20I5)') qc_flags( &
              index_vals(start_point+ipoint:min(start_point+ipoint+19, current_point-1)))
          call oops_log % debug(Message)
      end do
      do ipoint = 0, nobs_profile-1, 10
          write(Message,'(10E16.5)') obs_refractivity( &
              index_vals(start_point+ipoint:min(start_point+ipoint+9, current_point-1)))
          call oops_log % debug(Message)
      end do
    end if

    ! Save the diagnostic information
    allocate(indices(1:nobs_profile))
    do ipoint = 0, nobs_profile-1
      indices(ipoint+1) = index_vals(start_point+ipoint)
    end do
    niter(indices) = Ob % niter
    initial_cost(indices) = O_Bdiff
    final_cost(indices) = Ob % jcost
    dfs_list(indices) = DFS
    deallocate(indices)

    call deallocate_singlerefob(Ob)
  end do
  call obsspace_put_db(self % obsdb, "FortranQC", "atmosphericRefractivity", qc_flags)

  ! Save the diagnostics to the obs-space
  call obsspace_put_db(self % obsdb, "OneDVarDiags", "nIter", niter)
  call obsspace_put_db(self % obsdb, "OneDVarDiags", "initialCost", initial_cost)
  call obsspace_put_db(self % obsdb, "OneDVarDiags", "finalCost", final_cost)
  call obsspace_put_db(self % obsdb, "OneDVarDiags", "DFS", dfs_list)

  write(Message,*) "TRACE: ufo_refractivityonedvarcheck_refractivityonedvarcheck_apply: begin"
  call oops_log % debug(Message)

end subroutine ufo_refractivityonedvarcheck_apply

end module ufo_refractivityonedvarcheck_mod
