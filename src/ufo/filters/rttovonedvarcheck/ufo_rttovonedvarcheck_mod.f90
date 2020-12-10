! (C) Copyright 2017-2020 Met Office
! 
! this software is licensed under the terms of the apache licence version 2.0
! which can be obtained at http://www.apache.org/licenses/license-2.0.

!> Thae main Fortran module for implementing the rttov onedvar check

module ufo_rttovonedvarcheck_mod

use fckit_configuration_module, only: fckit_configuration
use fckit_log_module, only : fckit_log
use iso_c_binding
use kinds
use missing_values_mod
use obsspace_mod
use oops_variables_mod
use ufo_geovals_mod
use ufo_radiancerttov_mod
use ufo_rttovonedvarcheck_bmatrix_mod
use ufo_rttovonedvarcheck_constants_mod
use ufo_rttovonedvarcheck_minimize_utils_mod
use ufo_rttovonedvarcheck_minimize_newton_mod
use ufo_rttovonedvarcheck_minimize_ml_mod
use ufo_rttovonedvarcheck_ob_mod
use ufo_rttovonedvarcheck_obs_mod
use ufo_rttovonedvarcheck_pcemis_mod
use ufo_rttovonedvarcheck_profindex_mod
use ufo_rttovonedvarcheck_rmatrix_mod
use ufo_rttovonedvarcheck_rsubmatrix_mod
use ufo_rttovonedvarcheck_utils_mod
use ufo_vars_mod

implicit none
private
public :: ufo_rttovonedvarcheck_create
public :: ufo_rttovonedvarcheck_delete
public :: ufo_rttovonedvarcheck_apply

! ------------------------------------------------------------------------------
contains
! ------------------------------------------------------------------------------
!> Setup the main rttov onedvar object in Fortran
!!
!! \details Makes a call to the main setup routine.
!!
!! \author Met Office
!!
!! \date 09/06/2020: Created
!!
subroutine ufo_rttovonedvarcheck_create(self, obsspace, f_conf, channels, &
                                        onedvarflag, passflag)

  implicit none
  type(ufo_rttovonedvarcheck), intent(inout) :: self         !< rttovonedvarcheck main object
  type(c_ptr), value, intent(in)             :: obsspace     !< observation database pointer
  type(fckit_configuration), intent(in)      :: f_conf       !< yaml file contents
  integer(c_int), intent(in)                 :: channels(:)  !< all channels that can be used in 1D-Var
  integer(c_int), intent(in)                 :: onedvarflag  !< flag from qc flags
  integer(c_int), intent(in)                 :: passflag     !< pass flag from qc flags

  self % obsdb = obsspace
  self % conf = f_conf
  self % onedvarflag = onedvarflag
  self % passflag = passflag

  call ufo_rttovonedvarcheck_setup(self, channels) ! from init

end subroutine ufo_rttovonedvarcheck_create

! ------------------------------------------------------------------------------
!> Delete the main rttov onedvar object in Fortran
!!
!! \author Met Office
!!
!! \date 09/06/2020: Created
!!
subroutine ufo_rttovonedvarcheck_delete(self)

  implicit none
  type(ufo_rttovonedvarcheck), intent(inout) :: self !< rttovonedvarcheck main object

  if (allocated(self % retrieval_variables)) deallocate(self % retrieval_variables)
  if (allocated(self % channels))            deallocate(self % channels)

end subroutine ufo_rttovonedvarcheck_delete

! ------------------------------------------------------------------------------
!> The main routine that applys the rttov onedvar filter
!!
!! \details Heritage : Ops_SatRad_Do1DVar_RTTOV12.f90
!!
!! This routine is called from the c++ apply method.  The filter performs 
!! a 1D-Var minimization using rttov
!!
!! \author Met Office
!!
!! \date 09/06/2020: Created
!!
subroutine ufo_rttovonedvarcheck_apply(self, vars, retrieval_vars, geovals, apply)

  implicit none
  type(ufo_rttovonedvarcheck), intent(inout) :: self     !< rttovonedvarcheck main object
  type(oops_variables), intent(in)           :: vars     !< channels for 1D-Var
  type(oops_variables), intent(in)           :: retrieval_vars !< retrieval variables for 1D-Var
  type(ufo_geovals), intent(in)              :: geovals  !< model values at observation space
  logical, intent(in)                        :: apply(:) !< qc manager flags

  type(ufo_rttovonedvarcheck_obs)        :: obs            ! data for all observations read from db
  type(ufo_rttovonedvarcheck_bmatrix)    :: full_bmatrix   ! full bmatrix read from file
  type(ufo_geovals)                      :: local_geovals  ! geoval for one observation
  type(ufo_rttovonedvarcheck_ob)         :: ob             ! observation data for a single observation
  type(ufo_rttovonedvarcheck_profindex)  :: prof_index     ! index for mapping geovals to 1d-var state profile
  type(ufo_rttovonedvarcheck_rmatrix)    :: full_rmatrix   ! full r_matrix read from file
  type(ufo_rttovonedvarcheck_rsubmatrix) :: r_submatrix    ! r_submatrix object
  type(ufo_geovals)                      :: hofxdiags      ! hofxdiags containing jacobian
  type(ufo_rttovonedvarcheck_pcemis), target :: IR_pcemis  ! Infrared principal components object
  character(len=max_string)          :: sensor_id
  character(len=max_string)          :: var
  character(len=max_string)          :: varname
  character(len=max_string)          :: message
  integer                            :: jvar, ivar, jobs, band, ii ! counters
  integer                            :: nchans_used      ! counter for number of channels used for an ob
  integer                            :: jchans_used
  integer                            :: fileunit        ! unit number for reading in files
  integer                            :: apply_count
  integer                            :: nprofelements   ! number of elements in 1d-var state profile
  integer, allocatable               :: fields_in(:)
  real(kind_real)                    :: missing         ! missing value
  real(kind_real)                    :: t1, t2          ! timing
  real(kind_real), allocatable       :: b_matrix(:,:)   ! 1d-var profile b matrix
  real(kind_real), allocatable       :: b_inverse(:,:)  ! inverse for each 1d-var profile b matrix
  real(kind_real), allocatable       :: b_sigma(:)      ! b_matrix diagonal error
  logical                            :: file_exists     ! check if a file exists logical
  logical                            :: onedvar_success
  logical                            :: cloud_retrieval = .false.
  type(ufo_radiancerttov)            :: rttov_simobs

  ! ------------------------------------------
  ! 1. Setup
  ! ------------------------------------------
  missing = missing_value(missing)

  ! Setup rttov simobs
  call rttov_simobs % setup(self % conf, self % channels)

  ! Setup IR emissivity - if needed
  if (self % pcemiss) then
    if (len(self % EmisAtlas) > 4) then
      call IR_pcemis % setup(self % EmisEigVecPath, self % EmisAtlas)
    else
      call IR_pcemis % setup(self % EmisEigVecPath)
    end if
  end if

  ! Setup full B matrix object
  call full_bmatrix % setup(self % retrieval_variables, self % b_matrix_path, &
                            self % qtotal)
  
  ! Setup full R matrix object
  call full_rmatrix % setup(self % r_matrix_path)

  ! Check if cloud retrievals needed
  do ii = 1, size(self % retrieval_variables)
    if (trim(self % retrieval_variables(ii)) == "cloud_top_pressure") then
      write(*,*) "Simple cloud is part of the state vector"
      cloud_retrieval = .true.
    end if
  end do

  ! Create profile index for mapping 1d-var profile to b-matrix
  call prof_index % setup(full_bmatrix, self%nlevels)

  ! Read in observation data from obsspace
  call obs % setup(self, prof_index % nprofelements, geovals, vars, IR_pcemis)

  ! Initialize data arrays
  allocate(b_matrix(prof_index % nprofelements,prof_index % nprofelements))
  allocate(b_inverse(prof_index % nprofelements,prof_index % nprofelements))
  allocate(b_sigma(prof_index % nprofelements))

  ! ------------------------------------------
  ! 2. Beginning main observation loop
  ! ------------------------------------------
  write(*,*) "Beginning loop over observations: ",self%qcname
  apply_count = 0
  obs_loop: do jobs = 1, obs % iloc
    if (apply(jobs)) then

      apply_count = apply_count + 1
      write(*,*) "starting obs number    ",jobs
      !---------------------------------------------------
      ! 2.1 Setup Jb terms
      !---------------------------------------------------
      ! create one ob geovals from full all obs geovals
      call ufo_geovals_copy_one(local_geovals, geovals, jobs)
      call ufo_rttovonedvarcheck_check_geovals(local_geovals, &
                          prof_index, obs % surface_type(jobs))

      ! create b matrix arrays for this single observation location
      call full_bmatrix % reset( obs % lat(jobs), & ! in
                    b_matrix, b_inverse, b_sigma  ) ! out

      !---------------------------------------------------
      ! 2.2 Setup Jo terms
      !---------------------------------------------------
      ! Channel selection based on previous filters flags
      nchans_used = 0
      do jvar = 1, self%nchans
        if( obs % QCflags(jvar,jobs) == self % passflag ) then
          nchans_used = nchans_used + 1
        end if
      end do
      if (nchans_used == 0) then
        write(message, *) "No channels selected for observation number ", &
               jobs, " : skipping"
        call fckit_log % info(message)
        cycle obs_loop
      end if

      ! setup ob data for this observation
      call ob % setup(nchans_used, prof_index % nprofelements, self % nchans)
      ob % forward_mod_name = self % forward_mod_name
      ob % latitude = obs % lat(jobs)
      ob % longitude = obs % lon(jobs)
      ob % elevation = obs % elevation(jobs)
      ob % sensor_zenith_angle = obs % sat_zen(jobs)
      ob % sensor_azimuth_angle = obs % sat_azi(jobs)
      ob % solar_zenith_angle = obs % sol_zen(jobs)
      ob % solar_azimuth_angle = obs % sol_azi(jobs)
      ob % channels_all = self % channels
      ob % surface_type = obs % surface_type(jobs)
      ob % retrievecloud = cloud_retrieval
      ob % pcemis => IR_pcemis
      ob % calc_emiss = obs % calc_emiss(jobs)
      if(self % RTTOV_mwscattSwitch) ob % mwscatt = .true.
      if(self % RTTOV_usetotalice) ob % mwscatt_totalice = .true.

      ! Create obs vector and r matrix
      jchans_used = 0
      do jvar = 1, self%nchans
        if( obs % QCflags(jvar,jobs) == self % passflag ) then
          jchans_used = jchans_used + 1
          ob % yobs(jchans_used) = obs % yobs(jvar, jobs)
          ob % channels_used(jchans_used) = self % channels(jvar)
          ob % emiss(jchans_used) = obs % emiss(jvar, jobs)
        end if
      end do
      call r_submatrix % setup(nchans_used, ob % channels_used, full_rmatrix=full_rmatrix)

      ! Setup hofxdiags for this retrieval
      call ufo_geovals_setup(hofxdiags, retrieval_vars, 1)

      if (self % FullDiagnostics) then
        call ob % info()
        call r_submatrix % info()
        write(*, *) "Observations used = ",ob % yobs(:)
        write(*,*) "ob % emiss = ",ob % emiss
        write(*,*) "ob % calc_emiss = ",ob % calc_emiss
        write(*,*) "Channel selection = "
        write(*,'(15I5)') ob % channels_used
        write(*,*) "All Channels = "
        write(*,'(15I5)') ob % channels_all
      end if

      !---------------------------------------------------
      ! 2.3 Call minimization
      !---------------------------------------------------
      if (self % UseMLMinimization) then
        call ufo_rttovonedvarcheck_minimize_ml(self, ob, &
                                      r_submatrix, b_matrix, b_inverse, b_sigma, &
                                      local_geovals, hofxdiags, rttov_simobs, &
                                      prof_index, onedvar_success)
      else
        call ufo_rttovonedvarcheck_minimize_newton(self, ob, &
                                      r_submatrix, b_matrix, b_inverse, b_sigma, &
                                      local_geovals, hofxdiags, rttov_simobs, &
                                      prof_index, onedvar_success)
      end if

      obs % output_BT(:, jobs) = ob % output_BT(:)
      obs % output_profile(:,jobs) = ob % output_profile(:)
      obs % final_cost(jobs) = ob % final_cost

      ! Set QCflags based on output from minimization
      if (.NOT. onedvar_success) then
        do jvar = 1, self%nchans
          if( obs % QCflags(jvar,jobs) == 0 ) then
            obs % QCflags(jvar,jobs) = self % onedvarflag
          end if
        end do
      end if

      ! Tidy up memory specific to a single observation
      call ufo_geovals_delete(local_geovals)
      call ufo_geovals_delete(hofxdiags)
      call ob % delete()
      call r_submatrix % delete()

    else
      call fckit_log % info("Final 1Dvar cost, apply = F")

    endif
  end do obs_loop

  !---------------------------------------------------
  ! 3.0 Return variables and tidy up
  !---------------------------------------------------

  write(message, *) "Total number of observations = ", obs % iloc
  call fckit_log % info(message)
  write(message, *) "Number tested by 1dvar = ", apply_count
  call fckit_log % info(message)

  ! Put qcflags and output variables into observation space
  call obs % output(self % obsdb, prof_index, vars, self % nchans)

  ! Tidy up memory used for all observations
  call full_bmatrix % delete()
  call full_rmatrix % delete()
  call obs % delete()
  if (self % pcemiss) call IR_pcemis % delete()
  if (allocated(b_matrix))  deallocate(b_matrix)
  if (allocated(b_inverse)) deallocate(b_inverse)
  if (allocated(b_sigma))   deallocate(b_sigma)
  call rttov_simobs % delete()

end subroutine ufo_rttovonedvarcheck_apply

end module ufo_rttovonedvarcheck_mod
