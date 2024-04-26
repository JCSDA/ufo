! (C) Copyright 2017-2020 Met Office
! 
! this software is licensed under the terms of the apache licence version 2.0
! which can be obtained at http://www.apache.org/licenses/license-2.0.

!> The main Fortran module for implementing the rttov onedvar check

module ufo_rttovonedvarcheck_mod

use ufo_constants_mod, only: zero
use fckit_configuration_module, only: fckit_configuration
use fckit_log_module, only : fckit_log
use iso_c_binding
use kinds
use obsspace_mod
use oops_variables_mod
use ufo_geovals_mod
use ufo_metoffice_bmatrixstatic_mod
use ufo_metoffice_rmatrixradiance_mod
use ufo_radiancerttov_mod
use ufo_rttovonedvarcheck_constants_mod
use ufo_rttovonedvarcheck_minimize_utils_mod
use ufo_rttovonedvarcheck_minimize_newton_mod
use ufo_rttovonedvarcheck_minimize_ml_mod
use ufo_rttovonedvarcheck_ob_mod
use ufo_rttovonedvarcheck_obs_mod
use ufo_rttovonedvarcheck_profindex_mod
use ufo_rttovonedvarcheck_rsubmatrix_mod
use ufo_rttovonedvarcheck_setup_mod
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
  self % onedvarflag = onedvarflag
  self % passflag = passflag

  call ufo_rttovonedvarcheck_setup(self, f_conf, channels) ! from init

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
subroutine ufo_rttovonedvarcheck_apply(self, f_conf, vars, hofxdiags_vars, geovals, apply)
  use ufo_utils_mod, only: cmp_strings

  implicit none
  type(ufo_rttovonedvarcheck), intent(inout) :: self     !< rttovonedvarcheck main object
  type(fckit_configuration), intent(in)      :: f_conf       !< yaml file contents
  type(oops_variables), intent(in)           :: vars     !< channels for 1D-Var
  type(oops_variables), intent(in)           :: hofxdiags_vars !< retrieval variables for 1D-Var
  type(ufo_geovals), intent(in)              :: geovals  !< model values at observation space
  logical, intent(in)                        :: apply(:) !< qc manager flags

  type(ufo_rttovonedvarcheck_obs)        :: obs             ! data for all observations read from db
  type(ufo_metoffice_bmatrixstatic)      :: full_bmatrix    ! full bmatrix read from file
  type(ufo_geovals)                      :: firstguess_geovals ! geoval for one observation
  type(ufo_rttovonedvarcheck_ob)         :: ob              ! observation data for a single observation
  type(ufo_rttovonedvarcheck_profindex)  :: prof_index      ! index for mapping geovals to 1d-var state profile
  type(ufo_rttovonedvarcheck_profindex)  :: local_profindex ! local copy of prof_index needed since the intro. of mwemiss
  type(ufo_metoffice_rmatrixradiance)    :: full_rmatrix    ! full r_matrix read from file
  type(ufo_rttovonedvarcheck_rsubmatrix) :: r_submatrix     ! r_submatrix object
  type(ufo_geovals)                      :: hofxdiags       ! hofxdiags containing jacobian
  type(ufo_geoval), pointer          :: geoval
  character(len=max_string)          :: var
  character(len=max_string)          :: varname
  character(len=max_string)          :: message
  integer                            :: jvar, ivar, jobs, band, ii, irej, jnew ! counters
  integer                            :: igval, gv_index ! counters for geoval
  integer                            :: nchans_used      ! counter for number of channels used for an ob
  integer                            :: jchans_used
  integer                            :: fileunit        ! unit number for reading in files
  integer                            :: apply_count ! number of profiles that the 1dvar has been applied to
  integer                            :: failed_1dvar_count ! number of profiles that failed to converge
  integer                            :: failed_retrievedBTcheck_count ! number of profiles with retrieved BTs outside error
  integer                            :: nprofelements   ! number of elements in 1d-var state profile
  integer, allocatable               :: fields_in(:)
  real(kind_real)                    :: t1, t2          ! timing
  real(kind_real), allocatable       :: b_matrix(:,:)   ! 1d-var profile b matrix
  real(kind_real), allocatable       :: b_inverse(:,:)  ! inverse for each 1d-var profile b matrix
  real(kind_real), allocatable       :: b_sigma(:)      ! b_matrix diagonal error
  real(kind_real), allocatable       :: max_error(:)    ! max_error = error(stdev) * factor
  logical                            :: file_exists     ! check if a file exists logical
  logical                            :: onedvar_success
  logical                            :: reject_profile
  type(ufo_radiancerttov)            :: rttov_simobs
  integer(c_size_t), allocatable     :: ret_nlevs(:)
  integer(c_size_t)                  :: npaths_by_method(1)
  integer(c_size_t), allocatable     :: sampling_method_by_var(:)
  type(oops_variables)               :: reduced_vars
  integer(c_size_t)                  :: nreduced_vals(0)
  logical(c_bool), parameter         :: is_sampling_method_trivial(1) = (/ .true. /)

  ! ------------------------------------------
  ! 1. Setup
  ! ------------------------------------------
  ! Setup rttov simobs
  call rttov_simobs % setup(f_conf, self % channels)

  ! Setup full B matrix object
  call full_bmatrix % setup(self % retrieval_variables, self % b_matrix_path, &
                            self % qtotal)
  
  ! Setup full R matrix object
  call full_rmatrix % setup(self % r_matrix_path)

  ! Create profile index for mapping 1d-var profile to b-matrix
  call prof_index % setup(full_bmatrix, self % nlevels)
  if (prof_index % mwemiss(1) > 0 .and. .not. self % mwEmissRetrieval) then
    write(message, *) "MWemiss is in the b-matrix but the mw retrieval flag, ", &
                      "retrieve mw emissivity, is false.  Set this to true in ", &
                      "the surface emissivty configuration to allow the correct arrays ", &
                      "to be loaded => aborting."
    call abor1_ftn(message)
  end if

  ! sanity check on clw storage - only store to obs space
  ! if qtotal or ql is present in retrieval vector
  if (self % Store1DVarCLW) then
    if ((prof_index % qt(1) == 0) .and. (prof_index % ql(1) == 0)) then
      write(message,*) "Info: as qtotal or ql are not part of the state vector resetting Store1DVarCLW to false"
      self % Store1DVarCLW = .false.
    end if
  end if
  
  ! Read in observation data from obsspace
  ! geovals are only providing surface information
  call obs % setup(self, prof_index, geovals, vars)

  ! Initialize data arrays
  allocate(b_matrix(prof_index % nprofelements,prof_index % nprofelements))
  allocate(b_inverse(prof_index % nprofelements,prof_index % nprofelements))
  allocate(b_sigma(prof_index % nprofelements))
  allocate(ret_nlevs(hofxdiags_vars % nvars()))

  ! Decide on loop parameters - testing
  if (self % StartOb == 0) self % StartOb = 1
  if (self % FinishOb == 0) self % FinishOb = obs % iloc
  if (self % StartOb > self % FinishOb) then
    write(message,*) "start loop ",self % StartOb," is greater than finish loop ",self % FinishOb
    call abor1_ftn(message)
  end if

  ! Calculate hofxdiags levels for each variable
  call ufo_rttovonedvarcheck_hofxdiags_levels(hofxdiags_vars, self % nlevels, ret_nlevs)

  ! Define a set of interpolation paths containing just a single path
  npaths_by_method(1) = 1
  allocate(sampling_method_by_var(hofxdiags_vars%nvars()))
  sampling_method_by_var(:) = 1

  ! ------------------------------------------
  ! 2. Beginning main observation loop
  ! ------------------------------------------
  write(*,*) "Beginning loop over observations: ",trim(self%qcname)
  apply_count = 0
  failed_1dvar_count = 0
  failed_retrievedBTcheck_count = 0
  obs_loop: do jobs = self % StartOb, self % FinishOb
    if (apply(jobs)) then

      obs % output_to_db(jobs) = .true.
      apply_count = apply_count + 1
      write(message, *) "starting obs number    ",jobs
      call fckit_log % debug(message)

      ! Needs copying for each ob since mwemiss introduced
      call local_profindex % copy(prof_index)
      !---------------------------------------------------
      ! 2.1 Setup Jb terms
      !---------------------------------------------------
      ! create one ob geovals from full obs geovals and check
      ! to make sure the values are within sensible bounds.
      call ufo_geovals_copy_one(firstguess_geovals, geovals, jobs)
      call ufo_rttovonedvarcheck_check_geovals(self, firstguess_geovals, &
              local_profindex, obs % surface_type(jobs))

      ! create b matrix arrays for this single observation location
      call full_bmatrix % reset( obs % lat(jobs), & ! in
                    b_matrix, b_inverse, b_sigma  ) ! out

      ! adjust b matrix based on tskin error and if mwemiss is in local_profindex
      ! check that we are over land
      call ufo_rttovonedvarcheck_adjust_bmatrix(local_profindex, & ! inout
           obs, jobs, self,                                      & ! in
           b_matrix, b_inverse, b_sigma)                           ! inout

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
        call fckit_log % debug(message)
        cycle obs_loop
      end if

      ! setup ob data for this observation
      call ob % setup(nchans_used, self %  nlevels, local_profindex % nprofelements, self % nchans, &
           self % Store1DVarCLW, self % Store1DVarTransmittance)
      
      ob % forward_mod_name = self % forward_mod_name
      ob % latitude = obs % lat(jobs)
      ob % longitude = obs % lon(jobs)
      ob % date = obs % date(jobs)
      ob % elevation = obs % elevation(jobs)
      ob % sensor_zenith_angle = obs % sat_zen(jobs)
      ob % sensor_azimuth_angle = obs % sat_azi(jobs)
      ob % solar_zenith_angle = obs % sol_zen(jobs)
      ob % solar_azimuth_angle = obs % sol_azi(jobs)
      ob % channels_all = self % channels
      ob % surface_type = obs % surface_type(jobs)
      ob % satellite_identifier = obs % satellite_identifier(jobs)
      ob % calc_emiss = obs % calc_emiss(jobs)
      ob % emiss(:) = obs % emiss(:, jobs)
      if(self % cloud_retrieval) ob % retrievecloud = .true.
      if(self % cloud_retrieval) ob % cloudtopp = obs % cloudtopp(jobs)
      if(self % cloud_retrieval) ob % cloudfrac = obs % cloudfrac(jobs)
      if(self % RTTOV_mwscattSwitch) ob % mwscatt = .true.
      if(self % RTTOV_usetotalice) ob % mwscatt_totalice = .true.
      if(associated(obs % pcemiss_object)) then
        ob % pcemiss_object => obs % pcemiss_object
        allocate(ob % pcemiss(size(obs % pcemiss, 1)))
        ob % pcemiss(:) = obs % pcemiss(:, jobs)
      end if

      ! Use the skinTemperature from the obsspace if it was provided
      if (self % skinTemperatureFromObsSpace) then
        gv_index = 0
        do igval = 1, geovals % nvar
          if (cmp_strings(var_sfc_tskin, firstguess_geovals % variables(igval))) gv_index = igval
        end do
        firstguess_geovals % geovals(gv_index) % vals(1,1) = obs % skinTemperature(jobs)
      end if

      ! Check if ctp very close to model pressure level.  If so
      ! make them exactly equal to match OPS behaviour for RTTOV
      ! jacobian calculation.
      if(self % cloud_retrieval) then
        call ufo_rttovonedvarcheck_check_ctp(ob % cloudtopp, firstguess_geovals, self %  nlevels)
      end if

      ! Store background T in ob data space
      call ufo_geovals_get_var(firstguess_geovals, var_ts, geoval)
      ob % background_T(:) = geoval%vals(:, 1) ! K

      ! Create ob vector and r matrix
      jchans_used = 0
      do jvar = 1, self%nchans
        if( obs % QCflags(jvar,jobs) == self % passflag ) then
          jchans_used = jchans_used + 1
          ob % yobs(jchans_used) = obs % yobs(jvar, jobs)
          ob % channels_used(jchans_used) = self % channels(jvar)
        end if
      end do
      call r_submatrix % setup(nchans_used, ob % channels_used, full_rmatrix=full_rmatrix)

      ! Setup hofxdiags for this retrieval
      reduced_vars = oops_variables()
      call ufo_geovals_setup(hofxdiags, 1, hofxdiags_vars, hofxdiags_vars % nvars(), ret_nlevs, &
                             size(npaths_by_method), npaths_by_method, sampling_method_by_var, &
                             reduced_vars, reduced_vars % nvars(), nreduced_vals, &
                             is_sampling_method_trivial)
      call reduced_vars % destruct()
      call ufo_geovals_setup_trivial_sampling_method(hofxdiags, sampling_method = 1)

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
        call local_profindex % info()
      end if

      !---------------------------------------------------
      ! 2.3 Call minimization
      !---------------------------------------------------
      if (self % UseMLMinimization) then
        call ufo_rttovonedvarcheck_minimize_ml(self, ob, &
                                      r_submatrix, b_matrix, b_inverse, b_sigma, &
                                      firstguess_geovals, hofxdiags, rttov_simobs, &
                                      local_profindex, onedvar_success)
      else
        call ufo_rttovonedvarcheck_minimize_newton(self, ob, &
                                      r_submatrix, b_matrix, b_inverse, b_sigma, &
                                      firstguess_geovals, hofxdiags, rttov_simobs, &
                                      local_profindex, onedvar_success)
      end if

      obs % output_BT(:, jobs) = ob % output_BT(:)
      obs % background_BT(:, jobs) = ob % background_BT(:)
      obs % output_profile(:,jobs) = ob % output_profile(:)
      obs % emiss(:, jobs) = ob % emiss(:)
      obs % final_cost(jobs) = ob % final_cost
      obs % LWP(jobs) = ob % LWP
      obs % IWP(jobs) = ob % IWP
      if (self % store1dvarclw) obs % CLW(:,jobs) = ob % CLW(:)
      if (self % store1dvartransmittance) obs % transmittance(:, jobs) = ob % transmittance(:)
      if (self % cloud_retrieval) obs % cloudtopp(jobs) = ob % cloudtopp
      if (self % cloud_retrieval) obs % cloudfrac(jobs) = ob % cloudfrac
      if (self % cloud_retrieval) obs % cloudtopp_error(jobs) = ob % cloudtopp_error
      if (self % RecalculateBT) obs % recalc_BT(:, jobs) = ob % recalc_BT(:)
      obs % niter(jobs) = ob % niter

      ! Set QCflags based on output from minimization
      if (.NOT. onedvar_success) then
        failed_1dvar_count = failed_1dvar_count + 1
        do jvar = 1, self%nchans
          if( obs % QCflags(jvar,jobs) == 0 ) then
            obs % QCflags(jvar,jobs) = self % onedvarflag
          end if
        end do
      end if

      ! Remove channels that have been removed because of slow convergence
      if (ob % QC_SlowConvChans) then
        do jvar = 1, self % nchans
          if( obs % QCflags(jvar,jobs) == 0 .and. &
              any( self % ConvergeCheckChans == obs % channels(jvar) ) ) then
            obs % QCflags(jvar,jobs) = self % onedvarflag
          end if
        end do
      end if

      ! Reject channels that have failed the ctp check
      if (allocated(ob % rejected_channels_ctp)) then
        jnew = 1
        rejected: do irej = 1, size(ob % rejected_channels_ctp)
          jvar = jnew
          do while ( jvar <= size(obs % channels) )
            if (ob % rejected_channels_ctp(irej) == obs % channels(jvar)) then
              obs % QCflags(jvar, jobs) = self % onedvarflag
              cycle rejected
            end if
            jvar = jvar + 1
          end do
        end do rejected
      end if

      ! Check the BTs are within a factor of the error.  This only applies to channels that are still
      ! active.
      if (self % RetrievedErrorFactor > zero .and. any(obs % QCflags(:,jobs) == self % passflag)) then
        allocate(max_error(size(ob % channels_used)))
        call r_submatrix % multiply_factor_by_stdev(self % RetrievedErrorFactor, max_error)
        reject_profile = .false.
        chanloop: do jvar = 1, size(ob % channels_used)
          if (allocated(ob % rejected_channels_ctp)) then
            if(any(ob % rejected_channels_ctp == ob % channels_used(jvar))) cycle chanloop
          end if
          if (abs(ob % final_bt_diff(jvar)) > max_error(jvar)) reject_profile = .true.
        end do chanloop
        if (reject_profile) then
          failed_retrievedBTcheck_count = failed_retrievedBTcheck_count + 1
          do jvar = 1, self % nchans
            if( obs % QCflags(jvar,jobs) == self % passflag ) then
              obs % QCflags(jvar,jobs) = self % onedvarflag
            end if
          end do
        end if
        deallocate(max_error)
      end if

      ! Tidy up memory specific to a single observation
      call ufo_geovals_delete(firstguess_geovals)
      call ufo_geovals_delete(hofxdiags)
      call ob % delete()
      call r_submatrix % delete()

    else
      call fckit_log % debug("Final 1Dvar cost, apply = F")

    endif
  end do obs_loop

  !---------------------------------------------------
  ! 3.0 Return variables and tidy up
  !---------------------------------------------------

  write(message, *) "Total number of observations = ", obs % iloc
  call fckit_log % info(message)
  write(message, *) "Number tested by 1dvar = ", apply_count
  call fckit_log % info(message)
  write(message, *) "Number that failed to converge = ", failed_1dvar_count
  call fckit_log % info(message)
  write(message, *) "Number that failed the retrieved BT check = ", failed_retrievedBTcheck_count
  call fckit_log % info(message)

  ! Put qcflags and output variables into observation space
  call obs % output(self % obsdb, prof_index, vars, self % nchans)

  ! Tidy up memory used for all observations
  call full_bmatrix % delete()
  call full_rmatrix % delete()
  call obs % delete()
  if (allocated(b_matrix))  deallocate(b_matrix)
  if (allocated(b_inverse)) deallocate(b_inverse)
  if (allocated(b_sigma))   deallocate(b_sigma)
  if (allocated(ret_nlevs)) deallocate(ret_nlevs)
  call rttov_simobs % delete()

end subroutine ufo_rttovonedvarcheck_apply

end module ufo_rttovonedvarcheck_mod
