! (C) Copyright 2017-2018 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

!> Fortran module for radiancerttov observation operator

module ufo_radiancerttov_mod

  use fckit_configuration_module, only: fckit_configuration
  use fckit_log_module, only : fckit_log
  use iso_c_binding
  use kinds
  use missing_values_mod

  use obsspace_mod

  use ufo_geovals_mod, only: ufo_geovals, ufo_geoval, ufo_geovals_get_var
  use ufo_vars_mod
  use ufo_radiancerttov_utils_mod

  use rttov_types
  use rttov_const, only : errorstatus_fatal, errorstatus_success
  use rttov_unix_env

  implicit none
  private

  !> Fortran derived type for the observation type
  type, public :: ufo_radiancerttov
    private
    character(len=MAXVARLEN), public, allocatable :: varin(:)      ! variables requested from the model.
    integer, allocatable                          :: channels(:)   ! list of instrument channels to simulate.
    integer, allocatable                          :: qc_passed(:)  ! list of values indicating passed qc
    integer, allocatable                          :: coefindex(:)  ! list of the coefindex for the channels to simulate.
    type(rttov_conf)                              :: conf
    type(ufo_rttov_io)                            :: RTProf
  contains
    procedure :: setup  => ufo_radiancerttov_setup
    procedure :: delete => ufo_radiancerttov_delete
    procedure :: simobs => ufo_radiancerttov_simobs
  end type ufo_radiancerttov

contains

  ! ------------------------------------------------------------------------------
  subroutine ufo_radiancerttov_setup(self, f_confOper, channels, qc_passed)

    implicit none
    class(ufo_radiancerttov), intent(inout) :: self
    type(fckit_configuration), intent(in)   :: f_confOper
    integer(c_int),            intent(in)   :: channels(:)  !List of channels to use
    integer(c_int), optional,  intent(in)   :: qc_passed(:) !List of values indicating passed qc

    type(fckit_configuration)               :: f_confOpts ! RTcontrol
    integer                                 :: ind, j, jspec, ii, jj, jnew
    integer                                 :: nvars_in
    character(len=800)                      :: message

    call f_confOper % get_or_die("obs options",f_confOpts)

! Begin RTTOV configuration and determine ngas (Absorbers)
    call rttov_conf_setup(self % conf, f_confOpts, f_confOper)

! Count mandatory inputs and additional gases
    nvars_in = size(varin_default) + self%conf%ngas

! Add inputs required for RTTOV-SCATT
    if (self % conf % do_mw_scatt) then
      do j = 1, size(varin_scatt)
        if (.not. any(varin_scatt(j) == self%conf%Absorbers)) then
          nvars_in = nvars_in + 1
        end if
      end do
    end if

! Allocate and store names of requested inputs in self%varin for geovals request
    allocate(self%varin(nvars_in))
    self%varin(1:size(varin_default)) = varin_default
    ind = size(varin_default) + 1

    !Use list of Absorbers from conf
    do j = 1, self%conf%ngas
      self%varin(ind) = self%conf%Absorbers(j)
      ind = ind + 1
    end do

    ! Number of channels to be simulated for this instrument (from the configuration, not necessarily the full instrument complement)
    self % RTprof % nchan_inst = size(channels)
    allocate(self % channels(self % RTprof % nchan_inst))
    allocate(self % coefindex(self % RTprof % nchan_inst))
    self % coefindex(:) = 0
    self % channels(:) = channels

    ! Copy flags to fortran object
    if (present(qc_passed)) then
      allocate(self % qc_passed(size(qc_passed)))
      self % qc_passed(:) = qc_passed(:)
    end if

    jnew = 1
    coefloop: do ii = 1, self % RTprof % nchan_inst
      jj = jnew
      do while ( jj <= self % conf % rttov_coef_array(1) % coef % fmv_chn )
        if (channels(ii) == self % conf % rttov_coef_array(1) % coef % ff_ori_chn(jj)) then
          self % coefindex(ii) = jj
          cycle coefloop
        end if
        jj = jj + 1
      end do
    end do coefloop

    if ( any(self % coefindex == 0) ) then
      message = 'ufo_radiancerttov_setup error: input channels not in the coefficient file'
      call abor1_ftn(message)
    end if

    !Add RTTOV-SCATT inputs ensuring not to double count duplicates
    if (self % conf % do_mw_scatt) then
      do j = 1, size(varin_scatt)
        if (.not. any(varin_scatt(j) == self%conf%Absorbers)) then
          self%varin(ind) = varin_scatt(j)
          ind = ind + 1
        end if
      enddo
    end if

    if (self % conf % debug) then
      do j=1,size(self%varin)
        write(message,'(I4,1x,A)') j, trim(self%varin(j))
        call fckit_log % debug(message)
      enddo
    end if

    message = 'Finished setting up rttov'
    call fckit_log%info(message)

  end subroutine ufo_radiancerttov_setup

  ! ------------------------------------------------------------------------------
  subroutine ufo_radiancerttov_delete(self)
    implicit none
    class(ufo_radiancerttov), intent(inout) :: self

    call rttov_conf_delete(self%conf)
    if (allocated(self % varin)) deallocate(self % varin)
    if (allocated(self % channels)) deallocate(self % channels)
    if (allocated(self % coefindex)) deallocate(self % coefindex)

  end subroutine ufo_radiancerttov_delete

  ! ------------------------------------------------------------------------------
  subroutine ufo_radiancerttov_simobs(self, geovals, obss, nvars, nlocs,      &
                                      hofx, hofxdiags, qcf_p, ob_info)
    use fckit_mpi_module,   only: fckit_mpi_comm
    use obsdatavector_mod,  only: obsdatavector_int
    use iso_fortran_env,    only: int64
    use ufo_rttovonedvarcheck_ob_mod

    implicit none

    class(ufo_radiancerttov), intent(inout)   :: self
    type(ufo_geovals),        intent(in)      :: geovals
    type(c_ptr), value,       intent(in)      :: obss
    integer,                  intent(in)      :: nvars, nlocs

    real(c_double),        intent(inout)      :: hofx(nvars,nlocs)
    type(ufo_geovals),     intent(inout)      :: hofxdiags    !non-h(x) diagnostics
    type(c_ptr), value, optional, intent(in)  :: qcf_p
    type(ufo_rttovonedvarcheck_ob), optional, intent(inout) :: ob_info

    real(c_double)                          :: missing
    type(fckit_mpi_comm)                    :: f_comm

    ! Local Variables
    character(*), parameter                 :: routine_name = 'ufo_radiancerttov_simobs'
    character(len=800)                      :: message
    type(rttov_chanprof), allocatable       :: chanprof(:)
    type(ufo_geoval), pointer               :: geoval_temp
    type(rttov_hofxdiags)                   :: hofxdiags_methods

    integer                                 :: nprofiles, nlevels
    integer(kind=jpim)                      :: errorstatus  ! Error status of RTTOV subroutine calls

    integer                                 :: iprof_rttov, iprof, ichan, ichan_sim, jchan
    integer                                 :: nprof_sim, nprof_max_sim, nchan_total
    integer                                 :: nchan_sim, qc_flag
    integer(int64)                          :: ii, pp
    integer                                 :: prof_start, prof_end
    integer                                 :: sensor_index
    integer, allocatable                    :: prof_list(:,:)  ! store list of 'good' profiles

    logical                                 :: jacobian_needed, skip_profile
    real(kind_real), allocatable            :: sfc_emiss(:,:)
    type(obsdatavector_int)                 :: qc_flags

    include 'rttov_direct.interface'
    include 'rttov_scatt.interface'
    include 'rttov_k.interface'
    include 'rttov_scatt_ad.interface'

    message = trim(routine_name) // ': Simulating observations'
    call fckit_log%debug(message)

    !Initialisations
    missing = missing_value(missing)
    hofx    = missing

    !Return the name and name length of obsspace communicator (from ioda)
    call obsspace_get_comm(obss, f_comm)

    !! Parse hofxdiags%variables into independent/dependent variables and channel assumed formats:
    ! Note this sets the jacobian_needed flag
    !!   jacobian var -->     <ystr>_jacobian_<xstr>_<chstr>
    !!   non-jacobian var --> <ystr>_<chstr>
    call hofxdiags_methods % reset()
    call hofxdiags_methods % parse(hofxdiags, jacobian_needed)

    ! Get number of profiles and levels from geovals
    nprofiles = geovals % nlocs
    call ufo_geovals_get_var(geovals, var_ts, geoval_temp)
    nlevels = geoval_temp % nval
    nullify(geoval_temp)

    ! Setup qc flags will be run profile by profile
    if (present(qcf_p)) then
      if (self % conf % UseQCFlagsToSkipHofX) qc_flags % data_ptr = qcf_p
    end if

    ! Sanity checks
    if (nprofiles == 0) return

    ! Allocate RTTOV profiles for ALL geovals for the direct calculation
    write(message,'(2A, I0, A, I0, A)') &
      trim(routine_name), ': Allocating ', nprofiles, ' profiles with ', nlevels, ' levels'
    call fckit_log%debug(message)

    call self % RTprof % alloc_profiles(errorstatus, self % conf, nprofiles, nlevels, init=.true., asw=1)

    ! Assign the atmospheric and surface data from the GeoVaLs
    message = trim(routine_name) // ': Creating RTTOV profiles from geovals'
    call fckit_log%debug(message)
    if(present(ob_info)) then
      call self % RTprof % setup_rtprof(geovals,obss,self % conf,ob_info=ob_info)
    else
      call self % RTprof % setup_rtprof(geovals,obss,self % conf)
    end if

    ! Read emissivity from obs space if its requested
    if (self % conf % surface_emissivity_group /= "") then
      allocate(sfc_emiss(self % RTprof % nchan_inst, nprofiles)) ! nchans, nprofiles
      call rttov_read_emissivity_from_obsspace(obss, self % conf % surface_emissivity_group, &
                                               self % channels, sfc_emiss)
    end if

    ! Maximum number of profiles to be processed by RTTOV per pass
    if(self % conf % prof_by_prof) then
      nprof_max_sim = 1
    else
      nprof_max_sim = max(1,self % conf % nchan_max_sim / self % RTprof % nchan_inst)
    end if
    nprof_sim = min(nprof_max_sim, nprofiles)

    ! Determine the total number of radiances to simulate (nchan_sim).
    nchan_sim = nprof_sim * size(self%channels)

    ! Allocate structures for RTTOV direct code (and, if needed, K code and any structures for RTTOV-SCATT)
    write(message,'(2A,I0,A,I0,A)')                                                   &
      trim(routine_name), ': Allocating resources for RTTOV direct code: ', nprof_sim, ' and ', nchan_sim, ' channels'
    call fckit_log%debug(message)
    call self % RTprof % alloc_direct(errorstatus, self % conf, nprof_sim, nchan_sim, nlevels, init=.true., asw=1)

    if (jacobian_needed) then
      write(message,'(2A,I0,A,I0,A)')                                                 &
        trim(routine_name), ': Allocating resources for RTTOV K code: ', nprof_sim, ' and ', nchan_sim, ' channels'
      call fckit_log%debug(message)

      call self % RTprof % alloc_profiles_k(errorstatus, self % conf, nchan_sim, nlevels, init=.true., asw=1)
      call self % RTprof % alloc_k(errorstatus, self % conf, nprof_sim, nchan_sim, nlevels, init=.true., asw=1)
    end if

    ! Used for keeping track of profiles for setting emissivity
    allocate(self % RTprof % chanprof ( nprofiles * self % RTprof % nchan_inst ))

    prof_start = 1; prof_end = nprofiles
    nchan_total = 0

    RTTOV_loop : do while (prof_start <= prof_end)
      ! Check qc flags if requested and skip if no active channels
      ! note prof_by_prof is set to true if UseQCFlagsToSkipHofX is true
      if (present(qcf_p) .and. allocated(self % qc_passed)) then
        if (self % conf % UseQCFlagsToSkipHofX) then
          skip_profile = .true.
          do ichan = 1, self % RTprof % nchan_inst
            ii = ichan ! conversion to int64
            pp = prof_start ! conversion to int64
            qc_flag = qc_flags % get(ii, pp)
            if (any(self % qc_passed == qc_flag)) then
              skip_profile = .false.
              exit
            end if
          end do
          if (skip_profile) then
            ! increment profile and channel counters and skip this one
            nchan_total = nchan_total + nchan_sim
            prof_start = prof_start + nprof_sim
            cycle
          end if
        end if
      end if

      ! Zero all k code variables.  These arrays are of size prof_end-prof_start
      if (jacobian_needed) call self % RTprof % zero_k(self % conf, reset_profiles_k=.false.)

      ! Reduce number of simulated profiles/channel if at end of the of profiles to be processed
      nprof_sim = min(nprof_sim, prof_end - prof_start + 1)
      nchan_sim = nprof_sim * size(self%channels)

      ! allocate and initialise local chanprof structure
      allocate(chanprof(nchan_sim))
      chanprof(:) % prof = 0
      chanprof(:) % chan = 0

      ! index for simulated channel
      ichan_sim = 0_jpim
      nchan_sim = 0_jpim

      !allocate list used to store 'good' profiles
      !initialise to -1, so no bad profile is given an emissivity
      allocate(prof_list(nprof_sim,2))
      prof_list = -1

      ! Build the list of profile/channel indices in chanprof
      do iprof_rttov = 1, nprof_sim
        errorstatus = errorstatus_success

        ! iprof is the index for the full set of RTTOV profiles
        iprof = prof_start + iprof_rttov - 1

        ! print profile information if requested
        if(any(self % conf % inspect == iprof)) call self % RTprof % print_rtprof(self % conf, iprof)

        ! check RTTOV profile will be valid for RTTOV and flag it if it fails the check
        call self % RTprof % check_rtprof(self % conf, iprof, errorstatus)

        ! check sfc_emiss valid if read in
        if(errorstatus == errorstatus_success) then
          if (allocated(sfc_emiss)) then
            do ichan = 1, self % RTprof % nchan_inst
              if ((sfc_emiss(ichan,iprof) > 1.0) .or. (sfc_emiss(ichan,iprof) < 0.0)) then
                errorstatus = errorstatus_fatal
              end if
            end do
          end if

          prof_list(iprof_rttov,1) = iprof_rttov ! chunk index
          prof_list(iprof_rttov,2) = iprof       ! all-obs index
          do ichan = 1, self % RTprof % nchan_inst
            ichan_sim = ichan_sim + 1_jpim
            chanprof(ichan_sim) % prof = iprof_rttov ! this refers to the slice of the RTprofile array passed to RTTOV
            chanprof(ichan_sim) % chan = self % coefindex(ichan)
            self % RTprof % chanprof(nchan_total + ichan_sim) % prof = iprof ! this refers to the index of the profile from the geoval
            self % RTprof % chanprof(nchan_total + ichan_sim) % chan = self % coefindex(ichan)
          end do
          nchan_sim = ichan_sim

          ! Pick the last valid (unskipped) observation to get the sensor index
          ! but it will always be the same as the first because prof_by_prof is true for nsensors > 1
          ! This is to guard against a bad satellite identifier sneaking in when only processing one instrument.
          sensor_index = self % RTProf % sensor_index_array(iprof)
        end if

      end do ! loop over profiles in chunk

      ! If there has been an rttov check error then nchan_sim will equal zero for
      ! the 1DVar and the flag needs setting.
      if (present(ob_info)) then
        if (nchan_sim == 0) then
          ob_info % rterror = .true.
        end if
      end if

      ! Set surface emissivity
      if (present(ob_info)) then
        if (.not. ob_info % rterror) then
          self % RTprof % calcemis(1:nchan_sim) = ob_info % calc_emiss(:)
          self % RTprof % emissivity(1:nchan_sim) % emis_in = ob_info % emiss(:)
        end if
      else
        if (allocated(sfc_emiss)) then
          self % RTprof % calcemis(:) = .false.
          outerloop: do ichan = 1, ichan_sim  ! list of channels*profiles
            do jchan = 1, self % RTprof % nchan_inst  ! list of self % channels
              ! if the channel number for this channel * profile == channel number needed
              ! chanprof(ichan) % chan refers to the index in the coefficient file
              if (self % conf % rttov_coef_array(1) % coef % ff_ori_chn(chanprof(ichan) % chan) == self % channels(jchan)) then
                iprof = prof_start + chanprof(ichan) % prof - 1
                self % RTprof % emissivity(ichan) % emis_in = sfc_emiss(jchan, iprof)
                self % RTprof % calcemis(ichan) = .false.
                if (self % RTprof % emissivity(ichan) % emis_in == 0.0) then
                  self % RTprof % calcemis(ichan) = .true.
                end if
                cycle outerloop
              end if
            end do
          end do outerloop
        else
          call self % RTProf % init_default_emissivity(self % conf, prof_list)
        end if
      end if

      ! Write out emissivity if checking profile
      if (size(self % conf % inspect) > 0) then
        do ichan = 1, ichan_sim, self % RTprof % nchan_inst
          iprof = prof_start + chanprof(ichan) % prof - 1
          if(any(self % conf % inspect == iprof)) then
            write(*,*) "profile ", iprof
            write(*,*) "calcemiss = ",self % RTprof % calcemis(ichan:ichan+self%RTprof%nchan_inst-1)
            write(*,*) "emissivity in = ",self % RTprof % emissivity(ichan:ichan+self%RTprof%nchan_inst-1) % emis_in
          end if
        end do
      end if

      deallocate(prof_list)

      ! --------------------------------------------------------------------------
      ! Call RTTOV model
      ! --------------------------------------------------------------------------
      !N.B. different from TL/AD as we don't need to retain the full profiles_k so data can be overwritten
      if (nchan_sim > 0) then
        if (jacobian_needed) then
          if (self % conf % do_mw_scatt) then
            call rttov_scatt_ad(                                                         &
              errorstatus,                                                               &! out   error flag
              self % conf % mw_scatt % opts,                                             &! in    options structure
              int (nlevels, KIND = jpim),                                                &
              chanprof(1:nchan_sim),                                                     &! in    LOCAL channel and profile index structure
              self % RTprof % mw_scatt % freq_indices(1:nchan_sim),                      &! in    frequency indices
              self % RTProf % profiles(prof_start:prof_start + nprof_sim -1),            &! in    profile array
              self % RTProf % mw_scatt % profiles(prof_start:prof_start + nprof_sim -1), &! in    scattering profile array
              self % conf % rttov_coef_array(sensor_index),                                &! in    coefficients structure
              self % conf % mw_scatt % coef,                                             &! in    scatt coefficients structure
              self % RTProf % calcemis(1:nchan_sim),                                     &! in    flag for internal emissivity calcs
              self % RTProf % emissivity(1:nchan_sim),                                   &! inout input/output emissivities per channel
              self % RTProf % profiles_k(1:nchan_sim),                                   &! inout 
              self % RTProf % mw_scatt % profiles_k(1:nchan_sim),                        &! inout 
              self % RTProf % emissivity_k(1:nchan_sim),                                 &! inout input/output emissivity jacs per channel
              self % RTProf % radiance,                                                  &! inout computed radiances
              self % RTProf % radiance_k)                                                 ! inout computed radiance jacobians
          else
            call rttov_k(                                                     &
              errorstatus,                                                    &! out   error flag
              chanprof(1:nchan_sim),                                          &! in    LOCAL channel and profile index structure
              self % conf % rttov_opts,                                       &! in    options structure
              self % RTProf % profiles(prof_start:prof_start + nprof_sim -1), &! in    profile array
              self % RTProf % profiles_k(1:nchan_sim),                        &! inout profile array
              self % conf % rttov_coef_array(sensor_index),                   &! in    coefficients structure
              self % RTProf % transmission,                                   &! inout computed transmittances
              self % RTProf % transmission_k,                                 &! inout computed transmittances
              self % RTProf % radiance,                                       &! inout computed radiances
              self % RTProf % radiance_k,                                     &! inout computed radiances
              calcemis    = self % RTProf % calcemis(1:nchan_sim),            &! in    flag for internal emissivity calcs
              emissivity  = self % RTProf % emissivity(1:nchan_sim),          &! inout input/output emissivities per channel
              emissivity_k = self % RTProf % emissivity_k(1:nchan_sim))        ! inout input/output emissivities per channel
          end if
          
          if ( errorstatus /= errorstatus_success ) then
            write(message,'(2A, I6, A, I6, A, I6)') trim(routine_name), 'after rttov_k: error ', errorstatus, &
              ' skipping profiles ', prof_start, ' -- ', prof_start + nprof_sim - 1
            call fckit_log%info(message)
          end if
        else ! direct
          if (self % conf % do_mw_scatt) then
            call rttov_scatt (                                                           &
              errorstatus,                                                               &! out   error flag
              self % conf % mw_scatt % opts,                                             &! in    options structure
              int (nlevels, KIND = jpim),                                                &
              chanprof(1:nchan_sim),                                                     &! in    channel and profile index structure
              self % RTprof % mw_scatt % freq_indices(1:nchan_sim),                      &
              self % RTProf % profiles(prof_start:prof_start + nprof_sim -1),            &! in    profile array
              self % RTProf % mw_scatt % profiles(prof_start:prof_start + nprof_sim -1), &
              self % conf % rttov_coef_array(sensor_index),                                &! in    coefficients structure
              self % conf % mw_scatt % coef,                                             &! in    scatt coefficients structure
              self % RTProf % calcemis(1:nchan_sim),                                     &! in    flag for internal emissivity calcs
              self % RTProf % emissivity(1:nchan_sim),                                   &! inout input/output emissivities per channel
              self % RTProf % radiance,                                                  &! inout computed radiances
              emis_retrieval_terms = self % RTprof % mw_scatt % emis_retrieval)           !             
          else        
            call rttov_direct(                                                &
              errorstatus,                                                    &! out   error flag
              chanprof(1:nchan_sim),                                          &! in    channel and profile index structure
              self % conf % rttov_opts,                                       &! in    options structure
              self % RTProf % profiles(prof_start:prof_start + nprof_sim -1), &! in    profile array
              self % conf % rttov_coef_array(sensor_index),                     &! in    coefficients structure
              self % RTProf % transmission,                                   &! inout computed transmittances
              self % RTProf % radiance,                                       &! inout computed radiances
              calcemis    = self % RTProf % calcemis(1:nchan_sim),            &! in    flag for internal emissivity calcs
              emissivity  = self % RTProf % emissivity(1:nchan_sim))!,        &! inout input/output emissivities per channel          
          end if
          
          if ( errorstatus /= errorstatus_success ) then
            write(message,'(2A, I6, A, I6, A, I6)') trim(routine_name), 'after rttov_direct: error ', errorstatus, &
                                         ' skipping profiles ', prof_start, ' -- ', prof_start + nprof_sim - 1
            call fckit_log%info(message)
            if (present(ob_info)) ob_info % rterror = .true.
          end if

        end if

        ! Put simulated brightness temperature into hofx
        if ( errorstatus == errorstatus_success ) then
          do ichan = 1, nchan_sim, size(self%channels)
            iprof = self % RTProf % chanprof(nchan_total + ichan)%prof
            hofx(1:size(self%channels),iprof) = self % RTprof % radiance % bt(ichan:ichan+size(self%channels)-1)
            
            !store transmittance if ob_info present in call and transmittance part of structure
            if (present(ob_info)) then
              if (allocated(ob_info % transmittance)) then
                if (self % conf % do_mw_scatt) then
                  ob_info % transmittance(1:size(self%channels)) = self % RTprof % mw_scatt % emis_retrieval % tau_clr(ichan:ichan+size(self%channels)-1)
                else
                  ob_info % transmittance(1:size(self%channels)) = self % RTProf % transmission % tau_total(ichan:ichan+size(self%channels)-1)
                end if
              end if
            end if
          enddo

          ! Write out emissivity out and hofx
          if (size(self % conf % inspect) > 0) then
            do ichan = 1, ichan_sim, self % RTprof % nchan_inst
              iprof = prof_start + chanprof(ichan) % prof - 1
              if(any(self % conf % inspect == iprof)) then
                write(*,*) "profile ", iprof
                write(*,*) "emissivity out = ",self % RTprof % emissivity(ichan:ichan+self%RTprof%nchan_inst-1) % emis_out
                write(*,*) "hofx out = ",hofx(1:size(self%channels),iprof)
                
                if (self % conf % do_mw_scatt) then
                  write(*,*) "emis_retrieval cfrac = ",self % RTprof % mw_scatt % emis_retrieval % cfrac(ichan:ichan+self%RTprof%nchan_inst-1)
                  write(*,*) "emis_retrieval tau_clr = ",self % RTprof % mw_scatt % emis_retrieval % tau_clr(ichan:ichan+self%RTprof%nchan_inst-1)
                  write(*,*) "emis_retrieval tau_cld = ",self % RTprof % mw_scatt % emis_retrieval % tau_cld(ichan:ichan+self%RTprof%nchan_inst-1)
                else
                  write(*,*) "tau_total out = ",self % RTProf % transmission % tau_total(ichan:ichan+self%RTprof%nchan_inst-1)
                end if
              end if
            end do
          end if

          ! Put simulated diagnostics into hofxdiags
          if(hofxdiags % nvar > 0) &
            call hofxdiags_methods % populate(self % RTProf, chanprof, self % conf, prof_start, hofxdiags)
        end if
      end if ! nchan_sim > 0

      ! increment profile and channel counters
      nchan_total = nchan_total + nchan_sim
      prof_start = prof_start + nprof_sim

      ! deallocate local chanprof so it can be re-allocated with a different number of channels if reqd.
      deallocate(chanprof)
      if (jacobian_needed) call self % RTprof % zero_k(self % conf, reset_profiles_k=.true.)
      
    end do RTTOV_loop

    message = 'Deallocating resource for RTTOV...'
    call fckit_log%debug(message)

    ! Deallocate structures for rttov_direct
    if(jacobian_needed) then
      call self % RTprof % alloc_k(errorstatus, self % conf, -1, -1, -1, asw=0)
      call self % RTprof % alloc_profiles_k(errorstatus, self % conf, size(self % RTprof % profiles_k), -1, asw=0)
      ! deallocation of profiles_k isn't done by default in alloc_profs_k because it can contain the trajectory 
      ! which is currently used for the TL/AD but the 1dvar doesn't use it so it can be safely done here
      deallocate (self % RTprof % profiles_k)
      if (self % conf % do_mw_scatt) deallocate(self % RTProf % mw_scatt % profiles_k)
    endif
    call self % RTprof % alloc_direct(errorstatus, self % conf, -1, -1, -1, asw=0)
    call self % RTprof % alloc_profiles(errorstatus, self % conf, size(self % RTprof % profiles), -1, asw=0)

    deallocate(self % RTprof % chanprof)
    deallocate(self % RTprof % sensor_index_array)
    
    if (errorstatus /= errorstatus_success) then
      message = 'after rttov_alloc_direct (deallocation)'
      call abor1_ftn(message)
    else
      message = 'Done. Returning'
      call fckit_log%debug(message)
    end if

  end subroutine ufo_radiancerttov_simobs

  ! ------------------------------------------------------------------------------

end module ufo_radiancerttov_mod
