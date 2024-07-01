! (C) Copyright 2017-2018 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

!> Fortran module for radiancerttov tl/ad observation operator

module ufo_radiancerttov_tlad_mod

  use fckit_configuration_module, only: fckit_configuration
  use fckit_log_module, only : fckit_log
  use iso_c_binding
  use kinds
  use missing_values_mod

  use obsspace_mod

  use ufo_constants_mod, only : zero, g_to_kg
  use ufo_geovals_mod, only: ufo_geovals, ufo_geoval, ufo_geovals_get_var, ufo_geovals_print
  use ufo_vars_mod
  use ufo_radiancerttov_utils_mod

  use rttov_types
  use rttov_const ! errorstatus and gas_id
  use rttov_unix_env

  implicit none
  private

  !> Fortran derived type for radiancerttov trajectory
  type, public :: ufo_radiancerttov_tlad
    private
    character(len=MAXVARLEN), public, allocatable :: varin(:)      ! variables which will be part of the analysis.
    integer, allocatable                          :: channels(:)
    integer, allocatable                          :: coefindex(:)  ! list of the coefindex for the channels to simulate.
    type(rttov_conf)                              :: conf
    type(rttov_conf)                              :: conf_traj
    type(ufo_rttov_io)                            :: RTProf_K

    integer                                       :: nprofiles
    integer                                       :: nchan_total
    integer                                       :: nlevels

    logical                                       :: ltraj

  contains
    procedure :: setup  => ufo_radiancerttov_tlad_setup
    procedure :: delete  => ufo_radiancerttov_tlad_delete
    procedure :: settraj => ufo_radiancerttov_tlad_settraj
    procedure :: simobs_tl  => ufo_radiancerttov_simobs_tl
    procedure :: simobs_ad  => ufo_radiancerttov_simobs_ad
  end type ufo_radiancerttov_tlad

contains

  ! ------------------------------------------------------------------------------
  subroutine ufo_radiancerttov_tlad_setup(self, f_confOper, channels)
    implicit none

    class(ufo_radiancerttov_tlad), intent(inout) :: self
    type(fckit_configuration), intent(in)        :: f_confOper
    integer(c_int),            intent(in)        :: channels(:)  !List of channels to use

    type(fckit_configuration)                    :: f_confOpts ! RTcontrol
    type(fckit_configuration)                    :: f_confLinOper
    integer                                      :: ind, jspec, ii, jj, jnew, nvars
    character(len=:), allocatable                :: str_array(:)
    character(len=800)                           :: message

    call f_confOper % get_or_die("obs options", f_confOpts)
    call f_confOper % get_or_die("linear obs operator", f_confLinOper)

    ! Last argument is false because its setting up the forward model configuration
    ! This is not used at the moment so I have removed the setup.
    ! call rttov_conf_setup(self % conf_traj, f_confOpts, f_confOper, .false.)

    call rttov_conf_setup(self%conf, f_confOpts, f_confOper)

    !DAR what is the RTTOV equivalant of making sure that humidity and ozone data are present
    if ( ufo_vars_getindex(self%conf%Absorbers, var_mixr) < 1 .and. &
      ufo_vars_getindex(self%conf%Absorbers, var_q)    < 1 ) then
      message = 'ufo_radiancerttov_setup error: H2O must be included in RTTOV Absorbers'
      call abor1_ftn(message)
    end if

    ! Read in list of incremented variables from the yaml
    nvars = f_confLinOper % get_size("increment variables")
    allocate(self % varin(nvars))
    call f_confLinOper % get_or_die("increment variables", str_array)
    self % varin(1:nvars) = str_array

    ! channels contains a list of instrument channels. From this we need to work out the coefindex
    ! which RTTOV needs to index the entry in the coefficient file.  This allows cut down
    ! coefficient files to be used.

    ! Number of channels to be simulated for this instrument (from the configuration, not necessarily the full instrument complement)
    self % RTprof_k % nchan_inst = size(channels)
    allocate(self % channels(self % RTprof_k % nchan_inst))
    allocate(self % coefindex(self % RTprof_k % nchan_inst))
    self % coefindex(:) = 0
    self % channels(:) = channels

    jnew = 1
    coefloop: do ii = 1, self % RTprof_k % nchan_inst
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

  end subroutine ufo_radiancerttov_tlad_setup

  ! ------------------------------------------------------------------------------
  subroutine ufo_radiancerttov_tlad_delete(self)
    implicit none

    class(ufo_radiancerttov_tlad), intent(inout) :: self
    integer(kind=jpim)                           :: errorstatus

    self % ltraj = .false.

    call rttov_conf_delete(self % conf)
    if (allocated(self % varin)) deallocate(self % varin)
    if (allocated(self % channels)) deallocate(self % channels)
    if (allocated(self % coefindex)) deallocate(self % coefindex)
    !call rttov_conf_delete(self % conf_traj)

  end subroutine ufo_radiancerttov_tlad_delete

  ! ------------------------------------------------------------------------------
  subroutine ufo_radiancerttov_tlad_settraj(self, geovals, obss, hofxdiags)

    use fckit_mpi_module,   only: fckit_mpi_comm

    implicit none

    class(ufo_radiancerttov_tlad), intent(inout) :: self
    type(ufo_geovals),             intent(in)    :: geovals
    type(c_ptr), value,            intent(in)    :: obss
    type(ufo_geovals),             intent(inout) :: hofxdiags    !non-h(x) diagnostics

    real(c_double)                               :: missing
    type(fckit_mpi_comm)                         :: f_comm

    ! Local Variables
    character(*), parameter                      :: routine_name = 'ufo_radiancerttov_tlad_settraj'
    character(len=800)                           :: message
    type(rttov_chanprof), allocatable            :: chanprof(:)
    type(ufo_geoval), pointer                    :: geoval_temp
    type(rttov_hofxdiags)                        :: hofxdiags_methods

    integer(kind=jpim)                           :: errorstatus ! Return error status of RTTOV subroutine calls

    integer                                      :: iprof_rttov, iprof, ichan, ichan_sim, jchan, ilev, localchan
    integer                                      :: nprof_sim, nprof_max_sim, nchan_total
    integer                                      :: nchan_sim, ivar
    integer                                      :: prof_start, prof_end
    integer                                      :: sensor_index
    integer, allocatable                         :: prof_list(:,:)  ! store list of 'good' profiles

    logical                                      :: jacobian_needed
    real(kind_real), allocatable                 :: sfc_emiss(:,:)
    character(len = MAXVARLEN)                   :: varname

    include 'rttov_k.interface'
    include 'rttov_scatt_ad.interface'

    !Initialisations
    missing = missing_value(missing)

    !Return the name and name length of obsspace communicator (from ioda)
    call obsspace_get_comm(obss, f_comm)

    !! Parse hofxdiags%variables into independent/dependent variables and channel assumed formats:
    ! Note this sets the jacobian_needed flag
    !!   jacobian var -->     <ystr>_jacobian_<xstr>_<chstr>
    !!   non-jacobian var --> <ystr>_<chstr>
    call hofxdiags_methods % parse(hofxdiags, jacobian_needed)

    ! Get number of profiles and levels from geovals
    self % nprofiles = geovals % nlocs
    call ufo_geovals_get_var(geovals, var_ts, geoval_temp)
    self % nlevels = geoval_temp % nval
    nullify(geoval_temp)

    ! Sanity checks
    if (self % nprofiles == 0) return

    ! Allocate RTTOV profiles for ALL geovals for the direct calculation
    write(message,'(2A, I0, A, I0, A)') &
      trim(routine_name), ': Allocating ', self % nprofiles, ' profiles with ', self % nlevels, ' levels'
    call fckit_log%debug(message)
    call self % RTprof_K % alloc_profiles(errorstatus, self % conf, self % nprofiles, self % nlevels, init=.true., asw=1)

    !Assign the atmospheric and surface data from the GeoVaLs
    message = trim(routine_name) // ': Creating RTTOV profiles from geovals'
    call fckit_log%debug(message)
    call self % RTprof_K % setup_rtprof(geovals,obss,self % conf)

    !DAR: Removing sensor_loop until it's demonstrated to be needed and properly tested
    ! at the moment self % channels is a single 1D array so cannot adequately contain more than one set of channels

    ! Read emissivity from obs space if its requested
    if (self % conf % surface_emissivity_group /= "") then
      allocate(sfc_emiss(self % RTprof_K % nchan_inst, self % nprofiles)) ! nchans, nprofiles
      call rttov_read_emissivity_from_obsspace(obss, self % conf % surface_emissivity_group, &
                                               self % channels, sfc_emiss)
    end if

    ! Allocate memory for *ALL* RTTOV_K channels
    ! If RTTOV-SCATT is being used as the obs operator then memory is allocated for profiles_k in mw_scatt too
    write(message,'(2A,I0,A)') &
      trim(routine_name), ': Allocating Trajectory resources for RTTOV K: ', self % nprofiles * self % RTprof_K % nchan_inst, ' total channels'
    call self % RTprof_K % alloc_profiles_k(errorstatus, self % conf, self % nprofiles * self % RTprof_K % nchan_inst, self % nlevels, init=.true., asw=1)

    ! Used for keeping track of profiles for setting emissivity
    allocate(self % RTprof_K % chanprof ( self % nprofiles * self % RTprof_K % nchan_inst ))

    ! Maximum number of profiles to be processed by RTTOV per pass
    if(self % conf % prof_by_prof) then
      nprof_max_sim = 1
    else
      nprof_max_sim = max(1,self % conf % nchan_max_sim / self % RTprof_K % nchan_inst)
    endif
    nprof_sim = min(nprof_max_sim, self % nprofiles)

    ! Determine the total number of radiances to simulate (nchan_sim).
    nchan_sim = nprof_sim * size(self%channels)

    ! Allocate structures for RTTOV direct and K code
    write(message,'(2A,I0,A,I0,A)') &
      trim(routine_name), ': Allocating resources for RTTOV direct (K): ', nprof_sim, ' and ', nchan_sim, ' channels'
    call fckit_log%debug(message)
    call self % RTprof_K % alloc_direct(errorstatus, self % conf, nprof_sim, nchan_sim, self % nlevels, init=.true., asw=1)

    write(message,'(2A,I0,A,I0,A)') &
      trim(routine_name), ': Allocating resources for RTTOV K code: ', nprof_sim, ' and ', nchan_sim, ' channels'
    call fckit_log%debug(message)
    call self % RTprof_K % alloc_k(errorstatus, self % conf, nprof_sim, nchan_sim, self % nlevels, init=.true., asw=1)

    prof_start = 1
    prof_end = self % nprofiles
    nchan_total = 0

    RTTOV_loop : do while (prof_start <= prof_end)

      ! Zero all k code variables.  These arrays are of size prof_end-prof_start
      call self % RTprof_K % zero_k(self % conf, reset_profiles_k=.false.)

      ! Reduce number of simulated profiles/channel if at end of the of profiles to be processed
      nprof_sim = min(nprof_sim, prof_end - prof_start + 1)
      nchan_sim = nprof_sim * size(self%channels)

      ! allocate and initialise local chanprof structure
      allocate(chanprof ( nchan_sim ))
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
        if(any(self % conf % inspect == iprof)) then
          write(*,*) "tlad profile start"
          call self % RTprof_K % print_rtprof(self % conf, iprof)
          write(*,*) "tlad profile done"
        end if

        ! check RTTOV profile will be valid for RTTOV and flag it if it fails the check
        call self % RTprof_K % check_rtprof(self % conf, iprof, errorstatus)

        if (errorstatus == errorstatus_success) then
          ! check sfc_emiss valid if read in
          if (allocated(sfc_emiss)) then
            do ichan = 1, self % RTprof_K % nchan_inst
              if ((sfc_emiss(ichan,iprof) > 1.0) .or. (sfc_emiss(ichan,iprof) < 0.0)) then
                errorstatus = errorstatus_fatal
              end if
            end do
          end if

          prof_list(iprof_rttov,1) = iprof_rttov ! chunk index
          prof_list(iprof_rttov,2) = iprof       ! all-obs index
          do ichan = 1, self % RTprof_K % nchan_inst
            ichan_sim = ichan_sim + 1_jpim
            chanprof(ichan_sim) % prof = iprof_rttov ! this refers to the slice of the RTprofile array passed to RTTOV
            chanprof(ichan_sim) % chan = self % coefindex(ichan)
            self % RTprof_K % chanprof(nchan_total + ichan_sim) % prof = iprof ! this refers to the index of the profile from the geoval
            self % RTprof_K % chanprof(nchan_total + ichan_sim) % chan = self % coefindex(ichan)
          end do
          nchan_sim = ichan_sim

          ! Pick the last valid (unskipped) observation to get the sensor index
          ! but it will always be the same as the first because prof_by_prof is true for nsensors > 1
          ! This is to guard against a bad satellite identifier sneaking in when only processing one instrument.
          sensor_index = self % RTProf_K % sensor_index_array(iprof)
        endif

      end do

      ! Set surface emissivity
      if (allocated(sfc_emiss)) then
        outerloop: do ichan = 1, ichan_sim  ! list of channels*profiles
          do jchan = 1, self % RTprof_K % nchan_inst  ! list of self % channels
            ! if the channel number for this channel * profile == channel number needed
            ! chanprof(ichan) % chan refers to the index in the coefficient file
            if (self % conf % rttov_coef_array(1) % coef % ff_ori_chn(chanprof(ichan) % chan) == self % channels(jchan)) then
              iprof = prof_start + chanprof(ichan) % prof - 1
              self % RTprof_K % emissivity(ichan) % emis_in = sfc_emiss(jchan, iprof)
              self % RTprof_K % calcemis(ichan) = .false.
              if (self % RTprof_K % emissivity(ichan) % emis_in == 0.0) then
                self % RTprof_K % calcemis(ichan) = .true.
              end if
              cycle outerloop
            end if
          end do
        end do outerloop
      else
        call self % RTProf_K % init_default_emissivity(self % conf, prof_list)
      end if

      ! Write out emissivity if checking profile
      if(size(self % conf % inspect) > 0) then
        do ichan = 1, ichan_sim, self % RTprof_K % nchan_inst
          iprof = prof_start + chanprof(ichan) % prof - 1
          if(any(self % conf % inspect == iprof)) then
            write(*,*) "tlad settraj profile ", iprof
            write(*,*) "tlad settraj calcemiss = ",self % RTprof_K % calcemis(ichan:ichan+self%RTprof_K%nchan_inst-1)
            write(*,*) "tlad settraj emissivity in = ",self % RTprof_K % emissivity(ichan:ichan+self%RTprof_K%nchan_inst-1) % emis_in
          end if
        end do
      end if

      deallocate(prof_list)

      ! --------------------------------------------------------------------------
      ! Call RTTOV K model
      ! --------------------------------------------------------------------------

      if (self % conf % do_mw_scatt) then
        call rttov_scatt_ad(                                                                  &
          errorstatus,                                                                        &! out   error flag
          self % conf % mw_scatt % opts,                                                      &! in    options structure
          int (self % nlevels, KIND = jpim),                                                  &
          chanprof(1:nchan_sim),                                                              &! in    LOCAL channel and profile index structure
          self % RTprof_K % mw_scatt % freq_indices(1:nchan_sim),                             &! in    frequency indices
          self % RTProf_K % profiles(prof_start:prof_start + nprof_sim -1),                   &! in    profile array
          self % RTProf_K % mw_scatt % profiles(prof_start:prof_start + nprof_sim -1),        &! in    scattering profile array
          self % conf % rttov_coef_array(sensor_index),                                       &! in    coefficients structure
          self % conf % mw_scatt % coef,                                                      &! in    scatt coefficients structure
          self % RTProf_K % calcemis(1:nchan_sim),                                            &! in    flag for internal emissivity calcs
          self % RTProf_K % emissivity(1:nchan_sim),                                          &! inout input/output emissivities per channel
          self % RTProf_K % profiles_k(nchan_total + 1 : nchan_total + nchan_sim),            &! inout
          self % RTProf_K % mw_scatt % profiles_k(nchan_total + 1 : nchan_total + nchan_sim), &! inout
          self % RTProf_K % emissivity_k(1:nchan_sim),                                        &! inout input/output emissivity jacs per channel
          self % RTProf_K % radiance,                                                         &! inout computed radiances
          self % RTProf_K % radiance_k)                                                        ! inout computed radiance jacobians
      else
        call rttov_k(                                                                   &
          errorstatus,                                                                  &! out   error flag
          chanprof(1:nchan_sim),                                                        &! in channel and profile index structure
          self % conf % rttov_opts,                                                     &! in    options structure
          self % RTprof_K % profiles(prof_start:prof_start + nprof_sim - 1),            &! in    profile array
          self % RTprof_K % profiles_k(nchan_total + 1 : nchan_total + nchan_sim),      &! in    profile array
          self % conf % rttov_coef_array(sensor_index),                                 &! in    coefficients structure
          self % RTprof_K % transmission,                                               &! inout computed transmittances
          self % RTprof_K % transmission_k,                                             &! inout computed transmittances
          self % RTprof_K % radiance,                                                   &! inout computed radiances
          self % RTprof_K % radiance_k,                                                 &! inout computed radiances
          calcemis     = self % RTprof_K % calcemis(1:nchan_sim),                       &! in    flag for internal emissivity calcs
          emissivity   = self % RTprof_K % emissivity(1:nchan_sim),                     &! inout input/output emissivities per channel
          emissivity_k = self % RTprof_K % emissivity_k(1:nchan_sim))                    ! inout input/output emissivities gradients per channel
      end if

      if(self % conf % SatRad_compatibility) then
        ! Zero jacobians if the humidity had been set to the minimum threshold
        if (self % conf % UseMinimumQ) then
          do ichan = 1, nchan_sim, self % RTprof_K % nchan_inst
            iprof = prof_start + chanprof(ichan) % prof - 1
            ! Humidity profile
            do ilev = 1, self % nlevels
              if (self % RTProf_K % q_profile_reset(iprof, ilev)) then
                do localchan = 1, self % RTprof_K % nchan_inst
                  self % RTprof_K % profiles_k(nchan_total + ichan + localchan - 1) % q(ilev) = zero
                end do
              end if
            end do
          end do
        end if

        ! Zero jacobians if selected variables had been set to the minimum threshold
        do ivar = 1, size(self % varin)
          varname = self % varin(ivar)
          select case (trim(varname))

            ! cloud liquid water
            case(var_clw)
              if (self % conf % UseMinimumClw .or. self % conf % MWScattZeroJacPress > 0.0) then
                if (self % conf % do_mw_scatt) then
                  do ichan = 1, nchan_sim, self % RTprof_K % nchan_inst
                    iprof = prof_start + chanprof(ichan) % prof - 1
                    do ilev = 1, self % nlevels
                      if (self % RTProf_K % clw_profile_reset(iprof, ilev)) then
                        do localchan = 1, self % RTprof_K % nchan_inst
                          self % RTprof_K % mw_scatt &
                            % profiles_k(nchan_total + ichan + localchan - 1) &
                            % clw(ilev) = zero
                        end do
                      end if
                    end do
                  end do
                else if (self % conf % rttov_opts % rt_mw % clw_data) then
                  do ichan = 1, nchan_sim, self % RTprof_K % nchan_inst
                    iprof = prof_start + chanprof(ichan) % prof - 1
                    do ilev = 1, self % nlevels
                      if (self % RTProf_K % clw_profile_reset(iprof, ilev)) then
                        do localchan = 1, self % RTprof_K % nchan_inst
                          self % RTprof_K % profiles_k(nchan_total + ichan + localchan - 1) &
                            % clw(ilev) = zero
                        end do
                      end if
                    end do
                  end do
                end if ! self % conf % rttov_opts % rt_ir % addclouds option not included yet
              end if

            ! cloud ice
            case(var_cli)
              if (self % conf % UseMinimumCiw .or. self % conf % MWScattZeroJacPress > 0.0) then
                if (self % conf % do_mw_scatt) then
                  if (self % conf % mw_scatt % use_totalice) then
                    do ichan = 1, nchan_sim, self % RTprof_K % nchan_inst
                      iprof = prof_start + chanprof(ichan) % prof - 1
                      do ilev = 1, self % nlevels
                        if (self % RTProf_K % ciw_profile_reset(iprof, ilev)) then
                          do localchan = 1, self % RTprof_K % nchan_inst
                            self % RTprof_K % mw_scatt &
                              % profiles_k(nchan_total + ichan + localchan - 1) &
                              % totalice(ilev) = zero
                          end do
                        end if
                      end do
                    end do
                  else
                    do ichan = 1, nchan_sim, self % RTprof_K % nchan_inst
                      iprof = prof_start + chanprof(ichan) % prof - 1
                      do ilev = 1, self % nlevels
                        if (self % RTProf_K % ciw_profile_reset(iprof, ilev)) then
                          do localchan = 1, self % RTprof_K % nchan_inst
                            self % RTprof_K % mw_scatt &
                              % profiles_k(nchan_total + ichan + localchan - 1) &
                              % ciw(ilev) = zero
                          end do
                        end if
                      end do
                    end do
                  end if
                else
                  message = &
                    'ufo_radiancerttov_tlad_settraj: Cloud Ice Water only supported for RTTOV-SCATT'
                  call abor1_ftn(message)
                end if
              end if

          end select
        end do
      end if

      if(size(self % conf % inspect) > 0) then
        do ichan = 1, ichan_sim, self % RTprof_K % nchan_inst
          iprof = prof_start + chanprof(ichan) % prof - 1
          if(any(self % conf % inspect == iprof)) then
          write(*,*) "tlad profile ", iprof
          write(*,*) "tlad emissivity out = ",self % RTprof_K % emissivity(ichan:ichan+self%RTprof_K%nchan_inst-1) % emis_out
          write(*,*) "tlad hofx out = ", self % RTprof_K % radiance % bt(ichan:ichan+self%RTprof_K%nchan_inst-1)
          end if
        end do
      end if

      if ( errorstatus /= errorstatus_success ) then
        write(message,'(3A, I6, A, I6)') trim(routine_name), 'after rttov_k: error\n', 'skipping profiles ', &
                                         prof_start, ' -- ', prof_start + nprof_sim - 1
        call fckit_log%info(message)
      else
        ! Put simulated diagnostics into hofxdiags
        ! ----------------------------------------------
        if(hofxdiags % nvar > 0) &
          call hofxdiags_methods % populate(self % RTprof_K, chanprof, self % conf, prof_start, hofxdiags)
      end if

      ! increment profile and channel counters
      nchan_total = nchan_total + nchan_sim
      prof_start = prof_start + nprof_sim

      self % nchan_total = nchan_total
      deallocate (chanprof)
    end do RTTOV_loop

    !    end do Sensor_Loop
    ! Deallocate structures for rttov_direct
    call self % RTprof_K % alloc_k(errorstatus, self % conf, -1, -1, -1, asw=0)
    call self % RTprof_K % alloc_direct(errorstatus, self % conf, -1, -1, -1, asw=0)
    call self % RTprof_K % alloc_profiles(errorstatus, self % conf, -1, -1, asw=0)


    ! Set flag that the tracectory was set
    ! ------------------------------------
    self % ltraj = .true.

  end subroutine ufo_radiancerttov_tlad_settraj

  ! ------------------------------------------------------------------------------
  subroutine ufo_radiancerttov_simobs_tl(self, geovals, obss, nvars, nlocs, hofx)

    implicit none

  class(ufo_radiancerttov_tlad), intent(in)    :: self
    type(ufo_geovals),           intent(in)    :: geovals
    type(c_ptr), value,          intent(in)    :: obss
    integer,                     intent(in)    :: nvars, nlocs
    real(c_double),              intent(inout) :: hofx(nvars, nlocs)

    type(ufo_geoval), pointer                  :: geoval_d
    integer                                    :: ichan, jchan, prof, jspec, ivar
    character(len = MAXVARLEN)                 :: varname
    character(len = *), parameter              :: myname_="ufo_radiancerttov_simobs_tl"
    character(len=800)                         :: message

    ! Initial checks
    ! --------------

    ! Check if trajectory was set
    if (.not. self % ltraj) then
      message = myname_ // ' trajectory wasnt set!'
      call abor1_ftn(message)
    end if

    ! Check if nlocs is consistent in geovals & hofx
    if (geovals % nlocs /= self % nprofiles) then
      message = myname_ // ' error: nlocs inconsistent!'
      call abor1_ftn(message)
    end if

    ! Initialize hofx
    ! ---------------
    hofx(:,:) = zero

    do ivar = 1, size(self % varin)
      varname = self % varin(ivar)
      select case (trim(varname))
        ! Variables with nlevels
        case (var_ts, var_q, var_mixr, var_clw, var_cli)
          call ufo_geovals_get_var(geovals, trim(varname), geoval_d)

          ! Check model levels is consistent in geovals
          if (geoval_d % nval /= self % nlevels) then
            message = myname_ // ' error: layers inconsistent!'
            call abor1_ftn(message)
          end if

          do ichan = 1, self % nchan_total, size(self % channels)
            prof = self % RTprof_K % chanprof(ichan) % prof
            do jchan = 1, size(self % channels)
              if (trim(varname) == var_ts) then
                hofx(jchan, prof) = hofx(jchan, prof) + &
                  sum(self % RTprof_K % profiles_k(ichan+jchan-1) % t(:) * &
                      geoval_d % vals(1:geoval_d % nval, prof))
              else if (trim(varname) == var_q) then
                ! scale_fac = 1 if Jacobain is in kg/kg
                ! scale_fac converts from ppmv to kg/kg if Jacobain wrt ppmv
                hofx(jchan, prof) = hofx(jchan, prof) + &
                  sum(self % RTprof_K % profiles_k(ichan+jchan-1) % q(:) * &
                      geoval_d % vals(1:geoval_d % nval, prof)) * self%conf%scale_fac(gas_id_watervapour)
              else if (trim(varname) == var_mixr) then
                ! scale_fac = 1 if Jacobain is in kg/kg
                ! scale_fac converts from ppmv to kg/kg if Jacobain wrt ppmv
                ! humidity_mixing_ratio (var_mixr) is in g/kg therefore need additional conversion factor
                hofx(jchan, prof) = hofx(jchan, prof) + &
                  sum(self % RTprof_K % profiles_k(ichan+jchan-1) % q(:) * &
                      geoval_d % vals(1:geoval_d % nval, prof)) * self%conf%scale_fac(gas_id_watervapour) / g_to_kg
              else if (trim(varname) == var_clw) then
                if (self % conf % do_mw_scatt) then
                  hofx(jchan, prof) = hofx(jchan, prof) + &
                    sum(self % RTprof_K % mw_scatt % profiles_k(ichan+jchan-1) % clw(:) * &
                        geoval_d % vals(1:geoval_d % nval, prof))
                else if (self % conf % rttov_opts % rt_mw % clw_data) then
                  hofx(jchan, prof) = hofx(jchan, prof) + &
                    sum(self % RTprof_K % profiles_k(ichan+jchan-1) % clw(:) * &
                        geoval_d % vals(1:geoval_d % nval, prof))
                end if
              else if (trim(varname) == var_cli) then
                if (self % conf % do_mw_scatt) then
                  if (self % conf % mw_scatt % use_totalice) then
                    hofx(jchan, prof) = hofx(jchan, prof) + &
                      sum(self % RTprof_K % mw_scatt % profiles_k(ichan+jchan-1) % totalice(:) * &
                          geoval_d % vals(1:geoval_d % nval, prof))
                  else
                    hofx(jchan, prof) = hofx(jchan, prof) + &
                      sum(self % RTprof_K % mw_scatt % profiles_k(ichan+jchan-1) % ciw(:) * &
                          geoval_d % vals(1:geoval_d % nval, prof))
                  end if
                else
                  message = 'ufo_radiancerttov_simobs_tl: Cloud Ice Water only supported for RTTOV-SCATT'
                  call abor1_ftn(message)
                end if
              end if
            end do
          end do
        ! Variables with 1 level - surface
        case (var_sfc_t2m, var_sfc_q2m, var_sfc_u10, var_sfc_v10, var_sfc_tskin)
          call ufo_geovals_get_var(geovals, trim(varname), geoval_d)

          do ichan = 1, self % nchan_total, size(self % channels)
            prof = self % RTprof_K % chanprof(ichan) % prof
            do jchan = 1, size(self%channels)
              if (trim(varname) == var_sfc_t2m) then
                hofx(jchan, prof) = hofx(jchan, prof) + &
                  self % RTprof_K % profiles_k(ichan+jchan-1) % s2m % t * geoval_d % vals(1, prof)
              else if (trim(varname) == var_sfc_q2m) then
                ! scale_fac = 1 if Jacobain is in kg/kg
                ! scale_fac converts from ppmv to kg/kg if Jacobain wrt ppmv
                hofx(jchan, prof) = hofx(jchan, prof) + &
                  self % RTprof_K % profiles_k(ichan+jchan-1) % s2m % q * &
                  geoval_d % vals(1,prof) * self%conf%scale_fac(gas_id_watervapour)
              else if (trim(varname) == var_sfc_u10) then
                hofx(jchan, prof) = hofx(jchan, prof) + &
                  self % RTprof_K % profiles_k(ichan+jchan-1) % s2m % u * geoval_d % vals(1, prof)
              else if (trim(varname) == var_sfc_v10) then
                hofx(jchan, prof) = hofx(jchan, prof) + &
                  self % RTprof_K % profiles_k(ichan+jchan-1) % s2m % v * geoval_d % vals(1, prof)
              else if (trim(varname) == var_sfc_tskin) then
                hofx(jchan, prof) = hofx(jchan, prof) + &
                  self % RTprof_K % profiles_k(ichan+jchan-1) % skin % t * geoval_d % vals(1, prof)
              end if
            end do
          end do

        case default
          message = trim(varname) // ' in increment list but not setup in ' // myname_
          call abor1_ftn(message)
      end select

    end do ! main variable loop

  end subroutine ufo_radiancerttov_simobs_tl

  ! ------------------------------------------------------------------------------
  subroutine ufo_radiancerttov_simobs_ad(self, geovals, obss, nvars, nlocs, hofx)

    implicit none

    class(ufo_radiancerttov_tlad), intent(in)    :: self
    type(ufo_geovals),             intent(inout) :: geovals
    type(c_ptr), value,            intent(in)    :: obss
    integer,                       intent(in)    :: nvars, nlocs
    real(c_double),                intent(in)    :: hofx(nvars, nlocs)

    type(ufo_geoval), pointer                    :: geoval_d
    real(c_double)                               :: missing
    integer                                      :: ichan, jchan, prof, jspec, ivar
    character(len = MAXVARLEN)                   :: varname
    character(len = *), parameter                :: myname_ = "ufo_radiancerttov_simobs_ad"
    character(len=800)                           :: message

    ! Set missing value
    missing = missing_value(missing)

    ! Initial checks
    ! --------------

    ! Check if trajectory was set
    if (.not. self % ltraj) then
      message = myname_ // ' trajectory was not set!'
      call abor1_ftn(message)
    end if

    ! Check if nlocs is consistent in geovals & hofx
    if (geovals % nlocs /= self % nprofiles) then
      message = myname_ // ' error: nlocs inconsistent!'
      call abor1_ftn(message)
    end if

    do ivar = 1, size(self % varin)
      varname = self % varin(ivar)
      select case (trim(varname))
        ! Variables with nlevels
        case (var_ts, var_q, var_mixr, var_clw, var_cli)
          call ufo_geovals_get_var(geovals, trim(varname), geoval_d)

          do ichan = 1, self % nchan_total, size(self % channels)
            prof = self % RTprof_K % chanprof(ichan) % prof
            do jchan = 1, size(self % channels)
              if (hofx(jchan, prof) /= missing) then
                if (trim(varname) == var_ts) then
                  geoval_d % vals(:, prof) = geoval_d % vals(:, prof) + &
                    self % RTprof_K % profiles_k(ichan+jchan-1) % t(:) * hofx(jchan, prof)
                else if (trim(varname) == var_q) then
                  ! scale_fac = 1 if Jacobain is in kg/kg
                  ! scale_fac converts from ppmv to kg/kg if Jacobain wrt ppmv
                  geoval_d % vals(:, prof) = geoval_d % vals(:, prof) + &
                    self % RTprof_K % profiles_k(ichan+jchan-1) % q(:) * hofx(jchan, prof) * &
                    self%conf%scale_fac(gas_id_watervapour)
                else if (trim(varname) == var_mixr) then
                  ! scale_fac = 1 if Jacobain is in kg/kg
                  ! scale_fac converts from ppmv to kg/kg if Jacobain wrt ppmv
                  ! humidity_mixing_ratio (var_mixr) is in g/kg therefore need additional conversion factor
                  geoval_d % vals(:, prof) = geoval_d % vals(:, prof) + &
                    (self % RTprof_K % profiles_k(ichan+jchan-1) % q(:) * hofx(jchan, prof)) * &
                    self%conf%scale_fac(gas_id_watervapour) / g_to_kg
                else if (trim(varname) == var_clw) then
                  if (self % conf % do_mw_scatt) then
                    geoval_d % vals(:, prof) = geoval_d % vals(:, prof) + &
                      self % RTprof_K % mw_scatt % profiles_k(ichan+jchan-1) % clw(:) * hofx(jchan, prof)
                  else if (self % conf % rttov_opts % rt_mw % clw_data) then
                    geoval_d % vals(:, prof) = geoval_d % vals(:, prof) + &
                      self % RTprof_K % profiles_k(ichan+jchan-1) % clw(:) * hofx(jchan, prof)
                  end if
                else if (trim(varname) == var_cli) then
                  if (self % conf % do_mw_scatt) then
                    if (self % conf % mw_scatt % use_totalice) then
                      geoval_d % vals(:, prof) = geoval_d % vals(:, prof) + &
                        self % RTprof_K % mw_scatt % profiles_k(ichan+jchan-1) % totalice(:) * hofx(jchan, prof)
                    else
                      geoval_d % vals(:, prof) = geoval_d % vals(:, prof) + &
                        self % RTprof_K % mw_scatt % profiles_k(ichan+jchan-1) % ciw(:) * hofx(jchan, prof)
                    end if
                  else
                    message = 'ufo_radiancerttov_simobs_ad: Cloud Ice Water only supported for RTTOV-SCATT'
                    call abor1_ftn(message)
                  end if
                end if
              endif
            enddo
          end do
        ! Variables with 1 level - surface
        case (var_sfc_t2m, var_sfc_q2m, var_sfc_u10, var_sfc_v10, var_sfc_tskin)
          call ufo_geovals_get_var(geovals, trim(varname), geoval_d)

          do ichan = 1, self % nchan_total, size(self % channels)
            prof = self % RTprof_K % chanprof(ichan) % prof
            do jchan = 1, size(self%channels)
              if (hofx(jchan, prof) /= missing) then
                if (trim(varname) == var_sfc_t2m) then
                  geoval_d % vals(1, prof) = geoval_d % vals(1, prof) + &
                    self % RTprof_K % profiles_k(ichan+jchan-1) % s2m % t * hofx(jchan, prof)
                else if (trim(varname) == var_sfc_q2m) then
                  ! scale_fac = 1 if Jacobain is in kg/kg
                  ! scale_fac converts from ppmv to kg/kg if Jacobain wrt ppmv
                  geoval_d % vals(1, prof) = geoval_d % vals(1, prof) + &
                    self % RTprof_K % profiles_k(ichan+jchan-1) % s2m % q * &
                    hofx(jchan, prof) * self%conf%scale_fac(gas_id_watervapour)
                else if (trim(varname) == var_sfc_u10) then
                  geoval_d % vals(1, prof) = geoval_d % vals(1, prof) + &
                    self % RTprof_K % profiles_k(ichan+jchan-1) % s2m % u * hofx(jchan, prof)
                else if (trim(varname) == var_sfc_v10) then
                  geoval_d % vals(1, prof) = geoval_d % vals(1, prof) + &
                    self % RTprof_K % profiles_k(ichan+jchan-1) % s2m % v * hofx(jchan, prof)
                else if (trim(varname) == var_sfc_tskin) then
                  geoval_d % vals(1, prof) = geoval_d % vals(1, prof) + &
                    self % RTprof_K % profiles_k(ichan+jchan-1) % skin % t * hofx(jchan, prof)
                end if
              end if
            end do
          end do

        case default
          message = trim(varname) // ' in increment list but not setup in ' // myname_
          call abor1_ftn(message)
      end select

    end do ! main variable loop

    ! Once all geovals set replace flag
    ! ---------------------------------
    if (.not. geovals % linit ) geovals % linit=.true.

  end subroutine ufo_radiancerttov_simobs_ad

  ! ------------------------------------------------------------------------------

end module ufo_radiancerttov_tlad_mod
