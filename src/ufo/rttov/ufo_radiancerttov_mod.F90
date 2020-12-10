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
  use rttov_const, only : errorstatus_success
  use rttov_unix_env

  implicit none
  private

  !> Fortran derived type for the observation type
  type, public :: ufo_radiancerttov
    private
    character(len=MAXVARLEN), public, allocatable :: varin(:)  ! variables requested from the model
    integer, allocatable                          :: channels(:)
    type(rttov_conf)                              :: conf
    type(ufo_rttov_io)                            :: RTProf
  contains
    procedure :: setup  => ufo_radiancerttov_setup
    procedure :: delete => ufo_radiancerttov_delete
    procedure :: simobs => ufo_radiancerttov_simobs
  end type ufo_radiancerttov

contains

  ! ------------------------------------------------------------------------------
  subroutine ufo_radiancerttov_setup(self, f_confOper, channels)

    implicit none
    class(ufo_radiancerttov), intent(inout) :: self
    type(fckit_configuration), intent(in)   :: f_confOper
    integer(c_int),            intent(in)   :: channels(:)  !List of channels to use

    type(fckit_configuration)               :: f_confOpts ! RTcontrol
    integer                                 :: ind, jspec

    call f_confOper % get_or_die("obs options",f_confOpts)

    call rttov_conf_setup(self % conf,f_confOpts,f_confOper)

    if ( ufo_vars_getindex(self%conf%Absorbers, var_mixr) < 1 .and. &
      ufo_vars_getindex(self%conf%Absorbers, var_q)    < 1 ) then
      write(message,*) 'ufo_radiancerttov_setup error: H2O must be included in RTTOV Absorbers'
      call abor1_ftn(message)
    end if

    nvars_in = size(varin_default) + self%conf%ngas

    allocate(self%varin(nvars_in))
    self%varin(1:size(varin_default)) = varin_default
    ind = size(varin_default) + 1

    !Use list of Absorbers from conf
    do jspec = 1, self%conf%ngas
      self%varin(ind) = self%conf%Absorbers(jspec)
      ind = ind + 1
    end do

    ! save channels
    allocate(self%channels(size(channels)))
    self%channels(:) = channels(:)

    write(message,'(A, 2I6)') 'Finished setting up rttov'
    call fckit_log%info(message)

  end subroutine ufo_radiancerttov_setup

  ! ------------------------------------------------------------------------------
  subroutine ufo_radiancerttov_delete(self)
    implicit none
    class(ufo_radiancerttov), intent(inout) :: self

    call rttov_conf_delete(self%conf)

  end subroutine ufo_radiancerttov_delete

  ! ------------------------------------------------------------------------------
  subroutine ufo_radiancerttov_simobs(self, geovals, obss, nvars, nlocs, &
                                      hofx, hofxdiags, ob_info)
    use fckit_mpi_module,   only: fckit_mpi_comm
    use ufo_rttovonedvarcheck_ob_mod

    implicit none

    class(ufo_radiancerttov), intent(inout) :: self
    type(ufo_geovals),        intent(in)    :: geovals
    type(c_ptr), value,       intent(in)    :: obss
    integer,                  intent(in)    :: nvars, nlocs

    real(c_double),        intent(inout)    :: hofx(nvars,nlocs)
    type(ufo_geovals),     intent(inout)    :: hofxdiags    !non-h(x) diagnostics
    type(ufo_rttovonedvarcheck_ob), optional, intent(inout) :: ob_info

    real(c_double)                          :: missing
    type(fckit_mpi_comm)                    :: f_comm

    ! Local Variables
    character(*), parameter                 :: ROUTINE_NAME = 'ufo_radiancerttov_simobs'
    type(ufo_geoval), pointer               :: temp

    integer                                 :: nprofiles
    integer(kind=jpim)                      :: errorstatus              ! Return error status of RTTOV subroutine calls

    integer                                 :: i_inst, nlevels, nchan_total, ichan, iprof, prof
    integer                                 :: nprof_sim, nprof_max_sim
    integer                                 :: prof_start, prof_end

    logical                                 :: jacobian_needed
    logical, allocatable                    :: skip_profiles(:)

    logical                                 :: layer_quantities

    include 'rttov_direct.interface'
    include 'rttov_k.interface'
    include 'rttov_print_profile.interface'
    include 'rttov_user_profile_checkinput.interface'

    !DAR: What is this?
    call obsspace_get_comm(obss, f_comm)

    ! Get number of profile and layers from geovals
    ! ---------------------------------------------

    nprofiles = geovals % nlocs
    if (ufo_vars_getindex(geovals%variables, var_ts) > 0) then
      call ufo_geovals_get_var(geovals, var_ts, temp)
      nlevels = temp % nval
      layer_quantities = .false.
    endif

    nullify(temp)

    missing = missing_value(missing)
    hofx    = missing

    errorstatus = 0_jpim
    nchan_total = 0

    !DARFIX: This isn't ideal because it's not going to work for multiple instruments but we'll deal with that later
    nchan_inst = size(self % channels)

    !! Parse hofxdiags%variables into independent/dependent variables and channel
    !! assumed formats:
    ! Note this sets jacobian_needed
    !!   jacobian var -->     <ystr>_jacobian_<xstr>_<chstr>
    !!   non-jacobian var --> <ystr>_<chstr>
    call parse_hofxdiags(hofxdiags, jacobian_needed)

    Sensor_Loop:do i_inst = 1, self % conf % nSensors

      ! Ensure the options and coefficients are consistent
      call rttov_user_options_checkinput(errorstatus, self % conf % rttov_opts, &
        self % conf % rttov_coef_array(i_inst))

      if (errorstatus /= errorstatus_success) then
        write(message,'(A, I6)') 'after rttov_user_options_checkinput: error = ',&
          errorstatus
        call fckit_log%info(message)
      end if

      ! keep journal of which profiles have no obs data these will be skipped
      allocate(Skip_Profiles(nprofiles))
      if (present(ob_info)) then
        Skip_Profiles(:) = .false.
      else
        ! keep journal of which profiles have no obs data these will be skipped
        Skip_Profiles(:) = .false.
        call ufo_rttov_skip_profiles(nProfiles,nchan_inst,self%channels,obss,Skip_Profiles)
      end if

      ! Determine the total number of radiances to simulate (nchan_sim).
      ! In this example we simulate all specified channels for each profile, but
      ! in general one can simulate a different number of channels for each profile.

      nprof_max_sim = self % conf % nchan_max_sim / nchan_inst
      nprof_sim = min(nprof_max_sim, nprofiles)

      prof_start = 1
      prof_end = nprofiles

      !DARFIX should actually be packed count of skip_profiles * SIZE(channels)
      nprof_sim = min(nprof_sim, prof_end - prof_start + 1)
      nchan_sim = nprof_sim * size(self%channels)

      ! --------------------------------------------------------------------------
      ! Allocate RTTOV input and output structures
      ! --------------------------------------------------------------------------
      if (.not. jacobian_needed) then
        ! allocate RTTOV resources
        call self % RTprof % alloc(errorstatus, self % conf, nprof_sim, nchan_sim, nlevels, init=.true., asw=1)
      else
        call self % RTprof % alloc(errorstatus, self % conf, nprof_sim, nchan_sim, nlevels, init=.true., asw=2)
      endif

      self % RTProf % profiles(:) % skin % surftype = -1_jpim

      do while (prof_start <= prof_end)
        self % RTprof % chanprof(:) % prof = 0 
        self % RTprof % chanprof(:) % chan = 0 

        nchan_sim = 0_jpim
        nprof_sim = min(nprof_sim, prof_end - prof_start + 1)

        ! --------------------------------------------------------------------------
        ! Build the list of profile/channel indices in chanprof
        ! --------------------------------------------------------------------------

        do iprof = 1, min(nprof_sim, prof_end - prof_start + 1)
          if(.not. Skip_Profiles(prof_start + iprof - 1)) then
            do ichan = 1, nchan_inst

              nchan_sim = nchan_sim + 1_jpim
              self % RTprof % chanprof(nchan_sim) % prof = iprof
              self % RTprof % chanprof(nchan_sim) % chan = self % channels(ichan)
            end do
          else
            if (debug) write(*,*) 'skipping ', iprof, prof_start + iprof - 1
          end if
        end do

        !Assign the data from the GeoVaLs
        !--------------------------------
        
        if(present(ob_info)) then
          call load_atm_data_rttov(geovals,obss,self % RTprof % profiles,prof_start,self % conf,layer_quantities,ob_info=ob_info)
          call load_geom_data_rttov(obss,self % RTprof % profiles,prof_start,ob_info=ob_info)
        else
          call load_atm_data_rttov(geovals,obss,self % RTprof % profiles,prof_start,self%conf,layer_quantities)
          call load_geom_data_rttov(obss,self % RTprof % profiles,prof_start)
        end if

        ! --------------------------------------------------------------------------
        ! Set surface emissivity
        ! --------------------------------------------------------------------------

        call self % RTProf % init_emissivity(self % conf)

        if(self % conf % RTTOV_profile_checkinput) then
          if (self % conf % inspect > 0) then
            call rttov_print_profile(self % RTprof % profiles(self % conf % inspect))
          endif

          ! no error checking, could check multiple profiles in loop but what would be the point
          call rttov_user_profile_checkinput(rttov_errorstatus, &
            self % conf % rttov_opts, &
            self % conf % rttov_coef_array(i_inst), &
            self % RTprof % profiles(self % conf % inspect))
        endif


        ! --------------------------------------------------------------------------
        ! Call RTTOV model
        ! --------------------------------------------------------------------------

        if (jacobian_needed) then

          ! Inintialize the K-matrix INPUT so that the results are dTb/dx
          ! -------------------------------------------------------------
          self % RTprof % emissivity_k(:) % emis_out = 0
          self % RTprof % emissivity_k(:) % emis_in = 0
          self % RTprof % emissivity(:) % emis_out = 0
          self % RTprof % radiance_k % bt(:) = 1
          self % RTprof % radiance_k % total(:) = 1

          call rttov_k(                              &
            errorstatus,                             &! out   error flag
            self % RTProf % chanprof(nchan_total + 1:nchan_total + nchan_sim), &! in LOCAL channel and profile index structure
            self % conf % rttov_opts,                     &! in    options structure
            self % RTProf % profiles,                                &! in    profile array
            self % RTProf % profiles_k(nchan_total + 1 : nchan_total + nchan_sim), &! in    profile array
            self % conf % rttov_coef_array(i_inst), &! in    coefficients structure
            self % RTProf % transmission,                            &! inout computed transmittances
            self % RTProf % transmission_k,                          &! inout computed transmittances
            self % RTProf % radiance,                                &! inout computed radiances
            self % RTProf % radiance_k,                              &! inout computed radiances
            calcemis    = self % RTProf % calcemis,                  &! in    flag for internal emissivity calcs
            emissivity  = self % RTProf % emissivity,                &!, &! inout input/output emissivities per channel
            emissivity_k = self % RTProf % emissivity_k)!,           &! inout input/output emissivities per channel      

          if ( errorstatus /= errorstatus_success ) then
            write(message,'(A, 2I6)') 'after rttov_k: error ', errorstatus, i_inst
            call abor1_ftn(message)
          end if

        else
          call rttov_direct(                         &
            errorstatus,                             &! out   error flag
            self % RTprof % chanprof(1:nchan_sim),        &! in    channel and profile index structure
            self % conf % rttov_opts,                     &! in    options structure
            self % RTProf % profiles(1:nprof_sim),                                &! in    profile array
            self % conf % rttov_coef_array(i_inst), &! in    coefficients structure
            self % RTProf % transmission,                            &! inout computed transmittances
            self % RTProf % radiance,                                &! inout computed radiances
            calcemis    = self % RTProf % calcemis(1:nchan_sim),                  &! in    flag for internal emissivity calcs
            emissivity  = self % RTProf % emissivity(1:nchan_sim))!,              &! inout input/output emissivities per channel

          if ( errorstatus /= errorstatus_success ) then
            write(message,'(A, 2I6)') 'after rttov_direct: error ', errorstatus, i_inst
            call abor1_ftn(message)
          end if
        endif ! jacobian_needed

        !DARFIX: This should be available only when debugging
        ! Write out progress through batch
        if (debug) write(*,'(A1, i0, A, i0)',ADVANCE="NO") achar(13), prof_start+nprof_sim-1, ' locations processed out of ', geovals%nlocs

        ! Put simulated brightness temperature into hofx
        ! ----------------------------------------------
        do ichan=1, nchan_sim, size(self%channels)
          prof = self % RTProf % chanprof(ichan)%prof
          hofx(1:size(self%channels),prof_start + prof - 1) = self % RTprof % radiance % bt(ichan:ichan+size(self%channels)-1)
        enddo

        ! Put simulated diagnostics into hofxdiags
        ! ----------------------------------------------
        if(hofxdiags%nvar > 0) call populate_hofxdiags(self % RTProf, self % RTProf % chanprof, hofxdiags)

        nchan_total = nchan_total + nchan_sim
        prof_start = prof_start + nprof_sim

      end do

      ! Deallocate structures for rttov_direct

      call self % RTprof % alloc(errorstatus, self % conf, nprof_sim, nchan_sim, nlevels, asw=0)

      if (errorstatus /= errorstatus_success) then
        write(message,'(A, 2I6)') &
          'after rttov_alloc_direct (deallocation): errorstatus, i_inst =', &
          errorstatus, i_inst
        call abor1_ftn(message)
      end if

    end do Sensor_Loop
    write(*,*)

  end subroutine ufo_radiancerttov_simobs

  ! ------------------------------------------------------------------------------

end module ufo_radiancerttov_mod
