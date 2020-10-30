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

  use ufo_geovals_mod, only: ufo_geovals, ufo_geoval, ufo_geovals_get_var
  use ufo_vars_mod
  use ufo_radiancerttov_utils_mod

  use rttov_types
  use rttov_const, only : errorstatus_success
  use rttov_unix_env

  implicit none
  private

  !> Fortran derived type for radiancerttov trajectory
  type, public :: ufo_radiancerttov_tlad
    private
    character(len=MAXVARLEN), public, allocatable :: varin(:)  ! variables requested from the model
    integer, allocatable                          :: channels(:)
    type(rttov_conf)                              :: conf
    type(rttov_conf)                              :: conf_traj
    type(ufo_rttov_io)                            :: RTProf_K

    integer                                       :: nprofiles
    integer                                       :: nchan_total
    integer                                       :: nlevels

    logical                                       :: ltraj
    logical, allocatable                          :: Skip_Profiles(:)

  contains
    procedure :: setup  => ufo_radiancerttov_tlad_setup
    procedure :: delete  => ufo_radiancerttov_tlad_delete
    procedure :: settraj => ufo_radiancerttov_tlad_settraj
    procedure :: simobs_tl  => ufo_radiancerttov_simobs_tl
    procedure :: simobs_ad  => ufo_radiancerttov_simobs_ad
  end type ufo_radiancerttov_tlad

  character(len=maxvarlen), dimension(1), parameter :: varin_default_tlad = &
    (/var_ts/)

contains

  ! ------------------------------------------------------------------------------
  subroutine ufo_radiancerttov_tlad_setup(self, f_confOper, channels)
    implicit none

    class(ufo_radiancerttov_tlad), intent(inout) :: self
    type(fckit_configuration), intent(in)        :: f_confOper
    integer(c_int),            intent(in)        :: channels(:)  !List of channels to use

    type(fckit_configuration)                    :: f_confOpts ! RTcontrol
    type(fckit_configuration)                    :: f_confLinOper
    integer                                      :: nvars_in
    integer                                      :: ind, jspec

    call f_confOper % get_or_die("obs options",f_confOpts)

    call rttov_conf_setup(self % conf_traj, f_confOpts,f_confOper)

    if ( f_confOper%has("linear obs operator") ) then
      call f_confOper%get_or_die("linear obs operator",f_confLinOper)
      call rttov_conf_setup(self%conf, f_confOpts, f_confLinOper)
    else
      call rttov_conf_setup(self%conf, f_confOpts, f_confOper)
    end if

    !DAR what is the RTTOV equivalant of making sure that humidity and ozone data are present
    if ( ufo_vars_getindex(self%conf%Absorbers, var_mixr) < 1 .and. &
      ufo_vars_getindex(self%conf%Absorbers, var_q)    < 1 ) then
      write(message,*) 'ufo_radiancerttov_setup error: H2O must be included in RTTOV Absorbers'
      call abor1_ftn(message)
    end if

    nvars_in = size(varin_default_tlad) + self%conf%ngas + 5 ! 5 near-surface parameters

    allocate(self%varin(nvars_in))
    self%varin(1:size(varin_default_tlad)) = varin_default_tlad
    ind = size(varin_default_tlad) + 1

    do jspec = 1, self%conf%ngas
      self%varin(ind) = self%conf%Absorbers(jspec)
      ind = ind + 1
    end do

    self%varin(ind) = var_sfc_t2m
    ind = ind + 1
    self%varin(ind) = var_sfc_q2m
    ind = ind + 1
    self%varin(ind) = var_sfc_tskin
    ind = ind + 1
    self%varin(ind) = var_u
    ind = ind + 1
    self%varin(ind) = var_v

    ! save channels
    allocate(self%channels(size(channels)))
    self%channels(:) = channels(:)

  end subroutine ufo_radiancerttov_tlad_setup

  ! ------------------------------------------------------------------------------
  subroutine ufo_radiancerttov_tlad_delete(self)
    implicit none

    class(ufo_radiancerttov_tlad), intent(inout) :: self
    integer(kind=jpim)                           :: errorstatus

    self % ltraj = .false.

    call rttov_conf_delete(self % conf)
    call rttov_conf_delete(self % conf_traj)

    if (allocated(self%Skip_Profiles)) deallocate(self%Skip_Profiles)

  end subroutine ufo_radiancerttov_tlad_delete

  ! ------------------------------------------------------------------------------
  subroutine ufo_radiancerttov_tlad_settraj(self, geovals, obss, hofxdiags)

    use fckit_mpi_module,   only: fckit_mpi_comm

    implicit none

    class(ufo_radiancerttov_tlad), intent(inout) :: self
    type(ufo_geovals),             intent(in)    :: geovals
    type(c_ptr), value,            intent(in)    :: obss
    type(ufo_geovals),             intent(inout) :: hofxdiags    !non-h(x) diagnostics

    type(fckit_mpi_comm)                         :: f_comm

    ! Local Variables
    character(*), parameter                      :: ROUTINE_NAME = 'ufo_radiancerttov_tlad_settraj'
    type(ufo_geoval), pointer                    :: temp

    integer(kind=jpim)                           :: errorstatus ! Return error status of RTTOV subroutine calls

    integer                                      :: nchan_max_sim, nchan_count, nchan_total
    integer                                      :: nprof_sim, nprof_max_sim

    integer                                      :: iprof, ichan, i_inst
    integer                                      :: prof_start, prof_end

    logical                                      :: layer_quantities
    logical                                      :: jacobian_needed

    include 'rttov_k.interface'
    include 'rttov_print_profile.interface'
    include 'rttov_user_profile_checkinput.interface'

    !DAR: What is this?
    call obsspace_get_comm(obss, f_comm)

    ! Get number of profile and layers from geovals
    ! ---------------------------------------------
    self % nprofiles = geovals % nlocs
    if (ufo_vars_getindex(geovals%variables, var_ts) > 0) then
      call ufo_geovals_get_var(geovals, var_ts, temp)
      self % nlevels = temp % nval ! lfric passing nlevels
      layer_quantities = .false.
    endif

    nullify(temp)

    errorstatus = 0_jpim
    nchan_count = 0
    nchan_total = 0

    nchan_inst = size(self % channels)

    nchan_max_sim = self % nprofiles * nchan_inst ! Maximum number of channels to pass to RTTOV to simulate

    !! Parse hofxdiags%variables into independent/dependent variables and channel
    !! assumed formats:
    !!   jacobian var -->     <ystr>_jacobian_<xstr>_<chstr>
    !!   non-jacobian var --> <ystr>_<chstr>
    call parse_hofxdiags(hofxdiags, jacobian_needed)

    Sensor_Loop:do i_inst = 1, self % conf % nSensors

      ! Ensure the options and coefficients are consistent
      call rttov_user_options_checkinput(errorstatus, self % conf % rttov_opts, &
        self % conf % rttov_coef_array(i_inst))

      if (errorstatus /= errorstatus_success) then
        write(message,'(A, A,I6)') trim(ROUTINE_NAME), 'after rttov_user_options_checkinput: error = ',&
          errorstatus
        call fckit_log%info(message)
      end if

      ! keep journal of which profiles have no obs data these will be skipped
      allocate(self % Skip_Profiles(self % nprofiles))
      call ufo_rttov_skip_profiles(self%nProfiles,nchan_inst,self%channels,obss,self%Skip_Profiles)

      ! Determine the total number of radiances to simulate (nchanprof).
      ! In this example we simulate all specified channels for each profile, but
      ! in general one can simulate a different number of channels for each profile.

      nprof_max_sim = nchan_max_sim / nchan_inst
      nprof_sim = min(nprof_max_sim, self % nprofiles)

      prof_start = 1
      prof_end = self % nprofiles

      !DARFIX should actually be packed count of skip_profiles * SIZE(channels)
      nprof_sim = min(nprof_sim, prof_end - prof_start + 1)
      nchan_sim = nprof_sim * size(self%channels)

      ! --------------------------------------------------------------------------
      ! Allocate RTTOV input and output structures
      ! --------------------------------------------------------------------------

      call self % RTprof_K % alloc(errorstatus, self % conf, nprof_sim, nchan_sim, self % nlevels, init=.true., asw=2)

      do while ( prof_start <= prof_end)

        ! --------------------------------------------------------------------------
        ! Build the list of profile/channel indices in chanprof
        ! --------------------------------------------------------------------------

        nchan_sim = 0_jpim
        nprof_sim = min(nprof_sim, prof_end - prof_start + 1)

        do iprof = 1, min(nprof_sim, prof_end - prof_start + 1)
          if(.not. self % Skip_Profiles(prof_start + iprof - 1)) then
            do ichan = 1, nchan_inst
              nchan_sim = nchan_sim + 1_jpim

              self % RTProf_K % chanprof(nchan_total + nchan_sim) % prof = iprof
              self % RTprof_K % chanprof(nchan_total + nchan_sim) % chan = self % channels(ichan) 
            end do
          endif
        end do

        !Assign the data from the GeoVaLs
        !--------------------------------

        call load_atm_data_rttov(geovals,obss,self % RTprof_K % profiles,prof_start,self % conf,layer_quantities)
        call load_geom_data_rttov(obss,self % RTprof_K % profiles,prof_start)

        ! --------------------------------------------------------------------------
        ! Set surface emissivity
        ! --------------------------------------------------------------------------
        call self % RTProf_K % init_emissivity(self % conf)

        ! Inintialize the K-matrix INPUT so that the results are dTb/dx
        ! -------------------------------------------------------------

        self % RTprof_K % emissivity_k(:) % emis_out = 0
        self % RTprof_K % emissivity_k(:) % emis_in = 0
        self % RTprof_K % emissivity(:) % emis_out = 0
        self % RTprof_K % radiance_k % bt(:) = 1
        self % RTprof_K % radiance_k % total(:) = 1

        if(self % conf % RTTOV_profile_checkinput) then
          if (self % conf % inspect > 0) then
            ! no error checking, could check multiple profiles in loop but what would be the point
            call rttov_print_profile(self % RTprof_K % profiles(self % conf % inspect))

            call rttov_user_profile_checkinput(rttov_errorstatus, &
              self % conf % rttov_opts, &
              self % conf % rttov_coef_array(i_inst), &
              self % RTprof_K % profiles(self % conf % inspect))
          endif
        endif

        ! --------------------------------------------------------------------------
        ! Call RTTOV K model
        ! --------------------------------------------------------------------------

        call rttov_k(                              &
          errorstatus,                             &! out   error flag
          self % RTprof_K % chanprof(nchan_total + 1:nchan_total + nchan_sim), &! in channel and profile index structure
          self % conf % rttov_opts,                     &! in    options structure
          self % RTprof_K % profiles,                                &! in    profile array
          self % RTprof_K % profiles_k(nchan_total + 1 : nchan_total + nchan_sim), &! in    profile array
          self % conf % rttov_coef_array(i_inst), &! in    coefficients structure
          self % RTprof_K % transmission,                            &! inout computed transmittances
          self % RTprof_K % transmission_k,                          &! inout computed transmittances
          self % RTprof_K % radiance,                                &! inout computed radiances
          self % RTprof_K % radiance_k,                              &! inout computed radiances
          calcemis    = self % RTprof_K % calcemis,                  &! in    flag for internal emissivity calcs
          emissivity  = self % RTprof_K % emissivity,                &!, &! inout input/output emissivities per channel
          emissivity_k = self % RTprof_K % emissivity_k)!,           &! inout input/output emissivities per channel      

        if (self % conf % inspect > 0) then
          ! no error checking, could check multiple channels here using inspect_k array (NOT IMPLEMENTED)
          call rttov_print_profile(self % RTprof_K % profiles_k(self % conf % inspect))
        endif

        if ( errorstatus /= errorstatus_success ) then
          write(message,'(A, A, 2I6)') trim(ROUTINE_NAME), 'after rttov_k: error ', errorstatus, i_inst
          call abor1_ftn(message)
        end if

        write(*,'(A1, i0, A, i0)',ADVANCE="NO") achar(13), prof_start+nprof_sim-1, ' locations processed out of ', geovals%nlocs

        ! Put simulated diagnostics into hofxdiags
        ! ----------------------------------------------
        if(hofxdiags%nvar > 0)     call populate_hofxdiags(self % RTprof_K, self % RTprof_K % chanprof, hofxdiags)

        prof_start = prof_start + nprof_sim
        nchan_total = nchan_total + nchan_sim

        self % nchan_total = nchan_total
      end do
      ! Deallocate structures for rttov_direct
    end do Sensor_Loop

    write(*,*)

    ! Set flag that the tracectory was set
    ! ------------------------------------
    self % ltraj = .true.

  end subroutine ufo_radiancerttov_tlad_settraj

  ! ------------------------------------------------------------------------------
  subroutine ufo_radiancerttov_simobs_tl(self, geovals, obss, nvars, nlocs, hofx)
    
    use ufo_constants_mod, only : zero, g_to_kg

    implicit none
  
  class(ufo_radiancerttov_tlad), intent(in)    :: self
    type(ufo_geovals),           intent(in)    :: geovals
    type(c_ptr), value,          intent(in)    :: obss
    integer,                     intent(in)    :: nvars, nlocs
    real(c_double),              intent(inout) :: hofx(nvars, nlocs)

    character(len=*), parameter                :: myname_="ufo_radiancerttov_simobs_tl"
    integer                                    :: ichan, jchan, prof, jspec

    type(ufo_geoval), pointer                  :: geoval_d, geoval_d2

    ! Initial checks
    ! --------------

    ! Check if trajectory was set
    if (.not. self % ltraj) then
      write(message,*) myname_, ' trajectory wasnt set!'
      call abor1_ftn(message)
    end if

    ! Check if nlocs is consistent in geovals & hofx
    if (geovals % nlocs /= self % nprofiles) then
      write(message,*) myname_, ' error: nlocs inconsistent!'
      call abor1_ftn(message)
    end if

    ! Initialize hofx
    ! ---------------
    hofx(:,:) = zero

    ! Temperature
    ! -----------
    call ufo_geovals_get_var(geovals, var_ts, geoval_d) ! var_ts = air_temperature

    ! Check model levels is consistent in geovals
    if (geoval_d % nval /= self % nlevels) then
      write(message,*) myname_, ' error: layers inconsistent!'
      call abor1_ftn(message)

    end if

    do ichan = 1, self % nchan_total, size(self%channels)
      prof = self % RTprof_K % chanprof(ichan) % prof
      if (.not. self % Skip_Profiles(prof)) then
        do jchan = 1, size(self%channels)
          hofx(jchan,prof) = hofx(jchan,prof) + &
            sum(self % RTprof_K % profiles_k(ichan+jchan-1) % t(self % nlevels:1:-1) * geoval_d % vals(1:geoval_d % nval,prof))
        enddo
      endif
    end do

    do jspec = 1, self%conf%ngas
      call ufo_geovals_get_var(geovals, self%conf%Absorbers(jspec), geoval_d)

      ! Check model levels is consistent in geovals
      if (geoval_d % nval /= self % nlevels) then
        write(message,*) myname_, ' error: layers inconsistent!'
        call abor1_ftn(message)
      end if

      ! Absorbers
      ! ---------
      ! This is where CO2 and friends will live as well as CLW
      do ichan = 1, self % nchan_total, size(self%channels)
        prof = self % RTprof_K % chanprof(ichan) % prof
        if (.not. self % Skip_Profiles(prof)) then
          do jchan = 1, size(self%channels)
            if(self%conf%Absorbers(jspec) == var_q) then
              hofx(jchan,prof) = hofx(jchan,prof) + &
                sum(self % RTprof_K % profiles_k(ichan+jchan-1) % q(self % nlevels:1:-1) * geoval_d % vals(1:geoval_d % nval,prof))
            elseif(self%conf%Absorbers(jspec) == var_mixr) then
              hofx(jchan,prof) = hofx(jchan,prof) + &
                sum(self % RTprof_K % profiles_k(ichan+jchan-1) % q(self % nlevels:1:-1) * geoval_d % vals(1:geoval_d % nval,prof)) / &
                g_to_kg
            elseif(self%conf%Absorbers(jspec) == var_clw) then
              hofx(jchan,prof) = hofx(jchan,prof) + &
                sum(self % RTprof_K % profiles_k(ichan+jchan-1) % clw(self % nlevels:1:-1) * geoval_d % vals(1:geoval_d % nval,prof))
            endif
          enddo
        endif
      end do
    enddo

    ! Cloud
    ! --------------------------
    !IR 


    ! Surface + Single-valued Variables
    ! --------------------------
    !CTP/cloudfrac
    !O3total
    !LWP

    !T2m
    call ufo_geovals_get_var(geovals, var_sfc_t2m, geoval_d)

    do ichan = 1, self % nchan_total, size(self%channels)
      prof = self % RTprof_K % chanprof(ichan) % prof
      if (.not. self % Skip_Profiles(prof)) then
        do jchan = 1, size(self%channels)
          hofx(jchan,prof) = hofx(jchan,prof) + &
            self % RTprof_K % profiles_k(ichan+jchan-1) % s2m % t * geoval_d % vals(1,prof)
        enddo
      endif
    end do

    !q2m
    call ufo_geovals_get_var(geovals, var_sfc_q2m, geoval_d)

    do ichan = 1, self % nchan_total, size(self%channels)
      prof = self % RTprof_K % chanprof(ichan) % prof
      if (.not. self % Skip_Profiles(prof)) then
        do jchan = 1, size(self%channels)
          hofx(jchan,prof) = hofx(jchan,prof) + &
            self % RTprof_K % profiles_k(ichan+jchan-1) % s2m % q * geoval_d % vals(1,prof)
        enddo
      endif
    end do

    !windspeed
    call ufo_geovals_get_var(geovals, var_u, geoval_d)
    call ufo_geovals_get_var(geovals, var_v, geoval_d2)

    do ichan = 1, self % nchan_total, size(self%channels)
      prof = self % RTprof_K % chanprof(ichan) % prof
      if (.not. self % Skip_Profiles(prof)) then
        do jchan = 1, size(self%channels)
          hofx(jchan,prof) = hofx(jchan,prof) + &
            self % RTprof_K % profiles_k(ichan+jchan-1) % s2m % u * geoval_d % vals(1,prof) + &
            self % RTprof_K % profiles_k(ichan+jchan-1) % s2m % v * geoval_d2 % vals(1,prof)
        enddo
      endif
    end do

    !Tskin
    call ufo_geovals_get_var(geovals, var_sfc_tskin, geoval_d)

    do ichan = 1, self % nchan_total, size(self%channels)
      prof = self % RTprof_K % chanprof(ichan) % prof
      if (.not. self % Skip_Profiles(prof)) then
        do jchan = 1, size(self%channels)
          hofx(jchan,prof) = hofx(jchan,prof) + &
            self % RTprof_K % profiles_k(ichan+jchan-1) % skin % t * geoval_d % vals(1,prof)
        enddo
      endif
    end do

  end subroutine ufo_radiancerttov_simobs_tl

  ! ------------------------------------------------------------------------------
  subroutine ufo_radiancerttov_simobs_ad(self, geovals, obss, nvars, nlocs, hofx)

    use ufo_constants_mod, only : zero, g_to_kg

    implicit none

    class(ufo_radiancerttov_tlad), intent(in)    :: self
    type(ufo_geovals),             intent(inout) :: geovals
    type(c_ptr), value,            intent(in)    :: obss
    integer,                       intent(in)    :: nvars, nlocs
    real(c_double),                intent(in)    :: hofx(nvars, nlocs)

    type(ufo_geoval), pointer                    :: geoval_d, geoval_d2

    real(c_double)                               :: missing
    integer                                      :: ichan, jchan, prof, jspec

    character(len=*), parameter                  :: myname_ = "ufo_radiancerttov_simobs_ad"

    ! Set missing value
    missing = missing_value(missing)

    ! Initial checks
    ! --------------

    ! Check if trajectory was set
    if (.not. self % ltraj) then
      write(message,*) myname_, ' trajectory wasnt set!'
      call abor1_ftn(message)
    end if

    ! Check if nlocs is consistent in geovals & hofx
    if (geovals % nlocs /= self % nprofiles) then
      write(message,*) myname_, ' error: nlocs inconsistent!'
      call abor1_ftn(message)
    end if

    ! Temperature
    ! -----------
    call ufo_geovals_get_var(geovals, var_ts, geoval_d) ! var_ts = air_temperature

    ! allocate if not yet allocated
    if (.not. allocated(geoval_d % vals)) then
      geoval_d % nlocs = self % nprofiles
      geoval_d % nval = self % nlevels
      allocate(geoval_d % vals(geoval_d % nval,geoval_d % nlocs))
      geoval_d % vals = zero
    end if

    do ichan = 1, self % nchan_total, size(self%channels)
      prof = self % RTprof_K % chanprof(ichan) % prof
      if (.not. self % Skip_Profiles(prof)) then
        do jchan = 1, size(self%channels)
          if (hofx(jchan, prof) /= missing) then
            geoval_d % vals(:,prof) = geoval_d % vals(:,prof) + &
              self % RTprof_K % profiles_k(ichan+jchan-1) % t(self % nlevels:1:-1) * hofx(jchan,prof)
          endif
        enddo
      endif
    end do
    
    ! Absorbers
    ! ---------
    ! This is where CO2 and friends will live as well as CLW

    do jspec = 1, self%conf%ngas
      call ufo_geovals_get_var(geovals, self%conf%Absorbers(jspec), geoval_d)

      ! allocate if not yet allocated
      if (.not. allocated(geoval_d % vals)) then
        geoval_d % nlocs = self % nprofiles
        geoval_d % nval = self % nlevels
        allocate(geoval_d % vals(geoval_d % nval,geoval_d % nlocs))
        geoval_d % vals = zero
      end if

      do ichan = 1, self % nchan_total, size(self%channels)
        prof = self % RTprof_K % chanprof(ichan) % prof
        if (.not. self % Skip_Profiles(prof)) then
          do jchan = 1, size(self%channels)
            if (hofx(jchan, prof) /= missing) then
              
              if(self%conf%Absorbers(jspec) == var_q) then
                geoval_d % vals(:,prof) = geoval_d % vals(:,prof) + &
                  self % RTprof_K % profiles_k(ichan+jchan-1) % q(self % nlevels:1:-1) * hofx(jchan,prof)
              elseif(self%conf%Absorbers(jspec) == var_mixr) then
                geoval_d % vals(:,prof) = geoval_d % vals(:,prof) + &
                  (self % RTprof_K % profiles_k(ichan+jchan-1) % q(self % nlevels:1:-1) / g_to_kg) * hofx(jchan,prof)
              elseif(self%conf%Absorbers(jspec) == var_clw) then
                geoval_d % vals(:,prof) = geoval_d % vals(:,prof) + &
                  self % RTprof_K % profiles_k(ichan+jchan-1) % clw(self % nlevels:1:-1) * hofx(jchan,prof)
              endif
            endif
          enddo
        endif
      enddo
    enddo

    ! Cloud
    ! --------------------------
    !IR 

    ! Surface + Single-valued Variables
    ! --------------------------
    !CTP/cloudfrac
    !O3total
    !LWP
    !T2m

    call ufo_geovals_get_var(geovals, var_sfc_t2m, geoval_d) 

    ! allocate if not yet allocated
    if (.not. allocated(geoval_d % vals)) then
      geoval_d % nlocs = self % nprofiles
      geoval_d % nval = self % nlevels
      allocate(geoval_d % vals(geoval_d % nval,geoval_d % nlocs)) ! DARFIX try setting to 1?
      geoval_d % vals = zero
    end if

    do ichan = 1, self % nchan_total, size(self%channels)
      prof = self % RTprof_K % chanprof(ichan) % prof
      if (.not. self % Skip_Profiles(prof)) then
        do jchan = 1, size(self%channels)
          if (hofx(jchan, prof) /= missing) then
            geoval_d % vals(1,prof) = geoval_d % vals(1,prof) + &
              self % RTprof_K % profiles_k(ichan+jchan-1) % s2m % t * hofx(jchan,prof)
          endif
        enddo
      endif
    end do

    !q2m
    call ufo_geovals_get_var(geovals, var_sfc_q2m, geoval_d) 
    ! allocate if not yet allocated
    if (.not. allocated(geoval_d % vals)) then
      geoval_d % nlocs = self % nprofiles
      geoval_d % nval = self % nlevels
      allocate(geoval_d % vals(geoval_d % nval,geoval_d % nlocs)) ! DARFIX try setting to 1?
      geoval_d % vals = zero
    end if

    do ichan = 1, self % nchan_total, size(self%channels)
      prof = self % RTprof_K % chanprof(ichan) % prof
      if (.not. self % Skip_Profiles(prof)) then
        do jchan = 1, size(self%channels)
          if (hofx(jchan, prof) /= missing) then
            geoval_d % vals(1,prof) = geoval_d % vals(1,prof) + &
              self % RTprof_K % profiles_k(ichan+jchan-1) % s2m % q * hofx(jchan,prof)
          endif
        enddo
      endif
    end do
      
    !windspeed
    call ufo_geovals_get_var(geovals, var_u, geoval_d)
    call ufo_geovals_get_var(geovals, var_v, geoval_d2)

    ! allocate if not yet allocated
    if (.not. allocated(geoval_d % vals)) then
      geoval_d % nlocs = self % nprofiles
      geoval_d % nval = self % nlevels
      geoval_d2 % nlocs = self % nprofiles
      geoval_d2 % nval = self % nlevels
      allocate(geoval_d % vals(geoval_d % nval,geoval_d % nlocs), &
        geoval_d2 % vals(geoval_d % nval,geoval_d % nlocs)) ! DARFIX try setting to 1?
      geoval_d % vals = zero
      geoval_d2 % vals = zero
    end if

    do ichan = 1, self % nchan_total, size(self%channels)
      prof = self % RTprof_K % chanprof(ichan) % prof
      if (.not. self % Skip_Profiles(prof)) then
        do jchan = 1, size(self%channels)
          if (hofx(jchan, prof) /= missing) then
            geoval_d % vals(1,prof) = geoval_d % vals(1,prof) + &
              self % RTprof_K % profiles_k(ichan+jchan-1) % s2m % u * hofx(jchan,prof)
            
            geoval_d2 % vals(1,prof) = geoval_d2 % vals(1,prof) + &
              self % RTprof_K % profiles_k(ichan+jchan-1) % s2m % v * hofx(jchan,prof)
          endif
        enddo
      endif
    end do

    !Tskin
    call ufo_geovals_get_var(geovals, var_sfc_tskin, geoval_d)

    ! allocate if not yet allocated
    if (.not. allocated(geoval_d % vals)) then
      geoval_d % nlocs = self % nprofiles
      geoval_d % nval = self % nlevels
      allocate(geoval_d % vals(geoval_d % nval,geoval_d % nlocs)) ! DARFIX try setting to 1?
      geoval_d % vals = zero
    end if

    do ichan = 1, self % nchan_total, size(self%channels)
      prof = self % RTprof_K % chanprof(ichan) % prof
      if (.not. self % Skip_Profiles(prof)) then
        do jchan = 1, size(self%channels)
          if (hofx(jchan, prof) /= missing) then
            geoval_d % vals(1,prof) = geoval_d % vals(1,prof) + &
              self % RTprof_K % profiles_k(ichan+jchan-1) % skin % t * hofx(jchan,prof)
          endif
        enddo
      endif
    end do
      
    ! Once all geovals set replace flag
    ! ---------------------------------
    if (.not. geovals % linit ) geovals % linit=.true.

  end subroutine ufo_radiancerttov_simobs_ad

  ! ------------------------------------------------------------------------------

end module ufo_radiancerttov_tlad_mod
