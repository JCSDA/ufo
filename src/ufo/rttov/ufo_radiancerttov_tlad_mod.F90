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
  use rttov_const ! errorstatus and gas_id
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
    type(rttov_chanprof), allocatable            :: chanprof(:)
    type(ufo_geoval), pointer                    :: geoval_temp

    integer(kind=jpim)                           :: errorstatus ! Return error status of RTTOV subroutine calls

    integer                                      :: i_inst, nchan_total, ichan, iprof, prof, iprof_rttov
    integer                                      :: nprof_sim, nprof_max_sim, ichan_sim
    integer                                      :: prof_start, prof_end

    logical                                      :: jacobian_needed

    include 'rttov_k.interface'

    !Initialisations
    missing = missing_value(missing)

    !Return the name and name length of obsspace communicator (from ioda)
    call obsspace_get_comm(obss, f_comm)
    
    !! Parse hofxdiags%variables into independent/dependent variables and channel assumed formats:
    ! Note this sets the jacobian_needed flag
    !!   jacobian var -->     <ystr>_jacobian_<xstr>_<chstr>
    !!   non-jacobian var --> <ystr>_<chstr>
    call parse_hofxdiags(hofxdiags, jacobian_needed)

    ! Get number of profiles and levels from geovals
    self % nprofiles = geovals % nlocs
    call ufo_geovals_get_var(geovals, var_ts, geoval_temp)
    self % nlevels = geoval_temp % nval
    nullify(geoval_temp)

    ! Allocate RTTOV profiles for ALL geovals for the direct calculation
    write(message,'(A, A, I0, A, I0, A)') &
      trim(routine_name), ': Allocating ', self % nprofiles, ' profiles with ', self % nlevels, ' levels'
    call fckit_log%debug(message)
    call self % RTprof_K % alloc_profs(errorstatus, self % conf, self % nprofiles, self % nlevels, init=.true., asw=1)

    !Assign the atmospheric and surface data from the GeoVaLs
    write(message,'(A, A, I0, A, I0, A)') trim(routine_name), ': Creating RTTOV profiles from geovals'
    call fckit_log%debug(message)
    call self % RTprof_K % setup(geovals,obss,self % conf)

    !DAR: Removing sensor_loop until it's demonstrated to be needed and properly tested
    ! at the moment self % channels is a single 1D array so cannot adequately contain more than one set of channels
!    Sensor_Loop:do i_inst = 1, self % conf % nSensors
    i_inst = 1

    ! Number of channels to be simulated for this instrument (from the configuration, not necessarily the full instrument complement)
    nchan_inst = size(self % channels)

    ! Allocate memory for *ALL* RTTOV_K channels
    write(message,'(A,A,I0,A)') &
      trim(routine_name), ': Allocating Trajectory resources for RTTOV K: ', self % nprofiles * nchan_inst, ' total channels'
    call self % RTprof_K % alloc_profs_K(errorstatus, self % conf, self % nprofiles * nchan_inst, self % nlevels, init=.true., asw=1)

    ! Used for keeping track of profiles for setting emissivity
    allocate(self % RTprof_K % chanprof ( self % nprofiles * nchan_inst )) 

    ! Maximum number of profiles to be processed by RTTOV per pass
    if(self % conf % prof_by_prof) then
      nprof_max_sim = 1
    else
      nprof_max_sim = max(1,self % conf % nchan_max_sim / nchan_inst)
    endif
    nprof_sim = min(nprof_max_sim, self % nprofiles)

    ! Determine the total number of radiances to simulate (nchan_sim).
    nchan_sim = nprof_sim * size(self%channels)

    ! Allocate structures for RTTOV direct and K code
    write(message,'(A,A,I0,A,I0,A)') &
      trim(routine_name), ': Allocating resources for RTTOV direct (K): ', nprof_sim, ' and ', nchan_sim, ' channels'
    call fckit_log%debug(message)
    call self % RTprof_K % alloc_direct(errorstatus, self % conf, nprof_sim, nchan_sim, self % nlevels, init=.true., asw=1)

    write(message,'(A,A,I0,A,I0,A)') &
      trim(routine_name), ': Allocating resources for RTTOV K code: ', nprof_sim, ' and ', nchan_sim, ' channels'
    call fckit_log%debug(message)
    call self % RTprof_K % alloc_k(errorstatus, self % conf, nprof_sim, nchan_sim, self % nlevels, init=.true., asw=1)

    prof_start = 1; prof_end = self % nprofiles
    nchan_total = 0

    RTTOV_loop : do while (prof_start <= prof_end)

      ! Reduce number of simulated profiles/channel if at end of the of profiles to be processed
      nprof_sim = min(nprof_sim, prof_end - prof_start + 1)
      nchan_sim = nprof_sim * size(self%channels)

      ! allocate and initialise local chanprof structure
      allocate(chanprof ( nchan_sim ))
      chanprof(:) % prof = 0
      chanprof(:) % chan = 0

      ! index for simulated channel
      ichan_sim = 0_jpim

      ! Build the list of profile/channel indices in chanprof
      do iprof_rttov = 1, nprof_sim

        ! iprof is the index for the full set of RTTOV profiles
        iprof = prof_start + iprof_rttov - 1

        do ichan = 1, nchan_inst
          ichan_sim = ichan_sim + 1_jpim
          chanprof(ichan_sim) % prof = iprof_rttov ! this refers to the slice of the RTprofile array passed to RTTOV
          chanprof(ichan_sim) % chan = self % channels(ichan)
          self % RTprof_K % chanprof(nchan_total + ichan_sim) % prof = iprof
          self % RTprof_K % chanprof(nchan_total + ichan_sim) % chan = self % channels(ichan)
        end do
        nchan_sim = ichan_sim

        if(self % conf % RTTOV_profile_checkinput) call self % RTprof_K % check(self % conf, iprof, i_inst)
        if(any(self % conf % inspect == iprof)) call self % RTprof_K % print(self % conf, iprof, i_inst)
            
      end do

      ! Set surface emissivity
      call self % RTProf_K % init_emissivity(self % conf, prof_start)

      ! --------------------------------------------------------------------------
      ! Call RTTOV K model
      ! --------------------------------------------------------------------------
    
      call rttov_k(                              &
        errorstatus,                             &! out   error flag
        chanprof(1:nchan_sim), &! in channel and profile index structure
        self % conf % rttov_opts,                     &! in    options structure
        self % RTprof_K % profiles(prof_start:prof_start + nprof_sim - 1), &! in    profile array
        self % RTprof_K % profiles_k(nchan_total + 1 : nchan_total + nchan_sim), &! in    profile array
        self % conf % rttov_coef_array(i_inst), &! in    coefficients structure
        self % RTprof_K % transmission,                            &! inout computed transmittances
        self % RTprof_K % transmission_k,                          &! inout computed transmittances
        self % RTprof_K % radiance,                                &! inout computed radiances
        self % RTprof_K % radiance_k,                              &! inout computed radiances
        calcemis    = self % RTprof_K % calcemis(1:nchan_sim),                  &! in    flag for internal emissivity calcs
        emissivity  = self % RTprof_K % emissivity(1:nchan_sim),                &!, &! inout input/output emissivities per channel
        emissivity_k = self % RTprof_K % emissivity_k(1:nchan_sim))!,           &! inout input/output emissivities per channel      
      
      if ( errorstatus /= errorstatus_success ) then
        write(message,'(A, A, 2I6)') trim(routine_name), 'after rttov_k: error ', errorstatus, i_inst
        call abor1_ftn(message)
      end if
      
      ! Put simulated diagnostics into hofxdiags
      ! ----------------------------------------------
      if(hofxdiags%nvar > 0)     call populate_hofxdiags(self % RTprof_K, chanprof, self % conf, hofxdiags)

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
    call self % RTprof_K % alloc_profs(errorstatus, self % conf, -1, -1, asw=0)
    
 
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
      do jchan = 1, size(self%channels)
        hofx(jchan,prof) = hofx(jchan,prof) + &
          sum(self % RTprof_K % profiles_k(ichan+jchan-1) % t(self % nlevels:1:-1) * geoval_d % vals(1:geoval_d % nval,prof))
      enddo
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
          do jchan = 1, size(self%channels)
            if(self%conf%Absorbers(jspec) == var_q) then
              hofx(jchan,prof) = hofx(jchan,prof) + &
                sum(self % RTprof_K % profiles_k(ichan+jchan-1) % q(self % nlevels:1:-1) * &
                    geoval_d % vals(1:geoval_d % nval,prof)) * self%conf%scale_fac(gas_id_watervapour)
            elseif(self%conf%Absorbers(jspec) == var_mixr) then
              hofx(jchan,prof) = hofx(jchan,prof) + &
                sum(self % RTprof_K % profiles_k(ichan+jchan-1) % q(self % nlevels:1:-1) * &
                    geoval_d % vals(1:geoval_d % nval,prof)) * self%conf%scale_fac(gas_id_watervapour) / g_to_kg
            elseif(self%conf%Absorbers(jspec) == var_clw) then
              hofx(jchan,prof) = hofx(jchan,prof) + &
                sum(self % RTprof_K % profiles_k(ichan+jchan-1) % clw(self % nlevels:1:-1) * &
                    geoval_d % vals(1:geoval_d % nval,prof))
            endif
          enddo
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
      do jchan = 1, size(self%channels)
        hofx(jchan,prof) = hofx(jchan,prof) + &
          self % RTprof_K % profiles_k(ichan+jchan-1) % s2m % t * geoval_d % vals(1,prof)
      enddo
    end do

    !q2m
    call ufo_geovals_get_var(geovals, var_sfc_q2m, geoval_d)

    do ichan = 1, self % nchan_total, size(self%channels)
      prof = self % RTprof_K % chanprof(ichan) % prof
      do jchan = 1, size(self%channels)
        hofx(jchan,prof) = hofx(jchan,prof) + &
          self % RTprof_K % profiles_k(ichan+jchan-1) % s2m % q * &
          geoval_d % vals(1,prof) * self%conf%scale_fac(gas_id_watervapour)
      enddo
    end do

    !windspeed
    call ufo_geovals_get_var(geovals, var_u, geoval_d)
    call ufo_geovals_get_var(geovals, var_v, geoval_d2)

    do ichan = 1, self % nchan_total, size(self%channels)
      prof = self % RTprof_K % chanprof(ichan) % prof
      do jchan = 1, size(self%channels)
        hofx(jchan,prof) = hofx(jchan,prof) + &
          self % RTprof_K % profiles_k(ichan+jchan-1) % s2m % u * geoval_d % vals(1,prof) + &
          self % RTprof_K % profiles_k(ichan+jchan-1) % s2m % v * geoval_d2 % vals(1,prof)
      enddo
    end do

    !Tskin
    call ufo_geovals_get_var(geovals, var_sfc_tskin, geoval_d)

    do ichan = 1, self % nchan_total, size(self%channels)
      prof = self % RTprof_K % chanprof(ichan) % prof
      do jchan = 1, size(self%channels)
        hofx(jchan,prof) = hofx(jchan,prof) + &
          self % RTprof_K % profiles_k(ichan+jchan-1) % skin % t * geoval_d % vals(1,prof)
      enddo
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
      do jchan = 1, size(self%channels)
        if (hofx(jchan, prof) /= missing) then
          geoval_d % vals(:,prof) = geoval_d % vals(:,prof) + &
            self % RTprof_K % profiles_k(ichan+jchan-1) % t(self % nlevels:1:-1) * hofx(jchan,prof)
        endif
      enddo
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
        do jchan = 1, size(self%channels)
          if (hofx(jchan, prof) /= missing) then
            
            if(self%conf%Absorbers(jspec) == var_q) then
              geoval_d % vals(:,prof) = geoval_d % vals(:,prof) + &
                self % RTprof_K % profiles_k(ichan+jchan-1) % q(self % nlevels:1:-1) * hofx(jchan,prof) * &
                self%conf%scale_fac(gas_id_watervapour)
            elseif(self%conf%Absorbers(jspec) == var_mixr) then
              geoval_d % vals(:,prof) = geoval_d % vals(:,prof) + &
                (self % RTprof_K % profiles_k(ichan+jchan-1) % q(self % nlevels:1:-1) * hofx(jchan,prof)) * &
                self%conf%scale_fac(gas_id_watervapour) / g_to_kg
            elseif(self%conf%Absorbers(jspec) == var_clw) then
              geoval_d % vals(:,prof) = geoval_d % vals(:,prof) + &
                self % RTprof_K % profiles_k(ichan+jchan-1) % clw(self % nlevels:1:-1) * hofx(jchan,prof)
            endif
          endif
        enddo
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
      do jchan = 1, size(self%channels)
        if (hofx(jchan, prof) /= missing) then
          geoval_d % vals(1,prof) = geoval_d % vals(1,prof) + &
            self % RTprof_K % profiles_k(ichan+jchan-1) % s2m % t * hofx(jchan,prof)
        endif
      enddo
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
      do jchan = 1, size(self%channels)
        if (hofx(jchan, prof) /= missing) then
          geoval_d % vals(1,prof) = geoval_d % vals(1,prof) + &
            self % RTprof_K % profiles_k(ichan+jchan-1) % s2m % q * &
            hofx(jchan,prof) * self%conf%scale_fac(gas_id_watervapour)
        endif
      enddo
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
      do jchan = 1, size(self%channels)
        if (hofx(jchan, prof) /= missing) then
          geoval_d % vals(1,prof) = geoval_d % vals(1,prof) + &
            self % RTprof_K % profiles_k(ichan+jchan-1) % s2m % u * hofx(jchan,prof)
          
          geoval_d2 % vals(1,prof) = geoval_d2 % vals(1,prof) + &
            self % RTprof_K % profiles_k(ichan+jchan-1) % s2m % v * hofx(jchan,prof)
        endif
      enddo
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
      do jchan = 1, size(self%channels)
        if (hofx(jchan, prof) /= missing) then
          geoval_d % vals(1,prof) = geoval_d % vals(1,prof) + &
            self % RTprof_K % profiles_k(ichan+jchan-1) % skin % t * hofx(jchan,prof)
        endif
      enddo
    end do
      
    ! Once all geovals set replace flag
    ! ---------------------------------
    if (.not. geovals % linit ) geovals % linit=.true.

  end subroutine ufo_radiancerttov_simobs_ad

  ! ------------------------------------------------------------------------------

end module ufo_radiancerttov_tlad_mod
