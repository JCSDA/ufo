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
use ufo_rttovonedvarcheck_bmatrix_mod
use ufo_rttovonedvarcheck_constants_mod
use ufo_rttovonedvarcheck_minimize_utils_mod
use ufo_rttovonedvarcheck_minimize_newton_mod
use ufo_rttovonedvarcheck_minimize_ml_mod
use ufo_rttovonedvarcheck_obinfo_mod
use ufo_rttovonedvarcheck_profindex_mod
use ufo_rttovonedvarcheck_rmatrix_mod
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
!! \author M. Cooke (Met Office)
!!
!! \date 09/06/2020: Created
!!
subroutine ufo_rttovonedvarcheck_create(self, obspace, f_conf, channels, &
                                        onedvarflag)

  implicit none
  type(ufo_rttovonedvarcheck), intent(inout) :: self         !< rttovonedvarcheck main object
  type(c_ptr), value, intent(in)             :: obspace      !< observation database pointer
  type(fckit_configuration), intent(in)      :: f_conf       !< yaml file contents
  integer(c_int), intent(in)                 :: channels(:)  !< all channels that can be used in 1D-Var
  integer(c_int), intent(in)                 :: onedvarflag  !< flag for qc manager

  self % obsdb = obspace
  self % conf = f_conf
  self % onedvarflag = onedvarflag

  call ufo_rttovonedvarcheck_setup(self, channels) ! from init

end subroutine ufo_rttovonedvarcheck_create

! ------------------------------------------------------------------------------
!> Delete the main rttov onedvar object in Fortran
!!
!! \author M. Cooke (Met Office)
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
!! \author M. Cooke (Met Office)
!!
!! \date 09/06/2020: Created
!!
subroutine ufo_rttovonedvarcheck_apply(self, vars, geovals, apply)

  implicit none
  type(ufo_rttovonedvarcheck), intent(inout) :: self     !< rttovonedvarcheck main object
  type(oops_variables), intent(in)           :: vars     !< channels for 1D-Var
  type(ufo_geovals), intent(in)              :: geovals  !< model values at observation space
  logical, intent(in)                        :: apply(:) !< qc manager flags

  type(bmatrix_type)   :: full_bmatrix
  type(ufo_geovals)    :: local_geovals  ! geoval for one observation
  type(obinfo_type)    :: ob_info
  type(profindex_type) :: prof_index     ! index for mapping geovals to 1d-var state profile
  type(rmatrix_type)   :: r_submatrix    ! r_submatrix object
  character(len=max_string)          :: sensor_id
  character(len=max_string)          :: var
  character(len=max_string)          :: varname
  character(len=max_string)          :: message
  integer                            :: iloc, jvar, jobs, ivar, band ! counters
  integer                            :: ii
  integer                            :: chans_used      ! counter for number of channels used for an ob
  integer                            :: jchans_used
  integer                            :: fileunit        ! unit number for reading in files
  integer                            :: apply_count
  integer                            :: nprofelements   ! number of elements in 1d-var state profile
  integer(c_int32_t), allocatable    :: flags(:,:)      ! qc flag for return to var obs file
  integer, allocatable               :: channels_used(:)
  integer, allocatable               :: fields_in(:)
  integer, allocatable               :: QCflags(:,:)    ! current qc flags needed for channel selection
  real(kind_real)                    :: missing         ! missing value
  real(kind_real)                    :: t1, t2          ! timing
  real(kind_real), allocatable       :: b_matrix(:,:)   ! 1d-var profile b matrix
  real(kind_real), allocatable       :: b_inverse(:,:)  ! inverse for each 1d-var profile b matrix
  real(kind_real), allocatable       :: b_sigma(:)      ! b_matrix diagonal error
  real(kind_real), allocatable       :: obs_error(:)
  real(kind_real), allocatable       :: yobs(:,:)       ! observation value from obs files
  real(kind_real), allocatable       :: yerr(:,:)       ! observation error from obs files
  real(kind_real), allocatable       :: lat(:)          ! observation latitude
  real(kind_real), allocatable       :: lon(:)          ! observation longitude
  real(kind_real), allocatable       :: elevation(:)    ! observation elevation
  real(kind_real), allocatable       :: sat_zen(:)      ! observation satellite zenith angle
  real(kind_real), allocatable       :: sat_azi(:)      ! observation satellite azimuth angle
  real(kind_real), allocatable       :: sol_zen(:)      ! observation solar zenith angle
  real(kind_real), allocatable       :: sol_azi(:)      ! observation solar azimuth angle
  real(kind_real), allocatable       :: emissivity(:,:) ! initial surface emissivity
  logical                            :: file_exists     ! check if a file exists logical
  logical                            :: onedvar_success
  logical                            :: cloud_retrieval = .false.

  ! ------------------------------------------
  ! Setup
  ! ------------------------------------------
  missing = missing_value(missing)
  iloc = obsspace_get_nlocs(self%obsdb)

  ! allocate arrays
  allocate(yobs(self%nchans,iloc))
  allocate(yerr(self%nchans,iloc))
  allocate(QCflags(self%nchans,iloc))
  allocate(lat(iloc))
  allocate(lon(iloc))
  allocate(elevation(iloc))
  allocate(sat_zen(iloc))
  allocate(sat_azi(iloc))
  allocate(sol_zen(iloc))
  allocate(sol_azi(iloc))
  allocate(flags(self%nchans,iloc))

  ! initialize arrays
  yobs(:,:) = 0.0
  yerr(:,:) = 0.0
  lat(:) = 0.0
  lon(:) = 0.0
  elevation(:) = 0.0
  sat_zen(:) = 0.0
  sat_azi(:) = 0.0
  sol_zen(:) = 0.0
  sol_azi(:) = 0.0
  flags(:,:) = 0

  ! Setup optional arrays
  if (self % ReadMWemiss .or. self % ReadIRemiss) then
    allocate(emissivity(self%nchans,iloc))
    emissivity(:,:) = 0.0
  end if

  ! read in observations and associated errors for full ObsSpace
  do jvar = 1, self%nchans
    var = vars%variable(jvar)
    call obsspace_get_db(self%obsdb, "ObsValue",  trim(var), yobs(jvar,:))
    call obsspace_get_db(self%obsdb, "ObsError",  trim(var), yerr(jvar,:))
    call obsspace_get_db(self%obsdb, "FortranQC", trim(var), QCflags(jvar,:))
    if (self % ReadMWemiss) then
      call obsspace_get_db(self%obsdb, "MwEmiss",   trim(var), emissivity(jvar,:))
    else if (self % ReadIRemiss) then
      call obsspace_get_db(self%obsdb, "IREmiss",   trim(var), emissivity(jvar,:))
    end if
  end do

  call obsspace_get_db(self%obsdb, "MetaData", "latitude", lat(:))
  call obsspace_get_db(self%obsdb, "MetaData", "longitude", lon(:))
  call obsspace_get_db(self%obsdb, "MetaData", "elevation", elevation(:))
  call obsspace_get_db(self%obsdb, "MetaData", "sensor_zenith_angle", sat_zen(:))
  call obsspace_get_db(self%obsdb, "MetaData", "sensor_azimuth_angle", sat_azi(:))
  call obsspace_get_db(self%obsdb, "MetaData", "solar_zenith_angle", sol_zen(:))
  call obsspace_get_db(self%obsdb, "MetaData", "solar_azimuth_angle", sol_azi(:))

  ! Setup full B matrix object
  call full_bmatrix % setup(self % retrieval_variables, self % b_matrix_path, self % qtotal)
  
  ! Check if cloud retrievals needed
  do ii = 1, size(self % retrieval_variables)
    if (trim(self % retrieval_variables(ii)) == "cloud_top_pressure") then
      write(*,*) "Simple cloud is part of the state vector"
      cloud_retrieval = .true.
    end if
  end do

  ! Create profile index for mapping 1d-var profile to b-matrix
  call prof_index % setup(full_bmatrix)

  ! Initialize data arrays
  allocate(b_matrix(prof_index % nprofelements,prof_index % nprofelements))
  allocate(b_inverse(prof_index % nprofelements,prof_index % nprofelements))
  allocate(b_sigma(prof_index % nprofelements))

  ! ------------------------------------------
  ! Beginning mains observations loop
  ! ------------------------------------------
  write(*,*) "Beginning observations loop: ",self%qcname
  apply_count = 0
  obs_loop: do jobs = 1, iloc
    if (apply(jobs)) then

      apply_count = apply_count + 1
      !---------------------------------------------------
      ! Setup Jb terms
      !---------------------------------------------------
      ! create one ob geovals from full all obs geovals
      call ufo_geovals_copy_one(geovals, local_geovals, jobs)
      call ufo_rttovonedvarcheck_check_geovals(local_geovals, prof_index)
      if (self % FullDiagnostics) call ufo_geovals_print(local_geovals, 1)

      ! select appropriate b matrix for latitude of observation
      b_matrix(:,:) = 0.0
      b_inverse(:,:) = 0.0
      b_sigma(:) = 0.0
      do band = 1, full_bmatrix % nbands
        if (lat(jobs) <  full_bmatrix % north(band)) exit
      end do
      ! Call to Ops_SatRad_ResetCovariances
      b_matrix(:,:) = full_bmatrix % store(:,:,band)
      b_inverse(:,:) = full_bmatrix % inverse(:,:,band)
      b_sigma(:) = full_bmatrix % sigma(:,band)

      !---------------------------------------------------
      ! Setup Jo terms
      !---------------------------------------------------
      ! Channel selection based on previous filter flags
      chans_used = 0
      do jvar = 1, self%nchans
        if( QCflags(jvar,jobs) == 0 ) then
          chans_used = chans_used + 1
        end if
      end do
      if (chans_used == 0) then
        write(message, *) "No channels selected for observation number ", &
                    jobs, " : skipping"
        call fckit_log % info(message)
        cycle obs_loop
      end if

      ! setup ob data for this observation
      call ob_info % setup(chans_used)
      call ob_info % init_emiss(self, local_geovals)
      ob_info % forward_mod_name = self % forward_mod_name
      ob_info % latitude = lat(jobs)
      ob_info % longitude = lon(jobs)
      ob_info % elevation = elevation(jobs)
      ob_info % sensor_zenith_angle = sat_zen(jobs)
      ob_info % sensor_azimuth_angle = sat_azi(jobs)
      ob_info % solar_zenith_angle = sol_zen(jobs)
      ob_info % solar_azimuth_angle = sol_azi(jobs)
      ob_info % retrievecloud = cloud_retrieval

      ! allocate arrays
      allocate(channels_used(chans_used))
      allocate(obs_error(chans_used))

      ! create obs vector and r matrix
      ! if emissivity read in use that
      jchans_used = 0
      do jvar = 1, self%nchans
        if( QCflags(jvar,jobs) == 0 ) then
          jchans_used = jchans_used + 1
          obs_error(jchans_used) = yerr(jvar, jobs)
          write(*,*) "jchans_used, ob number = ",trim(vars%variable(jvar)),jobs
          ob_info % yobs(jchans_used) = yobs(jvar, jobs)
          channels_used(jchans_used) = self%channels(jvar)
          if (self % ReadMWemiss .OR. self % ReadIRemiss) then
            write(*,*) "Copying emissivity from db"
            ob_info % emiss(jchans_used) = emissivity(jvar, jobs)
            where(ob_info % emiss == 0.0_kind_real) 
              ob_info % calc_emiss = .true.
            else where
              ob_info % calc_emiss = .false.
            end where
          end if
        end if
      end do

      write(*,*) "ob_info % emiss = ",ob_info % emiss
      write(*,*) "ob_info % calc_emiss = ",ob_info % calc_emiss

      call r_submatrix % setup(self % rtype, chans_used, obs_error)
      if (self % FullDiagnostics) then
        call r_submatrix % info()
        write(*, *) "Observations used = ",ob_info % yobs(:)
      end if
      write(*,*) "Channels used = ",channels_used

      !---------------------------------------------------
      ! Call minimization
      !---------------------------------------------------
      if (self % UseMLMinimization) then
        call ufo_rttovonedvarcheck_minimize_ml(self, ob_info, &
                                      r_submatrix, b_matrix, b_inverse, b_sigma, &
                                      local_geovals, prof_index,           &
                                      channels_used, onedvar_success)
      else
        call ufo_rttovonedvarcheck_minimize_newton(self, ob_info, &
                                      r_submatrix, b_matrix, b_inverse, b_sigma, &
                                      local_geovals, prof_index,           &
                                      channels_used, onedvar_success)
      end if

      ! Set QCflags based on output from minimization
      if (.NOT. onedvar_success) then
        do jvar = 1, self%nchans
          if( QCflags(jvar,jobs) == 0 ) then
            QCflags(jvar,jobs) = self % onedvarflag
          end if
        end do
      end if

      ! Tidy up memory specific to a single observation
      call ob_info % delete()
      call r_submatrix % delete()
      call ufo_geovals_delete(local_geovals)
      if (allocated(channels_used)) deallocate(channels_used)
      if (allocated(obs_error))     deallocate(obs_error)

    else
      call fckit_log % info("Final 1Dvar cost, apply = F")

    endif
  end do obs_loop

  write(message, *) "Total number of observations = ",iloc
  call fckit_log % info(message)
  write(message, *) "Number tested by 1dvar = ",apply_count
  call fckit_log % info(message)

  ! Put QC flags back in database
  do jvar = 1, self%nchans
    var = vars%variable(jvar)
    call obsspace_put_db(self%obsdb, "FortranQC", trim(var), QCflags(jvar,:))
  end do

  ! Tidy up memory used for all observations
  call full_bmatrix % delete()
  if (allocated(yobs))       deallocate(yobs)
  if (allocated(yerr))       deallocate(yerr)
  if (allocated(lat))        deallocate(lat)
  if (allocated(lon))        deallocate(lon)
  if (allocated(elevation))  deallocate(elevation)
  if (allocated(sat_zen))    deallocate(sat_zen)
  if (allocated(sat_azi))    deallocate(sat_azi)
  if (allocated(sol_zen))    deallocate(sol_zen)
  if (allocated(sol_azi))    deallocate(sol_azi)
  if (allocated(flags))      deallocate(flags)
  if (allocated(b_matrix))   deallocate(b_matrix)
  if (allocated(b_inverse))  deallocate(b_inverse)
  if (allocated(b_sigma))    deallocate(b_sigma)
  if (allocated(emissivity)) deallocate(emissivity)

end subroutine ufo_rttovonedvarcheck_apply

end module ufo_rttovonedvarcheck_mod
