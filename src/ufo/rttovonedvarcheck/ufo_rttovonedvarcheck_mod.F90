! (C) Copyright 2017-2020 Met Office
! 
! this software is licensed under the terms of the apache licence version 2.0
! which can be obtained at http://www.apache.org/licenses/license-2.0.

!> fortran module to implement onedvar fortran check

module ufo_rttovonedvarcheck_mod

use fckit_configuration_module, only: fckit_configuration
use iso_c_binding
use kinds
use ufo_geovals_mod
use ufo_vars_mod
use oops_variables_mod
use obsspace_mod
use ufo_rttovonedvarcheck_utils_mod
use ufo_rttovonedvarcheck_minimize_utils_mod

implicit none
private
public :: ufo_rttovonedvarcheck
public :: ufo_rttovonedvarcheck_create
public :: ufo_rttovonedvarcheck_delete
public :: ufo_rttovonedvarcheck_apply

! ------------------------------------------------------------------------------
contains
! ------------------------------------------------------------------------------

subroutine ufo_rttovonedvarcheck_create(self, obspace, f_conf, channels)

  implicit none
  type(ufo_rttovonedvarcheck), intent(inout) :: self
  type(c_ptr), value, intent(in)             :: obspace
  type(fckit_configuration), intent(in)      :: f_conf
  integer(c_int), intent(in)                 :: channels(:)

  self % obsdb = obspace
  self % conf = f_conf

  call ufo_rttovonedvarcheck_setup(self, channels) ! from init

end subroutine ufo_rttovonedvarcheck_create

! ------------------------------------------------------------------------------

subroutine ufo_rttovonedvarcheck_delete(self)

  implicit none
  type(ufo_rttovonedvarcheck), intent(inout) :: self

  if (allocated(self % model_variables)) deallocate(self % model_variables)
  if (allocated(self % channels))        deallocate(self % channels)

end subroutine ufo_rttovonedvarcheck_delete

! ------------------------------------------------------------------------------

subroutine ufo_rttovonedvarcheck_apply(self, vars, geovals, apply)

  ! ------------------------------------------
  ! load modules only used in this subroutine
  ! ------------------------------------------
  use missing_values_mod
  use ufo_rttovonedvarcheck_minimize_newton_mod, only: &
          ufo_rttovonedvarcheck_minimize_newton
  use ufo_rttovonedvarcheck_minimize_ml_mod, only: &
          ufo_rttovonedvarcheck_minimize_ml
  use ufo_rttovonedvarcheck_rmatrix_mod, only: rmatrix_type
  use ufo_rttovonedvarcheck_bmatrix_mod, only: bmatrix_type
  use ufo_rttovonedvarcheck_profindex_mod, only: profindex_type

  implicit none
  type(ufo_rttovonedvarcheck), intent(inout) :: self     ! one d var check setup info
  type(oops_variables), intent(in)           :: vars
  type(ufo_geovals), intent(in)              :: geovals  ! model values at observation space
  logical, intent(in)                        :: apply(:)

  type(bmatrix_type)   :: full_bmatrix
  type(ufo_geovals)    :: local_geovals  !< geoval for one observation
  type(obinfo_type)    :: ob_info
  type(profindex_type) :: prof_index     !< index for mapping geovals to 1d-var state profile
  type(rmatrix_type)   :: r_submatrix    !< r_submatrix object
  character(len=max_string)          :: sensor_id
  character(len=max_string)          :: var
  character(len=max_string)          :: varname
  integer                            :: iloc, jvar, jobs, ivar, band ! counters
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
  logical                            :: file_exists     ! check if a file exists logical
  logical                            :: onedvar_success

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

  ! read in observations and associated errors for full ObsSpace
  do jvar = 1, self%nchans
    var = vars%variable(jvar)
    call obsspace_get_db(self%obsdb, "ObsValue",  trim(var), yobs(jvar,:))
    call obsspace_get_db(self%obsdb, "ObsError",  trim(var), yerr(jvar,:))
    call obsspace_get_db(self%obsdb, "FortranQC", trim(var), QCflags(jvar,:))
  end do

  call obsspace_get_db(self%obsdb,  "MetaData", "latitude", lat(:))
  call obsspace_get_db(self%obsdb,  "MetaData", "longitude", lon(:))
  call obsspace_get_db(self%obsdb,  "MetaData", "elevation", elevation(:))
  call obsspace_get_db(self%obsdb,  "MetaData", "sensor_zenith_angle", sat_zen(:))
  call obsspace_get_db(self%obsdb,  "MetaData", "sensor_azimuth_angle", sat_azi(:))
  call obsspace_get_db(self%obsdb,  "MetaData", "solar_zenith_angle", sol_zen(:))
  call obsspace_get_db(self%obsdb,  "MetaData", "solar_azimuth_angle", sol_azi(:))

  ! Setup full B matrix object
  call full_bmatrix % setup(self % model_variables, self % b_matrix_path, self % qtotal)

  ! Create profile index for mapping 1d-var profile to b-matrix
  call prof_index % setup(full_bmatrix)

  ! Initialize data arrays
  call ufo_rttovonedvarcheck_InitObInfo(ob_info, self%nchans)
  allocate(b_matrix(prof_index % nprofelements,prof_index % nprofelements))
  allocate(b_inverse(prof_index % nprofelements,prof_index % nprofelements))
  allocate(b_sigma(prof_index % nprofelements))

  ! ------------------------------------------
  ! Beginning mains observations loop
  ! ------------------------------------------
  print *,"beginning observations loop: ",self%qcname
  apply_count = 0
  do jobs = 1, iloc
    write(*,*) "Apply for ob number = ",jobs,apply(jobs)
    if (apply(jobs)) then

      apply_count = apply_count + 1
      !---------------------------------------------------
      ! Setup Jb terms
      !---------------------------------------------------
      ! create one ob geovals from full all obs geovals
      call ufo_geovals_copy_one(geovals, local_geovals, jobs)
      call ufo_geovals_print(local_geovals, 1)
  
      ! select appropriate b matrix for latitude of observation
      b_matrix(:,:) = 0.0
      b_inverse(:,:) = 0.0
      b_sigma(:) = 0.0
      do band = 1, full_bmatrix % nbands
        if (lat(jobs) <  full_bmatrix % north(band)) exit
      end do
      b_matrix(:,:) = full_bmatrix % store(:,:,band)
      b_inverse(:,:) = full_bmatrix % inverse(:,:,band)
      b_sigma(:) = full_bmatrix % sigma(:,band)

      !---------------------------------------------------
      ! Setup Jo terms
      !---------------------------------------------------
      ! setup ob data for this observation
      ob_info%forward_mod_name=self%forward_mod_name
      ob_info%latitude=lat(jobs)
      ob_info%longitude=lon(jobs)
      ob_info%elevation=elevation(jobs)
      ob_info%sensor_zenith_angle=sat_zen(jobs)
      ob_info%sensor_azimuth_angle=sat_azi(jobs)
      ob_info%solar_zenith_angle=sol_zen(jobs)
      ob_info%solar_azimuth_angle=sol_azi(jobs)

      ! Channel selection based on previous filter flags
      chans_used = 0
      do jvar = 1, self%nchans
        if( QCflags(jvar,jobs) == 0 ) then
          chans_used = chans_used + 1
        end if
      end do

      ! allocate arrays
      allocate(ob_info%yobs(chans_used))
      allocate(channels_used(chans_used))
      allocate(obs_error(chans_used))

      ! create obs vector and r matrix
      jchans_used = 0
      do jvar = 1, self%nchans
        if( QCflags(jvar,jobs) == 0 ) then
          jchans_used = jchans_used + 1
          obs_error(jchans_used) = yerr(jvar, jobs)
          write(*,*) "jchans_used err = ",jchans_used,obs_error(jchans_used)
          ob_info % yobs(jchans_used) = yobs(jvar, jobs)
          channels_used(jchans_used) = self%channels(jvar)
        end if
      end do
      call r_submatrix % setup(self % rtype, chans_used, obs_error)
      call r_submatrix % info()

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
        QCflags(:,jobs) = 1
      end if

      ! Tidy up memory specific to a single observation
      if (allocated(ob_info%yobs))  deallocate(ob_info%yobs)
      if (allocated(channels_used)) deallocate(channels_used)
      if (allocated(obs_error))     deallocate(obs_error)
      call r_submatrix % delete()
      call ufo_geovals_delete(local_geovals)

    endif
  end do

  write(*,*) "Number being tested by 1dvar = ",apply_count

  ! Put QC flags back in database
  do jvar = 1, self%nchans
    var = vars%variable(jvar)
    call obsspace_put_db(self%obsdb, "FortranQC", trim(var), QCflags(jvar,:))
  end do

  ! tidy up
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

end subroutine ufo_rttovonedvarcheck_apply

end module ufo_rttovonedvarcheck_mod
