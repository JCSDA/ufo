! (C) Copyright 2017-2018 ucar
! 
! this software is licensed under the terms of the apache licence version 2.0
! which can be obtained at http://www.apache.org/licenses/license-2.0.

!> fortran module to implement onedvar fortran check

module ufo_rttovonedvarcheck_mod

use iso_c_binding
use kinds
use ufo_geovals_mod
use ufo_vars_mod
use oops_variables_mod
use obsspace_mod
use config_mod
use ufo_rttovonedvarcheck_utils_mod
use ufo_rttovonedvarcheck_init_mod
use ufo_rttovonedvarcheck_process_mod

implicit none
public :: ufo_rttovonedvarcheck
public :: ufo_rttovonedvarcheck_create, ufo_rttovonedvarcheck_delete
public :: ufo_rttovonedvarcheck_prior, ufo_rttovonedvarcheck_post

! ------------------------------------------------------------------------------
type :: ufo_rttovonedvarcheck
  character(len=max_string_length) :: qcname
  character(len=max_string_length) :: b_matrix_path
  character(len=max_string_length) :: forward_mod_name
  integer     :: qtotal
  type(c_ptr) :: obsdb
  type(c_ptr) :: conf
  integer :: nmvars
  character(len=max_string_length), allocatable :: model_variables(:)
  integer :: nchans
  integer(c_int), allocatable        :: channels(:)
end type ufo_rttovonedvarcheck

! ------------------------------------------------------------------------------
contains
! ------------------------------------------------------------------------------

subroutine ufo_rttovonedvarcheck_create(self, obspace, conf, channels)

  implicit none
  type(ufo_rttovonedvarcheck), intent(inout) :: self
  type(c_ptr), value, intent(in)          :: obspace
  type(c_ptr), value, intent(in)          :: conf
  integer(c_int), intent(in)              :: channels(:)

  self%qcname = "rttovonedvarcheck"
  self%b_matrix_path = config_get_string(conf, max_string_length, "BMatrix")
  self%forward_mod_name = config_get_string(conf, max_string_length, "ModName")
  self%qtotal = config_get_int(conf, "qtotal")
  self%obsdb = obspace
  self%conf = conf

  self%nmvars = size(config_get_string_vector(conf, max_string_length, "model_variables"))
  allocate(self%model_variables(self%nmvars))
  self%model_variables = config_get_string_vector(conf, max_string_length, "model_variables")

  self%nchans = size(channels)
  allocate(self%channels(self%nchans))
  self%channels(:) = channels(:)

  write(*,*) "nchans setup = ",self%nchans
  write(*,*) "channels setup = ",self%channels

end subroutine ufo_rttovonedvarcheck_create

! ------------------------------------------------------------------------------

subroutine ufo_rttovonedvarcheck_delete(self)

  implicit none
  type(ufo_rttovonedvarcheck), intent(inout) :: self

  if (allocated(self % model_variables)) deallocate(self % model_variables)
  if (allocated(self % channels))        deallocate(self % channels)

end subroutine ufo_rttovonedvarcheck_delete

! ------------------------------------------------------------------------------

subroutine ufo_rttovonedvarcheck_prior(self, geovals)

  implicit none
  type(ufo_rttovonedvarcheck), intent(in) :: self
  type(ufo_geovals), intent(in) :: geovals

end subroutine ufo_rttovonedvarcheck_prior

! ------------------------------------------------------------------------------

subroutine ufo_rttovonedvarcheck_post(self, vars, geovals, apply)

  ! ------------------------------------------
  ! load modules only used in this subroutine
  ! ------------------------------------------
  use missing_values_mod
  use ufo_rttovonedvarcheck_minimize_newton_mod, only: &
                    Ops_SatRad_MinimizeNewton_RTTOV12
  use ufo_rttovonedvarcheck_minimize_ml_mod, only: &
                    ufo_rttovonedvarcheck_MinimizeML

  implicit none
  type(ufo_rttovonedvarcheck), intent(inout) :: self        ! one d var check setup info
  type(oops_variables), intent(in)        :: vars
  type(ufo_geovals), intent(in)           :: geovals     ! model values at observation space
  logical, intent(in)                     :: apply(:)

  integer, allocatable               :: fields_in(:)
  integer                            :: iloc, jvar, jobs, ivar, band ! counters
  integer                            :: chans_used      ! counter for number of channels used for an ob
  integer                            :: jchans_used
  integer                            :: fileunit        ! unit number for reading in files
  integer(c_int32_t), allocatable    :: flags(:,:)      ! qc flag for return to var obs file
  real(kind_real)                    :: missing         ! missing value
  character(len=max_string_length)   :: var
  character(len=100)                 :: varname
  real(kind_real), allocatable       :: r_matrix(:,:)   ! r-matrix = obs + forward model error
  real(kind_real), allocatable       :: r_inverse(:,:)  ! inverse of the r-matrix 
  type(bmatrix_type)                 :: full_b_matrix   ! full b matrix
  real(kind_real), allocatable       :: b_matrix(:,:)   ! 1d-var profile b matrix
  real(kind_real), allocatable       :: b_inverse(:,:)  ! inverse for each 1d-var profile b matrix
  type(profileinfo_type)             :: profile_index   ! index for mapping geovals to 1d-var state profile
  integer                            :: nprofelements   ! number of elements in 1d-var state profile
  type(ufo_geovals)                  :: local_geovals   ! geoval for one observation
  real(c_double), allocatable        :: iter_hofx(:)    ! model equivalent of observations during 1d-var
  real(kind_real)                    :: obs_error
  integer, allocatable               :: channels_used(:)

  ! variables from ioda
  real(kind_real), allocatable       :: yobs(:,:)       ! observation value from obs files
  real(kind_real), allocatable       :: yerr(:,:)       ! observation error from obs files
  integer, allocatable               :: QCflags(:,:)    ! current qc flags needed for channel selection
  real(kind_real), allocatable       :: lat(:)          ! observation latitude
  real(kind_real), allocatable       :: lon(:)          ! observation longitude
  real(kind_real), allocatable       :: elevation(:)    ! observation elevation
  real(kind_real), allocatable       :: sat_zen(:)      ! observation satellite zenith angle
  real(kind_real), allocatable       :: sat_azi(:)      ! observation satellite azimuth angle
  real(kind_real), allocatable       :: sol_zen(:)      ! observation solar zenith angle
  real(kind_real), allocatable       :: sol_azi(:)      ! observation solar azimuth angle

  real                               :: t2,t1           ! timing

  logical                            :: file_exists     ! check if a file exists logical
  logical                            :: onedvar_success
  character(len=255)                 :: sensor_id
  integer                            :: apply_count

  type(obinfo_type)                  :: ob_info

  ! ------------------------------------------
  ! 1. setup
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
  allocate(iter_hofx(self%nchans))

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
  iter_hofx(:) = 0.0

  ! read in observations and associated errors
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

  !QCflags(1,1) = 1

  call cpu_time(t1)

  ! Select fields for B matrix
  allocate(fields_in(7))
  fields_in(:) = 0
  do jvar = 1, self % nmvars
    varname = self % model_variables(jvar)
    select case (trim(varname))

      case ("air_temperature")
        fields_in(1) = 1 ! air_temperature

      case ("specific_humidity")
        if (self%qtotal == 1) then
          fields_in(2) = 10 ! water profile in bmatrix - specific humidity in geovals!?!
        else
          fields_in(2) = 2 ! water profile in bmatrix - specific humidity in geovals!?!
        end if

      case("air_temperature_at_two_meters_above_surface")
        fields_in(3) = 3 ! 2m air_temperature

      case("specific_humidity_at_two_meters_above_surface")
        fields_in(4) = 4 ! 2m specific_humidity

      case("skin_temperature")
        fields_in(5) = 5 ! surface skin temperature

      case("surface_air_pressure")
        fields_in(6) = 6 ! surface air pressure

      case default
        write(*,*) "Variable not implemented yet in OneDVarCheck Covariance"

    end select
  end do
  write(*,*) "fields in = ",fields_in(:)

  ! setup and read in background covariance
  inquire(file=trim(self%b_matrix_path), exist=file_exists)
  if (file_exists) then
    call ufo_rttovonedvarcheck_IOGetFreeUnit(fileunit)
    open(unit = fileunit, file = trim(self%b_matrix_path))
    call ufo_rttovonedvarcheck_InitBmatrix(full_b_matrix)
    !call ufo_rttovonedvarcheck_GetBmatrix(fileunit, full_b_matrix)
    call ufo_rttovonedvarcheck_GetBmatrix(fileunit, full_b_matrix, fieldlist=fields_in)
    close(unit = fileunit)
  else
    write(*,*) "rttovonedvarcheck bmatrix file not found"
  end if

  ! get mapping for 1d-var profile from b-matrix
  call ufo_rttovonedvarcheck_InitProfInfo(profile_index)
  call ufo_rttovonedvarcheck_MapProfileToB(full_b_matrix, profile_index, nprofelements)
  write(*,*) "1DVar number of profile elements = ",nprofelements
  allocate(b_matrix(nprofelements,nprofelements))
  allocate(b_inverse(nprofelements,nprofelements))

  call cpu_time(t2)

  write(*,*) "time to read in b matrix and map profile = ",(t2-t1)

  ! initialize ob info type
  call ufo_rttovonedvarcheck_InitObInfo(ob_info,self%nchans)

  ! print geovals infor
  !call ufo_geovals_print(geovals, 1)

  ! ------------------------------------------
  ! 2. beginning mains observations loop
  ! ------------------------------------------
  print *,"beginning observations loop: ",self%qcname
  apply_count = 0
  do jobs = 1, iloc
  !do jobs = 1, 1
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
      do band = 1, full_b_matrix % nbands
        if (lat(jobs) <  full_b_matrix % north(band)) exit
      end do
      write(*,*) "band = ",band
      b_matrix(:,:) = full_b_matrix % store(:,:,band)
      b_inverse(:,:) = full_b_matrix % inverse(:,:,band)
      write(*,'(144E13.5)') b_matrix(:,:)

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
      allocate(r_matrix(chans_used,chans_used))
      allocate(r_inverse(chans_used,chans_used))
      allocate(channels_used(chans_used))

      ! create obs vector and r matrix
      r_matrix(:,:) = 0.0_kind_real
      r_inverse(:,:) = 0.0_kind_real
      jchans_used = 0
      do jvar = 1, self%nchans
        if( QCflags(jvar,jobs) == 0 ) then
          jchans_used = jchans_used + 1
          obs_error = yerr(jvar, jobs)
          write(*,*) "jchans_used err = ",jchans_used,obs_error
          r_matrix(jchans_used,jchans_used) = obs_error * obs_error
          r_inverse(jchans_used,jchans_used) = 1.0_kind_real / (obs_error * obs_error)
          ob_info%yobs(jchans_used) = yobs(jvar, jobs)
          channels_used(jchans_used) = self%channels(jvar)
        end if
      end do

      write(*,*) "Ob number = ",jobs
      write(*,*) "channels used = ",channels_used(:)
      write(*,*) "channels used number = ",chans_used
      write(*,*) "r_inverse = ",r_inverse
      write(*,*) "r_matrix = ",r_matrix

      !---------------------------------------------------
      ! Call minimization
      !---------------------------------------------------
      ! do 1d-var using marquardt-levenberg
      !call ufo_rttovonedvarcheck_MinimizeML(ob_info, r_matrix, r_inverse, b_matrix, &
      !                                   b_inverse, local_geovals, & 
      !                                   profile_index, nprofelements, self%conf, &
      !                                   self%obsdb, channels_used, onedvar_success)
      call Ops_SatRad_MinimizeNewton_RTTOV12(ob_info, r_matrix, r_inverse, b_matrix, &
                                         b_inverse, local_geovals, & 
                                         profile_index, nprofelements, self%conf, &
                                         self%obsdb, channels_used, onedvar_success)
    
      ! Set QCflags based on output from minimization
      if (.NOT. onedvar_success) then
        do jvar = 1, self%nchans
          QCflags(jvar,jobs) = 1
        end do
      end if

      ! Deallocate arrays allocated in loop
      if (allocated(ob_info%yobs))  deallocate(ob_info%yobs)
      if (allocated(r_matrix))      deallocate(r_matrix)
      if (allocated(r_inverse))     deallocate(r_inverse)
      if (allocated(channels_used)) deallocate(channels_used)

      ! tidy up
      call ufo_geovals_delete(local_geovals)

      write(*,*) "iteration = ",jobs
      write(*,*) "QCflags = ",QCflags(:,jobs)

    endif
  end do

  write(*,*) "Number being tested by 1dvar = ",apply_count

  ! Put QC flags back in database
  do jvar = 1, self%nchans
    var = vars%variable(jvar)
    call obsspace_put_db(self%obsdb, "FortranQC", trim(var), QCflags(jvar,:))
  end do

  ! Put QC flags back in database
  !do jvar = 1, self%nchans
  !  var = vars%variable(jvar)
  !  write(*,*) "QCflags(jvar,:) = ",QCflags(jvar,:)
  !  call obsspace_get_db(self%obsdb, "FortranQC", trim(var), QCflags(jvar,:))
  !end do

  ! tidy up
  if (allocated(yobs))            deallocate(yobs)
  if (allocated(yerr))            deallocate(yerr)
  if (allocated(lat))             deallocate(lat)
  if (allocated(lon))             deallocate(lon)
  if (allocated(elevation))       deallocate(elevation)
  if (allocated(sat_zen))         deallocate(sat_zen)
  if (allocated(sat_azi))         deallocate(sat_azi)
  if (allocated(sol_zen))         deallocate(sol_zen)
  if (allocated(sol_azi))         deallocate(sol_azi)
  if (allocated(flags))           deallocate(flags)
  if (allocated(r_matrix))        deallocate(r_matrix)
  if (allocated(r_inverse))       deallocate(r_inverse)
  if (allocated(b_matrix))        deallocate(b_matrix)
  if (allocated(b_inverse))       deallocate(b_inverse)
  if (allocated(iter_hofx))       deallocate(iter_hofx)

end subroutine ufo_rttovonedvarcheck_post

end module ufo_rttovonedvarcheck_mod
