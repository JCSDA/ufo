! (C) Copyright 2017-2018 ucar
! 
! this software is licensed under the terms of the apache licence version 2.0
! which can be obtained at http://www.apache.org/licenses/license-2.0.

!> fortran module to implement onedvar fortran check

module ufo_onedvarfortran_mod

use iso_c_binding
use kinds
use ufo_geovals_mod
use ufo_vars_mod
use obsspace_mod
use config_mod
use ufo_onedvarfortran_utils_mod
use ufo_onedvarfortran_init_mod
use ufo_onedvarfortran_process_mod

implicit none
public  :: ufo_onedvarfortran, ufo_onedvarfortran_create, ufo_onedvarfortran_delete, ufo_onedvarfortran_prior, ufo_onedvarfortran_post

! ------------------------------------------------------------------------------
type :: ufo_onedvarfortran
  character(len=max_string_length) :: qcname
  character(len=max_string_length) :: b_matrix_path
  character(len=max_string_length) :: forward_mod_name
  type(c_ptr) :: obsdb
  type(c_ptr) :: conf
  integer :: nvars
  character(len=max_string_length), allocatable :: variables(:)
  integer :: nmvars
  character(len=max_string_length), allocatable :: model_variables(:)
  integer :: nchans
  integer(c_int), allocatable        :: channels(:)
end type ufo_onedvarfortran

! ------------------------------------------------------------------------------
contains
! ------------------------------------------------------------------------------

subroutine ufo_onedvarfortran_create(self, obspace, conf, channels)

  implicit none
  type(ufo_onedvarfortran), intent(inout) :: self
  type(c_ptr), value, intent(in)        :: obspace
  type(c_ptr), value, intent(in)        :: conf
  integer(c_int), intent(in)            :: channels(:)

  self%qcname = "onedvarfortran"
  self%b_matrix_path = config_get_string(conf, max_string_length, "BMatrix")
  self%forward_mod_name = config_get_string(conf, max_string_length, "ModName")
  self%obsdb = obspace
  self%conf = conf
  
  self%nvars = size(config_get_string_vector(conf, max_string_length, "variables"))
  allocate(self%variables(self%nvars))
  self%variables = config_get_string_vector(conf, max_string_length, "variables")

  self%nmvars = size(config_get_string_vector(conf, max_string_length, "model_variables"))
  allocate(self%model_variables(self%nmvars))
  self%model_variables = config_get_string_vector(conf, max_string_length, "model_variables")

  self%nchans = size(channels)
  allocate(self%channels(self%nchans))
  self%channels(:) = channels(:)

  write(*,*) "nchans setup = ",self%nchans
  write(*,*) "channels setup = ",self%channels

end subroutine ufo_onedvarfortran_create

! ------------------------------------------------------------------------------

subroutine ufo_onedvarfortran_delete(self)

  implicit none
  type(ufo_onedvarfortran), intent(inout) :: self

  if (allocated(self % variables))       deallocate(self % variables)
  if (allocated(self % model_variables)) deallocate(self % model_variables)
  if (allocated(self % channels))        deallocate(self % channels)

end subroutine ufo_onedvarfortran_delete

! ------------------------------------------------------------------------------

subroutine ufo_onedvarfortran_prior(self, geovals)

  implicit none
  type(ufo_onedvarfortran), intent(in) :: self
  type(ufo_geovals), intent(in) :: geovals

end subroutine ufo_onedvarfortran_prior

! ------------------------------------------------------------------------------

subroutine ufo_onedvarfortran_post(self, hofx, hofxvars, geovals, conf)

  ! ------------------------------------------
  ! load modules only used in this subroutine
  ! ------------------------------------------
  use missing_values_mod

  implicit none
  type(ufo_onedvarfortran), intent(inout) :: self        ! one d var check setup info
  real(c_double),  intent(in)             :: hofx(:,:)   ! model equivalent of observations
  character(len=max_string_length)        :: hofxvars(:) ! channels names
  type(ufo_geovals), intent(in)           :: geovals     ! model values at observation space
  type(c_ptr), intent(in)                 :: conf        ! pointer to conf

  integer, allocatable               :: fields_in(:)
  integer                            :: iloc, jvar, jobs, ivar, band ! counters
  integer                            :: fileunit        ! unit number for reading in files
  integer(c_int32_t), allocatable    :: flags(:,:)      ! qc flag for return to var obs file
  real(kind_real)                    :: missing         ! missing value
  character(len=max_string_length)   :: var
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

  ! variables from ioda
  real(kind_real), allocatable       :: yobs(:,:)       ! observation value from obs files
  real(kind_real), allocatable       :: yerr(:,:)       ! observation error from obs files
  real(kind_real), allocatable       :: lat(:)          ! observation latitude
  real(kind_real), allocatable       :: lon(:)          ! observation longitude
  real(kind_real), allocatable       :: elevation(:)    ! observation elevation
  real(kind_real), allocatable       :: sat_zen(:)      ! observation satellite zenith angle
  real(kind_real), allocatable       :: sat_azi(:)      ! observation satellite azimuth angle
  real(kind_real), allocatable       :: sol_zen(:)      ! observation solar zenith angle
  real(kind_real), allocatable       :: sol_azi(:)      ! observation solar azimuth angle

  real                               :: t2,t1           ! timing

  logical                            :: file_exists     ! check if a file exists logical
  character(len=255)                 :: sensor_id

  type(obinfo_type)                  :: ob_info

  ! ------------------------------------------
  ! 1. setup
  ! ------------------------------------------
  missing = missing_value(missing)
  iloc = obsspace_get_nlocs(self%obsdb)

  ! allocate arrays
  allocate(yobs(self%nvars,iloc))
  allocate(yerr(self%nvars,iloc))
  allocate(lat(iloc))
  allocate(lon(iloc))
  allocate(elevation(iloc))
  allocate(sat_zen(iloc))
  allocate(sat_azi(iloc))
  allocate(sol_zen(iloc))
  allocate(sol_azi(iloc))
  allocate(flags(self%nvars,iloc))
  allocate(r_matrix(self%nvars,self%nvars))
  allocate(r_inverse(self%nvars,self%nvars))
  allocate(iter_hofx(self%nvars))

  ! initialize arrays
  yobs(:,:) = 0.0
  yerr(:,:) = 0.0
  lat(:) = 0.0
  lon(:) = 0.0
  elevation(iloc) = 0.0
  sat_zen(:) = 0.0
  sat_azi(:) = 0.0
  sol_zen(:) = 0.0
  sol_azi(:) = 0.0
  flags(:,:) = 0
  iter_hofx(:) = 0.0

  ! read in observations and associated errors
  do jvar = 1, self%nvars
    var = trim(self%variables(jvar))
    ivar = find_index(hofxvars, var)

    ! read in observations and associated errors
    call obsspace_get_db(self%obsdb,  "ObsValue", var       , yobs(ivar,:))
    call obsspace_get_db(self%obsdb,  "ObsError", var       , yerr(ivar,:))
  end do

  call obsspace_get_db(self%obsdb,  "MetaData", "latitude", lat(:))
  call obsspace_get_db(self%obsdb,  "MetaData", "longitude", lon(:))
  call obsspace_get_db(self%obsdb,  "MetaData", "elevation", elevation(:))
  call obsspace_get_db(self%obsdb,  "MetaData", "sensor_zenith_angle", sat_zen(:))
  call obsspace_get_db(self%obsdb,  "MetaData", "sensor_azimuth_angle", sat_azi(:))
  call obsspace_get_db(self%obsdb,  "MetaData", "solar_zenith_angle", sol_zen(:))
  call obsspace_get_db(self%obsdb,  "MetaData", "solar_azimuth_angle", sol_azi(:))

  call cpu_time(t1)

  ! Just the air temp for now
  allocate(fields_in(7))
  fields_in(:) = 0
  fields_in(1) = 1 ! just air_temperature for now

  ! setup and read in background covariance
  inquire(file=trim(self%b_matrix_path), exist=file_exists)
  if (file_exists) then
    call ufo_onedvarfortran_IOGetFreeUnit(fileunit)
    open(unit = fileunit, file = trim(self%b_matrix_path))
    call ufo_onedvarfortran_InitBmatrix(full_b_matrix)
    !call ufo_onedvarfortran_GetBmatrix(fileunit, full_b_matrix)
    call ufo_onedvarfortran_GetBmatrix(fileunit, full_b_matrix, fieldlist=fields_in)
    close(unit = fileunit)
  else
    write(*,*) "onedvarfortran bmatrix file not found"
  end if

  ! get mapping for 1d-var profile from b-matrix
  call ufo_onedvarfortran_InitProfInfo(profile_index)
  call ufo_onedvarfortran_MapProfileToB(full_b_matrix, profile_index, nprofelements)
  allocate(b_matrix(nprofelements,nprofelements))
  allocate(b_inverse(nprofelements,nprofelements))

  call cpu_time(t2)

  write(*,*) "time to read in b matrix and map profile = ",(t2-t1)

  ! initialize ob info type
  call ufo_onedvarfortran_InitObInfo(ob_info,self%nvars)

  ! print geovals infor
  !call ufo_geovals_print(geovals, 1)

  ! ------------------------------------------
  ! 2. beginning mains observations loop
  ! ------------------------------------------
  print *,"beginning observations loop: ",self%qcname
  do jobs = 1, iloc

    ! reset loop variables
    b_matrix(:,:) = 0.0
    b_inverse(:,:) = 0.0
    r_matrix(:,:) = 0.0
    r_inverse(:,:) = 0.0

    ! create one ob geovals from full all obs geovals
    call ufo_onedvarfortran_LocalGeovals(geovals,jobs,local_geovals)
    call ufo_geovals_print(geovals, 1)

    ! create r matrix
    do jvar = 1, self%nvars
      obs_error = yerr(jvar, jobs)
      r_matrix(jvar,jvar) = obs_error * obs_error
      r_inverse(jvar,jvar) = 1.0 / (obs_error * obs_error)
    end do

    ! select appropriate b matrix for latitude of observation
    do band = 1, full_b_matrix % nbands
      if (lat(jobs) <  full_b_matrix % north(band)) exit
    end do
    write(*,*) "band = ",band
    b_matrix(:,:) = full_b_matrix % store(:,:,band)
    b_inverse(:,:) = full_b_matrix % inverse(:,:,band)
    
    ! setup ob data for this observation
    ob_info%forward_mod_name=self%forward_mod_name
    ob_info%latitude=lat(jobs)
    ob_info%longitude=lon(jobs)
    ob_info%elevation=elevation(jobs)
    ob_info%sensor_zenith_angle=sat_zen(jobs)
    ob_info%sensor_azimuth_angle=sat_azi(jobs)
    ob_info%solar_zenith_angle=sol_zen(jobs)
    ob_info%solar_azimuth_angle=sol_azi(jobs)
    ob_info%hofx=hofx(:, jobs)
    ob_info%yobs=yobs(:, jobs)

    ! do 1d-var using marquardt-levenberg
    call ufo_onedvarfortran_Minimize(ob_info, r_matrix, r_inverse, b_matrix, &
                                   b_inverse, local_geovals, & 
                                   profile_index, nprofelements, conf, &
                                   self%obsdb, self%channels)

    ! tidy up
    call ufo_geovals_delete(local_geovals)
    
  end do

  ! tidy up
  if (allocated(yobs))      deallocate(yobs)
  if (allocated(yerr))      deallocate(yerr)
  if (allocated(lat))       deallocate(lat)
  if (allocated(lon))       deallocate(lon)
  if (allocated(elevation)) deallocate(elevation)
  if (allocated(sat_zen))   deallocate(sat_zen)
  if (allocated(sat_azi))   deallocate(sat_azi)
  if (allocated(sol_zen))   deallocate(sol_zen)
  if (allocated(sol_azi))   deallocate(sol_azi)
  if (allocated(flags))     deallocate(flags)
  if (allocated(r_matrix))  deallocate(r_matrix)
  if (allocated(r_inverse)) deallocate(r_inverse)
  if (allocated(b_matrix))  deallocate(b_matrix)
  if (allocated(b_inverse)) deallocate(b_inverse)
  if (allocated(iter_hofx)) deallocate(iter_hofx)

end subroutine ufo_onedvarfortran_post

end module ufo_onedvarfortran_mod
