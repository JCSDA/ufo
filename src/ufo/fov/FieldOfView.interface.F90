! (C) Copyright 2024 UCAR.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

!> Fortran module providing a C interface to the Fortran implementation of satellite FOV

module ufo_fov_mod_c

  use iso_c_binding
  use ufo_fov_mod, only: ufo_fov

  implicit none
  private

  ! ------------------------------------------------------------------------------
#define LISTED_TYPE ufo_fov

  !> Linked list interface - defines registry_t type
#include "oops/util/linkedList_i.f"

  !> Global registry
  type(registry_t) :: ufo_fov_registry

  ! ------------------------------------------------------------------------------

contains

  ! ------------------------------------------------------------------------------
  !> Linked list implementation
#include "oops/util/linkedList_c.f"

! ------------------------------------------------------------------------------

subroutine ufo_fov_setup_c(c_key_self, sensor_len, sensor_cstr, platform_len, platform_cstr, &
                           valid, npoly) bind(c, name='ufo_fov_setup_f90')
  use string_f_c_mod, only: c_f_string

  integer(c_int), intent(inout) :: c_key_self
  integer(c_int), intent(in) :: sensor_len
  character(kind=c_char, len=1), intent(in) :: sensor_cstr(sensor_len + 1)
  integer(c_int), intent(in) :: platform_len
  character(kind=c_char, len=1), intent(in) :: platform_cstr(platform_len + 1)
  logical(c_bool), intent(out) :: valid
  integer(c_int), intent(out) :: npoly

  type(ufo_fov), pointer :: self
  character(len=sensor_len) :: sensor
  character(len=platform_len) :: platform
  logical :: valid_local

  ! Copy C char* into Fortran char array
  call c_f_string(sensor_cstr, sensor)
  call c_f_string(platform_cstr, platform)

  call ufo_fov_registry%setup(c_key_self, self)
  call self%setup(sensor, platform, valid_local, npoly)

  valid = valid_local

end subroutine ufo_fov_setup_c

! ------------------------------------------------------------------------------

subroutine ufo_fov_delete_c(c_key_self) bind(c, name='ufo_fov_delete_f90')
  integer(c_int), intent(inout) :: c_key_self

  type(ufo_fov), pointer :: self

  ! even though ufo_fov type itself has no allocatable data, it wraps GSI code
  ! that does have allocated arrays, so we must delete it
  call ufo_fov_registry%get(c_key_self, self)
  call self%delete
  call ufo_fov_registry%remove(c_key_self)

end subroutine ufo_fov_delete_c

! ------------------------------------------------------------------------------

subroutine ufo_fov_ellipse_c(c_key_self, sensor_len, sensor_cstr, scan_position, &
                             sat_azimuth_angle, fov_center_lon, fov_center_lat, npoly, &
                             fov_ellipse_lons, fov_ellipse_lats) bind(c, name='ufo_fov_ellipse_f90')
  use string_f_c_mod, only: c_f_string

  integer(c_int), intent(inout) :: c_key_self
  integer(c_int), intent(in) :: sensor_len
  character(kind=c_char, len=1), intent(in) :: sensor_cstr(sensor_len + 1)
  integer(c_int), intent(in) :: scan_position
  real(c_double), intent(in) :: sat_azimuth_angle
  real(c_double), intent(in) :: fov_center_lon
  real(c_double), intent(in) :: fov_center_lat
  integer(c_int), intent(in) :: npoly
  real(c_double), intent(out) :: fov_ellipse_lons(npoly)
  real(c_double), intent(out) :: fov_ellipse_lats(npoly)

  type(ufo_fov), pointer :: self

  character(len=sensor_len) :: sensor

  ! Copy C char* into Fortran char array
  call c_f_string(sensor_cstr, sensor)

  call ufo_fov_registry%get(c_key_self, self)
  call self%fov_ellipse(sensor, scan_position, sat_azimuth_angle, fov_center_lon, fov_center_lat, &
                        fov_ellipse_lons, fov_ellipse_lats)

end subroutine ufo_fov_ellipse_c

! ------------------------------------------------------------------------------

subroutine ufo_antenna_power_within_fov_c(c_key_self, sensor_len, sensor_cstr, scan_position, &
                                          sat_azimuth_angle, fov_center_lon, fov_center_lat, &
                                          test_lon, test_lat, antenna_power) &
                                          bind(c, name='ufo_antenna_power_within_fov_f90')
  use string_f_c_mod, only: c_f_string

  integer(c_int), intent(inout) :: c_key_self
  integer(c_int), intent(in) :: sensor_len
  character(kind=c_char, len=1), intent(in) :: sensor_cstr(sensor_len + 1)
  integer(c_int), intent(in) :: scan_position
  real(c_double), intent(in) :: sat_azimuth_angle
  real(c_double), intent(in) :: fov_center_lon
  real(c_double), intent(in) :: fov_center_lat
  real(c_double), intent(in) :: test_lon
  real(c_double), intent(in) :: test_lat
  real(c_double), intent(out) :: antenna_power

  type(ufo_fov), pointer :: self
  character(len=sensor_len) :: sensor

  ! Copy C char* into Fortran char array
  call c_f_string(sensor_cstr, sensor)

  call ufo_fov_registry%get(c_key_self, self)
  call self%antenna_power_within_fov(sensor, scan_position, sat_azimuth_angle, fov_center_lon, &
                                     fov_center_lat, test_lon, test_lat, antenna_power)

end subroutine ufo_antenna_power_within_fov_c

! ------------------------------------------------------------------------------

end module ufo_fov_mod_c
