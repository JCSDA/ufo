!
! (C) Copyright 2020 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
!
module test_sampled_locations

use iso_c_binding

implicit none
private

contains

! ------------------------------------------------------------------------------
!> Tests whether latitudes and longitudes from ufo::SampledLocations object
!! are matching references in yaml file
integer function test_sampled_locations_lonslats_c(c_locs, c_conf) &
  bind(c,name='test_sampled_locations_lonslats_f90')
use fckit_configuration_module, only: fckit_configuration
use fckit_log_module, only: fckit_log
use ufo_sampled_locations_mod
implicit none
type(c_ptr), value, intent(in) :: c_locs  !< ufo::SampledLocations object
type(c_ptr), value, intent(in) :: c_conf  !< test configuration (described in SampledLocationsTestParameters)

!> local variables
type(ufo_sampled_locations) :: locs
type(fckit_configuration) :: f_conf
integer :: npaths, ipath
real    :: tolerance = 1.e-8
real(c_double), dimension(:), allocatable :: lons_ref, lats_ref, lons, lats
character(len=200) :: logmessage

!> default value: test passed
test_sampled_locations_lonslats_c = 1

!> read reference longitudes and latitudes
f_conf = fckit_configuration(c_conf)
npaths = f_conf%get_size("reference lons")
call f_conf%get_or_die("reference lons", lons_ref)
call f_conf%get_or_die("reference lats", lats_ref)

!> get longitudes and latitudes from sampled_locations object
locs = ufo_sampled_locations(c_locs)
allocate(lons(npaths), lats(npaths))
call locs%get_lons(lons)
call locs%get_lats(lats)

do ipath = 1, npaths
  write(logmessage, *) "lons: path and ref: ", lons(ipath), ", ", lons_ref(ipath)
  call fckit_log%debug(logmessage)
  write(logmessage, *) "lats: path and ref: ", lats(ipath), ", ", lats_ref(ipath)
  call fckit_log%debug(logmessage)
enddo

!> compare to reference
if(any(abs(lons-lons_ref) > tolerance)) test_sampled_locations_lonslats_c = 0
if(any(abs(lats-lats_ref) > tolerance)) test_sampled_locations_lonslats_c = 0

deallocate(lons, lats)

end function test_sampled_locations_lonslats_c

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------
!> Tests whether ufo::SampledLocations time mask method results
!! are matching references in yaml file
integer function test_sampled_locations_timemask_c(c_locs, c_conf) &
  bind(c,name='test_sampled_locations_timemask_f90')
use fckit_configuration_module, only: fckit_configuration
use fckit_log_module, only: fckit_log
use datetime_mod
use ufo_sampled_locations_mod
implicit none
type(c_ptr), value, intent(in) :: c_locs  !< ufo::SampledLocations object
type(c_ptr), value, intent(in) :: c_conf  !< test configuration (described in SampledLocationsTestParameters

!> local variables
type(ufo_sampled_locations) :: locs
type(fckit_configuration) :: f_conf
type(fckit_configuration), allocatable :: testconfigs(:)
character(kind=c_char,len=:), allocatable :: t1str, t2str
logical(kind=c_bool), allocatable :: mask(:)
logical, allocatable :: mask_ref(:)
type(datetime) :: t1, t2
integer :: npaths, ipath, itest
character(len=200) :: logmessage

!> default value: test passed
test_sampled_locations_timemask_c = 1

locs = ufo_sampled_locations(c_locs)
f_conf = fckit_configuration(c_conf)
call f_conf%get_or_die("time mask tests", testconfigs)
!> loop over all the tests for time masks
do itest = 1, size(testconfigs)
  call testconfigs(itest)%get_or_die("t1", t1str)
  call testconfigs(itest)%get_or_die("t2", t2str)
  call datetime_create(t1str, t1)
  call datetime_create(t2str, t2)
  call testconfigs(itest)%get_or_die("reference mask", mask_ref)

  !> call ufo::SampledLocations method to compute time mask
  npaths = locs%npaths()
  allocate(mask(npaths))
  call locs%get_timemask(t1, t2, mask)

  write(logmessage, *) "Test ", itest
  call fckit_log%debug(logmessage)
  do ipath = 1, npaths
    write(logmessage, *) "mask: path and ref: ", mask(ipath), ", ", mask_ref(ipath)
    call fckit_log%debug(logmessage)
  enddo

  !> compare to reference
  if(any(mask .neqv. mask_ref)) test_sampled_locations_timemask_c = 0

  deallocate(mask)

enddo

end function test_sampled_locations_timemask_c

end module test_sampled_locations
