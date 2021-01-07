!
! (C) Copyright 2020 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
!
module test_locations

use iso_c_binding

implicit none
private

contains

! ------------------------------------------------------------------------------
!> Tests whether latitudes and longitudes from ufo::Locations object
!! are matching references in yaml file
integer function test_locations_lonslats_c(c_locs, c_conf) bind(c,name='test_locations_lonslats_f90')
use fckit_configuration_module, only: fckit_configuration
use fckit_log_module, only: fckit_log
use ufo_locations_mod
implicit none
type(c_ptr), value, intent(in) :: c_locs  !< ufo::Locations object
type(c_ptr), value, intent(in) :: c_conf  !< test configuration (described in LocationsTestParameters

!> local variables
type(ufo_locations) :: locs
type(fckit_configuration) :: f_conf
integer :: nlocs, iloc
real    :: tolerance = 1.e-8
real(c_double), dimension(:), allocatable :: lons_ref, lats_ref, lons, lats
character(len=200) :: logmessage

!> default value: test passed
test_locations_lonslats_c = 1

!> read reference longitudes and latitudes
f_conf = fckit_configuration(c_conf)
nlocs = f_conf%get_size("reference lons")
call f_conf%get_or_die("reference lons", lons_ref)
call f_conf%get_or_die("reference lats", lats_ref)

!> get longitudes and latitudes from locations object
locs = ufo_locations(c_locs)
allocate(lons(nlocs), lats(nlocs))
call locs%get_lons(lons)
call locs%get_lats(lats)

do iloc = 1, nlocs
  write(logmessage, *) "lons: loc and ref: ", lons(iloc), ", ", lons_ref(iloc)
  call fckit_log%debug(logmessage)
  write(logmessage, *) "lats: loc and ref: ", lats(iloc), ", ", lats_ref(iloc)
  call fckit_log%debug(logmessage)
enddo

!> compare to reference
if(any(abs(lons-lons_ref) > tolerance)) test_locations_lonslats_c = 0
if(any(abs(lats-lats_ref) > tolerance)) test_locations_lonslats_c = 0

deallocate(lons, lats)

end function test_locations_lonslats_c

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------
!> Tests whether ufo::Locations time mask method results
!! are matching references in yaml file
integer function test_locations_timemask_c(c_locs, c_conf) bind(c,name='test_locations_timemask_f90')
use fckit_configuration_module, only: fckit_configuration
use fckit_log_module, only: fckit_log
use datetime_mod
use ufo_locations_mod
implicit none
type(c_ptr), value, intent(in) :: c_locs  !< ufo::Locations object
type(c_ptr), value, intent(in) :: c_conf  !< test configuration (described in LocationsTestParameters

!> local variables
type(ufo_locations) :: locs
type(fckit_configuration) :: f_conf
type(fckit_configuration), allocatable :: testconfigs(:)
character(kind=c_char,len=:), allocatable :: t1str, t2str
logical(kind=c_bool), allocatable :: mask(:)
logical, allocatable :: mask_ref(:)
type(datetime) :: t1, t2
integer :: nlocs, iloc, itest
character(len=200) :: logmessage

!> default value: test passed
test_locations_timemask_c = 1

locs = ufo_locations(c_locs)
f_conf = fckit_configuration(c_conf)
call f_conf%get_or_die("time mask tests", testconfigs)
!> loop over all the tests for time masks
do itest = 1, size(testconfigs)
  call testconfigs(itest)%get_or_die("t1", t1str)
  call testconfigs(itest)%get_or_die("t2", t2str)
  call datetime_create(t1str, t1)
  call datetime_create(t2str, t2)
  call testconfigs(itest)%get_or_die("reference mask", mask_ref)

  !> call ufo::Locations method to compute time mask
  nlocs = locs%nlocs()
  allocate(mask(nlocs))
  call locs%get_timemask(t1, t2, mask)

  write(logmessage, *) "Test ", itest
  call fckit_log%debug(logmessage)
  do iloc = 1, nlocs
    write(logmessage, *) "mask: loc and ref: ", mask(iloc), ", ", mask_ref(iloc)
    call fckit_log%debug(logmessage)
  enddo

  !> compare to reference
  if(any(mask .neqv. mask_ref)) test_locations_timemask_c = 0

  deallocate(mask)

enddo

end function test_locations_timemask_c

end module test_locations
