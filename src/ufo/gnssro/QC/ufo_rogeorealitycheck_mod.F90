! (C) Copyright 2017-2018 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.  

!> Fortran module to implement RO geophysical reality check

module ufo_rogeorealitycheck_mod
use iso_c_binding
use kinds
use ufo_geovals_mod
use obsspace_mod
use config_mod
use missing_values_mod

implicit none
public :: ufo_rogeorealitycheck, ufo_rogeorealitycheck_create, ufo_rogeorealitycheck_delete
public :: ufo_rogeorealitycheck_prior, ufo_rogeorealitycheck_post
public :: max_string
private
integer, parameter :: max_string=99 

! ------------------------------------------------------------------------------
type :: ufo_rogeorealitycheck
  character(len=max_string) :: variable
  integer(c_int)  :: ro_top_meter
  integer(c_int)  :: ro_bot_meter
  type(c_ptr)     :: obsdb
end type ufo_rogeorealitycheck

! ------------------------------------------------------------------------------
contains
! ------------------------------------------------------------------------------

subroutine ufo_rogeorealitycheck_create(self, obspace, conf)
implicit none
type(ufo_rogeorealitycheck), intent(inout) :: self
type(c_ptr),  value,       intent(in)    :: obspace
type(c_ptr),               intent(in)    :: conf

self%variable      = config_get_string(conf, max_string, "variable", "refractivity")
self%ro_bot_meter  = config_get_int(conf, "ro_bot_meter", 0)
if (trim(self%variable) .eq. "refractivity" ) then
  self%ro_top_meter  = config_get_int(conf, "ro_top_meter", 30000)
else
  self%ro_top_meter  = config_get_int(conf, "ro_top_meter", 50000)
endif
if (self%ro_top_meter<=self%ro_bot_meter .or. self%ro_bot_meter<0.0_kind_real)  then
   call abor1_ftn("ufo_rogeorealitycheck_create: Error RO RANGE threshold")
endif 
self%obsdb = obspace

end subroutine ufo_rogeorealitycheck_create

! ------------------------------------------------------------------------------

subroutine ufo_rogeorealitycheck_delete(self)
implicit none
type(ufo_rogeorealitycheck), intent(inout) :: self
end subroutine ufo_rogeorealitycheck_delete

! ------------------------------------------------------------------------------

subroutine ufo_rogeorealitycheck_prior(self, geovals)
use fckit_log_module, only : fckit_log
implicit none
type(ufo_rogeorealitycheck), intent(in) :: self
type(ufo_geovals),      intent(in)    :: geovals
integer                           :: iobs, jobs,ireject
integer, parameter                :: ro_geocheck_flag   = 76
integer, parameter                :: ro_rangecheck_flag = 77
real(kind_real),    allocatable   :: yalt(:), yimpar(:), ygeoid(:), yearthr(:)
real(kind_real)                   :: yimph !! impact height
integer(c_int32_t), allocatable   :: flags(:)
real(kind_real)                   :: missing
character(max_string)             :: err_msg

missing = missing_value(missing)
iobs    = obsspace_get_nlocs(self%obsdb)
allocate(yalt(iobs))
allocate(ygeoid(iobs))
allocate(yearthr(iobs))
allocate(yimpar(iobs))
allocate(flags(iobs))
flags(:)  = 0
ireject   = 0

call obsspace_get_db(self%obsdb, "", "impact_parameter", yimpar)
call obsspace_get_db(self%obsdb, "", "geoid_height_above_reference_ellipsoid",ygeoid)
call obsspace_get_db(self%obsdb, "", "altitude",  yalt)
call obsspace_get_db(self%obsdb, "", "earth_radius_of_curvature", yearthr)
call obsspace_get_db(self%obsdb, "FortranQC", trim(self%variable),flags )

do jobs = 1, iobs

   yimph = yimpar(jobs) - ygeoid(jobs) - yearthr(jobs)
!   if (    yearthr(jobs)<=6250000 .or. yearthr(jobs)>=6450000  &
!      .or. ygeoid(jobs)<=-200 .or. ygeoid(jobs)>=200           &
!      .or. yimph<=0 .or. yimph>=90000 .and. flags(jobs)<=0.01 ) then
   if (     yimph<=0 .or. yimph>=90000 .and. flags(jobs)<=0.01 ) then
      flags(jobs) = ro_geocheck_flag
      ireject     = ireject + 1
   endif 

!   if ( (yalt(jobs)>self%ro_top_meter .or. yalt(jobs)<self%ro_bot_meter)   &
!        .and. flags(jobs)<=0.01 ) then
!      flags(jobs) = ro_rangecheck_flag
!      ireject     = ireject + 1
!   endif

enddo

write(err_msg,*)'ufo_ufo_rogeorealitycheck_prior: ',ireject, " out of ", iobs," are rejected"
call fckit_log%info(err_msg)

call obsspace_put_db(self%obsdb, "FortranQC", trim(self%variable), flags)

deallocate(yalt)
deallocate(ygeoid)
deallocate(yearthr)
deallocate(yimpar)
deallocate(flags)

end subroutine ufo_rogeorealitycheck_prior

! ------------------------------------------------------------------------------
subroutine ufo_rogeorealitycheck_post(self, hofx)
implicit none
type(ufo_rogeorealitycheck), intent(in) :: self
real(c_double),  intent(in)           :: hofx(:)
end subroutine ufo_rogeorealitycheck_post

end module ufo_rogeorealitycheck_mod
