! (C) Copyright 2017-2018 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.  

!> Fortran module to implement RO observational error

module ufo_roobserror_mod
use iso_c_binding
use kinds
use ufo_geovals_mod
use obsspace_mod
use config_mod
use missing_values_mod
use gnssro_mod_obserror

implicit none
public :: ufo_roobserror, ufo_roobserror_create, ufo_roobserror_delete
public :: ufo_roobserror_prior, ufo_roobserror_post
public :: max_string
private
integer, parameter :: max_string=99 

! ------------------------------------------------------------------------------
type :: ufo_roobserror
  character(len=max_string) :: variable
  character(len=max_string) :: errmodel
  character(len=max_string) :: qcname
  type(c_ptr)               :: obsdb
end type ufo_roobserror

! ------------------------------------------------------------------------------
contains
! ------------------------------------------------------------------------------

subroutine ufo_roobserror_create(self, obspace, conf)
implicit none
type(ufo_roobserror), intent(inout) :: self
type(c_ptr),  value,       intent(in)    :: obspace
type(c_ptr),               intent(in)    :: conf

self%variable   = config_get_string(conf, max_string, "variable", "bending_angle")
self%errmodel   = config_get_string(conf, max_string, "errmodel", "GSI")
self%obsdb      = obspace
self%qcname     = config_get_string(conf, max_string, "QCname")
end subroutine ufo_roobserror_create

! ------------------------------------------------------------------------------

subroutine ufo_roobserror_delete(self)
implicit none
type(ufo_roobserror), intent(inout) :: self
end subroutine ufo_roobserror_delete

! ------------------------------------------------------------------------------

subroutine ufo_roobserror_prior(self)
use fckit_log_module, only : fckit_log
implicit none
type(ufo_roobserror), intent(in) :: self
integer                          :: nobs
real(kind_real),    allocatable   :: obsZ(:), obsLat(:)
real(kind_real),    allocatable   :: obsImpH(:),obsImpP(:),obsGeoid(:),obsLocR(:)
real(kind_real),    allocatable   :: obsValue(:)
real(kind_real),    allocatable   :: obsErr(:)
integer(c_int), allocatable   :: obsSaid(:)
integer(c_int), allocatable   :: QCflags(:)
real(kind_real)                   :: missing
character(max_string)             :: err_msg

missing = missing_value(missing)
nobs    = obsspace_get_nlocs(self%obsdb)
allocate(QCflags(nobs))
allocate(obsErr(nobs))
QCflags(:)  = 0

! read QC flags
call obsspace_get_db(self%obsdb, self%qcname, trim(self%variable),QCflags )

!-------------------------------
select case (trim(self%variable))

!-------------------------------
case ("bending_angle")

  allocate(obsImpP(nobs))
  allocate(obsGeoid(nobs))
  allocate(obsLocR(nobs))
  allocate(obsImpH(nobs))
  call obsspace_get_db(self%obsdb, "", "impact_parameter", obsImpP)
  call obsspace_get_db(self%obsdb, "", "geoid_height_above_reference_ellipsoid",obsGeoid)
  call obsspace_get_db(self%obsdb, "", "earth_radius_of_curvature", obsLocR)
  obsImpH(:) = obsImpP(:) - obsGeoid(:) - obsLocR(:)

  select case (trim(self%errmodel))
  case ("GSI")
    allocate(obsSaid(nobs))
    call obsspace_get_db(self%obsdb, " ", "occulting_sat_id", obsSaid)
    call bending_angle_obserr_GSI(obsLat, obsImpH, obsSaid, nobs, obsErr, QCflags, missing)
    write(err_msg,*) "ufo_roobserror_mod: setting up bending_angle obs error with GSI method"
    deallocate(obsSaid)

  case ("ROPP")
    allocate(obsValue(nobs))
    call obsspace_get_db(self%obsdb, "ObsValue", "bending_angle", obsValue)
    call bending_angle_obserr_ROPP(obsImpH, obsValue, nobs, obsErr, QCflags, missing)
    write(err_msg,*) "ufo_roobserror_mod: setting up bending_angle obs error with ROPP method"
    deallocate(obsValue)

  case default
    write(err_msg,*) "ufo_roobserror_mod: bending_angle error model must be GSI or ROPP"
    call fckit_log%info(err_msg) 

  end select
  deallocate(obsImpP)
  deallocate(obsGeoid)
  deallocate(obsLocR)
  deallocate(obsImpH)

!-------------------------------
case ("refractivity")

  select case (trim(self%errmodel))

  case ("GSI")
    allocate(obsZ(nobs))
    allocate(obsLat(nobs))
    call obsspace_get_db(self%obsdb, "", "altitude",  obsZ) 
    call refractivity_obserr_GSI(obsLat, obsZ, nobs, obsErr, QCflags, missing)
    write(err_msg,*) "ufo_roobserror_mod: setting up refractivity obs error with GSI method" 
    deallocate(obsZ)
    deallocate(obsLat)  

  case ("ROPP")
     write(err_msg,*) "ufo_roobserror_mod: ROPP refractivity error model is not avaiable now"
     call fckit_log%info(err_msg)

  case default
    write(err_msg,*) "ufo_roobserror_mod: only GSI refractivity model is avaiable now"
    call fckit_log%info(err_msg)
  end select

case default
  call abor1_ftn("ufo_roobserror_prior: variable has to be bending_angle or refractivity")
end select

! up date obs error
call obsspace_put_db(self%obsdb, "EffectiveError", trim(self%variable), obsErr)

deallocate(QCflags)
deallocate(obsErr)

end subroutine ufo_roobserror_prior

! ------------------------------------------------------------------------------
subroutine ufo_roobserror_post(self)
implicit none
type(ufo_roobserror), intent(in) :: self
end subroutine ufo_roobserror_post

end module ufo_roobserror_mod
