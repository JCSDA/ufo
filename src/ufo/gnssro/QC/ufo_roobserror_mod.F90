! (C) Copyright 2017-2018 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.  

!> Fortran module to implement RO observational error

module ufo_roobserror_mod
use fckit_configuration_module, only: fckit_configuration 
use kinds
use ufo_geovals_mod
use obsspace_mod
use oops_variables_mod
use missing_values_mod
use gnssro_mod_obserror
use fckit_log_module, only : fckit_log

implicit none
public :: ufo_roobserror, ufo_roobserror_create, ufo_roobserror_delete
public :: ufo_roobserror_prior
public :: max_string
private

! ------------------------------------------------------------------------------
type :: ufo_roobserror
  character(len=max_string)    :: variable
  character(len=max_string)    :: errmodel
  character(len=max_string)    :: rmatrix_filename
  character(len=max_string)    :: err_variable
  type(oops_variables), public :: obsvar
  type(c_ptr)                  :: obsdb
  integer                      :: n_horiz            ! For ROPP-2D multiplier of the number of geovals
end type ufo_roobserror

! ------------------------------------------------------------------------------
contains
! ------------------------------------------------------------------------------

subroutine ufo_roobserror_create(self, obspace, f_conf)
use iso_c_binding
use oops_variables_mod
implicit none
type(ufo_roobserror), intent(inout)   :: self
type(c_ptr),  value,       intent(in) :: obspace
type(fckit_configuration), intent(in) :: f_conf
character(len=:), allocatable         :: str

self%errmodel = "NBAM"
if (f_conf%has("errmodel")) then
   call f_conf%get_or_die("errmodel",str)
   self%errmodel = str
end if

self % rmatrix_filename = ""
if (f_conf % has("rmatrix_filename")) then
   call f_conf % get_or_die("rmatrix_filename", str)
   self % rmatrix_filename = str
end if

self % err_variable = ""
if (f_conf % has("err_variable")) then
   call f_conf % get_or_die("err_variable", str)
   self % err_variable = str
end if

! Get the number of extra geovals as a multiplier, default to 1
self % n_horiz = 1
if (f_conf % has("n_horiz")) then
   call f_conf % get_or_die("n_horiz", self % n_horiz)
end if

self%obsdb      = obspace

end subroutine ufo_roobserror_create

! ------------------------------------------------------------------------------

subroutine ufo_roobserror_delete(self)
implicit none
type(ufo_roobserror), intent(inout) :: self
end subroutine ufo_roobserror_delete

! ------------------------------------------------------------------------------

subroutine ufo_roobserror_prior(self, model_nobs, model_nlevs, air_temperature, &
    geopotential_height)

use fckit_log_module, only : fckit_log
use ufo_utils_mod, only: cmp_strings

implicit none

type(ufo_roobserror), intent(in) :: self
integer, intent(in)              :: model_nobs
integer, intent(in)              :: model_nlevs
real, intent(in)                 :: air_temperature(:,:)
real, intent(in)                 :: geopotential_height(:,:)

integer                        :: nobs
real(kind_real), allocatable   :: obsZ(:), obsLat(:)
integer,         allocatable   :: obsSatid(:)             ! Satellite identifier
integer,         allocatable   :: obsOrigC(:)             ! Originating centre number for each ob
real(kind_real), allocatable   :: obsImpA(:)              ! The observation impact alitude
real(kind_real), allocatable   :: obsImpH(:)              ! The observation impact height
real(kind_real), allocatable   :: obsImpP(:)              ! The observation impact parameter
real(kind_real), allocatable   :: obsGeoid(:)             ! The geoid undulation at the observation location
real(kind_real), allocatable   :: obsLocR(:)              ! The local radius of curvature at the observation location
real(kind_real), allocatable   :: obsValue(:)
real(kind_real), allocatable   :: obsErr(:)
integer(c_int),  allocatable   :: obsSaid(:)
integer(c_int),  allocatable   :: QCflags(:)
real(kind_real)                :: missing
character(max_string)          :: err_msg

missing = missing_value(missing)
nobs    = obsspace_get_nlocs(self%obsdb)
allocate(QCflags(nobs))
allocate(obsErr(nobs))
QCflags(:)  = 0

if (model_nobs /= nobs * self % n_horiz) then
  write(err_msg, '(A,2I8)') 'nobs from model and observations must be equal', nobs, model_nobs
  call abor1_ftn(err_msg)
end if

! read QC flags
call obsspace_get_db(self%obsdb, "FortranQC", trim(self%variable),QCflags )

!-------------------------------
select case (trim(self%variable))

!-------------------------------
case ("bending_angle")

  allocate(obsImpP(nobs))
  allocate(obsGeoid(nobs))
  allocate(obsLocR(nobs))
  allocate(obsImpH(nobs))
  allocate(obsImpA(nobs))
  call obsspace_get_db(self%obsdb, "MetaData", "impact_parameter", obsImpP)
  call obsspace_get_db(self%obsdb, "MetaData", "geoid_height_above_reference_ellipsoid",obsGeoid)
  call obsspace_get_db(self%obsdb, "MetaData", "earth_radius_of_curvature", obsLocR)
  obsImpH(:) = obsImpP(:) - obsLocR(:)
  obsImpA(:) = obsImpP(:) - obsGeoid(:) - obsLocR(:)

  select case (trim(self%errmodel))
  case ("NBAM")
    allocate(obsSaid(nobs))
    allocate(obsLat(nobs))
    call obsspace_get_db(self%obsdb, "MetaData", "occulting_sat_id", obsSaid)
    call obsspace_get_db(self%obsdb, "MetaData", "latitude", obsLat)
    call bending_angle_obserr_NBAM(obsLat, obsImpA, obsSaid, nobs, obsErr, QCflags, missing)
    write(err_msg,*) "ufo_roobserror_mod: setting up bending_angle obs error with NBAM method"
    call fckit_log%info(err_msg)
    deallocate(obsSaid)
    deallocate(obsLat)
    ! update obs error
    call obsspace_put_db(self%obsdb, "FortranERR", trim(self%variable), obsErr)

  case ("ECMWF")
    allocate(obsValue(nobs))
    call obsspace_get_db(self%obsdb, "ObsValue", "bending_angle", obsValue)
    call bending_angle_obserr_ECMWF(obsImpA, obsValue, nobs, obsErr, QCflags, missing)
    write(err_msg,*) "ufo_roobserror_mod: setting up bending_angle obs error with ECMWF method"
    call fckit_log%info(err_msg)
    deallocate(obsValue)
    ! update obs error
    call obsspace_put_db(self%obsdb, "FortranERR", trim(self%variable), obsErr)
  case ("NRL")
    allocate(obsValue(nobs))
    allocate(obsLat(nobs))
    call obsspace_get_db(self%obsdb, "ObsValue", "bending_angle", obsValue)
    call obsspace_get_db(self%obsdb, "MetaData", "latitude", obsLat)
    call bending_angle_obserr_NRL(obsLat, obsImpA, obsValue, nobs, obsErr, QCflags, missing)
    write(err_msg,*) "ufo_roobserror_mod: setting up bending_angle obs error with NRL method"
    call fckit_log%info(err_msg)
    deallocate(obsValue)
    deallocate(obsLat)
    ! update obs error
    call obsspace_put_db(self%obsdb, "FortranERR", trim(self%variable), obsErr)

  case ("MetOffice")
    write(err_msg,*) "ufo_roobserror_mod: setting up bending_angle obs error with the Met Office method"
    call fckit_log%info(err_msg)
    if (cmp_strings(self % rmatrix_filename, "")) then
      err_msg = "If you choose the Met Office method, then you must specify rmatrix_filename"
      call abor1_ftn(err_msg)
    end if
    allocate(obsSatid(nobs))
    allocate(obsOrigC(nobs))
    allocate(obsValue(nobs))
    call obsspace_get_db(self%obsdb, "MetaData", "occulting_sat_id", obsSatid)
    call obsspace_get_db(self%obsdb, "MetaData", "originating_center", obsOrigC)
    call obsspace_get_db(self%obsdb, "ObsValue", "bending_angle", obsValue)
    call obsspace_get_db(self%obsdb, "ObsError", trim(self%variable), obsErr)
    if (self % err_variable == "latitude") then
      allocate(obsLat(nobs))
      call obsspace_get_db(self%obsdb, "MetaData", "latitude", obsLat)
      call gnssro_obserr_latitude(nobs, self % rmatrix_filename, obsSatid, obsOrigC, obsLat, obsImpH, &
                                  obsValue, obsErr, QCflags, missing)
      deallocate(obsLat)
    else if (self % err_variable == "average_temperature") then
      call gnssro_obserr_avtemp(nobs, self % n_horiz, self % rmatrix_filename, obsSatid, obsOrigC, &
                                model_nlevs, air_temperature, geopotential_height, obsImpH, obsValue, &
                                obsErr, QCflags, missing)
    else
      err_msg = "The error variable should be either 'latitude' or 'average_temperature', but you gave " // &
                trim(self % err_variable)
      call fckit_log%info(err_msg)
    end if
    deallocate(obsValue)
    deallocate(obsSatid)
    deallocate(obsOrigC)
    ! up date obs error
    call obsspace_put_db(self%obsdb, "FortranERR", trim(self%variable), obsErr)

  case default
    write(err_msg,*) "ufo_roobserror_mod: bending_angle error model must be NBAM, ECMWF, NRL or MetOffice"
    call fckit_log%info(err_msg)
  end select
  deallocate(obsImpP)
  deallocate(obsGeoid)
  deallocate(obsLocR)
  deallocate(obsImpH)
  deallocate(obsImpA)

!-------------------------------
case ("refractivity")

  select case (trim(self%errmodel))

  case ("NBAM")

    allocate(obsZ(nobs))
    allocate(obsLat(nobs))
    call obsspace_get_db(self%obsdb, "MetaData", "altitude",  obsZ) 
    call obsspace_get_db(self%obsdb, "MetaData", "latitude", obsLat)
    call refractivity_obserr_NBAM(obsLat, obsZ, nobs, obsErr, QCflags, missing)
    write(err_msg,*) "ufo_roobserror_mod: setting up refractivity obs error with NBAM method" 
    call fckit_log%info(err_msg)
    deallocate(obsZ)
    deallocate(obsLat)  
    ! up date obs error
    call obsspace_put_db(self%obsdb, "FortranERR", trim(self%variable), obsErr)

  case ("ECMWF")
     write(err_msg,*) "ufo_roobserror_mod: ECMWF refractivity error model is not available now"
     call fckit_log%info(err_msg)

  case default
    write(err_msg,*) "ufo_roobserror_mod: only NBAM refractivity model is available now"
    call fckit_log%info(err_msg)
  end select

case default
  call abor1_ftn("ufo_roobserror_prior: variable has to be bending_angle or refractivity")
end select

deallocate(QCflags)
deallocate(obsErr)

end subroutine ufo_roobserror_prior

end module ufo_roobserror_mod
