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
use obs_variables_mod
use missing_values_mod
use gnssro_mod_obserror
use fckit_log_module, only : fckit_log
use fckit_exception_module, only: fckit_exception
use gnssro_mod_constants, only: max_string
use iso_c_binding, only: c_ptr, c_int, c_size_t

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
  character(len=max_string)    :: average_temperature_name
  logical                      :: verbose_output       ! Whether to give extra output messages
  type(obs_variables), public :: obsvar
  type(c_ptr)                  :: obsdb
  integer                      :: n_horiz              ! For ROPP-2D multiplier of the number of geovals
  logical                      :: allow_extrapolation  ! Allow errors to be extrapolated outside of range?
  logical                      :: use_profile          ! Use a single profile to give errors?
end type ufo_roobserror

! ------------------------------------------------------------------------------
contains
! ------------------------------------------------------------------------------

subroutine ufo_roobserror_create(self, obspace, f_conf)
use iso_c_binding
use obs_variables_mod
implicit none
type(ufo_roobserror), intent(inout)   :: self
type(c_ptr),  value,       intent(in) :: obspace
type(fckit_configuration), intent(in) :: f_conf
character(len=:), allocatable         :: str
character(len=800)                    :: message    ! Message to be output

self%errmodel = "NBAM"
if (f_conf%has("errmodel")) then
   call f_conf%get_or_die("errmodel",str)
   self%errmodel = str
end if
write(message,*) 'errmodel = ', trim(self % errmodel)
call fckit_log%debug(message)

self % rmatrix_filename = ""
if (f_conf % has("rmatrix_filename")) then
   call f_conf % get_or_die("rmatrix_filename", str)
   self % rmatrix_filename = str
end if
write(message,*) 'rmatrix_filename = ', trim(self % rmatrix_filename)
call fckit_log%debug(message)

self % err_variable = ""
if (f_conf % has("err_variable")) then
   call f_conf % get_or_die("err_variable", str)
   self % err_variable = str
end if
write(message,*) 'err_variable = ', trim(self % err_variable)
call fckit_log%debug(message)

! Get the number of extra geovals as a multiplier, default to 1
self % n_horiz = 1
if (f_conf % has("n_horiz")) then
   call f_conf % get_or_die("n_horiz", self % n_horiz)
end if
write(message,*) 'n_horiz = ', self % n_horiz
call fckit_log%debug(message)

self % allow_extrapolation = .false.
if (f_conf % has("allow extrapolation")) then
   call f_conf % get_or_die("allow extrapolation", self % allow_extrapolation)
end if
write(message,*) 'allow_extrapolation = ', self % allow_extrapolation
call fckit_log%debug(message)

self % use_profile = .false.
if (f_conf % has("use profile")) then
   call f_conf % get_or_die("use profile", self % use_profile)
end if
write(message,*) 'use_profile = ', self % use_profile
call fckit_log%debug(message)

self % verbose_output = .false.
if (f_conf % has("verbose output")) then
   call f_conf % get_or_die("verbose output", self % verbose_output)
end if
write(message,*) 'verbose_output = ', self % verbose_output
call fckit_log%debug(message)

self % average_temperature_name = ""
if (f_conf % has("average temperature name")) then
   call f_conf % get_or_die("average temperature name", str)
   self % average_temperature_name = str
end if
write(message,*) 'average_temperature_name = ', trim(self % average_temperature_name)
call fckit_log%debug(message)

self%obsdb      = obspace

end subroutine ufo_roobserror_create

! ------------------------------------------------------------------------------

subroutine ufo_roobserror_delete(self)
implicit none
type(ufo_roobserror), intent(inout) :: self
end subroutine ufo_roobserror_delete

! ------------------------------------------------------------------------------

subroutine ufo_roobserror_prior(self)

use fckit_log_module, only : fckit_log
use ufo_utils_mod, only: cmp_strings, sort_and_unique

implicit none

type(ufo_roobserror), intent(in) :: self

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
real(kind_real), allocatable   :: averageTemp(:)          ! The average model temperature
integer(c_int),  allocatable   :: obsSaid(:)
integer(c_int),  allocatable   :: QCflags(:)
real(kind_real)                :: missing
character(max_string)          :: err_msg
integer                        :: iob                     ! Loop variable, observation number
integer(c_size_t), allocatable :: record_number(:)        ! Number used to identify unique profiles in the data
integer, allocatable           :: sort_order(:)           ! Sorting order for the profile numbers
integer, allocatable           :: unique(:)               ! Set of unique profile numbers

missing = missing_value(missing)
nobs    = obsspace_get_nlocs(self%obsdb)
allocate(QCflags(nobs))
allocate(obsErr(nobs))
QCflags(:)  = 0

! read QC flags
call obsspace_get_db(self%obsdb, "FortranQC", trim(self%variable),QCflags )

!-------------------------------
select case (trim(self%variable))

!-------------------------------
case ("bendingAngle")

  allocate(obsImpP(nobs))
  allocate(obsGeoid(nobs))
  allocate(obsLocR(nobs))
  allocate(obsImpH(nobs))
  allocate(obsImpA(nobs))
  call obsspace_get_db(self%obsdb, "MetaData", "impactParameterRO", obsImpP)
  call obsspace_get_db(self%obsdb, "MetaData", "geoidUndulation",obsGeoid)
  call obsspace_get_db(self%obsdb, "MetaData", "earthRadiusCurvature", obsLocR)
  obsImpH(:) = obsImpP(:) - obsLocR(:)
  obsImpA(:) = obsImpP(:) - obsGeoid(:) - obsLocR(:)

  allocate(record_number(nobs))
  if (self % use_profile) then
    ! Get the record numbers and the set of unique profile numbers
    call obsspace_get_recnum(self % obsdb, record_number)
    call sort_and_unique(nobs, obsImpP, record_number, sort_order, unique)
  else
    allocate(unique(1:nobs), sort_order(1:nobs))
    do iob = 1, nobs
      record_number(iob) = iob
      sort_order(iob) = iob
      unique(iob) = iob
    end do
  end if

  select case (trim(self%errmodel))
  case ("NBAM")
    allocate(obsSaid(nobs))
    allocate(obsLat(nobs))
    call obsspace_get_db(self%obsdb, "MetaData", "satelliteIdentifier", obsSaid)
    call obsspace_get_db(self%obsdb, "MetaData", "latitude", obsLat)
    call bending_angle_obserr_NBAM(obsLat, obsImpH, obsSaid, nobs, obsErr, QCflags, missing)
    write(err_msg,*) "ufo_roobserror_mod: setting up bending_angle obs error with NBAM method"
    call fckit_log%debug(err_msg)
    deallocate(obsSaid)
    deallocate(obsLat)
    ! update obs error
    call obsspace_put_db(self%obsdb, "FortranERR", trim(self%variable), obsErr)

  case ("ECMWF")
    allocate(obsValue(nobs))
    call obsspace_get_db(self%obsdb, "ObsValue", "bendingAngle", obsValue)
    call bending_angle_obserr_ECMWF(obsImpA, obsValue, nobs, obsErr, QCflags, missing)
    write(err_msg,*) "ufo_roobserror_mod: setting up bending_angle obs error with ECMWF method"
    call fckit_log%debug(err_msg)
    deallocate(obsValue)
    ! update obs error
    call obsspace_put_db(self%obsdb, "FortranERR", trim(self%variable), obsErr)
  case ("NRL")
    allocate(obsValue(nobs))
    allocate(obsLat(nobs))
    call obsspace_get_db(self%obsdb, "ObsValue", "bendingAngle", obsValue)
    call obsspace_get_db(self%obsdb, "MetaData", "latitude", obsLat)
    call bending_angle_obserr_NRL(obsLat, obsImpA, obsValue, nobs, obsErr, QCflags, missing)
    write(err_msg,*) "ufo_roobserror_mod: setting up bending_angle obs error with NRL method"
    call fckit_log%debug(err_msg)
    deallocate(obsValue)
    deallocate(obsLat)
    ! update obs error
    call obsspace_put_db(self%obsdb, "FortranERR", trim(self%variable), obsErr)

  case ("MetOffice")
    write(err_msg,*) "ufo_roobserror_mod: setting up bending_angle obs error with the Met Office method"
    call fckit_log%debug(err_msg)
    if (cmp_strings(self % rmatrix_filename, "")) then
      err_msg = "If you choose the Met Office method, then you must specify rmatrix_filename"
      call abor1_ftn(err_msg)
    end if
    allocate(obsSatid(nobs))
    allocate(obsOrigC(nobs))
    allocate(obsValue(nobs))
    call obsspace_get_db(self%obsdb, "MetaData", "satelliteIdentifier", obsSatid)
    call obsspace_get_db(self%obsdb, "MetaData", "dataProviderOrigin", obsOrigC)
    call obsspace_get_db(self%obsdb, "ObsValue", "bendingAngle", obsValue)
    call obsspace_get_db(self%obsdb, "ObsError", trim(self%variable), obsErr)
    if (self % err_variable == "latitude") then
      allocate(obsLat(nobs))
      call obsspace_get_db(self%obsdb, "MetaData", "latitude", obsLat)
      call gnssro_obserr_latitude(nobs, self % rmatrix_filename, obsSatid, obsOrigC, obsLat, obsImpH, &
                                  obsValue, obsErr, QCflags, missing, self % allow_extrapolation, &
                                  record_number, sort_order, unique, self % verbose_output)
      deallocate(obsLat)
    else if (self % err_variable == "average_temperature") then
      allocate(averageTemp(nobs))
      if (self % average_temperature_name == "") then
        call fckit_exception % throw("When using average_temperature error model," &
          // " you must specify an average_temperature_name")
      else
        call obsspace_get_db(self%obsdb, "MetaData", self % average_temperature_name, averageTemp)
      end if
      call gnssro_obserr_avtemp(nobs, self % n_horiz, self % rmatrix_filename, obsSatid, obsOrigC, &
                                obsImpH, obsValue, averageTemp, &
                                obsErr, QCflags, missing, self % allow_extrapolation, record_number, &
                                sort_order, unique, self % verbose_output)
    else
      err_msg = "The error variable should be either 'latitude' or 'average_temperature', but you gave " // &
                trim(self % err_variable)
      call abor1_ftn(err_msg)
    end if
    deallocate(obsValue)
    deallocate(obsSatid)
    deallocate(obsOrigC)
    ! up date obs error
    call obsspace_put_db(self%obsdb, "FortranERR", trim(self%variable), obsErr)

  case default
    write(err_msg,*) "ufo_roobserror_mod: bending_angle error model must be NBAM, NBAM_GEOS, ECMWF, NRL or MetOffice"
    call abor1_ftn(err_msg)
  end select
  deallocate(obsImpP)
  deallocate(obsGeoid)
  deallocate(obsLocR)
  deallocate(obsImpH)
  deallocate(obsImpA)

!-------------------------------
case ("atmosphericRefractivity")

  select case (trim(self%errmodel))

  case ("NCEP")

    allocate(obsZ(nobs))
    allocate(obsLat(nobs))
    call obsspace_get_db(self%obsdb, "MetaData", "height",  obsZ) 
    call obsspace_get_db(self%obsdb, "MetaData", "latitude", obsLat)
    call refractivity_obserr_NCEP(obsLat, obsZ, nobs, obsErr, QCflags, missing)
    write(err_msg,*) "ufo_roobserror_mod: setting up refractivity obs error with NCEP method" 
    call fckit_log%debug(err_msg)
    deallocate(obsZ)
    deallocate(obsLat)  
    ! up date obs error
    call obsspace_put_db(self%obsdb, "FortranERR", trim(self%variable), obsErr)

  case ("ECMWF")
     write(err_msg,*) "ufo_roobserror_mod: ECMWF refractivity error model is not available now"
     call fckit_log%info(err_msg)

  case default
     write(err_msg,*) "ufo_roobserror_mod: only NCEP refractivity model is available now"
     call abor1_ftn(err_msg)
  end select

case default
  call abor1_ftn("ufo_roobserror_prior: variable has to be bending_angle or refractivity")
end select

deallocate(QCflags)
deallocate(obsErr)

end subroutine ufo_roobserror_prior

end module ufo_roobserror_mod
