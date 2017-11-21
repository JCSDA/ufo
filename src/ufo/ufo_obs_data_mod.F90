! (C) Copyright 2017 UCAR
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 

!> Fortran module to handle radiosondeentional observations

module ufo_obs_data_mod
  use iso_c_binding
  use fckit_log_module, only : fckit_log

  use ufo_obs_data_basis_mod

  use radDiag_mod,  only: RadDiag
  use radDiag_mod,  only: RadDiag_read
  use raobDiag_mod, only: RaobDiag
  use raobDiag_mod, only: RaobDiag_read

  implicit none

  type, extends(BasisObsData) :: Obs_Data
  contains
    ! Implementation for the deferred procedure in Basis
    procedure :: Setup
    procedure :: Delete
  end type Obs_Data

  type(RadDiag),  pointer, save::  Radiance
  type(RaobDiag), pointer, save::  Radiosonde
contains

  subroutine Setup(self, filein,obstype,nobs,nlocs)
    class(Obs_Data), intent(inout) :: self
    character(len=*), intent(in)    :: filein
    character(len=*), intent(in)    :: obstype
    integer(c_int),   intent(inout) :: nobs
    integer(c_int),   intent(inout) :: nlocs

    select case(trim(obstype))
      case("Radiance")
        allocate(Radiance)
        self%Obspoint => Radiance
        call SetupRadiance(self, Radiance, filein,nobs,nlocs)
      case("Radiosonde")
        allocate(Radiosonde)
        self%Obspoint => Radiosonde
        call SetupRaob(self, Radiosonde, filein,nobs,nlocs)
    end select

  end subroutine Setup

  subroutine SetupRadiance(self, mytype, filein,nobs,nlocs)
    class(Obs_Data), intent(inout) :: self
    type(RadDiag), pointer, intent(inout)    :: mytype
    character(len=*), intent(in)    :: filein
    integer(c_int),   intent(inout) :: nobs
    integer(c_int),   intent(inout) :: nlocs

    character(len=*),parameter:: myname_=myname//"::SetupRadiance"
    character(len=255) :: record
    integer :: failed

    call radDiag_read(mytype,filein,'Radiance',nobs,nlocs)
    self%Nobs = nobs
    self%Nlocs= nlocs

    failed=0
    if(failed==0 .and. size(mytype%datafix(:)%Lat)==nlocs) then
       allocate(self%lat(nlocs))
       self%lat(:) = mytype%datafix(:)%Lat
    else
       failed=1
    endif
    if(failed==0 .and. size(mytype%datafix(:)%Lon)==nlocs) then
       allocate(self%lon(nlocs))
       self%lon(:) = mytype%datafix(:)%Lon
    else
       failed=2
    endif
    if(failed==0 .and. size(mytype%datafix(:)%obstime)==nlocs) then
       allocate(self%time(nlocs))
       self%time(:) = mytype%datafix(:)%obstime
    else
       failed=3
    endif
    if(failed==0)then
      write(record,*)myname_,': allocated/assinged obs-data'
      call fckit_log%info(record)
    else
      write(record,*)myname_,': failed allocation/assignment of obs-data, ier: ', failed
      call fckit_log%info(record)
      ! should exit in error here
    endif

  end subroutine SetupRadiance

  subroutine SetupRaob(self, mytype, filein,nobs,nlocs)
    class(Obs_Data), intent(inout) :: self
    type(RaobDiag), pointer, intent(inout)   :: mytype
    character(len=*), intent(in)    :: filein
    integer(c_int),   intent(inout) :: nobs
    integer(c_int),   intent(inout) :: nlocs

    character(len=*),parameter:: myname_=myname//"::SetupRaob"
    character(len=255) :: record
    integer :: failed 

    call raobDiag_read(mytype,filein,'Radiosonde',nobs,nlocs)
    self%Nobs = nobs
    self%Nlocs= nlocs

    failed=0
    if(failed==0 .and. size(mytype%mass(:)%Longitude)==nlocs) then
       allocate(self%lon(nlocs))
       self%lon(:) = mytype%mass(:)%Longitude
    else
       failed=1
    endif
    if(failed==0 .and. size(mytype%mass(:)%Latitude) ==nlocs) then
       allocate(self%lat(nlocs))
       self%lat(:) = mytype%mass(:)%Latitude
    else
       failed=2
    endif
    if(failed==0 .and. size(mytype%mass(:)%Pressure) ==nlocs) then
       allocate(self%lev(nlocs))
       self%lev(:) = mytype%mass(:)%Pressure
    else
       failed=3
    endif
    if(failed==0 .and. size(mytype%mass(:)%Time) ==nlocs) then
       allocate(self%time(nlocs))
       self%time(:) = mytype%mass(:)%Time
    else
       failed=4
    endif
    if(failed==0)then
      write(record,*)myname_,': allocated/assinged obs-data'
      call fckit_log%info(record)
    else
      write(record,*)myname_,': failed allocation/assignment of obs-data, ier: ', failed
      call fckit_log%info(record)
      ! should exit in error here
    endif

  end subroutine SetupRaob

  subroutine Delete(self)
    class(Obs_Data), intent(inout) :: self
  end subroutine Delete

  function vname2vmold_(vname) result(obsmold_)
    implicit none
    class(BasisObsData),pointer:: obsmold_
    character(len=*),intent(in):: vname

    character(len=*),parameter:: myname_=myname//"::vname2vmold_"

    obsmold_=>NUll()
    select case(trim(vname))
      case("Radiance");   obsmold_%Obspoint => Radiance
      case("Radiosonde"); obsmold_%Obspoint => Radiosonde
    end select
  end function vname2vmold_

end module ufo_obs_data_mod
