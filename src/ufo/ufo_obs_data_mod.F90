! (C) Copyright 2017 UCAR
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 

!> Fortran module to handle radiosondeentional observations

module ufo_obs_data_mod
  use iso_c_binding

  use ufo_obs_data_basis_mod

  use radDiag_mod,  only: RadDiag
  use radDiag_mod,  only: RadDiag_read
  use raobDiag_mod, only: RaobDiag
  use raobDiag_mod, only: RaobDiag_read
  use aodDiag_mod,  only: AodDiag
  use aodDiag_mod,  only: AodDiag_read
  

  implicit none

  type, extends(BasisObsData) :: Obs_Data
  contains
    ! Implementation for the deferred procedure in Basis
    procedure :: Setup
    procedure :: Delete
    procedure :: GetLocs
  end type Obs_Data

  type(RadDiag),  pointer, save::  Radiance
  type(RaobDiag), pointer, save::  Radiosonde
  type(AodDiag),  pointer, save::  Aod

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
      case("Aod")
        allocate(Aod)
        self%Obspoint => Aod
        call SetupAod(self, Aod, filein,nobs,nlocs)
    end select

  end subroutine Setup

  subroutine SetupRadiance(self, mytype, filein,nobs,nlocs)
    class(Obs_Data), intent(inout) :: self
    type(RadDiag), pointer, intent(inout)    :: mytype
    character(len=*), intent(in)    :: filein
    integer(c_int),   intent(inout) :: nobs
    integer(c_int),   intent(inout) :: nlocs

    character(len=*),parameter:: myname_=myname//"::SetupRadiance"

    call radDiag_read(mytype,filein,'Radiance',nobs,nlocs)
    self%Nobs = nobs
    self%Nlocs= nlocs

  end subroutine SetupRadiance

  subroutine SetupRaob(self, mytype, filein,nobs,nlocs)
    class(Obs_Data), intent(inout) :: self
    type(RaobDiag), pointer, intent(inout)   :: mytype
    character(len=*), intent(in)    :: filein
    integer(c_int),   intent(inout) :: nobs
    integer(c_int),   intent(inout) :: nlocs

    character(len=*),parameter:: myname_=myname//"::SetupRaob"

    call raobDiag_read(mytype,filein,'Radiosonde',nobs,nlocs)
    self%Nobs = nobs
    self%Nlocs= nlocs

  end subroutine SetupRaob

  subroutine SetupAod(self, mytype, filein,nobs,nlocs)
    class(Obs_Data), intent(inout) :: self
    type(AodDiag), pointer, intent(inout)    :: mytype
    character(len=*), intent(in)    :: filein
    integer(c_int),   intent(inout) :: nobs
    integer(c_int),   intent(inout) :: nlocs

    character(len=*),parameter:: myname_=myname//"::SetupAod"

    call aodDiag_read(mytype,filein,'Aod',nobs,nlocs)
    self%Nobs = nobs
    self%Nlocs= nlocs

  end subroutine SetupAod

  subroutine Delete(self)
    class(Obs_Data), intent(inout) :: self
  end subroutine Delete

  subroutine GetLocs(self, nlocs, locs)
    class(Obs_Data), intent(in) :: self
    integer, intent(in) :: nlocs
    type(ufo_locs), intent(inout) :: locs
  end subroutine GetLocs

  function vname2vmold_(vname) result(obsmold_)
    implicit none
    class(BasisObsData),pointer:: obsmold_
    character(len=*),intent(in):: vname

    character(len=*),parameter:: myname_=myname//"::vname2vmold_"

    obsmold_=>NUll()
    select case(trim(vname))
      case("Radiance");   obsmold_%Obspoint => Radiance
      case("Radiosonde"); obsmold_%Obspoint => Radiosonde
      case("Aod");   obsmold_%Obspoint => Aod
    end select
  end function vname2vmold_

end module ufo_obs_data_mod
