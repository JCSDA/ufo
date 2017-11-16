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

  implicit none

  type, extends(BasisObsData) :: Obs_Data
  contains
    ! Extending the mapping SetupBasis => SetupRaob
    ! the parent type.
    generic :: SetupBasis => SetupRaob
    generic :: SetupBasis => SetupRadiance
    procedure :: SetupRaob
    procedure :: SetupRadiance
    ! Implementation for the deferred procedure in Basis
    procedure :: Setup
    procedure :: Delete
  end type Obs_Data

  type(RadDiag),  pointer, save::   Rad
  type(RaobDiag), pointer, save::  Raob
contains

  subroutine Setup(self, filein,obstype,nobs)
    class(Obs_Data), intent(inout) :: self
    character(len=*), intent(in)    :: filein
    character(len=*), intent(in)    :: obstype
    integer(c_int),   intent(inout) :: nobs

    select case(trim(obstype))
      case("Radiance")
        allocate(Rad)
        self%Obspoint => Rad
        call SetupRadiance(self, Rad, filein,nobs)
      case("Radiosonde")
        allocate(Raob)
        self%Obspoint => Raob
        call SetupRaob(self, Raob, filein,nobs)
    end select

  end subroutine Setup

  subroutine Delete(self)
    class(Obs_Data), intent(inout) :: self
  end subroutine Delete

  subroutine SetupRadiance(self, mytype, filein,nobs)
    class(Obs_Data), intent(inout) :: self
    type(RadDiag), pointer, intent(inout)    :: mytype
    character(len=*), intent(in)    :: filein
    integer(c_int),   intent(inout) :: nobs

    character(len=*),parameter:: myname_=myname//"::SetupRadiance"

    print *, trim(myname_)
    call radDiag_read(mytype,filein,'Radiance',nobs)
    self%Nobs = nobs

  end subroutine SetupRadiance

  subroutine SetupRaob(self, mytype, filein,nobs)
    class(Obs_Data), intent(inout) :: self
    type(RaobDiag), pointer, intent(inout)   :: mytype
    character(len=*), intent(in)    :: filein
    integer(c_int),   intent(inout) :: nobs

    character(len=*),parameter:: myname_=myname//"::SetupRaob"

    print *, trim(myname_)
    call raobDiag_read(mytype,filein,'Radiosonde',nobs)
    self%Nobs = nobs

  end subroutine SetupRaob

  function vname2vmold_(vname) result(obsmold_)
    implicit none
    class(BasisObsData),pointer:: obsmold_
    character(len=*),intent(in):: vname

    character(len=*),parameter:: myname_=myname//"::vname2vmold_"

    obsmold_=>NUll()
    select case(trim(vname))
      case("Radiance");   obsmold_%Obspoint => Rad
      case("Radiosonde"); obsmold_%Obspoint => Raob
    end select
  end function vname2vmold_

end module ufo_obs_data_mod
