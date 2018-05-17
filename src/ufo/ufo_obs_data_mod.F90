! (C) Copyright 2017 UCAR
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 

!> Fortran module to handle radiosondeentional observations

module ufo_obs_data_mod
  use iso_c_binding

  use ufo_obs_data_basis_mod

  implicit none

  type, extends(BasisObsData) :: Obs_Data
  contains
    ! Implementation for the deferred procedure in Basis
    procedure :: Setup
    procedure :: Delete
    procedure :: GetLocs
  end type Obs_Data

contains

  subroutine Setup(self, filein,obstype,nobs,nlocs)
    class(Obs_Data), intent(inout) :: self
    character(len=*), intent(in)    :: filein
    character(len=*), intent(in)    :: obstype
    integer(c_int),   intent(inout) :: nobs
    integer(c_int),   intent(inout) :: nlocs

  end subroutine Setup

  subroutine Delete(self)
    class(Obs_Data), intent(inout) :: self
  end subroutine Delete

  subroutine GetLocs(self, nlocs, locs)
    class(Obs_Data), intent(in) :: self
    integer, intent(in) :: nlocs
    type(ioda_locs), intent(inout) :: locs
  end subroutine GetLocs

  function vname2vmold_(vname) result(obsmold_)
    implicit none
    class(BasisObsData),pointer:: obsmold_
    character(len=*),intent(in):: vname

    character(len=*),parameter:: myname_=myname//"::vname2vmold_"

    obsmold_=>NUll()
  end function vname2vmold_

end module ufo_obs_data_mod
