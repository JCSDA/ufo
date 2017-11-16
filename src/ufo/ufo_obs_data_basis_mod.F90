! (C) Copyright 2017 UCAR
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 

!> Fortran module to handle radiosondeentional observations

module ufo_obs_data_basis_mod
  use iso_c_binding
  implicit none

  character(len=*), parameter :: MyName='basis_obs_data_mod'
  type, abstract :: BasisObsData
    integer :: nobs =0
    integer :: nlocs=0
    character(len=800) :: filein, fileout
    class(BasisObsData),pointer :: Obspoint => NULL()
    contains
      procedure, nopass :: echoMyname
      procedure(Setup_),  deferred :: Setup
      procedure(Delete_), deferred :: Delete
      generic :: SetupBasis => Setup
      generic :: DeleteBasis => Delete
  end type BasisObsData

  abstract interface 
    ! Interface for setup
    subroutine Setup_(self, filein,obstype,nobs,nlocs)
      import
      class(BasisObsData), intent(inout) :: self
      character(len=*), intent(in)    :: filein
      character(len=*), intent(in)    :: obstype
      integer(c_int),   intent(inout) :: nobs
      integer(c_int),   intent(inout) :: nlocs
    end subroutine Setup_
    ! Interface for setup
    subroutine Delete_(self)
      import
      class(BasisObsData), intent(inout) :: self
    end subroutine Delete_
  end interface

contains
    !  echo module name
    subroutine echoMyname(self)
      class(BasisObsData), intent(inout) :: self
      print *, 'basis module name: ', MyName
    end subroutine echoMyname
end module ufo_obs_data_basis_mod
