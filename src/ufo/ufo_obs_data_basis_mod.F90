module ufo_obs_data_basis_mod
  use iso_c_binding
  implicit none

  character(len=*), parameter :: MyName='basis_obs_data_mod'
  type, abstract :: BasisObsData
    integer :: nobs=0
    character(len=800) :: filein, fileout
    contains
      procedure, nopass :: echoMyname
      procedure(SetupSomethingInterface), deferred :: SetupSomething
      procedure(DeleteInterface), deferred :: Delete
      generic :: SetupBasis => SetupSomething
      generic :: DeleteBasis => Delete
  end type BasisObsData

  interface 
    ! Interface for setup
    subroutine SetupSomethingInterface(self, filein,obstype,nobs)
      import
      class(BasisObsData), intent(inout) :: self
      character(len=*), intent(in)    :: filein
      character(len=*), intent(in)    :: obstype
      integer(c_int),   intent(inout) :: nobs
    end subroutine SetupSomethingInterface
    ! Interface for setup
    subroutine DeleteInterface(self)
      import
      class(BasisObsData), intent(inout) :: self
    end subroutine DeleteInterface
  end interface

contains
    !  echo module name
    subroutine echoMyname(self)
      class(BasisObsData), intent(inout) :: self
      print *, 'basis module name: ', MyName
    end subroutine echoMyname
end module ufo_obs_data_basis_mod
