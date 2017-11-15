module ufo_obs_data_mod
  use iso_c_binding

  use ufo_obs_data_basis_mod

  use radDiag_mod, only: RadDiag
  use radDiag_mod, only: RadDiag_read
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
    procedure :: SetupSomething
    procedure :: Delete
  end type Obs_Data

  type Radiance_Obs
     character(len=80) :: myname
     integer :: ndim=0
     real, pointer :: lat(:) 
     real, pointer :: lon(:)
     real, pointer :: lev(:)
     real, pointer :: obs(:,:) 
  end type Radiance_Obs
  type Radiosonde_Obs
     character(len=80) :: myname
     integer :: ndim=0
     real, pointer :: lat(:)
     real, pointer :: lon(:)
     real, pointer :: lev(:)
     real, pointer :: obs(:) 
  end type Radiosonde_Obs
contains

  subroutine SetupSomething(self, filein,obstype,nobs)
    class(Obs_Data), intent(inout) :: self
    character(len=*), intent(in)    :: filein
    character(len=*), intent(in)    :: obstype
    integer(c_int),   intent(inout) :: nobs

    type(Radiance_Obs),  pointer :: Radiance
    type(Radiosonde_Obs),pointer :: Radiosonde
    type(RadDiag),       pointer :: Rad
    type(RaobDiag),      pointer :: Raob

    if (trim(obstype) == "Radiance"   ) then
        call SetupRadiance(self, Radiance, filein,nobs)
    endif
    if (trim(obstype) == "Radiosonde" ) then
        call SetupRaob(self, Radiosonde, filein,nobs)
    endif

  end subroutine SetupSomething

  subroutine Delete(self)
    class(Obs_Data), intent(inout) :: self
  end subroutine Delete

  subroutine SetupRadiance(self, mytype, filein,nobs)
    class(Obs_Data), intent(inout) :: self
    type(Radiance_Obs), intent(inout)  :: mytype
    character(len=*), intent(in)    :: filein
    integer(c_int),   intent(inout) :: nobs

    type(RadDiag), pointer :: Rad

    print *, "Obs_Data:ReadRadiance"
    allocate(Rad)
    call radDiag_read(Rad,filein,nobs)
    self%Nobs = nobs

    ! NOTE: at this point this routine (or another one) can fill 
    !       in the entries in the general vector type

  end subroutine SetupRadiance

  subroutine SetupRaob(self, mytype, filein,nobs)
    class(Obs_Data), intent(inout) :: self
    type(Radiosonde_Obs), intent(inout) :: mytype
    character(len=*), intent(in)    :: filein
    integer(c_int),   intent(inout) :: nobs

    type(RaobDiag),pointer :: Raob

    print *, "Obs_Data:ReadRaob"
    allocate(Raob)
    call raobDiag_read(Raob,filein,nobs)
    self%Nobs = nobs

    ! NOTE: at this point this routine (or another one) can fill 
    !       in the entries in the general vector type

  end subroutine SetupRaob

end module ufo_obs_data_mod
