! (C) Copyright 2017-2018 UCAR
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 

!> Fortran module to handle radiosonde observations

module ufo_radiosonde_mod
  
  use ufo_conventional_profile_mod, only: ufo_conventional_profile

  implicit none
  public :: ufo_radiosonde
  private

  !> Fortran derived type for radiosonde_t trajectory
  type, extends(ufo_conventional_profile) :: ufo_radiosonde
  contains
  end type ufo_radiosonde

contains
    
! ------------------------------------------------------------------------------
! ------------------------------------------------------------------------------

end module ufo_radiosonde_mod
