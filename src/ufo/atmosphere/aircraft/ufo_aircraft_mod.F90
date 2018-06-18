! (C) Copyright 2017 UCAR
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 

!> Fortran module to handle aircraft observations

module ufo_aircraft_mod
  
  use ufo_conventional_profile_mod, only: ufo_conventional_profile

  implicit none
  public :: ufo_aircraft
  private

  !> Fortran derived type for aircraft_t trajectory
  type, extends(ufo_conventional_profile) :: ufo_aircraft
  contains
  end type ufo_aircraft

contains
    
! ------------------------------------------------------------------------------


! ------------------------------------------------------------------------------

end module ufo_aircraft_mod
