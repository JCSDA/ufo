! (C) Copyright 2017-2018 UCAR
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 

!> Fortran module to handle radiosonde observations

module ufo_radiosonde_tlad_mod
  
  use ufo_conventional_profile_tlad_mod, only: ufo_conventional_profile_tlad

  implicit none
  public :: ufo_radiosonde_tlad
  private

!> Fortran derived type for radiosonde_t trajectory
type, extends(ufo_conventional_profile_tlad) :: ufo_radiosonde_tlad
contains
end type ufo_radiosonde_tlad

contains

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------

end module ufo_radiosonde_tlad_mod
