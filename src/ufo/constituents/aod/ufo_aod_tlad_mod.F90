! (C) Copyright 2017 UCAR
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 

!> Fortran module to handle aod observations

MODULE ufo_aod_tlad_mod
implicit none

public :: ufo_aod_tlad
  
!> Fortran derived type for aod trajectory
type :: ufo_aod_tlad
   logical :: ltraj = .false. !< trajectory set?
end type ufo_aod_tlad

! ------------------------------------------------------------------------------

contains
  

END MODULE ufo_aod_tlad_mod
