! (C) Copyright 2023 NOAA, UCAR, and Met Office
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

module ufo_interp_param_mod

  implicit none

  public

  integer, parameter :: UNSPECIFIED_INTERP = -1
  integer, parameter :: LINEAR_INTERP = 1
  integer, parameter :: LOG_LINEAR_INTERP = 2
  integer, parameter :: NEAREST_NEIGHBOR_INTERP = 3

end module ufo_interp_param_mod

