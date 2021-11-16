! (C) Copyright 2021
!
! This software is developed by NOAA/NWS/EMC under the Apache 2.0 license
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

module aero_kinds_mod
  use, intrinsic :: iso_c_binding
  implicit none

  private
  public kind_real

  integer, parameter :: kind_real=c_double
end module aero_kinds_mod
