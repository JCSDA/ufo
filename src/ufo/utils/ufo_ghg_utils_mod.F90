! (C) Copyright 2018-2025 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

module ufo_ghg_utils_mod

use kinds, only: kind_real

implicit none
private

public co2_mauna_loa_rise

contains

!----------------------------------------------------------------------------
! Carbon Dioxide Options ----------------------------------------------------
!----------------------------------------------------------------------------

subroutine co2_mauna_loa_rise(julday,co2)

! note: the CO2 value that comes out of this routine might be
!       sligthly different than the value returned by GSI because
!       in some cases, the time defined in "today" is the start
!       of the DA cycle, which could be offset from the center
!       synoptic time used in GSI. That should not be a big deal.

  implicit none
  real(kind=kind_real), intent(in)    :: julday
  real(kind=kind_real), intent(inout) :: co2     !Carbon dioxide (ppmv)

  co2 = 0.00602410_kind_real * (julday - 2455563.0_kind_real) + 389.5_kind_real

end subroutine co2_mauna_loa_rise

!----------------------------------------------------------------------------

end module ufo_ghg_utils_mod
