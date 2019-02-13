! (C) Copyright 2018 UCAR
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 

!> Fortran module to perform linear interpolation

module twindow_utils_mod

use kinds, only: kind_real
use datetime_mod
use duration_mod

implicit none
private

public gen_twindow_index

contains

! ------------------------------------------------------------------------------
subroutine gen_twindow_index(refdate, t1, t2, nlocs, time_offset, tw_indx, tw_nlocs)
  implicit none

  type(datetime)                            :: refdate
  type(datetime)                            :: t1
  type(datetime)                            :: t2
  integer, intent(in)                       :: nlocs
  real(kind_real), dimension(:), intent(in) :: time_offset
  integer, dimension(:), intent(out)        :: tw_indx
  integer, intent(out)                      :: tw_nlocs

  integer :: i

  type(duration), dimension(:), allocatable :: dt
  type(datetime), dimension(:), allocatable :: t

  ! Convert the refdate, time offset pairs to an absolute time that
  ! can be compared to the timing window edges.
  allocate(dt(nlocs))
  allocate(t(nlocs))

  do i = 1, nlocs
    dt(i) = int(3600*time_offset(i))
    t(i) = refdate
    call datetime_update(t(i), dt(i))
  enddo

  ! Find number of locations in this timeframe
  tw_nlocs = 0
  do i = 1, nlocs
    if (t(i) > t1 .and. t(i) <= t2) then
      tw_nlocs = tw_nlocs + 1
      tw_indx(tw_nlocs) = i
    endif
  enddo

  deallocate(dt)
  deallocate(t)

end subroutine gen_twindow_index

! ------------------------------------------------------------------------------

end module twindow_utils_mod
