module ufo_timeoper_locs_mod

use iso_c_binding
use fckit_log_module, only : fckit_log
use kinds,            only : kind_real
use ufo_locs_mod
use ufo_timeoper_mod


public:: ufo_timeoper_locs_init

contains

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
subroutine ufo_timeoper_locs_init(self, locs, obss, t0, t1, t2, t3, st)

use kinds
use datetime_mod
use duration_mod
use twindow_utils_mod
use obsspace_mod

implicit none

type(ufo_locs),      intent(inout) :: locs
class(ufo_timeoper), intent(inout) :: self
type(c_ptr),  value, intent(in)    :: obss
type(datetime),      intent(in)    :: t0, t1, t2, t3
type(datetime),      intent(in)    :: st ! stateTime

! local variables
character(len=*),parameter    :: myname = "ufo_timeoper_locs_init"

integer :: i, j, k, tw_nlocs, nlocs
integer,         dimension(:), allocatable :: tw_indx
real(kind_real), dimension(:), allocatable :: lon, lat, lon_in, lat_in
type(datetime),  dimension(:), allocatable :: date_time, date_time_in
type(duration)               :: dt,dt2
character(len=12), parameter :: timeop_interp_str = 'TimeOpInterp'

logical                      :: l_full_window

! Local copies pre binning
nlocs = obsspace_get_nlocs(obss)

allocate(date_time_in(nlocs), lon_in(nlocs), lat_in(nlocs))
allocate(date_time(self%time_stencil*nlocs), lon(self%time_stencil*nlocs), &
         lat(self%time_stencil*nlocs))

call obsspace_get_db(obss, "MetaData", "datetime", date_time_in)
call obsspace_get_db(obss, "MetaData", "longitude", lon_in)
call obsspace_get_db(obss, "MetaData", "latitude", lat_in)

k = 0
do i = 1, nlocs
  do j = 1, self%time_stencil
    k = k + 1
    date_time(k) = date_time_in(i)
    lon(k) = lon_in(i)
    lat(k) = lat_in(i)
  end do
end do

! Generate the timing window indices
call datetime_diff(t0, t1, dt)
call datetime_diff(t2, t3, dt2)
l_full_window = ((duration_seconds(dt) == 0) .and. (duration_seconds(dt2) == 0))

if (l_full_window) then
  tw_nlocs = nlocs * self%time_stencil
  allocate(tw_indx(tw_nlocs))
  do i = 1, tw_nlocs
    tw_indx(i) = i
  end do
else
  allocate(tw_indx(nlocs))
  tw_nlocs = 0
  do i = 1, nlocs
    if (date_time_in(i) >= t0 .and. date_time_in(i) < t3 &
       .and. date_time_in(i) >= st) then
      tw_nlocs = tw_nlocs + 1
      tw_indx(tw_nlocs) = 2 * i - 1
    elseif  (date_time_in(i) >= t0 .and. date_time_in(i) < t3 &
      .and. date_time_in(i) < st) then
      tw_nlocs = tw_nlocs + 1
      tw_indx(tw_nlocs) = 2 * i
    endif
  enddo
end if

call ufo_locs_setup(locs, tw_nlocs)

do i = 1, tw_nlocs
  locs%lon(i) = lon(tw_indx(i))
  locs%lat(i) = lat(tw_indx(i))
  locs%time(i) = date_time(tw_indx(i))
  locs%indx(i) = tw_indx(i)
end do

do i = 1, nlocs
  call datetime_delete(date_time_in(i))
enddo

do i = 1, nlocs * self%time_stencil
  call datetime_delete(date_time(i))
enddo

deallocate(date_time, date_time_in, lon, lon_in, lat, lat_in, tw_indx)

end subroutine ufo_timeoper_locs_init

!------------------------------------------------------------------------------
end module ufo_timeoper_locs_mod
