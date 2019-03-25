module ufo_gnssro_2d_locs_mod

use iso_c_binding
use fckit_log_module, only : fckit_log
use kinds,            only: kind_real
use ufo_locs_mod

public:: ufo_gnssro_2d_locs_init

contains

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
subroutine ufo_gnssro_2d_locs_init(self, obss, t1, t2)
  use kinds
  use datetime_mod
  use twindow_utils_mod
  use fckit_log_module, only : fckit_log
  use obsspace_mod

  implicit none

  type(ufo_locs), intent(inout) :: self
  type(c_ptr), value, intent(in)              :: obss
  type(datetime), intent(in)                  :: t1, t2

  integer :: nlocs

  character(len=*),parameter:: &
     myname = "ufo_gnssro_2d_locs_init"
  integer :: i
  integer :: tw_nlocs
  integer, dimension(:), allocatable :: tw_indx
  real(kind_real), dimension(:), allocatable :: lon, lat
  type(datetime), dimension(:), allocatable :: date_time

! gnss ro data 2d location  
  integer, parameter          :: n_horiz = 31           ! from ropp2d lib
  real*8,  parameter          :: dtheta = 40.0/6371.0   ! from ropp2d lib
  real*8,  dimension(n_horiz) :: plat_2d, plon_2d
  integer                     ::  kerror
  real(kind_real), dimension(:), allocatable :: obsAzim

 ! Local copies pre binning
  nlocs = obsspace_get_nlocs(obss)

  allocate(date_time(nlocs), lon(nlocs), lat(nlocs))

!TODO(JG): Add "MetaData" or similar group attribute to all ioda ObsSpace objects
  if (obsspace_has(obss,"MetaData", "time")) then
    call obsspace_get_db(obss, "MetaData", "datetime", date_time)
  else
    call obsspace_get_db(obss, "", "datetime", date_time)
  endif

  ! Generate the timing window indices
  allocate(tw_indx(nlocs))
  tw_nlocs = 0
  do i = 1, nlocs
    if (date_time(i) > t1 .and. date_time(i) <= t2) then
      tw_nlocs = tw_nlocs + 1
      tw_indx(tw_nlocs) = i
    endif
  enddo

!TODO(JG): Add "MetaData" or similar group attribute to all ioda ObsSpace objects
  if (obsspace_has(obss,"MetaData", "longitude")) then
    call obsspace_get_db(obss, "MetaData", "longitude", lon)
    call obsspace_get_db(obss, "MetaData", "latitude", lat)
  else
    call obsspace_get_db(obss, "", "longitude", lon)
    call obsspace_get_db(obss, "", "latitude", lat)
  endif

  allocate(obsAzim(nlocs))
  if (obsspace_has(obss,"ObsValue","bending_angle")) then
     call obsspace_get_db(obss, " ", "sensor_azimuth_angle", obsAzim)
  endif

  !Setup ufo locations
  call ufo_locs_setup(self, tw_nlocs*n_horiz)
  do i = 1, tw_nlocs
    call ropp_fm_2d_plane(lon(tw_indx(i)),lat(tw_indx(i)),obsAzim(tw_indx(i)),dtheta, n_horiz,plat_2d,plon_2d,kerror)
    self%lon( (i-1)*n_horiz+1:i*n_horiz )  = plon_2d
    self%lat( (i-1)*n_horiz+1:i*n_horiz )  = plat_2d
    self%time((i-1)*n_horiz+1:i*n_horiz )  = date_time(tw_indx(i))
  enddo
  self%indx = tw_indx(1:tw_nlocs)
 do i = 1, nlocs
    call datetime_delete(date_time(i))
  enddo
  deallocate(date_time, lon, lat, tw_indx, obsAzim)

end subroutine ufo_gnssro_2d_locs_init

!----------------------------
end module ufo_gnssro_2d_locs_mod
