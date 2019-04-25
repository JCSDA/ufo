module ufo_gnssro_2d_locs_mod

use iso_c_binding
use fckit_log_module, only : fckit_log
use kinds,            only : kind_real
use ufo_locs_mod
public:: ufo_gnssro_2d_locs_init

contains

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
subroutine ufo_gnssro_2d_locs_init(self, locs, obss, t1, t2)
  use kinds
  use datetime_mod
  use twindow_utils_mod
  use fckit_log_module, only : fckit_log
  use obsspace_mod
  use ufo_gnssro_bndropp2d_mod

  implicit none

  type(ufo_locs),              intent(inout) :: locs
  class(ufo_gnssro_BndROPP2D), intent(inout) :: self
  type(c_ptr),  value, intent(in)    :: obss
  type(datetime),      intent(in)    :: t1, t2

  character(len=*),parameter    :: myname = "ufo_gnssro_2d_locs_init"
  integer,         parameter    :: max_string = 800
  character(max_string)         :: err_msg

  integer :: i, j, tw_nlocs,nlocs
  integer,         dimension(:), allocatable :: tw_indx
  real(kind_real), dimension(:), allocatable :: lon, lat
  type(datetime),  dimension(:), allocatable :: date_time

! gnss ro data 2d location  
  real(kind_real), dimension(:), allocatable      :: obsAzim
  real(kind_real), dimension(self%roconf%n_horiz) :: plat_2d, plon_2d
  integer         :: kerror, n_horiz
  real(kind_real) :: dtheta

  dtheta  = self%roconf%dtheta
  n_horiz = self%roconf%n_horiz
 ! Local copies pre binning
  nlocs = obsspace_get_nlocs(obss)

  allocate(date_time(nlocs), lon(nlocs), lat(nlocs))

  call obsspace_get_db(obss, "", "datetime", date_time)
  call obsspace_get_db(obss, "", "longitude", lon)
  call obsspace_get_db(obss, "", "latitude", lat)

  ! Generate the timing window indices
  allocate(tw_indx(nlocs))
  tw_nlocs = 0
  do i = 1, nlocs
    if (date_time(i) > t1 .and. date_time(i) <= t2) then
      tw_nlocs = tw_nlocs + 1
      tw_indx(tw_nlocs) = i
    endif
  enddo

  allocate(obsAzim(nlocs))
  if (obsspace_has(obss,"ObsValue","bending_angle")) then
     if (obsspace_has(obss, "GroupUndefined", "sensor_azimuth_angle")) then
       call obsspace_get_db(obss, " ", "sensor_azimuth_angle", obsAzim)
     else
       write(err_msg,*) myname, ' error: sensor_azimuth_angle not found'
       call abor1_ftn(err_msg)
     endif
  endif

  !Setup ufo 2d locations 
  call ufo_locs_setup(locs, tw_nlocs*n_horiz)
  do i = 1, tw_nlocs
    call ropp_fm_2d_plane(lat(tw_indx(i)),lon(tw_indx(i)),obsAzim(tw_indx(i)),dtheta,n_horiz,plat_2d,plon_2d,kerror)
    locs%lon( (i-1)*n_horiz+1 : i*n_horiz) =  plon_2d
    locs%lat( (i-1)*n_horiz+1 : i*n_horiz) =  plat_2d
    locs%time((i-1)*n_horiz+1 : i*n_horiz) =  date_time(tw_indx(i))
  ! save ufo_locs to self
    self%obsLat2d( (i-1)*n_horiz+1 : i*n_horiz)  = locs%lat( (i-1)*n_horiz+1 : i*n_horiz) 
    self%obsLon2d( (i-1)*n_horiz+1 : i*n_horiz)  = locs%lon( (i-1)*n_horiz+1 : i*n_horiz) 
    do j = 1, n_horiz
      locs%indx((i-1)*n_horiz+j) =  (tw_indx(i)-1)*n_horiz+j
    end do
  end do

  do i = 1, nlocs
    call datetime_delete(date_time(i))
  enddo
  deallocate(date_time, lon, lat, tw_indx, obsAzim)

end subroutine ufo_gnssro_2d_locs_init

!----------------------------
end module ufo_gnssro_2d_locs_mod
