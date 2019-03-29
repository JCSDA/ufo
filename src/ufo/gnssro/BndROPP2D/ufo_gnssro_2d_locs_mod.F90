module ufo_gnssro_2d_locs_mod

use iso_c_binding
use fckit_log_module, only : fckit_log
use kinds,            only: kind_real
use ufo_locs_mod
  use gnssro_mod_conf
public:: ufo_gnssro_2d_locs_init

contains

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
subroutine ufo_gnssro_2d_locs_init(self, obss, t1, t2, loc2dconf)
  use kinds
  use datetime_mod
  use twindow_utils_mod
  use fckit_log_module, only : fckit_log
  use obsspace_mod

  implicit none

  type(ufo_locs),      intent(inout) :: self
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
  real(kind_real), dimension(:), allocatable     :: obsAzim
  type(conf2d),                  intent(in)      :: loc2dconf
  real(kind_real), dimension(loc2dconf%n_horiz)  :: plat_2d, plon_2d
  integer         :: kerror
  real(kind_real) :: dtheta

  dtheta=loc2dconf%res/6371.0

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
     if (obsspace_has(obss, "GroupUndefined", "sensor_azimuth_angle")) then
       call obsspace_get_db(obss, " ", "sensor_azimuth_angle", obsAzim)
     else
       write(err_msg,*) myname, ' error: sensor_azimuth_angle not found'
       call abor1_ftn(err_msg)
     endif
  endif

  !Setup ufo 2d locations 
  call ufo_locs_setup(self, tw_nlocs*loc2dconf%n_horiz)
  do i = 1, tw_nlocs
    call ropp_fm_2d_plane(lon(tw_indx(i)),lat(tw_indx(i)),obsAzim(tw_indx(i)),dtheta, loc2dconf%n_horiz,plat_2d,plon_2d,kerror)
    do j = 1, loc2dconf%n_horiz
      self%lon( (j-1)*tw_nlocs+i ) = plon_2d(j)
      self%lat( (j-1)*tw_nlocs+i ) = plat_2d(j)
      self%time((j-1)*tw_nlocs+i ) = date_time(tw_indx(i))
      self%indx((j-1)*tw_nlocs+i ) = tw_indx(i)+(j-1)*tw_nlocs
    enddo
  end do

 do i = 1, nlocs
    call datetime_delete(date_time(i))
  enddo
  deallocate(date_time, lon, lat, tw_indx, obsAzim)

end subroutine ufo_gnssro_2d_locs_init

!----------------------------
end module ufo_gnssro_2d_locs_mod
