module ufo_gnssro_2d_locs_mod

use iso_c_binding
use fckit_log_module, only : fckit_log
use kinds,            only : kind_real
use ufo_locs_mod
public:: ufo_gnssro_2d_locs_init

contains

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
subroutine ufo_gnssro_2d_locs_init(self, locs, obss)
  use kinds
  use datetime_mod
  use fckit_log_module, only : fckit_log
  use obsspace_mod
  use ufo_gnssro_bndropp2d_mod

  implicit none

  type(ufo_locs),              intent(inout) :: locs
  class(ufo_gnssro_BndROPP2D), intent(inout) :: self
  type(c_ptr),  value, intent(in)    :: obss

  character(len=*),parameter    :: myname = "ufo_gnssro_2d_locs_init"
  integer,         parameter    :: max_string = 800
  character(max_string)         :: err_msg

  integer :: i, j, nlocs
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

  call obsspace_get_db(obss, "MetaData", "datetime", date_time)
  call obsspace_get_db(obss, "MetaData", "longitude", lon)
  call obsspace_get_db(obss, "MetaData", "latitude", lat)

  allocate(obsAzim(nlocs))
  if (obsspace_has(obss,"ObsValue","bending_angle")) then
     if (obsspace_has(obss, "MetaData", "sensor_azimuth_angle")) then
       call obsspace_get_db(obss, "MetaData", "sensor_azimuth_angle", obsAzim)
     else
       write(err_msg,*) myname, ' error: sensor_azimuth_angle not found'
       call abor1_ftn(err_msg)
     endif
  endif

  !Setup ufo 2d locations 
  call ufo_locs_setup(locs, nlocs*n_horiz)
  do i = 1, nlocs
    locs%lon( (i-1)*n_horiz+1 : i*n_horiz) =  lon(i)
    locs%lat( (i-1)*n_horiz+1 : i*n_horiz) =  lat(i)
    do j = 1, n_horiz
      locs%indx((i-1)*n_horiz+j) =  (i-1)*n_horiz+j
      locs%time((i-1)*n_horiz+j) =  date_time(i)
    end do
  end do

  ! save ufo_locs to self
  self%obsLat2d = locs%lat
  self%obsLon2d = locs%lon

  do i = 1, nlocs
    call datetime_delete(date_time(i))
  enddo
  deallocate(date_time, lon, lat, obsAzim)

end subroutine ufo_gnssro_2d_locs_init

!----------------------------
end module ufo_gnssro_2d_locs_mod
