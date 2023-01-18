module ufo_gnssro_2d_locs_mod

use iso_c_binding
use fckit_log_module, only : fckit_log
use kinds,            only : kind_real

private
public:: ufo_gnssro_2d_locs_init

contains

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
subroutine ufo_gnssro_2d_locs_init(self, obss, nlocs_ext, lons, lats)
  use kinds
  use datetime_mod
  use obsspace_mod
  use ufo_gnssro_bndropp2d_mod

  implicit none

  class(ufo_gnssro_BndROPP2D), intent(inout) :: self
  type(c_ptr),  value, intent(in)  :: obss
  integer, intent(in) :: nlocs_ext
  real(c_float), dimension(nlocs_ext), intent(inout)  :: lons, lats

  character(len=*),parameter    :: myname = "ufo_gnssro_2d_locs_init"
  integer,         parameter    :: max_string = 800
  character(max_string)         :: err_msg

  integer :: i, j,nlocs
  real(kind_real), dimension(:), allocatable :: lon, lat

! gnss ro data 2d location  
  real(kind_real), dimension(:), allocatable      :: obsAzim
  real(kind_real), dimension(self%roconf%n_horiz) :: plat_2d, plon_2d
  integer         :: kerror, n_horiz
  real(kind_real) :: dtheta

  dtheta  = self%roconf%dtheta
  n_horiz = self%roconf%n_horiz
  nlocs = nlocs_ext / n_horiz

  allocate(lon(nlocs), lat(nlocs))
  call obsspace_get_db(obss, "MetaData", "longitude", lon)
  call obsspace_get_db(obss, "MetaData", "latitude", lat)

  allocate(obsAzim(nlocs))
  call obsspace_get_db(obss, "MetaData", "sensorAzimuthAngle", obsAzim)

  !Setup ufo 2d locations 
  do i = 1, nlocs
    call ropp_fm_2d_plane(lat(i),lon(i),obsAzim(i),dtheta,n_horiz,plat_2d,plon_2d,kerror)
    lons( (i-1)*n_horiz+1 : i*n_horiz) =  plon_2d
    lats( (i-1)*n_horiz+1 : i*n_horiz) =  plat_2d
  ! save ufo_locs to self
    self%obsLat2d( (i-1)*n_horiz+1 : i*n_horiz)  = lats( (i-1)*n_horiz+1 : i*n_horiz) 
    self%obsLon2d( (i-1)*n_horiz+1 : i*n_horiz)  = lons( (i-1)*n_horiz+1 : i*n_horiz) 
  end do

  deallocate(lon, lat, obsAzim)

end subroutine ufo_gnssro_2d_locs_init

!----------------------------
end module ufo_gnssro_2d_locs_mod
