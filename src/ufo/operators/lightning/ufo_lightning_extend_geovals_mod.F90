module ufo_lightning_extend_geovals_mod

  use iso_c_binding, only: c_ptr, c_float
  use fckit_log_module, only: fckit_log
  use kinds, only: kind_real

  implicit none

  private
  public :: ufo_lightning_extend_geovals
contains

!-------------------------------------------------------------------------
! This subroutine extends geovals data sets by calculating
! 2D locations for a predefined area around each observation point.
! We referenced code from another UFO directory: BndROPP2D/ufo_gnssro_2d_locs_mod.F90
subroutine ufo_lightning_extend_geovals(self, obss, nlocs_ext, lons, lats)
  use obsspace_mod, only: obsspace_get_db
  use ufo_lightning_mod, only: ufo_lightning

  implicit none

  class(ufo_lightning), intent(inout) :: self
  type(c_ptr), value, intent(in) :: obss
  integer, intent(in) :: nlocs_ext
  real(c_float), dimension(nlocs_ext), intent(inout) :: lons, lats

  character(len=*), parameter :: myname = "ufo_lightning_extend_geovals"
  integer :: i, nlocs, n_horiz
  real(kind_real), dimension(:), allocatable :: lon, lat, plat_2d, plon_2d
  integer :: kerror
  integer, parameter    :: max_string = 800
  character(max_string) :: err_msg         ! Error message for output

  kerror = 1
  n_horiz = self%n_horiz
  nlocs = nlocs_ext / n_horiz
  allocate(lon(nlocs), lat(nlocs))
  allocate(plat_2d(n_horiz), plon_2d(n_horiz))
  call obsspace_get_db(obss, "MetaData", "longitude", lon)
  call obsspace_get_db(obss, "MetaData", "latitude", lat)

  ! Populate 2D locations for the area around each observation
  do i = 1, nlocs
    call find_partner_points(lat(i), lon(i), n_horiz, plat_2d, plon_2d, kerror)
    if (kerror > 0) then
        write(err_msg,*) myname, ' error: find_partner_points does not run successfully!' 
    endif
    lons((i-1)*n_horiz+1 : i*n_horiz) = plon_2d
    lats((i-1)*n_horiz+1 : i*n_horiz) = plat_2d
    ! Save extended locations to the ufo_lightning object
    !self%obsLat2d((i-1)*n_horiz+1 : i*n_horiz) = plat_2d
    !self%obsLon2d((i-1)*n_horiz+1 : i*n_horiz) = plon_2d
  end do

  deallocate(lon, lat, plat_2d, plon_2d)

end subroutine ufo_lightning_extend_geovals
!-------------------------------------------------------------------------

! This subroutine calculates partner points around a central observation point,
! simulating a 15x15 km^2 grid area. These points are used for ordering extended 
! geovals and facilitating horizontal integration within the observation operators.
subroutine find_partner_points(center_lat, center_lon, n_horiz, plat_2d, plon_2d, kerror)
  use ufo_constants_mod, only: mean_earth_rad, pi
  implicit none

  real(kind_real), intent(in) :: center_lat, center_lon
  integer, intent(in) :: n_horiz
  real(kind_real), dimension(n_horiz), intent(out) :: plat_2d, plon_2d
  integer, intent(inout) :: kerror

  integer :: i, j, grid_size, indx
  real(kind_real) :: lat, lon, delta_lat, delta_lon, distance

  distance = 15000.0 / sqrt(real(n_horiz, kind_real))  ! Calculate spacing based on number of partner points
  grid_size = nint(sqrt(real(n_horiz, kind_real)))  ! Determine the grid size from n_horiz
  ! Calculate latitude increment (delta_lat) for grid points
  ! This approximation works well for small areas where Earth's curvature is minimal.
  ! It converts a distance on Earth's surface to an angular change in latitude.
  delta_lat = (distance /(1000 * mean_earth_rad)) * (180.0 / pi)

  ! Calculate longitude increment (delta_lon) for grid points
  ! Adjusts for Earth's shape by accounting for the cosine of the latitude,
  ! ensuring the longitude increment matches the specified distance at this latitude.
  ! This approximation assumes a spherical Earth and is accurate for small distances (< 15 km in our case).
  delta_lon = delta_lat / cos(center_lat * pi / 180.0)

  ! Generate and print grid points around the center
  if (mod(grid_size, 2) == 0) then
    ! For even grid_size, adjust bounds to ensure symmetric distribution around the center
    do i = -(grid_size / 2), grid_size / 2 - 1
       lat = center_lat + i * delta_lat
       do j = -(grid_size / 2), grid_size / 2 - 1
          lon = center_lon + j * delta_lon
          ! Adjust longitude to ensure it's within the 0 to 360 range
          !if (lon > 180.0) lon = lon - 360.0
          if (lon < -180.0) lon = lon + 360.0
          indx = (i + grid_size / 2) * grid_size + (j + grid_size / 2) + 1
          plon_2d(indx) = lon
          plat_2d(indx) = lat
       end do
    end do
  else
    ! For odd grid_size, loop bounds are centered around 0
    do i = -(grid_size - 1) / 2, (grid_size - 1) / 2
       lat = center_lat + i * delta_lat
       do j = -(grid_size - 1) / 2, (grid_size - 1) / 2
          lon = center_lon + j * delta_lon
          ! Ensure longitude is within the 0 to 360 range
          !if (lon > 180.0) lon = lon - 360.0
          if (lon < -180.0) lon = lon + 360.0
          indx = (i + (grid_size - 1) / 2) * grid_size + (j + (grid_size - 1) / 2) + 1
          plon_2d(indx) = lon
          plat_2d(indx) = lat
       end do
    end do
  end if
  kerror = 0  ! Assume successful operation
end subroutine find_partner_points

end module ufo_lightning_extend_geovals_mod
