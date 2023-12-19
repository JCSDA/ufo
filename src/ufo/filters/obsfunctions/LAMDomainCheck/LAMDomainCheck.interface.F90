! (C) Copyright 2020 NOAA NWS NCEP EMC
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

!> Fortran module to handle verifying that an observation is within a limited area model domain

module ufo_lamdomaincheck_mod_c

  use iso_c_binding
  use kinds
  use ufo_constants_mod, only: pi, deg2rad, rad2deg, two, half, mean_earth_rad

  implicit none

  private

contains

! -----------------------------------------------------------------------------
!> \brief subroutine lam_domaincheck_esg_c
!!
!! \details **lam_domaincheck_esg_c()** is a subroutine that for a given input defined ESG regional grid
!! and a given input lat/lon point, determines if the point is within or outside the regional domain
!! It takes the following arguments as input:
!! * float c_a - alpha parameter for ESG grid definition
!! * float c_k - kappa paremeter for ESG grid definition
!! * float c_plat - center point latitude of ESG grid (degrees)
!! * float c_plon - center point longitude of ESG grid (degrees)
!! * float c_pazi - azimuth angle for ESG grid definition (radians)
!! * int c_npx - number of grid points in x direction
!! * int c_npy - number of grid points in y direction
!! * float c_dx - grid spacing in degrees
!! * float c_dy - grid spacing in degrees
!! * float c_lat - input latitude (degrees)
!! * float c_lon - input longitude (degrees)
!!
!! and returns c_mask, an integer of 1 (inside the domain) or 0 (outside the domain)
!! The above input arguments, c_lat and c_lon are independent for each observation.
!! The other input arguments are available as global attributes in the FV3 regional grid netCDF file.
!!

subroutine lam_domaincheck_esg_c(c_a, c_k, c_plat, c_plon, c_pazi, c_npx, c_npy,&
                                 c_dx, c_dy, c_nbdy, c_lat, c_lon, c_mask) &
                                 bind(c, name='lam_domaincheck_esg_f90')
  use esg_grid_mod, only: gtoxm_ak_rr
  implicit none
  real(c_float),  intent(in   ) :: c_a, c_k, c_plat, c_plon, c_pazi, c_dx, c_dy
  real(c_float),  intent(in   ) :: c_lat, c_lon
  integer(c_int), intent(in   ) :: c_npx, c_npy, c_nbdy
  integer(c_int), intent(inout) :: c_mask
  real(kind_real), dimension(2) :: xm
  logical :: failure

  real(kind_real) :: a, k, plat, plon, pazi, dx, dy, lat, lon
  integer :: npx, npy, nbdy

  !! convert integers
  npx = int(c_npx)
  npy = int(c_npy)
  nbdy = int(c_nbdy)
  !! convert from C to kind_real
  a = real(c_a, kind_real)
  k = real(c_k, kind_real)
  ! some need converted from degrees to radians
  plat = real(c_plat, kind_real)*deg2rad
  plon = real(c_plon, kind_real)*deg2rad
  pazi = real(c_pazi, kind_real)
  ! dx and dy are on the supergrid, for actual grid resolution is half
  dx = real(c_dx, kind_real)*deg2rad*two
  dy = real(c_dy, kind_real)*deg2rad*two
  lat = real(c_lat, kind_real)*deg2rad
  lon = real(c_lon, kind_real)*deg2rad

  ! call highest level subroutine from Jim Purser's pesg.f90 code
  call gtoxm_ak_rr(a, k, plat, plon, pazi,&
                   dx, dy, lat, lon, xm, failure)

  ! use xm to determine if mask is 1 (good) or 0 (bad)
  c_mask = 0
!  if ((abs(xm(1)) < npx/2) .and. (abs(xm(2)) < npy/2) .and. (.not. failure)) then
! Add a parameter "nbdy" to control the depth of model grid-points from the lateral boundary where observations are rejected.
   if ((abs(xm(1)) < npx/2-nbdy) .and. (abs(xm(2)) < npy/2-nbdy) .and. (.not. failure)) then
    c_mask = 1
  end if

end subroutine lam_domaincheck_esg_c

! -----------------------------------------------------------------------------
!> \brief subroutine lam_domaincheck_circle_c
!!
!! \details **lam_domaincheck_circle_c()** is a subroutine that for a given input defined circle regional grid
!! and a given input central lat/lon point and radius of circle, determines if the point is within or outside
!! the regional domain. Circle domain could be typical for regional MPAS applications.
!! It takes the following arguments as input:
!! * float c_cenlat - center point latitude of circle domain (degrees)
!! * float c_cenlon - center point longitude of circle domain (degrees)
!! * float c_radius - radius of circle domain (km)
!! * float c_lat - input latitude (degrees)
!! * float c_lon - input longitude (degrees)
!!
!! and returns c_mask, an integer of 1 (inside the domain) or 0 (outside the domain)
!! The above input arguments, c_lat and c_lon are independent for each observation.
!!

subroutine lam_domaincheck_circle_c(c_cenlat, c_cenlon, c_radius, &
                                    c_lat, c_lon, c_mask) &
                                 bind(c, name='lam_domaincheck_circle_f90')
  implicit none
  real(c_float),  intent(in   ) :: c_cenlat, c_cenlon, c_radius
  real(c_float),  intent(in   ) :: c_lat, c_lon
  integer(c_int), intent(inout) :: c_mask

  real(kind_real) :: dlat, dlon, rr
  real(kind_real) :: radius, cenlat, cenlon, lat, lon

  ! some need converted from degrees to radians
  cenlat = real(c_cenlat, kind_real)*deg2rad
  cenlon = real(c_cenlon, kind_real)*deg2rad
  radius = real(c_radius, kind_real)    ! in km
  lat = real(c_lat, kind_real)*deg2rad
  lon = real(c_lon, kind_real)*deg2rad

  ! calculate great-circle distance using haversine formula
  dlat = half*abs(lat - cenlat)
  dlon = half*abs(lon - cenlon)
  rr = sqrt( sin(dlat)**2 + cos(lat)*cos(cenlat)*sin(dlon)**2 )
  rr = two*asin(rr)*mean_earth_rad

  ! use rr to determine if mask is 1 (good) or 0 (bad)
  c_mask = 0  ! outside domain
  if (rr < radius) then
    c_mask = 1
  end if

end subroutine lam_domaincheck_circle_c

end module ufo_lamdomaincheck_mod_c
