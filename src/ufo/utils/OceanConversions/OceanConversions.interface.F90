! (C) Copyright 2022 UK Met Office
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

!> Fortran module to handle calculation of ocean variables

module ufo_oceanconversions_mod_c

  use iso_c_binding
  use kinds

  implicit none

  private

contains

! -----------------------------------------------------------------------------
!> \brief real(c_float) function gsw_rho_t_exact_c
!!
!! \details  real(c_float) function **gsw_rho_t_exact_c**(c_sal, c_temp, c_pressure)
!! calculates density as a function of salinity, temperature and pressure.
!! It takes the following arguments as input:
!! * float c_sal - salinity (g/kg)
!! * float c_temp - temperature (deg.C)
!! * float c_pressure - pressure (dbar)
!!
!! and returns the density (kg/m^3) according to TEOS-10 function gsw_rho_t_exact.
!!

real(c_float) function gsw_rho_t_exact_c(c_sal, c_temp, c_pressure) &
                                 bind(c, name='gsw_rho_t_exact_f90')
  use gsw_mod_toolbox, only: gsw_rho_t_exact
  implicit none
  real(c_float),  intent(in   ) :: c_sal, c_temp, c_pressure
  real(kind_real) :: sal, temp, pressure

  !! convert reals
  sal = real(c_sal, kind_real)
  temp = real(c_temp, kind_real)
  pressure = real(c_pressure, kind_real)

  ! call GSW fortran function
  gsw_rho_t_exact_c = gsw_rho_t_exact(sal, temp, pressure)
  return

end function gsw_rho_t_exact_c


! -----------------------------------------------------------------------------
!> \brief real(c_float) function gsw_p_from_z_c
!!
!! \details real(c_float) function **gsw_p_from_z_c**(c_depth, c_lat)
!! calculates pressure (dbar) using TEOS-10 function gsw_p_from_z.
!! It takes the following arguments as input:
!! * float c_depth - depth (m) - negative down from 0 at surface
!! * float c_lat - input latitude (degrees)
!!
!! and returns pressure (dbar).
!!

real(c_float) function gsw_p_from_z_c(c_depth, c_lat) &
                                 bind(c, name='gsw_p_from_z_f90')
  use gsw_mod_toolbox, only: gsw_p_from_z
  implicit none
  real(c_float),  intent(in   ) :: c_depth, c_lat
  real(kind_real) :: depth, lat, pressure, geo_strf_dyn_height, sea_surface_geopotential

  !! convert reals
  depth = real(c_depth, kind_real)
  lat = real(c_lat, kind_real)
  geo_strf_dyn_height = 0.0
  sea_surface_geopotential = 0.0

  ! call GSW fortran function
  gsw_p_from_z_c = gsw_p_from_z(depth, lat, geo_strf_dyn_height, sea_surface_geopotential)
  return

end function gsw_p_from_z_c


! -----------------------------------------------------------------------------
!> \brief real(c_float) function gsw_pt_from_t_c
!!
!! \details  real(c_float) function **gsw_pt_from_t_c**(c_sal, c_temp, c_pressure)
!! calculates potential temperature as a function of salinity, temperature and pressure.
!! It takes the following arguments as input:
!! * float c_sal - salinity (g/kg)
!! * float c_temp - temperature (deg.C)
!! * float c_pressure - pressure (dbar)
!!
!! and returns potential temperature (deg.C) according to TEOS-10 function gsw_pt_from_t.
!!

real(c_float) function gsw_pt_from_t_c(c_sal, c_temp, c_pressure) &
                                 bind(c, name='gsw_pt_from_t_f90')
  use gsw_mod_toolbox, only: gsw_pt_from_t
  implicit none
  real(c_float),  intent(in   ) :: c_sal, c_temp, c_pressure
  real(kind_real) :: sal, temp, pressure, ref_pressure

  !! convert reals
  sal = real(c_sal, kind_real)
  temp = real(c_temp, kind_real)
  pressure = real(c_pressure, kind_real)
  ref_pressure = 0.0

  ! call GSW fortran function
  gsw_pt_from_t_c = gsw_pt_from_t(sal, temp, pressure, ref_pressure)
  return

end function gsw_pt_from_t_c

! -----------------------------------------------------------------------------
!> \brief real(c_float) function gsw_ct_from_t_c
!!
!! \details  real(c_float) function **gsw_ct_from_t_c**(c_sal, c_temp, c_pressure)
!! calculates conservative temperature as a function of salinity, temperature and pressure.
!! It takes the following arguments as input:
!! * float c_sal - absolute salinity (g/kg)
!! * float c_temp - in situ temperature (deg.C)
!! * float c_pressure - pressure (dbar)
!!
!! and returns conservative temperature (deg.C) according to TEOS-10 function gsw_CT_from_t.
!!

real(c_float) function gsw_ct_from_t_c(c_sal, c_temp, c_pressure) &
                                 bind(c, name='gsw_ct_from_t_f90')
  use gsw_mod_toolbox, only: gsw_CT_from_t
  implicit none
  real(c_float),  intent(in   ) :: c_sal, c_temp, c_pressure
  real(kind_real) :: sal, temp, pressure

  !! convert reals
  sal = real(c_sal, kind_real)
  temp = real(c_temp, kind_real)
  pressure = real(c_pressure, kind_real)

  ! call GSW fortran function
  gsw_ct_from_t_c = gsw_CT_from_t(sal, temp, pressure)
  return

end function gsw_ct_from_t_c

! -----------------------------------------------------------------------------
!> \brief real(c_float) function gsw_sa_from_sp_c
!!
!! \details  real(c_float) function **gsw_sa_from_sp_c**(c_sal, c_pressure, c_longitude, c_latitude)
!! calculates absolute salinity as a function of practical salinity, pressure, longitude and latitude.
!! It takes the following arguments as input:
!! * float c_sal - absolute salinity (g/kg)
!! * float c_pressure - pressure (dbar)
!! * float c_longitude - longitude (degrees)
!! * float c_latitude - latitude (degrees)
!!
!! and returns absolute salinity (g/kg) according to TEOS-10 function gsw_SA_from_SP.
!!

real(c_float) function gsw_sa_from_sp_c(c_sal, c_pressure, c_longitude, c_latitude) &
                                 bind(c, name='gsw_sa_from_sp_f90')
  use gsw_mod_toolbox, only: gsw_SA_from_SP
  implicit none
  real(c_float),  intent(in   ) :: c_sal, c_pressure, c_longitude, c_latitude
  real(kind_real) :: sal, pressure, longitude, latitude

  !! convert reals
  sal = real(c_sal, kind_real)
  pressure = real(c_pressure, kind_real)
  longitude = real(c_longitude, kind_real)
  latitude = real(c_latitude, kind_real)

  ! call GSW fortran function
  gsw_sa_from_sp_c = gsw_SA_from_SP(sal, pressure, longitude, latitude)
  return

end function gsw_sa_from_sp_c

end module ufo_oceanconversions_mod_c
