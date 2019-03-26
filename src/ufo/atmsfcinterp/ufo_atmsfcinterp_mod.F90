! (C) Copyright 2017-2018 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

!> Fortran module for atmsfcinterp observation operator

module ufo_atmsfcinterp_mod

  use iso_c_binding
  use config_mod
  use kinds

  use ufo_geovals_mod, only: ufo_geovals, ufo_geoval, ufo_geovals_get_var
  use ufo_basis_mod, only: ufo_basis
  use ufo_vars_mod
  use obsspace_mod

  implicit none
  private

  integer, parameter :: max_string = 800
  real(kind_real), parameter :: grav = 9.80665e+0_kind_real
  !> Fortran derived type for the observation type
  type, extends(ufo_basis), public :: ufo_atmsfcinterp
  private
    integer :: nvars
    real(kind_real) :: magl
    character(len=max_string), allocatable :: varin(:)
    character(len=max_string), allocatable :: varout(:)
  end type

  contains
    procedure :: setup  => ufo_atmsfcinterp_setup
    procedure :: delete => ufo_atmsfcinterp_delete
    procedure :: simobs => ufo_atmsfcinterp_simobs
  end type ufo_atmsfcinterp

contains

! ------------------------------------------------------------------------------
subroutine ufo_atmsfcinterp_setup(self, c_conf)
  implicit none
  class(ufo_atmsfcinterp), intent(inout) :: self
  type(c_ptr),        intent(in)    :: c_conf

  !> Size of variables
  self%nvars = size(config_get_string_vector(c_conf, max_string, "variables"))
  !> Allocate varout: variables in the observation vector
  allocate(self%varout(self%nvars))
  !> Read variable list and store in varout
  self%varout = config_get_string_vector(c_conf, max_string, "variables")
  !> Allocate varin: variables we need from the model
  allocate(self%varin(self%nvars))
  !> Set vars_in based on vars_out
  do ii = 1, self%nvars
    self%varin(ii) = self%varout(ii)
  enddo

end subroutine ufo_atmsfcinterp_setup

! ------------------------------------------------------------------------------
! TODO: add cleanup of your observation operator (optional)
subroutine ufo_atmsfcinterp_delete(self)
implicit none
class(ufo_atmsfcinterp), intent(inout) :: self

end subroutine ufo_atmsfcinterp_delete

! ------------------------------------------------------------------------------
subroutine ufo_atmsfcinterp_simobs(self, geovals, hofx, obss)
  implicit none
  class(ufo_atmsfcinterp), intent(in)    :: self
  type(ufo_geovals),  intent(in)    :: geovals
  real(c_double),     intent(inout) :: hofx(:)
  type(c_ptr), value, intent(in)    :: obss
  type(ufo_geoval), pointer :: phi, hgt
  integer :: nlocs,i
  real(kind_real), allocatable :: obselev(:), obshgt(:)

  ! for low altitudes, we can assume: phi = z
  ! get geopotential height profile
  call ufo_geovals_get_var(geovals, var_z, phi)
  ! get surface geopotential height
  call ufo_geovals_get_var(geovals, var_sfc_z, hgt)

  ! get number of obs
  nlocs = obsspace_get_nobs(obss)

  ! get station elevation from obs
  allocate(obselev(nlocs))
  call obsspace_get_db(obss, "MetaData", "station_elevation",obselev)

  ! get observation height (above sea level)
  allocate(obshgt(nlocs))
  call obsspace_get_db(obss, "MetaData", "height",obshgt)

  print *, shape(phi), shape(hgt)
  do i=1,nlocs
   print *, obselev(i),obshgt(i)
  end do



end subroutine ufo_atmsfcinterp_simobs

! ------------------------------------------------------------------------------

end module ufo_atmsfcinterp_mod
