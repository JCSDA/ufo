! (C) Copyright 2017-2018 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

module ufo_conventional_profile_mod

  use iso_c_binding
  use kinds
  use ufo_vars_mod
  use ufo_geovals_mod
  use ufo_geovals_mod_c,   only: ufo_geovals_registry
  use vert_interp_mod
  use ufo_basis_mod, only: ufo_basis
  use obsspace_mod

  integer, parameter :: max_string=800

  type, extends(ufo_basis) :: ufo_conventional_profile
   private
     integer :: nval, nobs
     real(kind_real), allocatable :: wf(:)
     integer, allocatable :: wi(:)
  contains
    procedure :: simobs    => conventional_profile_simobs_
  end type ufo_conventional_profile
contains

! ------------------------------------------------------------------------------

    subroutine conventional_profile_simobs_(self, geovals, hofx, obss)

      implicit none
      class(ufo_conventional_profile), intent(in)  :: self
      type(ufo_geovals), intent(in)                :: geovals
      real(c_double),  intent(inout)               :: hofx(:)
      type(c_ptr), value, intent(in)               :: obss

      character(len=*), parameter :: myname_="ufo_conventional_profile_simobs"
      character(max_string) :: err_msg

      integer :: iobs
      real(kind_real) :: wf
      integer :: wi, nobs
      real(kind_real), allocatable :: pressure(:)
      type(ufo_geoval), pointer :: prsl, tv

      ! get pressure from geovals
      call ufo_geovals_get_var(geovals, var_prsl, prsl)

      ! get tv from geovals
      call ufo_geovals_get_var(geovals, var_tv, tv)

      ! observation of pressure (for vertical interpolation)
      nobs = obsspace_get_nobs(obss)
      allocate(pressure(nobs))
      call obsspace_get_db(obss, "MetaData", "air_pressure", pressure)

      ! obs operator
      do iobs = 1, geovals%nobs
        call vert_interp_weights(prsl%nval,log(pressure(iobs)/10.),prsl%vals(:,iobs),wi,wf)
        call vert_interp_apply(tv%nval, tv%vals(:,iobs), hofx(iobs), wi, wf)
      enddo

      ! cleanup
      deallocate(pressure)
    end subroutine conventional_profile_simobs_

! ------------------------------------------------------------------------------

end module ufo_conventional_profile_mod
