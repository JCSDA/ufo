! (C) Copyright 2017-2022 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

!> Fortran insitutemperature module for observation operator

module ufo_insitutemperature_mod

use gsw_pot_to_insitu
use iso_c_binding
use kinds
use obsspace_mod
use oops_variables_mod
use obs_variables_mod
use ufo_geovals_mod, only: ufo_geovals, ufo_geoval, ufo_geovals_get_var
use ufo_tpsp2ti_mod
use ufo_vars_mod
use vert_interp_mod

implicit none
private

!> Fortran derived type for the observation type
type, public :: ufo_insitutemperature
   type(oops_variables), public :: geovars
   type(obs_variables), public :: obsvars ! Variables to be simulated
   integer,              public :: obsvaridx
contains
   procedure :: setup  => ufo_insitutemperature_setup
   procedure :: simobs => ufo_insitutemperature_simobs
end type ufo_insitutemperature

contains

! ------------------------------------------------------------------------------
subroutine ufo_insitutemperature_setup(self)
   class(ufo_insitutemperature), intent(inout) :: self

   call self%geovars%push_back(var_ocn_pot_temp)
   call self%geovars%push_back(var_ocn_salt)
   call self%geovars%push_back(var_ocn_depth)

end subroutine ufo_insitutemperature_setup


! ------------------------------------------------------------------------------
subroutine ufo_insitutemperature_simobs(self, geovals, obss, nvars, nlocs, hofx)
   use, intrinsic :: iso_c_binding, only: c_double, c_ptr
   implicit none

   class(ufo_insitutemperature), intent(in) :: self
   type(ufo_geovals),            intent(in) :: geovals
   type(c_ptr),           value, intent(in) :: obss
   integer,                      intent(in) :: nvars, nlocs
   real(c_double),            intent(inout) :: hofx(nvars, nlocs)

   integer :: iobs, ilev
   type(ufo_geoval), pointer :: temp, salt, depth
   real(kind_real), allocatable :: obs_lon(:)
   real(kind_real), allocatable :: obs_lat(:)
   real(kind_real), allocatable :: obs_depth(:)
   real(kind_real) :: tp, sp
   integer, allocatable :: wi(:)
   real(kind_real), allocatable :: wf(:)
   real(c_double) :: missing

   missing = missing_value(missing)

   ! Associate geoval pointers
    call ufo_geovals_get_var(geovals, var_ocn_pot_temp, temp)
    call ufo_geovals_get_var(geovals, var_ocn_salt, salt)
    call ufo_geovals_get_var(geovals, var_ocn_depth, depth)

    ! Read in obs data
    allocate(obs_lon(nlocs))
    allocate(obs_lat(nlocs))
    allocate(obs_depth(nlocs))
    call obsspace_get_db(obss, "MetaData", "longitude", obs_lon)
    call obsspace_get_db(obss, "MetaData", "latitude", obs_lat)
    call obsspace_get_db(obss, "MetaData", "depth", obs_depth)

   ! calculate interpolation weights
   allocate(wi(nlocs))
   allocate(wf(nlocs))
   do iobs = 1, nlocs
      call vert_interp_weights(depth%nval, obs_depth(iobs), depth%vals(:,iobs), wi(iobs), wf(iobs))
   end do

   outer: do iobs = 1, nlocs
      ! if any values in the the geovals profile are missing, skip
      ! TODO: be less restrictive if there is the possibility of partially 
      !  missing geoval profiles? (probably would never happen)
      do ilev = 1, temp%nval
         if (temp%vals(ilev,iobs) == missing .or. &
             salt%vals(ilev,iobs) == missing) then
            hofx(self%obsvaridx, iobs) = missing
            cycle outer
         end if
      end do

       ! Interpolate temp_p, salt_p to deptho
       call vert_interp_apply(temp%nval, temp%vals(:,iobs), tp, wi(iobs), wf(iobs))
       call vert_interp_apply(salt%nval, salt%vals(:,iobs), sp, wi(iobs), wf(iobs))

       call insitu_t_nl(hofx(self%obsvaridx, iobs), tp, sp, obs_lon(iobs), obs_lat(iobs), obs_depth(iobs))
   end do outer

   ! done, cleanup

   deallocate(obs_lon)
   deallocate(obs_lat)
   deallocate(obs_depth)
   deallocate(wi)
   deallocate(wf)

end subroutine ufo_insitutemperature_simobs

end module ufo_insitutemperature_mod
