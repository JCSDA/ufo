! (C) Copyright 2017-2022 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

!> Fortran insitutemperature module for tl/ad observation operator

module ufo_insitutemperature_tlad_mod

use iso_c_binding
use kinds
use missing_values_mod
use obsspace_mod
use oops_variables_mod
use obs_variables_mod
use ufo_basis_tlad_mod, only: ufo_basis_tlad
use ufo_geovals_mod, only: ufo_geovals, ufo_geoval, ufo_geovals_get_var
use ufo_tpsp2ti_mod
use ufo_vars_mod
use vert_interp_mod

implicit none
private


 !> Fortran derived type for the tl/ad observation operator
 type, public :: ufo_insitutemperature_tlad
 private
   type(obs_variables), public :: obsvars ! Variables to be simulated
   integer,              public :: obsvaridx
   type(oops_variables), public :: geovars
   integer                            :: nlocs       !< Number of observations
   real (kind=kind_real), allocatable :: tempo(:)   !< temp interpolated at observation location
   real (kind=kind_real), allocatable :: salto(:)   !< salt interpolated at observation location
   real(kind_real), allocatable       :: wf(:)      !< Vertical interpolation weights
   integer, allocatable               :: wi(:)      !< Vertical interpolation indices
   real (kind=kind_real), allocatable :: jac(:,:)   !< Jacobian     [2 x nlocs]
 contains
  procedure :: setup  => ufo_insitutemperature_tlad_setup
  procedure :: delete  => ufo_insitutemperature_tlad_delete
  procedure :: settraj => ufo_insitutemperature_tlad_settraj
  procedure :: simobs_tl  => ufo_insitutemperature_simobs_tl
  procedure :: simobs_ad  => ufo_insitutemperature_simobs_ad
 end type ufo_insitutemperature_tlad

contains

! ------------------------------------------------------------------------------
subroutine ufo_insitutemperature_tlad_setup(self)
   class(ufo_insitutemperature_tlad), intent(inout) :: self

   call self%geovars%push_back("sea_water_potential_temperature")
   call self%geovars%push_back("sea_water_salinity")

end subroutine ufo_insitutemperature_tlad_setup

! ------------------------------------------------------------------------------
subroutine ufo_insitutemperature_tlad_delete(self)
   class(ufo_insitutemperature_tlad), intent(inout) :: self

   if (allocated(self%jac)) deallocate(self%jac)
   if (allocated(self%wi)) deallocate(self%wi)
   if (allocated(self%wf)) deallocate(self%wf)
   if (allocated(self%tempo)) deallocate(self%tempo)
   if (allocated(self%salto)) deallocate(self%salto)

end subroutine ufo_insitutemperature_tlad_delete

! ------------------------------------------------------------------------------
subroutine ufo_insitutemperature_tlad_settraj(self, geovals, obss)
   class(ufo_insitutemperature_tlad), intent(inout) :: self    !< Complete trajectory needed by the operator
   type(ufo_geovals),                 intent(in)    :: geovals !< Model background
   type(c_ptr), value,                intent(in)    :: obss    !< Insitu temperature observations

   type(ufo_geoval), pointer :: temp, salt, depth
   real(kind_real), allocatable :: obs_lon(:)
   real(kind_real), allocatable :: obs_lat(:)
   real(kind_real), allocatable :: obs_depth(:)
   integer :: ilev, iobs
   real(c_double) :: missing

   missing = missing_value(missing)

   ! make sure nothing allocated
   call self%delete()

   ! get geovals
   call ufo_geovals_get_var(geovals, var_ocn_pot_temp, temp)
   call ufo_geovals_get_var(geovals, var_ocn_salt, salt)
   call ufo_geovals_get_var(geovals, var_ocn_depth, depth)
   self%nlocs = obsspace_get_nlocs(obss)

   ! Read in obs data
   allocate(obs_lon(self%nlocs))
   allocate(obs_lat(self%nlocs))
   allocate(obs_depth(self%nlocs))
   call obsspace_get_db(obss, "MetaData", "longitude", obs_lon)
   call obsspace_get_db(obss, "MetaData", "latitude", obs_lat)
   call obsspace_get_db(obss, "MetaData", "depth", obs_depth)

   ! calculate interpolation weights
   allocate(self%wi(self%nlocs))
   allocate(self%wf(self%nlocs))
   do iobs = 1, self%nlocs
      call vert_interp_weights(depth%nval, obs_depth(iobs), depth%vals(:,iobs), self%wi(iobs), self%wf(iobs))
   end do

   ! Jacobian
   allocate(self%jac(2,self%nlocs))
   allocate(self%tempo(self%nlocs))
   allocate(self%salto(self%nlocs))
   outer: do iobs = 1, self%nlocs
      ! if any values in the the geovals profile are missing, skip
      ! TODO: be less restrictive if there is the possibility of partially 
      !  missing geoval profiles? (probably would never happen)
      do ilev = 1, temp%nval
         if (temp%vals(ilev,iobs) == missing .or. &
             salt%vals(ilev,iobs) == missing) then
            self%jac(:,iobs) = 0.0
            cycle outer
         end if
      end do

      ! Interpolate background to obs depth and save in traj
      call vert_interp_apply(temp%nval, temp%vals(:,iobs), self%tempo(iobs), self%wi(iobs), self%wf(iobs))
      call vert_interp_apply(salt%nval, salt%vals(:,iobs), self%salto(iobs), self%wi(iobs), self%wf(iobs))

      ! Compute jacobian
      call insitu_t_jac(self%jac(:,iobs), self%tempo(iobs), self%salto(iobs), obs_lon(iobs), obs_lat(iobs), obs_depth(iobs))
   end do outer

   deallocate(obs_lon, obs_lat, obs_depth)

end subroutine ufo_insitutemperature_tlad_settraj

! ------------------------------------------------------------------------------
subroutine ufo_insitutemperature_simobs_tl(self, geovals, obss, nvars, nlocs, hofx )
   use ufo_tpsp2ti_mod
   use gsw_pot_to_insitu
   use vert_interp_mod
   implicit none
   class(ufo_insitutemperature_tlad), intent(in)    :: self
   type(ufo_geovals),         intent(in) :: geovals
   integer,                   intent(in) :: nvars, nlocs
   real(c_double),         intent(inout) :: hofx(nvars, nlocs)
   type(c_ptr), value,        intent(in) :: obss

   type(ufo_geoval), pointer :: temp_d, salt_d !< Increments from geovals
   integer :: iobs
   real(kind_real) :: dtp, dsp
   real(kind_real), allocatable :: obs_lon(:), obs_lat(:), obs_depth(:)

   ! get geovals
   call ufo_geovals_get_var(geovals, var_ocn_pot_temp, temp_d)
   call ufo_geovals_get_var(geovals, var_ocn_salt, salt_d)

   ! Read in obs data
   allocate(obs_lon(self%nlocs))
   allocate(obs_lat(self%nlocs))
   allocate(obs_depth(self%nlocs))
   call obsspace_get_db(obss, "MetaData", "longitude", obs_lon)
   call obsspace_get_db(obss, "MetaData", "latitude", obs_lat)
   call obsspace_get_db(obss, "MetaData", "depth", obs_depth)

   ! interpolate and linear sea temperature profile obs operator
   do iobs = 1, self%nlocs
      call vert_interp_apply(temp_d%nval, temp_d%vals(:,iobs), dtp, self%wi(iobs), self%wf(iobs))
      call vert_interp_apply(salt_d%nval, salt_d%vals(:,iobs), dsp, self%wi(iobs), self%wf(iobs))
      call insitu_t_tl(hofx(self%obsvaridx, iobs), dtp, dsp, self%tempo(iobs), self%salto(iobs), obs_lon(iobs), obs_lat(iobs), obs_depth(iobs), self%jac(:,iobs))
   end do

end subroutine ufo_insitutemperature_simobs_tl

! ------------------------------------------------------------------------------
subroutine ufo_insitutemperature_simobs_ad(self, geovals, obss, nvars, nlocs, hofx )
   use ufo_tpsp2ti_mod
   use gsw_pot_to_insitu
   use vert_interp_mod
   implicit none
   class(ufo_insitutemperature_tlad), intent(in)    :: self
   type(ufo_geovals),         intent(inout) :: geovals
   integer,                   intent(in)    :: nvars, nlocs
   real(c_double),            intent(in)    :: hofx(nvars, nlocs)
   type(c_ptr), value,        intent(in)    :: obss

   integer:: iobs
   type(ufo_geoval), pointer :: dtemp, dsalt
   real (kind_real) :: dtp, dsp
   real(c_double) :: missing
   real(kind_real), allocatable :: obs_lon(:), obs_lat(:), obs_depth(:)

   !> Set missing value
   missing = missing_value(missing)

   ! Read in obs data
   allocate(obs_lon(self%nlocs))
   allocate(obs_lat(self%nlocs))
   allocate(obs_depth(self%nlocs))
   call obsspace_get_db(obss, "MetaData", "longitude", obs_lon)
   call obsspace_get_db(obss, "MetaData", "latitude", obs_lat)
   call obsspace_get_db(obss, "MetaData", "depth", obs_depth)

   ! get geovals
   call ufo_geovals_get_var(geovals, var_ocn_pot_temp, dtemp)
   call ufo_geovals_get_var(geovals, var_ocn_salt, dsalt)

   ! backward sea temperature profile obs operator
   do iobs = 1, self%nlocs
      if (hofx(self%obsvaridx, iobs) == missing) cycle

      ! adjoint obs operator
      dtp = 0.0
      dsp = 0.0
      call insitu_t_tlad(hofx(self%obsvaridx, iobs), dtp, dsp, self%tempo(iobs), self%salto(iobs), &
                         obs_lon(iobs), obs_lat(iobs), obs_depth(iobs), self%jac(:,iobs))

      ! backward interpolate
      call vert_interp_apply_ad(dtemp%nval, dtemp%vals(:,iobs), dtp, self%wi(iobs), self%wf(iobs))
      call vert_interp_apply_ad(dsalt%nval, dsalt%vals(:,iobs), dsp, self%wi(iobs), self%wf(iobs))

   end do
end subroutine ufo_insitutemperature_simobs_ad

end module ufo_insitutemperature_tlad_mod
