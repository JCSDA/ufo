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
use ufo_geovals_mod, only: ufo_geovals, ufo_geoval, ufo_geovals_get_var
use ufo_tpsp2ti_mod
use ufo_vars_mod
use vert_interp_mod

implicit none
private

!> Fortran derived type for the observation type
type, public :: ufo_insitutemperature
   type(oops_variables), public :: geovars
   type(oops_variables), public :: obsvars ! Variables to be simulated
   integer,              public :: obsvaridx
contains
   procedure :: setup  => ufo_insitutemperature_setup
   procedure :: simobs => ufo_insitutemperature_simobs
end type ufo_insitutemperature

contains

! ------------------------------------------------------------------------------
subroutine ufo_insitutemperature_setup(self)
   class(ufo_insitutemperature), intent(inout) :: self

   call self%geovars%push_back("sea_water_potential_temperature")
   call self%geovars%push_back("sea_water_salinity")
   call self%geovars%push_back("sea_water_cell_thickness")

end subroutine ufo_insitutemperature_setup


! ------------------------------------------------------------------------------
subroutine ufo_insitutemperature_simobs(self, geovals, obss, nvars, nlocs, hofx)
   class(ufo_insitutemperature), intent(in) :: self
   type(ufo_geovals),            intent(in) :: geovals
   type(c_ptr),           value, intent(in) :: obss
   integer,                      intent(in) :: nvars, nlocs
   real(c_double),            intent(inout) :: hofx(nvars, nlocs)

   integer :: iobs, ilev
   type(ufo_geoval), pointer :: temp, salt, h
   real (kind_real), allocatable :: depth(:,:)
   real(kind_real), allocatable :: obs_lon(:)
   real(kind_real), allocatable :: obs_lat(:)
   real(kind_real), allocatable :: obs_depth(:)
   real(kind_real) :: tp, sp
   integer, allocatable :: wi(:)
   real(kind_real), allocatable :: wf(:)

   ! Associate geoval pointers
    call ufo_geovals_get_var(geovals, var_ocn_pot_temp, temp)
    call ufo_geovals_get_var(geovals, var_ocn_salt, salt)
    call ufo_geovals_get_var(geovals, var_ocn_lay_thick, h)

    ! Read in obs data
    allocate(obs_lon(nlocs))
    allocate(obs_lat(nlocs))
    allocate(obs_depth(nlocs))
    call obsspace_get_db(obss, "MetaData", "longitude", obs_lon)
    call obsspace_get_db(obss, "MetaData", "latitude", obs_lat)
    call obsspace_get_db(obss, "MetaData", "depth", obs_depth)

   ! calculate depth from layer thickness, TODO just ask model for depth??
   allocate(depth(h%nval, nlocs))
   do iobs = 1, nlocs
      !< Depth from layer thickness
      depth(1,iobs)=0.5*h%vals(1,iobs)
      do ilev = 2, h%nval
         depth(ilev,iobs)=sum(h%vals(1:ilev-1,iobs))+0.5*h%vals(ilev,iobs)
      end do
   end do

   ! calculate interpolation weights
   allocate(wi(nlocs))
   allocate(wf(nlocs))
   do iobs = 1, nlocs
      call vert_interp_weights(h%nval, obs_depth(iobs), depth(:,iobs), wi(iobs), wf(iobs))
   end do

   deallocate(depth)


   do iobs = 1, nlocs
       ! Interpolate temp_p, salt_p to deptho
       call vert_interp_apply(temp%nval, temp%vals(:,iobs), tp, wi(iobs), wf(iobs))
       call vert_interp_apply(salt%nval, salt%vals(:,iobs), sp, wi(iobs), wf(iobs))

       call insitu_t_nl(hofx(self%obsvaridx, iobs), tp, sp, obs_lon(iobs), obs_lat(iobs), obs_depth(iobs))
   end do

   ! done, cleanup

   deallocate(obs_lon)
   deallocate(obs_lat)
   deallocate(obs_depth)
   deallocate(wi)
   deallocate(wf)

end subroutine ufo_insitutemperature_simobs

end module ufo_insitutemperature_mod
