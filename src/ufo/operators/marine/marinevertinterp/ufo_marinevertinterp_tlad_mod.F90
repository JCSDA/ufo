! (C) Copyright 2017-2022 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

!> Fortran marinevertinterp module for tl/ad observation operator

module ufo_marinevertinterp_tlad_mod

use gsw_pot_to_insitu
use iso_c_binding
use kinds
use missing_values_mod
use obsspace_mod
use oops_variables_mod
use ufo_geovals_mod, only: ufo_geovals, ufo_geoval, ufo_geovals_get_var
use ufo_tpsp2ti_mod
use ufo_vars_mod
use vert_interp_mod

implicit none
private

!> Fortran derived type for the tl/ad observation operator
type, public :: ufo_marinevertinterp_tlad
 private
   type(oops_variables), public :: obsvars ! Variables to be simulated
   integer, allocatable, public :: obsvarindices(:) ! Indices of obsvars in the list of all
                                                  ! simulated variables in the ObsSpace
   type(oops_variables), public :: geovars
   integer                      :: nlocs      !< Number of observations
   integer                      :: nval       !< Number of level in model's profiles
   real(kind_real), allocatable :: wf(:)      !< Vertical interpolation weights
   integer,         allocatable :: wi(:)      !< Vertical interpolation indices
 contains
   procedure :: setup  => ufo_marinevertinterp_tlad_setup
   procedure :: delete  => ufo_marinevertinterp_tlad_delete
   procedure :: settraj => ufo_marinevertinterp_tlad_settraj
   procedure :: simobs_tl  => ufo_marinevertinterp_simobs_tl
   procedure :: simobs_ad  => ufo_marinevertinterp_simobs_ad
end type ufo_marinevertinterp_tlad

contains

! ------------------------------------------------------------------------------
subroutine ufo_marinevertinterp_tlad_setup(self)
   class(ufo_marinevertinterp_tlad), intent(inout) :: self

end subroutine ufo_marinevertinterp_tlad_setup

! ------------------------------------------------------------------------------
subroutine ufo_marinevertinterp_tlad_delete(self)
   class(ufo_marinevertinterp_tlad), intent(inout) :: self

   if (allocated(self%wi)) deallocate(self%wi)
   if (allocated(self%wf)) deallocate(self%wf)

end subroutine ufo_marinevertinterp_tlad_delete

! ------------------------------------------------------------------------------
subroutine ufo_marinevertinterp_tlad_settraj(self, geovals, obss)

   class(ufo_marinevertinterp_tlad), intent(inout)  :: self    !< Complete trajectory needed by the operator
   type(ufo_geovals), intent(in)                    :: geovals !< Model background
   type(c_ptr), value, intent(in)                   :: obss    !<

   integer :: ilev, iobs
   real (kind_real), allocatable :: depth(:,:)
   real(kind_real), allocatable :: obs_depth(:)
   type(ufo_geoval), pointer :: var, h

   ! make sure nothing allocated
   call self%delete()

   ! get h geovals
   call ufo_geovals_get_var(geovals, var_ocn_lay_thick, h)
   self%nval = h%nval
   self%nlocs = obsspace_get_nlocs(obss)

   ! calculate depth from layer thickness, TODO just ask model for depth??
   allocate(depth(h%nval, self%nlocs))
   do iobs = 1, self%nlocs
      !< Depth from layer thickness
      depth(1,iobs)=0.5*h%vals(1,iobs)
      do ilev = 2, h%nval
         depth(ilev,iobs)=sum(h%vals(1:ilev-1,iobs))+0.5*h%vals(ilev,iobs)
      end do
   end do

   ! Read in obs depth
   allocate(obs_depth(self%nlocs))
   call obsspace_get_db(obss, "MetaData", "depth", obs_depth)

   ! calculate interpolation weights
   allocate(self%wi(self%nlocs))
   allocate(self%wf(self%nlocs))
   do iobs = 1, self%nlocs
      call vert_interp_weights(h%nval, obs_depth(iobs), depth(:,iobs), self%wi(iobs), self%wf(iobs))
   end do

   ! done cleanup
   deallocate(depth)
   deallocate(obs_depth)

end subroutine ufo_marinevertinterp_tlad_settraj

! ------------------------------------------------------------------------------
subroutine ufo_marinevertinterp_simobs_tl(self, geovals, obss, nvars, nlocs, hofx )
   class(ufo_marinevertinterp_tlad), intent(in) :: self
   type(ufo_geovals),                intent(in) :: geovals
   integer,                          intent(in) :: nvars, nlocs
   real(c_double),                intent(inout) :: hofx(nvars, nlocs)
   type(c_ptr), value,               intent(in) :: obss

   integer :: iobsvar, ivar, iobs
   character(len=MAXVARLEN) :: geovar
   type(ufo_geoval), pointer :: profile

   do iobsvar = 1, size(self%obsvarindices)
      ! Get the index of the row of hofx to fill
      ivar = self%obsvarindices(iobsvar)

      ! Get the name of input variable in geovals
      geovar = self%geovars%variable(iobsvar)

      ! Get profile for this variable from geovals
      call ufo_geovals_get_var(geovals, geovar, profile)

      ! Interpolate from geovals to observational location into hofx
      do iobs = 1, nlocs
         call vert_interp_apply_tl(profile%nval, profile%vals(:,iobs), &
                                & hofx(ivar,iobs), self%wi(iobs), self%wf(iobs))
      enddo
   end do
end subroutine ufo_marinevertinterp_simobs_tl

! ------------------------------------------------------------------------------
subroutine ufo_marinevertinterp_simobs_ad(self, geovals, obss, nvars, nlocs, hofx )
   class(ufo_marinevertinterp_tlad),intent(in)    :: self
   type(ufo_geovals),               intent(inout) :: geovals
   integer,                         intent(in)    :: nvars, nlocs
   real(c_double),                  intent(in)    :: hofx(nvars, nlocs)
   type(c_ptr), value,              intent(in)    :: obss

   integer :: iobs, iobsvar, ivar
   type(ufo_geoval), pointer :: profile
   character(len=MAXVARLEN) :: geovar
   real(c_double) :: missing

   missing = missing_value(missing)

   do iobsvar = 1, size(self%obsvarindices)
     ! Get the index of the row of hofx to fill
     ivar = self%obsvarindices(iobsvar)

     ! Get the name of input variable in geovals
     geovar = self%geovars%variable(iobsvar)

     ! Get pointer to profile for this variable in geovals
     call ufo_geovals_get_var(geovals, geovar, profile)

     ! Adjoint of interpolate, from hofx into geovals
     do iobs = 1, self%nlocs
       if (hofx(ivar,iobs) /= missing) then
         call vert_interp_apply_ad(profile%nval, profile%vals(:,iobs), &
                                 & hofx(ivar,iobs), self%wi(iobs), self%wf(iobs))
       endif
     enddo
   enddo

end subroutine ufo_marinevertinterp_simobs_ad

end module ufo_marinevertinterp_tlad_mod
