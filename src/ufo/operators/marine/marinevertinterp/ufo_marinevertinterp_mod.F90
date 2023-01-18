! (C) Copyright 2017-2022 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

!> Fortran marinevertinterp module for observation operator

module ufo_marinevertinterp_mod

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

 integer, parameter :: max_string=800

!> Fortran derived type for the observation type
 type, public :: ufo_marinevertinterp
   type(oops_variables), public :: geovars
   type(oops_variables), public :: obsvars ! Variables to be simulated
   integer, allocatable, public :: obsvarindices(:) ! Indices of obsvars in the list of all
                                                    ! simulated variables in the ObsSpace
 contains
   procedure :: setup  => ufo_marinevertinterp_setup
   procedure :: simobs => ufo_marinevertinterp_simobs
 end type ufo_marinevertinterp

contains

! ------------------------------------------------------------------------------
subroutine ufo_marinevertinterp_setup(self)
  class(ufo_marinevertinterp), intent(inout) :: self

  call self%geovars%push_back("sea_water_cell_thickness")

end subroutine ufo_marinevertinterp_setup

! ------------------------------------------------------------------------------
subroutine ufo_marinevertinterp_simobs(self, geovals, obss, nvars, nlocs, hofx)
  class(ufo_marinevertinterp), intent(in) :: self
  type(ufo_geovals),           intent(in) :: geovals
  type(c_ptr),          value, intent(in) :: obss
  integer,                     intent(in) :: nvars, nlocs
  real(c_double),           intent(inout) :: hofx(nvars, nlocs)

  integer :: iobs, ilev, iobsvar, ivar
  type(ufo_geoval), pointer :: profile, h
  real (kind_real), allocatable :: depth(:,:)
  real(kind_real), allocatable :: obs_depth(:)
  integer, allocatable :: wi(:)
  real(kind_real), allocatable :: wf(:)
  character(len=MAXVARLEN) :: geovar

  ! get thickness TODO, get model level instead?
  call ufo_geovals_get_var(geovals, var_ocn_lay_thick, h)

  ! Read in obs depth
  allocate(obs_depth(nlocs))
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

  ! depths are no longer needed after this point
  deallocate(depth)
  deallocate(obs_depth)

  ! Vertical interpolation
  do iobsvar = 1, size(self%obsvarindices)
    ! get the index of row of hofx to fill
    ivar = self%obsvarindices(iobsvar)

    ! get name of input variable in geovals
    geovar = self%geovars%variable(iobsvar)

    ! get profile for this variable from geovals
    call ufo_geovals_get_var(geovals, geovar, profile)

    ! interpolate
    do iobs = 1, nlocs
      call vert_interp_apply(profile%nval, profile%vals(:,iobs), &
                             hofx(ivar, iobs), wi(iobs), wf(iobs))
    end do
  end do

  ! done, cleanup
  deallocate(wi)
  deallocate(wf)

end subroutine ufo_marinevertinterp_simobs

end module ufo_marinevertinterp_mod
