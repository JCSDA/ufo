! (C) Copyright 2017-2020 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

!> Fortran module for avgkernel tl/ad observation operator

module ufo_avgkernel_tlad_mod
  use oops_variables_mod
  use ufo_vars_mod
  use ufo_geovals_mod
  use missing_values_mod
  use kinds
  use iso_c_binding

  implicit none
  private
  integer, parameter :: max_string=800

! ------------------------------------------------------------------------------

 type, public :: ufo_avgkernel_tlad
 private
   type(oops_variables), public :: obsvars
   type(oops_variables), public :: geovars
   integer :: nlayers_kernel
   integer :: nlocs, nvars, nval
   character(kind=c_char,len=:), allocatable :: obskernelvar, obspressurevar, tracervars(:)
   logical :: troposphere, totalcolumn
   real(kind_real) :: convert_factor_model
   real(kind_real), allocatable, dimension(:,:) :: avgkernel_obs, prsi_obs
   real(kind_real), allocatable, dimension(:,:) :: prsi
   real(kind_real), allocatable, dimension(:) :: airmass_tot, airmass_trop
   integer, allocatable, dimension(:) :: troplev_obs
   logical :: flip_it
 contains
   procedure :: setup  => avgkernel_tlad_setup_
   procedure :: cleanup  => avgkernel_tlad_cleanup_
   procedure :: settraj => avgkernel_tlad_settraj_
   procedure :: simobs_tl  => avgkernel_simobs_tl_
   procedure :: simobs_ad  => avgkernel_simobs_ad_
   final :: destructor
 end type ufo_avgkernel_tlad

contains

! ------------------------------------------------------------------------------
subroutine avgkernel_tlad_setup_(self, f_conf)
  use fckit_configuration_module, only: fckit_configuration
  use ufo_constants_mod, only: one
  implicit none
  class(ufo_avgkernel_tlad), intent(inout) :: self
  type(fckit_configuration), intent(in)  :: f_conf
  integer :: nlevs_yaml
  integer :: ivar, nvars
  character(len=max_string) :: err_msg
  character(len=:), allocatable :: str_array(:)

  ! get configuration for the averaging kernel operator
  call f_conf%get_or_die("nlayers_kernel", self%nlayers_kernel)

  ! get variable name from IODA for observation averaging kernel
  call f_conf%get_or_die("AvgKernelVar", self%obskernelvar)

  ! get vertical pressure grid 
  call f_conf%get_or_die("PresLevVar", self%obspressurevar)
 
  ! get name of geoval/tracer to use from the model
  nvars = self%obsvars%nvars()
  call f_conf%get_or_die("tracer variables", str_array)
  self%tracervars = str_array

  ! determine if this is a total column or troposphere calculation
  call f_conf%get_or_die("tropospheric column", self%troposphere)
  call f_conf%get_or_die("total column", self%totalcolumn)

  ! both of these cannot be true
  if (self%troposphere .and. self%totalcolumn) then
    write(err_msg, *) "ufo_avgkernel_tlad_setup error: both tropospheric and total column set to TRUE, only one can be TRUE"
    call abor1_ftn(err_msg)
  end if

  ! do we need a conversion factor, say between ppmv and unity?
  call f_conf%get_or_die("model units coeff", self%convert_factor_model)

  ! add variables to geovars that are needed
  ! specified tracers
  do ivar = 1, nvars
    call self%geovars%push_back(self%tracervars(ivar))
  end do

end subroutine avgkernel_tlad_setup_

! ------------------------------------------------------------------------------
subroutine destructor(self)
  implicit none
  type(ufo_avgkernel_tlad), intent(inout) :: self

  call self%cleanup()

end subroutine destructor

! ------------------------------------------------------------------------------
subroutine avgkernel_tlad_settraj_(self, geovals_in, obss)
  use iso_c_binding
  use ufo_geovals_mod, only: ufo_geovals, ufo_geoval, ufo_geovals_get_var
  use ufo_constants_mod, only: half, zero
  use obsspace_mod
  implicit none
  class(ufo_avgkernel_tlad), intent(inout) :: self
  type(ufo_geovals),       intent(in)    :: geovals_in
  type(c_ptr), value,      intent(in)    :: obss

  ! Local variables
  type(ufo_geoval), pointer :: prsi
  integer :: ivar, iobs, ilev
  character(len=MAXVARLEN) :: varstring
  character(len=4) :: levstr
  type(ufo_geovals) :: geovals

  ! get nlocs and nvars
  self%nlocs = obsspace_get_nlocs(obss)
  self%nvars = obsspace_get_nvars(obss)

  ! get geovals of atmospheric pressure
  call ufo_geovals_copy(geovals_in, geovals)  ! dont want to change geovals_in
  call ufo_geovals_get_var(geovals, var_prsi, prsi)

  ! get nval
  self%nval = prsi%nval

  allocate(self%prsi(self%nval, self%nlocs))
  do iobs = 1, self%nlocs

    self%prsi(:,iobs) = prsi%vals(:,iobs)
  end do

  ! grab necesary metadata from IODA
  ! get observation averaging kernel
  ! once 2D arrays are allowed, rewrite/simplify this part
  ! TEMPORARY: reverse do loops to make sure we follow the convention
  ! TEMPORARY: top->bottom; increasing pressure
  allocate(self%avgkernel_obs(self%nlayers_kernel, self%nlocs))
  do ilev = self%nlayers_kernel, 1, -1
    write(levstr, fmt = "(I3)") ilev
    levstr = adjustl(levstr)
    varstring = trim(self%obskernelvar)//"_"//trim(levstr)
    call obsspace_get_db(obss, "RtrvlAncData", trim(varstring), &
            self%avgkernel_obs(self%nlayers_kernel+1-ilev, :))
  end do

  ! get prsi_obs
  allocate(self%prsi_obs(self%nlayers_kernel+1, self%nlocs))
  do ilev = self%nlayers_kernel, 1, -1
    write(levstr, fmt = "(I3)") ilev
    levstr = adjustl(levstr)
    varstring = trim(self%obspressurevar)//"_"//trim(levstr)
    call obsspace_get_db(obss, "RtrvlAncData", trim(varstring), &
            self%prsi_obs(self%nlayers_kernel+2-ilev, :))
  end do
  !last vertices should be always TOA (0 hPa)
  self%prsi_obs(1,:) = zero
  self%prsi(1,:) = zero

  if (self%troposphere) then
    allocate(self%troplev_obs(self%nlocs))
    allocate(self%airmass_trop(self%nlocs))
    allocate(self%airmass_tot(self%nlocs))
    call obsspace_get_db(obss, "RtrvlAncData", "troposphere_layer_index", self%troplev_obs)
    call obsspace_get_db(obss, "RtrvlAncData", "air_mass_factor_troposphere", self%airmass_trop)
    call obsspace_get_db(obss, "RtrvlAncData", "air_mass_factor_total", self%airmass_tot)
  end if

end subroutine avgkernel_tlad_settraj_

! ------------------------------------------------------------------------------
subroutine avgkernel_simobs_tl_(self, geovals_in, obss, nvars, nlocs, hofx)
  use iso_c_binding
  use ufo_geovals_mod, only: ufo_geovals, ufo_geoval, ufo_geovals_get_var
  use obsspace_mod
  use satcolumn_mod, only: simulate_column_ob_tl
  implicit none
  class(ufo_avgkernel_tlad), intent(in)    :: self
  type(ufo_geovals),       intent(in)    :: geovals_in
  integer,                 intent(in)    :: nvars, nlocs
  real(c_double),          intent(inout) :: hofx(nvars, nlocs)
  type(c_ptr), value,      intent(in)    :: obss
  type(ufo_geoval), pointer :: tracer
  integer :: ivar, iobs
  character(len=MAXVARLEN) :: geovar
  real(kind_real) :: hofx_tmp
  type(ufo_geovals) :: geovals
  real(c_double) :: missing

  missing = missing_value(missing)

  call ufo_geovals_copy(geovals_in, geovals)  ! dont want to change geovals_in

  ! loop through all variables
  do ivar = 1, nvars
    geovar = self%tracervars(ivar)
    call ufo_geovals_get_var(geovals, geovar, tracer)
    do iobs = 1, nlocs
      if (self%avgkernel_obs(1,iobs) /= missing) then ! take care of missing obs
        if(self%flip_it) tracer%vals(1:tracer%nval,iobs) = tracer%vals(tracer%nval:1:-1,iobs)
        if (self%troposphere) then
          call simulate_column_ob_tl(self%nlayers_kernel, tracer%nval, self%avgkernel_obs(:,iobs), &
                                     self%prsi_obs(:,iobs), self%prsi(:,iobs),&
                                     tracer%vals(:,iobs)*self%convert_factor_model, &
                                     hofx_tmp, self%troplev_obs(iobs), self%airmass_tot(iobs), self%airmass_trop(iobs))
        else if (self%totalcolumn) then
          call simulate_column_ob_tl(self%nlayers_kernel, tracer%nval, self%avgkernel_obs(:,iobs), &
                                     self%prsi_obs(:,iobs), self%prsi(:,iobs),&
                                     tracer%vals(:,iobs)*self%convert_factor_model, &
                                     hofx_tmp)
        end if
        hofx(ivar,iobs) = hofx_tmp
      else
        hofx(ivar,iobs) = missing
      end if
    end do
  end do
  call ufo_geovals_delete(geovals)

end subroutine avgkernel_simobs_tl_

! ------------------------------------------------------------------------------
subroutine avgkernel_simobs_ad_(self, geovals_in, obss, nvars, nlocs, hofx)
  use iso_c_binding
  use ufo_geovals_mod, only: ufo_geovals, ufo_geoval, ufo_geovals_get_var
  use satcolumn_mod, only: simulate_column_ob_ad
  use obsspace_mod
  implicit none
  class(ufo_avgkernel_tlad), intent(in)    :: self
  type(ufo_geovals),       intent(inout) :: geovals_in
  integer,                 intent(in)    :: nvars, nlocs
  real(c_double),          intent(in)    :: hofx(nvars, nlocs)
  type(c_ptr), value,      intent(in)    :: obss
  type(ufo_geoval), pointer :: tracer
  character(len=MAXVARLEN) :: geovar
  type(ufo_geovals) :: geovals
  real(kind_real) :: hofx_tmp
  integer :: ivar, iobs
  real(c_double) :: missing

  missing = missing_value(missing)

  if (.not. geovals_in%linit ) geovals_in%linit=.true. ! need this for var exe

  ! loop through all variables
  do ivar = 1, nvars
    geovar = self%tracervars(ivar)
    call ufo_geovals_get_var(geovals_in, geovar, tracer)

    do iobs = 1, nlocs
      if (hofx(ivar,iobs) /= missing) then ! take care of missing obs
        hofx_tmp = hofx(ivar,iobs)
        if (self%troposphere) then
          call simulate_column_ob_ad(self%nlayers_kernel, tracer%nval, self%avgkernel_obs(:,iobs), &
                                     self%prsi_obs(:,iobs), self%prsi(:,iobs),&
                                     tracer%vals(:,iobs), &
                                     hofx_tmp, self%troplev_obs(iobs), self%airmass_tot(iobs), self%airmass_trop(iobs))
        else if (self%totalcolumn) then
          call simulate_column_ob_ad(self%nlayers_kernel, tracer%nval, self%avgkernel_obs(:,iobs), &
                                     self%prsi_obs(:,iobs), self%prsi(:,iobs),&
                                     tracer%vals(:,iobs), hofx_tmp)
        end if
        tracer%vals(:,iobs) = tracer%vals(:,iobs) * self%convert_factor_model
      end if
      if(self%flip_it) tracer%vals(1:tracer%nval,iobs) = tracer%vals(tracer%nval:1:-1,iobs)
    end do
  end do

end subroutine avgkernel_simobs_ad_

! ------------------------------------------------------------------------------
subroutine avgkernel_tlad_cleanup_(self)
  implicit none
  class(ufo_avgkernel_tlad), intent(inout) :: self
  self%nvars = 0
  self%nlocs = 0
  self%nval = 0
  ! deallocate things
  if (allocated(self%obskernelvar)) deallocate(self%obskernelvar)
  if (allocated(self%obspressurevar)) deallocate(self%obspressurevar)
  if (allocated(self%tracervars)) deallocate(self%tracervars)
  if (allocated(self%avgkernel_obs)) deallocate(self%avgkernel_obs)
  if (allocated(self%prsi_obs)) deallocate(self%prsi_obs)
  if (allocated(self%prsi)) deallocate(self%prsi)
  if (allocated(self%airmass_tot)) deallocate(self%airmass_tot)
  if (allocated(self%airmass_trop)) deallocate(self%airmass_trop)
  if (allocated(self%troplev_obs)) deallocate(self%troplev_obs)
end subroutine avgkernel_tlad_cleanup_

! ------------------------------------------------------------------------------

end module ufo_avgkernel_tlad_mod
