! (C) Copyright 2017-2020 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

!> Fortran module for averaging kernel observation operator

module ufo_avgkernel_mod

 use oops_variables_mod
 use ufo_vars_mod
 use missing_values_mod
 use kinds
 use iso_c_binding

 implicit none
 private
 integer, parameter :: max_string=800

!> Fortran derived type for the observation type
 type, public :: ufo_avgkernel
 private
   type(oops_variables), public :: obsvars
   type(oops_variables), public :: geovars
   integer :: nlayers_kernel
   real(kind_real), allocatable, dimension(:) :: ak_kernel, bk_kernel
   character(kind=c_char,len=:), allocatable :: obskernelvar, tracervars(:)
   logical :: troposphere, totalcolumn
   real(kind_real) :: convert_factor_model, convert_factor_hofx
 contains
   procedure :: setup  => ufo_avgkernel_setup
   procedure :: simobs => ufo_avgkernel_simobs
   final :: destructor
 end type ufo_avgkernel

contains

! ------------------------------------------------------------------------------
subroutine ufo_avgkernel_setup(self, f_conf)
  use fckit_configuration_module, only: fckit_configuration
  use ufo_constants_mod, only: one
  implicit none
  class(ufo_avgkernel), intent(inout)     :: self
  type(fckit_configuration), intent(in) :: f_conf
  type(fckit_configuration) :: f_confOpts
  integer :: nlevs_yaml
  integer :: ivar, nvars
  character(len=max_string) :: err_msg
  character(len=:), allocatable :: str_array(:)

  ! get configuration for the averaging kernel operator
  call f_conf%get_or_die("obs options",f_confOpts)

  call f_confOpts%get_or_die("nlayers_kernel", self%nlayers_kernel)
  nlevs_yaml = f_confOpts%get_size("ak")
  if (nlevs_yaml /= self%nlayers_kernel+1) then
    write(err_msg, *) "ufo_avgkernel_setup error: YAML nlayers_kernel != size of ak array"
    call abor1_ftn(err_msg)
  end if
  nlevs_yaml = f_confOpts%get_size("bk")
  if (nlevs_yaml /= self%nlayers_kernel+1) then
    write(err_msg, *) "ufo_avgkernel_setup error: YAML nlayers_kernel != size of bk array"
    call abor1_ftn(err_msg)
  end if

  ! allocate and read ak/bk for the averaging kernel
  allocate(self%ak_kernel(nlevs_yaml))
  allocate(self%bk_kernel(nlevs_yaml))
  call f_confOpts%get_or_die("ak", self%ak_kernel)
  call f_confOpts%get_or_die("bk", self%bk_kernel)

  ! get variable name from IODA for observation averaging kernel
  if (.not. f_confOpts%get("AvgKernelVar", self%obskernelvar)) then
    self%obskernelvar = "averaging_kernel" ! default option
  end if

  ! get name of geoval/tracer to use from the model
  nvars = self%obsvars%nvars()
  call f_confOpts%get_or_die("tracer variables", str_array)
  self%tracervars = str_array

  ! determine if this is a total column or troposphere calculation
  ! support stratosphere, etc. later?
  if (.not. f_confOpts%get("tropospheric column", self%troposphere)) self%troposphere = .false.
  if (.not. f_confOpts%get("total column", self%totalcolumn)) self%totalcolumn = .false.

  ! do we need a conversion factor, say between ppmv and unity?
  if (.not. f_confOpts%get("model units coeff", self%convert_factor_model)) self%convert_factor_model = one
  if (.not. f_confOpts%get("hofx units coeff", self%convert_factor_hofx)) self%convert_factor_hofx = one

  ! add variables to geovars that are needed
  ! specified tracers
  do ivar = 1, nvars
    call self%geovars%push_back(self%tracervars(ivar))
  end do
  ! surface pressure
  call self%geovars%push_back(var_ps)
  ! column pressure both layer and interface
  call self%geovars%push_back(var_prs)
  call self%geovars%push_back(var_prsi)
  ! need air temperature for number density conversion
  call self%geovars%push_back(var_ts)
  ! need geopotential height on levels to convert units
  call self%geovars%push_back(var_zi)

end subroutine ufo_avgkernel_setup

! ------------------------------------------------------------------------------
subroutine destructor(self)
  implicit none
  type(ufo_avgkernel), intent(inout) :: self

end subroutine destructor

! ------------------------------------------------------------------------------
! averaging kernel observation operator
subroutine ufo_avgkernel_simobs(self, geovals_in, obss, nvars, nlocs, hofx)
  use kinds
  use ufo_geovals_mod, only: ufo_geovals, ufo_geoval, ufo_geovals_get_var, &
                             ufo_geovals_reorderzdir, ufo_geovals_copy
  use ufo_constants_mod, only: half
  use satcolumn_mod, only: simulate_column_ob
  use iso_c_binding
  use obsspace_mod
  implicit none
  class(ufo_avgkernel), intent(in)    :: self
  integer, intent(in)               :: nvars, nlocs
  type(ufo_geovals),  intent(in)    :: geovals_in
  real(c_double),     intent(inout) :: hofx(nvars, nlocs)
  type(c_ptr), value, intent(in)    :: obss

  ! Local variables
  type(ufo_geoval), pointer :: prsi, prsl, psfc, temp, phii, tracer
  integer :: ivar, iobs, ilev
  character(len=MAXVARLEN) :: geovar, varstring
  character(len=4) :: levstr
  real(kind_real), allocatable, dimension(:,:) :: avgkernel_obs, prsl_obs, prsi_obs
  real(kind_real), allocatable, dimension(:) :: airmass_tot, airmass_trop
  integer, allocatable, dimension(:) :: troplev_obs
  real(kind_real) :: hofx_tmp
  type(ufo_geovals) :: geovals
  real(c_double) :: missing

  missing = missing_value(missing)

  ! get geovals of atmospheric pressure
  call ufo_geovals_copy(geovals_in, geovals)  ! dont want to change geovals_in
  call ufo_geovals_reorderzdir(geovals, self%geovars%variable(nvars+2), "bottom2top")
  call ufo_geovals_get_var(geovals, self%geovars%variable(nvars+1), psfc)
  call ufo_geovals_get_var(geovals, self%geovars%variable(nvars+2), prsl)
  call ufo_geovals_get_var(geovals, self%geovars%variable(nvars+3), prsi)
  call ufo_geovals_get_var(geovals, self%geovars%variable(nvars+4), temp)
  call ufo_geovals_get_var(geovals, self%geovars%variable(nvars+5), phii)

  ! grab necesary metadata from IODA
  ! get observation averaging kernel
  ! once 2D arrays are allowed, rewrite/simplify this part
  allocate(avgkernel_obs(self%nlayers_kernel, nlocs))
  do ilev = 1, self%nlayers_kernel
    write(levstr, fmt = "(I3)") ilev
    levstr = adjustl(levstr)
    varstring = trim(self%obskernelvar)//"_"//trim(levstr)
    call obsspace_get_db(obss, "MetaData", trim(varstring), avgkernel_obs(ilev, :))
  end do

  ! compute prsl_obs/prsi_obs from ak/bk/psfc
  allocate(prsl_obs(self%nlayers_kernel, nlocs))
  allocate(prsi_obs(self%nlayers_kernel+1, nlocs))
  ! prsi_obs calculation
  do ilev = 1, self%nlayers_kernel+1
    prsi_obs(ilev,:) = self%ak_kernel(ilev) + self%bk_kernel(ilev) * psfc%vals(1,:)
  end do
  ! using simple averaging for now for prsl, can use more complex way later
  do ilev = 1, self%nlayers_kernel
    prsl_obs(ilev,:) = (prsi_obs(ilev,:) + prsi_obs(ilev+1,:)) * half
  end do

  if (self%troposphere) then
    allocate(troplev_obs(nlocs))
    allocate(airmass_trop(nlocs))
    allocate(airmass_tot(nlocs))
    call obsspace_get_db(obss, "MetaData", "troposphere_layer_index", troplev_obs)
    call obsspace_get_db(obss, "MetaData", "air_mass_factor_troposphere", airmass_trop)
    call obsspace_get_db(obss, "MetaData", "air_mass_factor_total", airmass_tot)
  end if

  ! loop through all variables
  do ivar = 1, nvars
    geovar = self%tracervars(ivar)
    call ufo_geovals_get_var(geovals, geovar, tracer)
    do iobs = 1, nlocs
      if (avgkernel_obs(1,iobs) /= missing) then ! take care of missing obs
        if (self%troposphere) then
          call simulate_column_ob(self%nlayers_kernel, tracer%nval, avgkernel_obs(:,iobs), &
                                  prsl_obs(:,iobs), prsl%vals(:,iobs), temp%vals(:,iobs),&
                                  phii%vals(:,iobs), tracer%vals(:,iobs)*self%convert_factor_model, &
                                  hofx_tmp, troplev_obs(iobs), airmass_tot(iobs), airmass_trop(iobs))
          hofx(ivar,iobs) = hofx_tmp * self%convert_factor_hofx
        end if
      else
        hofx(ivar,iobs) = missing ! default if we are unable to compute averaging kernel
      end if
    end do
  end do


end subroutine ufo_avgkernel_simobs


! ------------------------------------------------------------------------------

end module ufo_avgkernel_mod
