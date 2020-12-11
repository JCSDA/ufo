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
   real(kind_real), allocatable, dimension(:) :: ak_kernel, bk_kernel
   character(kind=c_char,len=:), allocatable :: obskernelvar, tracervars(:)
   logical :: troposphere, totalcolumn
   real(kind_real) :: convert_factor_model, convert_factor_hofx
   real(kind_real), allocatable, dimension(:,:) :: avgkernel_obs, prsl_obs
   real(kind_real), allocatable, dimension(:,:) :: prsl, temp, phii
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
  use ufo_constants_mod, only: half
  use obsspace_mod
  implicit none
  class(ufo_avgkernel_tlad), intent(inout) :: self
  type(ufo_geovals),       intent(in)    :: geovals_in
  type(c_ptr), value,      intent(in)    :: obss

  ! Local variables
  type(ufo_geoval), pointer :: prsi, psfc
  integer :: ivar, iobs, ilev
  character(len=MAXVARLEN) :: varstring
  character(len=4) :: levstr
  real(kind_real), allocatable, dimension(:,:) :: prsi_obs
  type(ufo_geovals) :: geovals
  type(ufo_geoval), pointer :: prsl, temp, phii, p_temp

  ! get nlocs and nvars
  self%nlocs = obsspace_get_nlocs(obss)
  self%nvars = obsspace_get_nvars(obss)

  ! get geovals of atmospheric pressure
  call ufo_geovals_copy(geovals_in, geovals)  ! dont want to change geovals_in
  call ufo_geovals_get_var(geovals, var_prs, p_temp)
  if (p_temp%vals(1,1) < p_temp%vals(p_temp%nval,1) ) then
    self%flip_it = .true.
  else
    self%flip_it = .false.
  end if
  call ufo_geovals_reorderzdir(geovals, var_prs, "bottom2top")
  call ufo_geovals_get_var(geovals, var_ps, psfc)
  call ufo_geovals_get_var(geovals, var_prs, prsl)
  call ufo_geovals_get_var(geovals, var_prsi, prsi)
  call ufo_geovals_get_var(geovals, var_ts, temp)
  call ufo_geovals_get_var(geovals, var_zi, phii)

  ! get nval
  self%nval = prsl%nval

  allocate(self%prsl(self%nval, self%nlocs))
  allocate(self%temp(self%nval, self%nlocs))
  allocate(self%phii(self%nval+1, self%nlocs))
  do iobs = 1, self%nlocs
    self%prsl(:,iobs) = prsl%vals(:,iobs)
    self%temp(:,iobs) = temp%vals(:,iobs)
    self%phii(:,iobs) = phii%vals(:,iobs)
  end do

  ! grab necesary metadata from IODA
  ! get observation averaging kernel
  ! once 2D arrays are allowed, rewrite/simplify this part
  allocate(self%avgkernel_obs(self%nlayers_kernel, self%nlocs))
  do ilev = 1, self%nlayers_kernel
    write(levstr, fmt = "(I3)") ilev
    levstr = adjustl(levstr)
    varstring = trim(self%obskernelvar)//"_"//trim(levstr)
    call obsspace_get_db(obss, "MetaData", trim(varstring), self%avgkernel_obs(ilev, :))
  end do

  ! compute prsl_obs/prsi_obs from ak/bk/psfc
  allocate(self%prsl_obs(self%nlayers_kernel, self%nlocs))
  allocate(prsi_obs(self%nlayers_kernel+1, self%nlocs))
  ! prsi_obs calculation
  do ilev = 1, self%nlayers_kernel+1
    prsi_obs(ilev,:) = self%ak_kernel(ilev) + self%bk_kernel(ilev) * psfc%vals(1,:)
  end do
  ! using simple averaging for now for prsl, can use more complex way later
  do ilev = 1, self%nlayers_kernel
    self%prsl_obs(ilev,:) = (prsi_obs(ilev,:) + prsi_obs(ilev+1,:)) * half
  end do

  if (self%troposphere) then
    allocate(self%troplev_obs(self%nlocs))
    allocate(self%airmass_trop(self%nlocs))
    allocate(self%airmass_tot(self%nlocs))
    call obsspace_get_db(obss, "MetaData", "troposphere_layer_index", self%troplev_obs)
    call obsspace_get_db(obss, "MetaData", "air_mass_factor_troposphere", self%airmass_trop)
    call obsspace_get_db(obss, "MetaData", "air_mass_factor_total", self%airmass_tot)
  end if

  ! cleanup things we do not need now
  deallocate(prsi_obs)

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

  call ufo_geovals_copy(geovals_in, geovals)  ! dont want to change geovals_in

  ! loop through all variables
  do ivar = 1, nvars
    geovar = self%tracervars(ivar)
    call ufo_geovals_get_var(geovals, geovar, tracer)
    do iobs = 1, nlocs
        if(self%flip_it) tracer%vals(1:tracer%nval,iobs) = tracer%vals(tracer%nval:1:-1,iobs)
        if (self%troposphere) then
          call simulate_column_ob_tl(self%nlayers_kernel, tracer%nval, self%avgkernel_obs(:,iobs), &
                                     self%prsl_obs(:,iobs), self%prsl(:,iobs), self%temp(:,iobs),&
                                     self%phii(:,iobs), tracer%vals(:,iobs)*self%convert_factor_model, &
                                     hofx_tmp, self%troplev_obs(iobs), self%airmass_tot(iobs), self%airmass_trop(iobs))
          hofx(ivar,iobs) = hofx_tmp * self%convert_factor_hofx
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

    ! Allocate geovals profile if not yet allocated
    if (.not. allocated(tracer%vals)) then
       tracer%nlocs = self%nlocs
       tracer%nval  = self%nval
       allocate(tracer%vals(tracer%nval, tracer%nlocs))
       tracer%vals(:,:) = 0.0_kind_real
    endif

    do iobs = 1, nlocs
      if (hofx(ivar,iobs) /= missing) then ! take care of missing obs
        if (self%troposphere) then
          hofx_tmp = hofx(ivar,iobs) * self%convert_factor_hofx
          call simulate_column_ob_ad(self%nlayers_kernel, tracer%nval, self%avgkernel_obs(:,iobs), &
                                     self%prsl_obs(:,iobs), self%prsl(:,iobs), self%temp(:,iobs),&
                                     self%phii(:,iobs), tracer%vals(:,iobs), &
                                     hofx_tmp, self%troplev_obs(iobs), self%airmass_tot(iobs), self%airmass_trop(iobs))
          tracer%vals(:,iobs) = tracer%vals(:,iobs) * self%convert_factor_model
        end if
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
  if (allocated(self%ak_kernel)) deallocate(self%ak_kernel)
  if (allocated(self%bk_kernel)) deallocate(self%bk_kernel)
  if (allocated(self%obskernelvar)) deallocate(self%obskernelvar)
  if (allocated(self%tracervars)) deallocate(self%tracervars)
  if (allocated(self%avgkernel_obs)) deallocate(self%avgkernel_obs)
  if (allocated(self%prsl_obs)) deallocate(self%prsl_obs)
  if (allocated(self%prsl)) deallocate(self%prsl)
  if (allocated(self%temp)) deallocate(self%temp)
  if (allocated(self%phii)) deallocate(self%phii)
  if (allocated(self%airmass_tot)) deallocate(self%airmass_tot)
  if (allocated(self%airmass_trop)) deallocate(self%airmass_trop)
  if (allocated(self%troplev_obs)) deallocate(self%troplev_obs)
end subroutine avgkernel_tlad_cleanup_

! ------------------------------------------------------------------------------

end module ufo_avgkernel_tlad_mod
