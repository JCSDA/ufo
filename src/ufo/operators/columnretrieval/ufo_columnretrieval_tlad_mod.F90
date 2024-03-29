! (C) Copyright 2017-2022 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

!> Fortran module for columnretrieval tl/ad observation operator

module ufo_columnretrieval_tlad_mod
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

 type, public :: ufo_columnretrieval_tlad
 private
   type(oops_variables), public :: obsvars
   type(oops_variables), public :: geovars
   integer :: nlayers_retrieval
   integer :: nlocs, nvars, nval
   character(kind=c_char,len=:), allocatable :: obskernelvar, obspressurevar
   character(kind=c_char,len=:), allocatable :: tracervars(:), stretch
   logical :: isaveragingkernel
   real(kind_real) :: convert_factor_model
   real(kind_real), allocatable, dimension(:,:) :: avgkernel_obs, prsi_obs
   real(kind_real), allocatable, dimension(:,:) :: prsi
 contains
   procedure :: setup  => columnretrieval_tlad_setup_
   procedure :: cleanup  => columnretrieval_tlad_cleanup_
   procedure :: settraj => columnretrieval_tlad_settraj_
   procedure :: simobs_tl  => columnretrieval_simobs_tl_
   procedure :: simobs_ad  => columnretrieval_simobs_ad_
   final :: destructor
 end type ufo_columnretrieval_tlad

contains

! ------------------------------------------------------------------------------
subroutine columnretrieval_tlad_setup_(self, f_conf)
  use fckit_configuration_module, only: fckit_configuration
  use ufo_constants_mod, only: one
  implicit none
  class(ufo_columnretrieval_tlad), intent(inout) :: self
  type(fckit_configuration), intent(in)  :: f_conf
  integer :: nlevs_yaml
  integer :: ivar, nvars
  character(len=max_string) :: err_msg

  ! get configuration for the averaging kernel operator
  call f_conf%get_or_die("nlayers_retrieval", self%nlayers_retrieval)

  ! get variable name from IODA for observation averaging kernel
  call f_conf%get_or_die("AvgKernelVar", self%obskernelvar)

  ! get vertical pressure grid 
  call f_conf%get_or_die("PresLevVar", self%obspressurevar)
 
  ! get name of geoval/tracer to use from the model
  call f_conf%get_or_die("tracer variables", self%tracervars)

  ! determine if the averaging kernel is needed
  call f_conf%get_or_die("isAveragingKernel", self%isaveragingkernel)

  ! determine if the averaging kernel is needed
  call f_conf%get_or_die("stretchVertices", self%stretch)

  ! do we need a conversion factor, say between ppmv and unity?
  call f_conf%get_or_die("model units coeff", self%convert_factor_model)

  ! add variables to geovars that are needed
  ! specified tracers
  nvars = self%obsvars%nvars()
  do ivar = 1, nvars
    call self%geovars%push_back(self%tracervars(ivar))
  end do

end subroutine columnretrieval_tlad_setup_

! ------------------------------------------------------------------------------
subroutine destructor(self)
  implicit none
  type(ufo_columnretrieval_tlad), intent(inout) :: self

  call self%cleanup()

end subroutine destructor

! ------------------------------------------------------------------------------
subroutine columnretrieval_tlad_settraj_(self, geovals_in, obss)
  use iso_c_binding
  use ufo_geovals_mod, only: ufo_geovals, ufo_geoval, ufo_geovals_get_var
  use ufo_constants_mod, only: one, zero
  use obsspace_mod
  implicit none
  class(ufo_columnretrieval_tlad), intent(inout) :: self
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

  allocate(self%avgkernel_obs(self%nlayers_retrieval, self%nlocs))
  ! set to 1.0 if no ak provided
  self%avgkernel_obs = one
  if (self%isaveragingkernel) then
    do ilev = self%nlayers_retrieval, 1, -1
      write(levstr, fmt = "(I3)") ilev
      levstr = adjustl(levstr)
      varstring = trim(self%obskernelvar)//"_"//trim(levstr)
      call obsspace_get_db(obss, "RtrvlAncData", trim(varstring), &
            self%avgkernel_obs(self%nlayers_retrieval+1-ilev, :))
    end do
  end if

  ! get prsi_obs
  allocate(self%prsi_obs(self%nlayers_retrieval+1, self%nlocs))
  do ilev = self%nlayers_retrieval+1, 1, -1
    write(levstr, fmt = "(I3)") ilev
    levstr = adjustl(levstr)
    varstring = trim(self%obspressurevar)//"_"//trim(levstr)
    call obsspace_get_db(obss, "RtrvlAncData", trim(varstring), &
            self%prsi_obs(self%nlayers_retrieval+2-ilev, :))
  end do

end subroutine columnretrieval_tlad_settraj_

! ------------------------------------------------------------------------------
subroutine columnretrieval_simobs_tl_(self, geovals_in, obss, nvars, nlocs, hofx)
  use iso_c_binding
  use ufo_geovals_mod, only: ufo_geovals, ufo_geoval, ufo_geovals_get_var
  use obsspace_mod
  use satcolumn_mod, only: simulate_column_ob_tl
  implicit none
  class(ufo_columnretrieval_tlad), intent(in)    :: self
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
        call simulate_column_ob_tl(self%nlayers_retrieval, tracer%nval, & 
                                   self%avgkernel_obs(:,iobs), &
                                   self%prsi_obs(:,iobs), self%prsi(:,iobs),&
                                   tracer%vals(:,iobs)*self%convert_factor_model, &
                                   hofx_tmp, self%stretch)
        hofx(ivar,iobs) = hofx_tmp
      else
        hofx(ivar,iobs) = missing
      end if
    end do
  end do
  call ufo_geovals_delete(geovals)

end subroutine columnretrieval_simobs_tl_

! ------------------------------------------------------------------------------
subroutine columnretrieval_simobs_ad_(self, geovals_in, obss, nvars, nlocs, hofx)
  use iso_c_binding
  use ufo_geovals_mod, only: ufo_geovals, ufo_geoval, ufo_geovals_get_var
  use satcolumn_mod, only: simulate_column_ob_ad
  use obsspace_mod
  implicit none
  class(ufo_columnretrieval_tlad), intent(in)    :: self
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
        call simulate_column_ob_ad(self%nlayers_retrieval, tracer%nval, &
                                   self%avgkernel_obs(:,iobs), &
                                   self%prsi_obs(:,iobs), self%prsi(:,iobs),&
                                   tracer%vals(:,iobs), hofx_tmp, self%stretch)
        tracer%vals(:,iobs) = tracer%vals(:,iobs) * self%convert_factor_model
      end if
    end do
  end do

end subroutine columnretrieval_simobs_ad_

! ------------------------------------------------------------------------------
subroutine columnretrieval_tlad_cleanup_(self)
  implicit none
  class(ufo_columnretrieval_tlad), intent(inout) :: self
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
  if (allocated(self%stretch)) deallocate(self%stretch)
end subroutine columnretrieval_tlad_cleanup_

! ------------------------------------------------------------------------------

end module ufo_columnretrieval_tlad_mod
