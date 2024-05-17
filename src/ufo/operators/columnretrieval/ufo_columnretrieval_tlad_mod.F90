! (C) Copyright 2017-2022 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

!> Fortran module for columnretrieval tl/ad observation operator

module ufo_columnretrieval_tlad_mod
  use oops_variables_mod
  use obs_variables_mod
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
   type(obs_variables), public :: obsvars
   type(oops_variables), public :: geovars
   integer :: nlayers_retrieval
   integer :: nlocs, nvars, nval
   character(kind=c_char,len=:), allocatable :: obskernelvar, obspressurevar
   character(kind=c_char,len=:), allocatable :: tracervars, stretch
   logical :: isaveragingkernel, totalnovertice
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
  character(len=:), allocatable :: value(:)

  ! get configuration for the averaging kernel operator
  call f_conf%get_or_die("nlayers_retrieval", self%nlayers_retrieval)

  ! get name of geoval/tracer to use from the model
  call f_conf%get_or_die("tracer variables", self%tracervars)

  ! determine if the averaging kernel is needed
  call f_conf%get_or_die("isAveragingKernel", self%isaveragingkernel)

  ! determine if the averaging kernel is needed
  call f_conf%get_or_die("stretchVertices", self%stretch)

  ! do we need a conversion factor, say between ppmv and unity?
  call f_conf%get_or_die("model units coeff", self%convert_factor_model)

  ! perform a simple total column calculation using model whole profile
  call f_conf%get_or_die("totalNoVertice", self%totalnovertice)

  ! add variables to geovars that are needed
  ! specified tracers
  call self%geovars%push_back(self%tracervars)

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
  integer :: iobs, ilev
  character(len=MAXVARLEN) :: varstring
  character(len=4) :: levstr
  type(ufo_geovals) :: geovals
  character(len=max_string) :: err_msg

  !prevent impossible options
  if (self%totalnovertice .and. (self%nlayers_retrieval > 1 .or. &
                                   self%isaveragingkernel )) then
    write(err_msg, *) "Error: wrong combination of yaml options, &
                  & totalNoVertice, isAveragingKerne, nlayers_retrieval"
    call abor1_ftn(err_msg)
  end if

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

  allocate(self%avgkernel_obs(self%nlayers_retrieval, self%nlocs))
  ! set to 1.0 if no ak provided
  self%avgkernel_obs = one
  if (self%isaveragingkernel) then
    call obsspace_get_db_2d(obss, "RetrievalAncillaryData", &
                            "averagingKernel", self%avgkernel_obs)
  end if

  ! get prsi_obs
  allocate(self%prsi_obs(self%nlayers_retrieval+1, self%nlocs))
  if (.not. self%totalnovertice) then
    call obsspace_get_db_2d(obss, "RetrievalAncillaryData", &
                          "pressureVertice", self%prsi_obs)
  else
    self%prsi_obs(self%nlayers_retrieval+1, :) = self%prsi(self%nval, :)
    self%prsi_obs(1, :) = self%prsi(1,:)
  end if

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

  ! nvars is always 1 for this operator
  call ufo_geovals_get_var(geovals, self%tracervars, tracer)
  do iobs = 1, nlocs
    if (self%avgkernel_obs(1,iobs) /= missing) then ! take care of missing obs
      call simulate_column_ob_tl(self%nlayers_retrieval, tracer%nval, &
                                 self%avgkernel_obs(:,iobs), &
                                 self%prsi_obs(:,iobs), self%prsi(:,iobs),&
                                 tracer%vals(:,iobs)*self%convert_factor_model, &
                                 hofx_tmp, self%stretch)
      hofx(nvars,iobs) = hofx_tmp
    else
      hofx(nvars,iobs) = missing
    end if
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
  type(ufo_geovals) :: geovals
  real(kind_real) :: hofx_tmp
  integer :: iobs
  real(c_double) :: missing

  missing = missing_value(missing)

  if (.not. geovals_in%linit ) geovals_in%linit=.true. ! need this for var exe

  ! nvars is always 1 for this operator
  call ufo_geovals_get_var(geovals_in, self%tracervars, tracer)

  do iobs = 1, nlocs
    if (hofx(nvars,iobs) /= missing) then ! take care of missing obs
      hofx_tmp = hofx(nvars,iobs)
      call simulate_column_ob_ad(self%nlayers_retrieval, tracer%nval, &
                                 self%avgkernel_obs(:,iobs), &
                                 self%prsi_obs(:,iobs), self%prsi(:,iobs),&
                                 tracer%vals(:,iobs), hofx_tmp, self%stretch)
      tracer%vals(:,iobs) = tracer%vals(:,iobs) * self%convert_factor_model
    end if
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
