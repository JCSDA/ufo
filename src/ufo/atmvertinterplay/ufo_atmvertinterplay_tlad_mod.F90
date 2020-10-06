! (C) Copyright 2017-2019 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

module ufo_atmvertinterplay_tlad_mod

  use oops_variables_mod
  use ufo_vars_mod
  use ufo_geovals_mod 
  use vert_interp_lay_mod
  use missing_values_mod
  use, intrinsic :: iso_c_binding
  use kinds, only: kind_real
  use missing_values_mod

! ------------------------------------------------------------------------------

  type, public :: ufo_atmvertinterplay_tlad
  private
    type(oops_variables), public :: obsvars
    type(oops_variables), public :: geovars
    character(len=MAXVARLEN), public :: v_coord ! GeoVaL to use to interpolate in vertical)
    logical, public :: use_ln
    integer :: nval, nlocs
    real(kind_real), dimension(:), allocatable :: toppressure, botpressure, airpressure
    integer, allocatable :: nlevels(:)
    real, allocatable :: coefficients(:) ! unit conversion from geoval to obs
    real(kind_real),allocatable :: modelpressures(:,:)
    real(kind_real),allocatable :: modeloz(:,:)
    !real(kind_real),allocatable :: profiles(:,:)
  contains
    procedure :: setup => atmvertinterplay_tlad_setup_
    procedure :: cleanup => atmvertinterplay_tlad_cleanup_
    procedure :: settraj => atmvertinterplay_tlad_settraj_
    procedure :: simobs_tl => atmvertinterplay_simobs_tl_
    procedure :: simobs_ad => atmvertinterplay_simobs_ad_
    final :: destructor
  end type ufo_atmvertinterplay_tlad

! ------------------------------------------------------------------------------
contains
! ------------------------------------------------------------------------------
subroutine atmvertinterplay_tlad_setup_(self, grid_conf)
  use fckit_configuration_module, only: fckit_configuration
  implicit none
  class(ufo_atmvertinterplay_tlad), intent(inout) :: self
  type(fckit_configuration), intent(in) :: grid_conf

  character(kind=c_char,len=:), allocatable :: coord_name
  character(kind=c_char,len=:), allocatable :: gvars(:)
  real(kind=c_double), allocatable :: coefficients(:)
  integer(kind=c_int), allocatable :: nlevels(:)
  !Local Variables
  integer :: ivar, nlevs=0, nvars=0, ngvars=0, ncoefs=0
  ! Check configurations
  if (grid_conf%has("geovals")) then
    ngvars = grid_conf%get_size("geovals")
    call grid_conf%get_or_die("geovals", gvars)
    ! add to geovars list
    do ivar = 1, ngvars
      call self%geovars%push_back(gvars(ivar))
    enddo
  endif
  nvars = self%obsvars%nvars()
  if (ngvars == 0 .and. nvars > 0) then
    allocate(self%coefficients(nvars))
    do ivar = 1, nvars
      call self%geovars%push_back(self%obsvars%variable(ivar))
      self%coefficients(ivar) = 1.0
    enddo
  endif
  if (grid_conf%has("coefficients")) then
    ncoefs = grid_conf%get_size("coefficients")
    call grid_conf%get_or_die("coefficients", coefficients)
    allocate(self%coefficients(ncoefs))
    self%coefficients(1:ncoefs) = coefficients(1:ncoefs)
  endif
  if (grid_conf%has("nlevels")) then
    nlevs = grid_conf%get_size("nlevels")
    call grid_conf%get_or_die("nlevels", nlevels)
    allocate(self%nlevels(nlevs))
    self%nlevels(1:nlevs) = nlevels(1:nlevs)
  endif

  ! Put pressure to the geovars (vars from the model) list
  call self%geovars%push_back(var_prsi)

end subroutine atmvertinterplay_tlad_setup_

! ------------------------------------------------------------------------------

subroutine atmvertinterplay_tlad_settraj_(self, geovals, obss)
  use obsspace_mod
  implicit none
  class(ufo_atmvertinterplay_tlad), intent(inout) :: self
  type(ufo_geovals),         intent(in)    :: geovals
  type(c_ptr), value,        intent(in)    :: obss
  integer :: iobs,nlevs,nsig,ilev
  !type(ufo_geovals) :: geovals
  type(ufo_geoval),pointer :: modelpres
  type(ufo_geoval),pointer :: profile
  character(len=MAXVARLEN) :: var_zdir
  character(len=MAXVARLEN) :: geovar
  real(kind_real), dimension(:), allocatable :: airpressure
  ! Make sure nothing already allocated
  call self%cleanup()

  ! Get the observation vertical coordinates
  self%nlocs = obsspace_get_nlocs(obss)
  ! Allocate arrays for top and bottom pressures for integral
  allocate(self%toppressure(self%nlocs))
  allocate(self%botpressure(self%nlocs))
  allocate(airpressure(self%nlocs)) 
  print *, 'in traj'
  !call ufo_geovals_copy(geovals_in, geovals)  ! dont want to change geovals_in
  var_zdir = var_prsi                         ! vertical coordinate variable
  !call ufo_geovals_reorderzdir(geovals, var_zdir, "bottom2top")

  ! Get pressure profiles from geovals [Pa]
  call ufo_geovals_get_var(geovals, var_prsi, modelpres)
  nlevs = self%nlevels(1)
  nsig = modelpres%nval - 1
  print *, 'nsig traj', nsig
  call obsspace_get_db(obss, "MetaData", "air_pressure", airpressure)  

  allocate(self%modelpressures(modelpres%nval,self%nlocs))
  allocate(self%modeloz(modelpres%nval,self%nlocs))
  allocate(self%airpressure(self%nlocs))
  self%modelpressures(1:nsig+1,1:self%nlocs) = modelpres%vals(1:nsig+1,1:self%nlocs)
  self%airpressure(1:self%nlocs) = airpressure(1:self%nlocs)
  call get_integral_limits(airpressure(:), self%botpressure(:), self%toppressure(:), self%modelpressures(:,:), nlevs, self%nlocs, nsig) 

  geovar = self%geovars%variable(1)
  ! interplayolate from geovals to observational location into hofx
  var_zdir = var_prsi                         ! vertical coordinate variable
  call ufo_geovals_get_var(geovals, geovar, profile)
  self%modeloz(1:nsig+1,1:self%nlocs) = profile%vals(1:nsig+1,1:self%nlocs)

  ! Cleanup memory
  deallocate(airpressure)
  !call ufo_geovals_delete(geovals)
  print *, 'done traj'
end subroutine atmvertinterplay_tlad_settraj_

! ------------------------------------------------------------------------------

subroutine atmvertinterplay_simobs_tl_(self, geovals, obss, nvars, nlocs, hofx)
  implicit none
  class(ufo_atmvertinterplay_tlad), intent(in) :: self
  !type(ufo_geovals),         intent(in) :: geovals_in
  integer,                   intent(in) :: nvars, nlocs
  real(c_double),         intent(inout) :: hofx(nvars, nlocs)
  type(c_ptr), value,        intent(in) :: obss

  integer :: iobs, ivar
  type(ufo_geoval), pointer :: profile
  type(ufo_geoval), pointer :: pressure
  character(len=MAXVARLEN) :: geovar
  character(len=MAXVARLEN) :: var_zdir
  type(ufo_geovals) :: geovals
  integer :: nsig
  !call ufo_geovals_copy(geovals_in, geovals)  ! dont want to change geovals_in
  do ivar = 1, nvars
    ! Get the name of input variable in geovals
    geovar = self%geovars%variable(ivar)
    ! interplayolate from geovals to observational location into hofx
    var_zdir = var_prsi                         ! vertical coordinate variable
    !call ufo_geovals_reorderzdir(geovals, var_zdir, "bottom2top")
    call ufo_geovals_get_var(geovals, geovar, profile)
    nsig = profile%nval
    do iobs = 1, nlocs
       call vert_interp_lay_apply_tl(profile%vals(:,iobs), hofx(ivar,iobs), self%coefficients(ivar),  self%modelpressures(:,iobs), self%botpressure(iobs), self%toppressure(iobs), nsig)
       if(iobs>1 .and. iobs < 23) then 
         print *,'iobs,h*(x)',iobs,hofx(ivar,iobs)
         !print *, 'modelpres h*(x)', pressure%vals(1,iobs)
         !print *, 'end modelpres'
       endif
       !if(bp(iobs) - tp(iobs)< 0.0) print *,'pressurebroken',iobs
    enddo
  enddo
end subroutine atmvertinterplay_simobs_tl_

! ------------------------------------------------------------------------------

subroutine atmvertinterplay_simobs_ad_(self, geovals, obss, nvars, nlocs, hofx)
  implicit none
  class(ufo_atmvertinterplay_tlad), intent(in) :: self
  type(ufo_geovals),         intent(inout) :: geovals
  integer,                   intent(in)    :: nvars, nlocs
  real(c_double),            intent(in)    :: hofx(nvars, nlocs)
  type(c_ptr), value,        intent(in)    :: obss

  integer :: iobs, ivar
  type(ufo_geoval), pointer :: profile
  character(len=MAXVARLEN) :: geovar
  real(c_double) :: missing

  missing = missing_value(missing)

  do ivar = 1, nvars
    ! Get the name of input variable in geovals
    geovar = self%geovars%variable(ivar)

    ! Get pointer to profile for this variable in geovals
    call ufo_geovals_get_var(geovals, geovar, profile)

    ! Allocate geovals profile if not yet allocated
    if (.not. allocated(profile%vals)) then
       profile%nlocs = self%nlocs
       profile%nval  = self%nval
       allocate(profile%vals(profile%nval, profile%nlocs))
       profile%vals(:,:) = 0.0_kind_real
    endif
    if (.not. geovals%linit ) geovals%linit=.true.

    ! Adjoint of interplayolate, from hofx into geovals
    do iobs = 1, self%nlocs
      if (hofx(ivar,iobs) /= missing) then
        
        !call vert_interplay_apply_ad(profile%nval, profile%vals(:,iobs), &
        !                        & hofx(ivar,iobs), self%wi(iobs), self%wf(iobs))
      endif
    enddo
  enddo
end subroutine atmvertinterplay_simobs_ad_

! ------------------------------------------------------------------------------

subroutine atmvertinterplay_tlad_cleanup_(self)
  implicit none
  class(ufo_atmvertinterplay_tlad), intent(inout) :: self
  self%nval = 0
  self%nlocs = 0
  if (allocated(self%toppressure)) deallocate(self%toppressure)
  if (allocated(self%botpressure)) deallocate(self%botpressure)
  if (allocated(self%modelpressures)) deallocate(self%modelpressures)
end subroutine atmvertinterplay_tlad_cleanup_

! ------------------------------------------------------------------------------

subroutine  destructor(self)
  type(ufo_atmvertinterplay_tlad), intent(inout)  :: self

  call self%cleanup()

end subroutine destructor

! ------------------------------------------------------------------------------

end module ufo_atmvertinterplay_tlad_mod
