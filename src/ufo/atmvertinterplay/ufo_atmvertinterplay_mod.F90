! (C) Copyright 2017-2018 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

!> Fortran module for atmvertinterplay observation operator

module ufo_atmvertinterplay_mod

 use oops_variables_mod
 use ufo_vars_mod
 use vert_interp_lay_mod
 implicit none
 private

!> Fortran derived type for the observation type
 type, public :: ufo_atmvertinterplay
   type(oops_variables), public :: obsvars
   type(oops_variables), public :: geovars
   integer, public, allocatable :: nlevels(:)
   real, public, allocatable :: coefficients(:) ! unit conversion from geoval to obs
 contains
   procedure :: setup  => ufo_atmvertinterplay_setup
   procedure :: simobs => ufo_atmvertinterplay_simobs
 end type ufo_atmvertinterplay

contains

! ------------------------------------------------------------------------------
subroutine ufo_atmvertinterplay_setup(self, conf)
use iso_c_binding
use fckit_configuration_module, only: fckit_configuration
implicit none
class(ufo_atmvertinterplay), intent(inout) :: self
type(fckit_configuration), intent(in) :: conf

character(kind=c_char,len=:), allocatable :: coord_name
character(kind=c_char,len=:), allocatable :: gvars(:)
real(kind=c_double), allocatable :: coefficients(:)
integer(kind=c_int), allocatable :: nlevels(:)
!Local Variables
integer :: ivar, nlevs=0, nvars=0, ngvars=0, ncoefs=0

! Check configurations
if (conf%has("geovals")) then
  ngvars = conf%get_size("geovals")
  call conf%get_or_die("geovals", gvars)
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

if (conf%has("coefficients")) then
  ncoefs = conf%get_size("coefficients")
  call conf%get_or_die("coefficients", coefficients)
  allocate(self%coefficients(ncoefs))
  self%coefficients(1:ncoefs) = coefficients(1:ncoefs)
endif

if (conf%has("nlevels")) then
   nlevs = conf%get_size("nlevels")
   call conf%get_or_die("nlevels", nlevels)
   allocate(self%nlevels(nlevs))
   self%nlevels(1:nlevs) = nlevels(1:nlevs)
endif

! Put pressure to the geovars (vars from the model) list
call self%geovars%push_back(var_prsi)

end subroutine ufo_atmvertinterplay_setup


! ------------------------------------------------------------------------------
subroutine ufo_atmvertinterplay_simobs(self, geovals_in, obss, nvars, nlocs, hofx)
use ufo_geovals_mod, only: ufo_geovals, ufo_geoval, ufo_geovals_get_var
use ufo_geovals_mod, only: ufo_geovals_delete, ufo_geovals_copy, ufo_geovals_reorderzdir
use ufo_constants_mod
use obsspace_mod
implicit none
class(ufo_atmvertinterplay), intent(in)    :: self
integer, intent(in)               :: nvars, nlocs
type(ufo_geovals),  intent(in)    :: geovals_in
real(c_double),     intent(inout) :: hofx(nvars, nlocs)
type(c_ptr), value, intent(in)    :: obss

! Local variables
integer :: iobs, ivar 
integer :: nsig, nlevs
real(kind_real), dimension(:), allocatable :: toppressure,botpressure,airpressure
type(ufo_geovals) :: geovals
type(ufo_geoval), pointer :: modelpressures, modelozone
character(len=MAXVARLEN) :: geovar
character(len=MAXVARLEN) :: var_zdir
real :: pob,delp4,delz,dz1
real(kind_real) :: layer_oz

  ! Notes:
  ! (1) Set desired vertical coordinate direction (top2bottom or bottom2top) based
  !     on vertical coodinate variable and reload geoavls according to the set
  !     direction
  ! (2) This is done because this observation operator assums pressure levels
  !     are from bottom to top (bottom2top) with morelpressure(1) for surface and
  !     modelpressure(nsig+1) for model top

  call ufo_geovals_copy(geovals_in, geovals)  ! dont want to change geovals_in
  var_zdir = var_prsi                         ! vertical coordinate variable
  call ufo_geovals_reorderzdir(geovals, var_zdir, "bottom2top")

  ! Get pressure profiles from geovals [Pa]
  call ufo_geovals_get_var(geovals, var_prsi, modelpressures)
  nsig = modelpressures%nval - 1
  ! Allocate pressure limits and get air_pressure metadata from obs
  allocate(toppressure(nlocs))
  allocate(botpressure(nlocs))
  allocate(airpressure(nlocs))
  call obsspace_get_db(obss, "MetaData", "air_pressure", airpressure)  

  do ivar = 1, nvars
    write(6,*) 'ufo_atmvertinterplay_simobs: self%nlevels = ', self%nlevels
    nlevs = self%nlevels(ivar) 
    call get_integral_limits(airpressure, botpressure, toppressure, modelpressures%vals(:,:), nlevs, nlocs, nsig) 

    !Get the name of input variable in geovals
    geovar = self%geovars%variable(ivar)

    !Get model output
    call ufo_geovals_get_var(geovals, geovar, modelozone)

    do iobs = 1, nlocs
      call apply_layer_integral(self%coefficients(ivar), modelozone%vals(:,iobs), modelpressures%vals(:,iobs), botpressure(iobs), toppressure(iobs), nsig, layer_oz )
      hofx(ivar,iobs) = layer_oz
    enddo
  enddo
  deallocate(toppressure)
  deallocate(botpressure)
  call ufo_geovals_delete(geovals)
end subroutine ufo_atmvertinterplay_simobs


! ------------------------------------------------------------------------------

end module ufo_atmvertinterplay_mod
