! (C) Copyright 2017-2018 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

!> Fortran module for atmvertinterplay observation operator

module ufo_atmvertinterplay_mod

 use oops_variables_mod
 use ufo_vars_mod

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
integer :: iobs, ivar, iprof
integer :: iz1, iz2, kk
integer :: k1, k2
integer :: nsig, nprof, nlev
real(kind_real), dimension(:), allocatable :: toppressure,botpressure,airpressure
type(ufo_geovals) :: geovals
type(ufo_geoval), pointer :: modelpressures, modelozone
character(len=MAXVARLEN) :: geovar
character(len=MAXVARLEN) :: var_zdir
real :: pob,delp4,delz,dz1
real(kind_real) :: rozcon, g
real(kind_real) :: topozp, botozp
real :: pindex

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

  allocate(toppressure(nlocs))
  allocate(botpressure(nlocs))
  allocate(airpressure(nlocs))

  do ivar = 1, nvars
    write(6,*) 'ufo_atmvertinterplay_simobs: self%nlevels = ', self%nlevels
    if (self%nlevels(ivar) == 1) then ! total column ozone
       do iobs = 1, nlocs
          toppressure(iobs) = modelpressures%vals(nsig+1, iobs)
          botpressure(iobs) = modelpressures%vals(1, iobs)
       enddo
    else
      !Obs pressures read in as Pa
      call obsspace_get_db(obss, "MetaData", "air_pressure", airpressure)
      nlev = self%nlevels(ivar)
      nprof = nlocs/nlev
      iobs = 0
      do iprof = 1, nprof
        do kk = 1, nlev
          k1 = kk
          k2 = kk - 1
          if (k2 == 0) k2 = 1
          if (kk == nlev) then
            k1 = nlev - 1
            k2 = 1
          endif
          iobs = iobs+1
          toppressure(iobs) = airpressure(k2)
          botpressure(iobs) = airpressure(k1)
          if( kk== 1 ) then
             toppressure(iobs) =modelpressures%vals(nsig+1, iobs)
             botpressure(iobs) = airpressure(1)
          else if( kk == nlev) then
             toppressure(iobs) = modelpressures%vals(nsig+1, iobs)
             botpressure(iobs) = modelpressures%vals(1, iobs)
          endif
        enddo
      enddo
    endif

    !Get the name of input variable in geovals
    geovar = self%geovars%variable(ivar)

    !Get model output
    call ufo_geovals_get_var(geovals, geovar, modelozone)

    do iobs = 1, nlocs
      topozp = pindex(nsig+1, modelpressures%vals(1, iobs), toppressure(iobs))
      botozp = pindex(nsig+1, modelpressures%vals(1, iobs), botpressure(iobs))

      pob = botozp
      iz1 = topozp
      if (iz1>nsig) iz1=nsig
      iz2 = pob
      !For total column ozone
      if(iz1 .eq. nsig .and. iz2 .lt.7)iz2 = 1
      g = 0.
      dz1 = topozp
      do kk=iz1,iz2,-1
        delz = 1.
        if(kk==iz1)delz=dz1-iz1
        if (kk==iz2) delz=delz-pob+iz2
        !For total column ozone
        if(iz1 .eq. nsig .and. iz2 .eq. 1)delz = 1
        !Interpolate in cbars
        delp4 = (modelpressures%vals(kk,iobs)-modelpressures%vals(kk+1,iobs))  ! [Pa]
        g = g + modelozone%vals(kk,iobs)*self%coefficients(ivar)*(delz*delp4)
      enddo
      hofx(ivar,iobs) = g
      dz1 = pob
    enddo
  enddo
  deallocate(toppressure)
  deallocate(botpressure)
  call ufo_geovals_delete(geovals)

end subroutine ufo_atmvertinterplay_simobs


! ------------------------------------------------------------------------------

end module ufo_atmvertinterplay_mod
