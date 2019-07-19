! (C) Copyright 2017-2018 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

!> Fortran module for atmvertinterplay observation operator

module ufo_atmvertinterplay_mod

 use iso_c_binding
 use config_mod
 use kinds

 use ufo_geovals_mod, only: ufo_geovals, ufo_geoval, ufo_geovals_get_var
 use ufo_geovals_mod_c, only: ufo_geovals_registry
 use ufo_basis_mod, only: ufo_basis
 use ufo_vars_mod
 use ufo_constants_mod
 use obsspace_mod

 implicit none
 private
 integer, parameter :: max_string=800

!> Fortran derived type for the observation type
 type, public :: ufo_atmvertinterplay
 private
   integer :: nvars  ! number of variables to be interpolated
   character(len=max_string), public, allocatable :: varin(:)
   character(len=max_string), public, allocatable :: varout(:)
 contains
   procedure :: setup  => ufo_atmvertinterplay_setup
   procedure :: simobs => ufo_atmvertinterplay_simobs
   final :: destructor
 end type ufo_atmvertinterplay

contains

! ------------------------------------------------------------------------------
subroutine ufo_atmvertinterplay_setup(self, vars)
implicit none
class(ufo_atmvertinterplay), intent(inout) :: self
character(len=MAXVARLEN), dimension(:), intent(inout) :: vars
  
!Local Variables
integer :: i

  self%nvars = size(vars)
  allocate(self%varout(self%nvars))
  self%varout = vars

  ! Allocate varin: variables we need from the model
  !  need additional slot to hold vertical coord.
  allocate(self%varin(self%nvars+1))

  ! Set vars_in based on vars_out
  do i = 1, self%nvars
    self%varin(i) = self%varout(i)
  enddo

  ! Put pressure to the varin (vars from the model) list
  self%varin(self%nvars+1) = "air_pressure_levels"

end subroutine ufo_atmvertinterplay_setup

! ------------------------------------------------------------------------------
subroutine destructor(self)
implicit none
type(ufo_atmvertinterplay), intent(inout) :: self

  if (allocated(self%varout)) deallocate(self%varout)
  if (allocated(self%varin))  deallocate(self%varin)

end subroutine destructor

! ------------------------------------------------------------------------------
subroutine ufo_atmvertinterplay_simobs(self, geovals, obss, nvars, nlocs, hofx)
implicit none
class(ufo_atmvertinterplay), intent(in)    :: self
integer, intent(in)               :: nvars, nlocs
type(ufo_geovals),  intent(in)    :: geovals
real(c_double),     intent(inout) :: hofx(nvars, nlocs)
type(c_ptr), value, intent(in)    :: obss

! Local variables
type(ufo_geoval), pointer :: geoval
integer :: i,iobs, ivar
integer :: iz1,iz2,k,kk
integer :: nsig
real(kind_real), dimension(:), allocatable :: obss_metadata
real(kind_real), dimension(:), allocatable :: toppressure,botpressure
real(kind_real), dimension(:), allocatable :: lat, lon
type(ufo_geoval), pointer :: modelpressures, modelozone
character(len=MAXVARLEN) :: geovar
real :: pob,delp4,delz,dz1
real(kind_real) :: rozcon,g
real(kind_real) :: topozp,botozp
real :: pindex

  rozcon = 1./(1.e-3*DU*grav)

  ! Get pressure profiles from geovals log(cb)
  call ufo_geovals_get_var(geovals, var_prsi, modelpressures)
  nsig = modelpressures%nval - 1
  !geoval pressures read in as Pa

  allocate(toppressure(nlocs))
  allocate(botpressure(nlocs))
  allocate(lat(nlocs))
  allocate(lon(nlocs))

  !obs pressures read in as Pa
  call obsspace_get_db(obss, "MetaData", "top_level_pressure", toppressure)
  call obsspace_get_db(obss, "MetaData", "bottom_level_pressure", botpressure)
  call obsspace_get_db(obss, "MetaData", "latitude", lat)
  call obsspace_get_db(obss, "MetaData", "longitude", lon)

  do ivar = 1, self%nvars
    !get the name of input variable in geovals
    geovar = self%varin(ivar)

    !Get model output
    call ufo_geovals_get_var(geovals, geovar, modelozone)

     do iobs = 1, nlocs
      topozp = pindex(nsig, modelpressures%vals(1, iobs), toppressure(iobs))
      botozp = pindex(nsig, modelpressures%vals(1, iobs), botpressure(iobs))
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
        delp4=(modelpressures%vals(kk,iobs)-modelpressures%vals(kk+1,iobs))*.001
        g=g + &
             modelozone%vals(kk,iobs)*rozcon*delz*delp4
   print *,iobs,k,g,modelpressures%vals(kk,iobs),modelpressures%vals(kk+1,iobs)
   print *,lat(iobs),lon(iobs)
      enddo
      hofx(ivar,iobs) = g
      dz1 = pob
     enddo
     enddo
     deallocate(toppressure)
     deallocate(botpressure)

end subroutine ufo_atmvertinterplay_simobs


! ------------------------------------------------------------------------------

end module ufo_atmvertinterplay_mod
