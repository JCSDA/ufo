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
 contains
   procedure :: setup  => ufo_atmvertinterplay_setup
   procedure :: simobs => ufo_atmvertinterplay_simobs
 end type ufo_atmvertinterplay

contains

! ------------------------------------------------------------------------------
subroutine ufo_atmvertinterplay_setup(self)
implicit none
class(ufo_atmvertinterplay), intent(inout) :: self

!Local Variables
integer :: ivar, nvars

! Fill in geovars: variables we need from the model
!  need additional slot to hold vertical coord.
nvars = self%obsvars%nvars()
do ivar = 1, nvars
  call self%geovars%push_back(self%obsvars%variable(ivar))
enddo

! Put pressure to the geovars (vars from the model) list
call self%geovars%push_back(var_prsi)

end subroutine ufo_atmvertinterplay_setup

! ------------------------------------------------------------------------------
subroutine ufo_atmvertinterplay_simobs(self, geovals, obss, nvars, nlocs, hofx)
use ufo_geovals_mod, only: ufo_geovals, ufo_geoval, ufo_geovals_get_var
use ufo_constants_mod
use obsspace_mod
implicit none
class(ufo_atmvertinterplay), intent(in)    :: self
integer, intent(in)               :: nvars, nlocs
type(ufo_geovals),  intent(in)    :: geovals
real(c_double),     intent(inout) :: hofx(nvars, nlocs)
type(c_ptr), value, intent(in)    :: obss

! Local variables
integer :: iobs, ivar
integer :: iz1,iz2,kk
integer :: nsig
real(kind_real), dimension(:), allocatable :: toppressure,botpressure
type(ufo_geoval), pointer :: modelpressures, modelozone
character(len=MAXVARLEN) :: geovar
real :: pob,delp4,delz,dz1
real(kind_real) :: rozcon,g
real(kind_real) :: topozp,botozp
real :: pindex

  rozcon = 1.e3/(DU*grav)

  ! Get pressure profiles from geovals log(cb)
  call ufo_geovals_get_var(geovals, var_prsi, modelpressures)
  nsig = modelpressures%nval - 1
  !geoval pressures read in as Pa

  allocate(toppressure(nlocs))
  allocate(botpressure(nlocs))

  !obs pressures read in as Pa
  call obsspace_get_db(obss, "MetaData", "top_level_pressure", toppressure)
  call obsspace_get_db(obss, "MetaData", "bottom_level_pressure", botpressure)

  do ivar = 1, nvars
    !get the name of input variable in geovals
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
        delp4=(modelpressures%vals(kk,iobs)-modelpressures%vals(kk+1,iobs))*1.0e-3_kind_real
        g = g + modelozone%vals(kk,iobs)*rozcon*delz*delp4
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
