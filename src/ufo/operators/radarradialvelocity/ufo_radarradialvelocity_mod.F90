! (C) Copyright 2017-2018 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

!> Fortran module for radarradialvelocity observation operator

module ufo_radarradialvelocity_mod

 use oops_variables_mod
 use obs_variables_mod
 use ufo_vars_mod

 implicit none
 private

!> Fortran derived type for the observation type
 type, public :: ufo_radarradialvelocity
 private
   type(obs_variables), public :: obsvars
   type(oops_variables), public :: geovars
   character(len=MAXVARLEN), public :: v_coord ! GeoVaL to use to interpolate in vertical
 contains
   procedure :: setup  => ufo_radarradialvelocity_setup
   procedure :: simobs => ufo_radarradialvelocity_simobs
 end type ufo_radarradialvelocity

 character(len=maxvarlen), dimension(3), parameter :: geovars_default = (/var_u, &
                                                                        var_v, &
                                                                        var_w  /)


contains

! ------------------------------------------------------------------------------
subroutine ufo_radarradialvelocity_setup(self, yaml_conf)
use fckit_configuration_module, only: fckit_configuration
use iso_c_binding
implicit none
class(ufo_radarradialvelocity), intent(inout)     :: self
type(fckit_configuration), intent(in) :: yaml_conf

character(kind=c_char,len=:), allocatable :: coord_name

  call self%geovars%push_back(geovars_default)

  call yaml_conf%get_or_die("VertCoord",coord_name)
  self%v_coord = coord_name
  if( trim(self%v_coord) .ne. var_z .and. trim(self%v_coord) .ne. var_zm ) then
      call abor1_ftn("ufo_radarradialvelocity: incorrect vertical coordinate specified")
  endif

  call self%geovars%push_back(self%v_coord)

end subroutine ufo_radarradialvelocity_setup

! ------------------------------------------------------------------------------
! Code in this routine is for radar radialvelocity only
subroutine ufo_radarradialvelocity_simobs(self, geovals, obss, nvars, nlocs, hofx)
  use kinds
  use vert_interp_mod
  use ufo_geovals_mod, only: ufo_geovals, ufo_geoval, ufo_geovals_get_var
  use obsspace_mod
  implicit none
  class(ufo_radarradialvelocity), intent(in)    :: self
  integer, intent(in)               :: nvars, nlocs
  type(ufo_geovals),  intent(in)    :: geovals
  real(c_double),     intent(inout) :: hofx(nvars, nlocs)
  type(c_ptr), value, intent(in)    :: obss

  ! Local variables
  integer :: iobs, ivar, nvars_geovars
  real(kind_real),  dimension(:), allocatable :: obsvcoord
  real(kind_real),  dimension(:), allocatable :: cosazm_costilt, sinazm_costilt, sintilt, vterminal
  type(ufo_geoval), pointer :: vcoordprofile, profile
  real(kind_real),  allocatable :: wf(:)
  integer,          allocatable :: wi(:)

  character(len=MAXVARLEN) :: geovar

  real(kind_real), allocatable :: tmp(:)
  real(kind_real) :: tmp2

  real(kind_real), allocatable :: vfields(:,:)  ! background fields (u,v,w) interplated vertically to the observation height

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

! Get heigh profiles from geovals
  call ufo_geovals_get_var(geovals, self%v_coord, vcoordprofile)

! Get the observation vertical coordinates
  allocate(obsvcoord(nlocs))
  allocate(cosazm_costilt(nlocs))
  allocate(sinazm_costilt(nlocs))
  allocate(sintilt(nlocs))
  allocate(vterminal(nlocs))

  call obsspace_get_db(obss, "MetaData", "height", obsvcoord)
  call obsspace_get_db(obss, "MetaData", "cosazm_costilt", cosazm_costilt)
  call obsspace_get_db(obss, "MetaData", "sinazm_costilt", sinazm_costilt)
  call obsspace_get_db(obss, "MetaData", "sintilt", sintilt)

! put observation operator code here
! Allocate arrays for interpolation weights

  allocate(wi(nlocs))
  allocate(wf(nlocs))

! Calculate the interpolation weights

  allocate(tmp(vcoordprofile%nval))
  do iobs = 1, nlocs
    tmp = vcoordprofile%vals(:,iobs)
    tmp2 = obsvcoord(iobs)
    call vert_interp_weights(vcoordprofile%nval, tmp2, tmp, wi(iobs), wf(iobs))
  enddo

! Number of variables in geovars (without the vertical coordinate)
  nvars_geovars = self%geovars%nvars() - 1
  allocate(vfields(nvars_geovars,nlocs))

  do ivar = 1, nvars_geovars
! Get the name of input variable in geovals
    geovar = self%geovars%variable(ivar)

! Get profile for this variable from geovals
    call ufo_geovals_get_var(geovals, geovar, profile)

! Interpolate from geovals to observational location into hofx
    do iobs = 1, nlocs
      call vert_interp_apply(profile%nval, profile%vals(:,iobs), &
                             & vfields(ivar,iobs), wi(iobs), wf(iobs))
    enddo
  enddo

  vterminal=0.0
  do ivar = 1, nvars
    do iobs=1,nlocs
      hofx(ivar,iobs) = vfields(1,iobs)*cosazm_costilt(iobs) &
                      + vfields(2,iobs)*sinazm_costilt(iobs) &
                      + (vfields(3,iobs)-vterminal(iobs))*sintilt(iobs)
    enddo
  end do

! Cleanup memory
  deallocate(obsvcoord)
  deallocate(cosazm_costilt)
  deallocate(sinazm_costilt)
  deallocate(sintilt )
  deallocate(vterminal)
  deallocate(wi)
  deallocate(wf)

  deallocate(vfields)

end subroutine ufo_radarradialvelocity_simobs


! ------------------------------------------------------------------------------

end module ufo_radarradialvelocity_mod
