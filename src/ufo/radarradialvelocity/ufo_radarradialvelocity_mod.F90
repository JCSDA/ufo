! (C) Copyright 2017-2018 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

!> Fortran module for radarradialvelocity observation operator

module ufo_radarradialvelocity_mod

 use fckit_configuration_module, only: fckit_configuration
 use kinds
 use vert_interp_mod

 use ufo_geovals_mod, only: ufo_geovals, ufo_geoval, ufo_geovals_get_var
 use ufo_vars_mod
 use obsspace_mod

 implicit none
 private

!> Fortran derived type for the observation type
! DONE
 type, public :: ufo_radarradialvelocity
 private
   integer, public :: nvars_in, nvars_out
   character(len=MAXVARLEN), public, allocatable :: varin(:)
   character(len=MAXVARLEN), public, allocatable :: varout(:)
   character(len=MAXVARLEN), public :: v_coord ! GeoVaL to use to interpolate in vertical
 contains
   procedure :: setup  => ufo_radarradialvelocity_setup
   procedure :: simobs => ufo_radarradialvelocity_simobs
   final :: destructor
 end type ufo_radarradialvelocity

 character(len=maxvarlen), dimension(3), parameter :: varin_default = (/var_u, &
                                                                        var_v, &
                                                                        var_w  /)


contains

! ------------------------------------------------------------------------------
! Done
subroutine ufo_radarradialvelocity_setup(self, yaml_conf, vars)
implicit none
class(ufo_radarradialvelocity), intent(inout)     :: self
type(fckit_configuration), intent(in) :: yaml_conf
character(len=MAXVARLEN), dimension(:), intent(inout) :: vars
character(kind=c_char,len=:), allocatable :: coord_name

  self%nvars_out = size(vars)
  allocate(self%varout(self%nvars_out))
  self%varout = vars

  self%nvars_in  = size(varin_default)
  allocate(self%varin(self%nvars_in+1))
  self%varin(1:self%nvars_in) = varin_default

  if( yaml_conf%has("VertCoord") ) then
      call yaml_conf%get_or_die("VertCoord",coord_name)
      self%v_coord = coord_name
      if( trim(self%v_coord) .ne. var_z ) then
        call abor1_ftn("ufo_radarradialvelocity: incorrect vertical coordinate specified")
      endif
  else  ! default
      self%v_coord = var_z
  endif

  self%varin(self%nvars_in+1) = self%v_coord

end subroutine ufo_radarradialvelocity_setup

! ------------------------------------------------------------------------------
! Done
subroutine destructor(self)
implicit none
type(ufo_radarradialvelocity), intent(inout) :: self

  if (allocated(self%varout)) deallocate(self%varout)
  if (allocated(self%varin))  deallocate(self%varin)

end subroutine destructor

! ------------------------------------------------------------------------------
! TODO: put code for your nonlinear observation operator in this routine
! Code in this routine is for radar radialvelocity only

subroutine ufo_radarradialvelocity_simobs(self, geovals, obss, nvars, nlocs, hofx)

  implicit none
  class(ufo_radarradialvelocity), intent(in)    :: self
  integer, intent(in)               :: nvars, nlocs
  type(ufo_geovals),  intent(in)    :: geovals
  real(c_double),     intent(inout) :: hofx(nvars, nlocs)
  type(c_ptr), value, intent(in)    :: obss

  ! Local variables
  integer :: iobs, ivar
  real(kind_real),  dimension(:), allocatable :: obsvcoord
  real(kind_real),  dimension(:), allocatable :: radarazim, radartilt, radardir, vterminal
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
  allocate(radarazim(nlocs))
  allocate(radartilt(nlocs))
  allocate(radardir(nlocs))
  allocate(vterminal(nlocs))
  call obsspace_get_db(obss, "MetaData", "height", obsvcoord)
  call obsspace_get_db(obss, "MetaData", "radar_azimuth", radarazim)
  call obsspace_get_db(obss, "MetaData", "radar_tilt", radartilt)
  call obsspace_get_db(obss, "MetaData", "radar_dir3", radardir)
  call obsspace_get_db(obss, "MetaData", "vterminal", vterminal)

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

  allocate(vfields(self%nvars_in,nlocs))

  do ivar = 1, self%nvars_in
! Get the name of input variable in geovals
    geovar = self%varin(ivar)

! Get profile for this variable from geovals
    call ufo_geovals_get_var(geovals, geovar, profile)

! Interpolate from geovals to observational location into hofx
    do iobs = 1, nlocs
      call vert_interp_apply(profile%nval, profile%vals(:,iobs), &
                             & vfields(ivar,iobs), wi(iobs), wf(iobs))
    enddo
  enddo

  do ivar = 1, self%nvars_out
    do iobs=1,nlocs
      hofx(ivar,iobs) = vfields(1,iobs)*radarazim(iobs) &
                      + vfields(2,iobs)*radartilt(iobs) &
                      + (vfields(3,iobs)-vterminal(iobs))*radardir(iobs)
    enddo
  end do

! Cleanup memory
  deallocate(obsvcoord)
  deallocate(radarazim)
  deallocate(radartilt)
  deallocate(radardir )
  deallocate(vterminal)
  deallocate(wi)
  deallocate(wf)

  deallocate(vfields)

end subroutine ufo_radarradialvelocity_simobs


! ------------------------------------------------------------------------------

end module ufo_radarradialvelocity_mod
