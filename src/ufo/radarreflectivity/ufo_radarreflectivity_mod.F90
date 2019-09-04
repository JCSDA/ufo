! (C) Copyright 2017-2018 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

!> Fortran module for radarreflectivity observation operator

module ufo_radarreflectivity_mod

 use fckit_configuration_module, only: fckit_configuration
 use kinds
 use vert_interp_mod

 use ufo_geovals_mod, only: ufo_geovals, ufo_geoval, ufo_geovals_get_var
 use ufo_vars_mod
 use obsspace_mod

 implicit none
 private
 integer, parameter :: max_string=800

!> Fortran derived type for the observation type

 type, public :: ufo_radarreflectivity
 private
   integer, public :: nvars_in, nvars_out
   character(len=max_string), public, allocatable :: varin(:)
   character(len=max_string), public, allocatable :: varout(:)
   character(len=MAXVARLEN), public :: v_coord ! GeoVaL to use to interpolate in vertical
 contains
   procedure :: setup  => ufo_radarreflectivity_setup
   procedure :: simobs => ufo_radarreflectivity_simobs
   final :: destructor
 end type ufo_radarreflectivity

 character(len=maxvarlen), dimension(1), parameter :: varin_default = (/var_refl/)

contains

! ------------------------------------------------------------------------------
! Done
subroutine ufo_radarreflectivity_setup(self, yaml_conf, vars)

  implicit none
  class(ufo_radarreflectivity), intent(inout)     :: self
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
        call abor1_ftn("ufo_radarreflectivity: incorrect vertical coordinate specified")
      endif
  else  ! default
      self%v_coord = var_z
  endif

  self%varin(self%nvars_in+1) = self%v_coord

end subroutine ufo_radarreflectivity_setup

! ------------------------------------------------------------------------------
! Done
subroutine destructor(self)

  implicit none
  type(ufo_radarreflectivity), intent(inout) :: self

  if (allocated(self%varout)) deallocate(self%varout)
  if (allocated(self%varin))  deallocate(self%varin)

end subroutine destructor

! ------------------------------------------------------------------------------
! TODO: put code for your nonlinear observation operator in this routine
! Code in this routine is for radarreflectivity only

subroutine ufo_radarreflectivity_simobs(self, geovals, obss, nvars, nlocs, hofx)

  implicit none
  class(ufo_radarreflectivity), intent(in)    :: self
  integer, intent(in)               :: nvars, nlocs
  type(ufo_geovals),  intent(in)    :: geovals
  real(c_double),     intent(inout) :: hofx(nvars, nlocs)
  type(c_ptr), value, intent(in)    :: obss

  ! Local variables
  type(ufo_geoval), pointer :: geoval
  real(kind_real), dimension(:), allocatable :: obss_metadata

  integer :: iobs, ivar
  real(kind_real),  dimension(:), allocatable :: obsvcoord
  type(ufo_geoval), pointer :: vcoordprofile, profile
  real(kind_real),  allocatable :: wf(:)
  integer,          allocatable :: wi(:)

  character(len=MAXVARLEN) :: geovar

  real(kind_real), allocatable :: tmp(:)
  real(kind_real) :: tmp2

  !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  ! Get height profiles from geovals
  call ufo_geovals_get_var(geovals, self%v_coord, vcoordprofile)


  ! Get the observation vertical coordinates
  allocate(obsvcoord(nlocs))
  call obsspace_get_db(obss, "MetaData", "height", obsvcoord)

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

  do ivar = 1, self%nvars_in
    ! Get the name of input variable in geovals
    geovar = self%varin(ivar)

    ! Get profile for this variable from geovals
    call ufo_geovals_get_var(geovals, geovar, profile)

    ! Interpolate from geovals to observational location into hofx
    do iobs = 1, nlocs
      call vert_interp_apply(profile%nval, profile%vals(:,iobs), &
                             & hofx(ivar,iobs), wi(iobs), wf(iobs))
    enddo
  enddo

  ! Cleanup memory
  deallocate(obsvcoord)
  deallocate(wi)
  deallocate(wf)

end subroutine ufo_radarreflectivity_simobs


! ------------------------------------------------------------------------------

end module ufo_radarreflectivity_mod
