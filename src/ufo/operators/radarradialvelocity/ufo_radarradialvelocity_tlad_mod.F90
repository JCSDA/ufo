! (C) Copyright 2017-2019 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

module ufo_radarradialvelocity_tlad_mod

  use oops_variables_mod
  use obs_variables_mod
  use ufo_vars_mod
  use ufo_geovals_mod
  use vert_interp_mod
  use missing_values_mod


! ------------------------------------------------------------------------------

  type, public :: ufo_radarradialvelocity_tlad
  private
    type(obs_variables), public :: obsvars
    type(oops_variables), public :: geovars
    integer :: nval, nlocs
    character(len=MAXVARLEN), public :: v_coord ! GeoVaL to use to interpolate in vertical
    real(kind_real), allocatable :: wf(:)
    integer, allocatable :: wi(:)
    real(kind_real),  dimension(:), allocatable :: cosazm_costilt, sinazm_costilt, sintilt, vterminal
  contains
    procedure :: setup => radarradialvelocity_tlad_setup_
    procedure :: cleanup => radarradialvelocity_tlad_cleanup_
    procedure :: settraj => radarradialvelocity_tlad_settraj_
    procedure :: simobs_tl => radarradialvelocity_simobs_tl_
    procedure :: simobs_ad => radarradialvelocity_simobs_ad_
    final :: destructor
  end type ufo_radarradialvelocity_tlad

  character(len=maxvarlen), dimension(3), parameter :: geovars_default = (/var_u, &
                                                                        var_v, &
                                                                        var_w /)

! ------------------------------------------------------------------------------
contains
! ------------------------------------------------------------------------------

subroutine radarradialvelocity_tlad_setup_(self, yaml_conf)
  use fckit_configuration_module, only: fckit_configuration
  implicit none
  class(ufo_radarradialvelocity_tlad), intent(inout) :: self
  type(fckit_configuration), intent(in) :: yaml_conf

  character(kind=c_char,len=:), allocatable :: coord_name

  integer :: ivar, nvars
  call self%geovars%push_back(geovars_default)

  if( yaml_conf%has("VertCoord") ) then
      call yaml_conf%get_or_die("VertCoord",coord_name)
      self%v_coord = coord_name
      if( trim(self%v_coord) .ne. var_z .and. trim(self%v_coord) .ne. var_zm ) then
        call abor1_ftn("ufo_radarradialvelocity: incorrect vertical coordinate specified")
      endif
  else  ! default
      !self%v_coord = var_z
      self%v_coord = var_zm
  endif


end subroutine radarradialvelocity_tlad_setup_

! ------------------------------------------------------------------------------

subroutine radarradialvelocity_tlad_settraj_(self, geovals, obss)
  use obsspace_mod
  implicit none
  class(ufo_radarradialvelocity_tlad), intent(inout) :: self
  type(ufo_geovals),         intent(in)    :: geovals
  type(c_ptr), value,        intent(in)    :: obss

  !local variables
  integer :: iobs, ivar, nvars_geovars
  real(kind_real),  dimension(:), allocatable :: obsvcoord
  !real(kind_real),  dimension(:), allocatable :: cosazm_costilt, sinazm_costilt, sintilt, vterminal
  type(ufo_geoval), pointer :: vcoordprofile, profile  

  real(kind_real), allocatable :: tmp(:)
  real(kind_real) :: tmp2
  integer:: i

  real(kind_real), allocatable :: vfields(:,:)  ! background fields (u,v,w) interplated vertically to the observation height

  ! Make sure nothing already allocated
  call self%cleanup()

  ! Get height profiles from geovals
  call ufo_geovals_get_var(geovals, self%v_coord, vcoordprofile)
  self%nval = vcoordprofile%nval

  ! Get the observation vertical coordinates
  self%nlocs = obsspace_get_nlocs(obss)
  allocate(obsvcoord(self%nlocs))
  allocate(self%cosazm_costilt(self%nlocs))
  allocate(self%sinazm_costilt(self%nlocs))
  allocate(self%sintilt(self%nlocs))
  allocate(self%vterminal(self%nlocs))

  call obsspace_get_db(obss, "MetaData", "height", obsvcoord)
  call obsspace_get_db(obss, "MetaData", "cosAzimuthCosTilt", self%cosazm_costilt)
  call obsspace_get_db(obss, "MetaData", "sinAzimuthCosTilt", self%sinazm_costilt)
  call obsspace_get_db(obss, "MetaData", "sinTilt", self%sintilt)
! call obsspace_get_db(obss, "MetaData", "vterminal", self%vterminal)

  ! Allocate arrays for interpolation weights
  allocate(self%wi(self%nlocs))
  allocate(self%wf(self%nlocs))

  ! Calculate the interpolation weights
  allocate(tmp(vcoordprofile%nval))
  do iobs = 1, self%nlocs
    tmp = vcoordprofile%vals(:,iobs)
    tmp2 = obsvcoord(iobs)
    call vert_interp_weights(vcoordprofile%nval, tmp2, tmp, self%wi(iobs), self%wf(iobs))
  enddo

! Cleanup memory
  deallocate(obsvcoord)
  deallocate(tmp)

end subroutine radarradialvelocity_tlad_settraj_

! ------------------------------------------------------------------------------

subroutine radarradialvelocity_simobs_tl_(self, geovals, obss, nvars, nlocs, hofx)
  implicit none
  class(ufo_radarradialvelocity_tlad), intent(in) :: self
  type(ufo_geovals),         intent(in) :: geovals
  integer,                   intent(in) :: nvars, nlocs
  real(c_double),         intent(inout) :: hofx(nvars, nlocs)
  type(c_ptr), value,        intent(in) :: obss

  integer :: iobs, ivar, nvars_geovars
  type(ufo_geoval), pointer :: profile
  character(len=MAXVARLEN) :: geovar

  real(kind_real), allocatable :: vfields(:,:)  ! background fields (u,v,w) interplated vertically to the observation height


! Number of variables in geovars (without the vertical coordinate)
  nvars_geovars = self%geovars%nvars() - 1
  allocate(vfields(nvars_geovars,nlocs))
  vfields=0.0

  do ivar = 1, nvars_geovars
    geovar = self%geovars%variable(ivar)

! Get profile for this variable from geovals
    call ufo_geovals_get_var(geovals, geovar, profile)

! Interpolate from geovals to observational location into hofx
    do iobs = 1, nlocs
      call vert_interp_apply_tl(profile%nval, profile%vals(:,iobs), &
                             & vfields(ivar,iobs), self%wi(iobs), self%wf(iobs))
    enddo
  enddo

  do ivar = 1, nvars
    do iobs=1,nlocs
      hofx(ivar,iobs) = vfields(1,iobs)*self%cosazm_costilt(iobs) &
                      + vfields(2,iobs)*self%sinazm_costilt(iobs)
    enddo
  end do

  deallocate(vfields)

end subroutine radarradialvelocity_simobs_tl_

! ------------------------------------------------------------------------------

subroutine radarradialvelocity_simobs_ad_(self, geovals, obss, nvars, nlocs, hofx)
  implicit none
  class(ufo_radarradialvelocity_tlad), intent(in) :: self
  type(ufo_geovals),         intent(inout) :: geovals
  integer,                   intent(in)    :: nvars, nlocs
  real(c_double),            intent(in)    :: hofx(nvars, nlocs)
  type(c_ptr), value,        intent(in)    :: obss

  integer :: iobs, ivar, nvars_geovars
  type(ufo_geoval), pointer :: profile
  character(len=MAXVARLEN) :: geovar
  real(c_double) :: missing
  real(kind_real), allocatable :: vfields(:,:)  ! background fields (u,v,w) interplated vertically to the observation height

  missing = missing_value(missing)

! Number of variables in geovars (without the vertical coordinate)
  nvars_geovars = self%geovars%nvars() - 1
  allocate(vfields(nvars_geovars,nlocs))

  vfields=0.0
  do ivar = 1, nvars
    do iobs=1,nlocs
     ! no vertical velocity and terminal velocity in GSI rw observer, it can add
     ! in future after acceptance test
     if (hofx(ivar,iobs) .ne. missing) then
      vfields(1,iobs) = vfields(1,iobs) + hofx(ivar,iobs)*self%cosazm_costilt(iobs)
      vfields(2,iobs) = vfields(2,iobs) + hofx(ivar,iobs)*self%sinazm_costilt(iobs)
     end if
    enddo
  end do

  do ivar = 1, nvars_geovars
! Get the name of input variable in geovals
    geovar = self%geovars%variable(ivar)

! Get profile for this variable from geovals
    call ufo_geovals_get_var(geovals, geovar, profile)

! Interpolate from geovals to observational location into hofx
    do iobs = 1, nlocs
      call vert_interp_apply_ad(profile%nval, profile%vals(:,iobs), &
                             & vfields(ivar,iobs), self%wi(iobs), self%wf(iobs))
    enddo
  enddo


  deallocate(vfields)
end subroutine radarradialvelocity_simobs_ad_

! ------------------------------------------------------------------------------

subroutine radarradialvelocity_tlad_cleanup_(self)
  implicit none
  class(ufo_radarradialvelocity_tlad), intent(inout) :: self
  self%nval = 0
  self%nlocs = 0
  if (allocated(self%wi)) deallocate(self%wi)
  if (allocated(self%wf)) deallocate(self%wf)
  if (allocated(self%cosazm_costilt)) deallocate(self%cosazm_costilt)
  if (allocated(self%sinazm_costilt)) deallocate(self%sinazm_costilt)
  if (allocated(self%sintilt)) deallocate(self%sintilt)
  if (allocated(self%vterminal)) deallocate(self%vterminal)
end subroutine radarradialvelocity_tlad_cleanup_

! ------------------------------------------------------------------------------

subroutine  destructor(self)
  type(ufo_radarradialvelocity_tlad), intent(inout)  :: self

  call self%cleanup()

end subroutine destructor

! ------------------------------------------------------------------------------

end module ufo_radarradialvelocity_tlad_mod
