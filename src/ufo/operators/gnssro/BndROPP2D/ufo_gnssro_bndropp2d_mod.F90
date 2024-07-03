! (C) Copyright 2017-2018 UCAR
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 

!> Fortran module for gnssro bending angle ropp2d forward operator
!> following the ROPP (2018 Aug) implementation

module ufo_gnssro_bndropp2d_mod

use fckit_configuration_module, only: fckit_configuration 
!use iso_c_binding
use kinds
use ufo_vars_mod
use ufo_geovals_mod
use ufo_geovals_mod_c, only: ufo_geovals_registry
use ufo_basis_mod,     only: ufo_basis
use vert_interp_mod
use lag_interp_mod,    only: lag_interp_const, lag_interp_smthWeights
use obsspace_mod
use missing_values_mod
use ufo_gnssro_ropp2d_utils_mod
use ufo_gnssro_ropp1d_utils_mod

use gnssro_mod_conf
use fckit_log_module,  only : fckit_log

implicit none
public             :: ufo_gnssro_bndropp2d
private
  !> Fortran derived type for gnssro trajectory
type, extends(ufo_basis) :: ufo_gnssro_BndROPP2D
  type(gnssro_conf)  :: roconf
  real(kind_real), allocatable  :: obsLon2d(:), obsLat2d(:)  !2d location
  contains
    procedure :: setup     => ufo_gnssro_bndropp2d_setup
    procedure :: simobs    => ufo_gnssro_bndropp2d_simobs
end type ufo_gnssro_BndROPP2D

contains

! ------------------------------------------------------------------------------
subroutine ufo_gnssro_bndropp2d_setup(self, f_conf, c_size)
  implicit none
  class(ufo_gnssro_BndROPP2D), intent(inout) :: self
  type(fckit_configuration),   intent(in)    :: f_conf
  integer,                     intent(in)    :: c_size ! 1d obsspace vector length

  call gnssro_conf_setup(self%roconf,f_conf)

  allocate(self%obsLon2d(c_size*self%roconf%n_horiz))
  allocate(self%obsLat2d(c_size*self%roconf%n_horiz))

  self%obsLon2d = 0.0
  self%obsLat2d = 0.0

end subroutine ufo_gnssro_bndropp2d_setup

! ------------------------------------------------------------------------------
subroutine ufo_gnssro_bndropp2d_simobs(self, geovals, hofx, obss)
  use ropp_fm_types, only: State2dFM, State1dFM
  use ropp_fm_types, only: Obs1dBangle, Obs1dRefrac
  use geodesy,       only: geometric2geopotential
  use typesizes,     only: wp => EightByteReal
  use datetimetypes, only: dp

  implicit none
  class(ufo_gnssro_BndROPP2D), intent(in)    :: self
  type(ufo_geovals),           intent(in)    :: geovals
  real(kind_real),             intent(inout) :: hofx(:)
  type(c_ptr), value,          intent(in)    :: obss
  real(c_double)                             :: missing

  type(State2dFM)                    :: x
  type(State1dFM)                    :: x1d
  type(Obs1dBangle)                  :: y
  type(Obs1dRefrac)                  :: y2    ! Observation vector (levels required)

  character(len=*), parameter   :: myname_="ufo_gnssro_bndropp2d_simobs"
  integer, parameter            :: max_string = 800
  character(max_string)         :: err_msg
  integer                       :: nlev, nlocs, iobs, nvprof
  integer                       :: iflip, i
  type(ufo_geoval), pointer     :: t, q, prs, gph, gph_sfc
  real(kind_real), allocatable  :: obsImpP(:),obsLocR(:),obsGeoid(:),obsAzim(:) !nlocs
  real(kind_real), allocatable  :: obsLat(:),obsLon(:)                          !nlocs
  real(kind_real), allocatable  :: obsLonnh(:),obsLatnh(:)                      !n_horiz
  integer                       :: use_compress
  integer                       :: n_horiz
  real(kind_real)               :: dtheta
  real(kind_real)                    :: ob_time
  integer, allocatable, dimension(:) :: ichk
  character(len=20)             :: ro_type
  real(kind_real), allocatable  :: obsAlt(:),obsRef(:),geop(:)                  !nlocs

  n_horiz = self%roconf%n_horiz
  dtheta  = self%roconf%dtheta
  ro_type = self%roconf%ro_type
  use_compress = self%roconf%use_compress

  write(err_msg,*) "TRACE: ufo_gnssro_bndropp2d_simobs: begin"
  call fckit_log%debug(err_msg)

#ifndef ropp_aro
! airborne ro_type can only be used with the ROPP ropp_fm_bangle_2d_aro routine.
! This routine is only avaible in certain ROPP distributions. When this code
! is available, the compilation must include the -Dropp_aro argument
! in the command line, to link with the ropp_fm_bangle_2d_aro routine.
! Abort now if ro_type is airborne but the routine is not available
  if ( ro_type .eq. "airborne" ) then
     write(err_msg,*) myname_, ' ERROR: option "ro_type = airborne"',&
   ' requires compiling UFO with the "-Dropp_aro" cpp argument'
     call abor1_ftn(err_msg)
  endif
#endif

! get variables from geovals
  call ufo_geovals_get_var(geovals, var_ts,    t)         ! temperature
  call ufo_geovals_get_var(geovals, var_q,     q)         ! specific humidity
  call ufo_geovals_get_var(geovals, var_prs,   prs)       ! pressure
  call ufo_geovals_get_var(geovals, var_z,     gph)       ! geopotential height
  call ufo_geovals_get_var(geovals, var_sfc_geomz, gph_sfc)   ! surface geopotential height

! check if the number of geoval profiles is correct
  if (t%nprofiles /= size(hofx)*n_horiz .or. q%nprofiles /= size(hofx)*n_horiz .or. &
      prs%nprofiles /= size(hofx)*n_horiz .or. gph%nprofiles /= size(hofx)*n_horiz .or. &
      gph_sfc%nprofiles /= size(hofx)*n_horiz) then
     write(err_msg,*) myname_, ' error: npaths inconsistent!'
     call abor1_ftn(err_msg)
  endif

  missing = missing_value(missing)
  nlev    = t%nval ! number of model levels
  nlocs   = obsspace_get_nlocs(obss)

  iflip = 0
  if (prs%vals(1,1) .lt. prs%vals(prs%nval,1) ) then
    iflip = 1
    write(err_msg,'(a)') '  ufo_gnssro_bndropp2d_simobs:'//new_line('a')//                         &
                         '  Model vertical height profile is in descending order,'//new_line('a')// &
                         '  but ROPP requires it to be ascending order, need flip'
    call fckit_log%debug(err_msg)
  end if

! set obs space struture
  allocate(obsLon(nlocs))
  allocate(obsLat(nlocs))
  allocate(obsImpP(nlocs))
  allocate(obsLocR(nlocs))
  allocate(obsGeoid(nlocs))
  allocate(obsAzim(nlocs))
  allocate(obsLatnh(n_horiz))
  allocate(obsLonnh(n_horiz))
  allocate(obsAlt(nlocs))
  allocate(obsRef(nlocs))
  allocate(y2%refrac(nlocs))
  allocate(y2%geop(nlocs))
  allocate(geop(nlocs))

  call obsspace_get_db(obss, "MetaData", "longitude",               obsLon)
  call obsspace_get_db(obss, "MetaData", "latitude",                obsLat)
  call obsspace_get_db(obss, "MetaData", "impactParameterRO",       obsImpP)
  call obsspace_get_db(obss, "MetaData", "earthRadiusCurvature",    obsLocR)
  call obsspace_get_db(obss, "MetaData", "geoidUndulation",         obsGeoid)
  call obsspace_get_db(obss, "MetaData", "sensorAzimuthAngle",      obsAzim)
  call obsspace_get_db(obss, "MetaData", "height",                  obsAlt)
  call obsspace_get_db(obss, "ObsValue", "atmosphericRefractivity", obsRef)

  do i = 1, nlocs
    geop(i) = geometric2geopotential(obsLat(i), obsAlt(i))
  enddo

  y2%refrac = obsRef(:)
  y2%geop = geop(:)

  nvprof  = 1  ! no. of bending angles in profile
  ob_time = 0.0
  allocate(ichk(nvprof))
  ichk(:) = 0
  write(err_msg,*) "TRACE: ufo_gnssro_bndropp2d_simobs: begin observation loop, nlocs =  ", nlocs
  call fckit_log%debug(err_msg)

! loop through the obs
  obs_loop: do iobs = 1, nlocs

    if ( ( obsImpP(iobs)-obsLocR(iobs)-obsGeoid(iobs) ) <= self%roconf%top_2d .and. &
           obsAzim(iobs) /= missing ) then

      obsLatnh = self%obsLat2d( (iobs-1)*n_horiz+1:iobs*n_horiz )
      obsLonnh = self%obsLon2d( (iobs-1)*n_horiz+1:iobs*n_horiz )
      call init_ropp_2d_statevec(obsLonnh, obsLatnh,                  &
                      t%vals(:,(iobs-1)*n_horiz+1:iobs*n_horiz),      &
                      q%vals(:,(iobs-1)*n_horiz+1:iobs*n_horiz),      &
                    prs%vals(:,(iobs-1)*n_horiz+1:iobs*n_horiz),      &
                    gph%vals(:,(iobs-1)*n_horiz+1:iobs*n_horiz),      &
                    nlev,x,n_horiz,dtheta,iflip,use_compress)
      call init_ropp_2d_obvec(nvprof,      &
                    obsImpP(iobs),         &
                     obsLat(iobs),         &
                     obsLon(iobs),         &
                    obsLocR(iobs),         &
                   obsGeoid(iobs),         &
                             y)

      if ( ro_type .eq. "airborne" ) then
#ifdef ropp_aro
        call ropp_fm_bangle_2d_aro(x,y,y2)
#else
        write(err_msg,*) myname_, ' ERROR: option "ro_type = airborne"',& 
      ' requires compiling UFO with the "-Dropp_aro" cpp argument'
        call abor1_ftn(err_msg)
#endif
      else
        call ropp_fm_bangle_2d(x,y)
      end if

    else ! apply ropp1d above top_2d

      call init_ropp_1d_statevec(ob_time,            &
                               obsLon(iobs),         &
                               obsLat(iobs),         &
                               t%vals(:,(iobs-1)*n_horiz+1+(n_horiz-1)/2),       &
                               q%vals(:,(iobs-1)*n_horiz+1+(n_horiz-1)/2),       &
                               prs%vals(:,(iobs-1)*n_horiz+1+(n_horiz-1)/2),     &
                               gph%vals(:,(iobs-1)*n_horiz+1+(n_horiz-1)/2),     &
                               nlev,                                             &
                               gph_sfc%vals(1,(iobs-1)*n_horiz+1+(n_horiz-1)/2), &
                               x1d, iflip, use_compress)

      call init_ropp_1d_obvec(nvprof,          &
                              obsImpP(iobs),   &
                              ichk, ob_time,   &
                              obsLat(iobs),    &
                              obsLon(iobs),    &
                              obsLocR(iobs),   &
                              obsGeoid(iobs),  &
                              y)

      call ropp_fm_bangle_1d(x1d,y)

    end if

!   hack -- handling ropp missing value
    if (y%bangle(nvprof) .lt. -900.0_wp ) then
       hofx(iobs) = missing
       y%bangle(nvprof) = missing
    else
       hofx(iobs) = y%bangle(nvprof)  ! nvprof is just one point
    end if
!   hack -- handling ropp missing value

!   deallocate ropp structures
    if ( ( obsImpP(iobs)-obsLocR(iobs)-obsGeoid(iobs) ) <= self%roconf%top_2d .and. &
           obsAzim(iobs) /= missing ) then
      call ropp_tidy_up_2d(x,y)
    else
      call ropp_tidy_up_1d(x1d,y)
    end if

  end do obs_loop

  deallocate(obsLat)
  deallocate(obsLon)
  deallocate(obsImpP)
  deallocate(obsLocR)
  deallocate(obsGeoid)
  deallocate(obsAzim)
  deallocate(obsLatnh)
  deallocate(obsLonnh)
  deallocate(ichk)
  deallocate(obsAlt)
  deallocate(obsRef)
  deallocate(y2%refrac)
  deallocate(y2%geop)

  write(err_msg,*) "TRACE: ufo_gnssro_bndropp2d_simobs: complete"
  call fckit_log%debug(err_msg)

end subroutine ufo_gnssro_bndropp2d_simobs
! ------------------------------------------------------------------------------

end module ufo_gnssro_bndropp2d_mod
