! (C) Copyright 2017-2018 UCAR
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 

!> Fortran module for gnssro bending angle ropp2d forward operator
!> following the ROPP (2018 Aug) implementation

module ufo_gnssro_bndropp2d_mod

use iso_c_binding
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
use fckit_log_module,  only : fckit_log

implicit none
public             :: ufo_gnssro_bndropp2d
private

  !> Fortran derived type for gnssro trajectory
type, extends(ufo_basis) :: ufo_gnssro_BndROPP2D
  contains
    procedure :: simobs    => ufo_gnssro_bndropp2d_simobs
end type ufo_gnssro_BndROPP2D

contains

! ------------------------------------------------------------------------------
! ------------------------------------------------------------------------------
subroutine ufo_gnssro_bndropp2d_simobs(self, geovals, hofx, obss)
  use ropp_fm_types, only: State2dFM
  use ropp_fm_types, only: Obs1dBangle
  use typesizes,     only: wp => EightByteReal
  use datetimetypes, only: dp

  implicit none
  class(ufo_gnssro_BndROPP2D), intent(in)    :: self
  type(ufo_geovals),           intent(in)    :: geovals
  real(kind_real),             intent(inout) :: hofx(:)
  type(c_ptr), value,          intent(in)    :: obss
  real(c_double)                             :: missing

  type(State2dFM)                    :: x
  type(Obs1dBangle)                  :: y

  character(len=*), parameter     :: myname_="ufo_gnssro_bndropp2d_simobs"
  integer, parameter              :: max_string = 800
  character(max_string)           :: err_msg
  integer                         :: nlev, nobs, iobs, nvprof
  integer                         :: ierr
  real(kind=dp)                   :: ob_time
  type(ufo_geoval), pointer          :: t, q, prs, gph !, gph_sfc
  real(kind_real), allocatable       :: obsLat(:), obsLon(:), obsImpP(:), obsLocR(:), obsGeoid(:)
  integer                            :: iflip
  integer,                 parameter :: n_horiz=3

  write(err_msg,*) "TRACE: ufo_gnssro_bndropp2d_simobs: begin"
  call fckit_log%info(err_msg)

  print*,  "geovals%nobs & size(hofx)*n_horiz)", geovals%nobs, size(hofx)

! check if nobs is consistent in geovals & hofx
  if (geovals%nobs /= size(hofx)) then
      write(err_msg,*) myname_, ' error: nobs inconsistent!'
      call abor1_ftn(err_msg)
  endif

! get variables from geovals
  call ufo_geovals_get_var(geovals, var_t,     t)         ! temperature
  call ufo_geovals_get_var(geovals, var_q,     q)         ! specific humidity
  call ufo_geovals_get_var(geovals, var_prs,   prs)       ! pressure
  call ufo_geovals_get_var(geovals, var_z,     gph)       ! geopotential height
!  call ufo_geovals_get_var(geovals, var_sfc_z, gph_sfc)   ! surface geopotential height

  missing = missing_value(missing)

  nlev  = t%nval ! number of model levels
  nobs  = obsspace_get_nlocs(obss)

  iflip = 0
  if (prs%vals(1,1) .lt. prs%vals(prs%nval,1) ) then
    iflip = 1
    write(err_msg,'(a)') '  ufo_gnssro_bndropp2d_simobs:'//new_line('a')//                         &
                         '  Model vertical height profile is in descending order,'//new_line('a')// &
                         '  but ROPP requires it to be ascending order, need flip'
    call fckit_log%info(err_msg)
  end if
! set obs space struture
  allocate(obsLon(nobs))
  allocate(obsLat(nobs))
  allocate(obsImpP(nobs))
  allocate(obsLocR(nobs))
  allocate(obsGeoid(nobs))

  call obsspace_get_db(obss, " ", "longitude",        obsLon)
  call obsspace_get_db(obss, " ", "latitude",         obsLat)
  call obsspace_get_db(obss, " ", "impact_parameter", obsImpP)
  call obsspace_get_db(obss, " ", "earth_radius_of_curvature", obsLocR)
  call obsspace_get_db(obss, " ", "geoid_height_above_reference_ellipsoid", obsGeoid)

  nvprof=1  ! no. of bending angles in profile 

  write(err_msg,*) "TRACE: ufo_gnssro_bndropp2d_simobs: begin observation loop, nobs =  ", nobs
  call fckit_log%info(err_msg)

! loop through the obs
  obs_loop: do iobs = 1, nobs  

    call init_ropp_2d_statevec(          &
                      obsLon(iobs),      &
                      obsLat(iobs),      &
                    t%vals(:,iobs),      &
                    q%vals(:,iobs),      &
                  prs%vals(:,iobs),      &
                  gph%vals(:,iobs),      &
                     nlev,x, iflip)
     
    call init_ropp_2d_obvec(nvprof,      &
                  obsImpP(iobs),         &
                   obsLat(iobs),         &
                   obsLon(iobs),         &
                  obsLocR(iobs),         &
                 obsGeoid(iobs),         &
                             y)

    call ropp_fm_bangle_2d(x,y)

!   hack -- handling ropp missing value 
    if (y%bangle(nvprof) .lt. -900.0_wp ) then
       hofx(iobs) = missing
       y%bangle(nvprof) = missing
    else
       hofx(iobs) = y%bangle(nvprof)  ! nvprof is just one point
    end if
!   hack -- handling ropp missing value 

!   deallocate ropp structures  
    call ropp_tidy_up_2d(x,y)   

  end do obs_loop

  deallocate(obsLat)
  deallocate(obsLon)
  deallocate(obsImpP)
  deallocate(obsLocR)
  deallocate(obsGeoid)

  write(err_msg,*) "TRACE: ufo_gnssro_bndropp2d_simobs: completed"
  call fckit_log%info(err_msg)
     
  return

end subroutine ufo_gnssro_bndropp2d_simobs
! ------------------------------------------------------------------------------

end module ufo_gnssro_bndropp2d_mod
