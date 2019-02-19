! (C) Copyright 2017-2018 UCAR
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 

!> Fortran module for gnssro bending angle ropp1d forward operator
!> following the ROPP (2018 Aug) implementation

module ufo_gnssro_bndropp1d_mod

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
use ufo_gnssro_ropp1d_utils_mod
use fckit_log_module,  only : fckit_log

implicit none
public             :: ufo_gnssro_bndropp1d
private

  !> Fortran derived type for gnssro trajectory
type, extends(ufo_basis) :: ufo_gnssro_BndROPP1D
  contains
    procedure :: simobs    => ufo_gnssro_bndropp1d_simobs
end type ufo_gnssro_BndROPP1D

contains

! ------------------------------------------------------------------------------
! ------------------------------------------------------------------------------
subroutine ufo_gnssro_bndropp1d_simobs(self, geovals, hofx, obss)
  use ropp_fm_types, only: State1dFM
  use ropp_fm_types, only: Obs1dBangle
  use typesizes,     only: wp => EightByteReal
  use datetimetypes, only: dp

  implicit none
  class(ufo_gnssro_BndROPP1D), intent(in)    :: self
  type(ufo_geovals),           intent(in)    :: geovals
  real(kind_real),             intent(inout) :: hofx(:)
  type(c_ptr), value,          intent(in)    :: obss
  real(c_double)                             :: missing

  type(State1dFM)                 :: x
  type(Obs1dBangle)               :: y

  character(len=*), parameter     :: myname_="ufo_gnssro_bndropp1d_simobs"
  real(kind=dp)                   :: ob_time
  integer, parameter              :: max_string = 800

  character(max_string)              :: err_msg
  character(len=250)                 :: record
  integer                            :: nlev, nobs, iobs,nvprof, obss_nobs
  integer, allocatable, dimension(:) :: ichk
  type(ufo_geoval), pointer          :: t, q, prs, gph, gph_sfc
  real(kind_real), allocatable       :: obsLat(:), obsLon(:), obsImpP(:), obsLocR(:), obsGeoid(:)

  write(err_msg,*) "TRACE: ufo_gnssro_bndropp1d_simobs: begin"
  call fckit_log%info(err_msg)

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
  call ufo_geovals_get_var(geovals, var_sfc_z, gph_sfc)   ! surface geopotential height

  missing = missing_value(missing)

  nlev  = t%nval ! number of model levels
  nobs  = obsspace_get_nlocs(obss)

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

  nvprof = 1 ! number of vertical profiles (occultation points)
  allocate(ichk(nvprof))
  ichk(:) = 0   ! this will hold QC values for observation from QC flags

  write(err_msg,*) "TRACE: ufo_gnssro_bndropp1d_simobs: begin observation loop, nobs =  ", nobs
  call fckit_log%info(err_msg)

  obs_loop: do iobs = 1, nobs 

    ob_time = 0.0
    call init_ropp_1d_statevec(ob_time,              &
                               obsLon(iobs),         &
                               obsLat(iobs),         &
                               t%vals(:,iobs),       &
                               q%vals(:,iobs),       &
                               prs%vals(:,iobs),     &
                               gph%vals(:,iobs),     &
                               nlev,                 &
                               gph_sfc%vals(1,iobs), &
                               x)

    call init_ropp_1d_obvec(nvprof,          &
                            obsImpP(iobs),   &
                            ichk, ob_time,   &
                            obsLat(iobs),    &
                            obsLon(iobs),    &
                            obsLocR(iobs),   &  
                            obsGeoid(iobs),  &
                            y)

    call ropp_fm_bangle_1d(x,y)

!  hack -- handling ropp missing value 
    if (y%bangle(nvprof) .lt. -900.0_wp ) then
       hofx(iobs) = missing 
       y%bangle(nvprof) = missing
    else
       hofx(iobs) = y%bangle(nvprof)  ! nvprof is just one point
    end if
!  hack -- handling ropp missing value 
   call ropp_tidy_up_1d(x,y)

  end do obs_loop
      
  deallocate(obsLat) 
  deallocate(obsLon)
  deallocate(obsImpP)
  deallocate(obsLocR)
  deallocate(obsGeoid)

  write(err_msg,*) "TRACE: ufo_gnssro_bndropp1d_simobs: completed"
  call fckit_log%info(err_msg)

end subroutine ufo_gnssro_bndropp1d_simobs
! ------------------------------------------------------------------------------

end module ufo_gnssro_bndropp1d_mod
