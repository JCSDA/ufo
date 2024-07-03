! (C) Copyright 2017-2018 UCAR
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 

!> Fortran module for gnssro bending angle ropp1d forward operator
!> following the ROPP (2018 Aug) implementation

module ufo_gnssro_bndropp1d_mod

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
use ufo_gnssro_ropp1d_utils_mod
use fckit_log_module,  only : fckit_log
use gnssro_mod_conf

implicit none
public             :: ufo_gnssro_bndropp1d
private

  !> Fortran derived type for gnssro trajectory
type, extends(ufo_basis) :: ufo_gnssro_BndROPP1D
  type(gnssro_conf)  :: roconf
  contains
    procedure :: setup     => ufo_gnssro_bndropp1d_setup
    procedure :: simobs    => ufo_gnssro_bndropp1d_simobs
end type ufo_gnssro_BndROPP1D

contains

! ------------------------------------------------------------------------------
subroutine ufo_gnssro_bndropp1d_setup(self, f_conf)
  implicit none
  class(ufo_gnssro_BndROPP1D), intent(inout) :: self
  type(fckit_configuration),   intent(in)    :: f_conf

  call gnssro_conf_setup(self%roconf,f_conf)

end subroutine ufo_gnssro_bndropp1d_setup
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
  integer                            :: nlev, nobs, iobs,nvprof
  integer, allocatable, dimension(:) :: ichk
  type(ufo_geoval), pointer          :: t, q, prs, gph, gph_sfc
  real(kind_real), allocatable       :: obsLat(:), obsLon(:), obsImpP(:), obsLocR(:), obsGeoid(:)
  integer                            :: iflip 
  integer                            :: use_compress

  use_compress = self%roconf%use_compress

  write(err_msg,*) "TRACE: ufo_gnssro_bndropp1d_simobs: begin"
  call fckit_log%debug(err_msg)

! check if nlocs is consistent in geovals & hofx
  if (geovals%nlocs /= size(hofx)) then
      write(err_msg,*) myname_, ' error: nlocs inconsistent!'
      call abor1_ftn(err_msg)
  endif

! get variables from geovals
  call ufo_geovals_get_var(geovals, var_ts,    t)         ! temperature
  call ufo_geovals_get_var(geovals, var_q,     q)         ! specific humidity
  call ufo_geovals_get_var(geovals, var_prs,   prs)       ! pressure
  call ufo_geovals_get_var(geovals, var_z,     gph)       ! geopotential height
  call ufo_geovals_get_var(geovals, var_sfc_geomz, gph_sfc)   ! surface geopotential height

  missing = missing_value(missing)

  nlev  = t%nval ! number of model levels
  nobs  = obsspace_get_nlocs(obss)

  if (nobs > 0) then 
     iflip = 0
     if (prs%vals(1,1) .lt. prs%vals(prs%nval,1) ) then
       iflip = 1 
       write(err_msg,'(a)') '  ufo_gnssro_bndropp1d_simobs:'//new_line('a')//                         &
                            '  Model vertical height profile is in descending order,'//new_line('a')// &
                            '  but ROPP requires it to be ascending order, need flip'
       call fckit_log%debug(err_msg)
     end if

   ! set obs space struture
     allocate(obsLon(nobs))
     allocate(obsLat(nobs))
     allocate(obsImpP(nobs))
     allocate(obsLocR(nobs))
     allocate(obsGeoid(nobs))

     call obsspace_get_db(obss, "MetaData", "longitude",            obsLon) 
     call obsspace_get_db(obss, "MetaData", "latitude",             obsLat) 
     call obsspace_get_db(obss, "MetaData", "impactParameterRO",    obsImpP)
     call obsspace_get_db(obss, "MetaData", "earthRadiusCurvature", obsLocR) 
     call obsspace_get_db(obss, "MetaData", "geoidUndulation",      obsGeoid) 

     nvprof = 1 ! number of vertical profiles (occultation points)
     allocate(ichk(nvprof))
     ichk(:) = 0   ! this will hold QC values for observation from QC flags

     write(err_msg,*) "TRACE: ufo_gnssro_bndropp1d_simobs: begin observation loop, nobs =  ", nobs
     call fckit_log%debug(err_msg)

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
                                  x, iflip, use_compress)

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
 
     deallocate(ichk)
     deallocate(obsLat) 
     deallocate(obsLon)
     deallocate(obsImpP)
     deallocate(obsLocR)
     deallocate(obsGeoid)
  end if ! nobs > 0

  write(err_msg,*) "TRACE: ufo_gnssro_bndropp1d_simobs: complete"
  call fckit_log%debug(err_msg)

end subroutine ufo_gnssro_bndropp1d_simobs
! ------------------------------------------------------------------------------

end module ufo_gnssro_bndropp1d_mod
