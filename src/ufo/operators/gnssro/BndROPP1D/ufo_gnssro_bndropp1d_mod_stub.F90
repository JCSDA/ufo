! (C) Copyright 2017-2018 UCAR
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 

!>Stubbed Fortran module for gnssro bending angle ropp1d forward operator
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
subroutine ufo_gnssro_bndropp1d_setup(self, f_conf )
  implicit none
  class(ufo_gnssro_BndROPP1D), intent(inout) :: self
  type(fckit_configuration), intent(in)      :: f_conf

  call gnssro_conf_setup(self%roconf,f_conf)

end subroutine ufo_gnssro_bndropp1d_setup
! ------------------------------------------------------------------------------
! ------------------------------------------------------------------------------
subroutine ufo_gnssro_bndropp1d_simobs(self, geovals, hofx, obss)

  implicit none
  class(ufo_gnssro_BndROPP1D), intent(in)    :: self
  type(ufo_geovals),           intent(in)    :: geovals
  real(kind_real),             intent(inout) :: hofx(:)
  type(c_ptr), value,          intent(in)    :: obss
  real(c_double)                             :: missing

  character(len=*), parameter  :: myname_="ufo_gnssro_bndropp1d_simobs"
  integer, parameter           :: max_string = 800

  character(max_string)              :: err_msg
  integer                            :: nlev, nobs, iobs,nvprof, obss_nobs
  type(ufo_geoval), pointer          :: t, q, prs, gph, gph_sfc
  real(kind_real), allocatable       :: obsLat(:), obsLon(:), obsImpP(:), obsLocR(:), obsGeoid(:)

  write(err_msg,*) "TRACE: ufo_gnssro_bndropp1d_simobs: begin"
  call fckit_log%debug(err_msg)

! check if nobs is consistent in geovals & hofx
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

  deallocate(obsLat) 
  deallocate(obsLon)
  deallocate(obsImpP)
  deallocate(obsLocR)
  deallocate(obsGeoid)

  write(err_msg,*) "TRACE: ufo_gnssro_bndropp1d_simobs: complete"
  call fckit_log%debug(err_msg)

end subroutine ufo_gnssro_bndropp1d_simobs
! ------------------------------------------------------------------------------

end module ufo_gnssro_bndropp1d_mod
