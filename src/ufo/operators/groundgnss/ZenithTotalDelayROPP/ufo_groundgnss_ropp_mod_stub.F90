! (C) Copyright 2017-2018 UCAR
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 

!>Stubbed Fortran module for ground based gnss ropp1d forward operator
!> following the ROPP (2018 Aug) implementation

module ufo_groundgnss_ropp_mod

use iso_c_binding
use kinds
use ufo_vars_mod
use ufo_geovals_mod
use ufo_geovals_mod_c, only: ufo_geovals_registry
use ufo_basis_mod,     only: ufo_basis
use obsspace_mod   
use missing_values_mod
use fckit_log_module,  only : fckit_log

implicit none
public             :: ufo_groundgnss_ropp
private

  !> Fortran derived type for ground based gnss trajectory
type, extends(ufo_basis) :: ufo_groundgnss_ROPP
  contains
    procedure :: simobs    => ufo_groundgnss_ropp_simobs
end type ufo_groundgnss_ROPP

contains

! ------------------------------------------------------------------------------
! ------------------------------------------------------------------------------
subroutine ufo_groundgnss_ropp_simobs(self, geovals, hofx, obss)

  implicit none
  class(ufo_groundgnss_ROPP),  intent(in)    :: self
  type(ufo_geovals),           intent(in)    :: geovals
  real(kind_real),             intent(inout) :: hofx(:)
  type(c_ptr), value,          intent(in)    :: obss
  real(c_double)                             :: missingDouble

  character(len=*), parameter  :: myname_="ufo_groundgnss_ropp_simobs"
  integer, parameter           :: max_string = 800

  character(max_string)              :: err_msg

  write(err_msg,*) "TRACE: ufo_groundgnss_ropp_simobs_stub: begin"
  call fckit_log%info(err_msg)
  write(err_msg,*) "WARNING: GroundgnssROPP operator cannot run when ROPP code is not available"
  call fckit_log%info(err_msg)

! check if nobs is consistent in geovals & hofx
  if (geovals%nlocs /= size(hofx)) then
      write(err_msg,*) myname_, ' error: nlocs inconsistent!'
      call abor1_ftn(err_msg)
  endif

  missingDouble = missing_value(missingDouble)

! initialize HofX to missing
  hofx(:) = missingDouble

  write(err_msg,*) "TRACE: ufo_groundgnss_ropp_simobs_stub: completed"
  call fckit_log%info(err_msg)

end subroutine ufo_groundgnss_ropp_simobs
! ------------------------------------------------------------------------------

end module ufo_groundgnss_ropp_mod
