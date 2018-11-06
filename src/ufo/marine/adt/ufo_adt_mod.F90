! (C) Copyright 2017 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

!> Fortran module to handle adt observations

module ufo_adt_mod

  use iso_c_binding
  use ufo_vars_mod
  use ioda_locs_mod
  use ufo_geovals_mod
  use kinds
  use ncd_kinds, only:  i_kind,r_single,r_kind,r_double
  use obsspace_mod 

  implicit none
  public :: ufo_adt
  public :: ufo_adt_simobs
  private
  integer, parameter :: max_string=800

  !> Fortran derived type for adt observation operator
  type :: ufo_adt
  end type ufo_adt

  ! ------------------------------------------------------------------------------

contains

  ! ------------------------------------------------------------------------------

  subroutine ufo_adt_simobs(self, geovals, hofx, obss)

    use ufo_marine_ncutils
    use iso_c_binding
    
    implicit none
    
    type(ufo_adt),      intent(in) :: self
    type(ufo_geovals),  intent(in) :: geovals
    type(c_ptr), value, intent(in) :: obss     !< adt observations
    real(c_double),  intent(inout) :: hofx(:)

    character(len=*), parameter :: myname_="ufo_adt_simobs"
    character(max_string) :: err_msg

    integer :: iobs
    real(kind_real) :: offset_obs, offset_hofx
    type(ufo_geoval), pointer :: geoval_adt

    real(kind_real), allocatable :: obs_adt(:)
    integer :: nobs
    
    ! check if nobs is consistent in geovals & hofx
    nobs = size(hofx,1)
    if (geovals%nobs /= size(hofx,1)) then
       write(err_msg,*) myname_, ' error: nobs inconsistent!'
       call abor1_ftn(err_msg)
    endif

    ! check if adt variable is in geovals and get it
    call ufo_geovals_get_var(geovals, var_abs_topo, geoval_adt)

    hofx = 0.0

    !call obsspace_get_db(obss, "ObsValue", "adt", obs_adt)

    ! Compute offset
    offset_hofx=sum(geoval_adt%vals(1,:))/nobs
    offset_obs=0.0!sum(obs_adt(:))/nobs

    ! adt obs operator
    do iobs = 1, nobs
       ! remove offset from hofx
       hofx(iobs) = geoval_adt%vals(1,iobs)+(offset_obs-offset_hofx)
    enddo

    !deallocate(obs_adt)
    
  end subroutine ufo_adt_simobs

end module ufo_adt_mod
