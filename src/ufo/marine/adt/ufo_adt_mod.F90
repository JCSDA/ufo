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
  use ioda_obsdb_mod
  use fckit_mpi_module, only: fckit_mpi_comm, fckit_mpi_sum
  
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
    type(c_ptr), value, intent(in) :: obss
    real(c_double),  intent(inout) :: hofx(:)

    character(len=*), parameter :: myname_="ufo_adt_simobs"
    character(max_string) :: err_msg
    type(ufo_geoval), pointer :: geoval_adt
    real(kind_real), allocatable :: obs_adt(:)
    integer :: nobs, obss_nobs
    integer :: iobs, cnt, cnt_glb    
    real(kind_real) :: offset_hofx, pe_offset_hofx
    real(kind_real) :: offset_obs, pe_offset_obs
    type(fckit_mpi_comm) :: f_comm
    real(kind_real) :: r_miss_val = 9999.9
    
    f_comm = fckit_mpi_comm()
    
    ! check if nobs is consistent in geovals & hofx
    nobs = obsspace_get_nlocs(obss)

    !nobs = size(hofx,1)
    if (geovals%nobs /= size(hofx,1)) then
       write(err_msg,*) myname_, ' error: nobs inconsistent!'
       call abor1_ftn(err_msg)
    endif

    ! check if adt variable is in geovals and get it
    call ufo_geovals_get_var(geovals, var_abs_topo, geoval_adt)

    ! Read in obs data
    obss_nobs = obsspace_get_nobs(obss)
    allocate(obs_adt(obss_nobs))
    
    call obsspace_get_db(obss, "ObsValue", "obs_absolute_dynamic_topography", obs_adt)

    ! Local offset
    pe_offset_hofx = 0.0
    pe_offset_obs = 0.0    
    cnt = 0
    do iobs = 1, nobs
       if (abs(hofx(iobs)).lt.r_miss_val) then      
          pe_offset_hofx = pe_offset_hofx + geoval_adt%vals(1,iobs)
          pe_offset_obs = pe_offset_obs + obs_adt(iobs)          
          cnt = cnt + 1
       end if
    end do

    ! Global offsets
    call f_comm%allreduce(pe_offset_hofx, offset_hofx, fckit_mpi_sum())
    call f_comm%allreduce(pe_offset_obs, offset_obs, fckit_mpi_sum())     
    call f_comm%allreduce(cnt, cnt_glb, fckit_mpi_sum()) 
    offset_hofx = offset_hofx/cnt_glb
    offset_obs = offset_obs/cnt_glb    
    
    ! Adjust simulated obs to obs offset
    do iobs = 1, nobs
       hofx(iobs) = geoval_adt%vals(1,iobs) + (offset_obs-offset_hofx)
    enddo

    deallocate(obs_adt)
    
  end subroutine ufo_adt_simobs

end module ufo_adt_mod
