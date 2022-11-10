! (C) Copyright 2017-2018 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

!> Fortran adt module for observation operator

module ufo_adt_mod

 use fckit_configuration_module, only: fckit_configuration
 use iso_c_binding
 use kinds

 use ufo_geovals_mod, only: ufo_geovals, ufo_geoval, ufo_geovals_get_var
 use ufo_basis_mod, only: ufo_basis
 use ufo_vars_mod
 use obsspace_mod
 use missing_values_mod

 implicit none
 private

 integer, parameter :: max_string=800

!> Fortran derived type for the observation type
 type, extends(ufo_basis), public :: ufo_adt
 private
 contains
   procedure :: setup  => ufo_adt_setup
   procedure :: delete => ufo_adt_delete
   procedure :: simobs => ufo_adt_simobs
 end type ufo_adt

contains

! ------------------------------------------------------------------------------
subroutine ufo_adt_setup(self, f_conf)
implicit none
class(ufo_adt), intent(inout)         :: self
type(fckit_configuration), intent(in) :: f_conf

end subroutine ufo_adt_setup

! ------------------------------------------------------------------------------
subroutine ufo_adt_delete(self)
implicit none
class(ufo_adt), intent(inout) :: self

end subroutine ufo_adt_delete

! ------------------------------------------------------------------------------
subroutine ufo_adt_simobs(self, geovals, hofx, obss)
use fckit_mpi_module, only: fckit_mpi_comm, fckit_mpi_sum
implicit none
    class(ufo_adt), intent(in)    :: self
    type(ufo_geovals),  intent(in)    :: geovals
    real(c_double),     intent(inout) :: hofx(:)
    type(c_ptr), value, intent(in)    :: obss

    character(len=*), parameter :: myname_="ufo_adt_simobs"
    character(max_string) :: err_msg
    type(ufo_geoval), pointer :: geoval_adt
    real(kind_real), allocatable :: obs_adt(:)
    integer :: obss_nlocs
    integer :: iobs, cnt, cnt_glb
    real(kind_real) :: offset, pe_offset
    type(fckit_mpi_comm) :: f_comm
    real(c_double) :: missing

    call obsspace_get_comm(obss, f_comm)

    ! Set missing flag
    missing = missing_value(missing)

    ! check if nlocs is consistent in geovals & hofx
    obss_nlocs = obsspace_get_nlocs(obss)

    !nlocs = size(hofx,1)
    if (geovals%nlocs /= size(hofx,1)) then
       write(err_msg,*) myname_, ' error: nlocs inconsistent!'
       call abor1_ftn(err_msg)
    endif

    ! check if adt variable is in geovals and get it
    call ufo_geovals_get_var(geovals, var_abs_topo, geoval_adt)

    ! Read in obs data
    allocate(obs_adt(obss_nlocs))

    call obsspace_get_db(obss, "ObsValue", "absolute_dynamic_topography", obs_adt)

    ! Local offset
    pe_offset = 0.0
    cnt = 0
    do iobs = 1, obss_nlocs
      ! TODO this really should only be done on qc'd obs, otherwise bad obs
      ! that aren't assimilated will throw off this global offset
       if (geoval_adt%vals(1,iobs) /= missing .and. obs_adt(iobs) /= missing) then
          pe_offset = pe_offset + obs_adt(iobs) - geoval_adt%vals(1,iobs)
          cnt = cnt + 1
       end if
    end do

    ! Global offsets
    call f_comm%allreduce(pe_offset, offset, fckit_mpi_sum())
    call f_comm%allreduce(cnt, cnt_glb, fckit_mpi_sum())
    if (cnt_glb > 0) then
      offset = offset/cnt_glb
    end if

    ! Adjust simulated obs to obs offset
    do iobs = 1, obss_nlocs
      hofx(iobs) = geoval_adt%vals(1,iobs)
      ! only add offset if hofx is not a missing value
      if (hofx(iobs) /= missing) then
          hofx(iobs) = hofx(iobs) + offset
      end if
    enddo

    deallocate(obs_adt)

  end subroutine ufo_adt_simobs

end module ufo_adt_mod
