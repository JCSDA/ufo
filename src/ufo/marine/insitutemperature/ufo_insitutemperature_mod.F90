! (C) Copyright 2017-2018 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

!> Fortran insitutemperature module for observation operator

module ufo_insitutemperature_mod

 use iso_c_binding
 use config_mod
 use kinds

 use ufo_geovals_mod, only: ufo_geovals, ufo_geoval, ufo_geovals_get_var
 use ufo_basis_mod, only: ufo_basis
 use ufo_vars_mod
 use obsspace_mod

 implicit none
 private

 integer, parameter :: max_string=800

!> Fortran derived type for the observation type
 type, extends(ufo_basis), public :: ufo_insitutemperature
 private
 contains
   procedure :: setup  => ufo_insitutemperature_setup
   procedure :: delete => ufo_insitutemperature_delete
   procedure :: simobs => ufo_insitutemperature_simobs
 end type ufo_insitutemperature

contains

! ------------------------------------------------------------------------------
subroutine ufo_insitutemperature_setup(self, c_conf)
implicit none
class(ufo_insitutemperature), intent(inout) :: self
type(c_ptr),        intent(in)    :: c_conf

end subroutine ufo_insitutemperature_setup

! ------------------------------------------------------------------------------
subroutine ufo_insitutemperature_delete(self)
implicit none
class(ufo_insitutemperature), intent(inout) :: self

end subroutine ufo_insitutemperature_delete

! ------------------------------------------------------------------------------
subroutine ufo_insitutemperature_simobs(self, geovals, hofx, obss)
use gsw_pot_to_insitu
use vert_interp_mod
use ufo_tpsp2ti_mod
use ufo_marine_ncutils
implicit none
class(ufo_insitutemperature), intent(in)    :: self
type(ufo_geovals),  intent(in)    :: geovals
real(c_double),     intent(inout) :: hofx(:)
type(c_ptr), value, intent(in)    :: obss

    character(len=*), parameter :: myname_="ufo_insitutemperature_simobs"
    character(max_string)  :: err_msg

    integer :: iobs, ilev, nlev, nlocs
    type(ufo_geoval), pointer :: temp, salt, h
    real (kind_real), allocatable :: depth(:,:)
    real(kind_real) :: lono, lato, deptho
    real(kind_real), allocatable :: obs_lon(:)
    real(kind_real), allocatable :: obs_lat(:)
    real(kind_real), allocatable :: obs_depth(:)
    integer :: obss_nlocs
    real(kind_real) :: wf, tp, sp, prs
    integer :: wi
    
    ! check if nlocs is consistent in geovals & hofx
    if (geovals%nlocs /= size(hofx,1)) then
       write(err_msg,*) myname_, ' error: nlocs inconsistent!'
       call abor1_ftn(err_msg)
    endif

    ! Associate geoval pointers
    call ufo_geovals_get_var(geovals, var_ocn_pot_temp, temp)
    call ufo_geovals_get_var(geovals, var_ocn_salt, salt)
    call ufo_geovals_get_var(geovals, var_ocn_lay_thick, h)

    ! Read in obs data
    obss_nlocs = obsspace_get_nlocs(obss)
    allocate(obs_lon(obss_nlocs))
    allocate(obs_lat(obss_nlocs))
    allocate(obs_depth(obss_nlocs))
    call obsspace_get_db(obss, "MetaData", "longitude", obs_lon)
    call obsspace_get_db(obss, "MetaData", "latitude", obs_lat)
    call obsspace_get_db(obss, "MetaData", "depth", obs_depth)

    nlev = temp%nval
    nlocs = temp%nlocs        
    allocate(depth(nlev,nlocs))
    do iobs = 1,size(hofx,1)
       !< Depth from layer thickness
       depth(1,iobs)=0.5*h%vals(1,iobs)
       do ilev = 2, nlev
          depth(ilev,iobs)=sum(h%vals(1:ilev-1,iobs))+0.5*h%vals(ilev,iobs)
       end do          
    end do

    ! Information for temporary output file
    
    hofx = 0.0
    ! insitu temperature profile obs operator
    do iobs = 1,size(hofx,1)

       lono = obs_lon(iobs)
       lato = obs_lat(iobs)
       deptho = obs_depth(iobs)
    
       !< Interpolation weight
       call vert_interp_weights(nlev, deptho, depth(:,iobs), wi, wf)
       if (deptho.ge.maxval(depth)) then
          wi=nlev-1
          wf=0.0
       end if

       ! Interpolate temp_p, salt_p to deptho
       call vert_interp_apply(nlev, temp%vals(:,iobs), tp, wi, wf)
       call vert_interp_apply(nlev, salt%vals(:,iobs), sp, wi, wf)

       ! Get insitu temp at model levels and obs location (lono, lato, zo)
       call insitu_t_nl(hofx(iobs), tp, sp, lono, lato, deptho)

    enddo

    deallocate(depth)
    deallocate(obs_lon)
    deallocate(obs_lat)
    deallocate(obs_depth)
    
  end subroutine ufo_insitutemperature_simobs

end module ufo_insitutemperature_mod
