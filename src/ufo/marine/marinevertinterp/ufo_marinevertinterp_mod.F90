! (C) Copyright 2017-2019 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

!> Fortran marinevertinterp module for observation operator

module ufo_marinevertinterp_mod

 use iso_c_binding
 use config_mod
 use kinds
 use ufo_geovals_mod, only: ufo_geovals, ufo_geoval, ufo_geovals_get_var
 use ufo_vars_mod
 use obsspace_mod

 implicit none
 private

 integer, parameter :: max_string=800

!> Fortran derived type for the observation type
 type, public :: ufo_marinevertinterp    
    private
    character(len=max_string), public, allocatable :: varin(:)
    character(len=max_string), public, allocatable :: varout(:)    
 contains
   procedure :: setup  => ufo_marinevertinterp_setup
   procedure :: delete => ufo_marinevertinterp_delete
   procedure :: simobs => ufo_marinevertinterp_simobs
 end type ufo_marinevertinterp

contains

! ------------------------------------------------------------------------------
subroutine ufo_marinevertinterp_setup(self, c_conf)
implicit none
class(ufo_marinevertinterp), intent(inout) :: self
type(c_ptr),        intent(in)    :: c_conf

! Get output variable name (hard-coded to 1)
allocate(self%varout(1))
self%varout = config_get_string_vector(c_conf, max_string, "variable")

! Set input variable names (hard-coded to 2)
allocate(self%varin(2))
self%varin(1) = self%varout(1)
self%varin(2) = "sea_water_cell_thickness"


end subroutine ufo_marinevertinterp_setup

! ------------------------------------------------------------------------------
subroutine ufo_marinevertinterp_delete(self)
implicit none
class(ufo_marinevertinterp), intent(inout) :: self

end subroutine ufo_marinevertinterp_delete

! ------------------------------------------------------------------------------
subroutine ufo_marinevertinterp_simobs(self, geovals, hofx, obss)
use gsw_pot_to_insitu
use vert_interp_mod
use ufo_tpsp2ti_mod
use ufo_marine_ncutils
implicit none
class(ufo_marinevertinterp), intent(in)    :: self
type(ufo_geovals),  intent(in)    :: geovals
real(c_double),     intent(inout) :: hofx(:)
type(c_ptr), value, intent(in)    :: obss

    character(len=*), parameter :: myname_="ufo_marinevertinterp_simobs"
    character(max_string)  :: err_msg

    integer :: iobs, ilev, nlev, nobs
    type(ufo_geoval), pointer :: var, h
    real (kind_real), allocatable :: depth(:,:)
    real(kind_real) :: deptho
    real(kind_real), allocatable :: obs_depth(:)
    integer :: obss_nlocs
    real(kind_real) :: wf, sp, prs
    integer :: wi
    
    ! check if nobs is consistent in geovals & hofx
    if (geovals%nobs /= size(hofx,1)) then
       write(err_msg,*) myname_, ' error: nobs inconsistent!'
       call abor1_ftn(err_msg)
    endif

    ! Associate geoval pointers
    call ufo_geovals_get_var(geovals, self%varin(1), var)
    call ufo_geovals_get_var(geovals, var_ocn_lay_thick, h)

    ! Read in obs data
    obss_nlocs = obsspace_get_nlocs(obss)
    allocate(obs_depth(obss_nlocs))
    call obsspace_get_db(obss, "MetaData", "depth", obs_depth)

    nlev = var%nval
    nobs = var%nobs        
    allocate(depth(nlev,nobs))
    do iobs = 1,size(hofx,1)
       !< Depth from layer thickness
       depth(1,iobs)=0.5*h%vals(1,iobs)
       do ilev = 2, nlev
          depth(ilev,iobs)=sum(h%vals(1:ilev-1,iobs))+0.5*h%vals(ilev,iobs)
       end do          
    end do

    hofx = 0.0
    ! Vertical interpolation
    do iobs = 1,size(hofx,1)

       deptho = obs_depth(iobs)
    
       !< Interpolation weight
       call vert_interp_weights(nlev, deptho, depth(:,iobs), wi, wf)

       !Apply vertical interpolation
       call vert_interp_apply(nlev, var%vals(:,iobs), hofx(iobs), wi, wf)

    enddo

    deallocate(depth)
    deallocate(obs_depth)
    
  end subroutine ufo_marinevertinterp_simobs

end module ufo_marinevertinterp_mod
