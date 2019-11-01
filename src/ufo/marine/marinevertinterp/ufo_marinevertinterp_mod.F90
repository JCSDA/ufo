! (C) Copyright 2017-2019 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

!> Fortran marinevertinterp module for observation operator

module ufo_marinevertinterp_mod

 use ufo_vars_mod
 use oops_variables_mod

 implicit none
 private

 integer, parameter :: max_string=800

!> Fortran derived type for the observation type
 type, public :: ufo_marinevertinterp    
   type(oops_variables), public :: geovars
   type(oops_variables), public :: obsvars
 contains
   procedure :: setup  => ufo_marinevertinterp_setup
   procedure :: simobs => ufo_marinevertinterp_simobs
 end type ufo_marinevertinterp

contains

! ------------------------------------------------------------------------------
subroutine ufo_marinevertinterp_setup(self)
implicit none
class(ufo_marinevertinterp), intent(inout) :: self
character(max_string)  :: err_msg

integer :: ivar, nvars

nvars = self%obsvars%nvars()
if (nvars /= 1) then
  write(err_msg,*) 'ufo_marinevertinterp_setup error: only variables size 1 supported!'
  call abor1_ftn(err_msg)
endif

! Set variables requested from the model
do ivar = 1, nvars
  call self%geovars%push_back(self%obsvars%variable(ivar))
enddo
call self%geovars%push_back("sea_water_cell_thickness")

end subroutine ufo_marinevertinterp_setup

! ------------------------------------------------------------------------------
subroutine ufo_marinevertinterp_simobs(self, geovals, hofx, obss)
use gsw_pot_to_insitu
use vert_interp_mod
use ufo_tpsp2ti_mod
use iso_c_binding
use kinds
use ufo_geovals_mod, only: ufo_geovals, ufo_geoval, ufo_geovals_get_var
use obsspace_mod
implicit none
class(ufo_marinevertinterp), intent(in)    :: self
type(ufo_geovals),  intent(in)    :: geovals
real(c_double),     intent(inout) :: hofx(:)
type(c_ptr), value, intent(in)    :: obss

    character(len=*), parameter :: myname_="ufo_marinevertinterp_simobs"
    character(max_string)  :: err_msg

    integer :: iobs, ilev, nlev, nlocs
    type(ufo_geoval), pointer :: var, h
    real (kind_real), allocatable :: depth(:,:)
    real(kind_real) :: deptho
    real(kind_real), allocatable :: obs_depth(:)
    integer :: obss_nlocs
    real(kind_real) :: wf, sp, prs
    integer :: wi
    
    ! check if nlocs is consistent in geovals & hofx
    if (geovals%nlocs /= size(hofx,1)) then
       write(err_msg,*) myname_, ' error: nlocs inconsistent!'
       call abor1_ftn(err_msg)
    endif

    ! Associate geoval pointers
    call ufo_geovals_get_var(geovals, self%geovars%variable(1), var)
    call ufo_geovals_get_var(geovals, var_ocn_lay_thick, h)

    ! Read in obs data
    obss_nlocs = obsspace_get_nlocs(obss)
    allocate(obs_depth(obss_nlocs))
    call obsspace_get_db(obss, "MetaData", "depth", obs_depth)

    nlev = var%nval
    nlocs = var%nlocs        
    allocate(depth(nlev,nlocs))
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
