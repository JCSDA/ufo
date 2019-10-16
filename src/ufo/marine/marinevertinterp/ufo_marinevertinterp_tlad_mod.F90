! (C) Copyright 2017-2019 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

!> Fortran marinevertinterp module for tl/ad observation operator

module ufo_marinevertinterp_tlad_mod

 use oops_variables_mod
 use kinds
 use ufo_geovals_mod, only: ufo_geovals, ufo_geoval, ufo_geovals_get_var
 use ufo_vars_mod

 implicit none
 private

 integer, parameter :: max_string=800

 !> Fortran derived type for the tl/ad observation operator
 type, public :: ufo_marinevertinterp_tlad    
    private
    type(oops_variables), public :: geovars
    type(oops_variables), public :: obsvars
    logical,                    public :: ltraj = .false. !< trajectory set?    
    integer                            :: nlocs       !< Number of observations
    integer                            :: nval       !< Number of level in model's profiles 
    type(ufo_geoval)                   :: var        !< traj
    type(ufo_geoval)                   :: h          !< Layer thickness (traj) ] obs locations
    real (kind=kind_real), allocatable :: depth(:,:) !< Depth                     [nval x nlocs]     
    real (kind=kind_real), allocatable :: deptho(:)  !< Observation location
    real(kind_real), allocatable       :: wf(:)      !< Vertical interpolation weights
    integer, allocatable               :: wi(:)      !< Vertical interpolation indices
 contains
  procedure :: setup  => ufo_marinevertinterp_tlad_setup
  procedure :: delete  => ufo_marinevertinterp_tlad_delete
  procedure :: settraj => ufo_marinevertinterp_tlad_settraj
  procedure :: simobs_tl  => ufo_marinevertinterp_simobs_tl
  procedure :: simobs_ad  => ufo_marinevertinterp_simobs_ad
 end type ufo_marinevertinterp_tlad

contains

! ------------------------------------------------------------------------------
subroutine ufo_marinevertinterp_tlad_setup(self)
implicit none
class(ufo_marinevertinterp_tlad), intent(inout) :: self

integer :: ivar, nvars

nvars = self%obsvars%nvars()
do ivar = 1, nvars
  call self%geovars%push_back(self%obsvars%variable(ivar))
enddo

end subroutine ufo_marinevertinterp_tlad_setup

! ------------------------------------------------------------------------------
subroutine ufo_marinevertinterp_tlad_delete(self)
implicit none
class(ufo_marinevertinterp_tlad), intent(inout) :: self

    if (allocated(self%wi)) deallocate(self%wi)
    if (allocated(self%wf)) deallocate(self%wf)
    if (allocated(self%deptho)) deallocate(self%deptho)
    if (allocated(self%depth)) deallocate(self%depth)
    if (allocated(self%var%vals)) deallocate(self%var%vals)
    if (allocated(self%h%vals)) deallocate(self%h%vals)
    self%ltraj = .false.

end subroutine ufo_marinevertinterp_tlad_delete

! ------------------------------------------------------------------------------
subroutine ufo_marinevertinterp_tlad_settraj(self, geovals, obss)
    use vert_interp_mod
    use ufo_tpsp2ti_mod
    use obsspace_mod
    use iso_c_binding
    implicit none
    class(ufo_marinevertinterp_tlad), intent(inout)  :: self    !< Complete trajectory needed by the operator
    type(ufo_geovals), intent(in)                    :: geovals !< Model background
    type(c_ptr), value, intent(in)                   :: obss    !< 

    character(len=*), parameter :: myname_="ufo_marinevertinterp_tlad_settraj"
    character(max_string) :: err_msg

    type(ufo_geoval), pointer :: var, h
    integer :: nlocs, nlev, iobs, ilev

    real(kind_real), allocatable :: obs_depth(:)
    integer :: obss_nlocs

    ! Associate geovals pointers
    call ufo_geovals_get_var(geovals, self%geovars%variable(1), var)
    call ufo_geovals_get_var(geovals, var_ocn_lay_thick, h)

    call ufo_marinevertinterp_tlad_delete(self)

    nlocs = h%nlocs
    nlev = h%nval
    
    self%nlocs = nlocs
    self%nval = nlev
    
    self%var = var
    self%h    = h

    allocate(self%deptho(nlocs))

    obss_nlocs = obsspace_get_nlocs(obss)
    allocate(obs_depth(obss_nlocs))
    call obsspace_get_db(obss, "MetaData", "depth", obs_depth)

    self%deptho = obs_depth

    !< Depth from layer thickness
    allocate(self%depth(nlev,nlocs))
    do iobs = 1, nlocs
       self%depth(1,iobs)=0.5*self%h%vals(1,iobs)       
       do ilev = 2, nlev
          self%depth(ilev,iobs)=sum(self%h%vals(1:ilev-1,iobs))+0.5*self%h%vals(ilev,iobs)
       end do
    end do

    !< Interpolation weight
    allocate(self%wi(nlocs),self%wf(nlocs))
    do iobs = 1, nlocs    
       call vert_interp_weights(nlev,self%deptho(iobs),self%depth(:,iobs),self%wi(iobs),self%wf(iobs))
    end do
    self%ltraj    = .true.

    deallocate(obs_depth)

end subroutine ufo_marinevertinterp_tlad_settraj

! ------------------------------------------------------------------------------
subroutine ufo_marinevertinterp_simobs_tl(self, geovals, hofx, obss)
    use ufo_tpsp2ti_mod
    use gsw_pot_to_insitu
    use vert_interp_mod
    use obsspace_mod
    use missing_values_mod
    use iso_c_binding
    implicit none
    class(ufo_marinevertinterp_tlad), intent(in)    :: self
    type(ufo_geovals),       intent(in)    :: geovals
    real(c_double),          intent(inout) :: hofx(:)
    type(c_ptr), value,      intent(in)    :: obss

    character(len=*), parameter :: myname_="ufo_marinevertinterp_simobs_tl"
    character(max_string) :: err_msg
    integer :: iobs, ilev, nlev, nlocs
    type(ufo_geoval), pointer :: var_d !< Increments from geoval

    ! check if trajectory was set
    if (.not. self%ltraj) then
       write(err_msg,*) myname_, ' trajectory wasnt set!'
       call abor1_ftn(err_msg)
    endif

    ! check if nlocs is consistent in geovals & hofx
    if (geovals%nlocs /= size(hofx,1)) then
       write(err_msg,*) myname_, ' error: nlocs inconsistent!'
       call abor1_ftn(err_msg)
    endif

    ! Associate geovals
    call ufo_geovals_get_var(geovals, self%geovars%variable(1), var_d)

    nlev = var_d%nval
    nlocs = var_d%nlocs        

    ! linear vertical interp
    hofx = 0.0
    do iobs = 1,nlocs
       call vert_interp_apply(nlev, var_d%vals(:,iobs), hofx(iobs), self%wi(iobs), self%wf(iobs))
    enddo

end subroutine ufo_marinevertinterp_simobs_tl

! ------------------------------------------------------------------------------
subroutine ufo_marinevertinterp_simobs_ad(self, geovals, hofx, obss)
    use ufo_tpsp2ti_mod
    use gsw_pot_to_insitu
    use vert_interp_mod
    use obsspace_mod
    use missing_values_mod
    use iso_c_binding
    implicit none
    class(ufo_marinevertinterp_tlad), intent(in)    :: self
    type(ufo_geovals),       intent(inout) :: geovals
    real(c_double),          intent(in)    :: hofx(:)
    type(c_ptr), value,      intent(in)    :: obss

    character(len=*), parameter :: myname_="ufo_marinevertinterp_simobs_ad"
    character(max_string) :: err_msg

    real (kind=kind_real) :: deptho !< Observation location
        
    integer :: iobs, nlocs, ilev, nlev
    type(ufo_geoval), pointer :: dvar
    real(c_double) :: missing

    !> Set missing value
    missing = missing_value(missing)
    
    ! check if trajectory was set
    if (.not. self%ltraj) then
       write(err_msg,*) myname_, ' trajectory wasnt set!'
       call abor1_ftn(err_msg)
    endif

    ! check if nlocs is consistent in geovals & hofx
    if (geovals%nlocs /= size(hofx,1)) then
       write(err_msg,*) myname_, ' error: nlocs inconsistent!'
       call abor1_ftn(err_msg)
    endif

    if (.not. geovals%linit ) geovals%linit=.true.
    
    ! Associate geovals
    call ufo_geovals_get_var(geovals, var_ocn_salt, dvar)
    
    nlev = self%nval
    nlocs = self%nlocs
    
    if (.not. allocated(dvar%vals)) allocate(dvar%vals(nlev, size(hofx,1)))

    ! backward vertical interp
    dvar%vals = 0.0
    do iobs = 1, size(hofx,1)
       if (hofx(iobs) /= missing) then
          deptho = self%deptho(iobs)
      
          ! Adjoint obs operator
          call vert_interp_apply_ad(nlev, dvar%vals(:,iobs), hofx(iobs), self%wi(iobs), self%wf(iobs))
       end if
    enddo

end subroutine ufo_marinevertinterp_simobs_ad

end module ufo_marinevertinterp_tlad_mod
