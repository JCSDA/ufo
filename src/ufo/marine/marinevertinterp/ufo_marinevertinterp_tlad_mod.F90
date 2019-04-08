! (C) Copyright 2017-2019 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

!> Fortran marinevertinterp module for tl/ad observation operator

module ufo_marinevertinterp_tlad_mod

 use iso_c_binding
 use config_mod
 use kinds

 use ufo_geovals_mod, only: ufo_geovals, ufo_geoval, ufo_geovals_get_var
 use ufo_vars_mod
 use obsspace_mod
 use missing_values_mod

 implicit none
 private

 integer, parameter :: max_string=800

 !> Fortran derived type for the tl/ad observation operator
 type, public :: ufo_marinevertinterp_tlad    
    private
    logical,                    public :: ltraj = .false. !< trajectory set?    
    character(len=max_string), public, allocatable :: varin(:)
    character(len=max_string), public, allocatable :: varout(:)    

    integer                            :: nobs       !< Number of observations
    integer                            :: nval       !< Number of level in model's profiles 
    type(ufo_geoval)                   :: var        !< traj
    type(ufo_geoval)                   :: h          !< Layer thickness (traj) ] obs locations
    real (kind=kind_real), allocatable :: depth(:,:) !< Depth                     [nval x nobs]     
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
subroutine ufo_marinevertinterp_tlad_setup(self, c_conf)
implicit none
class(ufo_marinevertinterp_tlad), intent(inout) :: self
type(c_ptr),              intent(in)    :: c_conf

! Get output variable name (hard-coded to 1)
allocate(self%varout(1))
self%varout = config_get_string_vector(c_conf, max_string, "variable")

! Set input variable names (hard-coded to 2)
allocate(self%varin(2))
self%varin(1) = self%varout(1)
self%varin(2) = "sea_water_cell_thickness"

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
        
    implicit none
    class(ufo_marinevertinterp_tlad), intent(inout)  :: self    !< Complete trajectory needed by the operator
    type(ufo_geovals), intent(in)                    :: geovals !< Model background
    type(c_ptr), value, intent(in)                   :: obss    !< 

    character(len=*), parameter :: myname_="ufo_marinevertinterp_tlad_settraj"
    character(max_string) :: err_msg

    type(ufo_geoval), pointer :: var, h
    integer :: nobs, nlev, iobs, ilev

    real(kind_real), allocatable :: obs_depth(:)
    integer :: obss_nlocs

    ! Associate geovals pointers
    call ufo_geovals_get_var(geovals, self%varin(1), var)
    call ufo_geovals_get_var(geovals, var_ocn_lay_thick, h)

    call ufo_marinevertinterp_tlad_delete(self)

    nobs = h%nobs
    nlev = h%nval
    
    self%nobs = nobs
    self%nval = nlev
    
    self%var = var
    self%h    = h

    allocate(self%deptho(nobs))

    obss_nlocs = obsspace_get_nlocs(obss)
    allocate(obs_depth(obss_nlocs))
    call obsspace_get_db(obss, "MetaData", "depth", obs_depth)

    self%deptho = obs_depth

    !< Depth from layer thickness
    allocate(self%depth(nlev,nobs))
    do iobs = 1, nobs
       self%depth(1,iobs)=0.5*self%h%vals(1,iobs)       
       do ilev = 2, nlev
          self%depth(ilev,iobs)=sum(self%h%vals(1:ilev-1,iobs))+0.5*self%h%vals(ilev,iobs)
       end do
    end do

    !< Interpolation weight
    allocate(self%wi(nobs),self%wf(nobs))
    do iobs = 1, nobs    
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
    implicit none
    class(ufo_marinevertinterp_tlad), intent(in)    :: self
    type(ufo_geovals),       intent(in)    :: geovals
    real(c_double),          intent(inout) :: hofx(:)
    type(c_ptr), value,      intent(in)    :: obss

    character(len=*), parameter :: myname_="ufo_marinevertinterp_simobs_tl"
    character(max_string) :: err_msg

    integer :: iobs, ilev, nlev, nobs

    type(ufo_geoval), pointer :: var_d, dlayerthick !< Increments from geovals

    ! check if trajectory was set
    if (.not. self%ltraj) then
       write(err_msg,*) myname_, ' trajectory wasnt set!'
       call abor1_ftn(err_msg)
    endif

    ! check if nobs is consistent in geovals & hofx
    if (geovals%nobs /= size(hofx,1)) then
       write(err_msg,*) myname_, ' error: nobs inconsistent!'
       call abor1_ftn(err_msg)
    endif

    ! Associate geovals
    call ufo_geovals_get_var(geovals, self%varin(1), var_d)
    call ufo_geovals_get_var(geovals, var_ocn_lay_thick, dlayerthick)

    ! Make sure thickness is not perturbed
    dlayerthick%vals=0.0
    
    nlev = var_d%nval
    nobs = var_d%nobs        

    ! linear vertical interp
    hofx = 0.0
    do iobs = 1,nobs
       call vert_interp_apply(nlev, var_d%vals(:,iobs), hofx(iobs), self%wi(iobs), self%wf(iobs))
    enddo

end subroutine ufo_marinevertinterp_simobs_tl

! ------------------------------------------------------------------------------
subroutine ufo_marinevertinterp_simobs_ad(self, geovals, hofx, obss)
    use ufo_tpsp2ti_mod
    use gsw_pot_to_insitu
    use vert_interp_mod
    implicit none
    class(ufo_marinevertinterp_tlad), intent(in)    :: self
    type(ufo_geovals),       intent(inout) :: geovals
    real(c_double),          intent(in)    :: hofx(:)
    type(c_ptr), value,      intent(in)    :: obss

    character(len=*), parameter :: myname_="ufo_marinevertinterp_simobs_ad"
    character(max_string) :: err_msg

    real (kind=kind_real) :: deptho !< Observation location
        
    integer :: iobs, nobs, ilev, nlev
    type(ufo_geoval), pointer :: dvar, dlayerthick
    real(c_double) :: missing

    !> Set missing value
    missing = missing_value(missing)
    
    ! check if trajectory was set
    if (.not. self%ltraj) then
       write(err_msg,*) myname_, ' trajectory wasnt set!'
       call abor1_ftn(err_msg)
    endif

    ! check if nobs is consistent in geovals & hofx
    if (geovals%nobs /= size(hofx,1)) then
       write(err_msg,*) myname_, ' error: nobs inconsistent!'
       call abor1_ftn(err_msg)
    endif

    if (.not. geovals%linit ) geovals%linit=.true.
    
    ! Associate geovals
    call ufo_geovals_get_var(geovals, var_ocn_salt, dvar)
    call ufo_geovals_get_var(geovals, var_ocn_lay_thick, dlayerthick)
    
    nlev = self%nval
    nobs = self%nobs
    
    if (.not. allocated(dvar%vals)) allocate(dvar%vals(nlev, size(hofx,1)))
    if (.not. allocated(dlayerthick%vals)) allocate(dlayerthick%vals(nlev, size(hofx,1)))    

    ! backward vertical interp
    dvar%vals = 0.0
    do iobs = 1, size(hofx,1)

       if (hofx(iobs) /= missing) then

          deptho = self%deptho(iobs)
      
          ! Adjoint obs operator
          call vert_interp_apply_ad(nlev, dvar%vals(:,iobs), hofx(iobs), self%wi(iobs), self%wf(iobs))

          ! Layer thickness is not a control variable: zero it out!
          dlayerthick%vals=0.0
       end if
    enddo

end subroutine ufo_marinevertinterp_simobs_ad

end module ufo_marinevertinterp_tlad_mod
