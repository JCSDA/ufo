! (C) Copyright 2017-2018 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

!> Fortran insitutemperature module for tl/ad observation operator

module ufo_insitutemperature_tlad_mod

 use iso_c_binding
 use config_mod
 use kinds

 use ufo_geovals_mod, only: ufo_geovals, ufo_geoval, ufo_geovals_get_var
 use ufo_basis_tlad_mod, only: ufo_basis_tlad
 use ufo_vars_mod
 use obsspace_mod

 implicit none
 private

 integer, parameter :: max_string=800

 !> Fortran derived type for the tl/ad observation operator
 type, extends(ufo_basis_tlad), public :: ufo_insitutemperature_tlad
 private
     integer                            :: nobs       !< Number of observations
     integer                            :: nval       !< Number of level in model's profiles 
     type(ufo_geoval)                   :: temp       !< Temperature (traj)     ] Model vertical 
     type(ufo_geoval)                   :: salt       !< Salinity (traj)        ] profile at 
     type(ufo_geoval)                   :: h          !< Layer thickness (traj) ] obs locations
     real (kind=kind_real), allocatable :: depth(:,:) !< Depth                     [nval x nobs]     
     real (kind=kind_real), allocatable :: lono(:)    !< Observation location
     real (kind=kind_real), allocatable :: lato(:)    !< Observation location
     real (kind=kind_real), allocatable :: deptho(:)  !< Observation location
     real (kind=kind_real), allocatable :: tempo(:)   !< temp interpolated at observation location
     real (kind=kind_real), allocatable :: salto(:)   !< salt interpolated at observation location     
     real(kind_real), allocatable       :: wf(:)      !< Vertical interpolation weights
     integer, allocatable               :: wi(:)      !< Vertical interpolation indices
     real (kind=kind_real), allocatable :: jac(:,:)   !< Jacobian     [2 x nobs]
 contains
  procedure :: setup  => ufo_insitutemperature_tlad_setup
  procedure :: delete  => ufo_insitutemperature_tlad_delete
  procedure :: settraj => ufo_insitutemperature_tlad_settraj
  procedure :: simobs_tl  => ufo_insitutemperature_simobs_tl
  procedure :: simobs_ad  => ufo_insitutemperature_simobs_ad
 end type ufo_insitutemperature_tlad

contains

! ------------------------------------------------------------------------------
subroutine ufo_insitutemperature_tlad_setup(self, c_conf)
implicit none
class(ufo_insitutemperature_tlad), intent(inout) :: self
type(c_ptr),              intent(in)    :: c_conf

end subroutine ufo_insitutemperature_tlad_setup

! ------------------------------------------------------------------------------
subroutine ufo_insitutemperature_tlad_delete(self)
implicit none
class(ufo_insitutemperature_tlad), intent(inout) :: self

    if (allocated(self%jac)) deallocate(self%jac)
    if (allocated(self%wi)) deallocate(self%wi)
    if (allocated(self%wf)) deallocate(self%wf)
    if (allocated(self%deptho)) deallocate(self%deptho)
    if (allocated(self%lato)) deallocate(self%lato)
    if (allocated(self%lono)) deallocate(self%lono)
    if (allocated(self%depth)) deallocate(self%depth)
    if (allocated(self%temp%vals)) deallocate(self%temp%vals)
    if (allocated(self%salt%vals)) deallocate(self%salt%vals)
    if (allocated(self%h%vals)) deallocate(self%h%vals)
    if (allocated(self%tempo)) deallocate(self%tempo)
    if (allocated(self%salto)) deallocate(self%salto)    
    self%ltraj = .false.

end subroutine ufo_insitutemperature_tlad_delete

! ------------------------------------------------------------------------------
subroutine ufo_insitutemperature_tlad_settraj(self, geovals, obss)
    use vert_interp_mod
    use ufo_tpsp2ti_mod
        
    implicit none
    class(ufo_insitutemperature_tlad), intent(inout)  :: self    !< Complete trajectory needed by the operator
    type(ufo_geovals), intent(in)                    :: geovals !< Model background
    type(c_ptr), value, intent(in)                   :: obss    !< Insitu temperature observations

    character(len=*), parameter :: myname_="ufo_insitutemperature_tlad_settraj"
    character(max_string) :: err_msg

    type(ufo_geoval), pointer :: temp, salt, h
    integer :: nobs, nlev, iobs, ilev

    real(kind_real), allocatable :: obs_lat(:)
    real(kind_real), allocatable :: obs_lon(:)
    real(kind_real), allocatable :: obs_depth(:)
    integer :: obss_nlocs

    ! check if sea temperature profile variable is in geovals and get it
    call ufo_geovals_get_var(geovals, var_ocn_pot_temp, temp)

    ! check if sea salinity profile variable is in geovals and get it
    call ufo_geovals_get_var(geovals, var_ocn_salt, salt)

    ! check if ocean layer thickness variable is in geovals and get it
    call ufo_geovals_get_var(geovals, var_ocn_lay_thick, h)

    call ufo_insitutemperature_tlad_delete(self)

    nobs = h%nobs
    nlev = h%nval
    
    self%nobs = nobs
    self%nval = nlev
    
    self%temp = temp
    self%salt = salt
    self%h    = h

    allocate(self%lono(nobs))
    allocate(self%lato(nobs))
    allocate(self%deptho(nobs))

    obss_nlocs = obsspace_get_nlocs(obss)
    allocate(obs_lat(obss_nlocs))
    allocate(obs_lon(obss_nlocs))
    allocate(obs_depth(obss_nlocs))

    call obsspace_get_db(obss, "MetaData", "longitude", obs_lon)
    call obsspace_get_db(obss, "MetaData", "latitude", obs_lat)
    call obsspace_get_db(obss, "MetaData", "depth", obs_depth)

    self%lono = obs_lon
    self%lato = obs_lat
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
       if (self%deptho(iobs).ge.maxval(self%depth(:,iobs))) then
          self%wi(iobs)=nlev-1
          self%wf(iobs)=0.0
       end if
    end do
    self%ltraj    = .true.

    !< Jacobian
    allocate(self%jac(2,nobs),self%tempo(nobs),self%salto(nobs))
    do iobs = 1, nobs
       ! Interpolate background do obs depth and save in traj
       call vert_interp_apply(nlev, self%temp%vals(:,iobs), self%tempo(iobs), self%wi(iobs), self%wf(iobs))
       call vert_interp_apply(nlev, self%salt%vals(:,iobs), self%salto(iobs), self%wi(iobs), self%wf(iobs))
       
       ! Compute jacobian
       call insitu_t_jac(self%jac(:,iobs), self%tempo(iobs), self%salto(iobs), self%lono(iobs), self%lato(iobs), self%deptho(iobs))
    end do

    deallocate(obs_lat)
    deallocate(obs_lon)
    deallocate(obs_depth)

end subroutine ufo_insitutemperature_tlad_settraj

! ------------------------------------------------------------------------------
subroutine ufo_insitutemperature_simobs_tl(self, geovals, hofx, obss)
    use ufo_tpsp2ti_mod
    use gsw_pot_to_insitu
    use vert_interp_mod
    implicit none
    class(ufo_insitutemperature_tlad), intent(in)    :: self
    type(ufo_geovals),       intent(in)    :: geovals
    real(c_double),          intent(inout) :: hofx(:)
    type(c_ptr), value,      intent(in)    :: obss

    character(len=*), parameter :: myname_="ufo_insitutemperature_simobs_tl"
    character(max_string) :: err_msg

    integer :: iobs, ilev, nlev, nobs

    type(ufo_geoval), pointer :: temp_d, salt_d, dlayerthick !< Increments from geovals
    real (kind=kind_real) :: lono, lato, deptho !< Observation location

    ! Vertical interpolation
    real(kind_real) :: dtp, dsp

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

    ! check if sea temperature profile variable is in geovals and get it
    call ufo_geovals_get_var(geovals, var_ocn_pot_temp, temp_d)

    ! check if sea salinity profile variable is in geovals and get it
    call ufo_geovals_get_var(geovals, var_ocn_salt, salt_d)

    ! check if sea layer thickness variable is in geovals get it and zero it out
    call ufo_geovals_get_var(geovals, var_ocn_lay_thick, dlayerthick)
    ! Make sure thickness is not perturbed
    dlayerthick%vals=0.0
    
    nlev = temp_d%nval
    nobs = temp_d%nobs        

    ! linear sea temperature profile obs operator
    hofx = 0.0
    do iobs = 1,nobs

       lono = self%lono(iobs)
       lato = self%lato(iobs)
       deptho = self%deptho(iobs)

       !  Interpolate increment
       call vert_interp_apply(nlev, temp_d%vals(:,iobs), dtp, self%wi(iobs), self%wf(iobs))
       call vert_interp_apply(nlev, salt_d%vals(:,iobs), dsp, self%wi(iobs), self%wf(iobs))

       ! Get insitu temp at model levels and obs location (lono, lato, zo)
       call insitu_t_tl(hofx(iobs),dtp,dsp,self%tempo(iobs),self%salto(iobs),lono,lato,deptho,self%jac(:,iobs))

    enddo

end subroutine ufo_insitutemperature_simobs_tl

! ------------------------------------------------------------------------------
subroutine ufo_insitutemperature_simobs_ad(self, geovals, hofx, obss)
    use ufo_tpsp2ti_mod
    use gsw_pot_to_insitu
    use vert_interp_mod      
    implicit none
    class(ufo_insitutemperature_tlad), intent(in)    :: self
    type(ufo_geovals),       intent(inout) :: geovals
    real(c_double),          intent(in)    :: hofx(:)
    type(c_ptr), value,      intent(in)    :: obss

    character(len=*), parameter :: myname_="ufo_insitutemperature_simobs_ad"
    character(max_string) :: err_msg

    real (kind=kind_real) :: lono, lato, deptho !< Observation location
        
    integer :: iobs, nobs, ilev, nlev
    type(ufo_geoval), pointer :: dtemp, dsalt, dlayerthick
    real (kind_real) :: dtp, dsp
    
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
    
    ! check if sea temperature profile variable is in geovals and get it
    call ufo_geovals_get_var(geovals, var_ocn_pot_temp, dtemp)

    ! check if sea salinity profile variable is in geovals and get it
    call ufo_geovals_get_var(geovals, var_ocn_salt, dsalt)

    ! check if sea layer thickness variable is in geovals get it and zero it out
    call ufo_geovals_get_var(geovals, var_ocn_lay_thick, dlayerthick)
    
    nlev = self%nval
    nobs = self%nobs
    
    if (.not. allocated(dtemp%vals)) allocate(dtemp%vals(nlev, size(hofx,1)))
    if (.not. allocated(dsalt%vals)) allocate(dsalt%vals(nlev, size(hofx,1)))
    if (.not. allocated(dlayerthick%vals)) allocate(dlayerthick%vals(nlev, size(hofx,1)))    

    ! backward sea temperature profile obs operator
    dtemp%vals = 0.0
    dsalt%vals = 0.0
    do iobs = 1, size(hofx,1)

       lono = self%lono(iobs)
       lato = self%lato(iobs)
       deptho = self%deptho(iobs)
      
       ! Adjoint obs operator
       dtp = 0.0
       dsp = 0.0
       call insitu_t_tlad(hofx(iobs),dtp,dsp,self%tempo(iobs),self%salto(iobs),lono,lato,deptho,self%jac(:,iobs))

       ! Backward interpolate
       call vert_interp_apply_ad(nlev, dtemp%vals(:,iobs), dtp, self%wi(iobs), self%wf(iobs))
       call vert_interp_apply_ad(nlev, dsalt%vals(:,iobs), dsp, self%wi(iobs), self%wf(iobs))

       ! Layer thickness is not a control variable: zero it out!
       dlayerthick%vals=0.0
       
    enddo

end subroutine ufo_insitutemperature_simobs_ad

end module ufo_insitutemperature_tlad_mod
