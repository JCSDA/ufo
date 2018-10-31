! (C) Copyright 2017-2018 UCAR
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 

!> Fortran module to handle ice concentration observations

module ufo_insitutemperature_tlad_mod

  use iso_c_binding
  use obsspace_mod
  use ufo_vars_mod
  use ioda_locs_mod
  use ufo_geovals_mod
  use kinds
  
  implicit none
  public :: ufo_insitutemperature_tlad
  public :: ufo_insitutemperature_tlad_delete
  public :: ufo_insitutemperature_tlad_settraj
  public :: ufo_insitutemperature_simobs_tl
  public :: ufo_insitutemperature_simobs_ad
  private
  integer, parameter :: max_string=800

  !> Fortran derived type to hold trajectory
  !> for ocean insitu temperature observation operator
  type :: ufo_insitutemperature_tlad
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
     logical                            :: ltraj = .false. !< trajectory set?
  end type ufo_insitutemperature_tlad

  ! ------------------------------------------------------------------------------

contains

  ! ------------------------------------------------------------------------------

  subroutine ufo_insitutemperature_tlad_delete(self)
    implicit none
    type(ufo_insitutemperature_tlad), intent(inout) :: self

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

  subroutine ufo_insitutemperature_tlad_settraj(traj, geovals, obss)
    use vert_interp_mod
    use ufo_tpsp2ti_mod
        
    implicit none
    type(ufo_insitutemperature_tlad), intent(inout)  :: traj    !< Complete trajectory needed by the operator
    type(ufo_geovals), intent(in)                    :: geovals !< Model background
    type(c_ptr), value, intent(in)                   :: obss    !< Insitu temperature observations

    character(len=*), parameter :: myname_="ufo_insitutemperature_tlad_settraj"
    character(max_string) :: err_msg

    type(ufo_geoval), pointer :: temp, salt, h
    integer :: nobs, nlev, iobs, ilev

    real(kind_real), allocatable :: obs_lat(:)
    real(kind_real), allocatable :: obs_lon(:)
    real(kind_real), allocatable :: obs_depth(:)
    integer :: obss_nobs
    
    ! check if sea temperature profile variable is in geovals and get it
    call ufo_geovals_get_var(geovals, var_ocn_pot_temp, temp)

    ! check if sea salinity profile variable is in geovals and get it
    call ufo_geovals_get_var(geovals, var_ocn_salt, salt)

    ! check if ocean layer thickness variable is in geovals and get it
    call ufo_geovals_get_var(geovals, var_ocn_lay_thick, h)

    call ufo_insitutemperature_tlad_delete(traj)

    nobs = h%nobs
    nlev = h%nval
    
    traj%nobs = nobs
    traj%nval = nlev
    
    traj%temp = temp
    traj%salt = salt
    traj%h    = h

    allocate(traj%lono(nobs))
    allocate(traj%lato(nobs))
    allocate(traj%deptho(nobs))

    obss_nobs = obsspace_get_nobs(obss)
    allocate(obs_lat(obss_nobs))
    allocate(obs_lon(obss_nobs))
    allocate(obs_depth(obss_nobs))

    call obsspace_get_db(obss, "Metadata", "longitude", obss_nobs, obs_lon)
    call obsspace_get_db(obss, "Metadata", "latitude", obss_nobs, obs_lat)
    call obsspace_get_db(obss, "Metadata", "ocean_depth", obss_nobs, obs_depth)

    traj%lono = obs_lon
    traj%lato = obs_lat
    traj%deptho = obs_depth

    !< Depth from layer thickness
    allocate(traj%depth(nlev,nobs))
    do iobs = 1, nobs
       traj%depth(1,iobs)=0.5*traj%h%vals(1,iobs)       
       do ilev = 2, nlev
          traj%depth(ilev,iobs)=sum(traj%h%vals(1:ilev-1,iobs))+0.5*traj%h%vals(ilev,iobs)
       end do
    end do

    !< Interpolation weight
    allocate(traj%wi(nobs),traj%wf(nobs))
    do iobs = 1, nobs    
       call vert_interp_weights(nlev,traj%deptho(iobs),traj%depth(:,iobs),traj%wi(iobs),traj%wf(iobs))
       if (traj%deptho(iobs).ge.maxval(traj%depth(:,iobs))) then
          traj%wi(iobs)=nlev-1
          traj%wf(iobs)=0.0
       end if
    end do
    traj%ltraj    = .true.

    !< Jacobian
    allocate(traj%jac(2,nobs),traj%tempo(nobs),traj%salto(nobs))
    do iobs = 1, nobs
       ! Interpolate background do obs depth and save in traj
       call vert_interp_apply(nlev, traj%temp%vals(:,iobs), traj%tempo(iobs), traj%wi(iobs), traj%wf(iobs))
       call vert_interp_apply(nlev, traj%salt%vals(:,iobs), traj%salto(iobs), traj%wi(iobs), traj%wf(iobs))
       
       ! Compute jacobian
       call insitu_t_jac(traj%jac(:,iobs), traj%tempo(iobs), traj%salto(iobs), traj%lono(iobs), traj%lato(iobs), traj%deptho(iobs))
    end do
    
    deallocate(obs_lat)
    deallocate(obs_lon)
    deallocate(obs_depth)
  end subroutine ufo_insitutemperature_tlad_settraj

  ! ------------------------------------------------------------------------------

  subroutine ufo_insitutemperature_simobs_tl(traj, geovals, hofx)

    use ufo_tpsp2ti_mod
    use gsw_pot_to_insitu
    use vert_interp_mod    

    implicit none
    type(ufo_insitutemperature_tlad), intent(in) :: traj !< Trajectory
    type(ufo_geovals), intent(in)    :: geovals           !< Increments (dtp, dsp)
    real(c_double),  intent(inout) :: hofx(:)              !< dti

    character(len=*), parameter :: myname_="ufo_insitutemperature_simobs_tl"
    character(max_string) :: err_msg

    integer :: iobs, ilev, nlev, nobs

    type(ufo_geoval), pointer :: temp_d, salt_d, dlayerthick !< Increments from geovals
    real (kind=kind_real) :: lono, lato, deptho !< Observation location

    ! Vertical interpolation
    real(kind_real) :: dtp, dsp

    ! check if trajectory was set
    if (.not. traj%ltraj) then
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

       lono = traj%lono(iobs)
       lato = traj%lato(iobs)
       deptho = traj%deptho(iobs)

       !  Interpolate increment
       call vert_interp_apply(nlev, temp_d%vals(:,iobs), dtp, traj%wi(iobs), traj%wf(iobs))
       call vert_interp_apply(nlev, salt_d%vals(:,iobs), dsp, traj%wi(iobs), traj%wf(iobs))

       ! Get insitu temp at model levels and obs location (lono, lato, zo)
       call insitu_t_tl(hofx(iobs),dtp,dsp,traj%tempo(iobs),traj%salto(iobs),lono,lato,deptho,traj%jac(:,iobs))

    enddo

  end subroutine ufo_insitutemperature_simobs_tl

  ! ------------------------------------------------------------------------------

  subroutine ufo_insitutemperature_simobs_ad(traj, geovals, hofx)

    use ufo_tpsp2ti_mod
    use gsw_pot_to_insitu
    use vert_interp_mod    
    
    implicit none
    type(ufo_insitutemperature_tlad), intent(in)  :: traj
    type(ufo_geovals), intent(inout)              :: geovals
    real(c_double),  intent(in)                   :: hofx(:)

    character(len=*), parameter :: myname_="ufo_insitutemperature_simobs_ad"
    character(max_string) :: err_msg

    real (kind=kind_real) :: lono, lato, deptho !< Observation location
        
    integer :: iobs, nobs, ilev, nlev
    type(ufo_geoval), pointer :: dtemp, dsalt, dlayerthick
    real (kind_real) :: dtp, dsp
    
    ! check if trajectory was set
    if (.not. traj%ltraj) then
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
    
    nlev = traj%nval
    nobs = traj%nobs
    
    if (.not. allocated(dtemp%vals)) allocate(dtemp%vals(nlev, size(hofx,1)))
    if (.not. allocated(dsalt%vals)) allocate(dsalt%vals(nlev, size(hofx,1)))
    if (.not. allocated(dlayerthick%vals)) allocate(dlayerthick%vals(nlev, size(hofx,1)))    

    ! backward sea temperature profile obs operator
    dtemp%vals = 0.0
    dsalt%vals = 0.0
    do iobs = 1, size(hofx,1)

       lono = traj%lono(iobs)
       lato = traj%lato(iobs)
       deptho = traj%deptho(iobs)
      
       ! Adjoint obs operator
       dtp = 0.0
       dsp = 0.0
       call insitu_t_tlad(hofx(iobs),dtp,dsp,traj%tempo(iobs),traj%salto(iobs),lono,lato,deptho,traj%jac(:,iobs))

       ! Backward interpolate
       call vert_interp_apply_ad(nlev, dtemp%vals(:,iobs), dtp, traj%wi(iobs), traj%wf(iobs))
       call vert_interp_apply_ad(nlev, dsalt%vals(:,iobs), dsp, traj%wi(iobs), traj%wf(iobs))

       ! Layer thickness is not a control variable: zero it out!
       dlayerthick%vals=0.0
       
    enddo

  end subroutine ufo_insitutemperature_simobs_ad

end module ufo_insitutemperature_tlad_mod
