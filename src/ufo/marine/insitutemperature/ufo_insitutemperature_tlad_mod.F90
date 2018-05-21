! (C) Copyright 2017 UCAR
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 

!> Fortran module to handle ice concentration observations

module ufo_insitutemperature_tlad_mod

  use ufo_obs_insitutemperature_mod
  use ufo_obs_vectors
  use ufo_vars_mod
  use ufo_locs_mod
  use ufo_geovals_mod
  use kinds

  implicit none
  public :: ufo_insitutemperature_tlad
  public :: ufo_insitutemperature_tlad_delete
  public :: ufo_insitutemperature_tlad_settraj
  public :: ufo_insitutemperature_tlad_eqv_tl
  public :: ufo_insitutemperature_tlad_eqv_ad
  private
  integer, parameter :: max_string=800

  !> Fortran derived type for sea temperature profile observation operator
  type :: ufo_insitutemperature_tlad
     integer          :: nobs
     integer          :: nval      !< Number of level in model's profiles 
     type(ufo_geoval) :: temp      !< Temperature (traj)
     type(ufo_geoval) :: salt      !< Salinity (traj)
     type(ufo_geoval) :: h         !< Layer thickness (traj)
     real (kind=kind_real), allocatable :: lono(:), lato(:), deptho(:) !< Observation location   
     logical :: ltraj = .false.                  !< trajectory set?
  end type ufo_insitutemperature_tlad


  ! ------------------------------------------------------------------------------

contains

  ! ------------------------------------------------------------------------------

  subroutine ufo_insitutemperature_tlad_delete(self)
    implicit none
    type(ufo_insitutemperature_tlad), intent(inout) :: self

    self%ltraj = .false.

  end subroutine ufo_insitutemperature_tlad_delete

  ! ------------------------------------------------------------------------------

  subroutine ufo_insitutemperature_tlad_settraj(traj, geovals, obs_ti)

    implicit none
    type(ufo_insitutemperature_tlad), intent(inout) :: traj    !< Complete trajectory needed by the operator
    type(ufo_geovals), intent(in)                    :: geovals !< Model background
    type(ufo_obs_insitutemperature), intent(in)     :: obs_ti  !< Insitu temperature observations

    character(len=*), parameter :: myname_="ufo_insitutemperature_tlad_settraj"
    character(max_string) :: err_msg

    type(ufo_geoval), pointer :: temp, salt, h

    ! check if sea temperature profile variable is in geovals and get it
    if (.not. ufo_geovals_get_var(geovals, var_ocn_pot_temp, temp)) then
       write(err_msg,*) myname_, trim(var_ocn_pot_temp), ' doesnt exist'
       call abor1_ftn(err_msg)
    endif

    ! check if sea salinity profile variable is in geovals and get it
    if (.not. ufo_geovals_get_var(geovals, var_ocn_salt, salt)) then
       write(err_msg,*) myname_, trim(var_ocn_salt), ' doesnt exist'
       call abor1_ftn(err_msg)
    endif

    ! check if ocean layer thickness variable is in geovals and get it
    if (.not. ufo_geovals_get_var(geovals, var_ocn_lay_thick, h)) then
       write(err_msg,*) myname_, trim(var_ocn_lay_thick), ' does not exist'
       call abor1_ftn(err_msg)
    endif

    call ufo_insitutemperature_tlad_delete(traj)

    traj%nobs = h%nobs
    traj%nval = h%nval

    traj%temp = temp
    traj%salt = salt
    traj%h    = h

    allocate(traj%lono(obs_ti%nobs))
    allocate(traj%lato(obs_ti%nobs))
    allocate(traj%deptho(obs_ti%nobs))

    traj%lono = obs_ti%lon
    traj%lato = obs_ti%lat
    traj%deptho = obs_ti%depth
    traj%ltraj    = .true.

  end subroutine ufo_insitutemperature_tlad_settraj


  ! ------------------------------------------------------------------------------

  subroutine ufo_insitutemperature_tlad_eqv_tl(traj, geovals, hofx, obs_ti)

    use ufo_tpsp2ti_mod
    use gsw_pot_to_insitu
    use vert_interp_mod    

    implicit none
    type(ufo_insitutemperature_tlad), intent(in) :: traj !< Trajectory
    type(ufo_geovals), intent(in)    :: geovals           !< Increments (dtp, dsp)
    type(obs_vector),  intent(inout) :: hofx              !< dti
    type(ufo_obs_insitutemperature), intent(in) :: obs_ti     !< Insitu temperature observations

    character(len=*), parameter :: myname_="ufo_insitutemperature_tlad_eqv_tl"
    character(max_string) :: err_msg

    integer :: iobs, ilev, nlev, nobs

    type(ufo_geoval), pointer :: temp_d, salt_d !< Increments from geovals
    real (kind=kind_real) :: lono, lato, deptho !< Observation location
    real (kind_real), allocatable :: pressure(:,:), depth(:,:)

    ! Vertical interpolation
    real(kind_real) :: wf, tp, sp, prs, dtp, dsp
    integer :: wi

    ! check if trajectory was set
    if (.not. traj%ltraj) then
       write(err_msg,*) myname_, ' trajectory wasnt set!'
       call abor1_ftn(err_msg)
    endif

    ! check if nobs is consistent in geovals & hofx
    if (geovals%nobs /= hofx%nobs) then
       write(err_msg,*) myname_, ' error: nobs inconsistent!'
       call abor1_ftn(err_msg)
    endif

    ! check if sea temperature profile variable is in geovals and get it
    if (.not. ufo_geovals_get_var(geovals, var_ocn_pot_temp, temp_d)) then
       write(err_msg,*) myname_, trim(var_ocn_pot_temp), ' doesnt exist'
       call abor1_ftn(err_msg)
    endif

    ! check if sea salinity profile variable is in geovals and get it
    if (.not. ufo_geovals_get_var(geovals, var_ocn_salt, salt_d)) then
       write(err_msg,*) myname_, trim(var_ocn_salt), ' doesnt exist'
       call abor1_ftn(err_msg)
    endif

    nlev = temp_d%nval
    nobs = temp_d%nobs        

    allocate(pressure(nlev,nobs), depth(nlev,nobs))
    do iobs = 1,nobs

       lono = traj%lono(iobs)
       lato = traj%lato(iobs)
       deptho = traj%deptho(iobs)

       !< Depth from layer thickness
       depth(1,iobs)=0.5*traj%h%vals(1,iobs)
       do ilev = 2, nlev
          depth(ilev,iobs)=sum(traj%h%vals(1:ilev-1,iobs))+0.5*traj%h%vals(ilev,iobs)
       end do

       !< Pressure from depth          
       do ilev = 1, nlev
          pressure(ilev,iobs)=p_from_z(depth(ilev,iobs),lato)
       end do
       
    end do

    ! linear sea temperature profile obs operator
    hofx%values = 0.0
    do iobs = 1,nobs

       lono = traj%lono(iobs)
       lato = traj%lato(iobs)
       deptho = traj%deptho(iobs)

       !< Interpolation weight
       call vert_interp_weights(nlev,deptho,depth(:,iobs),wi,wf)
       if (deptho.ge.maxval(depth)) then
          wi=nlev-1
          wf=0.0
       end if

       ! Interpolate background
       call vert_interp_apply(nlev, traj%temp%vals(:,iobs), tp, wi, wf)
       call vert_interp_apply(nlev, traj%salt%vals(:,iobs), sp, wi, wf)
       call vert_interp_apply(nlev, pressure(:,iobs), prs, wi, wf)

       !  Interpolate increment
       call vert_interp_apply(nlev, temp_d%vals(:,iobs), dtp, wi, wf)
       call vert_interp_apply(nlev, salt_d%vals(:,iobs), dsp, wi, wf)

       ! Get insitu temp at model levels and obs location (lono, lato, zo)
       call insitu_t_tl(hofx%values(iobs),dtp,dsp,tp,sp,lono,lato,deptho)

       if (isnan(hofx%values(iobs))) then !!!!!! HACK !!!!!!!!!!!!!!!!!!!!!
          print *,'in tlm:',iobs,deptho
          hofx%values(iobs)=0.0 !!!! NEED TO QC OUT BAD OBS LOCATION !!!!!!
       end if

    enddo

  end subroutine ufo_insitutemperature_tlad_eqv_tl

  ! ------------------------------------------------------------------------------

  subroutine ufo_insitutemperature_tlad_eqv_ad(traj, geovals, hofx, obs_ti)

    use ufo_tpsp2ti_mod
    use gsw_pot_to_insitu
    use vert_interp_mod    

    
    implicit none
    type(ufo_insitutemperature_tlad), intent(in) :: traj
    type(ufo_geovals), intent(inout)              :: geovals
    type(obs_vector),  intent(in)                 :: hofx
    type(ufo_obs_insitutemperature), intent(in)  :: obs_ti     !< Insitu temperature observations

    character(len=*), parameter :: myname_="ufo_insitutemperature_tlad_eqv_ad"
    character(max_string) :: err_msg

    real (kind_real), allocatable :: pressure(:,:), depth(:,:)
    real (kind=kind_real) :: lono, lato, deptho !< Observation location
        
    integer :: iobs, nobs, ilev, nlev
    type(ufo_geoval), pointer :: dtemp, dsalt
    real (kind_real) :: dtp, dsp, tp, sp
    ! Vertical interpolation
    real(kind_real) :: wf
    integer :: wi
    
    print *,'--------------coucou0'
    ! check if trajectory was set
    if (.not. traj%ltraj) then
       write(err_msg,*) myname_, ' trajectory wasnt set!'
       call abor1_ftn(err_msg)
    endif

    ! check if nobs is consistent in geovals & hofx
    if (geovals%nobs /= hofx%nobs) then
       write(err_msg,*) myname_, ' error: nobs inconsistent!'
       call abor1_ftn(err_msg)
    endif

    ! check if sea temperature profile variable is in geovals and get it
    if (.not. ufo_geovals_get_var(geovals, var_ocn_pot_temp, dtemp)) then
       write(err_msg,*) myname_, trim(var_ocn_pot_temp), ' doesnt exist'
       call abor1_ftn(err_msg)
    endif

    ! check if sea salinity profile variable is in geovals and get it
    if (.not. ufo_geovals_get_var(geovals, var_ocn_salt, dsalt)) then
       write(err_msg,*) myname_, trim(var_ocn_salt), ' doesnt exist'
       call abor1_ftn(err_msg)
    endif
    print *,'--------------coucou1'
    if (.not. allocated(dtemp%vals)) allocate(dtemp%vals(nlev, hofx%nobs))
    if (.not. allocated(dsalt%vals)) allocate(dsalt%vals(nlev, hofx%nobs))

    nlev = dtemp%nval
    nobs = dtemp%nobs        

    allocate(depth(nlev,nobs))
    do iobs = 1,nobs

       lono = traj%lono(iobs)
       lato = traj%lato(iobs)
       deptho = traj%deptho(iobs)

       !< Depth from layer thickness
       depth(1,iobs)=0.5*traj%h%vals(1,iobs)
       do ilev = 2, nlev
          depth(ilev,iobs)=sum(traj%h%vals(1:ilev-1,iobs))+0.5*traj%h%vals(ilev,iobs)
       end do
    end do
    
    ! backward sea temperature profile obs operator
    dtemp%vals = 0.0
    dsalt%vals = 0.0
    do iobs = 1, hofx%nobs

       lono = traj%lono(iobs)
       lato = traj%lato(iobs)
       deptho = traj%deptho(iobs)
      
       !< Interpolation weight
       call vert_interp_weights(nlev,deptho,depth(:,iobs),wi,wf)
       if (deptho.ge.maxval(depth)) then
          wi=nlev-1
          wf=0.0
       end if

       ! Interpolate background
       call vert_interp_apply(nlev, traj%temp%vals(:,iobs), tp, wi, wf)
       call vert_interp_apply(nlev, traj%salt%vals(:,iobs), sp, wi, wf)
       
       ! Get insitu temp at model levels and obs location (lono, lato, zo)
       call insitu_t_tlad(hofx%values(iobs), dtp, dsp, tp, sp, lono, lato, deptho)

       !
       !  Backward interpolate
       call vert_interp_apply_ad(nlev, dtemp%vals(:,iobs), dtp, wi, wf)
       call vert_interp_apply_ad(nlev, dsalt%vals(:,iobs), dsp, wi, wf)       

!!$          print *,'adt=',dtemp%vals(:,iobs)
!!$          print *,'hofx%values(iobs)=',hofx%values(iobs)
!!$          print *,'dtp, dsp, tp, sp=',dtp, dsp, tp, sp
!!$          !print *,'traj temp=',traj%temp%vals(:,iobs)
!!$          read(*,*)          

       if (isnan(dtp*dsp)) then
          print *,'crap in adj:',iobs,sum(dtemp%vals(:,iobs))
          dtemp%vals(:,iobs)=0.0
          dsalt%vals(:,iobs)=0.0
       end if
       if (isnan(sum(dtemp%vals(1:nlev,iobs)))) then
          print *,'crap in adj:',iobs,sum(dtemp%vals(:,iobs))
          dtemp%vals(:,iobs)=0.0
          dsalt%vals(:,iobs)=0.0
          print *,'adt=',dtemp%vals(:,iobs)
          print *,'wi,wf=',wi,wf
          print *,'dtp, dsp, tp, sp=',dtp, dsp, tp, sp
          !print *,'traj temp=',traj%temp%vals(:,iobs)
          !read(*,*)          
       end if
       print *,'dtp, dsp, tp, sp=',dtp, dsp, tp, sp
       print *,'dtemp=',dtemp%vals(:,iobs)
       !read(*,*)          

    enddo

  end subroutine ufo_insitutemperature_tlad_eqv_ad

end module ufo_insitutemperature_tlad_mod
