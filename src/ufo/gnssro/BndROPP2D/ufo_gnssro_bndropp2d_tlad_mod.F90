! (C) Copyright 2017-2018 UCAR
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 

!> Fortran module for gnssro bending angle ropp2d tangent linear and adjoint
!> following the ROPP (2018 Aug) implementation

module ufo_gnssro_bndropp2d_tlad_mod

use iso_c_binding
use kinds
use ufo_vars_mod
use ufo_geovals_mod
use ufo_geovals_mod_c,   only: ufo_geovals_registry
use vert_interp_mod
use ufo_basis_tlad_mod,  only: ufo_basis_tlad
use obsspace_mod
use config_mod
use gnssro_mod_conf
use missing_values_mod
use ufo_gnssro_ropp2d_utils_mod
use fckit_log_module, only : fckit_log

integer, parameter         :: max_string=800

!> Fortran derived type for gnssro trajectory
type, extends(ufo_basis_tlad)   ::  ufo_gnssro_BndROPP2D_tlad
  private
  integer                       :: nval, nlocs
  real(kind_real), allocatable  :: prs(:,:), t(:,:), q(:,:), gph(:,:), gph_sfc(:,:)
  integer                       :: n_horiz      ! 2d points along ray path
  integer                       :: iflip        ! geoval ascending order flag
  type(gnssro_conf)             :: roconf       ! ro configuration
  real(kind_real), allocatable  :: obsLon2d(:), obsLat2d(:)  !2d locations - nlocs*n_horiz
  contains
    procedure :: setup      => ufo_gnssro_bndropp2d_tlad_setup
    procedure :: delete     => ufo_gnssro_bndropp2d_tlad_delete
    procedure :: settraj    => ufo_gnssro_bndropp2d_tlad_settraj
    procedure :: simobs_tl  => ufo_gnssro_bndropp2d_simobs_tl
    procedure :: simobs_ad  => ufo_gnssro_bndropp2d_simobs_ad
end type ufo_gnssro_bndropp2d_tlad

contains

! ------------------------------------------------------------------------------
subroutine ufo_gnssro_bndropp2d_tlad_setup(self, c_conf)
  implicit none
  class(ufo_gnssro_BndROPP2D_tlad), intent(inout) :: self
  type(c_ptr),                      intent(in)    :: c_conf

  call gnssro_conf_setup(self%roconf,c_conf)

end subroutine ufo_gnssro_bndropp2d_tlad_setup

! ------------------------------------------------------------------------------    
subroutine ufo_gnssro_bndropp2d_tlad_settraj(self, geovals, obss)
       
  implicit none
  class(ufo_gnssro_BndROPP2D_tlad), intent(inout) :: self
  type(ufo_geovals),                intent(in)    :: geovals
  type(c_ptr), value,               intent(in)    :: obss
  character(len=*), parameter :: myname_="ufo_gnssro_bndropp2d_tlad_settraj"
  character(max_string)       :: err_msg
  type(ufo_geoval), pointer   :: t, q, prs, gph, gph_sfc
  integer                     :: i, iobs, ierr, kerror
  real(kind_real), allocatable  :: obsAzim(:)                    ! nlocs
  real(kind_real), allocatable  :: obsLat(:), obsLon(:)          ! nlocs
  real(kind_real), allocatable  :: obsLon2d(:),obsLat2d(:)       ! nlocs * n_horiz
  real(kind_real), allocatable  :: obsLonnh(:),obsLatnh(:)       ! n_horiz
  integer                       :: n_horiz
  real(kind_real)               :: dtheta

  write(err_msg,*) "TRACE: ufo_gnssro_bndropp2d_tlad_settraj: begin"
  call fckit_log%info(err_msg) 

! get model state variables from geovals
  call ufo_geovals_get_var(geovals, var_ts,    t)         ! temperature
  call ufo_geovals_get_var(geovals, var_q,     q)         ! specific humidity
  call ufo_geovals_get_var(geovals, var_prs,   prs)       ! pressure
  call ufo_geovals_get_var(geovals, var_z,     gph)       ! geopotential height
  call self%delete()   

  self%nval    = prs%nval
  self%nlocs    = obsspace_get_nlocs(obss)
  self%iflip   = 0

  n_horiz = self%roconf%n_horiz
  dtheta  = self%roconf%dtheta

  if (prs%vals(1,1) .lt. prs%vals(prs%nval,1) ) then
    self%iflip = 1 
    write(err_msg,'(a)') '  ufo_gnssro_bndropp2d_tlad_settraj:'//new_line('a')//                   &
                         '  Model vertical height profile is in descending order,'//new_line('a')// &
                         '  but ROPP requires it to be ascending order, need flip'
    call fckit_log%info(err_msg)
  end if

  allocate(self%obsLat2d(self%nlocs*n_horiz))
  allocate(self%obsLon2d(self%nlocs*n_horiz))

  allocate(obsLon(self%nlocs))
  allocate(obsLat(self%nlocs))
  allocate(obsAzim(self%nlocs))

  call obsspace_get_db(obss, "MetaData", "longitude",        obsLon)
  call obsspace_get_db(obss, "MetaData", "latitude",         obsLat)
  call obsspace_get_db(obss, "MetaData", "sensor_azimuth_angle", obsAzim)

  allocate(obsLatnh(n_horiz))
  allocate(obsLonnh(n_horiz))

  do i = 1, self%nlocs
     call ropp_fm_2d_plane(obsLat(i),obsLon(i),obsAzim(i),dtheta,n_horiz,obsLatnh,obsLonnh,kerror)
     self%obsLon2d((i-1)*n_horiz+1:i*n_horiz) =  obsLonnh
     self%obsLat2d((i-1)*n_horiz+1:i*n_horiz) =  obsLatnh
  end do

  deallocate(obsLat)
  deallocate(obsLon)
  deallocate(obsLonnh)
  deallocate(obsLatnh)
  deallocate(obsAzim)

  allocate(self%t(self%nval,self%nlocs*n_horiz))
  allocate(self%q(self%nval,self%nlocs*n_horiz))
  allocate(self%prs(self%nval,self%nlocs*n_horiz))
  allocate(self%gph(self%nval,self%nlocs*n_horiz))

! allocate  
  self%gph     = gph%vals
  self%t       = t%vals
  self%q       = q%vals
  self%prs     = prs%vals

  self%ltraj   = .true.
       
end subroutine ufo_gnssro_bndropp2d_tlad_settraj
    
! ------------------------------------------------------------------------------
! ------------------------------------------------------------------------------    
subroutine ufo_gnssro_bndropp2d_simobs_tl(self, geovals, hofx, obss)

  use ropp_fm_types, only: State2dFM
  use ropp_fm_types, only: Obs1dBangle
  use datetimetypes, only: dp
  implicit none
  class(ufo_gnssro_BndROPP2D_tlad), intent(in)    :: self
  type(ufo_geovals),                intent(in)    :: geovals   ! perturbed quantities
  real(kind_real),                  intent(inout) :: hofx(:)
  type(c_ptr),   value,             intent(in)    :: obss

  type(State2dFM)                 :: x,x_tl
  type(Obs1dBangle)               :: y,y_tl

  integer                         :: iobs,nlev, nlocs,ierr,nvprof
    
  character(len=*), parameter  :: myname_="ufo_gnssro_bndropp2d_simobs_tl"
  character(max_string)        :: err_msg
  type(ufo_geoval), pointer    :: t_d, q_d, prs_d 

! hack - set local geopotential height to zero for ropp routines
  real(kind_real), allocatable  :: gph_d_zero(:,:)
! hack - set local geopotential height to zero for ropp routines
  real(kind_real), allocatable  :: obsImpP(:),obsLocR(:),obsGeoid(:),obsAzim(:) !nlocs
  real(kind_real), allocatable  :: obsLat(:),obsLon(:)                          !nlocs
  real(kind_real), allocatable  :: obsLon2d(:),obsLat2d(:)       ! nlocs * n_horiz
  real(kind_real), allocatable  :: obsLonnh(:),obsLatnh(:)       ! n_horiz
  integer                         :: n_horiz
  real(kind_real)                 :: dtheta

  n_horiz = self%roconf%n_horiz
  dtheta  = self%roconf%dtheta

  write(err_msg,*) "TRACE: ufo_gnssro_bndropp2d_simobs_tl: begin"
  call fckit_log%info(err_msg)

! check if trajectory was set
  if (.not. self%ltraj) then
     write(err_msg,*) myname_, ' trajectory wasnt set!'
     call abor1_ftn(err_msg)
  endif
      
! check if nlocs is consistent in geovals & hofx
  if (geovals%nlocs /= size(hofx)*n_horiz ) then
     write(err_msg,*) myname_, ' error: 2d nlocs inconsistent!'
     call abor1_ftn(err_msg)
  endif

! get variables from geovals
  call ufo_geovals_get_var(geovals, var_ts,     t_d)         ! temperature
  call ufo_geovals_get_var(geovals, var_q,     q_d)         ! specific humidity
  call ufo_geovals_get_var(geovals, var_prs,   prs_d)       ! pressure

  nlev    = self%nval
  nlocs   = self%nlocs

  allocate(gph_d_zero(nlev,nlocs*n_horiz))
  gph_d_zero     = 0.0

! set obs space struture
  allocate(obsLon(nlocs))
  allocate(obsLat(nlocs))
  allocate(obsImpP(nlocs))
  allocate(obsLocR(nlocs))
  allocate(obsGeoid(nlocs))

  call obsspace_get_db(obss, "MetaData", "longitude",        obsLon)
  call obsspace_get_db(obss, "MetaData", "latitude",         obsLat)
  call obsspace_get_db(obss, "MetaData", "impact_parameter", obsImpP)
  call obsspace_get_db(obss, "MetaData", "earth_radius_of_curvature", obsLocR)
  call obsspace_get_db(obss, "MetaData", "geoid_height_above_reference_ellipsoid", obsGeoid)

  nvprof = 1  ! no. of bending angles in profile 

! loop through the obs
  obs_loop: do iobs = 1, nlocs   ! order of loop doesn't matter

!   map the trajectory to ROPP structure x
    call init_ropp_2d_statevec(self%obsLon2d( (iobs-1)*n_horiz+1:iobs*n_horiz ), &
                               self%obsLat2d( (iobs-1)*n_horiz+1:iobs*n_horiz ), &
                               self%t(:,(iobs-1)*n_horiz+1:iobs*n_horiz),    &
                               self%q(:,(iobs-1)*n_horiz+1:iobs*n_horiz),    &
                               self%prs(:,(iobs-1)*n_horiz+1:iobs*n_horiz),  &
                               self%gph(:,(iobs-1)*n_horiz+1:iobs*n_horiz),  &
                               nlev, x, n_horiz, dtheta, self%iflip)

!  hack -- make non zero humidity to avoid zero denominator in tangent linear
!          see  ropp_fm/bangle_1d/ropp_fm_bangle_1d_tl.f90
    where(x%shum .le. 1e-8)        x%shum = 1e-8
!  hack -- make non zero humidity to avoid zero denominator in tangent linear

    call init_ropp_2d_statevec(self%obsLon2d( (iobs-1)*n_horiz+1:iobs*n_horiz ), &
                               self%obsLat2d( (iobs-1)*n_horiz+1:iobs*n_horiz ), &
                               t_d%vals(:,(iobs-1)*n_horiz+1:iobs*n_horiz),    &
                               q_d%vals(:,(iobs-1)*n_horiz+1:iobs*n_horiz),    &
                               prs_d%vals(:,(iobs-1)*n_horiz+1:iobs*n_horiz),  &
                               gph_d_zero(:,(iobs-1)*n_horiz+1:iobs*n_horiz),  &
                               nlev, x_tl, n_horiz, dtheta, self%iflip)

!   set both y and y_tl structures    
    call init_ropp_2d_obvec_tlad(iobs, nvprof, &
                      obsImpP(iobs),           &
                      obsLat(iobs),            &
                      obsLon(iobs),            &
                      obsLocR(iobs),           &
                      obsGeoid(iobs),          &
                             y,y_tl)

!   now call TL of forward model
    call ropp_fm_bangle_2d_tl(x,x_tl,y, y_tl)
    hofx(iobs) = y_tl%bangle(nvprof) ! this will need to change if profile is passed

!   tidy up -deallocate ropp structures 
    call ropp_tidy_up_tlad_2d(x,x_tl,y,y_tl)

  end do obs_loop

! tidy up - deallocate obsspace structures
  deallocate(obsLat) 
  deallocate(obsLon)
  deallocate(obsImpP)
  deallocate(obsLocR)
  deallocate(obsGeoid)

  write(err_msg,*) "TRACE: ufo_gnssro_bndropp2d_simobs_tl: complete"
  call fckit_log%info(err_msg)

  return
    
end subroutine ufo_gnssro_bndropp2d_simobs_tl
 
! ------------------------------------------------------------------------------
! ------------------------------------------------------------------------------
subroutine ufo_gnssro_bndropp2d_simobs_ad(self, geovals, hofx, obss)

  use ropp_fm_types, only: State2dFM
  use ropp_fm_types, only: Obs1dBangle
  use typesizes,     only: wp => EightByteReal
  use datetimetypes, only: dp

  implicit none
  class(ufo_gnssro_BndROPP2D_tlad), intent(in)    :: self
  type(ufo_geovals),                intent(inout) :: geovals   ! perturbed quantities
  real(kind_real),                  intent(in)    :: hofx(:)
  type(c_ptr),  value,              intent(in)    :: obss
  real(c_double)              :: missing

  type(ufo_geoval),     pointer   :: t_d, q_d, prs_d 

! set local geopotential height to zero for ropp routines
  real(kind_real),    allocatable :: gph_d_zero(:,:)

  real(kind_real),    allocatable :: obsLat(:), obsLon(:), obsImpP(:), obsLocR(:), obsGeoid(:)
  type(State2dFM)                 :: x,x_ad
  type(Obs1dBangle)               :: y,y_ad
  integer                         :: iobs,nlev,nlocs,ierr,nvprof
  character(len=*), parameter     :: myname_="ufo_gnssro_bndropp2d_simobs_ad"
  character(max_string)           :: err_msg
  integer                         :: n_horiz
  real(kind_real)                 :: dtheta

  n_horiz = self%roconf%n_horiz
  dtheta  = self%roconf%dtheta

  write(err_msg,*) "TRACE: ufo_gnssro_bndropp2d_simobs_ad: begin"
  call fckit_log%info(err_msg)

! check if trajectory was set
  if (.not. self%ltraj) then
     write(err_msg,*) myname_, ' trajectory wasnt set!'
     call abor1_ftn(err_msg)
  endif
! check if nlocs is consistent in geovals & hofx
  if (geovals%nlocs /= size(hofx)*n_horiz) then
     write(err_msg,*) myname_, ' error: 2d nlocs inconsistent!'
     call abor1_ftn(err_msg)
  endif
     
! get variables from geovals
  call ufo_geovals_get_var(geovals, var_ts,    t_d)         ! temperature
  call ufo_geovals_get_var(geovals, var_q,     q_d)         ! specific humidity
  call ufo_geovals_get_var(geovals, var_prs,   prs_d)       ! pressure

! allocate if not yet allocated   
  if (.not. allocated(t_d%vals)) then
      t_d%nlocs = self%nlocs*n_horiz
      t_d%nval = self%nval
      allocate(t_d%vals(t_d%nval,t_d%nlocs))
      t_d%vals = 0.0_kind_real
  endif

  if (.not. allocated(prs_d%vals)) then
      prs_d%nlocs = self%nlocs*n_horiz
      prs_d%nval = self%nval
      allocate(prs_d%vals(prs_d%nval,prs_d%nlocs))
      prs_d%vals = 0.0_kind_real
  endif

  if (.not. allocated(q_d%vals)) then
      q_d%nlocs = self%nlocs*n_horiz
      q_d%nval = self%nval
      allocate(q_d%vals(q_d%nval,q_d%nlocs))
      q_d%vals = 0.0_kind_real
  endif

  if (.not. geovals%linit ) geovals%linit=.true.

  nlev    = self%nval
  nlocs   = self%nlocs

  allocate(gph_d_zero(nlev,nlocs*n_horiz))
  gph_d_zero = 0.0

! set obs space struture
  allocate(obsLon(nlocs))
  allocate(obsLat(nlocs))
  allocate(obsImpP(nlocs))
  allocate(obsLocR(nlocs))
  allocate(obsGeoid(nlocs))

  call obsspace_get_db(obss, "MetaData", "longitude",        obsLon)
  call obsspace_get_db(obss, "MetaData", "latitude",         obsLat)
  call obsspace_get_db(obss, "MetaData", "impact_parameter", obsImpP)
  call obsspace_get_db(obss, "MetaData", "earth_radius_of_curvature", obsLocR)
  call obsspace_get_db(obss, "MetaData", "geoid_height_above_reference_ellipsoid", obsGeoid)

  missing = missing_value(missing)

! loop through the obs
  nvprof=1  ! no. of bending angles in profile 
  obs_loop: do iobs = 1, nlocs 

    if (hofx(iobs) .gt. missing) then

!       map the trajectory to ROPP structure x
        call init_ropp_2d_statevec(self%obsLon2d((iobs-1)*n_horiz+1:iobs*n_horiz), &
                                   self%obsLat2d((iobs-1)*n_horiz+1:iobs*n_horiz), &
                                   self%t(:,(iobs-1)*n_horiz+1:iobs*n_horiz),    &
                                   self%q(:,(iobs-1)*n_horiz+1:iobs*n_horiz),    &
                                   self%prs(:,(iobs-1)*n_horiz+1:iobs*n_horiz),  &
                                   self%gph(:,(iobs-1)*n_horiz+1:iobs*n_horiz),  &
                                   nlev, x, n_horiz, dtheta, self%iflip)

        call init_ropp_2d_statevec(self%obsLon2d( (iobs-1)*n_horiz+1:iobs*n_horiz ), &
                                   self%obsLat2d( (iobs-1)*n_horiz+1:iobs*n_horiz ), &
                                   t_d%vals(:,(iobs-1)*n_horiz+1:iobs*n_horiz),    &
                                   q_d%vals(:,(iobs-1)*n_horiz+1:iobs*n_horiz),    &
                                   prs_d%vals(:,(iobs-1)*n_horiz+1:iobs*n_horiz),  &
                                   gph_d_zero(:,(iobs-1)*n_horiz+1:iobs*n_horiz),  &
                                   nlev, x_ad, n_horiz, dtheta, self%iflip)

 !      x_ad is local so initialise to 0.0
        x_ad%temp(:,:) = 0.0_wp
        x_ad%pres(:,:) = 0.0_wp
        x_ad%shum(:,:) = 0.0_wp
        x_ad%geop(:,:) = 0.0_wp
 
 !      set both y and y_ad structures    
        call init_ropp_2d_obvec_tlad(iobs,  nvprof,  &
                         obsImpP(iobs),              &
                         obsLat(iobs),               &
                         obsLon(iobs),               &
                         obsLocR(iobs),              &
                         obsGeoid(iobs),             &
                                y,y_ad)


!       local variable initialise
        y_ad%bangle(:) = 0.0_wp

!       now call AD of forward model
        y_ad%bangle(nvprof)  = y_ad%bangle(nvprof) + hofx(iobs)
        call ropp_fm_bangle_2d_ad(x,x_ad,y,y_ad)

        call init_ropp_2d_statevec_ad(           &
                          t_d%vals(:,(iobs-1)*n_horiz+1:iobs*n_horiz),      &
                          q_d%vals(:,(iobs-1)*n_horiz+1:iobs*n_horiz),      &
                        prs_d%vals(:,(iobs-1)*n_horiz+1:iobs*n_horiz),      &
                        gph_d_zero(:,(iobs-1)*n_horiz+1:iobs*n_horiz),      &
                        nlev, x_ad, n_horiz,self%iflip)

!     tidy up - deallocate ropp structures  
      call ropp_tidy_up_tlad_2d(x,x_ad,y,y_ad)

    end if  ! end missing value check

  end do obs_loop

! tidy up - deallocate obsspace structures
  deallocate(obsLat) 
  deallocate(obsLon)
  deallocate(obsImpP)
  deallocate(obsLocR)
  deallocate(obsGeoid)
  deallocate(gph_d_zero)

  write(err_msg,*) "TRACE: ufo_gnssro_bndropp2d_simobs_ad: complete"
  call fckit_log%info(err_msg)

  return

end subroutine ufo_gnssro_bndropp2d_simobs_ad
    
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
subroutine ufo_gnssro_bndropp2d_tlad_delete(self)

  implicit none
  class(ufo_gnssro_BndROPP2D_tlad), intent(inout) :: self
  character(len=*), parameter :: myname_="ufo_gnssro_bndropp_tlad_delete"
      
  self%nval = 0
  if (allocated(self%prs)) deallocate(self%prs)
  if (allocated(self%t))   deallocate(self%t)
  if (allocated(self%q))   deallocate(self%q)
  if (allocated(self%gph)) deallocate(self%gph)

  self%ltraj = .false. 

end subroutine ufo_gnssro_bndropp2d_tlad_delete

!-------------------------------------------------------------------------
       
end module ufo_gnssro_bndropp2d_tlad_mod
