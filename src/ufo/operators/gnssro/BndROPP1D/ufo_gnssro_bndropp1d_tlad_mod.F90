! (C) Copyright 2017-2018 UCAR
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 

!> Fortran module for gnssro bending angle ropp1d tangent linear and adjoint
!> following the ROPP (2018 Aug) implementation

module ufo_gnssro_bndropp1d_tlad_mod

use fckit_configuration_module, only: fckit_configuration
!use iso_c_binding
use kinds
use ufo_vars_mod
use ufo_geovals_mod
use ufo_geovals_mod_c,   only: ufo_geovals_registry
use vert_interp_mod
use ufo_basis_tlad_mod,  only: ufo_basis_tlad
use obsspace_mod
use gnssro_mod_conf
use missing_values_mod
use ufo_gnssro_ropp1d_utils_mod
use fckit_log_module, only : fckit_log

private
public :: ufo_gnssro_BndROPP1D_tlad

integer, parameter         :: max_string=800

!> Fortran derived type for gnssro trajectory
type, extends(ufo_basis_tlad)   ::  ufo_gnssro_BndROPP1D_tlad
  private
  integer                       :: nval, nlocs, iflip
  type(gnssro_conf)             :: roconf       ! ro configuration
  real(kind_real), allocatable  :: prs(:,:), t(:,:), q(:,:), gph(:,:), gph_sfc(:,:)
  contains
    procedure :: setup      => ufo_gnssro_bndropp1d_tlad_setup
    procedure :: delete     => ufo_gnssro_bndropp1d_tlad_delete
    procedure :: settraj    => ufo_gnssro_bndropp1d_tlad_settraj
    procedure :: simobs_tl  => ufo_gnssro_bndropp1d_simobs_tl
    procedure :: simobs_ad  => ufo_gnssro_bndropp1d_simobs_ad
end type ufo_gnssro_bndropp1d_tlad

contains

! ------------------------------------------------------------------------------
subroutine ufo_gnssro_bndropp1d_tlad_setup(self, f_conf)
  implicit none
  class(ufo_gnssro_BndROPP1D_tlad), intent(inout) :: self
  type(fckit_configuration), intent(in)           :: f_conf

  call gnssro_conf_setup(self%roconf,f_conf)

end subroutine ufo_gnssro_bndropp1d_tlad_setup
! ------------------------------------------------------------------------------    
subroutine ufo_gnssro_bndropp1d_tlad_settraj(self, geovals, obss)
       
  implicit none
  class(ufo_gnssro_BndROPP1D_tlad), intent(inout) :: self
  type(ufo_geovals),                intent(in)    :: geovals
  type(c_ptr), value,               intent(in)    :: obss
  character(len=*), parameter :: myname_="ufo_gnssro_bndropp1d_tlad_settraj"
  character(max_string)       :: err_msg
  type(ufo_geoval), pointer   :: t, q, prs, gph, gph_sfc

  write(err_msg,*) "TRACE: ufo_gnssro_bndropp1d_tlad_settraj: begin"
  call fckit_log%debug(err_msg) 

! get model state variables from geovals
  call ufo_geovals_get_var(geovals, var_ts,    t)         ! temperature
  call ufo_geovals_get_var(geovals, var_q,     q)         ! specific humidity
  call ufo_geovals_get_var(geovals, var_prs,   prs)       ! pressure
  call ufo_geovals_get_var(geovals, var_z,     gph)       ! geopotential height
  call ufo_geovals_get_var(geovals, var_sfc_geomz, gph_sfc)   ! surface geopotential height
      
  call self%delete()   

! Keep copy of dimensions
  self%nval = prs%nval
  self%nlocs = obsspace_get_nlocs(obss)
  
  if (self%nlocs > 0) then
     self%iflip = 0
     if (prs%vals(1,1) .lt. prs%vals(prs%nval,1) ) then
       self%iflip = 1
       write(err_msg,'(a)') '  ufo_gnssro_bndropp1d_tlad_settraj:'//new_line('a')//                   &
                            '  Model vertical height profile is in descending order,'//new_line('a')// &
                            '  but ROPP requires it to be ascending order, need flip'
       call fckit_log%debug(err_msg)
     end if

     allocate(self%t(self%nval,self%nlocs))
     allocate(self%q(self%nval,self%nlocs))
     allocate(self%prs(self%nval,self%nlocs))
     allocate(self%gph(self%nval,self%nlocs))
     allocate(self%gph_sfc(1,self%nlocs))
     self%gph     = gph%vals
     self%t       = t%vals
     self%q       = q%vals
     self%prs     = prs%vals
     self%gph_sfc = gph_sfc%vals
  end if
  self%ltraj   = .true.
       
end subroutine ufo_gnssro_bndropp1d_tlad_settraj
    
! ------------------------------------------------------------------------------
! ------------------------------------------------------------------------------    
subroutine ufo_gnssro_bndropp1d_simobs_tl(self, geovals, hofx, obss)

  use ropp_fm_types, only: State1dFM
  use ropp_fm_types, only: Obs1dBangle
  use datetimetypes, only: dp
  implicit none
  class(ufo_gnssro_BndROPP1D_tlad), intent(in)    :: self
  type(ufo_geovals),                intent(in)    :: geovals   ! perturbed quantities
  real(kind_real),                  intent(inout) :: hofx(:)
  type(c_ptr),   value,             intent(in)    :: obss

  type(State1dFM)                 :: x,x_tl
  type(Obs1dBangle)               :: y,y_tl
 
  integer                         :: iobs,nlev, nlocs
  integer                         :: nvprof
    
  character(len=*), parameter  :: myname_="ufo_gnssro_bndropp1d_simobs_tl"
  character(max_string)        :: err_msg
  real(kind=dp)                :: ob_time
  type(ufo_geoval), pointer    :: t_d, q_d, prs_d 

! hack - set local geopotential height to zero for ropp routines
  real(kind_real), allocatable :: gph_d_zero(:)
  real(kind_real)              :: gph_sfc_d_zero 
  real(kind_real), allocatable :: obsLat(:), obsLon(:), obsImpP(:), obsLocR(:), obsGeoid(:)
  integer                      :: use_compress

  use_compress = self%roconf%use_compress

! hack - set local geopotential height to zero for ropp routines

  write(err_msg,*) "TRACE: ufo_gnssro_bndropp1d_simobs_tl: begin"
  call fckit_log%debug(err_msg)

! check if trajectory was set
  if (.not. self%ltraj) then
     write(err_msg,*) myname_, ' trajectory wasnt set!'
     call abor1_ftn(err_msg)
  endif
      
! check if nlocs is consistent in geovals & hofx
  if (geovals%nlocs /= size(hofx)) then
     write(err_msg,*) myname_, ' error: nlocs inconsistent!'
     call abor1_ftn(err_msg)
  endif

! get variables from geovals
  call ufo_geovals_get_var(geovals, var_ts,    t_d)         ! temperature
  call ufo_geovals_get_var(geovals, var_q,     q_d)         ! specific humidity
  call ufo_geovals_get_var(geovals, var_prs,   prs_d)       ! pressure

  nlev  = self%nval 
  nlocs  = self%nlocs ! number of observations
  if (nlocs > 0) then
     allocate(gph_d_zero(nlev))
     gph_d_zero     = 0.0
     gph_sfc_d_zero = 0.0

   ! set obs space struture
     allocate(obsLon(nlocs))
     allocate(obsLat(nlocs))
     allocate(obsImpP(nlocs))
     allocate(obsLocR(nlocs))
     allocate(obsGeoid(nlocs))
     call obsspace_get_db(obss, "MetaData", "longitude",            obsLon)
     call obsspace_get_db(obss, "MetaData", "latitude",             obsLat)
     call obsspace_get_db(obss, "MetaData", "impactParameterRO",    obsImpP)
     call obsspace_get_db(obss, "MetaData", "earthRadiusCurvature", obsLocR)
     call obsspace_get_db(obss, "MetaData", "geoidUndulation",      obsGeoid)

     nvprof = 1  ! no. of bending angles in profile 

   ! loop through the obs
     obs_loop: do iobs = 1, nlocs   ! order of loop doesn't matter

       ob_time = 0.0
   !   map the trajectory to ROPP structure x
       call init_ropp_1d_statevec( ob_time,       &
                           obsLon(iobs),          &
                           obsLat(iobs),          &
                           self%t(:,iobs),        &
                           self%q(:,iobs),        &
                           self%prs(:,iobs),      &
                           self%gph(:,iobs),      &
                                      nlev,       &
                           self%gph_sfc(1,iobs),  &
                                 x, self%iflip, use_compress)
   !  hack -- make non zero humidity to avoid zero denominator in tangent linear
   !          see  ropp_fm/bangle_1d/ropp_fm_bangle_1d_tl.f90
       where(x%shum .le. 1e-8)        x%shum = 1e-8
   !  hack -- make non zero humidity to avoid zero denominator in tangent linear

       call init_ropp_1d_statevec( ob_time,        &
                           obsLon(iobs),           &
                           obsLat(iobs),           &
                           t_d%vals(:,iobs),       &
                           q_d%vals(:,iobs),       &
                           prs_d%vals(:,iobs),     &
                           gph_d_zero(:),          &
                                      nlev,        &
                           gph_sfc_d_zero,         &
                           x_tl, self%iflip, use_compress)
   !   set both y and y_tl structures    
       call init_ropp_1d_obvec_tlad(iobs, nvprof, &
                         obsImpP(iobs),           &
                         obsLat(iobs),            &
                         obsLon(iobs),            &
                         obsLocR(iobs),           &
                         obsGeoid(iobs),          &
                                y,y_tl)

   !   now call TL of forward model
       call ropp_fm_bangle_1d_tl(x,x_tl,y, y_tl%bangle(nvprof:nvprof))
       hofx(iobs) = y_tl%bangle(nvprof) ! this will need to change if profile is passed

   !   tidy up -deallocate ropp structures 
       call ropp_tidy_up_tlad_1d(x,x_tl,y,y_tl)
     end do obs_loop
 
   ! tidy up - deallocate obsspace structures
     deallocate(obsLat) 
     deallocate(obsLon)
     deallocate(obsImpP)
     deallocate(obsLocR)
     deallocate(obsGeoid)
     deallocate(gph_d_zero)

  end if ! nlocs > 0

  write(err_msg,*) "TRACE: ufo_gnssro_bndropp1d_simobs_tl: complete"
  call fckit_log%debug(err_msg)

end subroutine ufo_gnssro_bndropp1d_simobs_tl
 
! ------------------------------------------------------------------------------
! ------------------------------------------------------------------------------
subroutine ufo_gnssro_bndropp1d_simobs_ad(self, geovals, hofx, obss)

  use ropp_fm_types, only: State1dFM
  use ropp_fm_types, only: Obs1dBangle
  use typesizes,     only: wp => EightByteReal
  use datetimetypes, only: dp
  use ropp_fm, only: ropp_fm_bangle_1d_ad

  implicit none
  class(ufo_gnssro_BndROPP1D_tlad), intent(in)    :: self
  type(ufo_geovals),                intent(inout) :: geovals   ! perturbed quantities
  real(kind_real),                  intent(in)    :: hofx(:)
  type(c_ptr),  value,              intent(in)    :: obss
  real(c_double)              :: missing

  type(ufo_geoval),     pointer   :: t_d, q_d, prs_d 
! set local geopotential height to zero for ropp routines
  real(kind_real),      parameter :: gph_sfc_d_zero = 0.0
  real(kind_real),    allocatable :: gph_d_zero(:)

  real(kind_real),    allocatable :: obsLat(:), obsLon(:), obsImpP(:), obsLocR(:), obsGeoid(:)
  type(State1dFM)                 :: x,x_ad
  type(Obs1dBangle)               :: y,y_ad
  integer                         :: iobs,nlev, nlocs
  integer                         :: nvprof
  real(kind=dp)                   :: ob_time
  character(len=*), parameter     :: myname_="ufo_gnssro_bndropp1d_simobs_ad"
  character(max_string)           :: err_msg
  integer                         :: use_compress

  use_compress = self%roconf%use_compress

  write(err_msg,*) "TRACE: ufo_gnssro_bndropp1d_simobs_ad: begin"
  call fckit_log%debug(err_msg)
  if (self%nlocs > 0) then
   ! check if trajectory was set
     if (.not. self%ltraj) then
        write(err_msg,*) myname_, ' trajectory wasnt set!'
        call abor1_ftn(err_msg)
     endif
   ! check if nlocs is consistent in geovals & hofx
     if (geovals%nlocs /= size(hofx)) then
        write(err_msg,*) myname_, ' error: nlocs inconsistent!'
        call abor1_ftn(err_msg)
     endif
 
   ! get variables from geovals
     call ufo_geovals_get_var(geovals, var_ts,    t_d)         ! temperature
     call ufo_geovals_get_var(geovals, var_q,     q_d)         ! specific humidity
     call ufo_geovals_get_var(geovals, var_prs,   prs_d)       ! pressure

     nlev  = self%nval 
     nlocs  = self%nlocs

     allocate(gph_d_zero(nlev))
     gph_d_zero = 0.0

   ! set obs space struture
     allocate(obsLon(nlocs))
     allocate(obsLat(nlocs))
     allocate(obsImpP(nlocs))
     allocate(obsLocR(nlocs))
     allocate(obsGeoid(nlocs))

     call obsspace_get_db(obss, "MetaData", "longitude",            obsLon)
     call obsspace_get_db(obss, "MetaData", "latitude",             obsLat) 
     call obsspace_get_db(obss, "MetaData", "impactParameterRO",    obsImpP)
     call obsspace_get_db(obss, "MetaData", "earthRadiusCurvature", obsLocR)
     call obsspace_get_db(obss, "MetaData", "geoidUndulation",      obsGeoid)

     missing = missing_value(missing)

   ! loop through the obs
     nvprof=1  ! no. of bending angles in profile 
     obs_loop: do iobs = 1, nlocs 

       if (hofx(iobs) .gt. missing) then
           ob_time = 0.0

   !       map the trajectory to ROPP structure x
           call init_ropp_1d_statevec( ob_time,   &
                             obsLon(iobs),        &
                             obsLat(iobs),        &
                             self%t(:,iobs),      &
                             self%q(:,iobs),      &
                             self%prs(:,iobs),    &
                             self%gph(:,iobs),    &
                                        nlev,     &
                             self%gph_sfc(1,iobs),&
                                   x, self%iflip, use_compress)

           call init_ropp_1d_statevec( ob_time,  &
                               obsLon(iobs),     &
                               obsLat(iobs),     &
                           t_d%vals(:,iobs),     &
                           q_d%vals(:,iobs),     &
                          prs_d%vals(:,iobs),    &
                               gph_d_zero(:),    &
                                        nlev,    &
                              gph_sfc_d_zero,    &
                             x_ad, self%iflip, use_compress)

    !      x_ad is local so initialise to 0.0
           x_ad%temp(:) = 0.0_wp
           x_ad%pres(:) = 0.0_wp
           x_ad%shum(:) = 0.0_wp
           x_ad%geop(:) = 0.0_wp
 
    !      set both y and y_ad structures    
           call init_ropp_1d_obvec_tlad(iobs,  nvprof,  &
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
           call ropp_fm_bangle_1d_ad(x,x_ad,y,y_ad%bangle(nvprof:nvprof))
           call init_ropp_1d_statevec_ad(           &
                             t_d%vals(:,iobs),      &
                             q_d%vals(:,iobs),      &
                           prs_d%vals(:,iobs),      &
                           gph_d_zero(:),           &
                           nlev, x_ad, self%iflip) 

   !     tidy up - deallocate ropp structures  
         call ropp_tidy_up_tlad_1d(x,x_ad,y,y_ad)

       end if  ! end missing value check

     end do obs_loop

   ! tidy up - deallocate obsspace structures
     deallocate(obsLat) 
     deallocate(obsLon)
     deallocate(obsImpP)
     deallocate(obsLocR)
     deallocate(obsGeoid)
     deallocate(gph_d_zero)
  end if ! nlocs > 0

  write(err_msg,*) "TRACE: ufo_gnssro_bndropp1d_simobs_ad: complete"
  call fckit_log%debug(err_msg)

end subroutine ufo_gnssro_bndropp1d_simobs_ad
    
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
subroutine ufo_gnssro_bndropp1d_tlad_delete(self)

  implicit none
  class(ufo_gnssro_BndROPP1D_tlad), intent(inout) :: self
  character(len=*), parameter :: myname_="ufo_gnssro_bndropp_tlad_delete"
      
  self%nval = 0
  self%nlocs = 0
  if (allocated(self%prs)) deallocate(self%prs)
  if (allocated(self%t))   deallocate(self%t)
  if (allocated(self%q))   deallocate(self%q)
  if (allocated(self%gph)) deallocate(self%gph)
  if (allocated(self%gph_sfc)) deallocate(self%gph_sfc)
  self%ltraj = .false. 

end subroutine ufo_gnssro_bndropp1d_tlad_delete

!-------------------------------------------------------------------------
       
end module ufo_gnssro_bndropp1d_tlad_mod
