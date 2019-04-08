! (C) Copyright 2017-2018 UCAR
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 

!> Stubbed Fortran module for gnssro bending angle ropp2d tangent linear and adjoint
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
use fckit_log_module, only : fckit_log

integer, parameter         :: max_string=800

!> Fortran derived type for gnssro trajectory
type, extends(ufo_basis_tlad)   ::  ufo_gnssro_BndROPP2D_tlad
  private
  integer                       :: nval, nobs
  real(kind_real), allocatable  :: prs(:,:), t(:,:), q(:,:), gph(:,:), gph_sfc(:,:)
  integer                       :: n_horiz      ! 2d points along ray path
  integer                       :: iflip        ! geoval ascending order flag
  type(gnssro_conf)             :: roconf       ! ro configuration
  real(kind_real), allocatable  :: obsLon2d(:), obsLat2d(:)  !2d locations - nobs*n_horiz
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
  integer                     :: iobs

  write(err_msg,*) "TRACE: ufo_gnssro_bndropp2d_tlad_settraj: begin"
  call fckit_log%info(err_msg) 

! get model state variables from geovals
  call ufo_geovals_get_var(geovals, var_t,     t)         ! temperature
  call ufo_geovals_get_var(geovals, var_q,     q)         ! specific humidity
  call ufo_geovals_get_var(geovals, var_prs,   prs)       ! pressure
  call ufo_geovals_get_var(geovals, var_z,     gph)       ! geopotential height
      
  call self%delete()   

! Keep copy of dimensions
  self%nval = prs%nval
  self%nobs = obsspace_get_nlocs(obss)

  allocate(self%t(self%nval,self%nobs))
  allocate(self%q(self%nval,self%nobs))
  allocate(self%prs(self%nval,self%nobs))
  allocate(self%gph(self%nval,self%nobs))

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

  implicit none
  class(ufo_gnssro_BndROPP2D_tlad), intent(in)    :: self
  type(ufo_geovals),                intent(in)    :: geovals   ! perturbed quantities
  real(kind_real),                  intent(inout) :: hofx(:)
  type(c_ptr),   value,             intent(in)    :: obss

  integer                         :: iobs,nlev, nobs, nvprof
    
  character(len=*), parameter  :: myname_="ufo_gnssro_bndropp2d_simobs_tl"
  character(max_string)        :: err_msg
  type(ufo_geoval), pointer    :: t_d, q_d, prs_d 

! hack - set local geopotential height to zero for ropp routines
  real(kind_real), allocatable :: gph_d_zero(:)
  real(kind_real)              :: gph_sfc_d_zero 
  real(kind_real), allocatable :: obsLat(:), obsLon(:), obsImpP(:), obsLocR(:), obsGeoid(:)
! hack - set local geopotential height to zero for ropp routines

  write(err_msg,*) "TRACE: ufo_gnssro_bndropp2d_simobs_tl: begin"
  call fckit_log%info(err_msg)

! check if trajectory was set
  if (.not. self%ltraj) then
     write(err_msg,*) myname_, ' trajectory wasnt set!'
     call abor1_ftn(err_msg)
  endif
      
! check if nobs is consistent in geovals & hofx
  if (geovals%nobs /= size(hofx)) then
     write(err_msg,*) myname_, ' error: nobs inconsistent!'
     call abor1_ftn(err_msg)
  endif

! get variables from geovals
  call ufo_geovals_get_var(geovals, var_t,     t_d)         ! temperature
  call ufo_geovals_get_var(geovals, var_q,     q_d)         ! specific humidity
  call ufo_geovals_get_var(geovals, var_prs,   prs_d)       ! pressure

  nlev  = self%nval 
  nobs  = self%nobs ! number of observations

  allocate(gph_d_zero(nlev))
  gph_d_zero     = 0.0
  gph_sfc_d_zero = 0.0

! set obs space struture
  allocate(obsLon(nobs))
  allocate(obsLat(nobs))
  allocate(obsImpP(nobs))
  allocate(obsLocR(nobs))
  allocate(obsGeoid(nobs))
  call obsspace_get_db(obss, " ", "longitude",        obsLon)
  call obsspace_get_db(obss, " ", "latitude",         obsLat)
  call obsspace_get_db(obss, " ", "impact_parameter", obsImpP)
  call obsspace_get_db(obss, " ", "earth_radius_of_curvature", obsLocR)
  call obsspace_get_db(obss, " ", "geoid_height_above_reference_ellipsoid", obsGeoid) 

  nvprof = 1  ! no. of bending angles in profile 

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

  implicit none
  class(ufo_gnssro_BndROPP2D_tlad), intent(in)    :: self
  type(ufo_geovals),                intent(inout) :: geovals   ! perturbed quantities
  real(kind_real),                  intent(in)    :: hofx(:)
  type(c_ptr),  value,              intent(in)    :: obss
  real(c_double)              :: missing

  type(ufo_geoval),     pointer   :: t_d, q_d, prs_d 
! set local geopotential height to zero for ropp routines
  real(kind_real),      parameter :: gph_sfc_d_zero = 0.0
  real(kind_real),    allocatable :: gph_d_zero(:)

  real(kind_real),    allocatable :: obsLat(:), obsLon(:), obsImpP(:), obsLocR(:), obsGeoid(:)
  integer                         :: iobs,nlev, nobs, nvprof
  character(len=*), parameter     :: myname_="ufo_gnssro_bndropp2d_simobs_ad"
  character(max_string)           :: err_msg

  write(err_msg,*) "TRACE: ufo_gnssro_bndropp2d_simobs_ad: begin"
  call fckit_log%info(err_msg)

! check if trajectory was set
  if (.not. self%ltraj) then
     write(err_msg,*) myname_, ' trajectory wasnt set!'
     call abor1_ftn(err_msg)
  endif
! check if nobs is consistent in geovals & hofx
  if (geovals%nobs /= size(hofx)) then
     write(err_msg,*) myname_, ' error: nobs inconsistent!'
     call abor1_ftn(err_msg)
  endif
     
! get variables from geovals
  call ufo_geovals_get_var(geovals, var_t,     t_d)         ! temperature
  call ufo_geovals_get_var(geovals, var_q,     q_d)         ! specific humidity
  call ufo_geovals_get_var(geovals, var_prs,   prs_d)       ! pressure

! allocate if not yet allocated   
  if (.not. allocated(t_d%vals)) then
      t_d%nobs = self%nobs
      t_d%nval = self%nval
      allocate(t_d%vals(t_d%nval,t_d%nobs))
      t_d%vals = 0.0_kind_real
  endif

  if (.not. allocated(prs_d%vals)) then
      prs_d%nobs = self%nobs
      prs_d%nval = self%nval
      allocate(prs_d%vals(prs_d%nval,prs_d%nobs))
      prs_d%vals = 0.0_kind_real
  endif

  if (.not. allocated(q_d%vals)) then
      q_d%nobs = self%nobs
      q_d%nval = self%nval
      allocate(q_d%vals(q_d%nval,q_d%nobs))
      q_d%vals = 0.0_kind_real
  endif

  if (.not. geovals%linit ) geovals%linit=.true.

  nlev  = self%nval 
  nobs  = self%nobs

  allocate(gph_d_zero(nlev))
  gph_d_zero = 0.0

! set obs space struture
  allocate(obsLon(nobs))
  allocate(obsLat(nobs))
  allocate(obsImpP(nobs))
  allocate(obsLocR(nobs))
  allocate(obsGeoid(nobs))

  call obsspace_get_db(obss, " ", "longitude", obsLon)
  call obsspace_get_db(obss, " ", "latitude", obsLat) 
  call obsspace_get_db(obss, " ", "impact_parameter", obsImpP)
  call obsspace_get_db(obss, " ", "earth_radius_of_curvature", obsLocR)
  call obsspace_get_db(obss, " ", "geoid_height_above_reference_ellipsoid", obsGeoid)
  
  missing = missing_value(missing)

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
!  if (allocated(self%gph_sfc)) deallocate(self%gph_sfc)

  self%ltraj = .false. 

end subroutine ufo_gnssro_bndropp2d_tlad_delete

!-------------------------------------------------------------------------
       
end module ufo_gnssro_bndropp2d_tlad_mod
