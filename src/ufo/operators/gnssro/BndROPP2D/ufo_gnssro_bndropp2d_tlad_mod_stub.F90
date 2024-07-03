! (C) Copyright 2017-2018 UCAR
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 

!> Stubbed Fortran module for gnssro bending angle ropp2d tangent linear and adjoint
!> following the ROPP (2018 Aug) implementation

module ufo_gnssro_bndropp2d_tlad_mod

use fckit_configuration_module, only: fckit_configuration 
use kinds
use ufo_vars_mod
use ufo_geovals_mod
use ufo_geovals_mod_c,   only: ufo_geovals_registry
use vert_interp_mod
use ufo_basis_tlad_mod,  only: ufo_basis_tlad
use obsspace_mod
use gnssro_mod_conf
use missing_values_mod
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
subroutine ufo_gnssro_bndropp2d_tlad_setup(self, f_conf)
  implicit none
  class(ufo_gnssro_BndROPP2D_tlad), intent(inout) :: self
  type(fckit_configuration), intent(in) :: f_conf

  call gnssro_conf_setup(self%roconf,f_conf)

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
  call fckit_log%debug(err_msg) 

! get model state variables from geovals
  call ufo_geovals_get_var(geovals, var_ts,    t)         ! temperature
  call ufo_geovals_get_var(geovals, var_q,     q)         ! specific humidity
  call ufo_geovals_get_var(geovals, var_prs,   prs)       ! pressure
  call ufo_geovals_get_var(geovals, var_z,     gph)       ! geopotential height
      
  call self%delete()   

! Keep copy of dimensions
  self%nval = prs%nval
  self%nlocs = obsspace_get_nlocs(obss)

  allocate(self%t(self%nval,self%nlocs))
  allocate(self%q(self%nval,self%nlocs))
  allocate(self%prs(self%nval,self%nlocs))
  allocate(self%gph(self%nval,self%nlocs))

! allocate  
  self%gph     = gph%vals
  self%t       = t%vals
  self%q       = q%vals
  self%prs     = prs%vals

  self%ltraj   = .true.

  write(err_msg,*) "TRACE: ufo_gnssro_bndropp2d_tlad_settraj: complete"
  call fckit_log%debug(err_msg)
    
end subroutine ufo_gnssro_bndropp2d_tlad_settraj
    
! ------------------------------------------------------------------------------
! ------------------------------------------------------------------------------    
subroutine ufo_gnssro_bndropp2d_simobs_tl(self, geovals, hofx, obss)

  implicit none
  class(ufo_gnssro_BndROPP2D_tlad), intent(in)    :: self
  type(ufo_geovals),                intent(in)    :: geovals   ! perturbed quantities
  real(kind_real),                  intent(inout) :: hofx(:)
  type(c_ptr),   value,             intent(in)    :: obss

  integer                         :: iobs,nlev, nlocs, nvprof
    
  character(len=*), parameter  :: myname_="ufo_gnssro_bndropp2d_simobs_tl"
  character(max_string)        :: err_msg
  type(ufo_geoval), pointer    :: t_d, q_d, prs_d 

! hack - set local geopotential height to zero for ropp routines
  real(kind_real), allocatable :: gph_d_zero(:)
  real(kind_real)              :: gph_sfc_d_zero 
  real(kind_real), allocatable :: obsLat(:), obsLon(:), obsImpP(:), obsLocR(:), obsGeoid(:)
! hack - set local geopotential height to zero for ropp routines
  integer                       :: n_horiz
  integer                       :: use_compress

  n_horiz = self%roconf%n_horiz
  use_compress = self%roconf%use_compress

  write(err_msg,*) "TRACE: ufo_gnssro_bndropp2d_simobs_tl: begin"
  call fckit_log%debug(err_msg)

! check if trajectory was set
  if (.not. self%ltraj) then
     write(err_msg,*) myname_, ' trajectory wasnt set!'
     call abor1_ftn(err_msg)
  endif

! get variables from geovals
  call ufo_geovals_get_var(geovals, var_ts,    t_d)         ! temperature
  call ufo_geovals_get_var(geovals, var_q,     q_d)         ! specific humidity
  call ufo_geovals_get_var(geovals, var_prs,   prs_d)       ! pressure

! check if the number of geoval profiles is correct
  if (t_d%nprofiles /= size(hofx)*n_horiz .or. q_d%nprofiles /= size(hofx)*n_horiz .or. &
      prs_d%nprofiles /= size(hofx)*n_horiz) then
     write(err_msg,*) myname_, ' error: npaths inconsistent!'
     call abor1_ftn(err_msg)
  endif

  nlev  = self%nval 
  nlocs  = self%nlocs ! number of observations

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

! tidy up - deallocate obsspace structures
  deallocate(obsLat) 
  deallocate(obsLon)
  deallocate(obsImpP)
  deallocate(obsLocR)
  deallocate(obsGeoid)

  write(err_msg,*) "TRACE: ufo_gnssro_bndropp2d_simobs_tl: complete"
  call fckit_log%debug(err_msg)

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
  integer                         :: iobs,nlev, nlocs, nvprof
  integer                         :: n_horiz
  character(len=*), parameter     :: myname_="ufo_gnssro_bndropp2d_simobs_ad"
  character(max_string)           :: err_msg


  write(err_msg,*) "TRACE: ufo_gnssro_bndropp2d_simobs_ad: begin"
  call fckit_log%debug(err_msg)

! check if trajectory was set
  if (.not. self%ltraj) then
     write(err_msg,*) myname_, ' trajectory wasnt set!'
     call abor1_ftn(err_msg)
  endif

  n_horiz = self%roconf%n_horiz

! get variables from geovals
  call ufo_geovals_get_var(geovals, var_ts,    t_d)         ! temperature
  call ufo_geovals_get_var(geovals, var_q,     q_d)         ! specific humidity
  call ufo_geovals_get_var(geovals, var_prs,   prs_d)       ! pressure

! check if the number of geoval profiles is correct
  if (t_d%nprofiles /= size(hofx)*n_horiz .or. q_d%nprofiles /= size(hofx)*n_horiz .or. &
      prs_d%nprofiles /= size(hofx)*n_horiz) then
     write(err_msg,*) myname_, ' error: npaths inconsistent!'
     call abor1_ftn(err_msg)
  endif

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

! tidy up - deallocate obsspace structures
  deallocate(obsLat) 
  deallocate(obsLon)
  deallocate(obsImpP)
  deallocate(obsLocR)
  deallocate(obsGeoid)
  deallocate(gph_d_zero)

  write(err_msg,*) "TRACE: ufo_gnssro_bndropp2d_simobs_ad: complete"
  call fckit_log%debug(err_msg)

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
