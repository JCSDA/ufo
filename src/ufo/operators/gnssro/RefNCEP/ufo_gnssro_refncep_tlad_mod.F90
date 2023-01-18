! (C) Copyright 2017-2018 UCAR
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 

!> Fortran module for gnssro refractivity tangent linear and adjoint
!> following NCEP's method

module ufo_gnssro_refncep_tlad_mod
  use fckit_configuration_module, only: fckit_configuration 
  use kinds
  use ufo_vars_mod
  use ufo_geovals_mod
  use ufo_geovals_mod_c,   only: ufo_geovals_registry
  use vert_interp_mod
  use ufo_basis_tlad_mod,  only: ufo_basis_tlad
  use obsspace_mod
  use gnssro_mod_constants
  use gnssro_mod_conf
  use missing_values_mod
  use ufo_constants_mod
  use iso_c_binding

  !> Fortran derived type for gnssro trajectory
  type, extends(ufo_basis_tlad) :: ufo_gnssro_RefNCEP_tlad
   private
     type(gnssro_conf)             :: roconf
     integer                       :: nval, nlocs
     real(kind_real), allocatable  :: wf(:)
     integer,         allocatable  :: wi(:)
     real(kind_real), allocatable  :: jac_t(:), jac_q(:), jac_prs(:)

  contains
    procedure :: setup      => ufo_gnssro_refncep_tlad_setup
    procedure :: delete     => ufo_gnssro_refncep_tlad_delete
    procedure :: settraj    => ufo_gnssro_refncep_tlad_settraj
    procedure :: simobs_tl  => ufo_gnssro_refncep_simobs_tl
    procedure :: simobs_ad  => ufo_gnssro_refncep_simobs_ad
  end type ufo_gnssro_RefNCEP_tlad

contains
! ------------------------------------------------------------------------------
subroutine ufo_gnssro_refncep_tlad_setup(self, f_conf)
  implicit none
  class(ufo_gnssro_RefNCEP_tlad), intent(inout) :: self
  type(fckit_configuration),     intent(in)     :: f_conf

  call gnssro_conf_setup(self%roconf,f_conf)

end subroutine ufo_gnssro_refncep_tlad_setup

! ------------------------------------------------------------------------------
subroutine ufo_gnssro_refncep_tlad_settraj(self, geovals, obss)
  use gnssro_mod_transform, only: geometric2geop
      
  implicit none
  class(ufo_gnssro_RefNCEP_tlad),intent(inout) :: self
  type(ufo_geovals),             intent(in)    :: geovals
  type(c_ptr), value,            intent(in)    :: obss
      
  character(len=*), parameter :: myname_="ufo_gnssro_refncep_tlad_settraj"
      
  type(ufo_geoval),    pointer :: t,q,prs,gph
  real(kind_real), allocatable :: obsZ(:), obsLat(:)  ! observation vector
  real(kind_real)  :: obsH,gesT,gesQ,gesP
  real(kind_real)  :: Tv, Tv0
  integer          :: wi0, iobs

! Get variables from geovals
  call ufo_geovals_get_var(geovals, var_prs, prs)
  call ufo_geovals_get_var(geovals, var_ts,t)
  call ufo_geovals_get_var(geovals, var_q, q)
  call ufo_geovals_get_var(geovals, var_z, gph)

! Make sure nothing already allocated
  call self%delete()

! Keep copy of dimensions
  self%nval  = prs%nval
  self%nlocs = obsspace_get_nlocs(obss)
 
  allocate(self%wi(self%nlocs))
  allocate(self%wf(self%nlocs))
  allocate(self%jac_t(self%nlocs))
  allocate(self%jac_q(self%nlocs))
  allocate(self%jac_prs(self%nlocs))
  allocate(obsZ(self%nlocs))
  allocate(obsLat(self%nlocs))

! get observation vectors
  call obsspace_get_db(obss, "MetaData", "height", obsZ)
  call obsspace_get_db(obss, "MetaData", "latitude", obsLat)
  call gnssro_ref_constants(self%roconf%use_compress)

  do iobs = 1, self%nlocs

!   calculate observation geopotential height using  MJ Mahoney's (2001)
    call geometric2geop(obsLat(iobs), obsZ(iobs), obsH)
    call vert_interp_weights(self%nval, obsH, gph%vals(:,iobs),self%wi(iobs),self%wf(iobs))
    wi0 = self%wi(iobs)
    call vert_interp_apply(t%nval, t%vals(:,iobs), gesT, self%wi(iobs),self%wf(iobs))
    call vert_interp_apply(q%nval, q%vals(:,iobs), gesQ, self%wi(iobs),self%wf(iobs))

!   use  hypsometric equation to calculate pressure
    Tv  = 0.0
    Tv0 = 0.0
    Tv  = gesT*( one + (rv_over_rd-one)*gesQ/(1.0-gesQ) )
    Tv0 = t%vals(wi0,iobs)*(one + (rv_over_rd-one)*q%vals(wi0,iobs)/(1.0-q%vals(wi0,iobs) ))
    gesP = prs%vals(wi0,iobs)/exp(two*grav*(obsH-gph%vals(wi0,iobs))/(rd*(Tv+Tv0)))

    self%jac_t(iobs)   = - n_a*gesP/gesT**2                                                 &
                         - n_b*two*gesP*gesQ/( ((1-rd_over_rv)*gesQ+rd_over_rv)*gesT**3  )  &
                         - n_c*gesP*gesQ/( ((1-rd_over_rv)*gesQ+rd_over_rv)*gesT**2 )

    self%jac_q(iobs)   =   n_b*gesP/( gesT**2*( (1-rd_over_rv)*gesQ+rd_over_rv)**2 )         &
                                   * rd_over_rv                                              &
                         + n_c*gesP/( gesT   *( (1-rd_over_rv)*gesQ+rd_over_rv)**2 )         &
                                   * rd_over_rv
    self%jac_prs(iobs) =   n_a/gesT                                                 &
                         + n_b*gesQ/ ( ((1-rd_over_rv)*gesQ+rd_over_rv)*gesT**2 )   &
                         + n_c*gesQ/ ( ((1-rd_over_rv)*gesQ+rd_over_rv)*gesT )

  enddo

  
  self%ltraj = .true.
! cleanup
  deallocate(obsZ)
  deallocate(obsLat)

end subroutine ufo_gnssro_refncep_tlad_settraj
    
! ------------------------------------------------------------------------------
    
subroutine ufo_gnssro_refncep_simobs_tl(self, geovals, hofx, obss)
  implicit none
  class(ufo_gnssro_RefNCEP_tlad), intent(in)     :: self
  type(ufo_geovals),              intent(in)     :: geovals
  real(kind_real),                intent(inout)  :: hofx(:)
  type(c_ptr), value,             intent(in)     :: obss
     
  character(len=*), parameter :: myname="ufo_gnssro_refncep_tlad_tl"
  character(max_string)       :: err_msg
      
  type(ufo_geoval), pointer :: t_tl, q_tl, prs_tl
  real(kind_real)           :: gesT_tl, gesQ_tl, gesP_tl
  integer                   :: iobs

! check if trajectory was set
  if (.not. self%ltraj) then
      write(err_msg,*) myname, ' trajectory wasnt set!'
      call abor1_ftn(err_msg)
  endif
      
! check if nlocs is consistent in geovals & hofx
  if (geovals%nlocs /= size(hofx)) then
    write(err_msg,*) myname, ' error: nlocs inconsistent!'
    call abor1_ftn(err_msg)
  endif
     
! get variables from geovals
  call ufo_geovals_get_var(geovals, var_ts,t_tl)
  call ufo_geovals_get_var(geovals, var_q, q_tl)
  call ufo_geovals_get_var(geovals, var_prs, prs_tl)
 
! tangent linear obs operator (linear)
  do iobs = 1, geovals%nlocs
     call vert_interp_apply_tl(  t_tl%nval,  t_tl%vals(:,iobs), gesT_tl, self%wi(iobs),self%wf(iobs))
     call vert_interp_apply_tl(  q_tl%nval,  q_tl%vals(:,iobs), gesQ_tl, self%wi(iobs),self%wf(iobs))
     call vert_interp_apply_tl(prs_tl%nval,prs_tl%vals(:,iobs), gesP_tl, self%wi(iobs),self%wf(iobs))
     hofx(iobs)  =  self%jac_t(iobs)*gesT_tl  + self%jac_q(iobs)*gesQ_tl + self%jac_prs(iobs)*gesP_tl
  enddo
    
end subroutine ufo_gnssro_refncep_simobs_tl
! ------------------------------------------------------------------------------
    
subroutine ufo_gnssro_refncep_simobs_ad(self, geovals, hofx, obss)
  implicit none
  class(ufo_gnssro_RefNCEP_tlad), intent(in)    :: self
  type(ufo_geovals),              intent(inout) :: geovals
  real(kind_real),                intent(in)    :: hofx(:)
  type(c_ptr), value,             intent(in)    :: obss

  character(len=*), parameter :: myname="ufo_gnssro_refncep_tlad_ad"
  character(max_string)       :: err_msg
  type(ufo_geoval), pointer :: t_d, q_d, prs_d
  real(kind_real)           :: gesT_d, gesQ_d, gesP_d
  real(c_double)            :: missing
  integer                   :: iobs

! check if trajectory was set
  if (.not. self%ltraj) then
     write(err_msg,*) myname, ' trajectory wasnt set!'
     call abor1_ftn(err_msg)
  endif
      
! check if nlocs is consistent in geovals & hofx
  if (geovals%nlocs /= size(hofx)) then
     write(err_msg,*) myname, ' error: nlocs inconsistent!'
     call abor1_ftn(err_msg)
  endif
      
! get variables from geovals
  call ufo_geovals_get_var(geovals, var_prs, prs_d)
  call ufo_geovals_get_var(geovals, var_ts,t_d)
  call ufo_geovals_get_var(geovals, var_q, q_d)

  missing = missing_value(missing)

  do iobs = 1, geovals%nlocs
    
    if (hofx(iobs) .ne. missing) then
      gesT_d = 0.0_kind_real
      gesQ_d = 0.0_kind_real
      gesP_d = 0.0_kind_real
      gesT_d = gesT_d + hofx(iobs)*self%jac_t(iobs)
      gesQ_d = gesQ_d + hofx(iobs)*self%jac_q(iobs)
      gesP_d = gesP_d + hofx(iobs)*self%jac_prs(iobs)
      call vert_interp_apply_ad(  t_d%nval,  t_d%vals(:,iobs), gesT_d, self%wi(iobs), self%wf(iobs))
      call vert_interp_apply_ad(  q_d%nval,  q_d%vals(:,iobs), gesQ_d, self%wi(iobs), self%wf(iobs))
      call vert_interp_apply_ad(prs_d%nval,prs_d%vals(:,iobs), gesP_d, self%wi(iobs), self%wf(iobs))
    endif

  enddo

end subroutine ufo_gnssro_refncep_simobs_ad
! ------------------------------------------------------------------------------

subroutine ufo_gnssro_refncep_tlad_delete(self)
  implicit none
  class(ufo_gnssro_RefNCEP_tlad), intent(inout) :: self
  character(len=*), parameter :: myname_="ufo_gnssro_refncep_tlad_delete"
      
  self%nval = 0
  if (allocated(self%wi))    deallocate(self%wi)
  if (allocated(self%wf))    deallocate(self%wf)
  if (allocated(self%jac_t)) deallocate(self%jac_t)
  if (allocated(self%jac_q)) deallocate(self%jac_q)
  if (allocated(self%jac_prs)) deallocate(self%jac_prs)
  self%ltraj = .false.
end subroutine ufo_gnssro_refncep_tlad_delete
    
! ------------------------------------------------------------------------------

end module ufo_gnssro_refncep_tlad_mod
