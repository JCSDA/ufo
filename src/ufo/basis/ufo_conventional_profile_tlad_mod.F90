! (C) Copyright 2017-2018 UCAR
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 

module ufo_conventional_profile_tlad_mod

  use iso_c_binding
  use kinds
  use ufo_vars_mod
  use ufo_geovals_mod
  use ufo_geovals_mod_c,   only: ufo_geovals_registry
  use ioda_obs_vectors
  use ioda_obsdb_mod
  use ioda_obsdb_mod_c, only: ioda_obsdb_registry
  use vert_interp_mod
  use ufo_basis_tlad_mod, only: ufo_basis_tlad

  integer, parameter :: max_string=800

  type, extends(ufo_basis_tlad) :: ufo_conventional_profile_tlad
   private
     integer :: nval, nobs
     real(kind_real), allocatable :: wf(:)
     integer, allocatable :: wi(:)
  contains
    procedure :: delete => conventional_profile_tlad_delete_
    procedure :: settraj => conventional_profile_tlad_settraj_
    procedure :: eqv_tl => conventional_profile_tlad_t_eqv_tl_
    procedure :: eqv_ad => conventional_profile_tlad_t_eqv_ad_
  end type ufo_conventional_profile_tlad
contains

! ------------------------------------------------------------------------------
    
    subroutine conventional_profile_tlad_settraj_(self, geovals, obss)
      implicit none
      class(ufo_conventional_profile_tlad), intent(inout) :: self
      type(ufo_geovals),         intent(in)    :: geovals
      type(ioda_obsdb),          intent(in)    :: obss
      
      character(len=*), parameter :: myname_="ufo_conventional_profile_tlad_settraj"
      character(max_string) :: err_msg
      
      type(obs_vector) :: pressure
      type(ufo_geoval), pointer :: prsl
      integer :: iobs, ierr
      
      !Check if conventional_profiles in geovals and get it
      call ufo_geovals_get_var(geovals, var_prsl, prsl, status=ierr)
      if (ierr/=0) then
        write(err_msg,*) myname_, trim(var_prsl), ' doesnt exist'
        call abor1_ftn(err_msg)
      endif
      
      !Make sure nothing already allocated
      call self%delete()
      
      !Keep copy of dimensions
      self%nval = prsl%nval
      self%nobs = obss%nobs
      
      allocate(self%wi(obss%nobs))
      allocate(self%wf(obss%nobs))
      
      ! observation of pressure (for vertical interpolation)
      call ioda_obsvec_setup(pressure, obss%nobs)
      call ioda_obsdb_var_to_ovec(obss, pressure, "Pressure")
      
      ! compute interpolation weights
      do iobs = 1, obss%nobs
        call vert_interp_weights(self%nval,log(pressure%values(iobs)/10.),prsl%vals(:,iobs),self%wi(iobs),self%wf(iobs))
      enddo
      
      self%ltraj = .true.
      ! cleanup
      call ioda_obsvec_delete(pressure)
    
    end subroutine conventional_profile_tlad_settraj_
    
! ------------------------------------------------------------------------------
    
    subroutine conventional_profile_tlad_t_eqv_tl_(self, geovals, hofx, obss)
      implicit none
      class(ufo_conventional_profile_tlad), intent(in)     :: self
      type(ufo_geovals),         intent(in)     :: geovals
      type(obs_vector),          intent(inout)  :: hofx
      type(ioda_obsdb),          intent(in)     :: obss
      
      character(len=*), parameter :: myname_="ufo_conventional_profile_tlad_t_eqv_tl"
      character(max_string) :: err_msg
      
      integer :: iobs,ierr
      type(ufo_geoval), pointer :: tv_d
      
      ! check if trajectory was set
      if (.not. self%ltraj) then
        write(err_msg,*) myname_, ' trajectory wasnt set!'
        call abor1_ftn(err_msg)
      endif
      
      ! check if nobs is consistent in geovals & hofx
      if (geovals%nobs /= hofx%nobs) then
        write(err_msg,*) myname_, ' error: nobs inconsistent!'
        call abor1_ftn(err_msg)
      endif
      
      ! check if tv variable is in geovals and get it
      call ufo_geovals_get_var(geovals, var_tv, tv_d, status=ierr )
      if (ierr/=0) then
        write(err_msg,*) myname_, trim(var_tv), ' doesnt exist'
        call abor1_ftn(err_msg)
      endif
      
      ! tangent linear obs operator (linear)
      do iobs = 1, hofx%nobs
        call vert_interp_apply_tl(tv_d%nval, tv_d%vals(:,iobs), hofx%values(iobs), self%wi(iobs), self%wf(iobs))
      enddo
    
    end subroutine conventional_profile_tlad_t_eqv_tl_
    
! ------------------------------------------------------------------------------
    
    subroutine conventional_profile_tlad_t_eqv_ad_(self, geovals, hofx, obss)
      implicit none
      class(ufo_conventional_profile_tlad), intent(in)     :: self
      type(ufo_geovals),         intent(inout)  :: geovals
      type(obs_vector),          intent(in)     :: hofx
      type(ioda_obsdb),          intent(in)     :: obss
      
      character(len=*), parameter :: myname_="ufo_conventional_profile_tlad_t_eqv_ad"
      character(max_string) :: err_msg
      
      integer :: iobs,ierr
      type(ufo_geoval), pointer :: tv_d, prsl_d
      
      ! check if trajectory was set
      if (.not. self%ltraj) then
        write(err_msg,*) myname_, ' trajectory wasnt set!'
        call abor1_ftn(err_msg)
      endif
      
      ! check if nobs is consistent in geovals & hofx
      if (geovals%nobs /= hofx%nobs) then
        write(err_msg,*) myname_, ' error: nobs inconsistent!'
        call abor1_ftn(err_msg)
      endif
      
      ! check if tv variable is in geovals and get it
      call ufo_geovals_get_var(geovals, var_prsl, prsl_d, status=ierr)
      if (ierr/=0) then
        write(err_msg,*) myname_, trim(var_prsl), ' doesnt exist'
        call abor1_ftn(err_msg)
      endif
      
      ! check if tv variable is in geovals and get it
      call ufo_geovals_get_var(geovals, var_tv, tv_d, status=ierr)
      if (ierr/=0) then
        write(err_msg,*) myname_, trim(var_tv), ' doesnt exist'
        call abor1_ftn(err_msg)
      endif
     
      ! allocate if not yet allocated	
      if (.not. allocated(tv_d%vals)) then
         tv_d%nobs = self%nobs
         tv_d%nval = self%nval
         allocate(tv_d%vals(tv_d%nval,tv_d%nobs))
         tv_d%vals = 0.0_kind_real
      endif
      if (.not. allocated(prsl_d%vals)) then
         prsl_d%nobs = self%nobs
         prsl_d%nval = self%nval
         allocate(prsl_d%vals(prsl_d%nval,prsl_d%nobs))
         prsl_d%vals = 0.0_kind_real
      endif
      if (.not. geovals%linit ) geovals%linit=.true.
 
      do iobs = 1, hofx%nobs
        call vert_interp_apply_ad(tv_d%nval, tv_d%vals(:,iobs), hofx%values(iobs), self%wi(iobs), self%wf(iobs))
      enddo
    
    end subroutine conventional_profile_tlad_t_eqv_ad_
    
! ------------------------------------------------------------------------------
    
    subroutine conventional_profile_tlad_delete_(self)
      implicit none
      class(ufo_conventional_profile_tlad), intent(inout) :: self
      
      character(len=*), parameter :: myname_="ufo_conventional_profile_tlad_delete"
      
      self%nval = 0
      if (allocated(self%wi)) deallocate(self%wi)
      if (allocated(self%wf)) deallocate(self%wf)
      self%ltraj = .false.
    
    end subroutine conventional_profile_tlad_delete_
    
! ------------------------------------------------------------------------------

end module ufo_conventional_profile_tlad_mod
