! (C) Copyright 2017-2018 UCAR
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 

!> Fortran module for gnssro refractivity tangent linear and adjoint

module ufo_gnssro_ref_tlad_mod
  use iso_c_binding
  use kinds
  use ufo_vars_mod
  use ufo_geovals_mod
  use ufo_geovals_mod_c,   only: ufo_geovals_registry
  use ioda_obs_vectors
  use ioda_obsdb_mod
  use ioda_obsdb_mod_c,    only: ioda_obsdb_registry
  use vert_interp_mod
  use ufo_basis_tlad_mod,  only: ufo_basis_tlad

  integer, parameter :: max_string=800

  !> Fortran derived type for gnssro trajectory
  type, extends(ufo_basis_tlad) :: ufo_gnssro_tlad
   private
     integer                       :: nval, nobs
     real(kind_real), allocatable  :: wf(:)
     integer,         allocatable  :: wi(:)
     real(kind_real), allocatable  :: prs_traj(:), t_traj(:), mixr_traj(:), gph_traj(:,:)
     real(kind_real), allocatable  :: obsH(:)
  contains
    procedure :: delete  => ufo_gnssro_ref_tlad_delete
    procedure :: settraj => ufo_gnssro_ref_tlad_settraj
    procedure :: eqv_tl  => ufo_gnssro_ref_tlad_tl
    procedure :: eqv_ad  => ufo_gnssro_ref_tlad_ad
  end type ufo_gnssro_tlad

contains
! ------------------------------------------------------------------------------
    
    subroutine ufo_gnssro_ref_tlad_settraj(self, geovals, obss)
      use gnssro_mod_constants
      use gnssro_mod_transform, only: geometric2geop
      
      implicit none
      class(ufo_gnssro_tlad), intent(inout)   :: self
      type(ufo_geovals),         intent(in)   :: geovals
      type(ioda_obsdb),          intent(in)   :: obss
      
      character(len=*), parameter :: myname_="ufo_gnssro_ref_tlad_settraj"
      character(max_string)       :: err_msg
      
      type(ufo_geoval), pointer :: t,mixr,prs,gph
      integer          :: iobs, ierr
      type(obs_vector) :: obsZ, obsLat  ! observation vector
      real(kind_real)  :: Tv_traj, Tv_traj0
      real(kind_real)  :: wf0
      integer          :: wi0

      !Check if variables are in geovals and get them
      call ufo_geovals_get_var(geovals, var_prs, prs, status=ierr)
      if (ierr/=0) then
        write(err_msg,*) myname_, trim(var_prs), ' doesnt exist'
        call abor1_ftn(err_msg)
      endif
      call ufo_geovals_get_var(geovals, var_t, t,status=ierr)
      if (ierr/=0) then
         write(err_msg,*) myname_, trim(var_t), ' doesnt exist'
         call abor1_ftn(err_msg)
      endif
      call ufo_geovals_get_var(geovals, var_mixr, mixr,status=ierr)
       if (ierr/=0) then
         write(err_msg,*) myname_, trim(var_mixr), ' doesnt exist'
         call abor1_ftn(err_msg)
      endif     
      call ufo_geovals_get_var(geovals, var_z, gph,status=ierr)
       if (ierr/=0) then
         write(err_msg,*) myname_, trim(var_z), ' doesnt exist'
         call abor1_ftn(err_msg)
      endif

      !Make sure nothing already allocated
      call self%delete()

      !Keep copy of dimensions
      self%nval = prs%nval
      self%nobs = obss%nobs
 
      allocate(self%wi(obss%nobs))
      allocate(self%wf(obss%nobs))
      allocate(self%t_traj(obss%nobs))
      allocate(self%mixr_traj(obss%nobs))
      allocate(self%prs_traj(obss%nobs))
      allocate(self%gph_traj(prs%nval, obss%nobs))
      allocate(self%obsH(obss%nobs))

      ! observation of altitude (MSL) (for vertical interpolation)
      call ioda_obsvec_setup(obsZ, obss%nobs)
      call ioda_obsdb_var_to_ovec(obss, obsZ, "HEIT")
      ! observation of Latitude (degree) (for geometric to geopotential height transform)
      call ioda_obsvec_setup(obsLat, obss%nobs)
      call ioda_obsdb_var_to_ovec(obss, obsLat, "Latitude")

      do iobs = 1, obss%nobs
        !  calculate observation geopotential height using  MJ Mahoney's (2001)
        call geometric2geop(obsLat%values(iobs), obsZ%values(iobs), self%obsH(iobs)  )
        call vert_interp_weights(self%nval, self%obsH(iobs), gph%vals(:,iobs),self%wi(iobs),self%wf(iobs))
        wi0 = self%wi(iobs)
        wf0 = self%wf(iobs)
        call vert_interp_apply(t%nval,   t%vals(:,iobs), self%t_traj(iobs),    self%wi(iobs),self%wf(iobs))   
        call vert_interp_apply(t%nval,mixr%vals(:,iobs), self%mixr_traj(iobs), self%wi(iobs),self%wf(iobs))
!        call vert_interp_apply(t%nval, prs%vals(:,iobs), self%prs_traj(iobs),  self%wi(iobs),self%wf(iobs))
       ! use  hypsometric equation to calculate pressure 
        Tv_traj  = 0.0
        Tv_traj0 = 0.0
        Tv_traj  = self%t_traj(iobs)*(one + (rv_over_rd-one)*self%mixr_traj(iobs) )
        Tv_traj0 = t%vals(wi0,iobs)*(one + (rv_over_rd-one)*mixr%vals(wi0,iobs) )
        self%prs_traj(iobs) = prs%vals(wi0,iobs)/exp(two*grav*(self%obsH(iobs)-gph%vals(wi0,iobs))/(rd*(Tv_traj+Tv_traj0)))

      enddo
  
       self%gph_traj=gph%vals
       self%ltraj = .true.
      ! cleanup
      call ioda_obsvec_delete(obsZ)
      call ioda_obsvec_delete(obsLat)
    end subroutine ufo_gnssro_ref_tlad_settraj
    
! ------------------------------------------------------------------------------
    
    subroutine ufo_gnssro_ref_tlad_tl(self, geovals, hofx, obss)
      use gnssro_mod_constants
      implicit none
      class(ufo_gnssro_tlad), intent(in)     :: self
      type(ufo_geovals),      intent(in)     :: geovals
      type(obs_vector),       intent(inout)  :: hofx
      type(ioda_obsdb),       intent(in)     :: obss
      logical,                parameter :: use_compress=.true.
     
      character(len=*), parameter :: myname_="ufo_gnssro_ref_tlad_tl"
      character(max_string)       :: err_msg
      
      integer                   :: iobs,ierr
      type(ufo_geoval), pointer :: t_d, mixr_d, prs_d
      real(kind_real)           :: t_coeff, prs_coeff, mixr_coeff
      real(kind_real)           :: gesT_d, gesQ_d, gesP_d,gesTv_d, gesTv0_d
      real(kind_real)           :: wf0
      integer                   :: wi0
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
     
      ! check if variables are in geovals and get them 
      call ufo_geovals_get_var(geovals, var_prs, prs_d, status=ierr )
      if (ierr/=0) then
        write(err_msg,*) myname_, trim(var_prs), ' doesnt exist'
        call abor1_ftn(err_msg)
      endif
      call ufo_geovals_get_var(geovals, var_t, t_d,status=ierr)
      if (ierr/=0) then
         write(err_msg,*) myname_, trim(var_t), ' doesnt exist'
         call abor1_ftn(err_msg)
      endif
      call ufo_geovals_get_var(geovals, var_mixr, mixr_d,status=ierr)
      if (ierr/=0) then
         write(err_msg,*) myname_, trim(var_mixr), ' doesnt exist'
         call abor1_ftn(err_msg)
      endif
  

      call gnssro_ref_constants(use_compress)

      ! tangent linear obs operator (linear)
      do iobs = 1, hofx%nobs
        wi0 = self%wi(iobs)
        wf0 = self%wf(iobs)
        call vert_interp_apply(t_d%nval,   t_d%vals(:,iobs), gesT_d, self%wi(iobs),self%wf(iobs))
        call vert_interp_apply(t_d%nval,mixr_d%vals(:,iobs), gesQ_d, self%wi(iobs),self%wf(iobs))
        call vert_interp_apply(t_d%nval,prs_d%vals(:,iobs), gesP_d, self%wi(iobs),self%wf(iobs))

        prs_coeff  =   n_a/self%t_traj(iobs)   &
                     + n_b*self%mixr_traj(iobs)/((self%mixr_traj(iobs)+rd_over_rv)*self%t_traj(iobs)**2)   &
                     + n_c*self%mixr_traj(iobs)/((self%mixr_traj(iobs)+rd_over_rv)*self%t_traj(iobs))
        t_coeff    = - n_a*self%prs_traj(iobs)/self%t_traj(iobs)**2           &
                     - n_b*two*self%prs_traj(iobs)*self%mixr_traj(iobs)/  &
                          ( (self%mixr_traj(iobs)+rd_over_rv)*self%t_traj(iobs)**3  )   &
                     - n_c*self%prs_traj(iobs)*self%mixr_traj(iobs)/((self%mixr_traj(iobs)+rd_over_rv)*self%t_traj(iobs)**2)
        mixr_coeff =   n_b*self%prs_traj(iobs)/( self%t_traj(iobs)**2*(self%mixr_traj(iobs)+rd_over_rv)**2 ) *  &
                       rd_over_rv                        &
                     + n_c*self%prs_traj(iobs)/( self%t_traj(iobs)**2*(self%mixr_traj(iobs)+rd_over_rv) ) *  &
                       rd_over_rv          
        hofx%values(iobs)   = prs_coeff*gesP_d  +                     &
                              t_coeff*gesT_d    +   mixr_coeff*gesQ_d

      enddo
    
    end subroutine ufo_gnssro_ref_tlad_tl
    
! ------------------------------------------------------------------------------
    
    subroutine ufo_gnssro_ref_tlad_ad(self, geovals, hofx, obss)
      use gnssro_mod_constants
      implicit none
      class(ufo_gnssro_tlad), intent(in)     :: self
      type(ufo_geovals),         intent(inout)  :: geovals
      type(obs_vector),          intent(in)     :: hofx
      type(ioda_obsdb),          intent(in)     :: obss
      logical,                   parameter :: use_compress=.true.

      character(len=*), parameter :: myname_="ufo_gnssro_ref_tlad_ad"
      character(max_string) :: err_msg
      
      integer :: iobs,ierr
      type(ufo_geoval), pointer :: t_d, mixr_d, prs_d
      real(kind_real) :: t_coeff, prs_coeff, mixr_coeff
      real(kind_real) :: gesT_d, gesQ_d, gesP_d,gesTv_d, gesTv0_d

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
      
      ! check if variables are in geovals and get them
      call ufo_geovals_get_var(geovals, var_prs, prs_d, status=ierr)
      if (ierr/=0) then
        write(err_msg,*) myname_, trim(var_prs), ' doesnt exist'
        call abor1_ftn(err_msg)
      endif
      
      call ufo_geovals_get_var(geovals, var_t, t_d, status=ierr)
      if (ierr/=0) then
        write(err_msg,*) myname_, trim(var_t), ' doesnt exist'
        call abor1_ftn(err_msg)
      endif 
      
      call ufo_geovals_get_var(geovals, var_mixr, mixr_d, status=ierr)
      if (ierr/=0) then
        write(err_msg,*) myname_, trim(var_mixr), ' doesnt exist'
        call abor1_ftn(err_msg)
      endif    
      ! allocate if not yet allocated	
      if (.not. allocated(t_d%vals)) then
         t_d%nobs = self%nobs
         t_d%nval = self%nval
         allocate(t_d%vals(t_d%nval,t_d%nobs))
      endif
      t_d%vals = 0.0_kind_real

      if (.not. allocated(prs_d%vals)) then
         prs_d%nobs = self%nobs
         prs_d%nval = self%nval
         allocate(prs_d%vals(prs_d%nval,prs_d%nobs))
      endif
      prs_d%vals = 0.0_kind_real

      if (.not. allocated(mixr_d%vals)) then
         mixr_d%nobs = self%nobs
         mixr_d%nval = self%nval
         allocate(mixr_d%vals(mixr_d%nval,mixr_d%nobs))
      endif
      mixr_d%vals = 0.0_kind_real

      if (.not. geovals%linit ) geovals%linit=.true.

      call gnssro_ref_constants(use_compress)

      do iobs = 1, hofx%nobs

        prs_coeff  =   n_a/self%t_traj(iobs)   &
                     + n_b*self%mixr_traj(iobs)/((self%mixr_traj(iobs)+rd_over_rv)*self%t_traj(iobs)**2)   &
                     + n_c*self%mixr_traj(iobs)/((self%mixr_traj(iobs)+rd_over_rv)*self%t_traj(iobs))
        t_coeff    = - n_a*self%prs_traj(iobs)/self%t_traj(iobs)**2           &
                     - n_b*two*self%prs_traj(iobs)*self%mixr_traj(iobs)/  &
                          ( (self%mixr_traj(iobs)+rd_over_rv)*self%t_traj(iobs)**3  )   &
                     - n_c*self%prs_traj(iobs)*self%mixr_traj(iobs)/((self%mixr_traj(iobs)+rd_over_rv)*self%t_traj(iobs)**2)
        mixr_coeff =   n_b*self%prs_traj(iobs)/( self%t_traj(iobs)**2*(self%mixr_traj(iobs)+rd_over_rv)**2 ) *  &
                       rd_over_rv                        &
                     + n_c*self%prs_traj(iobs)/( self%t_traj(iobs)**2*(self%mixr_traj(iobs)+rd_over_rv) ) *  &
                       rd_over_rv

        gesT_d = 0.0_kind_real
        gesQ_d = 0.0_kind_real
        gesP_d = 0.0_kind_real
        gesT_d = gesT_d + hofx%values(iobs)*t_coeff
        gesQ_d = gesQ_d + hofx%values(iobs)*mixr_coeff
        gesP_d = gesP_d + hofx%values(iobs)*prs_coeff
        call vert_interp_apply_ad(t_d%nval,    t_d%vals(:,iobs),    gesT_d, self%wi(iobs), self%wf(iobs))
        call vert_interp_apply_ad(mixr_d%nval, mixr_d%vals(:,iobs), gesQ_d, self%wi(iobs), self%wf(iobs))
        call vert_interp_apply_ad(prs_d%nval,  prs_d%vals(:,iobs),  gesP_d, self%wi(iobs), self%wf(iobs))

      enddo

    end subroutine ufo_gnssro_ref_tlad_ad
    
! ------------------------------------------------------------------------------
    
    subroutine ufo_gnssro_ref_tlad_delete(self)
      implicit none
      class(ufo_gnssro_tlad), intent(inout) :: self
      character(len=*), parameter :: myname_="ufo_gnssro_ref_tlad_delete"
      
      self%nval = 0
      if (allocated(self%wi)) deallocate(self%wi)
      if (allocated(self%wf)) deallocate(self%wf)
      if (allocated(self%prs_traj))  deallocate(self%prs_traj)
      if (allocated(self%t_traj))    deallocate(self%t_traj)
      if (allocated(self%mixr_traj)) deallocate(self%mixr_traj)
      if (allocated(self%gph_traj))    deallocate(self%gph_traj)
      if (allocated(self%obsH))    deallocate(self%obsH)
      self%ltraj = .false. 
    end subroutine ufo_gnssro_ref_tlad_delete
    
! ------------------------------------------------------------------------------

end module ufo_gnssro_ref_tlad_mod
