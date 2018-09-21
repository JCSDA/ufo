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
  type, extends(ufo_basis_tlad) :: ufo_gnssro_Ref_tlad
   private
     integer                       :: nval, nobs
     real(kind_real), allocatable  :: wf(:)
     integer,         allocatable  :: wi(:)
     real(kind_real), allocatable  :: prs(:), t(:), q(:)
     real(kind_real), allocatable  :: obsH(:)
  contains
    procedure :: delete  => ufo_gnssro_ref_tlad_delete
    procedure :: settraj => ufo_gnssro_ref_tlad_settraj
    procedure :: simobs_tl  => ufo_gnssro_ref_simobs_tl
    procedure :: simobs_ad  => ufo_gnssro_ref_simobs_ad
  end type ufo_gnssro_Ref_tlad

contains
! ------------------------------------------------------------------------------
    
    subroutine ufo_gnssro_ref_tlad_settraj(self, geovals, obss)
      use gnssro_mod_constants
      use gnssro_mod_transform, only: geometric2geop
      
      implicit none
      class(ufo_gnssro_Ref_tlad), intent(inout) :: self
      type(ufo_geovals),         intent(in)     :: geovals
      type(ioda_obsdb),          intent(in)     :: obss
      
      character(len=*), parameter :: myname_="ufo_gnssro_ref_tlad_settraj"
      character(max_string)       :: err_msg
      
      type(ufo_geoval), pointer :: t,q,prs,gph
      integer          :: iobs, ierr
      type(obs_vector) :: obsZ, obsLat  ! observation vector
      real(kind_real)  :: Tv, Tv0
      integer          :: wi0

      !Check if variables are in geovals and get them
      call ufo_geovals_get_var(geovals, var_prs, prs, status=ierr)
      if (ierr/=0) then
        write(err_msg,*) myname_, trim(var_prs), ' doesnt exist'
        call abor1_ftn(err_msg)
      endif
      call ufo_geovals_get_var(geovals, var_t, t, status=ierr)
      if (ierr/=0) then
         write(err_msg,*) myname_, trim(var_t), ' doesnt exist'
         call abor1_ftn(err_msg)
      endif
      call ufo_geovals_get_var(geovals, var_q, q, status=ierr)
       if (ierr/=0) then
         write(err_msg,*) myname_, trim(var_q), ' doesnt exist'
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
      allocate(self%t(obss%nobs))
      allocate(self%q(obss%nobs))
      allocate(self%prs(obss%nobs))
      allocate(self%obsH(obss%nobs))

      ! observation of altitude (MSL) (for vertical interpolation)
      call ioda_obsvec_setup(obsZ, obss%nobs)
      call ioda_obsdb_var_to_ovec(obss, obsZ, "MSL_ALT")
      ! observation of Latitude (degree) (for geometric to geopotential height transform)
      call ioda_obsvec_setup(obsLat, obss%nobs)
      call ioda_obsdb_var_to_ovec(obss, obsLat, "Latitude")

      do iobs = 1, obss%nobs

        !  calculate observation geopotential height using  MJ Mahoney's (2001)
        call geometric2geop(obsLat%values(iobs), obsZ%values(iobs), self%obsH(iobs))
        call vert_interp_weights(self%nval, self%obsH(iobs), gph%vals(:,iobs),self%wi(iobs),self%wf(iobs))
        wi0 = self%wi(iobs)
        call vert_interp_apply(t%nval, t%vals(:,iobs), self%t(iobs), self%wi(iobs),self%wf(iobs))
        call vert_interp_apply(q%nval, q%vals(:,iobs), self%q(iobs), self%wi(iobs),self%wf(iobs))
       ! use  hypsometric equation to calculate pressure 
        Tv  = 0.0
        Tv0 = 0.0
        Tv  = self%t(iobs)*(one + (rv_over_rd-one)*self%q(iobs)/(1.0-self%q(iobs)) )
        Tv0 = t%vals(wi0,iobs)*(one + (rv_over_rd-one)*q%vals(wi0,iobs)/(1.0-q%vals(wi0,iobs) ))
        self%prs(iobs) = prs%vals(wi0,iobs)/exp(two*grav*(self%obsH(iobs)-gph%vals(wi0,iobs))/(rd*(Tv+Tv0)))

      enddo
  
       self%ltraj = .true.
      ! cleanup
      call ioda_obsvec_delete(obsZ)
      call ioda_obsvec_delete(obsLat)
    end subroutine ufo_gnssro_ref_tlad_settraj
    
! ------------------------------------------------------------------------------
    
    subroutine ufo_gnssro_ref_simobs_tl(self, geovals, hofx, obss)
      use gnssro_mod_constants
      implicit none
      class(ufo_gnssro_Ref_tlad), intent(in) :: self
      type(ufo_geovals),      intent(in)     :: geovals
      type(obs_vector),       intent(inout)  :: hofx
      type(ioda_obsdb),       intent(in)     :: obss
      logical,                parameter      :: use_compress=.true.
     
      character(len=*), parameter :: myname_="ufo_gnssro_ref_tlad_tl"
      character(max_string)       :: err_msg
      
      integer                   :: iobs,ierr
      type(ufo_geoval), pointer :: t_d, q_d                          !, prs_d 
      real(kind_real)           :: t_coeff, q_coeff                  !, prs_coeff
      real(kind_real)           :: gesT_d, gesQ_d, gesTv_d, gesTv0_d !, gesP_d
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
      call ufo_geovals_get_var(geovals, var_t, t_d,status=ierr)
      if (ierr/=0) then
         write(err_msg,*) myname_, trim(var_t), ' doesnt exist'
         call abor1_ftn(err_msg)
      endif
      call ufo_geovals_get_var(geovals, var_q, q_d,status=ierr)
      if (ierr/=0) then
         write(err_msg,*) myname_, trim(var_q), ' doesnt exist'
         call abor1_ftn(err_msg)
      endif
  

      call gnssro_ref_constants(use_compress)

      ! tangent linear obs operator (linear)
      do iobs = 1, hofx%nobs
        wi0 = self%wi(iobs)
        call vert_interp_apply_tl(  t_d%nval,  t_d%vals(:,iobs), gesT_d, self%wi(iobs),self%wf(iobs)) 
        call vert_interp_apply_tl(  q_d%nval,  q_d%vals(:,iobs), gesQ_d, self%wi(iobs),self%wf(iobs))
!       call vert_interp_apply_tl(prs_d%nval,prs_d%vals(:,iobs), gesP_d, self%wi(iobs),self%wf(iobs))
!       pressure does not change during minimization

!       prs_coeff  =   n_a/self%t(iobs)   &
!                    + n_b*self%q(iobs)/ ( ((1-rd_over_rv)*self%q(iobs)+rd_over_rv)*self%t(iobs)**2 )   &
!                    + n_c*self%q(iobs)/ ( ((1-rd_over_rv)*self%q(iobs)+rd_over_rv)*self%t(iobs) )
        t_coeff    = - n_a*self%prs(iobs)/self%t(iobs)**2           &
                     - n_b*two*self%prs(iobs)*self%q(iobs)/  &
                           ( ((1-rd_over_rv)*self%q(iobs)+rd_over_rv)*self%t(iobs)**3  )   &
                     - n_c*self%prs(iobs)*self%q(iobs)/      &
                           ( ((1-rd_over_rv)*self%q(iobs)+rd_over_rv)*self%t(iobs)**2 )
        q_coeff =   n_b*self%prs(iobs)/( self%t(iobs)**2*( (1-rd_over_rv)*self%q(iobs)+rd_over_rv)**2 ) *  &
                       rd_over_rv                        &
                  + n_c*self%prs(iobs)/( self%t(iobs)   *( (1-rd_over_rv)*self%q(iobs)+rd_over_rv)**2 ) *  &
                       rd_over_rv          
        hofx%values(iobs)   =  t_coeff*gesT_d  + q_coeff*gesQ_d !+ prs_coeff*gesP_d

      enddo 
    
    end subroutine ufo_gnssro_ref_simobs_tl
    
! ------------------------------------------------------------------------------
    
    subroutine ufo_gnssro_ref_simobs_ad(self, geovals, hofx, obss)
      use gnssro_mod_constants
      implicit none
      class(ufo_gnssro_Ref_tlad), intent(in)   :: self
      type(ufo_geovals),         intent(inout) :: geovals
      type(obs_vector),          intent(in)    :: hofx
      type(ioda_obsdb),          intent(in)    :: obss
      logical,                   parameter     :: use_compress=.true.

      character(len=*), parameter :: myname_="ufo_gnssro_ref_tlad_ad"
      character(max_string) :: err_msg
      
      integer :: iobs,ierr
      type(ufo_geoval), pointer :: t_d, q_d, prs_d, gph_d
      real(kind_real)           :: t_coeff, q_coeff                    !, prs_coeff
      real(kind_real)           :: gesT_d, gesQ_d !, gesTv_d,gesTv0_d    !, gesP_d

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
      
      call ufo_geovals_get_var(geovals, var_q, q_d, status=ierr)
      if (ierr/=0) then
        write(err_msg,*) myname_, trim(var_q), ' doesnt exist'
        call abor1_ftn(err_msg)
      endif    
      call ufo_geovals_get_var(geovals, var_z, gph_d, status=ierr)
      if (ierr/=0) then
        write(err_msg,*) myname_, trim(var_z), ' doesnt exist'
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

      if (.not. allocated(q_d%vals)) then
         q_d%nobs = self%nobs
         q_d%nval = self%nval
         allocate(q_d%vals(q_d%nval,q_d%nobs))
      endif
      q_d%vals = 0.0_kind_real

      if (.not. allocated(gph_d%vals)) then
         gph_d%nobs = self%nobs
         gph_d%nval = self%nval 
         allocate(gph_d%vals(gph_d%nval,gph_d%nobs))
      endif
      gph_d%vals = 0.0_kind_real

      if (.not. geovals%linit ) geovals%linit=.true.

      call gnssro_ref_constants(use_compress)


      do iobs = 1, hofx%nobs

! zero impct on pressure during minimization
!       prs_coeff  =   n_a/self%t(iobs)   &
!                    + n_b*self%q(iobs)/ ( ((1-rd_over_rv)*self%q(iobs)+rd_over_rv)*self%t(iobs)**2 )   &
!                    + n_c*self%q(iobs)/ ( ((1-rd_over_rv)*self%q(iobs)+rd_over_rv)*self%t(iobs) )
        t_coeff    = - n_a*self%prs(iobs)/self%t(iobs)**2           &
                     - n_b*two*self%prs(iobs)*self%q(iobs)/  &
                           ( ((1-rd_over_rv)*self%q(iobs)+rd_over_rv)*self%t(iobs)**3  )   &
                     - n_c*self%prs(iobs)*self%q(iobs)/      &
                           ( ((1-rd_over_rv)*self%q(iobs)+rd_over_rv)*self%t(iobs)**2 )
        q_coeff =   n_b*self%prs(iobs)/( self%t(iobs)**2*( (1-rd_over_rv)*self%q(iobs)+rd_over_rv)**2 ) *  &
                       rd_over_rv                        &
                  + n_c*self%prs(iobs)/( self%t(iobs)   *( (1-rd_over_rv)*self%q(iobs)+rd_over_rv)**2 ) *  &
                       rd_over_rv
  
        gesT_d = 0.0_kind_real
        gesQ_d = 0.0_kind_real
!       gesP_d = 0.0_kind_real
        gesT_d = gesT_d + hofx%values(iobs)*t_coeff
        gesQ_d = gesQ_d + hofx%values(iobs)*q_coeff
!       gesP_d = gesP_d + hofx%values(iobs)*prs_coeff
        call vert_interp_apply_ad(  t_d%nval,  t_d%vals(:,iobs), gesT_d, self%wi(iobs), self%wf(iobs))
        call vert_interp_apply_ad(  q_d%nval,  q_d%vals(:,iobs), gesQ_d, self%wi(iobs), self%wf(iobs))
!       call vert_interp_apply_ad(prs_d%nval,prs_d%vals(:,iobs), gesP_d, self%wi(iobs), self%wf(iobs))

      enddo

    end subroutine ufo_gnssro_ref_simobs_ad
    
! ------------------------------------------------------------------------------
    
    subroutine ufo_gnssro_ref_tlad_delete(self)
      implicit none
      class(ufo_gnssro_Ref_tlad), intent(inout) :: self
      character(len=*), parameter :: myname_="ufo_gnssro_ref_tlad_delete"
      
      self%nval = 0
      if (allocated(self%wi))  deallocate(self%wi)
      if (allocated(self%wf))  deallocate(self%wf)
      if (allocated(self%prs)) deallocate(self%prs)
      if (allocated(self%t))   deallocate(self%t)
      if (allocated(self%q))   deallocate(self%q)
      if (allocated(self%obsH))deallocate(self%obsH)
      self%ltraj = .false. 
    end subroutine ufo_gnssro_ref_tlad_delete
    
! ------------------------------------------------------------------------------

end module ufo_gnssro_ref_tlad_mod
