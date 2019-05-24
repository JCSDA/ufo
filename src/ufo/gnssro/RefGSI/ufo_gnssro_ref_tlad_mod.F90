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
  use vert_interp_mod
  use ufo_basis_tlad_mod,  only: ufo_basis_tlad
  use obsspace_mod
  use config_mod
  use gnssro_mod_constants
  use gnssro_mod_conf
  use missing_values_mod

  !> Fortran derived type for gnssro trajectory
  type, extends(ufo_basis_tlad) :: ufo_gnssro_Ref_tlad
   private
  type(gnssro_conf) :: roconf
     integer                       :: nval, nlocs
     real(kind_real), allocatable  :: wf(:)
     integer,         allocatable  :: wi(:)
     real(kind_real), allocatable  :: prs(:), t(:), q(:)
     real(kind_real), allocatable  :: obsH(:)
  contains
    procedure :: setup      => ufo_gnssro_ref_tlad_setup
    procedure :: delete     => ufo_gnssro_ref_tlad_delete
    procedure :: settraj    => ufo_gnssro_ref_tlad_settraj
    procedure :: simobs_tl  => ufo_gnssro_ref_simobs_tl
    procedure :: simobs_ad  => ufo_gnssro_ref_simobs_ad
  end type ufo_gnssro_Ref_tlad

contains
! ------------------------------------------------------------------------------
   subroutine ufo_gnssro_ref_tlad_setup(self, c_conf)
    implicit none
    class(ufo_gnssro_Ref_tlad), intent(inout) :: self
    type(c_ptr),         intent(in)    :: c_conf

    call gnssro_conf_setup(self%roconf,c_conf)
   end subroutine ufo_gnssro_ref_tlad_setup
   
    subroutine ufo_gnssro_ref_tlad_settraj(self, geovals, obss)
      use gnssro_mod_transform, only: geometric2geop
      
      implicit none
      class(ufo_gnssro_Ref_tlad), intent(inout) :: self
      type(ufo_geovals),         intent(in)     :: geovals
      type(c_ptr), value,        intent(in)     :: obss
      
      character(len=*), parameter :: myname_="ufo_gnssro_ref_tlad_settraj"
      character(max_string)       :: err_msg
      
      type(ufo_geoval), pointer :: t,q,prs,gph
      integer          :: iobs
      real(kind_real), allocatable :: obsZ(:), obsLat(:)  ! observation vector
      real(kind_real)  :: Tv, Tv0
      integer          :: wi0

      !Get variables from geovals
      call ufo_geovals_get_var(geovals, var_prs, prs)
      call ufo_geovals_get_var(geovals, var_ts,t)
      call ufo_geovals_get_var(geovals, var_q, q)
      call ufo_geovals_get_var(geovals, var_z, gph)

      !Make sure nothing already allocated
      call self%delete()

      !Keep copy of dimensions
      self%nval = prs%nval
      self%nlocs = obsspace_get_nlocs(obss)
 
      allocate(self%wi(self%nlocs))
      allocate(self%wf(self%nlocs))
      allocate(self%t(self%nlocs))
      allocate(self%q(self%nlocs))
      allocate(self%prs(self%nlocs))
      allocate(self%obsH(self%nlocs))

      allocate(obsZ(self%nlocs))
      allocate(obsLat(self%nlocs))

      ! get observation vectors
      call obsspace_get_db(obss, "MetaData", "altitude", obsZ)
      call obsspace_get_db(obss, "MetaData", "latitude", obsLat)

      do iobs = 1, self%nlocs

        !  calculate observation geopotential height using  MJ Mahoney's (2001)
        call geometric2geop(obsLat(iobs), obsZ(iobs), self%obsH(iobs))
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
      deallocate(obsZ)
      deallocate(obsLat)
    end subroutine ufo_gnssro_ref_tlad_settraj
    
! ------------------------------------------------------------------------------
    
    subroutine ufo_gnssro_ref_simobs_tl(self, geovals, hofx, obss)
      implicit none
      class(ufo_gnssro_Ref_tlad), intent(in) :: self
      type(ufo_geovals),      intent(in)     :: geovals
      real(kind_real),        intent(inout)  :: hofx(:)
      type(c_ptr), value,     intent(in)     :: obss
     
      character(len=*), parameter :: myname_="ufo_gnssro_ref_tlad_tl"
      character(max_string)       :: err_msg
      
      integer                   :: iobs,ierr
      type(ufo_geoval), pointer :: t_d, q_d, prs_d 
      real(kind_real)           :: t_coeff, q_coeff, p_coeff
      real(kind_real)           :: gesT_d, gesQ_d, gesP_d
      integer                   :: wi0
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
      call ufo_geovals_get_var(geovals, var_ts,t_d)
      call ufo_geovals_get_var(geovals, var_q, q_d)
      call ufo_geovals_get_var(geovals, var_prs, prs_d)
 
      call gnssro_ref_constants(self%roconf%use_compress)

      ! tangent linear obs operator (linear)
      do iobs = 1, geovals%nlocs
        wi0 = self%wi(iobs)
        call vert_interp_apply_tl(  t_d%nval,  t_d%vals(:,iobs), gesT_d, self%wi(iobs),self%wf(iobs)) 
        call vert_interp_apply_tl(  q_d%nval,  q_d%vals(:,iobs), gesQ_d, self%wi(iobs),self%wf(iobs))
        call vert_interp_apply_tl(prs_d%nval,prs_d%vals(:,iobs), gesP_d, self%wi(iobs),self%wf(iobs))

!       pressure does not change during minimization

        t_coeff = - n_a*self%prs(iobs)/self%t(iobs)**2           &
                  - n_b*two*self%prs(iobs)*self%q(iobs)/  &
                       ( ((1-rd_over_rv)*self%q(iobs)+rd_over_rv)*self%t(iobs)**3  )   &
                  - n_c*self%prs(iobs)*self%q(iobs)/      &
                           ( ((1-rd_over_rv)*self%q(iobs)+rd_over_rv)*self%t(iobs)**2 )
        q_coeff =   n_b*self%prs(iobs)/( self%t(iobs)**2*( (1-rd_over_rv)*self%q(iobs)+rd_over_rv)**2 ) *  &
                       rd_over_rv                        &
                  + n_c*self%prs(iobs)/( self%t(iobs)   *( (1-rd_over_rv)*self%q(iobs)+rd_over_rv)**2 ) *  &
                       rd_over_rv          
        p_coeff =   n_a/self%t(iobs)   &
                  + n_b*self%q(iobs)/ ( ((1-rd_over_rv)*self%q(iobs)+rd_over_rv)*self%t(iobs)**2 )   &
                  + n_c*self%q(iobs)/ ( ((1-rd_over_rv)*self%q(iobs)+rd_over_rv)*self%t(iobs) )

     
        hofx(iobs)  =  t_coeff*gesT_d  + q_coeff*gesQ_d + p_coeff*gesP_d

      enddo 
    
    end subroutine ufo_gnssro_ref_simobs_tl
    
! ------------------------------------------------------------------------------
    
    subroutine ufo_gnssro_ref_simobs_ad(self, geovals, hofx, obss)
      implicit none
      class(ufo_gnssro_Ref_tlad), intent(in)   :: self
      type(ufo_geovals),         intent(inout) :: geovals
      real(kind_real),           intent(in)    :: hofx(:)
      type(c_ptr), value,        intent(in)    :: obss

      character(len=*), parameter :: myname_="ufo_gnssro_ref_tlad_ad"
      character(max_string)       :: err_msg
      real(c_double)              :: missing
      integer :: iobs,ierr
      type(ufo_geoval), pointer :: t_d, q_d, prs_d
      real(kind_real)           :: t_coeff, q_coeff, p_coeff
      real(kind_real)           :: gesT_d, gesQ_d, gesP_d


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
      call ufo_geovals_get_var(geovals, var_prs, prs_d)
      call ufo_geovals_get_var(geovals, var_ts,t_d)
      call ufo_geovals_get_var(geovals, var_q, q_d)

      ! allocate if not yet allocated	
      if (.not. allocated(t_d%vals)) then
         t_d%nlocs = self%nlocs
         t_d%nval = self%nval
         allocate(t_d%vals(t_d%nval,t_d%nlocs))
         t_d%vals = 0.0_kind_real
      endif

      if (.not. allocated(prs_d%vals)) then
         prs_d%nlocs = self%nlocs
         prs_d%nval = self%nval
         allocate(prs_d%vals(prs_d%nval,prs_d%nlocs))
         prs_d%vals = 0.0_kind_real
      endif

      if (.not. allocated(q_d%vals)) then
         q_d%nlocs = self%nlocs
         q_d%nval = self%nval
         allocate(q_d%vals(q_d%nval,q_d%nlocs))
         q_d%vals = 0.0_kind_real
      endif

      if (.not. geovals%linit ) geovals%linit=.true.

      call gnssro_ref_constants(self%roconf%use_compress)
      missing = missing_value(missing)

      do iobs = 1, geovals%nlocs
    
         if (hofx(iobs) .ne. missing) then
           t_coeff = - n_a*self%prs(iobs)/self%t(iobs)**2           &
                     - n_b*two*self%prs(iobs)*self%q(iobs)/  &
                           ( ((1-rd_over_rv)*self%q(iobs)+rd_over_rv)*self%t(iobs)**3  )   &
                     - n_c*self%prs(iobs)*self%q(iobs)/      &
                           ( ((1-rd_over_rv)*self%q(iobs)+rd_over_rv)*self%t(iobs)**2 )
           q_coeff =   n_b*self%prs(iobs)/( self%t(iobs)**2*( (1-rd_over_rv)*self%q(iobs)+rd_over_rv)**2 ) *  &
                           rd_over_rv                        &
                     + n_c*self%prs(iobs)/( self%t(iobs)   *( (1-rd_over_rv)*self%q(iobs)+rd_over_rv)**2 ) *  &
                           rd_over_rv
           p_coeff =   n_a/self%t(iobs)   &
                     + n_b*self%q(iobs)/ ( ((1-rd_over_rv)*self%q(iobs)+rd_over_rv)*self%t(iobs)**2 )   &
                     + n_c*self%q(iobs)/ ( ((1-rd_over_rv)*self%q(iobs)+rd_over_rv)*self%t(iobs) )

           gesT_d = 0.0_kind_real
           gesQ_d = 0.0_kind_real
           gesP_d = 0.0_kind_real
           gesT_d = gesT_d + hofx(iobs)*t_coeff
           gesQ_d = gesQ_d + hofx(iobs)*q_coeff
           gesP_d = gesP_d + hofx(iobs)*p_coeff
           call vert_interp_apply_ad(  t_d%nval,  t_d%vals(:,iobs), gesT_d, self%wi(iobs), self%wf(iobs))
           call vert_interp_apply_ad(  q_d%nval,  q_d%vals(:,iobs), gesQ_d, self%wi(iobs), self%wf(iobs))
           call vert_interp_apply_ad(prs_d%nval,prs_d%vals(:,iobs), gesP_d, self%wi(iobs), self%wf(iobs))
        endif
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
