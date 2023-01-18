! (C) Copyright 2017-2018 UCAR
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 

!> Fortran module to handle gnssro refractivity observations

module ufo_gnssro_refncep_mod
  use fckit_configuration_module, only: fckit_configuration 
  use kinds
  use ufo_vars_mod
  use ufo_geovals_mod
  use ufo_geovals_mod_c, only: ufo_geovals_registry
  use vert_interp_mod
  use ufo_basis_mod,     only: ufo_basis
  use obsspace_mod
  use gnssro_mod_conf
  use gnssro_mod_constants
  use ufo_constants_mod

  implicit none
  public             :: ufo_gnssro_RefNCEP
  private

  !> Fortran derived type for gnssro trajectory
  type, extends(ufo_basis) :: ufo_gnssro_RefNCEP
  private
  type(gnssro_conf) :: roconf
  contains
   procedure :: setup     => ufo_gnssro_refncep_setup
   procedure :: simobs    => ufo_gnssro_refncep_simobs
  end type ufo_gnssro_RefNCEP

contains
! ------------------------------------------------------------------------------
   subroutine ufo_gnssro_refncep_setup(self, f_conf)
    implicit none
    class(ufo_gnssro_RefNCEP), intent(inout)  :: self
    type(fckit_configuration), intent(in)     :: f_conf

    call gnssro_conf_setup(self%roconf,f_conf)

   end subroutine ufo_gnssro_refncep_setup

   subroutine ufo_gnssro_refncep_simobs(self, geovals, hofx, obss)
    use gnssro_mod_transform, only: geometric2geop
      implicit none
      class(ufo_gnssro_RefNCEP), intent(in)     :: self
      type(ufo_geovals),         intent(in)     :: geovals
      real(kind_real),           intent(inout)  :: hofx(:)
      type(c_ptr), value,        intent(in)     :: obss
      
      character(len=*), parameter :: myname_="ufo_gnssro_refncep_simobs"
      character(max_string) :: err_msg
      integer           :: iobs,nlocs
      real(kind_real)   :: wf 
      integer           :: wi
      type(ufo_geoval), pointer    :: t,q,prs,gph
      real(kind_real)              :: refr1, refr2,refr3
      real(kind_real), allocatable :: obsZ(:), obsLat(:)
      real(kind_real)  :: obsH, gesT,gesQ, gesTv, gesTv0,gesP

      ! check if nlocs is consistent in geovals & hofx
      if (geovals%nlocs /= size(hofx)) then
        write(err_msg,*) myname_, ' error: nlocs inconsistent!'
        call abor1_ftn(err_msg)
      endif
      ! get variables from geovals
      call ufo_geovals_get_var(geovals, var_prs, prs)
      call ufo_geovals_get_var(geovals, var_ts,t)
      call ufo_geovals_get_var(geovals, var_q, q)
      call ufo_geovals_get_var(geovals, var_z, gph)

      nlocs = obsspace_get_nlocs(obss)

      allocate(obsZ(nlocs))
      allocate(obsLat(nlocs))

      call obsspace_get_db(obss, "MetaData", "height", obsZ)
      call obsspace_get_db(obss, "MetaData", "latitude", obsLat)

      call gnssro_ref_constants(self%roconf%use_compress)


      ! obs operator
      do iobs = 1, geovals%nlocs

      !  Convert geometric height at observation to geopotential height 
           call geometric2geop(obsLat(iobs), obsZ(iobs), obsH)
           call vert_interp_weights(gph%nval,obsH, gph%vals(:,iobs),wi,wf)  ! calculate weights 
           call vert_interp_apply(t%nval,   t%vals(:,iobs), gesT, wi, wf)
           call vert_interp_apply(q%nval,   q%vals(:,iobs), gesQ, wi, wf)

      !    use  hypsometric equation to calculate pressure 
           gesTv  = 0.0
           gesTv0 = 0.0
           gesTv  = gesT*(one + (rv_over_rd-one)* (gesQ/(1-gesQ) ) )
           gesTv0 = t%vals(wi,iobs)*(one + (rv_over_rd-one) * (q%vals(wi,iobs)/(1-q%vals(wi,iobs)) ))
           gesP   = prs%vals(wi,iobs)/exp(two*grav*(obsH-gph%vals(wi,iobs))/(rd*(gesTv+gesTv0)))
           refr1  = n_a*gesP/gesT
           refr2  = n_b*gesP*gesQ/ ( gesT**2 * (rd_over_rv+(1-rd_over_rv)*gesQ) )
           refr3  = n_c*gesP*gesQ/ ( gesT    * (rd_over_rv+(1-rd_over_rv)*gesQ) )
           hofx(iobs)  = refr1 + refr2 + refr3
      enddo

      ! cleanup 
      deallocate(obsZ)
      deallocate(obsLat)
    end subroutine ufo_gnssro_refncep_simobs

! ------------------------------------------------------------------------------
end module ufo_gnssro_refncep_mod
