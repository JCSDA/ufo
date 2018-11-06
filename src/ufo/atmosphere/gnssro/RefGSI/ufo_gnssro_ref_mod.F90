! (C) Copyright 2017-2018 UCAR
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 

!> Fortran module to handle gnssro refractivity observations

module ufo_gnssro_ref_mod
  use iso_c_binding
  use kinds
  use ufo_vars_mod
  use ufo_geovals_mod
  use ufo_geovals_mod_c, only: ufo_geovals_registry
  use vert_interp_mod
  use ufo_basis_mod,     only: ufo_basis
  use obsspace_mod
  
  implicit none
  integer, parameter :: max_string=800
  public             :: ufo_gnssro_Ref
  private

  !> Fortran derived type for gnssro trajectory
  type, extends(ufo_basis) :: ufo_gnssro_Ref
  contains
   procedure :: simobs    => ufo_gnssro_ref_simobs
  end type ufo_gnssro_Ref

contains
! ------------------------------------------------------------------------------
   subroutine ufo_gnssro_ref_simobs(self, geovals, hofx, obss)
    use gnssro_mod_constants
    use gnssro_mod_transform, only: geometric2geop
      implicit none
      logical,parameter                          :: use_compress=.true.
      class(ufo_gnssro_Ref), intent(in)          :: self
      type(ufo_geovals), intent(in)              :: geovals
      real(kind_real),   intent(inout)           :: hofx(:)
      type(c_ptr), value,       intent(in)       :: obss

      character(len=*), parameter :: myname_="ufo_gnssro_ref_simobs"
      character(max_string) :: err_msg

      integer         :: iobs,k
      real(kind_real) :: wf 
      integer         :: wi,ierr
      integer         :: nobs
      type(ufo_geoval), pointer  :: t,q,prs,gph
      real(kind_real)            :: refr1, refr2,refr3
      real(kind_real), allocatable :: obsZ(:), obsLat(:)
      real(kind_real)  :: obsH, gesT,gesQ, gesTv, gesTv0,gesP
      ! check if nobs is consistent in geovals & hofx
      if (geovals%nobs /= size(hofx)) then
        write(err_msg,*) myname_, ' error: nobs inconsistent!'
        call abor1_ftn(err_msg)
      endif
      ! check if prs variable is in geovals and get it
      call ufo_geovals_get_var(geovals, var_prs, prs,status=ierr)
      if (ierr/=0) then
         write(err_msg,*) myname_, trim(var_prs), ' doesnt exist'
         call abor1_ftn(err_msg)
      endif
       ! check if t variable is in geovals and get it
      call ufo_geovals_get_var(geovals, var_t, t,status=ierr)
      if (ierr/=0) then
         write(err_msg,*) myname_, trim(var_t), ' doesnt exist'
         call abor1_ftn(err_msg)
      endif
      call ufo_geovals_get_var(geovals, var_q, q,status=ierr)
       if (ierr/=0) then
         write(err_msg,*) myname_, trim(var_q), ' doesnt exist'
         call abor1_ftn(err_msg)
      endif
      call ufo_geovals_get_var(geovals, var_z, gph,status=ierr)
      if (ierr/=0) then
         write(err_msg,*) myname_, trim(var_z), ' doesnt exist'
         call abor1_ftn(err_msg)
       endif

      nobs = obsspace_get_nobs(obss)
      allocate(obsZ(nobs))
      allocate(obsLat(nobs))

      call obsspace_get_db(obss, "Metadata", "altitude", obsZ)
      call obsspace_get_db(obss, "Metadata", "latitude", obsLat)

      call gnssro_ref_constants(use_compress)

      ! obs operator
      do iobs = 1, geovals%nobs
      ! Convert geometric height at observation to geopotential height 
        call geometric2geop(obsLat(iobs), obsZ(iobs), obsH)
        call vert_interp_weights(gph%nval,obsH, gph%vals(:,iobs),wi,wf)  ! calculate weights 
        call vert_interp_apply(t%nval,   t%vals(:,iobs), gesT, wi, wf)
        call vert_interp_apply(q%nval,   q%vals(:,iobs), gesQ, wi, wf)

      ! use  hypsometric equation to calculate pressure 
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
    end subroutine ufo_gnssro_ref_simobs

! ------------------------------------------------------------------------------
end module ufo_gnssro_ref_mod
