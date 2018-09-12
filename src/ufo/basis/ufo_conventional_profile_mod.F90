! (C) Copyright 2017-2018 UCAR
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 

module ufo_conventional_profile_mod

  use iso_c_binding
  use kinds
  use ufo_vars_mod
  use ufo_geovals_mod
  use ufo_geovals_mod_c,   only: ufo_geovals_registry
  use ioda_obs_vectors
  use ioda_obsdb_mod
  use ioda_obsdb_mod_c, only: ioda_obsdb_registry
  use vert_interp_mod
  use ufo_basis_mod, only: ufo_basis

  integer, parameter :: max_string=800

  type, extends(ufo_basis) :: ufo_conventional_profile
   private
     integer :: nval, nobs
     real(kind_real), allocatable :: wf(:)
     integer, allocatable :: wi(:)
  contains
    procedure :: eqv    => conventional_profile_t_eqv_
  end type ufo_conventional_profile
contains

! ------------------------------------------------------------------------------

    subroutine conventional_profile_t_eqv_(self, geovals, hofx, obss)
    
      implicit none
      class(ufo_conventional_profile), intent(in)  :: self
      type(ufo_geovals), intent(in)                :: geovals
      type(obs_vector),  intent(inout)             :: hofx
      type(ioda_obsdb), target, intent(in)         :: obss
      
      character(len=*), parameter :: myname_="ufo_conventional_profile_t_eqv"
      character(max_string) :: err_msg
      
      integer :: iobs
      real(kind_real) :: wf
      integer :: wi,ierr
      type(obs_vector) :: pressure
      type(ufo_geoval), pointer :: prsl, tv
      
      ! check if nobs is consistent in geovals & hofx
      if (geovals%nobs /= hofx%nobs) then
        write(err_msg,*) myname_, ' error: nobs inconsistent!'
        call abor1_ftn(err_msg)
      endif
      
      ! check if prsl variable is in geovals and get it
      call ufo_geovals_get_var(geovals, var_prsl, prsl,status=ierr)
      if (ierr/=0) then
         write(err_msg,*) myname_, trim(var_prsl), ' doesnt exist'
         call abor1_ftn(err_msg)
      endif

      ! check if tv variable is in geovals and get it
      call ufo_geovals_get_var(geovals, var_tv, tv,status=ierr)
      if (ierr/=0) then
         write(err_msg,*) myname_, trim(var_tv), ' doesnt exist'
         call abor1_ftn(err_msg)
      endif
      
      ! observation of pressure (for vertical interpolation)
      call ioda_obsvec_setup(pressure, obss%nobs)
      call ioda_obsdb_var_to_ovec(obss, pressure, "air_pressure")
      
      ! obs operator
      do iobs = 1, hofx%nobs
        call vert_interp_weights(prsl%nval,log(pressure%values(iobs)/10.),prsl%vals(:,iobs),wi,wf)
        call vert_interp_apply(tv%nval, tv%vals(:,iobs), hofx%values(iobs), wi, wf)
      enddo

      ! cleanup
      call ioda_obsvec_delete(pressure)
    
    end subroutine conventional_profile_t_eqv_

! ------------------------------------------------------------------------------

end module ufo_conventional_profile_mod
