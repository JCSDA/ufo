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
  use vert_interp_mod
  use ufo_basis_tlad_mod, only: ufo_basis_tlad
  use obsspace_mod

  integer, parameter :: max_string=800

  type, extends(ufo_basis_tlad) :: ufo_conventional_profile_tlad
   private
     integer :: nval, nobs
     real(kind_real), allocatable :: wf(:)
     integer, allocatable :: wi(:)
  contains
    procedure :: delete => conventional_profile_tlad_delete_
    procedure :: settraj => conventional_profile_tlad_settraj_
    procedure :: simobs_tl => conventional_profile_simobs_tl_
    procedure :: simobs_ad => conventional_profile_simobs_ad_
  end type ufo_conventional_profile_tlad
contains

! ------------------------------------------------------------------------------

    subroutine conventional_profile_tlad_settraj_(self, geovals, obss)
      implicit none
      class(ufo_conventional_profile_tlad), intent(inout) :: self
      type(ufo_geovals),         intent(in)    :: geovals
      type(c_ptr), value,        intent(in)    :: obss

      character(len=*), parameter :: myname_="ufo_conventional_profile_tlad_settraj"
      character(max_string) :: err_msg

      real(kind_real), allocatable :: pressure(:)
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
      self%nobs = obsspace_get_nobs(obss)

      allocate(self%wi(self%nobs))
      allocate(self%wf(self%nobs))

      ! observation of pressure (for vertical interpolation)
      allocate(pressure(self%nobs))
      call obsspace_get_db(obss, "MetaData", "air_pressure", pressure)

      ! compute interpolation weights
      do iobs = 1, self%nobs
        call vert_interp_weights(self%nval,log(pressure(iobs)/10.),prsl%vals(:,iobs),self%wi(iobs),self%wf(iobs))
      enddo

      self%ltraj = .true.
      ! cleanup
      deallocate(pressure)
    end subroutine conventional_profile_tlad_settraj_

! ------------------------------------------------------------------------------

    subroutine conventional_profile_simobs_tl_(self, geovals, hofx, obss)
      implicit none
      class(ufo_conventional_profile_tlad), intent(in)     :: self
      type(ufo_geovals),         intent(in) :: geovals
      real(c_double),         intent(inout) :: hofx(:)
      type(c_ptr), value,        intent(in) :: obss
      
      character(len=*), parameter :: myname_="ufo_conventional_profile_simobs_tl"
      character(max_string) :: err_msg

      integer :: iobs,ierr
      type(ufo_geoval), pointer :: tv_d

      ! check if trajectory was set
      if (.not. self%ltraj) then
        write(err_msg,*) myname_, ' trajectory wasnt set!'
        call abor1_ftn(err_msg)
      endif

      ! check if nobs is consistent in geovals & hofx
!      if (geovals%nobs /= hofx%nobs) then
!        write(err_msg,*) myname_, ' error: nobs inconsistent!'
!        call abor1_ftn(err_msg)
!      endif

      ! check if tv variable is in geovals and get it
      call ufo_geovals_get_var(geovals, var_tv, tv_d, status=ierr )
      if (ierr/=0) then
        write(err_msg,*) myname_, trim(var_tv), ' doesnt exist'
        call abor1_ftn(err_msg)
      endif

      ! tangent linear obs operator (linear)
      do iobs = 1, geovals%nobs
        call vert_interp_apply_tl(tv_d%nval, tv_d%vals(:,iobs), hofx(iobs), self%wi(iobs), self%wf(iobs))
      enddo

    end subroutine conventional_profile_simobs_tl_

! ------------------------------------------------------------------------------

    subroutine conventional_profile_simobs_ad_(self, geovals, hofx, obss)
      implicit none
      class(ufo_conventional_profile_tlad), intent(in)     :: self
      type(ufo_geovals),         intent(inout)  :: geovals
      real(c_double),            intent(in)     :: hofx(:)
      type(c_ptr), value,        intent(in)     :: obss
      
      character(len=*), parameter :: myname_="ufo_conventional_profile_simobs_ad"
      character(max_string) :: err_msg

      integer :: iobs,ierr
      type(ufo_geoval), pointer :: tv_d

      ! check if trajectory was set
      if (.not. self%ltraj) then
        write(err_msg,*) myname_, ' trajectory wasnt set!'
        call abor1_ftn(err_msg)
      endif

      ! check if nobs is consistent in geovals & hofx
!      if (geovals%nobs /= hofx%nobs) then
!        write(err_msg,*) myname_, ' error: nobs inconsistent!'
!        call abor1_ftn(err_msg)
!      endif

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
      if (.not. geovals%linit ) geovals%linit=.true.

      do iobs = 1, geovals%nobs
        call vert_interp_apply_ad(tv_d%nval, tv_d%vals(:,iobs), hofx(iobs), self%wi(iobs), self%wf(iobs))
      enddo

    end subroutine conventional_profile_simobs_ad_

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
