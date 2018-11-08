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
  use ufo_conventional_profile_mod, only: find_position

  integer, parameter :: max_string=800

  type, extends(ufo_basis_tlad) :: ufo_conventional_profile_tlad
   private
     integer :: nval, nlocs
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
      self%nlocs = obsspace_get_nlocs(obss)

      allocate(self%wi(self%nlocs))
      allocate(self%wf(self%nlocs))

      ! observation of pressure (for vertical interpolation)
      allocate(pressure(self%nlocs))
      call obsspace_get_db(obss, "4DLocation", "air_pressure", pressure)

      ! compute interpolation weights
      do iobs = 1, self%nlocs
        call vert_interp_weights(self%nval,log(DBLE(pressure(iobs))/10.),prsl%vals(:,iobs),self%wi(iobs),self%wf(iobs))
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

      integer :: iobs, ierr, ivar, geo_ivar, jj, nlocs, nvars

      type ufo_geoval_ptr
         type(ufo_geoval), pointer :: ptr
      end type ufo_geoval_ptr
      type(ufo_geoval_ptr), dimension(:), allocatable :: vals

      character(len=MAXVARLEN), allocatable :: geovnames(:)
      character(len=MAXVARLEN), allocatable :: obsvnames(:)


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

      ! **********************************************************
      !                           STEP 1
      ! **********************************************************

      ! Retrieving the required variables names for this ObsOperator
      geovnames = ufo_vars_vnames(geovals%variables)
      nvars = size(geovnames)
    
      ! Checking if all required model variables are in geovals and get its pointer.
      allocate(vals(nvars))
      do ivar = 1, nvars
         call ufo_geovals_get_var(geovals, geovnames(ivar), vals(ivar)%ptr, status=ierr)
         if (ierr/=0) then
            write(err_msg,*) myname_, " : ", trim(geovnames(ivar)), ' doesnt exist'
            call abor1_ftn(err_msg)
         endif
      enddo
      
      ! **********************************************************
      !                           STEP 2
      ! **********************************************************

      ! Get the variable names we wwant to caculate the hofx
      obsvnames = obsspace_get_vnames(obss, MAXVARLEN)
      nvars = size(obsvnames)

      nlocs = obsspace_get_nlocs(obss)

      jj = 1
      do ivar = 1, nvars
       ! Determine the location of this variable in geovals
        if (trim(obsvnames(ivar)) == "air_temperature") then ! not match, to be solved
          geo_ivar = find_position("virtual_temperature", size(geovnames), geovnames)
        else
          geo_ivar = find_position(obsvnames(ivar), size(geovnames), geovnames)
        endif
        if (geo_ivar == -999 ) then
          write(err_msg,*) myname_, " : ", trim(obsvnames(ivar)), ' is not in geovals'
          call abor1_ftn(err_msg)
        endif

        ! tangent linear obs operator (linear)
        do iobs = 1, nlocs
          call vert_interp_apply_tl(vals(geo_ivar)%ptr%nval, vals(geo_ivar)%ptr%vals(:,iobs), hofx(jj), self%wi(iobs), self%wf(iobs))
          jj = jj + 1
        enddo
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

      integer :: iobs, ierr, ivar, geo_ivar, jj, nlocs, nvars
      type(ufo_geoval), pointer :: tv_d

      type ufo_geoval_ptr
         type(ufo_geoval), pointer :: ptr
      end type ufo_geoval_ptr
      type(ufo_geoval_ptr), dimension(:), allocatable :: vals

      character(len=MAXVARLEN), allocatable :: geovnames(:)
      character(len=MAXVARLEN), allocatable :: obsvnames(:)

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

      ! **********************************************************
      !                           STEP 1
      ! **********************************************************

      ! Retrieving the required variables names for this ObsOperator
      geovnames = ufo_vars_vnames(geovals%variables)
      nvars = size(geovnames)
    
      ! Checking if all required model variables are in geovals and get its pointer.
      allocate(vals(nvars))
      do ivar = 1, nvars
         call ufo_geovals_get_var(geovals, geovnames(ivar), vals(ivar)%ptr, status=ierr)
         if (ierr/=0) then
            write(err_msg,*) myname_, " : ", trim(geovnames(ivar)), ' doesnt exist'
            call abor1_ftn(err_msg)
         endif
      enddo

      ! check if tv variable is in geovals and get it
      ! call ufo_geovals_get_var(geovals, var_tv, tv_d, status=ierr)
      ! if (ierr/=0) then
      !   write(err_msg,*) myname_, trim(var_tv), ' doesnt exist'
      !   call abor1_ftn(err_msg)
      ! endif

      ! **********************************************************
      !                           STEP 2
      ! **********************************************************

      ! Get the variable names we wwant to caculate the hofx
      obsvnames = obsspace_get_vnames(obss, MAXVARLEN)
      nvars = size(obsvnames)

      nlocs = obsspace_get_nlocs(obss)

      jj = 1
      do ivar = 1, nvars
       ! Determine the location of this variable in geovals
        if (trim(obsvnames(ivar)) == "air_temperature") then ! not match, to be solved
          geo_ivar = find_position("virtual_temperature", size(geovnames), geovnames)
        else
          geo_ivar = find_position(obsvnames(ivar), size(geovnames), geovnames)
        endif
        if (geo_ivar == -999 ) then
          write(err_msg,*) myname_, " : ", trim(obsvnames(ivar)), ' is not in geovals'
          call abor1_ftn(err_msg)
        endif
      
        ! allocate if not yet allocated
        if (.not. allocated(vals(geo_ivar)%ptr%vals)) then
           vals(geo_ivar)%ptr%nobs = nlocs
           vals(geo_ivar)%ptr%nval = nvars
           allocate(vals(geo_ivar)%ptr%vals(tv_d%nval, vals(geo_ivar)%ptr%nobs))
           vals(geo_ivar)%ptr%vals = 0.0_kind_real
        endif
        if (.not. geovals%linit ) geovals%linit=.true.

        do iobs = 1, nlocs
          call vert_interp_apply_ad(vals(geo_ivar)%ptr%nval, vals(geo_ivar)%ptr%vals(:,iobs), hofx(jj), self%wi(iobs), self%wf(iobs))
          jj = jj + 1
        enddo
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