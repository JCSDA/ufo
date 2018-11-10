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
  use vert_interp_mod
  use ufo_basis_mod, only: ufo_basis
  use obsspace_mod

  public find_position

  integer, parameter :: max_string=800

  type, extends(ufo_basis) :: ufo_conventional_profile
   private
     integer :: nval, nobs
     real(kind_real), allocatable :: wf(:)
     integer, allocatable :: wi(:)
     integer, public :: nvars
     character(len=max_string), public, allocatable :: varout(:)
  contains
    procedure :: simobs    => conventional_profile_simobs_
    final :: destructor
  end type ufo_conventional_profile
contains

! ------------------------------------------------------------------------------

    subroutine conventional_profile_simobs_(self, geovals, hofx, obss)

      implicit none
      class(ufo_conventional_profile), intent(in)  :: self
      type(ufo_geovals), intent(in)                :: geovals
      real(c_double),  intent(inout)               :: hofx(:)
      type(c_ptr), value, intent(in)               :: obss

      character(len=*), parameter :: myname_="ufo_conventional_profile_simobs"
      character(max_string) :: err_msg

      integer :: iobs, ivar, nvars, ivar_prsl, geo_ivar
      real(kind_real) :: wf
      integer :: wi, ierr, nlocs
      real(kind_real), dimension(:), allocatable :: pressure
      type(ufo_geoval), pointer :: prsl

      type ufo_geoval_ptr
         type(ufo_geoval), pointer :: ptr
      end type ufo_geoval_ptr
      type(ufo_geoval_ptr), dimension(:), allocatable :: vals

      character(len=MAXVARLEN), allocatable :: geovnames(:)

      ! check if nobs is consistent in geovals & hofx
      !if (geovals%nobs /= size(hofx)) then
      !  write(err_msg,*) myname_, ' error: nobs inconsistent!'
      !  call abor1_ftn(err_msg)
      !endif

      ! **********************************************************
      !                           STEP 1
      ! **********************************************************

      ! Retrieving the required variables names for this ObsOperator
      geovnames = ufo_vars_vnames(geovals%variables)
      nvars = size(geovnames)
    
      ! Checking if all required model variables are in geovals and get its pointer.
      ! also, locating the $var_prsl (vertical coordinates of model) from geovals for 
      ! vertical interpolation.
      ivar_prsl = -999
      allocate(vals(nvars))
      do ivar = 1, nvars
         call ufo_geovals_get_var(geovals, geovnames(ivar), vals(ivar)%ptr, status=ierr)
         if (ierr/=0) then
            write(err_msg,*) myname_, " : ", trim(geovnames(ivar)), ' doesnt exist'
            call abor1_ftn(err_msg)
         endif
         if (trim(geovnames(ivar)) == trim(var_prsl)) ivar_prsl = ivar
      enddo
      
      if (ivar_prsl == -999 ) then
        write(err_msg,*) myname_, " : ", trim(var_prsl), ' is not in geovals'
        call abor1_ftn(err_msg)
      endif

      ! **********************************************************
      !                           STEP 2
      ! **********************************************************

      ! Retrieving the observation vertical coordinate from ObsSpace.
      ! Different observation type variables may have different vertical 
      ! coordinate vector, because of the missing values layout.

      ! Because the vertical coordinate information is likely different 
      ! for different variabls, here, we use the maximum number to allocat
      ! the pressure 
      nlocs = obsspace_get_nlocs(obss)
      allocate(pressure(nlocs))

      do ivar = 1, self%nvars
        ! Get the vertical coordinate and its dimension for this variable
        ! To be revisited, should not use this hard-wired name
        call obsspace_get_db(obss, "MetaData", "air_pressure", pressure)

        ! Determine the location of this variable in geovals
        if (trim(self%varout(ivar)) == "air_temperature") then ! not match, to be solved
          geo_ivar = find_position("virtual_temperature", size(geovnames), geovnames)
        else
          geo_ivar = find_position(self%varout(ivar), size(geovnames), geovnames)
        endif
        if (geo_ivar == -999 ) then
          write(err_msg,*) myname_, " : ", trim(self%varout(ivar)), ' is not in geovals'
          call abor1_ftn(err_msg)
        endif

        ! Interpolation
        do iobs = 1, nlocs
          ! Calculate the interpolation weights 
          call vert_interp_weights(vals(ivar_prsl)%ptr%nval, log(DBLE(pressure(iobs))/10.), &
                                   vals(ivar_prsl)%ptr%vals(:,iobs), wi, wf)
          ! Interpolate from geovals to observational location.
          call vert_interp_apply(vals(geo_ivar)%ptr%nval, vals(geo_ivar)%ptr%vals(:,iobs), hofx(ivar+(iobs-1)*self%nvars), wi, wf)
        enddo

      enddo

      ! cleanup
      if (allocated(geovnames)) deallocate(geovnames)
      if (allocated(pressure)) deallocate(pressure)
      if (allocated(vals)) deallocate(vals)
    
    end subroutine conventional_profile_simobs_

! ------------------------------------------------------------------------------
    
    subroutine  destructor(self)
      type(ufo_conventional_profile), intent(inout)  :: self
      if (allocated(self%varout)) deallocate(self%varout)
    end subroutine destructor

! ------------------------------------------------------------------------------

    integer function find_position(str, length, strList)
      implicit none
      integer, intent(in) :: length
      character(len=*), intent(in) :: str
      character(len=*), dimension(length), intent(in) :: strList

      integer :: i

      find_position = -999
      do i = 1, length
        if (trim(str) == trim(strList(i))) find_position = i
      enddo
    end function find_position

! ------------------------------------------------------------------------------

end module ufo_conventional_profile_mod
