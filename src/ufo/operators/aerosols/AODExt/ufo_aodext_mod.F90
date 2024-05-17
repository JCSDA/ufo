! (C) Copyright 2021 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

!> Fortran module for aodext observation operator

module ufo_aodext_mod

 use iso_c_binding
 use kinds
 use oops_variables_mod
 use obs_variables_mod
 use fckit_log_module, only : fckit_log
 use ufo_vars_mod
 use missing_values_mod

 implicit none
 private

!> Fortran derived type for the observation type

 type, public :: ufo_aodext
 private
   type(obs_variables), public :: obsvars
   type(oops_variables), public :: geovars
   real(kind_real), public, allocatable :: wavelength(:)!(nprofiles)
   integer, public              :: nprofiles
 contains
   procedure :: setup  => ufo_aodext_setup
   procedure :: simobs => ufo_aodext_simobs
   final :: destructor
 end type ufo_aodext

!> Default variables required from model
 character(len=maxvarlen), dimension(2), parameter :: varindefault = (/var_delp, var_airdens/)
 character(len=maxvarlen), dimension(3), parameter :: extdefault = (/var_ext1, var_ext2, var_ext3/)
! needs at least 2 profiles of extinction in bkg for angstrom law
contains

!----------------
integer function b_channel( bracket, nprofiles, bkg_wavelengths, obs_wavelength ) ! given a wavelength, return the index of bracket wavelengths
implicit none
integer, intent(in)                   :: bracket         ! bracket index (1 or 2)
integer, intent(in)                   :: nprofiles       ! number of bkg profiles
real(kind_real), dimension(nprofiles) :: bkg_wavelengths ! array of bkg wavelengths
real(kind_real)                       :: obs_wavelength  ! observed wavelength

character(len=maxvarlen) :: err_msg
integer :: j
   j = 1
   do while(j < nprofiles)
       if(obs_wavelength >= bkg_wavelengths(j) .and.&
          obs_wavelength < bkg_wavelengths(j+1)) then
          if (bracket == 1) then
              b_channel = j
          else if (bracket == 2) then
              b_channel = j + 1
          else
              write(err_msg,*) "ufo_aodext_mod: function b_channel: bracket index should be 1 (lower) or 2 (upper)"
              call abor1_ftn(err_msg)
          endif
          j = j+ 1
       else if(obs_wavelength > bkg_wavelengths(j) .and.&
            obs_wavelength <= bkg_wavelengths(j+1)) then
          if (bracket == 1) then
              b_channel = j
          else if (bracket == 2) then
              b_channel = j + 1
          else
              write(err_msg,*) "ufo_aodext_mod: function b_channel: bracket index should be 1 (lower) or 2 (upper)"
              call abor1_ftn(err_msg)
          endif
          j = j +1
       endif
       j = j + 1
   enddo
   return
end function b_channel

! ------------------------------------------------------------------------------
subroutine ufo_aodext_setup(self, f_conf)
use fckit_configuration_module, only: fckit_configuration
implicit none
class(ufo_aodext), intent(inout)     :: self
type(fckit_configuration), intent(in) :: f_conf

!Locals
integer n
character(len=maxvarlen) :: err_msg

  ! Fill in geovars: input variables requested from the model
  ! Need slots for airdens and delp

   call f_conf%get_or_die("nprofiles", self%nprofiles)

   do n = 1, self%nprofiles
      call self%geovars%push_back(extdefault(n))
   enddo
   call self%geovars%push_back(varindefault)

  ! Wavelengths for extinction profiles, specified in yaml in croissant order
   allocate(self%wavelength(self%nprofiles))
   call f_conf%get_or_die("bkg_wavelengths", self%wavelength)
  !Check that the wavelengths are in an ascending order
   n = 1
   do while (n <  self%nprofiles)
      if(self%wavelength(n) > self%wavelength(n+1)) then
         write(err_msg,*) ' ufo_aodext_setup: bkg wavelengths should be in an ascending order'
         call abor1_ftn(err_msg)
      endif
      n = n + 1
   enddo

end subroutine ufo_aodext_setup

! ------------------------------------------------------------------------------
subroutine destructor(self)
implicit none
type(ufo_aodext), intent(inout) :: self

 if (allocated(self%wavelength)) deallocate(self%wavelength)

end subroutine destructor
! ------------------------------------------------------------------------------
subroutine ufo_aodext_simobs(self, geovals, obss, nvars, nlocs, hofx)
use kinds
use ufo_geovals_mod, only: ufo_geovals, ufo_geoval, ufo_geovals_get_var
use iso_c_binding
use obsspace_mod
use ufo_constants_mod, only: grav, zero

implicit none
class(ufo_aodext), intent(in)    :: self
integer, intent(in)               :: nvars, nlocs
type(ufo_geovals),  intent(in)    :: geovals
real(c_double),     intent(inout) :: hofx(nvars, nlocs)
type(c_ptr), value, intent(in)    :: obss

! Local variables
type(ufo_geoval), pointer :: ext_profile
type(ufo_geoval), pointer :: delp_profile
type(ufo_geoval), pointer :: airdens_profile

real(kind_real), dimension(:,:,:), allocatable :: ext      !(km, nlocs, nch) ext profiles interp at obs loc [km-1]
real(kind_real), dimension(:,:),   allocatable :: airdens  !(km, nlocs) airdens profiles interp at obs loc  [kg.m-3]
real(kind_real), dimension(:,:),   allocatable :: delp     !(km, nlocs) air pressure thickness profiles at obs loc[Pa]
real(kind_real), dimension(:,:),   allocatable :: aod_bkg  !(nlocs, nch) AOD computed from modeled ext profiles
real(kind_real), dimension(:), allocatable :: obss_wavelength ! (nvars) observed AOD wavelengths [nm]

real(kind_real) :: angstrom ! (Angstrom coefficient calculated from bkg AOD and wavelengths)
real(kind_real) :: logm

character(len=MAXVARLEN) :: geovar
real(c_double) :: missing

character(len=MAXVARLEN) :: message
integer :: nlayers
integer :: km, nobs, nch, ic, i, j, k

 ! Get airdens and delp and number of layers from geovals
 ! -----------------------
 geovar = self%geovars%variable(1)

 call ufo_geovals_get_var(geovals, var_delp, delp_profile)
 nlayers = delp_profile%nval   ! number of model layers

 allocate(delp(nlayers,nlocs))
 delp = delp_profile%vals

 call ufo_geovals_get_var(geovals, var_airdens, airdens_profile)
 allocate(airdens(nlayers,nlocs))
 airdens = airdens_profile%vals

 ! Get extinction profiles from geovals
 ! ---------------------
 allocate(ext(nlayers, nlocs, self%nprofiles))
 do nch = 1, self%nprofiles
    geovar = self%geovars%variable(nch)
    call ufo_geovals_get_var(geovals, geovar, ext_profile)
    ext(:,:,nch) = ext_profile%vals
 enddo

 ! Get some metadata from obsspace, observed AOD wavelengths
 ! -----------------------
 allocate(obss_wavelength(nvars))
 call obsspace_get_db(obss,"MetaData", "obs_wavelength", obss_wavelength)

 ! Check if observed wavelength AOD is within the range of bkg wavelength to apply angstrom law
 ! else hofx set to missing value
 do ic = 1, nvars
    if(obss_wavelength(ic) < self%wavelength(1) .or. obss_wavelength(ic) > self%wavelength(self%nprofiles)) then
       write(message,*) 'ufo_aodext_simobs: observed wavelength outside of bkg wavelengths range', obss_wavelength(ic)
       call fckit_log%info(message)
    endif
 enddo
 ! Observation operator
 ! ------------------

 ! Calculation of AOD from extinction profiles
 ! ------------------
 allocate(aod_bkg(nlocs, self%nprofiles))
 aod_bkg = zero
 do nch = 1, self%nprofiles
    do nobs = 1, nlocs
       do k =1, nlayers
        aod_bkg(nobs,nch) = aod_bkg(nobs, nch) + (ext(k,nobs,nch) * delp(k,nobs)/(airdens(k,nobs))/(grav*1000.0_kind_real))
       enddo
    enddo
 enddo

 ! hofx: angstrom law
 ! ------------------
 missing =missing_value(missing)
 hofx = zero
 do nobs = 1, nlocs
    do ic = 1, nvars
       if(obss_wavelength(ic) < self%wavelength(1) .or. obss_wavelength(ic) > self%wavelength(self%nprofiles)) then
           hofx(ic,nobs) = missing
       else
           i = b_channel(1, self%nprofiles, self%wavelength, obss_wavelength(ic))
           j = b_channel(2, self%nprofiles, self%wavelength, obss_wavelength(ic))
           logm = log(self%wavelength(i)/self%wavelength(j))
           angstrom = log(aod_bkg(nobs,i)/aod_bkg(nobs,j))/logm
           hofx(ic,nobs) = aod_bkg(nobs,i) * (obss_wavelength(ic)/self%wavelength(i))**angstrom
       endif
    enddo
 enddo
 deallocate(ext, airdens, delp, aod_bkg, obss_wavelength)

end subroutine ufo_aodext_simobs
! ------------------------------------------------------------------------------
end module ufo_aodext_mod
