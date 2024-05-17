! (C) Copyright 2021 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

!> Fortran module for aodext tl/ad observation operator

module ufo_aodext_tlad_mod

 use kinds
 use missing_values_mod
 use oops_variables_mod
 use obs_variables_mod
 use ufo_vars_mod

 implicit none
 private

 type, public :: ufo_aodext_tlad
 private
  type(obs_variables), public :: obsvars
  type(oops_variables), public :: geovars
  real(kind_real),    allocatable :: airdens(:,:)      !(km, nlocs) airdens profile interp at obs loc [kg.m-3]
  real(kind_real),    allocatable :: delp(:,:)         !(km, nlocs) air pressure thickness profiles at obs loc[Pa]
  real(kind_real),    allocatable :: ext(:,:,:)        !(km, nlocs, nch) extinction profiles at obs loc [km-1]
  real(kind_real),    allocatable :: obss_wavelength(:)!(nvars) observed AOD wavelengths [nm]
  real(kind_real), public, allocatable :: wavelength(:)!(nch) background extinction profile's wavelengths[nm]
  integer :: nlayers,  nprofiles

 contains
  procedure :: setup  => ufo_aodext_tlad_setup
  procedure :: settraj => ufo_aodext_tlad_settraj
  procedure :: simobs_tl  => ufo_aodext_simobs_tl
  procedure :: simobs_ad  => ufo_aodext_simobs_ad
  final :: destructor
 end type ufo_aodext_tlad

!> Default variables required from model
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
       if(obs_wavelength .GE. bkg_wavelengths(j) .and.&
          obs_wavelength .LT. bkg_wavelengths(j+1)) then
          if (bracket == 1) then
              b_channel = j
          else if (bracket == 2) then  
              b_channel = j + 1
          else
              write(err_msg,*) "ufo_aodext_mod: function b_channel: bracket index should be 1 (lower) or 2 (upper)"
              call abor1_ftn(err_msg)
          endif
          j = j+ 1
       else if(obs_wavelength .GT. bkg_wavelengths(j) .and.&
            obs_wavelength .LE. bkg_wavelengths(j+1)) then
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
!-----------------------------------
subroutine ufo_aodext_tlad_setup(self, f_conf)
use fckit_configuration_module, only: fckit_configuration
implicit none
class(ufo_aodext_tlad), intent(inout) :: self
type(fckit_configuration), intent(in)  :: f_conf

!Locals
integer n
character(len=maxvarlen) :: err_msg
  
  ! setup input variables varin (updated model variables)
   call f_conf%get_or_die("nprofiles", self%nprofiles)
   
   if (self%nprofiles < 2 .or. self%nprofiles > 3) then
      write(err_msg,*) 'ufo_aodext_tlad_setup: number of extinction profiles must be 2 or 3'
      call abor1_ftn(err_msg) 
   endif

   do n = 1, self%nprofiles
       call self%geovars%push_back(extdefault(n))
   enddo

  ! Wavelengths for bkg extinction profiles, specified in yaml
   allocate(self%wavelength(self%nprofiles))
   call f_conf%get_or_die("bkg_wavelengths", self%wavelength)      
  !Check that the wavelengths are in an ascending order
   n = 1
   do while (n <  self%nprofiles)
      if(self%wavelength(n) > self%wavelength(n+1)) then
         write(err_msg,*) ' ufo_aodext_tlad_setup: bkg wavelengths should be in an ascending order'
         call abor1_ftn(err_msg)
      endif
      n = n + 1
   enddo 
end subroutine ufo_aodext_tlad_setup

! ------------------------------------------------------------------------------
subroutine destructor(self)
implicit none
type(ufo_aodext_tlad), intent(inout) :: self
  
   if (allocated(self%airdens))         deallocate(self%airdens)
   if (allocated(self%delp))            deallocate(self%delp)
   if (allocated(self%ext))             deallocate(self%ext)
   if (allocated(self%obss_wavelength)) deallocate(self%obss_wavelength)
   if (allocated(self%wavelength))      deallocate(self%wavelength)

end subroutine destructor

! ------------------------------------------------------------------------------
subroutine ufo_aodext_tlad_settraj(self, geovals, obss)
use iso_c_binding
use ufo_geovals_mod, only: ufo_geovals, ufo_geoval, ufo_geovals_get_var
use obsspace_mod
implicit none
class(ufo_aodext_tlad), intent(inout) :: self
type(ufo_geovals),       intent(in)    :: geovals
type(c_ptr), value,      intent(in)    :: obss

!locals
 character(len=MAXVARLEN) :: geovar
 integer :: nch, nvars, nlocs

type(ufo_geoval), pointer :: ext_profile
type(ufo_geoval), pointer :: delp_profile
type(ufo_geoval), pointer :: airdens_profile

 ! Get number of locations
 nlocs = obsspace_get_nlocs(obss)

 ! Get the number of obs type, for AOD it is the number of wavelengths
 nvars = self%obsvars%nvars()

 ! Get airdens, delp, ext and number of layers from geovals
 call ufo_geovals_get_var(geovals, var_delp, delp_profile)
 self%nlayers = delp_profile%nval   ! number of model layers

 allocate(self%delp(self%nlayers,nlocs)) 
 self%delp = delp_profile%vals

 call ufo_geovals_get_var(geovals, var_airdens, airdens_profile)
 allocate(self%airdens(self%nlayers,nlocs)) 
 self%airdens = airdens_profile%vals
 
 allocate(self%ext(self%nlayers, nlocs, self%nprofiles))
 do nch = 1, self%nprofiles
    geovar = self%geovars%variable(nch)
    call ufo_geovals_get_var(geovals, geovar, ext_profile)
    self%ext(:,:,nch) = ext_profile%vals
 enddo    

 ! Get some metadata from obsspace, observed AOD wavelengths
 ! -----------------------
 allocate(self%obss_wavelength(nvars))  
 call obsspace_get_db(obss,"MetaData", "obs_wavelength", self%obss_wavelength)

end subroutine ufo_aodext_tlad_settraj

! ------------------------------------------------------------------------------
! Input geovals parameter represents dx for tangent linear model
subroutine ufo_aodext_simobs_tl(self, geovals, obss, nvars, nlocs, hofx)
use iso_c_binding
use kinds
use ufo_geovals_mod, only: ufo_geovals, ufo_geoval, ufo_geovals_get_var
use obsspace_mod
use ufo_constants_mod, only: grav, zero

implicit none
class(ufo_aodext_tlad), intent(in)     :: self
type(ufo_geovals),       intent(in)    :: geovals
integer,                 intent(in)    :: nvars, nlocs
real(c_double),          intent(inout) :: hofx(nvars, nlocs)
type(c_ptr), value,      intent(in)    :: obss
  
!locals
real(kind_real), dimension(:,:,:), allocatable :: ext_tl     !(km, nlocs, nch) extinction profiles perturbations
real(kind_real), dimension(:,:),   allocatable :: aod_bkg    !(nlocs, nch) AOD computed from modeled ext profiles
real(kind_real), dimension(:,:),   allocatable :: aod_bkg_tl !(nlocs, nch) AOD tangent linear 
real(kind_real), dimension(:,:), allocatable :: angstrom     !(nvars, nlocs) Angstrom coefficient calculated from bkg AOD and wavelength
real(kind_real), dimension(:,:), allocatable :: angstrom_tl  !(nvars, nlocs) Angstrom coefficient tangent linear
real(kind_real), dimension(:,:), allocatable :: logm         !(nvars, nlocs) denominator of Angstrom parameter    

real(kind_real) :: arg1, arg1_tl, tmp, coef, coef_tl

integer :: nch, nobs, ic, i, j, km

type(ufo_geoval), pointer :: ext_profile
character(len=MAXVARLEN) :: geovar
real(c_double) :: missing

 ! Get extinction profile perturbations interpolated at obs loc
 allocate(ext_tl(self%nlayers, nlocs, self%nprofiles))

 do nch = 1, self%nprofiles
    geovar = self%geovars%variable(nch)
    call ufo_geovals_get_var(geovals, geovar, ext_profile)
    ext_tl(:,:,nch) = ext_profile%vals
 enddo    

 allocate(aod_bkg(nlocs, self%nprofiles)) 
 allocate(aod_bkg_tl(nlocs, self%nprofiles)) 
 aod_bkg = zero
 aod_bkg_tl = zero

 do nch = 1, self%nprofiles
    do nobs = 1, nlocs
       do km = 1, self%nlayers
              
        aod_bkg(nobs,nch) = aod_bkg(nobs, nch) + (self%ext(km,nobs,nch) * self%delp(km,nobs)&
                            * 1./(self%airdens(km,nobs) * grav*1000.0_kind_real))
        aod_bkg_tl(nobs,nch) = aod_bkg_tl(nobs, nch) + (ext_tl(km,nobs,nch) * self%delp(km,nobs)&
                            * 1./(self%airdens(km,nobs) * grav*1000.0_kind_real))

       enddo
    enddo
 enddo
 
 allocate(angstrom_tl(nvars,nlocs))
 allocate(angstrom(nvars,nlocs))
 allocate(logm(nvars,nlocs))

 missing =missing_value(missing)
 hofx = zero
 angstrom_tl = zero

 do nobs = 1, nlocs
    do ic = 1, nvars
       
       if(self%obss_wavelength(ic) < self%wavelength(1) .or.&
          self%obss_wavelength(ic) > self%wavelength(self%nprofiles)) then
          hofx(ic, nobs) = missing
       else
          i = b_channel(1, self%nprofiles, self%wavelength, self%obss_wavelength(ic))
          j = b_channel(2, self%nprofiles, self%wavelength, self%obss_wavelength(ic))
     
          logm(ic, nobs) = log(self%wavelength(i)/self%wavelength(j))
          tmp = aod_bkg(nobs, i) / aod_bkg(nobs, j)
          arg1_tl = (aod_bkg_tl(nobs, i)-tmp*aod_bkg_tl(nobs,j))/aod_bkg(nobs,j) 
          arg1 = tmp
      
          angstrom_tl(ic, nobs) = arg1_tl/(logm(ic, nobs) * arg1)
          angstrom(ic, nobs) = log(arg1)/logm(ic,nobs)

          coef = (self%wavelength(i)/self%wavelength(j))**angstrom(ic, nobs)
          coef_tl = coef * log(self%wavelength(i)/self%wavelength(j))*angstrom_tl(ic,nobs)
          hofx(ic, nobs) = coef * aod_bkg_tl(nobs, i) + aod_bkg(nobs, i)* coef_tl
       endif 
   enddo
 enddo
 deallocate( angstrom_tl, angstrom, aod_bkg, aod_bkg_tl, ext_tl, logm)

end subroutine ufo_aodext_simobs_tl

! ------------------------------------------------------------------------------
subroutine ufo_aodext_simobs_ad(self, geovals, obss, nvars, nlocs, hofx)
use kinds
use iso_c_binding
use ufo_geovals_mod, only: ufo_geovals, ufo_geoval, ufo_geovals_get_var
use obsspace_mod
use ufo_constants_mod, only: grav, zero

implicit none
class(ufo_aodext_tlad), intent(in)    :: self
type(ufo_geovals),       intent(inout) :: geovals
integer,                 intent(in)    :: nvars, nlocs
real(c_double),          intent(in)    :: hofx(nvars, nlocs)
type(c_ptr), value,      intent(in)    :: obss


real(kind_real), dimension(:,:,:), allocatable :: ext_ad     !(km, nlocs, nch) Adjoint of ext profiles at obs loc
real(kind_real), dimension(:,:),   allocatable :: aod_bkg    !(nlocs, nch) AOD computed from modeled ext profiles
real(kind_real), dimension(:,:),   allocatable :: aod_bkg_ad !(nlocs, nch) Adjoint of AOD 
real(kind_real), dimension(:,:),   allocatable :: angstrom   !(nvars, nlocs) Angstrom coefficient 
real(kind_real), dimension(:,:),   allocatable :: angstrom_ad!(nvars, nlocs) Adjoint of Angstrom coefficient
real(kind_real), dimension(:,:),   allocatable :: logm       !(nvars, nlocs) Denominator of Angstrom coeff

real(kind_real) :: arg1, arg1_ad, tmp, tmp_ad, coef, coef_ad

integer :: nch, nobs, km, ic, i, j

character(len=MAXVARLEN) :: geovar
type(ufo_geoval), pointer :: ext_profile
real(c_double) :: missing

 allocate(aod_bkg(nlocs, self%nprofiles)) 
 aod_bkg = zero

 do nch = 1, self%nprofiles
    do nobs = 1, nlocs
       do km = 1, self%nlayers
              
        aod_bkg(nobs,nch) = aod_bkg(nobs, nch) + (self%ext(km,nobs,nch) * self%delp(km,nobs) &
                            * 1./(self%airdens(km,nobs)*grav*1000.0_kind_real))

       enddo
    enddo
 enddo
 
 missing = missing_value(missing)
 allocate(angstrom(nvars,nlocs))
 allocate(logm(nvars,nlocs))
 do nobs = 1, nlocs
    do ic = 1, nvars
       
       if(self%obss_wavelength(ic) < self%wavelength(1) .or.&
          self%obss_wavelength(ic) > self%wavelength(self%nprofiles)) then     
          angstrom(ic, nobs) = missing
       else
 
       i = b_channel(1, self%nprofiles, self%wavelength, self%obss_wavelength(ic))
       j = b_channel(2, self%nprofiles, self%wavelength, self%obss_wavelength(ic))

       logm(ic, nobs) = log(self%wavelength(i)/self%wavelength(j))

       angstrom(ic, nobs) = log(aod_bkg(nobs,i)/aod_bkg(nobs,j))/logm(ic, nobs)
       endif
    enddo
 enddo

 allocate(aod_bkg_ad(nlocs, self%nprofiles)) 
 allocate(angstrom_ad(nvars,nlocs))
 aod_bkg_ad = zero
 angstrom_ad = zero

 do nobs = nlocs, 1, -1
    do ic = nvars, 1, -1
      
       if( hofx(ic,nobs)/=missing.and.angstrom(ic,nobs)/=missing )then

       i = b_channel(1, self%nprofiles, self%wavelength, self%obss_wavelength(ic))
       j = b_channel(2, self%nprofiles, self%wavelength, self%obss_wavelength(ic))
       
       coef = (self%wavelength(i)/self%wavelength(j))**angstrom(ic, nobs)
       aod_bkg_ad(nobs, i) = aod_bkg_ad(nobs, i) + coef * hofx(ic, nobs)

       angstrom_ad(ic, nobs) = angstrom_ad(ic,nobs) + coef * log(self%wavelength(i)/self%wavelength(j)) &
                               * aod_bkg(nobs, i) * hofx(ic, nobs)
  
       j = b_channel(2, self%nprofiles, self%wavelength, self%obss_wavelength(ic))

       tmp = aod_bkg(nobs, i) / aod_bkg(nobs, j)
       tmp_ad = angstrom_ad(ic, nobs)/(aod_bkg(nobs,j)*tmp*logm(ic, nobs))
       angstrom_ad = zero
       aod_bkg_ad(nobs, i) = aod_bkg_ad(nobs,i) + tmp_ad
       aod_bkg_ad(nobs, j) = aod_bkg_ad(nobs,j) - tmp * tmp_ad
       endif
    enddo
 enddo

 allocate(ext_ad(self%nlayers, nlocs, self%nprofiles))
 ext_ad = zero
 
 do nch = self%nprofiles, 1, -1
    do nobs = nlocs, 1, -1
       do km = self%nlayers, 1, -1
          ext_ad(km, nobs, nch) = ext_ad(km, nobs, nch) + self%delp(km, nobs) &
                                  * aod_bkg_ad(nobs, nch) * 1./(self%airdens(km, nobs) * grav * 1000.)
       enddo
    enddo
 enddo

 !Get pointer to ext profiles in geovals and put adjoint of extinction profiles into geovals
 do nch = self%nprofiles, 1, -1
     
   geovar = self%geovars%variable(nch)
   call ufo_geovals_get_var(geovals, geovar, ext_profile)
   ext_profile%vals(:,:) = ext_ad(:,:,nch)
                         
 enddo
 deallocate(angstrom_ad, angstrom, aod_bkg_ad, aod_bkg, ext_ad)  

end subroutine ufo_aodext_simobs_ad

! ------------------------------------------------------------------------------

end module ufo_aodext_tlad_mod
