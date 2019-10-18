! (C) Copyright 2019 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

!> Fortran module for geosaod tl/ad observation operator

module ufo_geosaod_tlad_mod

 use iso_c_binding

 use kinds

 use ufo_geovals_mod, only: ufo_geovals, ufo_geoval, ufo_geovals_get_var
 use ufo_vars_mod
 use ufo_constants_mod, only: grav

 use GEOS_MieObs_mod
 use oops_variables_mod


 implicit none
 private
 integer, parameter :: max_string=800

 !> Fortran derived type for the tl/ad observation operator
 type, public :: ufo_geosaod_tlad
  integer :: nlocs, nlayers, ntracers, nvars
  type(oops_variables), public :: obsvars
  type(oops_variables), public :: geovars
  real(kind=kind_real), allocatable :: bext(:,:,:,:)
  real(kind=kind_real), public, allocatable :: wavelength(:)
  character(len=maxvarlen),public:: rcfile
  real(kind=kind_real), dimension(:,:), allocatable :: delp(:,:)
 contains
  procedure :: setup  => ufo_geosaod_tlad_setup
  procedure :: delete => ufo_geosaod_tlad_delete
  procedure :: settraj => ufo_geosaod_tlad_settraj
  procedure :: simobs_tl  => ufo_geosaod_simobs_tl
  procedure :: simobs_ad  => ufo_geosaod_simobs_ad
 end type ufo_geosaod_tlad

contains

! ------------------------------------------------------------------------------

subroutine ufo_geosaod_tlad_setup(self, f_conf)
use fckit_configuration_module, only: fckit_configuration
implicit none
class(ufo_geosaod_tlad),  intent(inout) :: self
type(fckit_configuration), intent(in)   :: f_conf

!Locals
integer :: iq
character(len=maxvarlen), allocatable :: tracer_variables(:)
integer(c_size_t),parameter :: csize = MAXVARLEN
character(len=:), allocatable :: str

 ! Let user choose specific aerosols needed.

 call f_conf%get_or_die("tracer_geovals",csize,tracer_variables)
 self%ntracers = f_conf%get_size("tracer_geovals")

 do iq = 1, self%ntracers
    call self%geovars%push_back(tracer_variables(iq))      ! aer MR
 enddo
 deallocate(tracer_variables)

 ! size of variables (number of obs type (wavelength for AOD))
 self%nvars = self%obsvars%nvars()

 allocate(self%wavelength(self%nvars))
 call f_conf%get_or_die("wavelengths", self%wavelength)

 ! RC File for ChemBase
 call f_conf%get_or_die("RCFile",str)
 self%rcfile = str

end subroutine ufo_geosaod_tlad_setup

! ------------------------------------------------------------------------------

subroutine ufo_geosaod_tlad_delete(self)

implicit none
class(ufo_geosaod_tlad), intent(inout) :: self
 if (allocated(self%bext))    deallocate(self%bext)
 if (allocated(self%delp))    deallocate(self%delp)
 if (allocated(self%wavelength)) deallocate(self%wavelength)
end subroutine ufo_geosaod_tlad_delete

! ------------------------------------------------------------------------------

subroutine ufo_geosaod_tlad_settraj(self, geovals, obss)
use obsspace_mod
implicit none
class(ufo_geosaod_tlad), intent(inout) :: self
type(ufo_geovals),       intent(in)    :: geovals
type(c_ptr), value,      intent(in)    :: obss

! Local variables
type(ufo_geoval), pointer :: aer_profile
type(ufo_geoval), pointer :: delp_profile
type(ufo_geoval), pointer :: rh_profile
integer :: rc, iq
character(len=MAXVARLEN) :: geovar
character(len=MAXVARLEN), dimension(:), allocatable:: tracer_name

real(kind=kind_real), dimension(:,:),   allocatable :: rh(:,:)
real(kind=kind_real), dimension(:,:,:), allocatable :: qm    ! aer concentration (kg/kg *delp/g) profiles at obs loc


 ! Get number of locations
 self%nlocs = obsspace_get_nlocs(obss)

 ! Get delp and rh from model interp at obs loc (from geovals)
 call ufo_geovals_get_var(geovals, var_delp, delp_profile)
 self%nlayers = delp_profile%nval                                          ! number of model layers
 allocate(self%delp(self%nlayers,self%nlocs))
 self%delp = delp_profile%vals

 ! Get RH from geovals
 allocate(rh(self%nlayers,self%nlocs))
 call ufo_geovals_get_var(geovals, var_RH, rh_profile)
 rh = rh_profile%vals

 ! Get Aer profiles interpolated at obs loc
 allocate(qm(self%ntracers, self%nlayers, self%nlocs))
 allocate(tracer_name(self%ntracers))
 do iq = 1, self%ntracers
    geovar = self%geovars%variable(iq)                   !self%geovars contains tracers 
    tracer_name(iq) = geovar
    call ufo_geovals_get_var(geovals, geovar, aer_profile)
    qm(iq,:,:) = aer_profile%vals
    qm(iq,:,:) = qm(iq,:,:) * self%delp / grav
 enddo

 allocate(self%bext(self%nlayers, self%nvars, self%ntracers, self%nlocs)) !mass extinction efficiency 
 call get_GEOS_AOD(self%nlayers, self%nlocs, self%nvars, self%ntracers, self%rcfile,  &
                   real(self%wavelength,4), tracer_name, qm, rh, ext=self%bext, rc = rc) 

 deallocate(rh)
 deallocate(qm)
 deallocate(self%wavelength)

end subroutine ufo_geosaod_tlad_settraj

! ------------------------------------------------------------------------------

subroutine ufo_geosaod_simobs_tl(self, geovals, obss, nvars, nlocs, hofx)

implicit none
class(ufo_geosaod_tlad), intent(in)    :: self
type(ufo_geovals),       intent(in)    :: geovals
type(c_ptr), value,      intent(in)    :: obss
integer,                 intent(in)    :: nvars, nlocs
real(c_double),          intent(inout) :: hofx(nvars, nlocs)

integer :: iq
real(kind_real), dimension(:,:,:), allocatable :: qm_tl
type(ufo_geoval), pointer :: aer_profile

character(len=MAXVARLEN) :: geovar

 ! Get Aer profiles interpolated at obs loc
 allocate(qm_tl(self%ntracers, self%nlayers, nlocs))

 do iq = 1, self%ntracers
    geovar = self%geovars%variable(iq)                      !self%geovars contains tracers 
    call ufo_geovals_get_var(geovals, geovar, aer_profile)
    qm_tl(iq,:,:) = aer_profile%vals                         ! aer mass mixing ratio
    qm_tl(iq,:,:) = qm_tl(iq,:,:) * self%delp / grav          ! aer concentration
 enddo

 call get_geos_aod_tl(self%nlayers, nlocs, nvars, self%ntracers, self%bext, qm_tl, aod_tot_tl=hofx)

 deallocate(qm_tl)

end subroutine ufo_geosaod_simobs_tl

! ------------------------------------------------------------------------------

subroutine ufo_geosaod_simobs_ad(self, geovals, obss, nvars, nlocs, hofx)

implicit none
class(ufo_geosaod_tlad), intent(in)    :: self
type(ufo_geovals),       intent(inout) :: geovals
type(c_ptr), value,      intent(in)    :: obss
integer,                 intent(in)    :: nvars, nlocs
real(c_double),          intent(in)    :: hofx(nvars, nlocs)

integer :: iq
real(kind_real), dimension(:,:,:), allocatable :: qm_ad
type(ufo_geoval), pointer :: aer_profile
character(len=MAXVARLEN) :: geovar

 allocate(qm_ad(self%ntracers, self%nlayers, nlocs))   

 call get_geos_aod_ad(self%nlayers, nlocs, nvars, self%ntracers, self%bext, hofx, qm_ad)
 
 do iq = self%ntracers,1,-1

   geovar = self%geovars%variable(iq)                   !self%geovars contains tracers 
   call ufo_geovals_get_var(geovals, geovar, aer_profile)
   if (.not. allocated(aer_profile%vals)) then
       aer_profile%nlocs = nlocs
       aer_profile%nval  = self%nlayers
       allocate(aer_profile%vals(aer_profile%nval, aer_profile%nlocs))
       aer_profile%vals(:,:) = 0.0_kind_real
   endif

   qm_ad(iq,:,:) = qm_ad(iq,:,:) * self%delp / grav
   aer_profile%vals = qm_ad(iq,:,:)

 enddo
 
 deallocate(qm_ad)

end subroutine ufo_geosaod_simobs_ad

! ------------------------------------------------------------------------------


end module ufo_geosaod_tlad_mod
