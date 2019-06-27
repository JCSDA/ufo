! (C) Copyright 2019 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

!> Fortran module for geosaod tl/ad observation operator

module ufo_geosaod_tlad_mod

 use iso_c_binding
 use config_mod
 use kinds

 use ufo_geovals_mod, only: ufo_geovals, ufo_geoval, ufo_geovals_get_var
 use ufo_vars_mod
 use obsspace_mod
 use ufo_constants_mod, only: grav

 use GEOS_MieObs_mod

 implicit none
 private
 integer, parameter :: max_string=800

 !> Fortran derived type for the tl/ad observation operator
 type, public :: ufo_geosaod_tlad
 private
  integer :: nvars_in, nvars_out, nlocs, nlayers, ntracers
  character(len=max_string), public, allocatable :: varin(:)
  real(4), allocatable :: bext(:,:,:,:)
  real(c_float), public, allocatable :: wavelength(:)
  character(len=10) :: rcfile
  real(4), dimension(:,:), allocatable :: delp(:,:)
 contains
  procedure :: setup  => ufo_geosaod_tlad_setup
  procedure :: delete  => ufo_geosaod_tlad_delete
  procedure :: settraj => ufo_geosaod_tlad_settraj
  procedure :: simobs_tl  => ufo_geosaod_simobs_tl
  procedure :: simobs_ad  => ufo_geosaod_simobs_ad
 end type ufo_geosaod_tlad

contains

! ------------------------------------------------------------------------------

subroutine ufo_geosaod_tlad_setup(self, c_conf, c_nvars_out)

implicit none
class(ufo_geosaod_tlad),  intent(inout) :: self
type(c_ptr),              intent(in)    :: c_conf
integer(c_int),           intent(in)    :: c_nvars_out

!Locals
integer :: iq
character(len=maxvarlen), allocatable :: tracer_variables(:)

 ! Let user choose specific aerosols needed.
 self%ntracers = size(config_get_string_vector(c_conf, max_string, "tracer_geovals"))
 allocate(tracer_variables(self%ntracers))
 tracer_variables = config_get_string_vector(c_conf, max_string, "tracer_geovals")

 self%nvars_in =  self%ntracers
 allocate(self%varin(self%nvars_in))
 do iq = 1, self%nvars_in
    self%varin(iq) = tracer_variables(iq)
 enddo

 deallocate(tracer_variables)

 ! List of wavelenths
 self%nvars_out = c_nvars_out
 allocate(self%wavelength(self%nvars_out))
 call config_get_float_vector(c_conf, "wavelengths", self%wavelength)

 ! RC File for ChemBase
 self%rcfile = config_get_string(c_conf,len(self%rcfile),"RCFile")

end subroutine ufo_geosaod_tlad_setup

! ------------------------------------------------------------------------------

subroutine ufo_geosaod_tlad_delete(self)

implicit none
class(ufo_geosaod_tlad), intent(inout) :: self

 if (allocated(self%varin))   deallocate(self%varin)
 if (allocated(self%bext))    deallocate(self%bext)
 if (allocated(self%delp))    deallocate(self%delp)

end subroutine ufo_geosaod_tlad_delete

! ------------------------------------------------------------------------------
subroutine ufo_geosaod_tlad_settraj(self, geovals, obss)

implicit none
class(ufo_geosaod_tlad), intent(inout) :: self
type(ufo_geovals),       intent(in)    :: geovals
type(c_ptr), value,      intent(in)    :: obss

! Local variables
type(ufo_geoval), pointer :: aer_profile
type(ufo_geoval), pointer :: delp_profile
type(ufo_geoval), pointer :: rh_profile
integer :: rc, n
character(len=MAXVARLEN) :: geovar

real(4), dimension(:,:),   allocatable :: rh(:,:)
real(4), dimension(:,:,:), allocatable :: qm    ! aer mass mix ratio (kg/kg *delp/g) prof at obs loc

 ! Get nlcos
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
 do n = 1, self%ntracers
    geovar = self%varin(n)                   !self%varin in setup contains tracers first then delp and rh
    call ufo_geovals_get_var(geovals, geovar, aer_profile)
    qm(n,:,:) = aer_profile%vals
    qm(n,:,:) = qm(n,:,:) * self%delp / grav
 enddo

 allocate(self%bext(self%nlayers, self%nvars_out, self%ntracers, self%nlocs))

 call get_GEOS_AOD(self%nlayers, self%nlocs, self%nvars_out, self%ntracers, self%rcfile,  &
                   real(self%wavelength,4), self%varin, qm, rh, bext=self%bext, rc = rc) 

 deallocate(rh)
 deallocate(qm)

end subroutine ufo_geosaod_tlad_settraj

! ------------------------------------------------------------------------------

subroutine ufo_geosaod_simobs_tl(self, geovals, obss, nvars, nlocs, hofx)

implicit none
class(ufo_geosaod_tlad), intent(in)    :: self
type(ufo_geovals),       intent(in)    :: geovals
type(c_ptr), value,      intent(in)    :: obss
integer,                 intent(in)    :: nvars, nlocs
real(c_double),          intent(inout) :: hofx(nvars, nlocs)

integer :: n
real(4) :: hofx4(nvars, nlocs)
real(4), dimension(:,:,:), allocatable :: qm_tl
type(ufo_geoval), pointer :: aer_profile

character(len=MAXVARLEN) :: geovar

 ! Get Aer profiles interpolated at obs loc
 allocate(qm_tl(self%ntracers, self%nlayers, nlocs))
 do n = 1, self%ntracers
    geovar = self%varin(n)
    call ufo_geovals_get_var(geovals, geovar, aer_profile)
    qm_tl(n,:,:) = aer_profile%vals
    qm_tl(n,:,:) = qm_tl(n,:,:) * self%delp / grav
 enddo

 call get_geos_aod_tl(self%nlayers, nlocs, self%nvars_out, self%ntracers, self%bext, qm_tl, aod_tot_tl=hofx4)

 ! Convert back to ufo precision
 ! -----------------------------
 hofx = real(hofx4,c_double)

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

real(4) :: hofx4(nvars, nlocs)
real(4), dimension(:,:,:), allocatable :: qm_ad


 allocate(qm_ad(self%ntracers, self%nlayers, nlocs))
 
 call get_geos_aod_ad(self%nlayers, nlocs, self%nvars_out, self%ntracers, self%bext, qm_ad, aod_tot_ad=hofx4)

 deallocate(qm_ad)

end subroutine ufo_geosaod_simobs_ad

! ------------------------------------------------------------------------------

end module ufo_geosaod_tlad_mod
