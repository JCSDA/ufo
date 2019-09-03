! (C) Copyright 2019 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

!> Fortran module for geosaod observation operator

module ufo_geosaod_mod

 use iso_c_binding
 use fckit_configuration_module, only: fckit_configuration
 use kinds

 use ufo_geovals_mod, only: ufo_geovals, ufo_geoval, ufo_geovals_get_var
 use ufo_basis_mod, only: ufo_basis
 use ufo_vars_mod
 use obsspace_mod

 use ufo_constants_mod, only: grav

 use GEOS_MieObs_mod

 implicit none
 private
 integer, parameter :: max_string=800

!> Fortran derived type for the observation type
 type, public :: ufo_geosaod
 private
   integer, public :: nvars_in, nvars_out, ntracers
   character(len=max_string), public, allocatable :: varin(:)
   character(len=max_string), public, allocatable :: varout(:)
   real(kind_real), public, allocatable :: wavelength(:)
   character(len=maxvarlen),public :: rcfile
 contains
   procedure :: setup  => ufo_geosaod_setup
   procedure :: simobs => ufo_geosaod_simobs
   final :: destructor
 end type ufo_geosaod

!> Default variables required from model
 character(len=maxvarlen), dimension(2), parameter :: varindefault = (/var_delp, var_rh/)

contains

! ------------------------------------------------------------------------------

subroutine ufo_geosaod_setup(self, f_conf, vars)
implicit none
class(ufo_geosaod), intent(inout) :: self
type(fckit_configuration), intent(in) :: f_conf
character(len=maxvarlen), dimension(:), intent(inout) :: vars

!Locals
integer :: iq
character(kind=c_char,len=MAXVARLEN), allocatable :: tracer_variables(:)
integer(c_size_t),parameter :: csize = MAXVARLEN
character(len=:), allocatable :: str

  self%nvars_out = size(vars)
  allocate(self%varout(self%nvars_out))
  self%varout = vars
 
  ! Let user choose specific aerosols needed in aod calculation.
  call f_conf%get_or_die("tracer_geovals",csize,tracer_variables)
  self%ntracers = f_conf%get_size("tracer_geovals")
   
  self%nvars_in =  size(varindefault) + self%ntracers
  allocate(self%varin(self%nvars_in))
  do iq = 1, self%ntracers
     self%varin(iq) = tracer_variables(iq)                       ! aer MR
  enddo
  self%varin(self%ntracers + 1 : self%nvars_in) = varindefault   ! delp and rh (for concentration)

  deallocate(tracer_variables)

  ! List of wavelengths for aod
  allocate(self%wavelength(self%nvars_out))
  call f_conf%get_or_die("wavelengths", self%wavelength)
  
  ! RC File for ChemBase
  call f_conf%get_or_die("RCFile",str)
  self%rcfile = str
  
end subroutine ufo_geosaod_setup

! ------------------------------------------------------------------------------

subroutine destructor(self)
implicit none
type(ufo_geosaod), intent(inout) :: self

  if (allocated(self%varout))     deallocate(self%varout)
  if (allocated(self%varin))      deallocate(self%varin)
  if (allocated(self%wavelength)) deallocate(self%wavelength)
 
end subroutine destructor

! ------------------------------------------------------------------------------

subroutine ufo_geosaod_simobs(self, geovals, obss, nvars, nlocs, hofx)
implicit none
class(ufo_geosaod), intent(in)    :: self
integer, intent(in)               :: nvars, nlocs
type(ufo_geovals),  intent(in)    :: geovals
real(c_double),     intent(inout) :: hofx(nvars, nlocs)
type(c_ptr), value, intent(in)    :: obss

! Local variables
type(ufo_geoval), pointer :: aer_profile
type(ufo_geoval), pointer :: delp_profile
type(ufo_geoval), pointer :: rh_profile
integer :: nlayers, rc, n

real(kind_real), dimension(:,:,:), allocatable :: qm    ! aer mass mix ratio(kg/kg) that becomes concentration (*delp/g) profiles at obs loc
real(kind_real), dimension(:,:),   allocatable :: rh    ! relative humidity profile interp at obs loc
real(kind_real), dimension(:,:),   allocatable :: delp  ! air pressure thickness profiles at obs loc

character(len=MAXVARLEN) :: geovar

  ! Get delp and rh from model interp at obs loc (from geovals)
  call ufo_geovals_get_var(geovals, var_delp, delp_profile)
  nlayers = delp_profile%nval                            ! number of model layers

  allocate(delp(nlayers,nlocs))
  delp = delp_profile%vals

  ! Get RH from geovals
  allocate(rh(nlayers,nlocs))
  call ufo_geovals_get_var(geovals, var_RH, rh_profile)
  rh = rh_profile%vals

  ! Get Aer profiles interpolated at obs loc
  allocate(qm(self%ntracers, nlayers, nlocs))
  do n = 1, self%ntracers
     geovar = self%varin(n)                   ! self%varin in setup contains aer tracers first then delp and rh
     call ufo_geovals_get_var(geovals, geovar, aer_profile)
     qm(n,:,:) = aer_profile%vals             ! mass mixing ratio
     qm(n,:,:) = qm(n,:,:) * delp / grav      ! aer concentration (kg/m2) 
  enddo

  ! call observation operator code
  ! -----------------------------
  hofx(:,:) = 0.0
  call get_GEOS_AOD(nlayers, nlocs, self%nvars_out, self%ntracers, self%rcfile,  &
                    real(self%wavelength,4), self%varin(1:self%ntracers), qm, rh,       &
                    aod_tot = hofx, rc = rc)  !self%varin includes rh and delp

  ! cleanup memory
  ! --------
  deallocate(qm, rh, delp)

end subroutine ufo_geosaod_simobs

! ------------------------------------------------------------------------------

end module ufo_geosaod_mod
