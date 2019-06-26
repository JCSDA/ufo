! (C) Copyright 2019 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

!> Fortran module for geosaod observation operator

module ufo_geosaod_mod

 use iso_c_binding
 use config_mod
 use kinds

 use ufo_geovals_mod, only: ufo_geovals, ufo_geoval, ufo_geovals_get_var
 use ufo_vars_mod
 use obsspace_mod

 use GEOS_MieObs_mod

 implicit none
 private
 integer, parameter :: max_string=800

 type, public :: ufo_geosaod
 private
   integer, public :: nvars_in, nwavelengths, ntracers
   character(len=max_string), public, allocatable :: varin(:)
   character(len=max_string), public, allocatable :: varout(:)
   real(c_float), public, allocatable :: wavelength(:)
 contains
   procedure :: setup  => ufo_geosaod_setup
   procedure :: delete => ufo_geosaod_delete
   procedure :: simobs => ufo_geosaod_simobs
 end type ufo_geosaod

 character(len=maxvarlen), dimension(2), parameter :: varindefault = &
                                                 (/var_delp, var_rh/)

contains

! ------------------------------------------------------------------------------
subroutine ufo_geosaod_setup(self, c_conf)
  implicit none
  class(ufo_geosaod), intent(inout) :: self
  type(c_ptr),        intent(in)    :: c_conf

  integer nvar_name
  character(len=*), allocatable :: var_name(:)
  character(len=3) :: wav
  character(len=*), allocatable :: tracer_variables(:)

  integer iq, n

  !  varin: aer tracer variables we need from the model (list in .yaml file)
  !  need also relative humidity and delp (in varindefault).
  !---------
  self%ntracers = size(config_get_string_vector(c_conf, max_string, "tracer_geovals"))
  allocate(tracer_variables(self%ntracers))
  tracer_variables = config_get_string_vector(c_conf, max_string, "tracer_geovals")
  self%nvars_in =  size(varindefault) + self%ntracers

  allocate(self%varin(self%nvars_in))

  do iq = 1, self%nvars_in
     self%varin(iq) = tracer_variables(iq)                       ! aer MR
  enddo
  self%varin(self%ntracers + 1 : self%nvars_in) = varindefault   ! delp and rh

  !varout: variables in the observation vector
  !------
  self%nwavelengths = config_get_int(c_conf, "n_wavelengths")
  allocate(self%wavelength(self%nwavelengths))
  call config_get_float_vector(c_conf, "wavelengths", self%wavelength)

  ! Read variable list and store in varout
  allocate(self%varout(self%nwavelengths))
  nvar_name = size( config_get_string_vector(c_conf, max_string, "variables"))  ! AOD for now
  allocate(var_name(nvar_name))
  var_name = config_get_string_vector(c_conf, max_string, "variables")
  do n = 1, self%nwavelengths
     write(wav, 'IO') int(self%wavelength(n))
     self%varout = var_name //'_'// trim(wav)   !name: aerosol_optical_depth_in_log_space_550 (for ex)
  enddo

  deallocate(var_name, tracer_variables)
end subroutine ufo_geosaod_setup

! ------------------------------------------------------------------------------
subroutine ufo_geosaod_delete(self)
implicit none
class(ufo_geosaod), intent(inout) :: self

  if (allocated(self%varout)) deallocate(self%varout)
  if (allocated(self%varin))  deallocate(self%varin)
  if (allocated(self%wavelength)) deallocate(self%wavelength)

end subroutine ufo_geosaod_delete

! ------------------------------------------------------------------------------
subroutine ufo_geosaod_simobs(self, geovals, obss, nvars, nlocs, hofx)

implicit none
class(ufo_geosaod), intent(in)    :: self
type(ufo_geovals),  intent(in)    :: geovals
integer,            intent(in)    :: nvars, nlocs
real(c_double),     intent(inout) :: hofx(nvars,nlocs)   ! nwavelength, nlocs
type(c_ptr), value, intent(in)    :: obss

! Local variables
character(*), parameter :: PROGRAM_NAME = ' ufo_geosaod_mod.F90'
type(ufo_geoval), pointer :: aer_profile
type(ufo_geoval), pointer :: delp_profile
type(ufo_geoval), pointer :: rh_profile
integer :: nlayers, rc, n

character(len=10), parameter :: rcfile = 'Aod_EOS.rc'   ! perhaps move it into geos-aero/
real(kind_real)  , parameter :: grav = 9.80616

real(c_double), dimension(:,:,:), allocatable :: qm  ! aer mass mix ratio (kg/kg *delp/g) prof at obs loc
real(c_double), dimension(:,:), allocatable   :: rh  ! relativ humidity prof interp at obs loc
real(c_double), dimension(:,:), allocatable   :: delp

character(len=MAXVARLEN) :: geovar


  ! Get delp and rh from model interp at obs loc (from geovals)
  call ufo_geovals_get_var(geovals, var_delp, delp_profile)
  nlayers = delp_profile%nval                                          ! number of model layers
  allocate(delp(nlayers,nlocs))
  delp = delp_profile%vals

  ! Get RH from geovals
  allocate(rh(nlayers,nlocs))
  call ufo_geovals_get_var(geovals, var_RH, rh_profile)
  rh = rh_profile%vals

  ! Get Aer profiles interpolated at obs loc
  allocate(qm(self%ntracers, nlayers, nlocs))
  do n = 1, self%ntracers
     geovar = self%varin(n)                   !self%varin in setup contains tracers first then delp and rh
     call ufo_geovals_get_var(geovals, geovar, aer_profile)
     qm(n,:,:) = aer_profile(n)%vals
     qm(n,:,:) = qm(n,:,:) * delp / grav
  enddo

  ! call observation operator code
  ! -----------------------------
  hofx(:,:) = 0.0_kind_real
  call get_GEOS_AOD(nlayers, nlocs, real(self%nwavelengths,4), self%ntracers, rcfile,  &
                    self%wavelength, self%varin(1:self%n_tracers), qm, rh,       &
                    aod_tot = hofx, rc = rc)  !self%varin includes rh and delp!!!!


  ! cleanup memory
  ! --------
  deallocate(qm, rh, delp)


end subroutine ufo_geosaod_simobs


! ------------------------------------------------------------------------------

end module ufo_geosaod_mod
