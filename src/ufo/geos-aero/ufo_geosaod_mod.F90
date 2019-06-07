! (C) Copyright 2017-2018 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

!> Fortran module for geosaod observation operator

module ufo_geosaod_mod

 use iso_c_binding
 use config_mod
 use kinds

 use ufo_geovals_mod, only: ufo_geovals, ufo_geoval, ufo_geovals_get_var
 use ufo_basis_mod, only: ufo_basis
 use ufo_vars_mod
 use obsspace_mod

 use GEOS_MieObs_mod
 

 implicit none
 private
 integer, parameter :: max_string=800

 type, extends(ufo_basis), public :: ufo_geosaod
 private
   integer, public :: nvars_in, n_wavelengths
   character(len=max_string), public, allocatable :: varin(:)
   character(len=max_string), public, allocatable :: varout(:)
   real(kind_real), public, allocatable :: wavelength(:)
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

character(len=*) :: var_name
character(len=3) :: wav
 
  !  varin: aer tracer variables we need from the model (list in .yaml file)
  !  need also relative humidity and delp (in varindefault).
  !---------
  tracer_variables = config_get_string_vector(c_conf, max_string, "tracer_geovals")
  self%n_tracers = size(tracer_variables)
  self%nvars_in =  size(varin_default) + self%n_tracers 

  allocate(self%varin(self%nvars_in))

  do ii = 1, self%nvars_in 
     self%varin(ii) = tracers_variables(ii)                      ! aer MR
  enddo
  self%varin(self%n_tracers + 1 : self%nvars_in) = varindefault  ! delp and rh
 ! I am not sure to keep delp and RH in self%varin ??? only tracers perhaps since we
 ! get delp and RH from geovals

  !varout: variables in the observation vector
  !------
  self%n_wavelengths = size(config_get_float_vector(c_conf, "wavelengths")) 
                        ! number of wavelengths at which AOD will be computed  

  self%wavelength = allocate(self%n_wavelengths)                 ! wavelength numb
  self%wavelength = config_get_float_vector(c_conf, "wavelengths")
 
  allocate(self%varout(self%n_wavelengths))
  ! Read variable list and store in varout                       ! AOD for now 
  !------- 
  var_name = config_get_string_vector(c_conf, max_string, "variables") 
  do jj = 1, self%n_wavelengths     
     write(wav, 'IO)') int(self%wavelength(jj))
     self%varout = var_name //'_'// trim(wav)   !name: aerosol_optical_depth_in_log_space_550 (for ex)
  enddo

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
subroutine ufo_geosaod_simobs(self, geovals, nvars, nlocs, hofx, obss)

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
integer :: nlayers, rc
real(kind_real), dimension(:), allocatable :: obss_metadata  ! check real kind

character(len=10), parameter :: rcfile = 'Aod_EOS.rc'   ! perhaps move it into geos-aero/
real(kind_real)  , parameter :: grav = 9.80616
type(Chem_Mie)               :: mieTables

real(c_double), dimension(:,:,:), allocatable :: qm  ! aer mass mix ratio (kg/kg *delp/g) prof at obs loc
real(c_double), dimension(:,:), allocatable   :: rh  ! relativ humidity prof interp at obs loc
real(c_double), dimension(:,:), allocatable   :: delp

character(len=MAXVARLEN) :: geovar

    
  ! Get delp and rh from model interp at obs loc (from geovals)
  call ufo_geovals_get_var(geovals, var_delp, delp_profile)
  nlayers = profile%nval                                          ! number of model layers
  allocate(delp(nlayers,nlocs))
  delp = profile%vals

  ! Get RH from geovals
  allocate(rh(nlayers,nlocs))
  call ufo_geovals_get_var(geovals, var_RH, rh_profile)
  rh = profile%vals

  ! Get Aer profiles interpolated at obs loc 
  allocate(qm(self%n_tracers, nlayers, nlocs))   
  do n = 1, self%n_tracers
     geovar = self%varin(n)                   !self%varin in setup contains tracers first then delp and rh
     call ufo_geovals_get_var(geovals, geovar, aer_profile)  
     qm(n,:,:) = aer_profile(n)%vals
     qm(n,:,:) = qm(n,:,:) * delp / grav
  enddo
   
  ! create Mie tables
  ! -----------------
  call get_Mie_Tables(mieTables, rcfile, rc)

  ! call observation operator code
  ! -----------------------------   
  hofx(:,:) = 0.0_kind_real
  call get_GEOS_AOD(nlayers, nlocs, self%n_wavelengths, self%n_tracers, rcfile, mieTables, &
                    self%wavelength, self%varin(1:self%n_tracers), qm, rh,  hofx, rc = rc)  !self%varin includes rh and delp!!!!
        
  ! delete the Mie tables
  ! --------------------
  call detete_Mie_Tables(mieTables, rc)

  ! cleanup memory
  ! --------
  deallocate(qm, rh, delp)


end subroutine ufo_geosaod_simobs


! ------------------------------------------------------------------------------

end module ufo_geosaod_mod
