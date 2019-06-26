! (C) Copyright 2017-2018 UCAR
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
 use GEOS_MieObs_mod    ! located in geos-aero/

 implicit none
 private
 integer, parameter :: max_string=800

 !> Fortran derived type for the tl/ad observation operator
 ! TODO: add to the below type what you need for your tl/ad observation operator
 !       this type can hold information on trajectory, for geosaod
 type, public :: ufo_geosaod_tlad
 private
  integer :: ntracers    ! nvars_in
  integer :: nlayers, nwavelengths, nlocs
  character(len=max_string), public, allocatable :: varin(:)   ! variables requested from the model
  real(kind_real), allocatable :: wavelength(:)   ! nobs=nlocs
  real(c_double), allocatable  :: delp(:,:)       ! nlayers, nobs
  real(c_double), allocatable  :: rh(:,:)         ! nlayers, nobs
!  real(c_double), allocatable  :: ext(:,:,:,:)    ! nlayers, nwavelengths, ntracers, nobs
  type(Chem_Mie) :: MieTables

 contains
  procedure :: setup  => ufo_geosaod_tlad_setup
  procedure :: delete  => ufo_geosaod_tlad_delete
  procedure :: settraj => ufo_geosaod_tlad_settraj
  procedure :: simobs_tl  => ufo_geosaod_simobs_tl
  procedure :: simobs_ad  => ufo_geosaod_simobs_ad
 end type ufo_geosaod_tlad


contains

! ------------------------------------------------------------------------------
subroutine ufo_geosaod_tlad_setup(self, c_conf)
implicit none
class(ufo_geosaod_tlad), intent(inout) :: self
type(c_ptr),              intent(in)    :: c_conf


  !  varin: aer tracer variables we need from the model (list in .yaml file)
  tracer_variables = config_get_string_vector(c_conf, max_string, "tracer_geovals")
  self%ntracers = size(tracer_variables)

  allocate(self%varin(self%ntracers))

  do ii = 1, self%ntracers
     self%varin(ii) = tracers_variables(ii)                      ! aer MR
  enddo

  ! Wavelengths for AOD
  self%nwavelengths = size(config_get_float_vector(c_conf, "wavelengths"))
                        ! number of wavelengths at which AOD will be computed

  self%wavelength = allocate(self%nwavelengths)                 ! wavelength number
  self%wavelength = config_get_float_vector(c_conf, "wavelengths")

end subroutine ufo_geosaod_tlad_setup

! ------------------------------------------------------------------------------
subroutine ufo_geosaod_tlad_delete(self)
implicit none
class(ufo_geosaod_tlad), intent(inout) :: self

  if (allocated(self%varin))   deallocate(self%varin)
  if (allocated(self%wavelength)) deallocate(self%wavelength)
  if (allocated(self%delp))   deallocate(self%delp)
  if (allocated(self%rh)) deallocate(self%rh)
!  if (allocated(self%ext)) deallocate(self%ext)

  ! delete Mie Tables
  ! --------
  if (allocated(self%mieTables)) then
     call delete_Mie_Tables(self%mieTables, rc)
  endif

end subroutine ufo_geosaod_tlad_delete

! ------------------------------------------------------------------------------

subroutine ufo_geosaod_tlad_settraj(self, geovals, obss)
implicit none
class(ufo_geosaod_tlad), intent(inout) :: self
type(ufo_geovals),       intent(in)    :: geovals
type(c_ptr), value,      intent(in)    :: obss

! local variables
type(ufo_geoval), pointer :: profile_delp
type(ufo_geoval), pointer :: profile_rh
type(ufo_geoval), pointer :: profile_aer

character(len=10), parameter :: rcfile = 'Aod.rc'   ! should be put in geos-aero perhaps??
real(kind_real), parameter ::  grav =9.80616

real(c_double), allocatable :: qm(:,:,:)

character(len=MAXVARLEN) :: geovar

integer :: n

   ! Make sure nothing already allocated
   call self%cleanup()     ! if I do that how does it know the number of tracers?

   ! Get the number of observations
   self%nlocs = obsspace_get_nlocs(obss)

   ! Get delp and rh profiles from geovals
   call ufo_geovals_get_var(geovals, var_delp, profile_delp)  ! delp
   self%nlayers = profile%nval                                ! number of layers
   allocate(self%delp(self%nlayers, self%nlocs))
   self%delp = profile%vals

   call ufo_geovals_get_var(geovals, var_rh, profile_rh)      ! rh
   allocate(self%rh(self%nlayers, self%nlocs))
   self%rh = profile%vals


   allocate(qm(self%ntracers, self%nlayers, self%nlocs))
   do n = 1, self%n_tracers
     geovar = self%varin(n)                           !self%varin contains tracers
     call ufo_geovals_get_var(geovals, geovar, profile_aer)
     qm(n,:,:) = profile_aer(n)%vals
     qm(n,:,:) = qm(n,:,:) * self%delp / grav
   enddo


   ! create Mie tables
   ! -----------------
   call get_Mie_Tables(self%mieTables, rcfile, rc)

   ! put extinction coefficient for each layer, tracers, wavelengths in traj
   ! -------
!   allocate(self%ext(self%nlayers,self%nwavelengths, self%ntracers, self%nlocs))
!   self%ext =  0.0_kind_real
!   call get_GEOS_AOD(self%nlayers, self%nlocs, self%nwavelengths, self%ntracers, rcfile, self%mieTables, &
!                    self%wavelength, self%varin(1:self%n_tracers),self%qm, self%rh, ext=self%ext, rc = rc)
   deallocate(qm)

end subroutine ufo_geosaod_tlad_settraj

! ------------------------------------------------------------------------------
subroutine ufo_geosaod_simobs_tl(self, geovals, obss, nlayers, nwavelengths, nlocs, hofx)
implicit none
class(ufo_geosaod_tlad), intent(in)    :: self
type(ufo_geovals),       intent(in)    :: geovals
type(c_ptr), value,      intent(in)    :: obss
integer,                 intent(in)    :: nlayers, nwavelengths, nlocs

real(c_double),          intent(inout) :: hofx(nlayers,nwavelengths,nlocs)

integer :: ivar
type(ufo_geoval), pointer :: profile_aer
character(len=MAXVARLEN)  :: geovar
real(kind_real), parameter ::  grav =9.80616

real(c_double), allocatable :: qm(:,:,:)


character(len=10), parameter :: rcfile = 'Aod.rc'   ! should be put in geos-aero perhaps??

   allocate(qm(self%ntracers, self%nlayers, self%nlocs))

   do ivar = 1, self%ntracers
   ! get the name of input var in geovals (aer profile and rh, delp)
   geovar = self%varin(ivar)

   ! Get profile for this var from geovals
   call ufo_geovals_get_var(geovals, geovar, profile_aer)
   qm(ivar,:,:) = profile_aer(ivar)%vals
   qm(ivar,:,:) = qm(ivar,:,:) * self%delp / grav
   enddo

   ! call observation operator code
   ! -----------------------------
   hofx(:,:,:) = 0.0_kind_real
   call get_GEOS_AOD(self%nlayers, self%nlocs, self%nwavelengths, self%ntracers, rcfile, self%mieTables, &
                    self%wavelength, self%varin, qm, self%rh,  aod=hofx, rc = rc)

   deallocate(qm)

end subroutine ufo_geosaod_simobs_tl

! ------------------------------------------------------------------------------
subroutine ufo_geosaod_simobs_ad(self, geovals,  ntracers, nlayers, nwavelengths, hofx, obss)
implicit none
class(ufo_geosaod_tlad), intent(in)    :: self
type(ufo_geovals),       intent(inout) :: geovals
real(c_double),          intent(in)    :: hofx(nlayers,nwavelengths,nlocs)
type(c_ptr), value,      intent(in)    :: obss
integer,                 intent(in)    :: ntracers, nlayers, nwavelengths

type(ufo_geoval), pointer :: profile_aer
character(len=MAXVARLEN)  :: geovar
integer :: ivar

real(c_double), allocatable :: qm_ad(:,:,:)
real(kind_real), parameter ::  grav =9.80616



    ! call adjoint
   ! -----------------------------
   allocate(qm_ad(nlayers,nwavelengths, nlocs))
   qm_ad(:, :,:) = 0.0_kind_real
   call get_GEOS_AOD_ad(self%nlayers, self%nlocs, self%nwavelengths, self%ntracers, rcfile, self%mieTables, &
                    self%wavelength, self%varin, self%rh,  hofx, qm_ad, rc = rc)


   do ivar = 1, self%ntracers
     ! get the name of input var in geovals (aer profile and rh, delp)
     geovar = self%varin(ivar)

     ! Get pointer to profile for this var in geovals
     call ufo_geovals_get_var(geovals, geovar, profile_aer)

     ! allocate geovals profile if not yet allocated

     if (.not. allocated(profile_aer%vals)) then
        profile_aer%nlocs = self%nlocs
        profile_aer%nval = self%nlayers
        allocate(profile_aer(ivar)%vals(profile_aer%nval, profile_aer%nlocs))
        profile_aer(ivar)%vals(:,:) =  0.0_kind_real
     endif

     profile_aer(ivar)%vals(:,:) = qm_ad(ivar,:,:) * grav/self%delp
   enddo

   deallocate(qm_ad)


end subroutine ufo_geosaod_simobs_ad

! ------------------------------------------------------------------------------

end module ufo_geosaod_tlad_mod
