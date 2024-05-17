! (C) Copyright 2017-2022 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

!> Fortran module for averaging kernel observation operator

module ufo_columnretrieval_mod

 use oops_variables_mod
 use obs_variables_mod
 use ufo_vars_mod
 use missing_values_mod
 use kinds
 use iso_c_binding

 implicit none
 private
 integer, parameter :: max_string=800

!> Fortran derived type for the observation type
 type, public :: ufo_columnretrieval
 private
   type(obs_variables), public :: obsvars
   type(oops_variables), public :: geovars
   integer :: nlayers_retrieval
   character(kind=c_char,len=:), allocatable :: obskernelvar, obspressurevar
   character(kind=c_char,len=:), allocatable :: tracervars, stretch
   logical :: isapriori, isaveragingkernel, totalnovertice
   real(kind_real) :: convert_factor_model
 contains
   procedure :: setup  => ufo_columnretrieval_setup
   procedure :: simobs => ufo_columnretrieval_simobs
   final :: destructor
 end type ufo_columnretrieval

contains

! ------------------------------------------------------------------------------
subroutine ufo_columnretrieval_setup(self, f_conf)
  use fckit_configuration_module, only: fckit_configuration
  implicit none
  class(ufo_columnretrieval), intent(inout)     :: self
  type(fckit_configuration), intent(in) :: f_conf
  integer :: nlevs_yaml
  integer :: ivar, nvars

  ! get configuration for the averaging kernel operator
  call f_conf%get_or_die("nlayers_retrieval", self%nlayers_retrieval)

  ! get name of geoval/tracer to use from the model
  call f_conf%get_or_die("tracer variables", self%tracervars)

  ! determine if the apriori term is needed
  call f_conf%get_or_die("isApriori", self%isapriori)

  ! determine if the averaging kernel is needed
  call f_conf%get_or_die("isAveragingKernel", self%isaveragingkernel)

  ! determine wht kind of vertice stretching we need
  call f_conf%get_or_die("stretchVertices", self%stretch)

  ! do we need a conversion factor, say between ppmv and unity?
  call f_conf%get_or_die("model units coeff", self%convert_factor_model)

  ! perform a simple total column calculation using model whole profile
  call f_conf%get_or_die("totalNoVertice", self%totalnovertice)

  ! add variables to geovars that are needed
  ! specified tracers
  call self%geovars%push_back(self%tracervars)
  ! column pressure at interface
  call self%geovars%push_back(var_prsi)

end subroutine ufo_columnretrieval_setup

! ------------------------------------------------------------------------------
subroutine destructor(self)
  implicit none
  type(ufo_columnretrieval), intent(inout) :: self

  if (allocated(self%obskernelvar)) deallocate(self%obskernelvar)
  if (allocated(self%obspressurevar)) deallocate(self%obspressurevar)
  if (allocated(self%tracervars)) deallocate(self%tracervars)
  if (allocated(self%stretch)) deallocate(self%stretch)

end subroutine destructor

! ------------------------------------------------------------------------------
! averaging kernel observation operator
subroutine ufo_columnretrieval_simobs(self, geovals_in, obss, nvars, nlocs, hofx)
  use kinds
  use ufo_geovals_mod, only: ufo_geovals, ufo_geoval, ufo_geovals_get_var, &
                             ufo_geovals_reorderzdir, ufo_geovals_copy
  use ufo_constants_mod, only: zero, one
  use satcolumn_mod, only: simulate_column_ob
  use iso_c_binding
  use obsspace_mod
  implicit none
  class(ufo_columnretrieval), intent(in)    :: self
  integer, intent(in)               :: nvars, nlocs
  type(ufo_geovals),  intent(in)    :: geovals_in
  real(c_double),     intent(inout) :: hofx(nvars, nlocs)
  type(c_ptr), value, intent(in)    :: obss

  ! Local variables
  type(ufo_geoval), pointer :: prsi, tracer
  integer :: iobs, ilev
  character(len=MAXVARLEN) :: varstring
  character(len=4) :: levstr
  real(kind_real), allocatable, dimension(:,:) :: avgkernel_obs, prsi_obs
  real(kind_real), allocatable, dimension(:) :: apriori_term
  real(kind_real) :: hofx_tmp
  type(ufo_geovals) :: geovals
  real(c_double) :: missing
  character(len=max_string) :: err_msg

  missing = missing_value(missing)

  !prevent impossible options
  if (self%totalnovertice .and. (self%nlayers_retrieval > 1 .or. &
                                   self%isaveragingkernel )) then
    write(err_msg, *) "Error: wrong combination of yaml options, &
                  & totalNoVertice, isAveragingKernel, nlayers_retrieval"
    call abor1_ftn(err_msg)
  end if


  ! get geovals of atmospheric pressure
  call ufo_geovals_copy(geovals_in, geovals)  ! dont want to change geovals_in
  call ufo_geovals_get_var(geovals, self%geovars%variable(nvars+1), prsi)

  ! getting the apriori term if applicable
  allocate(apriori_term(nlocs))
  ! set to zero
  apriori_term = zero
  if (self%isapriori) then
    call obsspace_get_db(obss, "RetrievalAncillaryData", &
                         "aprioriTerm", apriori_term)
  end if

  allocate(avgkernel_obs(self%nlayers_retrieval, nlocs))
  ! set to 1.0 if no ak provided
  avgkernel_obs = one
  if (self%isaveragingkernel) then
    call obsspace_get_db_2d(obss, "RetrievalAncillaryData", &
                            "averagingKernel", avgkernel_obs)
  end if

  ! get prsi_obs
  allocate(prsi_obs(self%nlayers_retrieval+1, nlocs))
  if (.not. self%totalnovertice) then
    call obsspace_get_db_2d(obss, "RetrievalAncillaryData", &
                          "pressureVertice", prsi_obs)
  else
    prsi_obs(1, :) = prsi%vals(1, :)
    prsi_obs(self%nlayers_retrieval+1, :) = prsi%vals(prsi%nval, :)
  end if

  ! nvars is always 1 for this operator
  call ufo_geovals_get_var(geovals, self%tracervars, tracer)


  do iobs = 1, nlocs
    if (avgkernel_obs(1,iobs) /= missing) then ! take care of missing obs
      call simulate_column_ob(self%nlayers_retrieval, tracer%nval, &
                              avgkernel_obs(:,iobs), &
                              prsi_obs(:,iobs), &
                              prsi%vals(:,iobs), &
                              tracer%vals(:,iobs)*self%convert_factor_model, &
                              hofx_tmp, self%stretch)
      hofx(nvars,iobs) = hofx_tmp + apriori_term(iobs)
    else
      hofx(nvars,iobs) = missing ! default if we are unable to compute averaging kernel
    end if
  end do

end subroutine ufo_columnretrieval_simobs

end module ufo_columnretrieval_mod
