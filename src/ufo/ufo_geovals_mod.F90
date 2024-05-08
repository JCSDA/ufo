!
! (C) Copyright 2017-2018 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
!
module ufo_geovals_mod

use iso_c_binding
use ufo_vars_mod
use kinds
use obsspace_mod
use missing_values_mod

use fckit_mpi_module, only: fckit_mpi_comm, fckit_mpi_sum

implicit none
private
integer, parameter :: max_string=800

public :: ufo_geovals, ufo_geoval, ufo_index_range
public :: ufo_geoval_sampled, ufo_geoval_reduced
public :: ufo_geovals_get_var, ufo_geovals_get_var_by_index
public :: ufo_geovals_default_constr, ufo_geovals_setup, ufo_geovals_partial_setup, ufo_geovals_delete
public :: ufo_geovals_setup_sampling_method, ufo_geovals_setup_trivial_sampling_method
public :: ufo_geovals_add_reduced_vars, ufo_geovals_get_vars
public :: ufo_geovals_get_default_format, ufo_geovals_set_default_format
public :: ufo_geovals_are_reduced_and_sampled_formats_aliased
public :: ufo_geovals_zero, ufo_geovals_random, ufo_geovals_scalmult
public :: ufo_geovals_allocate, ufo_geovals_print
public :: ufo_geovals_profmult
public :: ufo_geovals_reorderzdir
public :: ufo_geovals_assign, ufo_geovals_add, ufo_geovals_diff, ufo_geovals_abs
public :: ufo_geovals_minmaxavg, ufo_geovals_normalize, ufo_geovals_maxloc, ufo_geovals_schurmult
public :: ufo_geovals_read_netcdf, ufo_geovals_write_netcdf
public :: ufo_geovals_rms, ufo_geovals_copy, ufo_geovals_copy_one
public :: ufo_geovals_fill, ufo_geovals_fillad
public :: ufo_geovals_analytic_init

! ------------------------------------------------------------------------------

!> Half-open range of indices.
!>
!> The range consists of the indices i fulfilling the conditions begin <= i and i < end.
type, bind(c) :: ufo_index_range
  integer(c_size_t) :: begin
  integer(c_size_t) :: end
end type ufo_index_range

!> A method of sampling observation locations with a set of interpolation paths.
!>
!> Unlike in ufo_sampled_locations, here we aren't interested in the geometry of these paths
!> (their latitude, longitude etc.), but just in their mapping to obs locations.
type :: ufo_location_sampling_method
  !> Whether each location is sampled exactly once: by the interpolation path with the same index
  !> as the location.
  logical :: trivial = .false.
  !> Number of interpolation paths produced by this sampling method.
  integer :: npaths = 0
  !> Maps the index of each observation to the range of indices of the interpolation paths sampling
  !> its location (i.e. the region of space probed by that observation).
  type(ufo_index_range), allocatable :: paths_by_loc(:)
end type ufo_location_sampling_method

!> holds values of a model variable interpolated along a set of paths
type :: ufo_geoval
  real(kind_real), allocatable :: vals(:,:) !< values (nval, nprofiles)
  integer :: nval = 0               !< number of values in a profile (collected along a single path)
  integer :: nprofiles = 0          !< number of profiles
end type ufo_geoval

!> GeoVaL formats (consistent with the GeoVaLFormat C++ enum)
integer(c_int), parameter :: ufo_geoval_sampled = 1
integer(c_int), parameter :: ufo_geoval_reduced = 2

!> type to hold interpolated fields required by the obs operators
type :: ufo_geovals
  integer :: nlocs  = 0          !< number of observations
  integer :: nvar  = 0           !< number of elements of the `geovals` array
                                 !< (note: a given variable may be stored twice: once in the
                                 !< sampled and once in the reduced form)
  !> the variable format used when the optional parameter `format` of various `ufo_geovals_*`
  !> subroutines is not present
  integer(c_int) :: default_format = ufo_geoval_sampled

  type(ufo_geoval), allocatable :: geovals(:)  !< array of interpolated
                                               !< vertical profiles for all obs (nvar)

  ! The following two arrays are of length `nvar` and index elements of `geovals`.
  ! If `variables(i)` is nonempty and set to `varname`, `geovals(i)` stores the values of variable
  ! `varname` in the sampled format (with one profile per interpolation path).
  ! Similarly, if `reduced_variables(i)` is nonempty and set to `varname`, `geovals(i)` stores the
  ! values of variable `varname` in the reduced format (with one profile per location).
  ! Note that the sampled and reduced representations of a variable may be identical and stored
  ! in the same element of `geovals`.
  character(len=MAXVARLEN), allocatable :: variables(:)
  character(len=MAXVARLEN), allocatable :: reduced_variables(:)

  !< number of available methods of sampling observation locations with sets of interpolation paths
  integer :: nsampling_methods = 0
  type(ufo_location_sampling_method), allocatable :: sampling_methods(:)
  !> sampling_method_by_var(i) is the index of the observation location sampling method that
  !> produced the set of paths along which the ith variable has been or will be interpolated
  integer, allocatable :: sampling_method_by_var(:)

  real(c_double) :: missing_value !< obsspace missing value mark

  logical :: linit = .false.     !< .true. if all the ufo_geoval arrays inside geovals
                                 !< were allocated and have data and if all sampling methods have
                                 !< been set up
end type ufo_geovals

interface move
  module procedure :: move_ufo_geoval, move_location_sampling_method
end interface move

interface resize
  module procedure :: resize_geovals, resize_location_sampling_methods, &
                      resize_varnames, resize_integers
end interface resize

! ------------------------------------------------------------------------------
contains
! ------------------------------------------------------------------------------

subroutine ufo_geovals_default_constr(self)
implicit none
type(ufo_geovals), intent(inout) :: self

self%nlocs = 0
self%missing_value = missing_value(0.0)
self%nvar = 0
self%nsampling_methods = 0
self%linit = .false.

end subroutine ufo_geovals_default_constr

!> ------------------------------------------------------------------------------
!> Sets self%linit to .true. if and only if all variables have been allocated and all location
!> sampling methods set up.
subroutine ufo_geovals_update_linit(self)
implicit none
type(ufo_geovals), intent(inout) :: self

logical :: linit
integer :: i

linit = .true.
do i = 1, self%nvar
  if (.not. allocated(self%geovals(i)%vals)) linit = .false.
enddo

do i = 1, self%nsampling_methods
  if (.not. allocated(self%sampling_methods(i)%paths_by_loc)) linit = .false.
enddo

self%linit = linit

end subroutine ufo_geovals_update_linit

! ------------------------------------------------------------------------------
!> Initializes GeoVaLs and allocates storage for variables \p vars.
!>
!> @param self
!>   The GeoVaLs object to initialize.
!> @param nlocs
!>   Number of observation locations.
!> @param vars
!>   Names of variables to be stored in the new GeoVaLs object in the sampled format.
!> @param nvars
!>   Number of variables.
!> @param nvals
!>   Array indicating how many values per interpolation path will be allocated for each variable.
!> @param nsampling_methods
!>   Number of distinct methods of sampling observation locations with sets of interpolation paths.
!> @param npaths_by_method
!>   An array of length `nsampling_methods` whose ith element is the number of paths in the set
!>   of interpolation paths obtained with the ith sampling method.
!> @param sampling_method_by_var
!>   An array of length `nvars` whose ith element is the index of the observation location
!>   sampling method producing the set of paths along which the ith variable will be
!>   interpolated. Valid values are integers from 1 to `nsampling_methods`.
!> @param reduced_vars
!>   Names of the variables to be stored in the new GeoVaLs in the reduced format.
!> @param nreduced_vars
!>   Number of variables in `reduced_vars`.
!> @param nreduced_vals
!>   Array of length `nreduced_vars` whose ith element indicates how many values per
!>   location will be stored in the GeoVaL corresponding to the ith variable.
!> @param is_sampling_method_trivial
!>   Array of length `nsampling_methods` whose ith element indicates whether the ith sampling method
!>   produces a set of interpolation paths sampling all locations exactly once in ascending order.
!>
!> This call must be followed by calls to ufo_geovals_setup_(trivial_)sampling_method() mapping
!> locations to interpolation paths.
subroutine ufo_geovals_setup(self, nlocs, vars, nvars, nvals, &
                             nsampling_methods, npaths_by_method, sampling_method_by_var, &
                             reduced_vars, nreduced_vars, nreduced_vals, &
                             is_sampling_method_trivial)
use oops_variables_mod
implicit none
type(ufo_geovals), intent(inout) :: self
type(oops_variables), intent(in) :: vars, reduced_vars
integer, intent(in) :: nlocs, nvars, nreduced_vars
integer(c_size_t), intent(in) :: nvals(nvars), nreduced_vals(nreduced_vars)
integer, intent(in) :: nsampling_methods
integer(c_size_t), intent(in) :: npaths_by_method(nsampling_methods)
integer(c_size_t), intent(in) :: sampling_method_by_var(nvars)
logical(c_bool), intent(in) :: is_sampling_method_trivial(nsampling_methods)

character(max_string) :: err_msg
integer :: ivar, ireduced_var, igeoval, imethod
integer(c_size_t) :: reduced_vars_sampling_method
integer, allocatable :: reduced_geovals(:)  ! Indices of geovals storing reduced variables

if (nvars /= vars%nvars()) then
  write(err_msg,*) "ufo_geovals_setup: mismatch in the number of variables (nvars = ", &
    nvars, ", vars%nvars() = ", vars%nvars(), ")"
  call abor1_ftn(err_msg)
endif
if (nreduced_vars /= reduced_vars%nvars()) then
  write(err_msg,*) "ufo_geovals_setup: mismatch in the number of reduced variables ", &
    "(nreduced_vars = ", nreduced_vars, ", reduced_vars%nvars() = ", reduced_vars%nvars(), ")"
  call abor1_ftn(err_msg)
endif
if (any(sampling_method_by_var < 1)) then
  write(err_msg,*) "ufo_geovals_setup: sampling method ", &
    minval(sampling_method_by_var), " does not exist"
  call abor1_ftn(err_msg)
endif
if (any(sampling_method_by_var > nsampling_methods)) then
  write(err_msg,*) "ufo_geovals_setup: sampling method ", &
    maxval(sampling_method_by_var), " does not exist"
  call abor1_ftn(err_msg)
endif
do imethod = 1, nsampling_methods
  if (is_sampling_method_trivial(imethod) .and. npaths_by_method(imethod) /= nlocs) then
    write(err_msg,*) "ufo_geovals_setup: trivial sampling methods must sample each location exactly once"
    call abor1_ftn(err_msg)
  endif
enddo

call ufo_geovals_delete(self)
self%nlocs = nlocs
self%missing_value = missing_value(self%missing_value)

self%nvar = vars%nvars()
allocate(reduced_geovals(reduced_vars%nvars()))
do ireduced_var = 1, reduced_vars%nvars()
  if (vars%has(reduced_vars%variable(ireduced_var))) then
    ivar = vars%find(reduced_vars%variable(ireduced_var))
    if (.not. is_sampling_method_trivial(sampling_method_by_var(ivar)) .or. &
        nvals(ivar) /= nreduced_vals(ireduced_var)) then
      ! The reduced version of this variable is different from the sampled version,
      ! so it needs to be stored separately
      self%nvar = self%nvar + 1
      reduced_geovals(ireduced_var) = self%nvar
    else
      ! The reduced and sampled versions of this variable are identical
      reduced_geovals(ireduced_var) = ivar
    endif
  else
    ! Only the reduced version of this variable will be stored; the sampled version won't
    self%nvar = self%nvar + 1
    reduced_geovals(ireduced_var) = self%nvar
  endif
enddo

reduced_vars_sampling_method = 0
if (reduced_vars%nvars() > 0) then
  do imethod = 1, nsampling_methods
    if (is_sampling_method_trivial(imethod)) then
      reduced_vars_sampling_method = imethod
      exit
    endif
  enddo
  if (reduced_vars_sampling_method == 0) then
    ! Reduced variables need to be mapped to a trivial sampling method, but none has been provided.
    ! We need to create one on our own.
    reduced_vars_sampling_method = nsampling_methods + 1
  endif
endif

allocate(self%geovals(self%nvar))
allocate(self%variables(self%nvar))
allocate(self%reduced_variables(self%nvar))
allocate(self%sampling_method_by_var(self%nvar))

self%variables(:) = ""
self%reduced_variables(:) = ""

! Handle the sampled versions of variables
do ivar = 1, vars%nvars()
  igeoval = ivar
  self%variables(igeoval) = vars%variable(ivar)
  call ufo_geovals_setup_geoval(self, igeoval, sampling_method_by_var(ivar), &
                                nvals(ivar), npaths_by_method(sampling_method_by_var(ivar)))
enddo
! Handle the reduced versions of variables non-aliased with the sampled versions
do ireduced_var = 1, reduced_vars%nvars()
  igeoval = reduced_geovals(ireduced_var)
  self%reduced_variables(igeoval) = reduced_vars%variable(ireduced_var)
  if (igeoval > vars%nvars()) then
    call ufo_geovals_setup_geoval(self, igeoval, reduced_vars_sampling_method, &
                                  nreduced_vals(ireduced_var), int(nlocs, kind=c_size_t))
  endif
enddo

if (reduced_vars_sampling_method == nsampling_methods + 1) then
  self%nsampling_methods = nsampling_methods + 1
else
  self%nsampling_methods = nsampling_methods
endif
allocate(self%sampling_methods(self%nsampling_methods))
do imethod = 1, nsampling_methods
  self%sampling_methods(imethod)%trivial = is_sampling_method_trivial(imethod)
  self%sampling_methods(imethod)%npaths = npaths_by_method(imethod)
enddo
if (reduced_vars_sampling_method == nsampling_methods + 1) then
  self%sampling_methods(nsampling_methods + 1)%trivial = .true.
  self%sampling_methods(nsampling_methods + 1)%npaths = nlocs
  call ufo_geovals_setup_trivial_sampling_method(self, nsampling_methods + 1)
endif

end subroutine ufo_geovals_setup

! ------------------------------------------------------------------------------
!> Deprecated, use ufo_geovals_setup instead.
!>
!> Initializes GeoVaLs for variables `vars` but does not allocate storage (the `vals` array in each
!> ufo_geoval).
!>
!> @param self
!>   The GeoVaLs object to initialize.
!> @param nlocs
!>   Number of observation locations.
!> @param vars
!>   Names of variables to be stored in the new GeoVaLs object in the sampled format.
!> @param nvars
!>   Number of variables.
!> @param nsampling_methods
!>   Number of distinct methods of sampling observation locations with sets of interpolation paths.
!> @param sampling_method_by_var
!>   An array of length `nvars` whose ith element is the index of the observation location
!>   sampling method producing the set of paths along which the ith variable will be
!>   interpolated. Valid values are integers from 1 to `nsampling_methods`.
!> @param reduced_vars
!>   Names of the variables to be stored in the new GeoVaLs in the reduced format.
!> @param is_sampling_method_trivial
!>   Array of length `nsampling_methods` whose ith element indicates whether the ith sampling method
!>   produces a set of interpolation paths sampling all locations exactly once in ascending order.
!>
!> This call must be followed by calls to ufo_geovals_setup_(trivial_)sampling_method() mapping
!> locations to interpolation paths and ufo_geovals_allocate() allocating specific GeoVaLs.
subroutine ufo_geovals_partial_setup(self, nlocs, vars, nvars, nsampling_methods, &
                                     sampling_method_by_var, &
                                     reduced_vars, is_sampling_method_trivial)
use oops_variables_mod
implicit none
type(ufo_geovals), intent(inout) :: self
type(oops_variables), intent(in) :: vars, reduced_vars
integer, intent(in) :: nlocs, nvars
integer, intent(in) :: nsampling_methods
integer(c_size_t), intent(in) :: sampling_method_by_var(nvars)
logical(c_bool), intent(in) :: is_sampling_method_trivial(nsampling_methods)

character(max_string) :: err_msg
integer :: igeoval, ivar, ireduced_var, imethod
integer(c_size_t) :: reduced_vars_sampling_method
integer, allocatable :: reduced_geovals(:)

if (nvars /= vars%nvars()) then
  write(err_msg,*) "ufo_geovals_partial_setup: mismatch in the number of variables (nvars = ", &
    nvars, ", vars%nvars() = ", vars%nvars(), ")"
  call abor1_ftn(err_msg)
endif
if (any(sampling_method_by_var < 1)) then
  write(err_msg,*) "ufo_geovals_partial_setup: sampling method ", &
    minval(sampling_method_by_var), " does not exist"
  call abor1_ftn(err_msg)
endif
if (any(sampling_method_by_var > nsampling_methods)) then
  write(err_msg,*) "ufo_geovals_partial_setup: sampling method ", &
    maxval(sampling_method_by_var), " does not exist"
  call abor1_ftn(err_msg)
endif

call ufo_geovals_delete(self)
self%nlocs = nlocs
self%missing_value = missing_value(self%missing_value)

self%nvar = vars%nvars()
allocate(reduced_geovals(reduced_vars%nvars()))
do ireduced_var = 1, reduced_vars%nvars()
  if (vars%has(reduced_vars%variable(ireduced_var))) then
    ivar = vars%find(reduced_vars%variable(ireduced_var))
    if (.not. is_sampling_method_trivial(sampling_method_by_var(ivar))) then
      ! The reduced version of this variable is different from the sampled version,
      ! so it needs to be stored separately
      self%nvar = self%nvar + 1
      reduced_geovals(ireduced_var) = self%nvar
    else
      ! The reduced and sampled versions of this variable are identical
      reduced_geovals(ireduced_var) = ivar
    endif
  else
    ! Only the reduced version of this variable will be stored; the sampled version won't
    self%nvar = self%nvar + 1
    reduced_geovals(ireduced_var) = self%nvar
  endif
enddo

reduced_vars_sampling_method = 0
if (reduced_vars%nvars() > 0) then
  do imethod = 1, nsampling_methods
    if (is_sampling_method_trivial(imethod)) then
      reduced_vars_sampling_method = imethod
      exit
    endif
  enddo
  if (reduced_vars_sampling_method == 0) then
    ! Reduced variables need to be mapped to a trivial sampling method, but none has been provided.
    ! We need to create one on our own.
    reduced_vars_sampling_method = nsampling_methods + 1
  endif
endif

allocate(self%geovals(self%nvar))
allocate(self%variables(self%nvar))
allocate(self%reduced_variables(self%nvar))
allocate(self%sampling_method_by_var(self%nvar))

self%variables(:) = ""
self%reduced_variables(:) = ""

! Handle the sampled versions of variables
do ivar = 1, vars%nvars()
  igeoval = ivar
  self%variables(igeoval) = vars%variable(ivar)
  call ufo_geovals_partial_setup_geoval(self, igeoval, sampling_method_by_var(ivar))
enddo
! Handle the reduced versions of variables non-aliased with the sampled versions
do ireduced_var = 1, reduced_vars%nvars()
  igeoval = reduced_geovals(ireduced_var)
  self%reduced_variables(igeoval) = reduced_vars%variable(ireduced_var)
  if (igeoval > vars%nvars()) then
    call ufo_geovals_partial_setup_geoval(self, igeoval, reduced_vars_sampling_method)
  endif
enddo

if (reduced_vars_sampling_method == nsampling_methods + 1) then
  self%nsampling_methods = nsampling_methods + 1
else
  self%nsampling_methods = nsampling_methods
endif
allocate(self%sampling_methods(self%nsampling_methods))
do imethod = 1, nsampling_methods
  self%sampling_methods(imethod)%trivial = is_sampling_method_trivial(imethod)
enddo
if (reduced_vars_sampling_method == nsampling_methods + 1) then
  self%sampling_methods(nsampling_methods + 1)%trivial = .true.
  self%sampling_methods(nsampling_methods + 1)%npaths = nlocs
  call ufo_geovals_setup_trivial_sampling_method(self, nsampling_methods + 1)
endif

end subroutine ufo_geovals_partial_setup

! ------------------------------------------------------------------------------
!> Deprecated. Rely on ufo_geovals_setup to allocate GeoVaLs instead.
!> Allocates GeoVaLs for \p vars variables with \p nlevels number of levels.
!> If the GeoVaLs for this variable were allocated before with different size,
!> aborts.
subroutine ufo_geovals_allocate(self, vars, nlevels)
use oops_variables_mod
implicit none
type(ufo_geovals), intent(inout) :: self
type(oops_variables), intent(in) :: vars
integer, intent(in) :: nlevels

integer :: ivar, ivar_gvals, ireducedvar_gvals
character(max_string) :: err_msg

do ivar = 1, vars%nvars()
  ! find the index of the geoval to be allocated
  ivar_gvals = ufo_vars_getindex(self%variables, vars%variable(ivar))
  if (ivar_gvals > 0) call allocate_geoval(ivar_gvals)

  ireducedvar_gvals = ufo_vars_getindex(self%reduced_variables, vars%variable(ivar))
  if (ireducedvar_gvals > 0) call allocate_geoval(ireducedvar_gvals)

  if (ivar_gvals <= 0 .and. ireducedvar_gvals <= 0) then
    write(err_msg,*) "ufo_geovals_allocate: ", trim(vars%variable(ivar)), " doesn't exist in geovals"
    call abor1_ftn(err_msg)
  endif
enddo

! check if all variables are now allocated and all sampling methods set up, and set self%linit
! accordingly
call ufo_geovals_update_linit(self)

contains

  subroutine allocate_geoval(igeoval)
    implicit none
    integer, intent(in) :: igeoval

    if (allocated(self%geovals(igeoval)%vals) .and. (self%geovals(igeoval)%nval /= nlevels)) then
      write(err_msg,*) "ufo_geovals_allocate: attempting to allocate already allocated geovals for ",       &
                       trim(vars%variable(ivar)), ". Previously allocated as ", self%geovals(igeoval)%nval, &
                       " levels; now trying to allocate as ", nlevels, " levels."
      call abor1_ftn(err_msg)
    ! only allocate if not already allocated
    elseif (.not. allocated(self%geovals(igeoval)%vals)) then
      self%geovals(igeoval)%nval  = nlevels
      allocate(self%geovals(igeoval)%vals(nlevels, self%geovals(igeoval)%nprofiles))
    endif
  end subroutine

end subroutine ufo_geovals_allocate

! ------------------------------------------------------------------------------
!> @brief Specify which interpolation paths produced by a given method sample which observation
!> locations.
!>
!> @param self
!>   The GeoVaLs object to modify.
!> @param sampling_method
!>   Index of the observation location sampling method to set up: an integer from 1 to
!>   `self%nsampling_methods`.
!> @param npaths
!>   Number of interpolation paths produced by this sampling method.
!> @param nlocs
!>   Number of observation locations. Must match self%nlocs.
!> @param paths_by_loc
!>   An array mapping the index of each observation location to the range of (1-based) indices of
!>   the paths sampling that location. Specifically, the ith location is deemed to be sampled by
!>   the paths with indices ranging from `paths_by_loc(i)%begin` up to but not including
!>   `paths_by_loc(i)%end`.
subroutine ufo_geovals_setup_sampling_method(self, sampling_method, npaths, nlocs, paths_by_loc)
implicit none
type(ufo_geovals), intent(inout) :: self
integer, intent(in) :: sampling_method, npaths, nlocs
type(ufo_index_range), intent(in) :: paths_by_loc(nlocs)

character(max_string) :: err_msg
integer :: ivar, iloc

if (sampling_method < 1 .or. sampling_method > self%nsampling_methods) then
  write(err_msg,*) "ufo_geovals_setup_sampling_method: sampling method ", sampling_method, &
    " does not exist"
  call abor1_ftn(err_msg)
endif
if (nlocs /= self%nlocs) then
  write(err_msg,*) "ufo_geovals_setup_sampling_method: mismatch in nlocs"
  call abor1_ftn(err_msg)
endif
if (self%sampling_methods(sampling_method)%trivial) then
  do iloc = 1, nlocs
    if (paths_by_loc(iloc)%begin /= iloc .or. paths_by_loc(iloc)%end /= iloc + 1) then
      write(err_msg,*) "ufo_geovals_setup_sampling_method: the method should be trivial but is not"
      call abor1_ftn(err_msg)
    endif
  enddo
endif

self%sampling_methods(sampling_method)%npaths = npaths
do ivar = 1, self%nvar
  if (self%sampling_method_by_var(ivar) == sampling_method) then
    self%geovals(ivar)%nprofiles = npaths
  endif
enddo

allocate(self%sampling_methods(sampling_method)%paths_by_loc(self%nlocs))
self%sampling_methods(sampling_method)%paths_by_loc(:) = paths_by_loc(:)

call ufo_geovals_update_linit(self)

end subroutine ufo_geovals_setup_sampling_method

! ------------------------------------------------------------------------------
!> @brief Designate an observation location sampling method as "trivial", i.e. one producing a set
!> of interpolation paths such that each location is sampled solely by the path with the same index.
!>
!> @param self
!>   The GeoVaLs object to modify.
!> @param sampling_method
!>   Index of the observation location sampling method to set up: an integer from 1 to
!>   `self%nsampling_methods`.
subroutine ufo_geovals_setup_trivial_sampling_method(self, sampling_method)
implicit none
type(ufo_geovals), intent(inout) :: self
integer, intent(in) :: sampling_method

character(max_string) :: err_msg
integer :: ivar, iloc

if (sampling_method < 1 .or. sampling_method > self%nsampling_methods) then
  write(err_msg,*) "ufo_geovals_setup_trivial_sampling_method: sampling method ", sampling_method, &
    " does not exist"
  call abor1_ftn(err_msg)
endif

self%sampling_methods(sampling_method)%npaths = self%nlocs
do ivar = 1, self%nvar
  if (self%sampling_method_by_var(ivar) == sampling_method) then
    self%geovals(ivar)%nprofiles = self%nlocs
  endif
enddo

allocate(self%sampling_methods(sampling_method)%paths_by_loc(self%nlocs))
do iloc = 1, self%nlocs
  self%sampling_methods(sampling_method)%paths_by_loc(iloc)%begin = iloc
  self%sampling_methods(sampling_method)%paths_by_loc(iloc)%end = iloc + 1
enddo

call ufo_geovals_update_linit(self)

end subroutine ufo_geovals_setup_trivial_sampling_method

! ------------------------------------------------------------------------------

subroutine ufo_geovals_delete(self)
implicit none
type(ufo_geovals), intent(inout) :: self

integer :: ivar

if (allocated(self%geovals)) then
  do ivar = 1, self%nvar
    if (allocated(self%geovals(ivar)%vals)) deallocate(self%geovals(ivar)%vals)
  enddo
  deallocate(self%geovals)
endif
if (allocated(self%variables)) deallocate(self%variables)
if (allocated(self%reduced_variables)) deallocate(self%reduced_variables)
if (allocated(self%sampling_methods)) deallocate(self%sampling_methods)
if (allocated(self%sampling_method_by_var)) deallocate(self%sampling_method_by_var)
self%nvar = 0
self%nlocs = 0
self%nsampling_methods = 0
self%linit = .false.

end subroutine ufo_geovals_delete

! ------------------------------------------------------------------------------

subroutine ufo_geovals_add_reduced_vars(self, vars, nvars, nvals)
use oops_variables_mod
implicit none
type(ufo_geovals), intent(inout) :: self
type(oops_variables), intent(in) :: vars
integer(c_int), intent(in)       :: nvars
integer(c_size_t), intent(in)    :: nvals(nvars)

character(max_string) :: err_msg
integer :: ivar, inew_var, new_nvars, igeoval, imethod
integer :: sampling_method
integer, allocatable :: geoval_indices(:)  ! Indices of geovals storing the newly added variables

if (nvars /= vars%nvars()) then
  write(err_msg,*) "ufo_geovals_add_reduced_vars: mismatch in the number of variables ", &
    "(nvars = ", nvars, ", vars%nvars() = ", vars%nvars(), ")"
  call abor1_ftn(err_msg)
endif

if (nvars == 0) return  ! Nothing to do

new_nvars = self%nvar
allocate(geoval_indices(vars%nvars()))
do inew_var = 1, vars%nvars()
  ivar = ufo_vars_getindex(self%reduced_variables, vars%variable(inew_var))
  if (ivar >= 1) then
    write(err_msg,*) "ufo_geovals_add_reduced_vars: variable ", vars%variable(inew_var), &
      " already exists in reduced format"
    call abor1_ftn(err_msg)
  endif
  ivar = ufo_vars_getindex(self%variables, vars%variable(inew_var))
  if (ivar >= 1) then
    if (.not. self%sampling_methods(self%sampling_method_by_var(ivar))%trivial .or. &
        nvals(inew_var) /= self%geovals(ivar)%nval) then
      ! The reduced version of this variable is different from the sampled version,
      ! so it needs to be stored separately
      new_nvars = new_nvars + 1
      geoval_indices(inew_var) = new_nvars
    else
      ! The reduced and sampled versions of this variable are identical
      geoval_indices(inew_var) = ivar
    endif
  else
    ! Only the reduced version of this variable will be stored; the sampled version won't
    new_nvars = new_nvars + 1
    geoval_indices(inew_var) = new_nvars
  endif
enddo

sampling_method = 0
do imethod = 1, self%nsampling_methods
  if (self%sampling_methods(imethod)%trivial) then
    sampling_method = imethod
    exit
  endif
enddo
if (sampling_method == 0) then
  ! Reduced variables need to be mapped to a trivial sampling method, but none exists.
  ! We need to create one on our own.
  call resize(self%sampling_methods, self%nsampling_methods + 1)
  self%nsampling_methods = self%nsampling_methods + 1
  sampling_method = self%nsampling_methods
  self%sampling_methods(sampling_method)%trivial = .true.
  self%sampling_methods(sampling_method)%npaths = self%nlocs
  call ufo_geovals_setup_trivial_sampling_method(self, sampling_method)
endif

call resize(self%geovals, new_nvars)
call resize(self%variables, new_nvars)
call resize(self%reduced_variables, new_nvars)
call resize(self%sampling_method_by_var, new_nvars)

self%variables(self%nvar + 1 :) = ""

do inew_var = 1, vars%nvars()
  igeoval = geoval_indices(inew_var)
  self%reduced_variables(igeoval) = vars%variable(inew_var)
  if (igeoval > self%nvar) then
    call ufo_geovals_setup_geoval(self, igeoval, int(sampling_method, kind=c_size_t), &
                                  nvals(inew_var), int(self%nlocs, kind=c_size_t))
  endif
enddo

self%nvar = new_nvars

end subroutine ufo_geovals_add_reduced_vars

! ------------------------------------------------------------------------------

!> Set `vars` to the list of variables stored in the format `format`.
subroutine ufo_geovals_get_vars(self, vars, format)
use oops_variables_mod
implicit none
type(ufo_geovals), intent(in), target :: self
type(oops_variables), intent(inout) :: vars
integer(c_int), intent(in), optional :: format

integer(c_int) :: actual_format
character(len=MAXVARLEN), pointer :: variables(:)  !< either the sampled or reduced variable list
integer :: i

if (present(format)) then
  actual_format = format
else
  actual_format = self%default_format
endif

if (actual_format == ufo_geoval_sampled) then
  variables => self%variables
else if (actual_format == ufo_geoval_reduced) then
  variables => self%reduced_variables
else
  call abor1_ftn('ufo_geovals_get_vars: format must be either ufo_geoval_sampled or ufo_geoval_reduced')
endif

call vars%clear()
do i = 1, size(variables)
  if (variables(i) /= "") call vars%push_back(variables(i))
enddo

end subroutine ufo_geovals_get_vars

! ------------------------------------------------------------------------------

integer(c_int) function ufo_geovals_get_default_format(self)
implicit none
type(ufo_geovals), intent(in) :: self

ufo_geovals_get_default_format = self%default_format

end function ufo_geovals_get_default_format

! ------------------------------------------------------------------------------

subroutine ufo_geovals_set_default_format(self, format)
implicit none
type(ufo_geovals), intent(inout) :: self
integer(c_int), intent(in) :: format

if (format /= ufo_geoval_sampled .and. format /= ufo_geoval_reduced) then
  call abor1_ftn('ufo_geovals_set_default_format: format must be either ufo_geoval_sampled or ufo_geoval_reduced')
endif

self%default_format = format

end subroutine ufo_geovals_set_default_format

! ------------------------------------------------------------------------------

logical(c_bool) function ufo_geovals_are_reduced_and_sampled_formats_aliased(self, varname)
implicit none
type(ufo_geovals), intent(in) :: self
character(len=*), intent(in) :: varname

integer :: ivar, ireduced_var

ivar = ufo_vars_getindex(self%variables, varname)
ireduced_var = ufo_vars_getindex(self%reduced_variables, varname)

ufo_geovals_are_reduced_and_sampled_formats_aliased = ivar > 0 .and. ivar == ireduced_var

end function ufo_geovals_are_reduced_and_sampled_formats_aliased

! ------------------------------------------------------------------------------

subroutine ufo_geovals_get_var(self, varname, geoval, format, must_be_found)
implicit none
type(ufo_geovals), target, intent(in)    :: self
character(len=*), intent(in) :: varname
type(ufo_geoval), pointer, intent(inout)    :: geoval
integer(c_int), intent(in), optional :: format
!> By default, the program will abort if the GeoVaL `varname` is not found. Set `must_be_found` to
!> `.false.` if in these circumstances the subroutine should set `geoval` to a null pointer instead.
logical, intent(in), optional :: must_be_found

character(len=*), parameter :: myname_="ufo_geovals_get_var"

character(max_string) :: err_msg
integer(c_int) :: actual_format
logical :: actual_must_be_found
character(len=MAXVARLEN), pointer :: variables(:)  !< either the sampled or reduced variable list
integer :: ivar, jv

geoval => NULL()

if (present(format)) then
  actual_format = format
else
  actual_format = self%default_format
endif

if (present(must_be_found)) then
  actual_must_be_found = must_be_found
else
  actual_must_be_found = .true.
endif

if (actual_format == ufo_geoval_sampled) then
  variables => self%variables
else if (actual_format == ufo_geoval_reduced) then
  variables => self%reduced_variables
else
  write(err_msg,*) myname_, ' format must be either ufo_geoval_sampled or ufo_geoval_reduced'
  call abor1_ftn(err_msg)
endif

ivar = ufo_vars_getindex(variables, varname)

if (ivar < 0) then
  if (actual_must_be_found) then
    write(0,*)'ufo_geovals_get_var looking for '
    if (actual_format == ufo_geoval_sampled) then
      write(0,*)'sampled '
    else
      write(0,*)'reduced '
    endif
    write(0,*)'variable ',trim(varname),' in:'
    do jv=1,self%nvar
      write(0,*)'ufo_geovals_get_var ',jv,trim(variables(jv))
    enddo
    write(err_msg,*) myname_, " ", trim(varname), ' doesnt exist'
    call abor1_ftn(err_msg)
  endif  ! if actual_must_be_found is .false., we leave geoval as a null pointer
else
  geoval => self%geovals(ivar)
endif

end subroutine ufo_geovals_get_var

! ------------------------------------------------------------------------------

subroutine ufo_geovals_get_var_by_index(self, ivar, geoval)
implicit none
type(ufo_geovals), target, intent(in)    :: self
integer, intent(in) :: ivar
type(ufo_geoval), pointer, intent(inout)    :: geoval

character(len=*), parameter :: myname_="ufo_geovals_get_var_by_index"

character(max_string) :: err_msg
integer :: jv

geoval => NULL()

if (ivar <= 0 .or. ivar > self % nvar) then
  write(err_msg,*) myname_, " index ", ivar, ' doesnt exist'
  call abor1_ftn(err_msg)
else
  geoval => self%geovals(ivar)
endif

end subroutine ufo_geovals_get_var_by_index

! ------------------------------------------------------------------------------

subroutine ufo_geovals_zero(self)
implicit none
type(ufo_geovals), intent(inout) :: self
integer :: ivar

if (.not. self%linit) then
  call abor1_ftn("ufo_geovals_zero: geovals not initialized")
endif
do ivar = 1, self%nvar
  self%geovals(ivar)%vals(:,:) = 0.0
enddo

end subroutine ufo_geovals_zero

! ------------------------------------------------------------------------------

subroutine ufo_geovals_abs(self)
implicit none
type(ufo_geovals), intent(inout) :: self
integer :: ivar

if (.not. self%linit) then
  call abor1_ftn("ufo_geovals_abs: geovals not initialized")
endif
do ivar = 1, self%nvar
  self%geovals(ivar)%vals = abs(self%geovals(ivar)%vals)
enddo

end subroutine ufo_geovals_abs

! ------------------------------------------------------------------------------

subroutine ufo_geovals_rms(self,vrms)
implicit none
type(ufo_geovals), intent(in) :: self
real(kind_real), intent(inout) :: vrms
integer :: jv, jo
real(kind_real) :: N

if (.not. self%linit) then
  call abor1_ftn("ufo_geovals_rms: geovals not initialized")
endif
vrms=0.0_kind_real
N=0.0_kind_real
do jv = 1, self%nvar
   do jo = 1, self%geovals(jv)%nprofiles
      vrms = vrms + Sum(self%geovals(jv)%vals(:,jo)**2)
      N=N+self%geovals(jv)%nval
   enddo
enddo

if ( N > 0) vrms = sqrt(vrms/N)

end subroutine ufo_geovals_rms

! ------------------------------------------------------------------------------

subroutine ufo_geovals_random(self)
use random_mod
implicit none
type(ufo_geovals), intent(inout) :: self
integer :: ivar
integer :: rseed = 7

if (.not. self%linit) then
  call abor1_ftn("ufo_geovals_random: geovals not initialized")
endif
do ivar = 1, self%nvar
  call normal_distribution(self%geovals(ivar)%vals, 0.0_kind_real, 1.0_kind_real, rseed)
enddo

end subroutine ufo_geovals_random

! ------------------------------------------------------------------------------

subroutine ufo_geovals_scalmult(self, zz)
implicit none
type(ufo_geovals), intent(inout) :: self
real(kind_real), intent(in) :: zz
integer :: jv, jo, jz

if (.not. self%linit) then
  call abor1_ftn("ufo_geovals_scalmult: geovals not allocated")
endif

do jv=1,self%nvar
  do jo=1,self%geovals(jv)%nprofiles
    do jz = 1, self%geovals(jv)%nval
      self%geovals(jv)%vals(jz,jo) = zz * self%geovals(jv)%vals(jz,jo)
    enddo
  enddo
enddo

end subroutine ufo_geovals_scalmult

! ------------------------------------------------------------------------------

subroutine ufo_geovals_profmult(self, nlocs, values)
implicit none
type(ufo_geovals), intent(inout) :: self
integer(c_int), intent(in) :: nlocs
real(c_float), intent(in) :: values(nlocs)
integer :: jv, jp, jl, jm

if (.not. self%linit) then
  call abor1_ftn("ufo_geovals_profmult: geovals not allocated")
endif

if (nlocs /= self%nlocs) then
  call abor1_ftn("ufo_geovals_profmult: nlocs mismatch")
endif

do jv = 1, self%nvar
  jm = self%sampling_method_by_var(jv)
  do jl = 1, self%nlocs
    do jp = self%sampling_methods(jm)%paths_by_loc(jl)%begin, &
            self%sampling_methods(jm)%paths_by_loc(jl)%end - 1
      self%geovals(jv)%vals(:,jp) = values(jl) * self%geovals(jv)%vals(:,jp)
    enddo
  enddo
enddo

end subroutine ufo_geovals_profmult

! ------------------------------------------------------------------------------


subroutine ufo_geovals_assign(self, rhs)
implicit none
type(ufo_geovals), intent(inout), target :: self
type(ufo_geovals), intent(in), target :: rhs
integer :: jv, jo, jz
integer :: iv
character(max_string) :: err_msg
character(len=MAXVARLEN), pointer :: self_variables(:), rhs_variables(:)

if (.not. self%linit) then
  call abor1_ftn("ufo_geovals_assign: geovals not allocated")
endif
if (.not. rhs%linit) then
  call abor1_ftn("ufo_geovals_assign: geovals not allocated")
endif

if (self%nlocs /= rhs%nlocs) then
  call abor1_ftn("ufo_geovals_assign: nlocs different between lhs and rhs")
endif

do jv=1,self%nvar
  ! Determine if this is an sampled or reduced variable
  if (self%variables(jv) /= "") then
    self_variables => self%variables
    rhs_variables => rhs%variables
  else if (self%reduced_variables(jv) /= "") then
    self_variables => self%reduced_variables
    rhs_variables => rhs%reduced_variables
  else
    continue
  endif

  iv = ufo_vars_getindex(rhs_variables, self_variables(jv))
  if (iv < 0) then
    write(err_msg,*) 'ufo_geovals_assign: var ', trim(self_variables(jv)), ' doesnt exist in rhs'
    call abor1_ftn(trim(err_msg))
  endif
  if (self%geovals(jv)%nval /= rhs%geovals(iv)%nval) then
    write(err_msg,*) 'ufo_geovals_assign: nvals for var ', trim(self_variables(jv)), ' are different in lhs and rhs'
    call abor1_ftn(trim(err_msg))
  endif
  if (self%geovals(jv)%nprofiles /= rhs%geovals(iv)%nprofiles) then
    write(err_msg,*) 'ufo_geovals_assign: nprofiles for var ', trim(self_variables(jv)), ' are different in lhs and rhs'
    call abor1_ftn(trim(err_msg))
  endif
  do jo=1,self%geovals(jv)%nprofiles
    do jz = 1, self%geovals(jv)%nval
      self%geovals(jv)%vals(jz,jo) = rhs%geovals(iv)%vals(jz,jo)
    enddo
  enddo
enddo

end subroutine ufo_geovals_assign

! ------------------------------------------------------------------------------
subroutine ufo_geovals_reorderzdir(self, varname, zdir)
implicit none
type(ufo_geovals),intent(inout) :: self
character(len=*), intent(in) :: varname
character(len=*), intent(in) :: zdir

type(ufo_geovals) :: selfclone
type(ufo_geoval), pointer :: geoval
character(max_string) :: err_msg
integer:: iobs, ivar, ival, kval
logical :: do_flip                   !< .true. if all the ufo_geoval arrays inside geovals

do_flip = .false.

if (.not. self%linit) then
  call abor1_ftn("ufo_geovals_reorderzdir: geovals not allocated")
endif

if (trim(zdir) /= "bottom2top" .and. trim(zdir) /= "top2bottom") then
  write(err_msg, *) 'ufo_geovals_reorderzdir: z-coordinate direction ', trim(zdir), ' not defined. ', &
                    'use either bottom2top or top2bottom'
  call abor1_ftn(err_msg)
end if

! Get vertical coordinate variable
call ufo_geovals_get_var(self, varname, geoval)
if (.not. associated(geoval)) then
  write(err_msg, *) 'ufo_geovals_reorderzdir: geoval vertical coordinate variable ', trim(varname), ' doesnt exist'
  call abor1_ftn(err_msg)
endif

! Check if reorder variables is necessary based on the direction defined by zdir
if ((trim(zdir) == "bottom2top" .and. geoval%vals(1,1) < geoval%vals(geoval%nval,1)) .or. &
    (trim(zdir) == "top2bottom" .and. geoval%vals(1,1) > geoval%vals(geoval%nval,1))) then
   do_flip = .true.
else
   return
endif

call ufo_geovals_copy(self, selfclone)

if (do_flip) then
  do ivar = 1, self%nvar
    do ival = 1, self%geovals(ivar)%nval
      kval = self%geovals(ivar)%nval - ival + 1
      self%geovals(ivar)%vals(ival,:) = selfclone%geovals(ivar)%vals(kval,:)
    enddo
  enddo
endif

call ufo_geovals_delete(selfclone)

end subroutine ufo_geovals_reorderzdir

! ------------------------------------------------------------------------------
!> Sum of two GeoVaLs objects

subroutine ufo_geovals_add(self, other)
implicit none
type(ufo_geovals), intent(inout), target :: self
type(ufo_geovals), intent(in), target :: other
integer :: jv, jo, jz
integer :: iv
character(max_string) :: err_msg
character(len=MAXVARLEN), pointer :: self_variables(:), other_variables(:)

if (.not. self%linit) then
  call abor1_ftn("ufo_geovals_add: geovals not allocated")
endif
if (.not. other%linit) then
  call abor1_ftn("ufo_geovals_add: geovals not allocated")
endif

do jv=1,self%nvar
  ! Determine if this is an sampled or reduced variable
  if (self%variables(jv) /= "") then
    self_variables => self%variables
    other_variables => other%variables
  else if (self%reduced_variables(jv) /= "") then
    self_variables => self%reduced_variables
    other_variables => other%reduced_variables
  else
    continue
  endif

  iv = ufo_vars_getindex(other_variables, self_variables(jv))
  if (iv .ne. -1) then !Only add if exists in RHS
    if (self%geovals(jv)%nval /= other%geovals(iv)%nval) then
      write(err_msg,*) 'ufo_geovals_add: nvals for var ', trim(self_variables(jv)), ' are different in lhs and rhs'
      call abor1_ftn(trim(err_msg))
    endif
    if (self%geovals(jv)%nprofiles /= other%geovals(iv)%nprofiles) then
      write(err_msg,*) 'ufo_geovals_add: nprofiles for var ', trim(self_variables(jv)), ' are different in lhs and rhs'
      call abor1_ftn(trim(err_msg))
    endif
    do jo=1,self%geovals(jv)%nprofiles
      do jz = 1, self%geovals(jv)%nval
        self%geovals(jv)%vals(jz,jo) = self%geovals(jv)%vals(jz,jo) + other%geovals(iv)%vals(jz,jo)
      enddo
    enddo
  endif
enddo

end subroutine ufo_geovals_add

! ------------------------------------------------------------------------------
!> Difference between two GeoVaLs objects

subroutine ufo_geovals_diff(self, other)
implicit none
type(ufo_geovals), intent(inout), target :: self
type(ufo_geovals), intent(in), target :: other
integer :: jv, jo, jz
integer :: iv
character(max_string) :: err_msg
character(len=MAXVARLEN), pointer :: self_variables(:), other_variables(:)

if (.not. self%linit) then
  call abor1_ftn("ufo_geovals_diff: geovals not allocated")
endif
if (.not. other%linit) then
  call abor1_ftn("ufo_geovals_diff: geovals not allocated")
endif

do jv=1,self%nvar
  ! Determine if this is an sampled or reduced variable
  if (self%variables(jv) /= "") then
    self_variables => self%variables
    other_variables => other%variables
  else if (self%reduced_variables(jv) /= "") then
    self_variables => self%reduced_variables
    other_variables => other%reduced_variables
  else
    continue
  endif

  iv = ufo_vars_getindex(other_variables, self_variables(jv))
  if (iv .ne. -1) then !Only subtract if exists in RHS
    if (self%geovals(jv)%nval /= other%geovals(iv)%nval) then
      write(err_msg,*) 'ufo_geovals_diff: nvals for var ', trim(self_variables(jv)), ' are different in lhs and rhs'
      call abor1_ftn(trim(err_msg))
    endif
    if (self%geovals(jv)%nprofiles /= other%geovals(iv)%nprofiles) then
      write(err_msg,*) 'ufo_geovals_diff: nprofiles for var ', trim(self_variables(jv)), ' are different in lhs and rhs'
      call abor1_ftn(trim(err_msg))
    endif
    do jo=1,self%geovals(jv)%nprofiles
      do jz = 1, self%geovals(jv)%nval
        self%geovals(jv)%vals(jz,jo) = self%geovals(jv)%vals(jz,jo) - other%geovals(iv)%vals(jz,jo)
      enddo
    enddo
  endif
enddo

end subroutine ufo_geovals_diff

! ------------------------------------------------------------------------------
!> Schur product of two GeoVaLs objects

subroutine ufo_geovals_schurmult(self, other)
implicit none
type(ufo_geovals), intent(inout), target :: self
type(ufo_geovals), intent(in), target :: other
integer :: jv, jo, jz
integer :: iv
character(max_string) :: err_msg
character(len=MAXVARLEN), pointer :: self_variables(:), other_variables(:)

if (.not. self%linit) then
  call abor1_ftn("ufo_geovals_schurmult: geovals not allocated")
endif
if (.not. other%linit) then
  call abor1_ftn("ufo_geovals_schurmult: geovals not allocated")
endif

do jv=1,self%nvar
  ! Determine if this is an sampled or reduced variable
  if (self%variables(jv) /= "") then
    self_variables => self%variables
    other_variables => other%variables
  else if (self%reduced_variables(jv) /= "") then
    self_variables => self%reduced_variables
    other_variables => other%reduced_variables
  else
    continue
  endif

  iv = ufo_vars_getindex(other_variables, self_variables(jv))
  if (iv .ne. -1) then !Only mult if exists in RHS
    if (self%geovals(jv)%nval /= other%geovals(iv)%nval) then
      write(err_msg,*) 'ufo_geovals_schurmult: nvals for var ', trim(self_variables(jv)), ' are different in lhs and rhs'
      call abor1_ftn(trim(err_msg))
    endif
    if (self%geovals(jv)%nprofiles /= other%geovals(iv)%nprofiles) then
      write(err_msg,*) 'ufo_geovals_schurmult: nprofiles for var ', trim(self_variables(jv)), ' are different in lhs and rhs'
      call abor1_ftn(trim(err_msg))
    endif
    do jo=1,self%geovals(jv)%nprofiles
      do jz = 1, self%geovals(jv)%nval
        self%geovals(jv)%vals(jz,jo) = self%geovals(jv)%vals(jz,jo) * other%geovals(iv)%vals(jz,jo)
      enddo
    enddo
  endif
enddo

end subroutine ufo_geovals_schurmult

! ------------------------------------------------------------------------------
!> Copy one GeoVaLs object into another
!!

subroutine ufo_geovals_copy(self, other)
implicit none
type(ufo_geovals), intent(in) :: self
type(ufo_geovals), intent(inout) :: other
integer :: jv, jc

if (.not. self%linit) then
  call abor1_ftn("ufo_geovals_copy: geovals not defined")
endif

call ufo_geovals_delete(other)

other%nlocs = self%nlocs
other%nvar = self%nvar
other%default_format = self%default_format
allocate(other%variables(other%nvar))
other%variables(:) = self%variables(:)
allocate(other%reduced_variables(other%nvar))
other%reduced_variables(:) = self%reduced_variables(:)

other%nsampling_methods = self%nsampling_methods
allocate(other%sampling_methods(size(self%sampling_methods)))
do jc = 1, size(other%sampling_methods)
  other%sampling_methods(jc)%npaths = self%sampling_methods(jc)%npaths
  allocate(other%sampling_methods(jc)%paths_by_loc(size(self%sampling_methods(jc)%paths_by_loc)))
  other%sampling_methods(jc)%paths_by_loc(:) = self%sampling_methods(jc)%paths_by_loc(:)
enddo
allocate(other%sampling_method_by_var(size(self%sampling_method_by_var)))
other%sampling_method_by_var(:) = self%sampling_method_by_var(:)

allocate(other%geovals(other%nvar))
do jv = 1, other%nvar
  other%geovals(jv)%nval = self%geovals(jv)%nval
  other%geovals(jv)%nprofiles = self%geovals(jv)%nprofiles
  allocate(other%geovals(jv)%vals(other%geovals(jv)%nval, other%geovals(jv)%nprofiles))
  other%geovals(jv)%vals(:,:) = self%geovals(jv)%vals(:,:)
enddo

other%missing_value = self%missing_value
other%linit = .true.

end subroutine ufo_geovals_copy

! ------------------------------------------------------------------------------
!> Copy one location from GeoVaLs into a new object
!!

subroutine ufo_geovals_copy_one(self, other, loc_index)
implicit none
type(ufo_geovals), intent(inout) :: self !> GeoVaLs for one location
type(ufo_geovals), intent(in) :: other   !> GeoVaLs for many location
integer, intent(in) :: loc_index !> Index of the location in the "other" geoval

integer :: jv, jc
type(ufo_index_range) :: path_index_range

if (.not. other%linit) then
  call abor1_ftn("ufo_geovals_copy_one: geovals not defined")
endif

call ufo_geovals_delete(self)

self%nlocs = 1
self%nvar = other%nvar
self%default_format = other%default_format
allocate(self%variables(self%nvar))
self%variables(:) = other%variables(:)
allocate(self%reduced_variables(self%nvar))
self%reduced_variables(:) = other%reduced_variables(:)

self%nsampling_methods = other%nsampling_methods
allocate(self%sampling_methods(size(other%sampling_methods)))
do jc = 1, size(self%sampling_methods)
  path_index_range = other%sampling_methods(jc)%paths_by_loc(loc_index)
  self%sampling_methods(jc)%npaths = path_index_range%end - path_index_range%begin
  allocate(self%sampling_methods(jc)%paths_by_loc(1))
  self%sampling_methods(jc)%paths_by_loc(1)%begin = 1
  self%sampling_methods(jc)%paths_by_loc(1)%end = 1 + self%sampling_methods(jc)%npaths
enddo
allocate(self%sampling_method_by_var(size(other%sampling_method_by_var)))
self%sampling_method_by_var(:) = other%sampling_method_by_var(:)

allocate(self%geovals(self%nvar))
do jv = 1, self%nvar
  jc = other%sampling_method_by_var(jv)
  self%geovals(jv)%nval = other%geovals(jv)%nval
  self%geovals(jv)%nprofiles = self%sampling_methods(jc)%npaths
  allocate(self%geovals(jv)%vals(self%geovals(jv)%nval, self%geovals(jv)%nprofiles))
  path_index_range = other%sampling_methods(jc)%paths_by_loc(loc_index)
  self%geovals(jv)%vals(:,:) = &
    other%geovals(jv)%vals(:, path_index_range%begin : path_index_range%end-1)
enddo

self%missing_value = other%missing_value
self%linit = .true.

end subroutine ufo_geovals_copy_one

! ------------------------------------------------------------------------------
!> Initialize a GeoVaLs object based on an analytic state
!!
!! \details **ufo_geovals_analytic_init_c()** takes an existing ufo::GeoVaLs object
!! and fills in values based on one of several analytic solutions.  This initialization
!! is intended to be used with the **TestStateInterpolation()** test; see there for
!! further information.
!!
!! Currently implemented options for analytic_init include:
!! * dcmip-test-1-1: 3D deformational flow
!! * dcmip-test-1-2: 3D Hadley-like meridional circulation
!! * dcmip-test-3-1: Non-orographic gravity waves on a small planet
!! * dcmip-test-4-0: Baroclinic instability
!!
!! \warning Currently only temperature is implemented.  For variables other than
!! temperature, the input GeoVaLs object is not changed.  This effectively
!! disables the interpolation test for that variable by setting the normalized
!! error to zero.
!!
!! \warning Currently there is no conversion between temperature and virtual
!! temperature
!!
!! \date May, 2018: Created by M. Miesch (JCSDA)
!! \date June, 2018: Added dcmip-test-4.0 (M. Miesch, JCSDA)
!!
!! \sa test::TestStateInterpolation()
!!

subroutine ufo_geovals_analytic_init(self, sampled_locations, ic)
use dcmip_initial_conditions_test_1_2_3, only : test1_advection_deformation, &
                                  test1_advection_hadley, test3_gravity_wave
use dcmip_initial_conditions_test_4, only : test4_baroclinic_wave
use ufo_sampled_locations_mod
use ufo_utils_mod, only: cmp_strings

implicit none
type(ufo_geovals), intent(inout) :: self
type(ufo_sampled_locations), intent(in)  :: sampled_locations
character(*), intent(in)         :: ic

real(kind_real) :: pi = acos(-1.0_kind_real)
real(kind_real) :: deg_to_rad,rlat, rlon
real(kind_real) :: p0, kz, u0, v0, w0, t0, phis0, ps0, rho0, hum0
real(kind_real) :: q1, q2, q3, q4
real(kind_real), allocatable, dimension(:) :: lons, lats
integer :: npaths, ivar, iprofile, ival

if (.not. self%linit) then
  call abor1_ftn("ufo_geovals_analytic_init: geovals not defined")
endif

! The last variable should be the ln pressure coordinate.  That's
! where we get the height information for the analytic init
if (self%variables(self%nvar) /= var_prs .and. &
    self%variables(self%nvar) /= var_prsi) then
  call abor1_ftn("ufo_geovals_analytic_init: pressure coordinate not defined")
endif

deg_to_rad = pi/180.0_kind_real

npaths = sampled_locations%npaths()
allocate(lons(npaths), lats(npaths))
call sampled_locations%get_lons(lons)
call sampled_locations%get_lats(lats)

do ivar = 1, self%nvar-1

   do iprofile = 1, self%geovals(ivar)%nprofiles

      ! convert lat and lon to radians
      rlat = deg_to_rad * lats(iprofile)
      rlon = deg_to_rad*modulo(lons(iprofile)+180.0_kind_real,360.0_kind_real) - pi

      do ival = 1, self%geovals(ivar)%nval

         ! obtain height from the existing GeoVaLs object, which should be an
         ! output of the State::getValues() method
         ! should be delivered in units of Pa
         p0 = self%geovals(self%nvar)%vals(ival,iprofile)

         init_option: select case (trim(ic))

         case ("invent_state")

           t0 = cos(deg_to_rad * lons(iprofile) ) * cos(rlat)

         case ("dcmip-test-1-1")

            call test1_advection_deformation(rlon,rlat,p0,kz,0,u0,v0,w0,&
                                             t0,phis0,ps0,rho0,hum0,q1,q2,q3,q4)

         case ("dcmip-test-1-2")

            call test1_advection_hadley(rlon,rlat,p0,kz,0,u0,v0,w0,&
                                        t0,phis0,ps0,rho0,hum0,q1)

         case ("dcmip-test-3-1")

            call test3_gravity_wave(rlon,rlat,p0,kz,0,u0,v0,w0,&
                                        t0,phis0,ps0,rho0,hum0)

         case ("dcmip-test-4-0")

            call test4_baroclinic_wave(0,1.0_kind_real,rlon,rlat,p0,kz,0,u0,v0,w0,&
                                        t0,phis0,ps0,rho0,hum0,q1,q2)

         case default

            call abor1_ftn("ufo_geovals_analytic_init: invalid analytic_init")

         end select init_option

         ! currently only temperature is implemented
         if (cmp_strings(self%variables(ivar), var_tv)) then
            ! Warning: we may need a conversion from temperature to
            ! virtual temperture here
            self%geovals(ivar)%vals(ival,iprofile) = t0
         endif

      enddo
   enddo
enddo

deallocate(lons, lats)

end subroutine ufo_geovals_analytic_init

! ------------------------------------------------------------------------------
!> Normalization of one GeoVaLs object by another
!!
!! \details This is a normalization operator that first computes the normalization
!! factor for each variable based on the rms amplitude of that variable across
!! all locations in the reference GeoVaLs object (other).  Then each element of
!! the input GeoVals object (self) is divided by these normalization factors.
!! The operation is done in place.  So, after execution, the input GeoVaLs
!! object will be nondimensional.
!!
!! \warning If the reference variable is identially zero across all
!! locations, then the result of this operatution is set to zero for that
!! variable.  This is to used to bypass variables that do not have a reference
!! value in the State interpolation test.
!!

subroutine ufo_geovals_normalize(self, other)
implicit none
type(ufo_geovals), intent(inout) :: self !> Input GeoVaLs object (LHS)
type(ufo_geovals), intent(in) :: other   !> Reference GeoVaLs object (RHS)
integer :: jv, jo, jz
real(kind_real) :: over_nloc, vrms, norm

if (.not. self%linit) then
  call abor1_ftn("ufo_geovals_normalize: geovals not allocated")
endif
if (.not. other%linit) then
  call abor1_ftn("ufo_geovals_normalize: geovals not allocated")
endif
if (self%nvar /= other%nvar) then
  call abor1_ftn("ufo_geovals_normalize: reference geovals object must have the same variables as the sampled")
endif


do jv=1,self%nvar

   !> Compute normalization factors for the errors based on the rms amplitude of
   !! each variable across all of the selected locations.  Use the "other" GeoVaLs
   !! object as a reference, since this may be the exact analytic answer

   over_nloc = 1.0_kind_real / &
        (real(other%geovals(jv)%nprofiles,kind_real)*real(other%geovals(jv)%nval,kind_real))

   vrms = 0.0_kind_real
   do jo = 1, other%geovals(jv)%nprofiles
      do jz = 1, other%geovals(jv)%nval
         vrms = vrms + other%geovals(jv)%vals(jz,jo)**2
      enddo
   enddo

   if (vrms > 0.0_kind_real) then
      norm = 1.0_kind_real / sqrt(vrms*over_nloc)
   else
      norm = 0.0_kind_real
   endif

   ! Now loop through the LHS locations to compute the normalized value
   do jo=1,self%geovals(jv)%nprofiles
      do jz = 1, self%geovals(jv)%nval
         self%geovals(jv)%vals(jz,jo) = norm*self%geovals(jv)%vals(jz,jo)
      enddo
   enddo
enddo

end subroutine ufo_geovals_normalize

! ------------------------------------------------------------------------------

subroutine ufo_geovals_minmaxavg(self, kobs, kvar, pmin, pmax, prms)
implicit none
integer, intent(inout) :: kobs
integer, intent(in) :: kvar
real(kind_real), intent(inout) :: pmin, pmax, prms
type(ufo_geovals), intent(in) :: self
integer :: jo, jz, jv

if (.not. self%linit) then
  call abor1_ftn("ufo_geovals_minmaxavg: geovals not initialized")
endif

jv = kvar+1
kobs = 0
pmin = huge(pmin)
pmax = -huge(pmax)
prms = 0.0_kind_real
do jo = 1, self%geovals(jv)%nprofiles
  do jz = 1, self%geovals(jv)%nval
    if (self%geovals(jv)%vals(jz,jo) .ne. self%missing_value) then
      kobs = kobs + 1
      if (self%geovals(jv)%vals(jz,jo) < pmin) pmin = self%geovals(jv)%vals(jz,jo)
      if (self%geovals(jv)%vals(jz,jo) > pmax) pmax = self%geovals(jv)%vals(jz,jo)
      prms = prms + self%geovals(jv)%vals(jz,jo) * self%geovals(jv)%vals(jz,jo)
    endif
  enddo
enddo
if (kobs > 0) prms = sqrt(prms/real(kobs,kind_real))

end subroutine ufo_geovals_minmaxavg

! ------------------------------------------------------------------------------
!> Location where the summed geovals value is maximum
!!
!! \details This routine computes the rms value over the vertical profile for
!! each location and observation then returns the location number and the
!! variable number where this rms value is maximum.  Intended for use with
!! the State interpotation test in which the input GeoVaLs object is a
!! nondimensional, positive-definite error measurement.

subroutine ufo_geovals_maxloc(self, mxval, iobs, ivar)
implicit none
real(kind_real), intent(inout) :: mxval
integer, intent(inout) :: iobs, ivar

type(ufo_geovals), intent(in) :: self
real(kind_real) :: vrms
integer :: jv, jo, jz

if (.not. self%linit) then
  call abor1_ftn("ufo_geovals_maxloc: geovals not allocated")
endif

mxval = 0.0_kind_real
iobs = 1
ivar = 1

do jv = 1,self%nvar
   do jo = 1, self%geovals(jv)%nprofiles

      vrms = 0.0_kind_real
      do jz = 1, self%geovals(jv)%nval
         vrms = vrms + self%geovals(jv)%vals(jz,jo)**2
      enddo

      if ( self%geovals(jv)%nval > 0 ) then
        vrms = sqrt(vrms/real(self%geovals(jv)%nval,kind_real))
      end if

      if (vrms > mxval) then
         mxval = vrms
         iobs = jo
         ivar = jv
      endif

   enddo
enddo

end subroutine ufo_geovals_maxloc

! ------------------------------------------------------------------------------

subroutine ufo_geovals_read_netcdf(self, filename, loc_multiplier, levels_are_top_down, c_obspace, vars)
use netcdf
use oops_variables_mod
implicit none
type(ufo_geovals), intent(inout), target :: self
character(max_string), intent(in) :: filename
integer, intent(in)               :: loc_multiplier
logical, intent(in)               :: levels_are_top_down
type(c_ptr), intent(in)           :: c_obspace
type(oops_variables), intent(in)  :: vars

integer :: global_npaths, var_global_npaths
integer :: nval
integer :: obs_nlocs
integer :: obs_all_nlocs
integer :: my_npaths, my_loc, global_loc, my_path, global_path, global_path_start, global_path_end
character(len=MAXVARLEN), pointer :: varname
type(oops_variables)  :: reduced_vars

integer :: ncid, dimid, varid, vartype, ndims
integer, dimension(3) :: dimids
integer :: ivar
integer :: ierr

character(max_string) :: err_msg
character(len=30) :: obs_nlocs_str
character(len=30) :: geo_nlocs_str

integer(c_size_t), allocatable, dimension(:) :: global_loc_by_my_loc
integer(c_size_t), allocatable, dimension(:) :: global_path_by_my_path
type(ufo_index_range), allocatable :: path_ranges_by_loc(:)
real, allocatable :: field2d(:,:), field1d(:)

! At the time being, the NetCDF GeoVaLs file format requires all GeoVaLs to be interpolated
! along the same set of paths composed of as many paths as there are observation locations.
! The format will be made more general in future.
integer, parameter :: nsampling_methods = 1
integer(c_size_t), allocatable :: sampling_method_by_var(:)
logical(c_bool) :: is_sampling_method_trivial(1)

! open netcdf file
call check('nf90_open', nf90_open(trim(filename),nf90_nowrite,ncid))

! find how many locs are in the file
ierr = nf90_inq_dimid(ncid, "nlocs", dimid)
if(ierr /= nf90_noerr) then
  write(err_msg,*) "Error: Dimension nlocs not found in ", trim(filename)
  call abor1_ftn(err_msg)
endif

call check('nf90_inq_dimid', nf90_inq_dimid(ncid, "nlocs", dimid))
call check('nf90_inquire_dimension', nf90_inquire_dimension(ncid, dimid, len = global_npaths))

!> round-robin distribute the observations to PEs
!> Calculate how many obs. on each PE
obs_all_nlocs = obsspace_get_gnlocs(c_obspace)
obs_nlocs = obsspace_get_nlocs(c_obspace)
allocate(global_loc_by_my_loc(obs_nlocs))
call obsspace_get_index(c_obspace, global_loc_by_my_loc)

! loc_multiplier specifies how many locations in the geovals file per
! single location in the obs file. There needs to be at least
! loc_multiplier * obs_all_nlocs locations in the geovals file.

if (global_npaths .lt. (loc_multiplier * obs_all_nlocs)) then
  write(obs_nlocs_str, *) loc_multiplier * obs_all_nlocs
  write(geo_nlocs_str, *) global_npaths
  write(err_msg,'(7a)') &
     "Error: Number of locations in the geovals file (", &
     trim(adjustl(geo_nlocs_str)), ") must be greater than or equal to ", &
     "the product of loc_multiplier and number of locations in the ", &
     "obs file (", trim(adjustl(obs_nlocs_str)), ")"
  call abor1_ftn(err_msg)
endif

! We have enough locations in the geovals file to cover the span of the
! number of locations in the obs file. Generate the global_path_by_my_path map according
! to the loc_multiplier and obs_nlocs values.
allocate(path_ranges_by_loc(obs_nlocs))
if (loc_multiplier > 0) then
  my_npaths = loc_multiplier * obs_nlocs
  allocate(global_path_by_my_path(my_npaths))
  my_path = 1
  do my_loc = 1,obs_nlocs
    global_path_start = ((global_loc_by_my_loc(my_loc) - 1) * loc_multiplier) + 1
    global_path_end = global_loc_by_my_loc(my_loc) * loc_multiplier
    path_ranges_by_loc(my_loc)%begin = my_path
    do global_path = global_path_start, global_path_end
      global_path_by_my_path(my_path) = global_path
      my_path = my_path + 1
    enddo
    path_ranges_by_loc(my_loc)%end = my_path
  enddo
else

  ! Negative loc_multipliers are no longer supported. They used to map each location to a set
  ! of paths with non-consecutive indices, which cannot be represented with the current data
  ! structures.
  write(err_msg, *) "Error: loc_multiplier must be positive"
  call abor1_ftn(err_msg)
end if

! define a common sampling method for all variables
allocate(sampling_method_by_var(vars%nvars()))
sampling_method_by_var(:) = 1
is_sampling_method_trivial(1) = (loc_multiplier == 1)

! Currently NetCDF files store variables only in the sampled format. But if loc_multiplier is
! equal to 1, there is no difference between the sampled and reduced format, so we can set up a
! GeoVaLs object storing variables in both formats aliased with each other. Otherwise we set up an
! VaLs object storing variables only in the sampled format.
!
! In future we may want to extend the NetCDF file layout to accommodate both formats at once.
if (loc_multiplier == 1) then
  reduced_vars = oops_variables(vars)
else
  reduced_vars = oops_variables()
endif

! allocate geovals structure
call ufo_geovals_partial_setup(self, obs_nlocs, vars, vars%nvars(), &
                               nsampling_methods, sampling_method_by_var, &
                               reduced_vars, is_sampling_method_trivial)

call reduced_vars % destruct()

! specify which paths sample which locations
call ufo_geovals_setup_sampling_method(self, 1, my_npaths, obs_nlocs, path_ranges_by_loc)

deallocate(path_ranges_by_loc)
deallocate(sampling_method_by_var)

do ivar = 1, self%nvar

  if (self%variables(ivar) /= "") then
    varname => self%variables(ivar)
  else
    varname => self%reduced_variables(ivar)
  endif

  ierr = nf90_inq_varid(ncid, varname, varid)
  if(ierr /= nf90_noerr) then
    write(err_msg,*) "Error: Variable ", trim(varname), " not found in ", trim(filename)
    call abor1_ftn(err_msg)
  endif

  call check('nf90_inquire_variable', nf90_inquire_variable(ncid, varid, xtype = vartype, &
                                         ndims = ndims, dimids = dimids))
  !> read 1d variable

  if (ndims == 1) then
    call check('nf90_inquire_dimension', nf90_inquire_dimension(ncid, dimids(1), len = var_global_npaths))
    if (var_global_npaths /= global_npaths) then
      call abor1_ftn('ufo_geovals_read_netcdf: var dim /= global_npaths')
    endif

    nval = 1

    !> allocate geoval for this variable
    self%geovals(ivar)%nval = nval
    allocate(self%geovals(ivar)%vals(nval,my_npaths))

    !> fill the geoval
    allocate(field1d(var_global_npaths))

    call check('nf90_get_var', nf90_get_var(ncid, varid, field1d))

    self%geovals(ivar)%vals(1,:) = field1d(global_path_by_my_path)

    deallocate(field1d)

  !> read 2d variable
  elseif (ndims == 2) then
    call check('nf90_inquire_dimension', nf90_inquire_dimension(ncid, dimids(1), len = nval))
    call check('nf90_inquire_dimension', nf90_inquire_dimension(ncid, dimids(2), len = var_global_npaths))
    if (var_global_npaths /= global_npaths) then
      call abor1_ftn('ufo_geovals_read_netcdf: var dim /= global_npaths')
    endif

    !> allocate geoval for this variable
    self%geovals(ivar)%nval = nval
    allocate(self%geovals(ivar)%vals(nval,my_npaths))
    allocate(field2d(nval, var_global_npaths))

    call check('nf90_get_var', nf90_get_var(ncid, varid, field2d))
    if (.not. levels_are_top_down) then
      self%geovals(ivar)%vals(:,:) = field2d(nval:1:-1,global_path_by_my_path)
    else
      self%geovals(ivar)%vals(:,:) = field2d(:,global_path_by_my_path)
    endif
    deallocate(field2d)
  !> only 1d & 2d vars
  else

    call abor1_ftn('ufo_geovals_read_netcdf: can only read 1d and 2d fields')
  endif

  ! set the missing value equal to IODA missing_value
  where (self%geovals(ivar)%vals > 1.0e08) self%geovals(ivar)%vals = self%missing_value

enddo

call check('nf90_close', nf90_close(ncid))

call ufo_geovals_update_linit(self)

if (.not. self%linit) then
  write(err_msg,*) "ufo_geovals_read_netcdf: internal error: not all data structures have been properly initialized"
  call abor1_ftn(err_msg)
endif

end subroutine ufo_geovals_read_netcdf

! ------------------------------------------------------------------------------
subroutine ufo_geovals_write_netcdf(self, filename)
use netcdf
implicit none
type(ufo_geovals), intent(inout)  :: self
character(max_string), intent(in) :: filename

integer :: i
integer :: ncid, dimid_nlocs, dimid_nval, dims(2)
integer, allocatable :: ncid_var(:)

allocate(ncid_var(self%nvar))

call check('nf90_create', nf90_create(trim(filename),nf90_hdf5,ncid))
! TODO(wsmigaj): define a new format with nlocs replaced by nprofiles and stored
! separately for each variable
!------------------------------------------------------------------------------
!call check('nf90_def_dim', nf90_def_dim(ncid,'nlocs',self%nlocs, dimid_nlocs))
! This modification accounts for two scenarios: 
! the number of observations (nlocs) for the trivial sampling method, where self%geovals(1)%nprofiles 
! is equal to nlocs, thus not affecting previous cases; and the number of extended geovals required for non-trivial cases. 
! Further refinement is needed for storing each variable separately.
!------------------------------------------------------------------------------
call check('nf90_def_dim', nf90_def_dim(ncid,'nlocs',self%geovals(1)%nprofiles, dimid_nlocs))

dims(2) = dimid_nlocs

do i = 1, self%nvar
  call check('nf90_def_dim', &
       nf90_def_dim(ncid,trim(self%variables(i))//"_nval",self%geovals(i)%nval, dimid_nval))
  dims(1) = dimid_nval
  call check('nf90_def_var',  &
       nf90_def_var(ncid,trim(self%variables(i)),nf90_float,dims,ncid_var(i)))
enddo

call check('nf90_enddef', nf90_enddef(ncid))

do i = 1, self%nvar
  call check('nf90_put_var', nf90_put_var(ncid,ncid_var(i),self%geovals(i)%vals(:,:)))
enddo

call check('nf90_close', nf90_close(ncid))
deallocate(ncid_var)

end subroutine ufo_geovals_write_netcdf

! ------------------------------------------------------------------------------

subroutine ufo_geovals_fill(self, varname, c_nprofiles, c_indx, c_nlev, c_vals, levels_top_down)
implicit none
type(ufo_geovals), intent(inout) :: self
character(len=*), intent(in) :: varname
integer(c_int), intent(in) :: c_nprofiles
integer(c_int), intent(in) :: c_indx(c_nprofiles)
integer(c_int), intent(in) :: c_nlev
real(c_double), intent(in) :: c_vals(c_nprofiles, c_nlev)
logical(c_bool), intent(in) :: levels_top_down

type(ufo_geoval), pointer :: geoval
integer :: jlev, jprofile, ilev, iprofile
integer :: lbgn, linc

if (.not.self%linit) call abor1_ftn("ufo_geovals_fill: geovals not initialized")

call ufo_geovals_get_var(self, varname, geoval, ufo_geoval_sampled)

if (geoval%nval /= c_nlev) call abor1_ftn("ufo_geovals_fill: incorrect number of levels")

! setting loop indices to ensure geovals are filled top to bottom
if (levels_top_down) then
  lbgn=1
  linc=1
else
  lbgn=geoval%nval
  linc=-1
endif

ilev = lbgn
do jlev=1, c_nlev
  do jprofile=1, c_nprofiles
    iprofile = c_indx(jprofile) + 1
    if (iprofile<1 .or. iprofile> geoval%nprofiles) call abor1_ftn("ufo_geovals_fill: error iprofile")
    geoval%vals(ilev,iprofile) = c_vals(jprofile,jlev)
  enddo
  ilev = ilev + linc
enddo

end subroutine ufo_geovals_fill

! ------------------------------------------------------------------------------

subroutine ufo_geovals_fillad(self, varname, c_nprofiles, c_indx, c_nlev, c_vals, levels_top_down)
implicit none
type(ufo_geovals), intent(in) :: self
character(len=*), intent(in) :: varname
integer(c_int), intent(in) :: c_nprofiles
integer(c_int), intent(in) :: c_indx(c_nprofiles)
integer(c_int), intent(in) :: c_nlev
real(c_double), intent(inout) :: c_vals(c_nprofiles, c_nlev)
logical(c_bool), intent(in) :: levels_top_down

type(ufo_geoval), pointer :: geoval
integer :: jlev, jprofile, ilev, iprofile
integer :: lbgn, linc

if (.not.self%linit) call abor1_ftn("ufo_geovals_fillad: geovals not initialized")

call ufo_geovals_get_var(self, varname, geoval, ufo_geoval_sampled)

if (geoval%nval /= c_nlev) call abor1_ftn("ufo_geovals_fillad: incorrect number of levels")

! setting loop indices to ensure geovals are filled top to bottom
if (levels_top_down) then
  lbgn=1
  linc=1
else
  lbgn=geoval%nval
  linc=-1
endif

ilev = lbgn
do jlev = 1, geoval%nval
  do jprofile=1, c_nprofiles
    iprofile = c_indx(jprofile) + 1
    if (iprofile<1 .or. iprofile>geoval%nprofiles) call abor1_ftn("ufo_geovals_fillad: error iprofile")
    c_vals(jprofile, jlev) = geoval%vals(ilev,iprofile)
  enddo
  ilev = ilev + linc
enddo

end subroutine ufo_geovals_fillad

! ------------------------------------------------------------------------------

subroutine check(action, status)
use netcdf, only: nf90_noerr, nf90_strerror
implicit none

integer, intent (in) :: status
character (len=*), intent (in) :: action
character(max_string) :: err_msg

if(status /= nf90_noerr) then
  write(err_msg,*) "During action: ", trim(action), ", received error: ", trim(nf90_strerror(status))
  call abor1_ftn(err_msg)
end if

end subroutine check

! ------------------------------------------------------------------------------

subroutine ufo_geovals_print(self, iobs, format)
implicit none
! TODO(wsmigaj): iobs is currently treated as a profile index, but in reality it's a location index,
! and each location may correspond to multiple profiles.
type(ufo_geovals), intent(in), target :: self
integer, intent(in) :: iobs
integer(c_int), intent(in), optional :: format

type(ufo_geoval), pointer :: geoval
character(MAXVARLEN) :: varname
integer(c_int) :: actual_format
character(len=MAXVARLEN), pointer :: variables(:)  !< either the sampled or reduced variable list
integer :: ivar

if (present(format)) then
  actual_format = format
else
  actual_format = self%default_format
endif

if (actual_format == ufo_geoval_sampled) then
  variables => self%variables
else if (actual_format == ufo_geoval_reduced) then
  variables => self%reduced_variables
else
  call abor1_ftn("ufo_geovals_print: format must be either ufo_geoval_sampled or ufo_geoval_reduced")
endif

do ivar = 1, self%nvar
  varname = variables(ivar)
  call ufo_geovals_get_var(self, varname, geoval, actual_format)
  if (associated(geoval)) then
    print *, 'geoval test: ', trim(varname), geoval%nval, geoval%vals(:,iobs)
  else
    print *, 'geoval test: ', trim(varname), ' doesnt exist'
  endif
enddo

end subroutine ufo_geovals_print

! ------------------------------------------------------------------------------

!> Initialize and allocate a single geoval.
subroutine ufo_geovals_setup_geoval(self, igeoval, isampling_method, nvals, nprofiles)
implicit none

type(ufo_geovals), intent(inout) :: self
integer, intent(in)              :: igeoval
integer(c_size_t), intent(in)    :: isampling_method
integer(c_size_t), intent(in)    :: nvals
integer(c_size_t), intent(in)    :: nprofiles

self%geovals(igeoval)%nval = nvals
self%geovals(igeoval)%nprofiles = nprofiles
allocate(self%geovals(igeoval)%vals(self%geovals(igeoval)%nval, self%geovals(igeoval)%nprofiles))
self%geovals(igeoval)%vals(:,:) = 0.0
self%sampling_method_by_var(igeoval) = isampling_method

end subroutine ufo_geovals_setup_geoval

! ------------------------------------------------------------------------------

!> Initialize a single geoval without allocating it.
subroutine ufo_geovals_partial_setup_geoval(self, igeoval, isampling_method)
implicit none

type(ufo_geovals), intent(inout) :: self
integer, intent(in)           :: igeoval
integer(c_size_t), intent(in) :: isampling_method

self%geovals(igeoval)%nval = 0
self%geovals(igeoval)%nprofiles = 0
self%sampling_method_by_var(igeoval) = isampling_method

end subroutine ufo_geovals_partial_setup_geoval

! ------------------------------------------------------------------------------

elemental subroutine move_ufo_geoval(from, to)
implicit none
type(ufo_geoval), intent(inout) :: from
type(ufo_geoval), intent(out) :: to

call move_alloc(from%vals, to%vals)
to%nval = from%nval
to%nprofiles = from%nprofiles

end subroutine move_ufo_geoval

! ------------------------------------------------------------------------------

elemental subroutine move_location_sampling_method(from, to)
implicit none
type(ufo_location_sampling_method), intent(inout) :: from
type(ufo_location_sampling_method), intent(out) :: to

call move_alloc(from%paths_by_loc, to%paths_by_loc)
to%trivial = from%trivial
to%npaths = from%npaths

end subroutine move_location_sampling_method

! ------------------------------------------------------------------------------

subroutine resize_geovals(array, new_size)
implicit none
type(ufo_geoval), allocatable, intent(inout) :: array(:)
integer, intent(in) :: new_size

type(ufo_geoval), allocatable :: tmp(:)

if (new_size == size(array)) return

allocate(tmp(new_size))
call move(array(:), tmp(:size(array)))
call move_alloc(tmp, array)

end subroutine resize_geovals

! ------------------------------------------------------------------------------

subroutine resize_location_sampling_methods(array, new_size)
implicit none
type(ufo_location_sampling_method), allocatable, intent(inout) :: array(:)
integer, intent(in) :: new_size

type(ufo_location_sampling_method), allocatable :: tmp(:)

if (new_size == size(array)) return

allocate(tmp(new_size))
call move(array(:), tmp(:size(array)))
call move_alloc(tmp, array)

end subroutine resize_location_sampling_methods

! ------------------------------------------------------------------------------

subroutine resize_varnames(array, new_size)
implicit none
character(len=MAXVARLEN), allocatable, intent(inout) :: array(:)
integer, intent(in) :: new_size

character(len=MAXVARLEN), allocatable :: tmp(:)

if (new_size == size(array)) return

allocate(tmp(new_size))
tmp(:size(array)) = array(:)
call move_alloc(tmp, array)

end subroutine resize_varnames

! ------------------------------------------------------------------------------

subroutine resize_integers(array, new_size)
implicit none
integer, allocatable, intent(inout) :: array(:)
integer, intent(in) :: new_size

integer, allocatable :: tmp(:)

if (new_size == size(array)) return

allocate(tmp(new_size))
tmp(:size(array)) = array(:)
call move_alloc(tmp, array)

end subroutine resize_integers

! ------------------------------------------------------------------------------

end module ufo_geovals_mod
