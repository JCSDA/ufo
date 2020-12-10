!
! (C) Copyright 2017-2018 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
!
module ufo_geovals_mod

use fckit_configuration_module, only: fckit_configuration
use iso_c_binding
use ufo_vars_mod
use kinds
use obsspace_mod
use missing_values_mod

use fckit_mpi_module, only: fckit_mpi_comm, fckit_mpi_sum

implicit none
private
integer, parameter :: max_string=800

public :: ufo_geovals, ufo_geoval
public :: ufo_geovals_get_var, ufo_geovals_put_var
public :: ufo_geovals_default_constr, ufo_geovals_setup, ufo_geovals_delete, ufo_geovals_print
public :: ufo_geovals_zero, ufo_geovals_random, ufo_geovals_dotprod, ufo_geovals_scalmult
public :: ufo_geovals_profmult
public :: ufo_geovals_reorderzdir
public :: ufo_geovals_assign, ufo_geovals_add, ufo_geovals_diff, ufo_geovals_abs
public :: ufo_geovals_split, ufo_geovals_merge
public :: ufo_geovals_minmaxavg, ufo_geovals_normalize, ufo_geovals_maxloc, ufo_geovals_schurmult
public :: ufo_geovals_read_netcdf, ufo_geovals_write_netcdf
public :: ufo_geovals_rms, ufo_geovals_copy, ufo_geovals_copy_one
public :: ufo_geovals_analytic_init

private :: ufo_geovals_reset_sec_arg

! ------------------------------------------------------------------------------

!> type to hold interpolated field for one variable, one observation
type :: ufo_geoval
  real(kind_real), allocatable :: vals(:,:) !< values (nval, nlocs)
  integer :: nval = 0                !< number of values in profile
  integer :: nlocs = 0               !< number of observations
end type ufo_geoval

!> type to hold interpolated fields required by the obs operators
type :: ufo_geovals
  integer :: nlocs  = 0          !< number of observations
  integer :: nvar  = 0           !< number of variables (supposed to be the
                                 !  same for same obs operator

  type(ufo_geoval), allocatable :: geovals(:)  !< array of interpolated
                                               !  vertical profiles for all obs (nvar)

  character(len=MAXVARLEN), allocatable :: variables(:)  !< variable list

  real(c_double) :: missing_value !< obsspace missing value mark

  logical :: linit = .false.     !< .true. if all the ufo_geoval arrays inside geovals
                                 !  were allocated and have data
end type ufo_geovals

! ------------------------------------------------------------------------------
contains
! ------------------------------------------------------------------------------

subroutine ufo_geovals_default_constr(self)
implicit none
type(ufo_geovals), intent(inout) :: self

self%nlocs = 0
self%missing_value = missing_value(0.0)
self%nvar = 0
self%linit = .false.

end subroutine ufo_geovals_default_constr


subroutine ufo_geovals_setup(self, vars, nlocs)
use oops_variables_mod
implicit none
type(ufo_geovals), intent(inout) :: self
type(oops_variables), intent(in) :: vars
integer, intent(in) :: nlocs

integer :: ivar
type(fckit_configuration) :: f_vars

call ufo_geovals_delete(self)
self%nlocs = nlocs
self%missing_value = missing_value(self%missing_value)

self%nvar = vars%nvars()
allocate(self%geovals(self%nvar))
allocate(self%variables(self%nvar))
do ivar = 1, self%nvar
  self%variables(ivar) = vars%variable(ivar)
  self%geovals(ivar)%nlocs = nlocs
  self%geovals(ivar)%nval = 0
enddo

end subroutine ufo_geovals_setup

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
self%nvar = 0
self%nlocs = 0
self%linit = .false.

end subroutine ufo_geovals_delete

! ------------------------------------------------------------------------------

subroutine ufo_geovals_get_var(self, varname, geoval)
implicit none
type(ufo_geovals), target, intent(in)    :: self
character(len=*), intent(in) :: varname
type(ufo_geoval), pointer, intent(inout)    :: geoval

character(len=*), parameter :: myname_="ufo_geovals_get_var"

character(max_string) :: err_msg
integer :: ivar, jv

geoval => NULL()
if (.not. self%linit) then
   !return
endif

ivar = ufo_vars_getindex(self%variables, varname)

if (ivar < 0) then
  write(0,*)'ufo_geovals_get_var looking for ',trim(varname),' in:'
  do jv=1,self%nvar
    write(0,*)'ufo_geovals_get_var ',jv,trim(self%variables(jv))
  enddo
  write(err_msg,*) myname_, " ", trim(varname), ' doesnt exist'
  call abor1_ftn(err_msg)
else
  geoval => self%geovals(ivar)
endif

end subroutine ufo_geovals_get_var

! ------------------------------------------------------------------------------

subroutine ufo_geovals_put_var(self, varname, geoval,k)
type(ufo_geovals),intent(inout) :: self
character(len=*),    intent(in) :: varname
type(ufo_geoval),    intent(in) :: geoval
integer,             intent(in) :: k

integer :: ivar

ivar = ufo_vars_getindex(self%variables, varname)
self%geovals(ivar)%vals(k,:)=geoval%vals(k,:)

end subroutine ufo_geovals_put_var

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
   do jo = 1, self%nlocs
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
  do jo=1,self%nlocs
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
integer :: jv, jo

if (.not. self%linit) then
  call abor1_ftn("ufo_geovals_profmult: geovals not allocated")
endif

do jv=1,self%nvar
  do jo=1,self%nlocs
     self%geovals(jv)%vals(:,jo) = values(jo) * self%geovals(jv)%vals(:,jo)
  enddo
enddo

end subroutine ufo_geovals_profmult

! ------------------------------------------------------------------------------


subroutine ufo_geovals_assign(self, rhs)
implicit none
type(ufo_geovals), intent(inout) :: self
type(ufo_geovals), intent(in) :: rhs
integer :: jv, jo, jz
integer :: iv
character(max_string) :: err_msg

if (.not. self%linit) then
  call abor1_ftn("ufo_geovals_scalmult: geovals not allocated")
endif
if (.not. rhs%linit) then
  call abor1_ftn("ufo_geovals_scalmult: geovals not allocated")
endif

if (self%nlocs /= rhs%nlocs) then
  call abor1_ftn("ufo_geovals_assign: nlocs different between lhs and rhs")
endif

do jv=1,self%nvar
  iv = ufo_vars_getindex(rhs%variables, self%variables(jv))
  if (iv < 0) then
    write(err_msg,*) 'ufo_geovals_assign: var ', trim(self%variables(jv)), ' doesnt exist in rhs'
    call abor1_ftn(trim(err_msg))
  endif
  if (self%geovals(jv)%nval /= rhs%geovals(iv)%nval) then
    write(err_msg,*) 'ufo_geovals_assign: nvals for var ', trim(self%variables(jv)), ' are different in lhs and rhs'
    call abor1_ftn(trim(err_msg))
  endif
  do jo=1,self%nlocs
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
logical :: do_flip = .false.     !< .true. if all the ufo_geoval arrays inside geovals

if (.not. self%linit) then
  call abor1_ftn("ufo_geovals_reorderzdir: geovals not allocated")
endif

! Get vertical coordinate variable
call ufo_geovals_get_var(self, varname, geoval)
if (.not. associated(geoval)) then
  write(err_msg, *) 'ufo_geovals_reorderzdir: geoval vertical coordinate variable ', trim(varname), ' doesnt exist'
endif

! Check if reorder variables is necessary based on the direction defined by zdir
if ((zdir == "bottom2top" .and. geoval%vals(1,1) < geoval%vals(geoval%nval,1)) .or. &
    (zdir == "top2bottom" .and. geoval%vals(1,1) > geoval%vals(geoval%nval,1))) then
   do_flip = .true.
else if (zdir /= "bottom2top" .or. zdir /= "top2bottom") then
  write(err_msg, *) 'ufo_geovals_reorderzdir: z-coordinate direction ', trim(zdir), ' not defined'
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
type(ufo_geovals), intent(inout) :: self
type(ufo_geovals), intent(in) :: other
integer :: jv, jo, jz
integer :: iv
character(max_string) :: err_msg

if (.not. self%linit) then
  call abor1_ftn("ufo_geovals_add: geovals not allocated")
endif
if (.not. other%linit) then
  call abor1_ftn("ufo_geovals_add: geovals not allocated")
endif

if (self%nlocs /= other%nlocs) then
  call abor1_ftn("ufo_geovals_add: nlocs different between lhs and rhs")
endif

do jv=1,self%nvar
  iv = ufo_vars_getindex(other%variables, self%variables(jv))
  if (iv .ne. -1) then !Only add if exists in RHS
    if (self%geovals(jv)%nval /= other%geovals(iv)%nval) then
      write(err_msg,*) 'ufo_geovals_add: nvals for var ', trim(self%variables(jv)), ' are different in lhs and rhs'
      call abor1_ftn(trim(err_msg))
    endif
    do jo=1,self%nlocs
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
type(ufo_geovals), intent(inout) :: self
type(ufo_geovals), intent(in) :: other
integer :: jv, jo, jz
integer :: iv
character(max_string) :: err_msg

if (.not. self%linit) then
  call abor1_ftn("ufo_geovals_diff: geovals not allocated")
endif
if (.not. other%linit) then
  call abor1_ftn("ufo_geovals_diff: geovals not allocated")
endif

if (self%nlocs /= other%nlocs) then
  call abor1_ftn("ufo_geovals_diff: nlocs different between lhs and rhs")
endif

do jv=1,self%nvar
  iv = ufo_vars_getindex(other%variables, self%variables(jv))
  if (iv .ne. -1) then !Only subtract if exists in RHS
    if (self%geovals(jv)%nval /= other%geovals(iv)%nval) then
      write(err_msg,*) 'ufo_geovals_diff: nvals for var ', trim(self%variables(jv)), ' are different in lhs and rhs'
      call abor1_ftn(trim(err_msg))
    endif
    do jo=1,self%nlocs
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
type(ufo_geovals), intent(inout) :: self
type(ufo_geovals), intent(in) :: other
integer :: jv, jo, jz
integer :: iv
character(max_string) :: err_msg

if (.not. self%linit) then
  call abor1_ftn("ufo_geovals_schurmult: geovals not allocated")
endif
if (.not. other%linit) then
  call abor1_ftn("ufo_geovals_schurmult: geovals not allocated")
endif

if (self%nlocs /= other%nlocs) then
  call abor1_ftn("ufo_geovals_schurmult: nlocs different between lhs and rhs")
endif

do jv=1,self%nvar
  iv = ufo_vars_getindex(other%variables, self%variables(jv))
  if (iv .ne. -1) then !Only mult if exists in RHS
    if (self%geovals(jv)%nval /= other%geovals(iv)%nval) then
      write(err_msg,*) 'ufo_geovals_schurmult: nvals for var ', trim(self%variables(jv)), ' are different in lhs and rhs'
      call abor1_ftn(trim(err_msg))
    endif
    do jo=1,self%nlocs
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
integer :: jv

if (.not. self%linit) then
  call abor1_ftn("ufo_geovals_copy: geovals not defined")
endif

call ufo_geovals_delete(other)

other%nlocs = self%nlocs
other%nvar = self%nvar
allocate(other%variables(other%nvar))
other%variables(:) = self%variables(:)

allocate(other%geovals(other%nvar))
do jv = 1, other%nvar
  other%geovals(jv)%nval = self%geovals(jv)%nval
  other%geovals(jv)%nlocs = self%geovals(jv)%nlocs
  allocate(other%geovals(jv)%vals(other%geovals(jv)%nval, other%geovals(jv)%nlocs))
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
integer :: jv

if (.not. other%linit) then
  call abor1_ftn("ufo_geovals_copy_one: geovals not defined")
endif

call ufo_geovals_delete(self)

self%nlocs = 1
self%nvar = other%nvar
allocate(self%variables(self%nvar))
self%variables(:) = other%variables(:)

allocate(self%geovals(self%nvar))
do jv = 1, self%nvar
  self%geovals(jv)%nval = other%geovals(jv)%nval
  self%geovals(jv)%nlocs = 1
  allocate(self%geovals(jv)%vals(self%geovals(jv)%nval, self%geovals(jv)%nlocs))
  self%geovals(jv)%vals(:,self%nlocs) = other%geovals(jv)%vals(:,loc_index)
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

subroutine ufo_geovals_analytic_init(self, locs, ic)
use ufo_locs_mod, only : ufo_locs
use dcmip_initial_conditions_test_1_2_3, only : test1_advection_deformation, &
                                  test1_advection_hadley, test3_gravity_wave
use dcmip_initial_conditions_test_4, only : test4_baroclinic_wave

implicit none
type(ufo_geovals), intent(inout) :: self
type(ufo_locs), intent(in)       :: locs
character(*), intent(in)         :: ic

real(kind_real) :: pi = acos(-1.0_kind_real)
real(kind_real) :: deg_to_rad,rlat, rlon
real(kind_real) :: p0, kz, u0, v0, w0, t0, phis0, ps0, rho0, hum0
real(kind_real) :: q1, q2, q3, q4
integer :: ivar, iloc, ival

if (.not. self%linit) then
  call abor1_ftn("ufo_geovals_analytic_init: geovals not defined")
endif

! The last variable should be the ln pressure coordinate.  That's
! where we get the height information for the analytic init
if (trim(self%variables(self%nvar)) /= trim(var_prs)) then
  call abor1_ftn("ufo_geovals_analytic_init: pressure coordinate not defined")
endif

deg_to_rad = pi/180.0_kind_real

do ivar = 1, self%nvar-1

   do iloc = 1, self%geovals(ivar)%nlocs

      ! convert lat and lon to radians
      rlat = deg_to_rad * locs%lat(iloc)
      rlon = deg_to_rad*modulo(locs%lon(iloc)+180.0_kind_real,360.0_kind_real) - pi

      do ival = 1, self%geovals(ivar)%nval

         ! obtain height from the existing GeoVaLs object, which should be an
         ! output of the State::getValues() method
         ! should be delivered in units of Pa
         p0 = self%geovals(self%nvar)%vals(ival,iloc)

         init_option: select case (trim(ic))

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

         ! currently only temperture is implemented
         if (trim(self%variables(ivar)) == trim(var_tv)) then
            ! Warning: we may need a conversion from temperature to
            ! virtual temperture here
            self%geovals(ivar)%vals(ival,iloc) = t0
         endif

      enddo
   enddo
enddo

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
  call abor1_ftn("ufo_geovals_normalize: reference geovals object must have the same variables as the original")
endif


do jv=1,self%nvar

   !> Compute normalization factors for the errors based on the rms amplitude of
   !! each variable across all of the selected locations.  Use the "other" GeoVaLs
   !! object as a reference, since this may be the exact analytic answer

   over_nloc = 1.0_kind_real / &
        (real(other%nlocs,kind_real)*real(other%geovals(jv)%nval,kind_real))

   vrms = 0.0_kind_real
   do jo = 1, other%nlocs
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
   do jo=1,self%nlocs
      do jz = 1, self%geovals(jv)%nval
         self%geovals(jv)%vals(jz,jo) = norm*self%geovals(jv)%vals(jz,jo)
      enddo
   enddo
enddo

end subroutine ufo_geovals_normalize

! ------------------------------------------------------------------------------

subroutine ufo_geovals_dotprod(self, other, gprod, f_comm)
implicit none
real(kind_real), intent(inout) :: gprod
type(ufo_geovals), intent(in) :: self, other
integer :: ivar, iobs, ival, nval
real(kind_real) :: prod

type(fckit_mpi_comm), intent(in) :: f_comm

if (.not. self%linit) then
  call abor1_ftn("ufo_geovals_dotprod: geovals not allocated")
endif

if (.not. other%linit) then
  call abor1_ftn("ufo_geovals_dotprod: geovals not allocated")
endif

! just something to put in (dot product of the 1st var and 1st element in the profile
prod=0.0
do ivar = 1, self%nvar
  nval = self%geovals(ivar)%nval
  do ival = 1, nval
     do iobs = 1, self%nlocs
      if ((self%geovals(ivar)%vals(ival,iobs) .ne. self%missing_value) .and. &
          (other%geovals(ivar)%vals(ival,iobs) .ne. self%missing_value)) then
        prod = prod + self%geovals(ivar)%vals(ival,iobs) * &
                      other%geovals(ivar)%vals(ival,iobs)
      endif
    enddo
  enddo
enddo

!Get global dot product
call f_comm%allreduce(prod,gprod,fckit_mpi_sum())

end subroutine ufo_geovals_dotprod

!-------------------------------------------------------------------------------
subroutine ufo_geovals_reset_sec_arg(self, other, nlocs)
implicit none
type(ufo_geovals), intent(in) :: self
type(ufo_geovals), intent(inout) :: other
integer, intent(in) :: nlocs
integer :: ivar

if (other%linit) call abor1_ftn("ufo_geovals_reset_sec_arg: other already have data")

other%nlocs = nlocs
other%nvar = self%nvar
other%missing_value = self%missing_value
allocate(other%variables(self%nvar))
allocate(other%geovals(self%nvar))
do ivar = 1, self%nvar
  other%variables(ivar) = self%variables(ivar)
  other%geovals(ivar)%nlocs = nlocs
  other%geovals(ivar)%nval = self%geovals(ivar)%nval
  allocate(other%geovals(ivar)%vals(self%geovals(ivar)%nval, nlocs))
  other%geovals(ivar)%vals(:,:) = 0.0
enddo
other%linit = .false.

end subroutine ufo_geovals_reset_sec_arg
! ------------------------------------------------------------------------------

subroutine ufo_geovals_split(self, other1, other2)
implicit none
type(ufo_geovals), intent(in) :: self
type(ufo_geovals), intent(inout) :: other1
type(ufo_geovals), intent(inout) :: other2

integer :: ivar, iobs

if (.not. self%linit) &
  call abor1_ftn("ufo_geovals_split: geovals self is not allocated or has no data")

if (other1%linit .or. other2%linit) &
  call abor1_ftn("ufo_geovals_split: geovals other1 or other2 already have data")

call ufo_geovals_delete(other1)
call ufo_geovals_delete(other2)
call ufo_geovals_reset_sec_arg(self, other1, self%nlocs/2)
call ufo_geovals_reset_sec_arg(self, other2, self%nlocs - self%nlocs/2)

do ivar = 1, self%nvar
  do iobs = 1, self%nlocs/2
    other1%geovals(ivar)%vals(:,iobs) = self%geovals(ivar)%vals(:,iobs)
  enddo
  do iobs = self%nlocs/2 + 1, self%nlocs
    other2%geovals(ivar)%vals(:,iobs - self%nlocs/2) = self%geovals(ivar)%vals(:,iobs)
  enddo
enddo
other1%linit = .true.
other2%linit = .true.

end subroutine ufo_geovals_split
! ------------------------------------------------------------------------------

subroutine ufo_geovals_merge(self, other1, other2)
implicit none
type(ufo_geovals), intent(inout) :: self
type(ufo_geovals), intent(in) :: other1
type(ufo_geovals), intent(in) :: other2

integer :: ivar, iobs

if ((.not. other1%linit) .or. (.not. other2%linit)) &
  call abor1_ftn("ufo_geovals_merge: geovals other1 or other2 is not allocated or has no data")

call ufo_geovals_delete(self)
call ufo_geovals_reset_sec_arg(other1, self, other1%nlocs + other2%nlocs)

do ivar = 1, self%nvar
  do iobs = 1, other1%nlocs
    self%geovals(ivar)%vals(:,iobs) = other1%geovals(ivar)%vals(:,iobs)
  enddo
  do iobs = other1%nlocs + 1, self%nlocs
    self%geovals(ivar)%vals(:,iobs) = &
      other2%geovals(ivar)%vals(:,iobs - other1%nlocs)
  enddo
enddo
self%linit = .true.

end subroutine ufo_geovals_merge
! ------------------------------------------------------------------------------

subroutine ufo_geovals_minmaxavg(self, kobs, kvar, pmin, pmax, prms)
implicit none
integer, intent(inout) :: kobs
integer, intent(in) :: kvar
real(kind_real), intent(inout) :: pmin, pmax, prms
type(ufo_geovals), intent(in) :: self
integer :: jo, jz, jv

jv = kvar+1
kobs = 0
pmin = huge(pmin)
pmax = -huge(pmax)
prms = 0.0_kind_real
do jo = 1, self%nlocs
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
   do jo = 1, self%nlocs

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

subroutine ufo_geovals_read_netcdf(self, filename, loc_multiplier, c_obspace, vars)
use netcdf
use oops_variables_mod
implicit none
type(ufo_geovals), intent(inout)  :: self
character(max_string), intent(in) :: filename
integer, intent(in)               :: loc_multiplier
type(c_ptr), intent(in)           :: c_obspace
type(oops_variables), intent(in)  :: vars

integer :: nlocs, gv_all_nlocs, nlocs_var
integer :: nval
integer :: obs_nlocs
integer :: obs_all_nlocs
integer :: iloc
integer :: jloc, jloc_start, jloc_end
integer :: iloc_new

integer :: ncid, dimid, varid, vartype, ndims
integer, dimension(3) :: dimids
integer :: ivar
integer :: ierr

character(max_string) :: err_msg
character(len=30) :: obs_nlocs_str
character(len=30) :: geo_nlocs_str

integer(c_size_t), allocatable, dimension(:) :: dist_indx
integer(c_size_t), allocatable, dimension(:) :: obs_dist_indx

real, allocatable :: field2d(:,:), field1d(:)

! open netcdf file
call check('nf90_open', nf90_open(trim(filename),nf90_nowrite,ncid))

! find how many locs are in the file
ierr = nf90_inq_dimid(ncid, "nlocs", dimid)
if(ierr /= nf90_noerr) then
  write(err_msg,*) "Error: Dimension nlocs not found in ", trim(filename)
  call abor1_ftn(err_msg)
endif
call check('nf90_inq_dimid', nf90_inq_dimid(ncid, "nlocs", dimid))
call check('nf90_inquire_dimension', nf90_inquire_dimension(ncid, dimid, len = gv_all_nlocs))

!> round-robin distribute the observations to PEs
!> Calculate how many obs. on each PE
obs_all_nlocs = obsspace_get_gnlocs(c_obspace)
obs_nlocs = obsspace_get_nlocs(c_obspace)
allocate(obs_dist_indx(obs_nlocs))
call obsspace_get_index(c_obspace, obs_dist_indx)

! loc_multiplier specifies how many locations in the geovals file per
! single location in the obs file. There needs to be at least
! loc_multiplier * obs_all_nlocs locations in the geovals file.

if (gv_all_nlocs .lt. (loc_multiplier * obs_all_nlocs)) then
  write(obs_nlocs_str, *) loc_multiplier * obs_all_nlocs
  write(geo_nlocs_str, *) gv_all_nlocs
  write(err_msg,'(7a)') &
     "Error: Number of locations in the geovals file (", &
     trim(adjustl(geo_nlocs_str)), ") must be greater than or equal to ", &
     "the product of loc_multiplier and number of locations in the ", &
     "obs file (", trim(adjustl(obs_nlocs_str)), ")"
  call abor1_ftn(err_msg)
endif

! We have enough locations in the geovals file to cover the span of the
! number of locations in the obs file. Generate the dist_indx according
! to the loc_multiplier and obs_nlocs values.
if (loc_multiplier >= 0) then
  nlocs = loc_multiplier * obs_nlocs
  allocate(dist_indx(nlocs))
  iloc_new = 1
  do iloc = 1,obs_nlocs
    jloc_start = ((obs_dist_indx(iloc) - 1) * loc_multiplier) + 1
    jloc_end = obs_dist_indx(iloc) * loc_multiplier
    do jloc = jloc_start, jloc_end
      dist_indx(iloc_new) = jloc
      iloc_new = iloc_new + 1
    enddo
  enddo
else
  nlocs = - loc_multiplier * obs_nlocs
  allocate(dist_indx(nlocs))
  iloc_new = 1
  do jloc = 1, - loc_multiplier
    do iloc = 1, obs_nlocs
      dist_indx(iloc_new) = obs_dist_indx(iloc) + (jloc - 1) * obs_all_nlocs
      iloc_new = iloc_new + 1
    enddo
  enddo
end if

! allocate geovals structure
call ufo_geovals_setup(self, vars, nlocs)

do ivar = 1, self%nvar

  ierr = nf90_inq_varid(ncid, self%variables(ivar), varid)
  if(ierr /= nf90_noerr) then
    write(err_msg,*) "Error: Variable ", trim(self%variables(ivar)), " not found in ", trim(filename)
    call abor1_ftn(err_msg)
  endif

  call check('nf90_inquire_variable', nf90_inquire_variable(ncid, varid, xtype = vartype, &
                                         ndims = ndims, dimids = dimids))
  !> read 1d variable
  if (ndims == 1) then
    call check('nf90_inquire_dimension', nf90_inquire_dimension(ncid, dimids(1), len = nlocs_var))
    if (nlocs_var /= gv_all_nlocs) then
      call abor1_ftn('ufo_geovals_read_netcdf: var dim /= gv_all_nlocs')
    endif
    nval = 1
    !> allocate geoval for this variable
    self%geovals(ivar)%nval = nval
    allocate(self%geovals(ivar)%vals(nval,nlocs))

    allocate(field1d(nlocs_var))
    call check('nf90_get_var', nf90_get_var(ncid, varid, field1d))
    self%geovals(ivar)%vals(1,:) = field1d(dist_indx)
    deallocate(field1d)
  !> read 2d variable
  elseif (ndims == 2) then
    call check('nf90_inquire_dimension', nf90_inquire_dimension(ncid, dimids(1), len = nval))
    call check('nf90_inquire_dimension', nf90_inquire_dimension(ncid, dimids(2), len = nlocs_var))
    if (nlocs_var /= gv_all_nlocs) then
      call abor1_ftn('ufo_geovals_read_netcdf: var dim /= gv_all_nlocs')
    endif
    !> allocate geoval for this variable
    self%geovals(ivar)%nval = nval
    allocate(self%geovals(ivar)%vals(nval,nlocs))
    allocate(field2d(nval, nlocs_var))
    call check('nf90_get_var', nf90_get_var(ncid, varid, field2d))
    self%geovals(ivar)%vals(:,:) = field2d(:,dist_indx)
    deallocate(field2d)
  !> only 1d & 2d vars
  else
    call abor1_ftn('ufo_geovals_read_netcdf: can only read 1d and 2d fields')
  endif

  ! set the missing value equal to IODA missing_value
  where (self%geovals(ivar)%vals > 1.0e08) self%geovals(ivar)%vals = self%missing_value

enddo

if (allocated(dist_indx)) deallocate(dist_indx)
if (allocated(obs_dist_indx)) deallocate(obs_dist_indx)

self%linit = .true.

call check('nf90_close', nf90_close(ncid))

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
call check('nf90_def_dim', nf90_def_dim(ncid,'nlocs',self%nlocs, dimid_nlocs))
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

subroutine ufo_geovals_print(self, iobs)
implicit none
type(ufo_geovals), intent(in) :: self
integer, intent(in) :: iobs

type(ufo_geoval), pointer :: geoval
character(MAXVARLEN) :: varname
integer :: ivar

do ivar = 1, self%nvar
  varname = self%variables(ivar)
  call ufo_geovals_get_var(self, varname, geoval)
  if (associated(geoval)) then
    print *, 'geoval test: ', trim(varname), geoval%nval, geoval%vals(:,iobs)
  else
    print *, 'geoval test: ', trim(varname), ' doesnt exist'
  endif
enddo

end subroutine ufo_geovals_print

! ------------------------------------------------------------------------------

end module ufo_geovals_mod
