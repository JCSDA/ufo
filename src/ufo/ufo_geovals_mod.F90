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
use type_distribution, only: random_distribution
use obsspace_mod

use fckit_mpi_module, only: fckit_mpi_comm, fckit_mpi_sum

implicit none
private
integer, parameter :: max_string=800

public :: ufo_geovals, ufo_geoval, ufo_geovals_get_var
public :: ufo_geovals_init, ufo_geovals_setup, ufo_geovals_delete, ufo_geovals_print
public :: ufo_geovals_zero, ufo_geovals_random, ufo_geovals_dotprod, ufo_geovals_scalmult
public :: ufo_geovals_assign, ufo_geovals_add, ufo_geovals_diff, ufo_geovals_abs
public :: ufo_geovals_minmaxavg, ufo_geovals_normalize, ufo_geovals_maxloc
public :: ufo_geovals_read_netcdf, ufo_geovals_rms, ufo_geovals_copy
public :: ufo_geovals_analytic_init
public :: ufo_geovals_allocone

! ------------------------------------------------------------------------------

!> type to hold interpolated field for one variable, one observation
type :: ufo_geoval
  real(kind_real), allocatable :: vals(:,:) !< values (nval, nobs)
  integer :: nval                !< number of values in profile
  integer :: nobs                !< number of observations
end type ufo_geoval

!> type to hold interpolated fields required by the obs operators
type :: ufo_geovals
  integer :: nobs                !< number of observations
  integer :: nvar                !< number of variables (supposed to be the
                                 !  same for same obs operator

  type(ufo_geoval), allocatable :: geovals(:)  !< array of interpolated
                                               !  vertical profiles for all obs (nvar)

  type(ufo_vars) :: variables    !< variables list

  real(c_double) :: missing_value !< obsspace missing value mark

  logical :: lalloc              !< .true. if type was initialized and allocated
                                 !  (only geovals are allocated, not the arrays
                                 !   inside of the ufo_geoval type)
  logical :: linit               !< .true. if all the ufo_geoval arrays inside geovals
                                 !  were allocated and have data
end type ufo_geovals

! ------------------------------------------------------------------------------
contains
! ------------------------------------------------------------------------------

subroutine ufo_geovals_init(self)
implicit none
type(ufo_geovals), intent(inout) :: self

self%lalloc = .false.
self%linit  = .false.
self%nvar = 0
self%nobs = 0

end subroutine ufo_geovals_init

! ------------------------------------------------------------------------------

subroutine ufo_geovals_setup(self, vars, nobs)
implicit none
type(ufo_geovals), intent(inout) :: self
type(ufo_vars), intent(in) :: vars
integer, intent(in) :: nobs

integer :: ivar

call ufo_geovals_delete(self)
self%nobs = nobs
self%nvar = vars%nv
self%missing_value = obspace_missing_value()
call ufo_vars_clone(vars, self%variables) 
allocate(self%geovals(self%nvar))
do ivar = 1, self%nvar
  self%geovals(ivar)%nobs = nobs
  self%geovals(ivar)%nval = 0
enddo
self%lalloc = .true.
end subroutine ufo_geovals_setup

! ------------------------------------------------------------------------------

subroutine ufo_geovals_delete(self)
implicit none
type(ufo_geovals), intent(inout) :: self

integer :: ivar

if (self%linit) then
  do ivar = 1, self%nvar
    deallocate(self%geovals(ivar)%vals)
  enddo
  self%linit = .false.
endif
if (self%lalloc) then
  deallocate(self%geovals)
  self%lalloc = .false.
  self%nvar = 0
  self%nobs = 0
endif

end subroutine ufo_geovals_delete

! ------------------------------------------------------------------------------

subroutine ufo_geovals_get_var(self, varname, geoval)
implicit none
type(ufo_geovals), target, intent(in)    :: self
character(MAXVARLEN), intent(in) :: varname
type(ufo_geoval), pointer, intent(inout)    :: geoval

character(len=*), parameter :: myname_="ufo_geovals_get_var"

character(max_string) :: err_msg
integer :: ivar

geoval => NULL()
if (.not. self%lalloc .or. .not. self%linit) then
   !return
endif

ivar = ufo_vars_getindex(self%variables, varname)

if (ivar < 0) then
  write(err_msg,*) myname_, trim(varname), ' doesnt exist'
  call abor1_ftn(err_msg)
else
  geoval => self%geovals(ivar)
endif

end subroutine ufo_geovals_get_var

! ------------------------------------------------------------------------------

subroutine ufo_geovals_allocone(self) 
implicit none
type(ufo_geovals), intent(inout) :: self
integer :: ivar

if (.not. self%lalloc) then
  call abor1_ftn("ufo_geovals_zero: geovals not allocated")
endif

do ivar = 1,self%nvar
  self%geovals(ivar)%nval = 1
  allocate(self%geovals(ivar)%vals(1,self%nobs))
enddo
self%linit = .true.

end subroutine ufo_geovals_allocone

! ------------------------------------------------------------------------------

subroutine ufo_geovals_zero(self) 
implicit none
type(ufo_geovals), intent(inout) :: self
integer :: ivar

if (.not. self%lalloc) then
  call abor1_ftn("ufo_geovals_zero: geovals not allocated")
endif
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

if (.not. self%lalloc) then
  call abor1_ftn("ufo_geovals_abs: geovals not allocated")
endif
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

if (.not. self%lalloc) then
  call abor1_ftn("ufo_geovals_rms: geovals not allocated")
endif
if (.not. self%linit) then
  call abor1_ftn("ufo_geovals_rms: geovals not initialized")
endif
vrms=0.0_kind_real
N=0.0_kind_real
do jv = 1, self%nvar
   do jo = 1, self%nobs
      vrms = vrms + Sum(self%geovals(jv)%vals(:,jo)**2)
      N=N+self%geovals(jv)%nval
   enddo   
enddo

vrms = sqrt(vrms/N)

end subroutine ufo_geovals_rms

! ------------------------------------------------------------------------------

subroutine ufo_geovals_random(self) 
use random_vectors_mod
implicit none
type(ufo_geovals), intent(inout) :: self
integer :: ivar

if (.not. self%lalloc) then
  call abor1_ftn("ufo_geovals_random: geovals not allocated")
endif
if (.not. self%linit) then
  call abor1_ftn("ufo_geovals_random: geovals not initialized")
endif
do ivar = 1, self%nvar
  call random_vector(self%geovals(ivar)%vals)
enddo

end subroutine ufo_geovals_random

! ------------------------------------------------------------------------------

subroutine ufo_geovals_scalmult(self, zz) 
implicit none
type(ufo_geovals), intent(inout) :: self
real(kind_real), intent(in) :: zz
integer :: jv, jo, jz

if (.not. self%lalloc .or. .not. self%linit) then
  call abor1_ftn("ufo_geovals_scalmult: geovals not allocated")
endif

do jv=1,self%nvar
  do jo=1,self%nobs
    do jz = 1, self%geovals(jv)%nval
      self%geovals(jv)%vals(jz,jo) = zz * self%geovals(jv)%vals(jz,jo)
    enddo
  enddo
enddo

end subroutine ufo_geovals_scalmult

! ------------------------------------------------------------------------------

subroutine ufo_geovals_assign(self, rhs) 
implicit none
type(ufo_geovals), intent(inout) :: self
type(ufo_geovals), intent(in) :: rhs
integer :: jv, jo, jz
integer :: iv
character(max_string) :: err_msg

if (.not. self%lalloc .or. .not. self%linit) then
  call abor1_ftn("ufo_geovals_scalmult: geovals not allocated")
endif
if (.not. rhs%lalloc .or. .not. rhs%linit) then
  call abor1_ftn("ufo_geovals_scalmult: geovals not allocated")
endif

if (self%nobs /= rhs%nobs) then
  call abor1_ftn("ufo_geovals_assign: nobs different between lhs and rhs")
endif

do jv=1,self%nvar
  iv = ufo_vars_getindex(rhs%variables, self%variables%fldnames(jv))
  if (iv < 0) then
    write(err_msg,*) 'ufo_geovals_assign: var ', trim(self%variables%fldnames(jv)), ' doesnt exist in rhs'
    call abor1_ftn(trim(err_msg))
  endif
  if (self%geovals(jv)%nval /= rhs%geovals(iv)%nval) then
    write(err_msg,*) 'ufo_geovals_assign: nvals for var ', trim(self%variables%fldnames(jv)), ' are different in lhs and rhs'
    call abor1_ftn(trim(err_msg))
  endif
  do jo=1,self%nobs
    do jz = 1, self%geovals(jv)%nval
      self%geovals(jv)%vals(jz,jo) = rhs%geovals(iv)%vals(jz,jo)
    enddo
  enddo
enddo

end subroutine ufo_geovals_assign

! ------------------------------------------------------------------------------
!> Sum of two GeoVaLs objects

subroutine ufo_geovals_add(self, other) 
implicit none
type(ufo_geovals), intent(inout) :: self
type(ufo_geovals), intent(in) :: other
integer :: jv, jo, jz
integer :: iv
character(max_string) :: err_msg

if (.not. self%lalloc .or. .not. self%linit) then
  call abor1_ftn("ufo_geovals_add: geovals not allocated")
endif
if (.not. other%lalloc .or. .not. other%linit) then
  call abor1_ftn("ufo_geovals_add: geovals not allocated")
endif

if (self%nobs /= other%nobs) then
  call abor1_ftn("ufo_geovals_assign: nobs different between lhs and rhs")
endif

do jv=1,self%nvar
  iv = ufo_vars_getindex(other%variables, self%variables%fldnames(jv))
  if (iv .ne. -1) then !Only add if exists in RHS
    if (self%geovals(jv)%nval /= other%geovals(iv)%nval) then
      write(err_msg,*) 'ufo_geovals_assign: nvals for var ', trim(self%variables%fldnames(jv)), ' are different in lhs and rhs'
      call abor1_ftn(trim(err_msg))
    endif
    do jo=1,self%nobs
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

if (.not. self%lalloc .or. .not. self%linit) then
  call abor1_ftn("ufo_geovals_diff: geovals not allocated")
endif
if (.not. other%lalloc .or. .not. other%linit) then
  call abor1_ftn("ufo_geovals_diff: geovals not allocated")
endif

if (self%nobs /= other%nobs) then
  call abor1_ftn("ufo_geovals_assign: nobs different between lhs and rhs")
endif

do jv=1,self%nvar
  iv = ufo_vars_getindex(other%variables, self%variables%fldnames(jv))
  if (iv .ne. -1) then !Only subtract if exists in RHS
    if (self%geovals(jv)%nval /= other%geovals(iv)%nval) then
      write(err_msg,*) 'ufo_geovals_assign: nvals for var ', trim(self%variables%fldnames(jv)), ' are different in lhs and rhs'
      call abor1_ftn(trim(err_msg))
    endif
    do jo=1,self%nobs
      do jz = 1, self%geovals(jv)%nval
        self%geovals(jv)%vals(jz,jo) = self%geovals(jv)%vals(jz,jo) - other%geovals(iv)%vals(jz,jo)
      enddo
    enddo
  endif
enddo

end subroutine ufo_geovals_diff

! ------------------------------------------------------------------------------
!> Copy one GeoVaLs object into another
!!

subroutine ufo_geovals_copy(self, other) 
implicit none
type(ufo_geovals), intent(in) :: self
type(ufo_geovals), intent(inout) :: other
integer :: jv

if (.not. self%lalloc .or. .not. self%linit) then
  call abor1_ftn("ufo_geovals_copy: geovals not defined")
endif

call ufo_geovals_delete(other)

call ufo_vars_clone(self%variables,other%variables)

other%nobs = self%nobs
other%nvar = self%nvar

allocate(other%geovals(other%nvar))
do jv = 1, other%nvar
  other%geovals(jv)%nval = self%geovals(jv)%nval
  other%geovals(jv)%nobs = self%geovals(jv)%nobs
  allocate(other%geovals(jv)%vals(other%geovals(jv)%nval, other%geovals(jv)%nobs))
  other%geovals(jv)%vals(:,:) = self%geovals(jv)%vals(:,:)
enddo

other%lalloc = .true.
other%linit = .true.

end subroutine ufo_geovals_copy

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
use ioda_locs_mod, only : ioda_locs
use dcmip_initial_conditions_test_1_2_3, only : test1_advection_deformation, &
                                  test1_advection_hadley, test3_gravity_wave  
use dcmip_initial_conditions_test_4, only : test4_baroclinic_wave

implicit none
type(ufo_geovals), intent(inout) :: self
type(ioda_locs), intent(in)      :: locs
character(*), intent(in)         :: ic

real(kind_real) :: pi = acos(-1.0_kind_real)
real(kind_real) :: deg_to_rad,rlat, rlon
real(kind_real) :: p0, kz, u0, v0, w0, t0, phis0, ps0, rho0, hum0
real(kind_real) :: q1, q2, q3, q4
integer :: ivar, iloc, ival

if (.not. self%lalloc .or. .not. self%linit) then
  call abor1_ftn("ufo_geovals_analytic_init: geovals not defined")
endif

! The last variable should be the ln pressure coordinate.  That's
! where we get the height information for the analytic init
if (trim(self%variables%fldnames(self%nvar)) /= trim(var_prsl)) then
  call abor1_ftn("ufo_geovals_analytic_init: pressure coordinate not defined")
endif

deg_to_rad = pi/180.0_kind_real

do ivar = 1, self%nvar-1

   do iloc = 1, self%geovals(ivar)%nobs

      ! convert lat and lon to radians
      rlat = deg_to_rad * locs%lat(iloc)
      rlon = deg_to_rad*modulo(locs%lon(iloc)+180.0_kind_real,360.0_kind_real) - pi
    
      do ival = 1, self%geovals(ivar)%nval                  

         ! obtain height from the existing GeoVaLs object, which should be an
         ! output of the State::getValues() method
         ! convert from KPa (ufo standard) to Pa (dcmip standard)
         p0 = exp(self%geovals(self%nvar)%vals(ival,iloc))*1.0e3_kind_real

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
         if (trim(self%variables%fldnames(ivar)) == trim(var_tv)) then
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

if (.not. self%lalloc .or. .not. self%linit) then
  call abor1_ftn("ufo_geovals_normalize: geovals not allocated")
endif
if (.not. other%lalloc .or. .not. other%linit) then
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
        (real(other%nobs,kind_real)*real(other%geovals(jv)%nval,kind_real))

   vrms = 0.0_kind_real
   do jo = 1, other%nobs
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
   do jo=1,self%nobs
      do jz = 1, self%geovals(jv)%nval
         self%geovals(jv)%vals(jz,jo) = norm*self%geovals(jv)%vals(jz,jo)
      enddo
   enddo
enddo

end subroutine ufo_geovals_normalize

! ------------------------------------------------------------------------------

subroutine ufo_geovals_dotprod(self, other, gprod) 
implicit none
real(kind_real), intent(inout) :: gprod
type(ufo_geovals), intent(in) :: self, other
integer :: ivar, iobs, ival, nval
real(kind_real) :: prod

type(fckit_mpi_comm) :: f_comm

f_comm = fckit_mpi_comm()

if (.not. self%lalloc .or. .not. self%linit) then
  call abor1_ftn("ufo_geovals_dotprod: geovals not allocated")
endif

if (.not. other%lalloc .or. .not. other%linit) then
  call abor1_ftn("ufo_geovals_dotprod: geovals not allocated")
endif

! just something to put in (dot product of the 1st var and 1st element in the profile
prod=0.0
do ivar = 1, self%nvar
  nval = self%geovals(ivar)%nval
  do ival = 1, nval
     do iobs = 1, self%nobs
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

! ------------------------------------------------------------------------------

subroutine ufo_geovals_minmaxavg(self, kobs, pmin, pmax, prms) 
implicit none
integer, intent(inout) :: kobs
real(kind_real), intent(inout) :: pmin, pmax, prms
type(ufo_geovals), intent(in) :: self

kobs = self%nobs
pmin=minval(self%geovals(1)%vals)
pmax=maxval(self%geovals(1)%vals)
prms=0. !sqrt(sum(self%values(:,:)**2)/real(self%nobs*self%nvar,kind_real))

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

if (.not. self%lalloc .or. .not. self%linit) then
  call abor1_ftn("ufo_geovals_maxloc: geovals not allocated")
endif

mxval = 0.0_kind_real
iobs = 1
ivar = 1

do jv = 1,self%nvar
   do jo = 1, self%nobs

      vrms = 0.0_kind_real
      do jz = 1, self%geovals(jv)%nval
         vrms = vrms + self%geovals(jv)%vals(jz,jo)**2
      enddo
      vrms = sqrt(vrms/real(self%geovals(jv)%nval,kind_real))
      if (vrms > mxval) then
         mxval = vrms
         iobs = jo
         ivar = jv
      endif

   enddo
enddo

end subroutine ufo_geovals_maxloc

! ------------------------------------------------------------------------------

subroutine ufo_geovals_read_netcdf(self, filename, vars, t1, t2)
use datetime_mod
use duration_mod

USE netcdf, ONLY: NF90_FLOAT, NF90_DOUBLE, NF90_INT
use nc_diag_read_mod, only: nc_diag_read_get_global_attr
use nc_diag_read_mod, only: nc_diag_read_get_var
use nc_diag_read_mod, only: nc_diag_read_get_dim
use nc_diag_read_mod, only: nc_diag_read_get_var_dims, nc_diag_read_check_var
use nc_diag_read_mod, only: nc_diag_read_get_var_type
use nc_diag_read_mod, only: nc_diag_read_init, nc_diag_read_close

implicit none
type(ufo_geovals), intent(inout)  :: self
character(max_string), intent(in) :: filename
type(ufo_vars), intent(in)        :: vars
type(datetime), intent(in)        :: t1
type(datetime), intent(in)        :: t2

integer :: iunit, ivar, nobs, nval, fvlen
integer :: nvardim, vartype
integer, allocatable, dimension(:) :: vardims

real(kind_real), allocatable :: fieldr2d(:,:), fieldr1d(:)
real, allocatable :: fieldf2d(:,:), fieldf1d(:)
integer, allocatable :: fieldi2d(:,:), fieldi1d(:)

character(max_string) :: err_msg

type(random_distribution) :: distribution
integer, allocatable, dimension(:) :: dist_indx

integer :: date_time_attr
type(datetime) :: refdate
integer :: tw_nobs
integer, allocatable :: tw_indx(:)
real(kind_real), allocatable :: time_offset(:)
type(duration), dimension(:), allocatable :: dt
type(datetime), dimension(:), allocatable :: t

integer :: i
integer :: j

!DEBUG:!
character(max_string) :: tstr

call datetime_to_string(t1, tstr)
print*, "DEBUG: read_file: t1: ", trim(tstr)
call datetime_to_string(t2, tstr)
print*, "DEBUG: read_file: t2: ", trim(tstr)
!DEBUG:!

! open netcdf file and read dimensions
call nc_diag_read_init(filename, iunit)
if (allocated(vardims)) deallocate(vardims)
call nc_diag_read_get_var_dims(iunit, vars%fldnames(1), nvardim, vardims)
if (nvardim .eq. 1) then
  fvlen = vardims(1)
else
  fvlen = vardims(2)
endif

!> round-robin distribute the observations to PEs
!> Calculate how many obs. on each PE
distribution=random_distribution(fvlen)
nobs=distribution%nobs_pe()
allocate(dist_indx(nobs))
dist_indx = distribution%indx

! Strip out obs that fall outside the timing window.

! Read in the date_time attribute and for a datetime object
call nc_diag_read_get_global_attr(iunit, "date_time", date_time_attr)

! Create the datetime object with a dummy date, then set it from the
! date_time attribute.
call datetime_create("1000-01-01T00:00:00Z", refdate)
call datetime_from_ifs(refdate, date_time_attr/100, 0)

!DEBUG:!
call datetime_to_string(refdate, tstr)
print*, "DEBUG: read_file: refdate: ", trim(tstr)
!DEBUG:!

! Read in the time variable
allocate(time_offset(nobs))
vartype = nc_diag_read_get_var_type(iunit, "time")
if (vartype == NF90_DOUBLE) then
   allocate(fieldr1d(fvlen))
   call nc_diag_read_get_var(iunit, "time", fieldr1d)
   time_offset(:) = fieldr1d(dist_indx)
   deallocate(fieldr1d)
elseif (vartype == NF90_FLOAT) then
   allocate(fieldf1d(fvlen))
   call nc_diag_read_get_var(iunit, "time", fieldf1d)  
   time_offset(:) = dble(fieldf1d(dist_indx))
   deallocate(fieldf1d)
endif
print*, "DEBUG: read_file: time: ", time_offset(1:3)

! Remove any obs that are outside the timing window.
allocate(tw_indx(nobs))
allocate(dt(nobs))
allocate(t(nobs))

do i = 1, nobs
  dt(i) = int(3600*time_offset(i))
  t(i) = refdate
  call datetime_update(t(i), dt(i))
enddo

! Find number of locations in this timeframe
tw_nobs = 0
do i = 1, nobs
  if (t(i) > t1 .and. t(i) <= t2) then
    tw_nobs = tw_nobs + 1
    tw_indx(tw_nobs) = i
  endif
enddo

print*, "DEBUG: nobs, tw_nobs: ", nobs, tw_nobs
print*, "DEBUG: dist_indx: ", dist_indx
print*, "DEBUG: tw_indx: ", tw_indx

! Adjust dist_indx if tw_nobs is different than original nobs
if (tw_nobs .ne. nobs) then
  nobs = tw_nobs
  if (allocated(dist_indx)) deallocate(dist_indx)
  allocate(dist_indx(nobs))
  do i = 1, nobs
    dist_indx(i) = distribution%indx(tw_indx(i))
  enddo
endif
print*, "DEBUG: nobs, size(dist_indx) (after): ", nobs, size(dist_indx)
print*, "DEBUG: dist_indx (after): ", dist_indx

deallocate(time_offset)
deallocate(tw_indx)
deallocate(dt)
deallocate(t)

! allocate geovals structure
call ufo_geovals_init(self)
call ufo_geovals_setup(self, vars, nobs)

do ivar = 1, vars%nv
  if (.not. nc_diag_read_check_var(iunit, vars%fldnames(ivar))) then
     write(err_msg,*) 'ufo_geovals_read_netcdf: var ', trim(vars%fldnames(ivar)), ' doesnt exist'
     call abor1_ftn(trim(err_msg))
  endif
  !> get dimensions of variable
  if (allocated(vardims)) deallocate(vardims)
  call nc_diag_read_get_var_dims(iunit, vars%fldnames(ivar), nvardim, vardims)
  !> get variable type
  vartype = nc_diag_read_get_var_type(iunit, vars%fldnames(ivar))
  !> read 1d vars (only double precision and integer for now)
  if (nvardim == 1) then
    if (vardims(1) /= fvlen) call abor1_ftn('ufo_geovals_read_netcdf: var dim /= fvlen')
    nval = 1

    !> allocate geoval for this variable
    self%geovals(ivar)%nval = nval
    allocate(self%geovals(ivar)%vals(nval,nobs))

    if (vartype == NF90_DOUBLE) then
       allocate(fieldr1d(vardims(1)))
       call nc_diag_read_get_var(iunit, vars%fldnames(ivar), fieldr1d)
       self%geovals(ivar)%vals(1,:) = fieldr1d(dist_indx)
       deallocate(fieldr1d)
    elseif (vartype == NF90_FLOAT) then
       allocate(fieldf1d(vardims(1)))
       call nc_diag_read_get_var(iunit, vars%fldnames(ivar), fieldf1d)  
       self%geovals(ivar)%vals(1,:) = dble(fieldf1d(dist_indx))
       deallocate(fieldf1d)
    elseif (vartype == NF90_INT) then
       allocate(fieldi1d(vardims(1)))
       call nc_diag_read_get_var(iunit, vars%fldnames(ivar), fieldi1d)
       self%geovals(ivar)%vals(1,:) = fieldi1d(dist_indx)
       deallocate(fieldi1d)
    else
       call abor1_ftn('ufo_geovals_read_netcdf: can only read double, float and int')
    endif

    ! set the missing value equal to IODA missing_value
    where (self%geovals(ivar)%vals(1,:) > 1.0e08) self%geovals(ivar)%vals(1,:) = self%missing_value

  !> read 2d vars (only double precision and integer for now)
  elseif (nvardim == 2) then
    if (vardims(2) /= fvlen) call abor1_ftn('ufo_geovals_read_netcdf: var dim /= fvlen')
    nval = vardims(1)

    !> allocate geoval for this variable
    self%geovals(ivar)%nval = nval
    allocate(self%geovals(ivar)%vals(nval,nobs))

    if (vartype == NF90_DOUBLE) then
       allocate(fieldr2d(vardims(1), vardims(2)))
       call nc_diag_read_get_var(iunit, vars%fldnames(ivar), fieldr2d)
       self%geovals(ivar)%vals = fieldr2d(:,dist_indx)
       deallocate(fieldr2d)
    elseif (vartype == NF90_FLOAT) then
       allocate(fieldf2d(vardims(1), vardims(2)))
       call nc_diag_read_get_var(iunit, vars%fldnames(ivar), fieldf2d)
       self%geovals(ivar)%vals = fieldf2d(:,dist_indx)
       deallocate(fieldf2d)
    elseif (vartype == NF90_INT) then
       allocate(fieldi2d(vardims(1), vardims(2)))
       call nc_diag_read_get_var(iunit, vars%fldnames(ivar), fieldi2d)
       self%geovals(ivar)%vals = fieldi2d(:,dist_indx)
       deallocate(fieldi2d)
    else
       call abor1_ftn('ufo_geovals_read_netcdf: can only read double, float and int')
    endif

    ! set the missing value equal to IODA missing_value
    where (self%geovals(ivar)%vals > 1.0e08) self%geovals(ivar)%vals = self%missing_value
    
  !> only 1d & 2d vars
  else
    call abor1_ftn('ufo_geovals_read_netcdf: can only read 1d and 2d fields')
  endif
enddo

self%linit = .true.

call nc_diag_read_close(filename)

end subroutine ufo_geovals_read_netcdf

! ------------------------------------------------------------------------------

subroutine ufo_geovals_print(self, iobs)
implicit none
type(ufo_geovals), intent(in) :: self
integer, intent(in) :: iobs

type(ufo_geoval), pointer :: geoval
character(MAXVARLEN) :: varname
integer :: ivar

do ivar = 1, self%nvar
  varname = self%variables%fldnames(ivar)
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
