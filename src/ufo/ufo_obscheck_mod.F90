!
! (C) Copyright 2017 UCAR
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
!
module ufo_obscheck_mod

use iso_c_binding
use ufo_vars_mod
use ufo_obs_vectors
use ufo_locs_mod
use ufo_geovals_mod
use ufo_obs_data_mod, only: obs_data
use kinds

implicit none
private
public :: ufo_obscheck_registry

! ------------------------------------------------------------------------------

!> Fortran derived type to hold interpolated fields required by the obs operators
type :: ufo_obscheck
  integer :: nobs
  integer :: nvar
end type ufo_obscheck

#define LISTED_TYPE ufo_obscheck

!> Linked list interface - defines registry_t type
#include "linkedList_i.f"

!> Global registry
type(registry_t) :: ufo_obscheck_registry

! ------------------------------------------------------------------------------
contains
! ------------------------------------------------------------------------------
!> Linked list implementation
#include "linkedList_c.f"

! ------------------------------------------------------------------------------

subroutine c_ufo_obscheck_setup(c_key_self, c_conf) bind(c,name='ufo_obscheck_setup_f90')
   implicit none
   integer(c_int), intent(in) :: c_key_self
   type(c_ptr), intent(in)    :: c_conf
   type(ufo_obscheck), pointer :: self

   call ufo_obscheck_registry%init()
   call ufo_obscheck_registry%add(c_key_self)
   call ufo_obscheck_registry%get(c_key_self, self)

end subroutine c_ufo_obscheck_setup

! ------------------------------------------------------------------------------

subroutine c_ufo_obscheck_delete(c_key_self) bind(c,name='ufo_obscheck_delete_f90')
   implicit none
   integer(c_int), intent(inout) :: c_key_self

   type(ufo_obscheck), pointer :: self

   call ufo_obscheck_registry%get(c_key_self, self)
   call ufo_obscheck_registry%remove(c_key_self)

end subroutine c_ufo_obscheck_delete

! ------------------------------------------------------------------------------

subroutine c_ufo_postFilter_f90(c_key_geovals, c_key_hofx,c_key_obsspace) bind(c,name='ufo_postFilter_f90')

   implicit none
   integer(c_int), intent(in) :: c_key_geovals
   integer(c_int), intent(in) :: c_key_hofx
   integer(c_int), intent(in) :: c_key_obsspace
   type(ufo_geovals), pointer  :: geovals
   type(obs_vector), pointer :: hofx
   type(obs_data) :: obsdata

   write(*,*) '=======Start Post Filter (observation QC)========='
! Get pointers to geovals, hofx, and obsdata
   call ufo_geovals_registry%get(c_key_geovals,geovals)
   call ufo_obs_vect_registry%get(c_key_hofx,hofx)
!   call obs_data_registry%get(c_key_obsspace,obsdata)
!
! working
!
   write(*,*) 'read obs_vector  =========================='
   write(*,*) 'hofx%nobs=',hofx%nobs
   write(*,*) 'hofx%nobs=',hofx%values(1:10)
   if(geovals%lalloc .and. geovals%linit) then
      write(*,*) 'read geovals  =========================='
      write(*,*) 'nobs=',geovals%nobs
      write(*,*) 'nvar=',geovals%nvar
      write(*,*) 'geovals nval=',geovals%geovals(1:geovals%nvar)%nval
      write(*,*) 'geovals nobs=',geovals%geovals(1:geovals%nvar)%nobs
      write(*,*) 'geovals vals=',geovals%geovals(1)%vals(1:geovals%geovals(1)%nval,1)
      write(*,*) 'varaibles nv=',geovals%variables%nv
      write(*,*) 'varaibles fldnames==',geovals%variables%fldnames(1:geovals%variables%nv)
   else
      write(*,*) 'GeoVals has not allocated or initialized'
      write(*,*) 'lalloc=',geovals%lalloc
      write(*,*) 'linit=',geovals%linit
   endif
   write(*,*) '=======End Post Filter (observation QC)========='

end subroutine c_ufo_postFilter_f90

! ------------------------------------------------------------------------------

subroutine c_ufo_priorFilter_f90(c_key_obsspace) bind(c,name='ufo_priorFilter_f90')

    implicit none
    integer(c_int), intent(in) :: c_key_obsspace
    
    type(obs_data) :: obsdata

    write(*,*) '=======Start Prior Filter (observation QC)========='
    write(*,*) 'read obs==',obsdata%nobs
    write(*,*) '=======End Proir Filter (observation QC)========='

end subroutine c_ufo_priorFilter_f90

! ------------------------------------------------------------------------------

end module ufo_obscheck_mod
