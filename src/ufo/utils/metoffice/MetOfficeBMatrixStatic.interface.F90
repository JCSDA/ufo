!-------------------------------------------------------------------------------
! (C) Crown Copyright 2021 Met Office
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
!-------------------------------------------------------------------------------

module ufo_metoffice_bmatrixstatic_mod_c

use fckit_configuration_module, only: fckit_configuration
use iso_c_binding
use kinds
use oops_variables_mod
use ufo_metoffice_bmatrixstatic_mod
use ufo_vars_mod

implicit none

private

#define LISTED_TYPE ufo_metoffice_bmatrixstatic

!> Linked list interface - defines registry_t type
#include "oops/util/linkedList_i.f"

!> Global registry
type(registry_t) :: ufo_metoffice_bmatrixstatic_registry

contains

!> Linked list implementation
#include "oops/util/linkedList_c.f"

!-------------------------------------------------------------------------------
subroutine ufo_metoffice_bmatrixstatic_setup_c(c_self, c_conf, &
                                               nbands, nelements) &
           bind(c, name='ufo_metoffice_bmatrixstatic_setup_f90')

implicit none
integer(c_int), intent(inout)            :: c_self
type(c_ptr), value, intent(in)           :: c_conf
integer(c_size_t), intent(inout)         :: nbands       ! number of latitude bands in B-matrix file
integer(c_size_t), intent(inout)         :: nelements    ! number of elements per B-matrix dimension

type(ufo_metoffice_bmatrixstatic), pointer :: self
type(fckit_configuration)                  :: f_conf
character(len=:), allocatable              :: str
character(len=:), allocatable              :: str_array(:)
logical                                    :: qtotal_flag
character(len=200)                         :: filepath
integer                                    :: varsize
character(len=200), allocatable            :: background_fields(:)

! Interface and setup
call ufo_metoffice_bmatrixstatic_registry % setup(c_self, self)

! Get filepath from configuration
f_conf = fckit_configuration(c_conf)
call f_conf % get_or_die("BMatrix", str)
filepath = str

! Get variables from configuration
varsize = f_conf % get_size("background fields")
allocate(background_fields(varsize))
call f_conf % get_or_die("background fields", str_array)
background_fields(1:varsize) = str_array

! Get qtotal from configuration
call f_conf % get_or_die("qtotal", qtotal_flag)

! Call Fortran
call self % setup(background_fields, trim(filepath), qtotal_flag)

! B-matrix has dimensions (nelements, nelements, nbands)
nbands = self % nbands
nelements = size(self % store, 1)

end subroutine ufo_metoffice_bmatrixstatic_setup_c

!-------------------------------------------------------------------------------
subroutine ufo_metoffice_bmatrixstatic_delete_c(c_self) &
           bind(c, name='ufo_metoffice_bmatrixstatic_delete_f90')

implicit none
integer(c_int), intent(inout) :: c_self

! Interface and setup
type(ufo_metoffice_bmatrixstatic), pointer :: self
call ufo_metoffice_bmatrixstatic_registry % get(c_self, self)

! Delete
call self % delete()
call ufo_metoffice_bmatrixstatic_registry % delete(c_self, self)

end subroutine ufo_metoffice_bmatrixstatic_delete_c

!-------------------------------------------------------------------------------
! Extract elements of B-matrix given its dimensions
subroutine ufo_metoffice_bmatrixstatic_getelements_c(c_self, nelements, nbands, south, north, &
           bmatrix_store) bind(C, name='ufo_metoffice_bmatrixstatic_getelements_f90')

implicit none
integer(c_int), intent(inout) :: c_self
integer(c_size_t), intent(in) :: nelements
integer(c_size_t), intent(in) :: nbands
real(c_float), intent(inout)  :: south(nbands)
real(c_float), intent(inout)  :: north(nbands)
real(c_float), intent(inout)  :: bmatrix_store(nelements,nelements,nbands)

type(ufo_metoffice_bmatrixstatic), pointer :: self
call ufo_metoffice_bmatrixstatic_registry % get(c_self, self)

south = self % south(1:nbands)
north = self % north(1:nbands)
bmatrix_store = real(self % store, kind=kind_single)

end subroutine ufo_metoffice_bmatrixstatic_getelements_c

!-------------------------------------------------------------------------------

end module ufo_metoffice_bmatrixstatic_mod_c
