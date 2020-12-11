! (C) Copyright 2017 UCAR
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 

!> Fortran module to handle gnssro observations-bending angle ROPP 1d operator

module ufo_gnssro_bndropp2d_mod_c
  
  use fckit_configuration_module, only: fckit_configuration 
  use ufo_gnssro_bndropp2d_mod
  use ufo_locs_mod
  use ufo_locs_mod_c
  use ufo_gnssro_2d_locs_mod

  implicit none
  private
  
#define LISTED_TYPE ufo_gnssro_BndROPP2D
  
  !> Linked list interface - defines registry_t type
#include "oops/util/linkedList_i.f"
  
  !> Global registry
  type(registry_t) :: ufo_gnssro_BndROPP2D_registry
  
  ! ------------------------------------------------------------------------------
contains
  ! ------------------------------------------------------------------------------
  !> Linked list implementation
#include "oops/util/linkedList_c.f"
  
! ------------------------------------------------------------------------------
  
subroutine ufo_gnssro_bndropp2d_setup_c(c_key_self, c_conf,c_size) bind(c,name='ufo_gnssro_bndropp2d_setup_f90')
implicit none
integer(c_int), intent(inout)  :: c_key_self
type(c_ptr), value, intent(in) :: c_conf
integer(c_int), intent(in)     :: c_size  ! obsspace vector length

type(ufo_gnssro_BndROPP2D), pointer :: self
type(fckit_configuration) :: f_conf

call ufo_gnssro_BndROPP2D_registry%setup(c_key_self, self)
f_conf = fckit_configuration(c_conf)

call self%setup(f_conf, c_size)
   
end subroutine ufo_gnssro_BndROPP2D_setup_c
  
! ------------------------------------------------------------------------------
  
subroutine ufo_gnssro_bndropp2d_delete_c(c_key_self) bind(c,name='ufo_gnssro_bndropp2d_delete_f90')
implicit none
integer(c_int), intent(inout) :: c_key_self
    
type(ufo_gnssro_BndROPP2D), pointer :: self

call ufo_gnssro_BndROPP2D_registry%delete(c_key_self,self)
    
end subroutine ufo_gnssro_bndropp2d_delete_c
  
! ------------------------------------------------------------------------------

subroutine ufo_gnssro_bndropp2d_simobs_c(c_key_self, c_key_geovals, c_obsspace, c_nobs, c_hofx) &
    bind(c,name='ufo_gnssro_bndropp2d_simobs_f90')

implicit none
integer(c_int),     intent(in)    :: c_key_self
integer(c_int),     intent(in)    :: c_key_geovals
type(c_ptr), value, intent(in)    :: c_obsspace
integer(c_int),     intent(in)    :: c_nobs
real(c_double),     intent(inout) :: c_hofx(c_nobs)
type(ufo_gnssro_BndROPP2D),  pointer :: self

character(len=*), parameter :: myname_="ufo_gnssro_bndropp2d_simobs_c"

call ufo_gnssro_BndROPP2D_registry%get(c_key_self, self)
call self%opr_simobs(c_key_geovals, c_obsspace, c_hofx)

end subroutine ufo_gnssro_bndropp2d_simobs_c

! ------------------------------------------------------------------------------
subroutine ufo_gnssro_2d_locs_init_c(c_key_self, c_key_locs, c_obsspace) bind(c,name='ufo_gnssro_2d_locs_init_f90')
implicit none
integer(c_int),     intent(in)     :: c_key_self  ! operator key
integer(c_int),     intent(inout)  :: c_key_locs  ! location key
type(c_ptr), value, intent(in)     :: c_obsspace

type(ufo_locs),              pointer :: locs
type(ufo_gnssro_BndROPP2D),  pointer :: self

integer, parameter            :: max_string = 800

call ufo_locs_registry%get(c_key_locs, locs)
call ufo_gnssro_BndROPP2D_registry%get(c_key_self, self)
call ufo_gnssro_2d_locs_init(self,locs, c_obsspace)

end subroutine ufo_gnssro_2d_locs_init_c

! ------------------------------------------------------------------------------

end module ufo_gnssro_bndropp2d_mod_c
