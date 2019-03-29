! (C) Copyright 2017 UCAR
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 

!> Fortran module to handle gnssro observations-bending angle ROPP 1d operator

module ufo_gnssro_bndropp2d_mod_c
  
  use iso_c_binding
  use config_mod
  use ufo_gnssro_bndropp2d_mod
  use ufo_locs_mod
  use ufo_locs_mod_c
  use ufo_gnssro_2d_locs_mod

  implicit none
  private
  
#define LISTED_TYPE ufo_gnssro_BndROPP2D
  
  !> Linked list interface - defines registry_t type
#include "../../linkedList_i.f"
  
  !> Global registry
  type(registry_t) :: ufo_gnssro_BndROPP2D_registry
  
  ! ------------------------------------------------------------------------------
contains
  ! ------------------------------------------------------------------------------
  !> Linked list implementation
#include "../../linkedList_c.f"
  
! ------------------------------------------------------------------------------
  
subroutine ufo_gnssro_bndropp2d_setup_c(c_key_self, c_conf) bind(c,name='ufo_gnssro_bndropp2d_setup_f90')
implicit none
integer(c_int), intent(inout) :: c_key_self
type(c_ptr),    intent(in)    :: c_conf
    
type(ufo_gnssro_BndROPP2D), pointer :: self

call ufo_gnssro_BndROPP2D_registry%setup(c_key_self, self)
call self%setup(c_conf)
   
end subroutine ufo_gnssro_BndROPP2D_setup_c
  
! ------------------------------------------------------------------------------
  
subroutine ufo_gnssro_bndropp2d_delete_c(c_key_self) bind(c,name='ufo_gnssro_bndropp2d_delete_f90')
implicit none
integer(c_int), intent(inout) :: c_key_self
    
type(ufo_gnssro_BndROPP2D), pointer :: self

call ufo_gnssro_BndROPP2D_registry%delete(c_key_self,self)
    
end subroutine ufo_gnssro_bndropp2d_delete_c
  
! ------------------------------------------------------------------------------

subroutine ufo_gnssro_bndropp2d_simobs_c(c_key_self, c_key_geovals, c_obsspace, c_nobs, c_hofx, c_bias) bind(c,name='ufo_gnssro_bndropp2d_simobs_f90')

implicit none
integer(c_int),     intent(in)    :: c_key_self
integer(c_int),     intent(in)    :: c_key_geovals
type(c_ptr), value, intent(in)    :: c_obsspace
integer(c_int),     intent(in)    :: c_nobs
real(c_double),     intent(inout) :: c_hofx(c_nobs)
integer(c_int),      intent(in)   :: c_bias

type(ufo_gnssro_BndROPP2D),  pointer :: self

character(len=*), parameter :: myname_="ufo_gnssro_bndropp2d_simobs_c"
call ufo_gnssro_BndROPP2D_registry%get(c_key_self, self)
call self%opr_simobs(c_key_geovals, c_obsspace, c_hofx)

end subroutine ufo_gnssro_bndropp2d_simobs_c

! ------------------------------------------------------------------------------
subroutine ufo_gnssro_2d_locs_init_c(c_key_self, c_obsspace, c_t1, c_t2, c_conf) bind(c,name='ufo_gnssro_2d_locs_init_f90')
use datetime_mod
use gnssro_mod_conf, only : n_horiz_2d, res_2d, conf2d
implicit none
integer(c_int),     intent(inout)  :: c_key_self
type(c_ptr), value, intent(in)     :: c_obsspace
type(c_ptr),        intent(in)     :: c_t1, c_t2
type(c_ptr),        intent(in)     :: c_conf
type(ufo_locs),     pointer        :: self
integer, parameter            :: max_string = 800
character(max_string)         :: err_msg
type(datetime)   :: t1, t2
type(conf2d)     :: loc2dconf

call c_f_datetime(c_t1, t1)
call c_f_datetime(c_t2, t2)

if (config_element_exists(c_conf,"n_horiz")) then
  loc2dconf%n_horiz = config_get_int(c_conf,"n_horiz")
  if ( mod(loc2dconf%n_horiz,2) .eq. 0 ) then
    write(err_msg,*) 'ufo_gnssro_2d_locs_init_c error: n_horiz must be a odd number'
    call abor1_ftn(err_msg)
  end if
else
  loc2dconf%n_horiz = n_horiz_2d
endif

if (config_element_exists(c_conf,"res")) then
  loc2dconf%res = config_get_real(c_conf,"res")
else
  loc2dconf%res = res_2d
endif 

call ufo_locs_registry%get(c_key_self, self)
call ufo_gnssro_2d_locs_init(self, c_obsspace, t1, t2, loc2dconf )

end subroutine ufo_gnssro_2d_locs_init_c

! ------------------------------------------------------------------------------

end module ufo_gnssro_bndropp2d_mod_c
