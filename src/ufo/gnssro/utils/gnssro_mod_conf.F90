!==========================================================================
module gnssro_mod_conf
!==========================================================================
use iso_c_binding
use config_mod
use kinds
use obsspace_mod

implicit none
private
public   :: gnssro_conf
public   :: gnssro_conf_setup

type gnssro_conf
  integer(c_int)     :: ro_top_meter
  integer(c_int)     :: use_compress
  character(len=255) :: obserr_method
  integer(c_int)     :: n_horiz
  real(kind_real)    :: res
end type gnssro_conf

!--------- ropp2d location default parameters-----------------
integer(c_int),  parameter, public :: n_horiz_2d = 31   !should be odd number
real(kind_real), parameter, public :: res_2d     = 40.0 !km

contains
!-------------------------------

subroutine gnssro_conf_setup(roconf, c_conf)
implicit none
type(gnssro_conf), intent(inout)  :: roconf
type(c_ptr),      intent(in)  :: c_conf

if (config_element_exists(c_conf,"ro_top_meter")) then
  roconf%ro_top_meter  = config_get_int(c_conf,"ro_top_meter" )
else
  roconf%ro_top_meter  = 30000 !meter
endif

if (config_element_exists(c_conf,"use_compress")) then
  roconf%use_compress  = config_get_int(c_conf,"use_compress" )
else
  roconf%use_compress  =  1
endif

if (config_element_exists(c_conf,"obserr_method")) then
  roconf%obserr_method = config_get_string(c_conf,len(roconf%obserr_method),"obserr_method")
else
  roconf%obserr_method = "FILE"
endif
if (config_element_exists(c_conf,"n_horiz")) then
  roconf%n_horiz = config_get_int(c_conf,"n_horiz")
else 
  roconf%n_horiz = n_horiz_2d
endif
if (config_element_exists(c_conf,"res")) then
  roconf%res = config_get_real(c_conf,"res")
else
  roconf%res = res_2d
endif

end subroutine gnssro_conf_setup


end module gnssro_mod_conf

