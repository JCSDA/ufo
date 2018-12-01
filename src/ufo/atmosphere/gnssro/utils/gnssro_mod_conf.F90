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
  integer(c_int) :: ro_top_meter
  integer(c_int) :: use_compress
  character(len=255) :: obs_err
end type gnssro_conf
contains
!-------------------------------

subroutine gnssro_conf_setup(roconf, c_conf)
implicit none
type(gnssro_conf), intent(inout) :: roconf
type(c_ptr),      intent(in)     :: c_conf

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

if (config_element_exists(c_conf,"obs_err")) then
  roconf%obs_err = config_get_string(c_conf,len(roconf%obs_err),"obs_err")
else
  roconf%obs_err = "FILE"
endif

end subroutine gnssro_conf_setup


end module gnssro_mod_conf

