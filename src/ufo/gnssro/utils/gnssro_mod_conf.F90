!==========================================================================
module gnssro_mod_conf
!==========================================================================
use iso_c_binding
use config_mod
use kinds
use obsspace_mod
use gnssro_mod_constants
implicit none
private
public   :: gnssro_conf
public   :: gnssro_conf_setup

type gnssro_conf
  integer(c_int)     :: ro_top_meter
  integer(c_int)     :: use_compress
  integer(c_int)     :: n_horiz
  real(kind_real)    :: res
  real(kind_real)    :: dtheta
  character(len=20)  :: vertlayer
end type gnssro_conf

!--------- ropp2d location default parameters-----------------
integer(c_int),  parameter, public :: n_horiz_2d = 31   !should be odd number
real(kind_real), parameter, public :: res_2d     = 40.0 !km

contains
!-------------------------------

subroutine gnssro_conf_setup(roconf, c_conf)
implicit none
type(gnssro_conf), intent(inout)  :: roconf
type(c_ptr),       intent(in)     :: c_conf

roconf%ro_top_meter  = config_get_int(c_conf,  "ro_top_meter",  30000 )
roconf%use_compress  = config_get_int(c_conf,  "use_compress",  1 )
roconf%n_horiz       = config_get_int(c_conf,  "n_horiz",       n_horiz_2d)
roconf%res           = config_get_real(c_conf, "res",           res_2d)
roconf%dtheta        = roconf%res/mean_earth_rad
roconf%vertlayer     = config_get_string(c_conf,len(roconf%vertlayer), "vertlayer", "full")

end subroutine gnssro_conf_setup


end module gnssro_mod_conf

