module ufo_gnssro_2d_locs_mod

use iso_c_binding
use fckit_log_module, only : fckit_log
use kinds,            only : kind_real
use ufo_locs_mod
use gnssro_mod_conf

public:: ufo_gnssro_2d_locs_init

contains

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
subroutine ufo_gnssro_2d_locs_init(self, obss, t1, t2, roconf)
  use kinds
  use datetime_mod
  use obsspace_mod

  implicit none

  type(ufo_locs),      intent(inout) :: self
  type(c_ptr),  value, intent(in)    :: obss
  type(datetime),      intent(in)    :: t1, t2

  integer :: i, j, tw_nlocs,nlocs
  integer,         dimension(:), allocatable :: tw_indx
  real(kind_real), dimension(:), allocatable :: lon, lat
  type(datetime),  dimension(:), allocatable :: date_time

! gnss ro data 2d location  
  type(gnssro_conf),          intent(in)     :: roconf
  real(kind_real), dimension(roconf%n_horiz) :: plat_2d, plon_2d
  integer         :: kerror, n_horiz
  real(kind_real) :: dtheta

  dtheta  = roconf%dtheta 
  n_horiz = roconf%n_horiz

 ! Local copies pre binning
  nlocs = obsspace_get_nlocs(obss)

end subroutine ufo_gnssro_2d_locs_init

!--------------------------------
end module ufo_gnssro_2d_locs_mod
