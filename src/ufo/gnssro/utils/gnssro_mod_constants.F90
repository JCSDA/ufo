!==========================================================================
module gnssro_mod_constants
!==========================================================================

use kinds
use iso_c_binding
use ufo_constants_mod
implicit none
public   :: gnssro_ref_constants
real(kind_real),            public :: n_a, n_b,n_c
integer, parameter,         public :: max_string    = 800

contains
subroutine gnssro_ref_constants(use_compress)
implicit none
integer(c_int),intent(in) :: use_compress

! cucurull 2010, Healy 2011
if (use_compress .eq. 1) then
       ! Constants for gpsro refractivity (Rueger 2002)
       n_a = 77.6890_kind_real    
       n_b = 3.75463e5_kind_real  
       n_c = 71.2952_kind_real     
else
       ! Constants for gpsro refractivity (Bevis et al 1994)
       n_a = 77.60_kind_real       
       n_b = 3.739e5_kind_real    
       n_c = 70.4_kind_real       
endif

n_c = n_c - n_a
return 

end subroutine gnssro_ref_constants

end module gnssro_mod_constants

