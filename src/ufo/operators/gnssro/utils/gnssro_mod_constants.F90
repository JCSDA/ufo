!==========================================================================
module gnssro_mod_constants
!==========================================================================

use kinds
use iso_c_binding
use ufo_constants_mod

implicit none

private
public   :: gnssro_ref_constants

real(kind_real),            public :: n_a, n_b, n_c
integer, parameter,         public :: max_string    = 800
integer, parameter,         public :: MAXVARLEN     = 20
real(kind_real), parameter, public :: r1em6 = 1.0e-6_kind_real
real(kind_real), parameter, public :: r1em3 = 1.0e-3_kind_real
real(kind_real), parameter, public :: ds    = 10000.0_kind_real
real(kind_real), parameter, public :: crit_gradRefr = 157.0_kind_real !criteria for the refractivity gradient

contains
subroutine gnssro_ref_constants(use_compress)
implicit none
integer(c_int),intent(in) :: use_compress

! cucurull 2010, Healy 2011
if (use_compress .eq. 1) then
       ! Constants for gpsro refractivity (Rueger 2002)
       n_a = 0.776890_kind_real    
       n_b = 3.75463e3_kind_real  
       n_c = 0.712952_kind_real     
else
       ! Constants for gpsro refractivity (Bevis et al 1994)
       n_a = 0.7760_kind_real       
       n_b = 3.739e3_kind_real    
       n_c = 0.704_kind_real       
endif

n_c = n_c - n_a
return 

end subroutine gnssro_ref_constants

end module gnssro_mod_constants

