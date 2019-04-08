!==========================================================================
module gnssro_mod_constants
!==========================================================================

use kinds
use iso_c_binding

implicit none
public   :: gnssro_ref_constants
real(kind_real), parameter, public :: zero    = 0.0_kind_real
real(kind_real), parameter, public :: one     = 1.0_kind_real
real(kind_real), parameter, public :: two     = 2.0_kind_real
real(kind_real), parameter, public :: half    = 0.5_kind_real
real(kind_real), parameter, public :: rad2deg = 57.29577954572      ! copy from fv3jedi_constants.f90
real(kind_real), parameter, public :: deg2rad =  0.01745329251
real(kind_real), parameter, public :: grav    = 9.80665e+0_kind_real
real(kind_real), parameter, public :: t0c     = 2.7315e+2_kind_real ! temperature at zero celsius     (K)
real(kind_real), parameter, public :: rd     = 2.8705e2_kind_real
real(kind_real), parameter, public :: rv     = 4.6150e2_kind_real
real(kind_real), parameter, public :: rd_over_rv = rd/rv
real(kind_real), parameter, public :: rv_over_rd = rv/rd
real(kind_real), parameter, public :: rd_over_g  = rd/grav
real(kind_real),            public :: n_a, n_b,n_c
real(kind_real), parameter, public :: mean_earth_rad = 6371.0
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

