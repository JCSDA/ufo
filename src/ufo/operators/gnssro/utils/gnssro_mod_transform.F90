!==========================================================================
module gnssro_mod_transform
!==========================================================================

use kinds
use, intrinsic :: iso_c_binding
use gnssro_mod_constants
use ufo_constants_mod
implicit none

private
public :: geometric2geop, geop2geometric, compute_refractivity


contains

! ------------------------------
! variable converting between geopotential and geometric heights using  MJ Mahoney's (2001)
! Parameters from WGS-84 model software inside GPS receivers.
! copy from GSI
subroutine  geometric2geop(Latitude,geometricZ, geopotentialH )
implicit none
real(kind_real), intent(in)  :: Latitude,   geometricZ
real(kind_real), intent(out) :: geopotentialH
real(kind_real)              :: sino, termg, termr ! local variables

sino          = sin(deg2rad*Latitude)**2
termg         = grav_equator*( (one+somigliana*sino)/sqrt(one-eccentricity*eccentricity*sino) )
termr         = semi_major_axis / (one + flattening + grav_ratio - two*flattening*sino)
geopotentialH = (termg/grav) * ((termr*geometricZ)/(termr+geometricZ))  ! meter

end subroutine  geometric2geop


subroutine geop2geometric(latitude, geopotentialH, geometricZ, dzdh_jac)
implicit none
! calculate observation geometric height using  MJ Mahoney's (2001), eq(23)
real(kind_real),intent(in)   :: latitude,   geopotentialH
real(kind_real),intent(out)  :: geometricZ
real(kind_real),intent(out), optional   ::dzdh_jac !dzdh's jacobian
real(kind_real)            :: sino
real(kind_real)            :: termg, termr, termrg

sino          = sin(deg2rad*latitude)**2
termg         = grav_equator*( (one+somigliana*sino)/sqrt(one-eccentricity*eccentricity*sino) )
termr         = semi_major_axis / (one + flattening + grav_ratio - two*flattening*sino)
termrg        = termg/grav*termr

geometricZ = termr*geopotentialH/(termrg-geopotentialH)

if ( present(dzdh_jac)) then
  dzdh_jac   = termr/(termrg-geopotentialH) + (termr*geopotentialH)/(termrg-geopotentialH)**2
end if

end subroutine geop2geometric
! ------------------------------

subroutine compute_refractivity(temperature, specH, pressure,refr, use_compress)
implicit none
real(kind_real), intent(in)  :: temperature, specH, pressure
real(kind_real), intent(out) :: refr
integer(c_int),  intent(in)  :: use_compress
real(kind_real) :: refr1,refr2,refr3, tfact

! constants needed to compute refractivity
  call gnssro_ref_constants(use_compress)

  tfact = (1-(rd_over_rv))*specH+(rd_over_rv)
  refr1 = n_a*pressure/temperature
  refr2 = n_b*specH*pressure/(temperature**2*tfact)
  refr3 = n_c*specH*pressure/(temperature*tfact)
  refr  = refr1 + refr2 + refr3

end subroutine compute_refractivity
! ------------------------------

end module gnssro_mod_transform

