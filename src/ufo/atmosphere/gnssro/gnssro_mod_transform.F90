!==========================================================================
module gnssro_mod_transform 
!==========================================================================

use kinds
use gnssro_mod_constants, only: one, two, deg2rad, grav
contains

subroutine  geometric2geop(obsLat,geometricZ, geopotentialH )
! calculate observation geopotential height using  MJ Mahoney's (2001)
! copy from GSI, to  Convert geometric height at observation to geopotential height 
! Parameters from WGS-84 model software inside GPS receivers.
real(kind_real), intent(in)  :: obsLat,   geometricZ
real(kind_real), intent(out) :: geopotentialH
real(kind_real)              :: sino, termg, termr ! local variables

real(kind_real), parameter ::  semi_major_axis = 6378.1370e3_kind_real     !                     (m)
real(kind_real), parameter ::  semi_minor_axis = 6356.7523142e3_kind_real  !                     (m)
real(kind_real), parameter ::  grav_polar      = 9.8321849378_kind_real    !                     (m/s2)
real(kind_real), parameter ::  grav_equator    = 9.7803253359_kind_real    !                     (m/s2) 
real(kind_real), parameter ::  earth_omega     = 7.292115e-5_kind_real     !                     (rad/s)
real(kind_real), parameter ::  grav_constant   = 3.986004418e14_kind_real  !      
real(kind_real), parameter ::  flattening = (semi_major_axis-semi_minor_axis)/semi_major_axis
real(kind_real), parameter ::  somigliana = (semi_minor_axis/semi_major_axis) * (grav_polar/grav_equator) - one
real(kind_real), parameter ::  grav_ratio = (earth_omega*earth_omega * &
                                           semi_major_axis*semi_major_axis * semi_minor_axis) / grav_constant
real(kind_real), parameter :: eccentricity = sqrt(semi_major_axis**2 - semi_minor_axis**2)/semi_major_axis


sino          = sin(deg2rad*obsLat)
termg         = grav_equator*( (one+somigliana*sino)/sqrt(one-eccentricity*eccentricity*sino) )
termr         = semi_major_axis / (one + flattening + grav_ratio - two*flattening*sino)
geopotentialH = (termg/grav) * ((termr*geometricZ)/(termr+geometricZ))  ! meter
end subroutine  geometric2geop

end module gnssro_mod_transform

