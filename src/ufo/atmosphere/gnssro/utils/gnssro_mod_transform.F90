!==========================================================================
module gnssro_mod_transform 
!==========================================================================

use kinds
use gnssro_mod_constants

real(kind_real), parameter ::  semi_major_axis = 6378.1370e3_kind_real     !                     (m)
real(kind_real), parameter ::  semi_minor_axis = 6356.7523142e3_kind_real  !                     (m)
real(kind_real), parameter ::  grav_polar      = 9.8321849378_kind_real    !                     (m/s2)
real(kind_real), parameter ::  grav_equator    = 9.7803253359_kind_real    !                     (m/s2) 
real(kind_real), parameter ::  earth_omega     = 7.292115e-5_kind_real     !                     (rad/s)
real(kind_real), parameter ::  grav_constant   = 3.986004418e14_kind_real  !      
real(kind_real), parameter ::  flattening  = (semi_major_axis-semi_minor_axis)/semi_major_axis
real(kind_real), parameter ::  somigliana  = (semi_minor_axis/semi_major_axis) * (grav_polar/grav_equator) - one
real(kind_real), parameter ::  grav_ratio  = (earth_omega*earth_omega * &
                                              semi_major_axis*semi_major_axis * semi_minor_axis) / grav_constant
real(kind_real), parameter ::  eccentricity = sqrt(semi_major_axis**2 - semi_minor_axis**2)/semi_major_axis

contains

! ------------------------------
! variable converting between geopotential and geometric heights using  MJ Mahoney's (2001)
! Parameters from WGS-84 model software inside GPS receivers.
! copy from GSI
subroutine  geometric2geop(Latitude,geometricZ, geopotentialH ) 

real(kind_real), intent(in)  :: Latitude,   geometricZ
real(kind_real), intent(out) :: geopotentialH
real(kind_real)              :: sino, termg, termr ! local variables

sino          = sin(deg2rad*Latitude)
termg         = grav_equator*( (one+somigliana*sino)/sqrt(one-eccentricity*eccentricity*sino) )
termr         = semi_major_axis / (one + flattening + grav_ratio - two*flattening*sino)
geopotentialH = (termg/grav) * ((termr*geometricZ)/(termr+geometricZ))  ! meter

end subroutine  geometric2geop


!subroutine geop2geometric(Latitude, geopotentialH, geometricZ)

!real(kind_real),intent(in)   :: Latitude,   geopotentialH
!real(kind_real),intent(out)  :: geometricZ
!real(kind_real)              :: sino, termg, termr ! local variables

!sino          = sin(deg2rad*Latitude)
!termg         = grav_equator*( (one+somigliana*sino)/sqrt(one-eccentricity*eccentricity*sino) )
!termr         = semi_major_axis / (one + flattening + grav_ratio - two*flattening*sino)
!geometricZ    = (termr*geopotentialH)/((termg/grav)*termr-geopotentialH)
!
!end subroutine geop2geometric


  subroutine geop2geometric(latitude, geopotentialH, geometricZ, gp2gm)
  ! calculate observation geometric height using  MJ Mahoney's (2001), eq(23)
    real(kind_real),intent(in)  :: latitude,   geopotentialH
    real(kind_real),intent(out)  :: geometricZ
    real(kind_real),intent(out)  :: gp2gm !jocabian 
    real(kind_real)            :: sino
    real(kind_real):: termg, termr, termrg

    sino          = sin(deg2rad*latitude)
    termg         = grav_equator*( (one+somigliana*sino)/sqrt(one-eccentricity*eccentricity*sino) )
    termr         = semi_major_axis / (one + flattening + grav_ratio - two*flattening*sino)
    termrg        = termg/grav*termr

    gp2gm = termr/(termrg-geopotentialH) + (termr*geopotentialH)/(termrg-geopotentialH)**2

    geometricZ = (termr*geopotentialH)/((termg/grav)*termr-geopotentialH)

  end subroutine geop2geometric
! ------------------------------

  subroutine compute_refractivity(temperature, specH, pressure,refr, use_compress)

    real(kind_real), intent(in) :: temperature, specH, pressure
    real(kind_real), intent(out) :: refr
    real(kind_real) :: fact,pw,refr1,refr2,refr3, tfact
    logical      ,intent(in)  :: use_compress       

   ! constants needed to compute refractivity
    call gnssro_ref_constants(use_compress)

    fact  = one+fv*specH
    tfact = (1-rd_over_rv)*specH+rd_over_rv
    pw    = rd_over_rv+specH*(one-rd_over_rv)
    refr1 = n_a*pressure/temperature
    refr2 = n_b*specH*pressure/(temperature**2*tfact)
    refr3 = n_c*specH*pressure/(temperature*tfact)
    refr  = refr1 + refr2 + refr3

  end subroutine compute_refractivity


  subroutine compute_refractivity_Tv(virT, specH, pressure,refr, use_compress)
    
    real(kind_real), intent(in) :: virT, specH, pressure 
    real(kind_real), intent(out) :: refr
    real(kind_real) :: fact,pw,refr1,refr2,refr3
    logical      ,intent(in)  :: use_compress       ! use computed compressibility or value as 1

    ! constants needed to compute refractivity
    call gnssro_ref_constants(use_compress)

    fact  = one+fv*specH
    pw    = rd_over_rv+specH*(one-rd_over_rv)
    refr1 = n_a*(pressure/virT)*fact
    refr2 = n_b*specH*pressure*fact**2/(virT**2*pw)
    refr3 = n_c*fact*specH*pressure/(virT*pw)
    refr  = refr1 + refr2 + refr3

  end subroutine compute_refractivity_Tv
end module gnssro_mod_transform

