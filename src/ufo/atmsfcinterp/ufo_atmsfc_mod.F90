module atmsfc_mod

contains

subroutine sfc_wtq_fwd_gsi(psfc_in,tsfc,prsl1_in,tsen1,q1,u1,v1,&
                           prsl2_in,tsen2,q2,phi1,roughlen,landmask,&
                           obshgt_in,outvar,varname)
  ! sfc_wtq_fwd_gsi
  ! based off of subroutines from GSI sfc_model.f90 file
  use kinds
  use ufo_vars_mod, only: MAXVARLEN
  implicit none
  real(kind_real), intent(in) :: psfc_in, tsfc, prsl1_in, tsen1, q1, u1, v1,&
                                 prsl2_in, tsen2, q2, phi1, roughlen, landmask, &
                                 obshgt_in
  character(len=MAXVARLEN), intent(in) :: varname
  real(kind_real), intent(out) :: outvar

  real(kind_real), parameter :: zint0 = 0.01_kind_real ! default roughness over land
  real(kind_real), parameter :: k_kar = 0.4_kind_real ! Von Karman constant
  real(kind_real), parameter :: rd_over_cp = 0.28573561616_kind_real ! rd/cp
  real(kind_real), parameter :: fv = 0.60773384427_kind_real ! cv/cp - 1
  real(kind_real), parameter :: ka = 2.4e-5_kind_real
  real(kind_real), parameter :: r16 = 16.0_kind_real
  real(kind_real), parameter :: r1_1 = 1.1_kind_real
  real(kind_real), parameter :: r10 = 10.0_kind_real
  real(kind_real), parameter :: r100 = 100.0_kind_real
  real(kind_real), parameter :: r1000 = 1000.0_kind_real
  real(kind_real), parameter :: r0_9 = 0.9_kind_real
  real(kind_real), parameter :: r0_2 = 0.2_kind_real
  real(kind_real), parameter :: zero = 0.0_kind_real
  real(kind_real), parameter :: one = 1.0_kind_real
  real(kind_real), parameter :: two = 2.0_kind_real
  real(kind_real), parameter :: five = 5.0_kind_real

  real(kind_real) :: psfc, prsl1, prsl2, obshgt
  real(kind_real) :: tvg, tv1, tv2 
  real(kind_real) :: z0,zq0
  real(kind_real) :: gzzoz0, gzsoz0
  real(kind_real) :: th1, thg, thv1, thv2, thvg, eg, qg
  real(kind_real) :: wspd2, Vc2, V2
  real(kind_real) :: rib
  real(kind_real) :: psim, psimz, psih, psihz
  real(kind_real) :: cc, ust, mol, hol, holz
  real(kind_real) :: xx, yy
  real(kind_real) :: psiw, psit, psiwz, psitz, psiq, psiqz
 
  ! convert pressures to hPa from Pa
  psfc = psfc_in / r100
  prsl1 = prsl1_in / r100
  prsl2 = prsl2_in / r100 

  ! 2mTemp etc. is not at 2m in GSI output (it's 0m apparently),
  ! but 10m winds are at 10m... make height agl at least 2m
  obshgt = max(obshgt_in,two)
  print *, obshgt_in,obshgt

  ! minimum roughness length (should be in meters)
  z0 = roughlen
  if (z0 < 0.0001_kind_real) z0 = 0.0001_kind_real
  ! roughness length for over water
  if ( landmask < 0.01 ) then
     zq0 = z0
  else
     zq0 = zint0
  end if

  ! constant variable for psi
  gzsoz0 = log(phi1/z0)
  gzzoz0 = log(obshgt/z0)

  ! virtual temperature from sensible temperature
  tv1 = tsen1 * (one + fv * q1)
  tv2 = tsen2 * (one + fv * q2)

  ! convert temperature of the ground to virtual temp assuming saturation
  call da_tp_to_qs( tsfc, psfc, eg, qg)
  tvg = tsfc * (one + fv * qg)

  ! potential temperature calculations
  thg = tsfc * (r1000 / psfc) ** rd_over_cp ! surface theta
  th1 = tsen1 * (r1000 / prsl1) ** rd_over_cp ! theta for lowest model layer

  ! virtual potential temperature
  thv1 = tv1 * (r1000 / prsl1) ** rd_over_cp ! surface theta
  thv2 = tv2 * (r1000 / prsl2) ** rd_over_cp ! surface theta
  thvg = tvg * (r1000 / psfc) ** rd_over_cp ! surface theta

  ! wind speed
  wspd2 = u1*u1 + v1*v1  

  ! convective velocity
  if (thvg >= thv1) then
    Vc2 = 4.0_kind_real * (thvg - thv1)
  else
    Vc2 = zero
  end if

  V2 = 0.000001_kind_real + wspd2 + Vc2

  ! bulk richardson number
  rib = (9.80665_kind_real * phi1 / th1) * (thv1 - thvg) / V2

  ! calculate psi based off of regime
  ! stable conditions
  if (rib >= r0_2) then
    psim = -r10*gzsoz0 
    psimz = -r10*gzzoz0
    psim = max(psim,-r10)
    psimz = max(psimz,-r10)
    psih = psim
    psihz = psimz

  ! mechanically driven turbulence
  else if ((rib < r0_2) .and. (rib > zero)) then
    psim = ( -five * rib) * gzsoz0 / (r1_1 - five*rib)  
    psimz = ( -five * rib) * gzzoz0 / (r1_1 - five*rib)  
    psim = max(psim,-r10)
    psimz = max(psimz,-r10)
    psih = psim
    psihz = psimz

  ! unstable forced convection
  else if ((rib == zero) .or. (rib < zero .and. thv2>thv1)) then
    psim = zero
    psimz = zero
    psih = psim
    psihz = psimz

  ! free convection
  else
    psim = zero
    psih = zero
    cc = two * atan(one)
    
    ! friction speed
    ust = k_kar * sqrt(V2) / (gzsoz0 - psim)
    ! heat flux factor
    mol = k_kar * (th1 - thg)/(gzsoz0 - psih)
    ! ratio of PBL height to Monin-Obukhov length
    if (ust < 0.01_kind_real) then
      hol = rib * gzsoz0
    else
      hol = k_kar * 9.80665_kind_real * phi1 * mol / (th1 * ust * ust)
    end if
    hol = min(hol,zero)
    hol = max(hol,-r10)
    holz = (obshgt / phi1) * hol 
    holz = min(holz,zero)
    holz = max(holz,-r10)

    xx = (one - r16 * hol) ** 0.25_kind_real 
    yy = log((one+xx*xx)/two) 
    psim = two * log((one+xx)/two) + yy - two * atan(xx) + cc
    psih = two * yy

    xx = (one - r16 * holz) ** 0.25_kind_real
    yy = log((one+xx*xx)/two) 
    psimz = two * log((one+xx)/two) + yy - two * atan(xx) + cc
    psihz = two * yy

    psim = min(psim,r0_9*gzsoz0)
    psimz = min(psimz, r0_9*gzzoz0)
    psih = min(psih,r0_9*gzsoz0)
    psihz = min(psihz,r0_9*gzzoz0)

  end if
  
  psiw = gzsoz0 - psim
  psit = gzsoz0 - psih 
  psiwz = gzzoz0 - psimz
  psitz = gzzoz0 - psihz

  ust = k_kar * sqrt(V2) / (gzsoz0 - psim)

  psiq = log(k_kar*ust*phi1/ka + phi1 / zq0) - psih
  psiqz = log(k_kar*ust*obshgt/ka + obshgt / zq0) - psihz

  select case(trim(varname))
    case("air_temperature")
      outvar = (thg + (th1 - thg)*psitz/psit)*(psfc/r1000)**rd_over_cp
    case("virtual_temperature")
      outvar = (thg + (th1 - thg)*psitz/psit)*(psfc/r1000)**rd_over_cp
      outvar = outvar * (one + fv * q1)  
    case("specific_humidity")
      outvar = qg + (q1 - qg)*psiqz/psiq
    case("eastward_wind")
      outvar = u1 * psiwz / psiw 
    case("northward_wind")
      outvar = v1 * psiwz / psiw
  end select
  return

end subroutine

subroutine DA_TP_To_Qs( t, p, es, qs )

!$$$ subprogram documentation block
!               .      .    .                                       .
! subprogram:   DA_TP_To_Qs
!
!   prgrmmr:
!
! abstract:  Convert T/p to saturation specific humidity.
!
!  METHOD: qs = es_alpha * es / ( p - ( 1 - rd_over_rv ) * es ).
!          Use Rogers & Yau (1989) formula: es = a exp( bTc / (T_c + c) ).
!
! program history log: 
!   2000-10-03  Barker  - Creation of F90 version.
!   2008-04-14  safford - added standard documentation block
!
!   input argument list:
!     t          - Temperature.
!     p          - Pressure.
!
!   output argument list:
!     es         - Sat. vapour pressure.
!     qs         - Sat. specific humidity.
!
! attributes:
!   language:  f90
!   machine:   ibm RS/6000 SP
!
!$$$ end documentation block
   use kinds

   implicit none

   real(kind_real), intent(in   ) :: t                ! Temperature.
   real(kind_real), intent(in   ) :: p                ! Pressure.
   real(kind_real), intent(  out) :: es               ! Sat. vapour pressure.
   real(kind_real), intent(  out) :: qs               ! Sat. specific humidity.

!  Saturation Vapour Pressure Constants(Rogers & Yau, 1989)
   real(kind_real), parameter    :: es_alpha = 611.2_kind_real
   real(kind_real), parameter    :: es_beta = 17.67_kind_real
   real(kind_real), parameter    :: es_gamma = 243.5_kind_real
  
   real(kind_real), parameter    :: eps = 0.62199349945_kind_real ! rd/rv 
   real(kind_real), parameter    :: t0c = 2.7315e+2_kind_real ! 0 deg C in K
   real(kind_real) :: omeps

   real(kind_real)                          :: t_c              ! T in degreesC.


!------------------------------------------------------------------------------
!  [1.0] Initialise:
!------------------------------------------------------------------------------
   omeps = 1.0_kind_real - eps
   t_c = t - t0c

!------------------------------------------------------------------------------
!  [2.0] Calculate saturation vapour pressure:
!------------------------------------------------------------------------------

   es = 0.01_kind_real * es_alpha * exp( es_beta * t_c / ( t_c + es_gamma ) ) 

!------------------------------------------------------------------------------
!  [3.0] Calculate saturation specific humidity:
!------------------------------------------------------------------------------

   qs = eps * es / ( p - omeps * es )


return
end subroutine DA_TP_To_Qs


end module atmsfc_mod
