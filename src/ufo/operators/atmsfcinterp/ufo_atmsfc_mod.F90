module atmsfc_mod

contains

subroutine sfc_wind_fact_gsi(u, v, tsen, q, psfc, prsi1, prsi2,&
                             skint, z0, lsmask, f10m) 
  ! compute wind reduction factor
  ! aka the coefficient to multiply u,v in lowest model level
  ! to get u,v at 'obshgt' aka 10m agl
  ! based off of compute_fact10 from GSI
  use kinds
  use ufo_constants_mod, only: rd_over_cp, rv_over_rd, rd, grav, &
                               zero, quarter, half, one, two, four, ten
  implicit none
  real(kind_real), intent(in) :: u, v, tsen, q, psfc, prsi1, prsi2,&
                                 skint, z0, lsmask 
  real(kind_real), intent(out) :: f10m
  real(kind_real), parameter :: alpha = 5.0_kind_real
  real(kind_real), parameter :: a0 = -3.975_kind_real
  real(kind_real), parameter :: a1 = 12.32_kind_real
  real(kind_real), parameter :: b1 = -7.755_kind_real
  real(kind_real), parameter :: b2 = 6.041_kind_real
  real(kind_real), parameter :: a0p = -7.941_kind_real
  real(kind_real), parameter :: a1p = 24.75_kind_real
  real(kind_real), parameter :: b1p = -8.705_kind_real
  real(kind_real), parameter :: b2p = 7.899_kind_real
  real(kind_real), parameter :: vis = 1.4e-5_kind_real
  real(kind_real), parameter :: fv = rv_over_rd - one
  real(kind_real),parameter::  charnok = 0.014_kind_real

  real(kind_real) :: rat, restar, ustar, del, tem, prki1, prki2, prkl, prsl
  real(kind_real) :: wspd, wind, q0, theta1, tv1, thv1, tvs, z1
  real(kind_real) :: z0max, ztmax, dtv, adtv, rb, fm, fh, hlinf, fm10
  real(kind_real) :: hl1, hl0inf, hltinf, aa, aa0, bb, bb0, pm, ph, fms, fhs
  real(kind_real) :: hl0, hlt, hl110, pm10, hl12, ph2, olinf 

   
  rat=zero
  restar=zero
  ustar=zero
  del = (prsi1-prsi2)*1.0e-3_kind_real
  tem = (one + rd_over_cp) * del
  prki1 = (prsi1*1.0e-5_kind_real)**rd_over_cp
  prki2 = (prsi2*1.0e-5_kind_real)**rd_over_cp
  prkl = (prki1*(prsi1*1.0e-3_kind_real)-prki2*(prsi2*1.0e-3_kind_real))/tem
  prsl = 1.0e5_kind_real*prkl**(one/rd_over_cp)
  wspd = sqrt( u*u + v*v )
  wind = max(wspd,one)
  q0 = max(q,1.e-8_kind_real)
  theta1 = tsen * (prki1/prkl)
  tv1 = tsen * (one+fv*q0)
  thv1 = theta1 * (one+fv*q0)
  tvs = skint * (one+fv*q0) ! fix this
  z1 = -rd*tv1*log(prsl/psfc)/grav

  ! compute stability dependent exchange coefficients
  if (lsmask < 0.01_kind_real) then
    ustar = sqrt(grav * z0 / charnok)
  end if 
  ! compute stability indices (rb and hlinf)
  z0max = min(z0,one*z1)
  ztmax = z0max
  if (lsmask < 0.01_kind_real) then
    restar = ustar * z0max / vis
    restar = max(restar, 1.e-6_kind_real)
    ! rat taken from Zeng, Zhao and Dickinson 1997
    rat = 2.67_kind_real * restar**quarter - 2.57_kind_real
    rat = min(rat, 7.0_kind_real)
    ztmax = z0max * exp(-rat)
  end if

  dtv = thv1 - tvs
  adtv = abs(dtv)
  adtv = max(adtv,1.e-3_kind_real)
  dtv = sign(one,dtv)*adtv
  rb = grav * dtv * z1 / (half * (thv1 + tvs) * wind * wind)
  rb = max(rb,-5.e3_kind_real)
  fm = log((z0max + z1) / z0max)
  fh = log((ztmax + z1) / ztmax)
  hlinf = rb * fm * fm / fh
  fm10 = log((z0max + 10.0_kind_real) / z0max)
    
  ! stable case
  if (dtv >= zero) then
    hl1 = hlinf
  end if
  if ((dtv >= zero) .and. (hlinf > quarter)) then
    hl0inf = z0max * hlinf / z1
    hltinf = ztmax * hlinf / z1
    aa = sqrt(one + four*alpha*hlinf)
    aa0 = sqrt(one + four*alpha*hl0inf)
    bb = aa
    bb0 = sqrt(one + four*alpha*hltinf)
    pm = aa0 - aa + log((aa+one)/(aa0+one))
    ph = bb0 - bb + log((bb+one)/(bb0+one))
    fms = fm - pm
    fhs = fh - ph
    hl1 = fms * fms * rb / fhs
  end if
  ! second iteration
  if (dtv >= zero) then
    hl0 = z0max * hl1 / z1
    hlt = ztmax * hl1 / z1
    aa = sqrt(one + four*alpha*hl1)
    aa0 = sqrt(one + four*alpha*hl0)
    bb = aa
    bb0 = sqrt(one + four*alpha*hlt)
    pm = aa0 - aa + log((aa+one)/(aa0+one))
    ph = bb0 - bb + log((bb+one)/(bb0+one))
    hl110 = hl1 * ten / z1
    aa = sqrt(one + four*alpha*hl110)
    pm10 = aa0 - aa + log((aa+one)/(aa0+one))
    hl12 = hl1 * two / z1
    bb = sqrt(one + four * alpha * hl12)
    ph2 = bb0 - bb + log((bb+one)/(bb0+one))
  end if

  ! unstable case
  ! check for unphysical obukhov length
  if (dtv < zero) then
    olinf = z1 / hlinf
    if ( abs(olinf) <= z0max * 50.0_kind_real ) then
      hlinf = -z1 / (50.0_kind_real * z0max)
    end if
  end if

  ! get pm and ph
  if (dtv < zero .and. hlinf >= (-half)) then
    hl1 = hlinf
    pm = (a0 + a1*hl1) * hl1 / (one + b1*hl1 + b2*hl1*hl1)
    ph = (a0p + a1p*hl1)*hl1/(one + b1p*hl1 + b2*hl1*hl1)
    hl110 = hl1 * ten / z1
    pm10 = (a0 + a1*hl110)*hl110/(one + b1*hl110 + b2*hl110*hl110)
    hl12 = hl1 * two / z1
    ph2 = (a0p + a1p*hl12)*hl12/(one + b1p*hl12 + b2p*hl12*hl12)
  end if
  if (dtv < zero .and. hlinf < (-half)) then
    hl1 = -hlinf
    pm = log(hl1) + two*hl1**(-quarter) - 0.8776_kind_real
    ph = log(hl1) + half*hl1**(-half) + 1.386_kind_real
    hl110 = hl1 * ten / z1
    pm10 = log(hl110) + two*hl110**(-quarter) - 0.8776_kind_real
    hl12 = hl1 * two / z1
    ph2 = log(hl12) + half*hl12**(-half) + 1.386_kind_real
  end if

  ! finish the exchange coefficient computation to provide fm and fh
  fm = fm - pm
  fh = fh - ph
  fm10 = fm10 - pm10
  f10m = max(zero, min(fm10/fm, one))

  return

end subroutine sfc_wind_fact_gsi

!--------------------------------------------------------------------------

subroutine calc_psi_vars_gsi(rib, gzsoz0, gzzoz0, thv1, thv2,&
                             V2, th1, thg, phi1, obshgt, & 
                             psim, psih, psimz, psihz)
  ! calculate psi based off of near-surface atmospheric regime
  use kinds
  use ufo_constants_mod, only: grav, von_karman, zero, quarter, &
                               one, two, five, ten
  implicit none
  real(kind_real), intent(in) :: rib, gzsoz0, gzzoz0, thv1, thv2, &
                                 V2, th1, thg, phi1, obshgt
  real(kind_real), intent(out) :: psim, psih, psimz, psihz 

  real(kind_real), parameter :: r0_2 = 0.2_kind_real
  real(kind_real), parameter :: r1_1 = 1.1_kind_real
  real(kind_real), parameter :: r0_9 = 0.9_kind_real
  real(kind_real), parameter :: r16 = 16.0_kind_real
  real(kind_real) :: cc, mol, ust, hol, holz, xx, yy 

  ! stable conditions
  if (rib >= r0_2) then
    psim = -ten*gzsoz0 
    psimz = -ten*gzzoz0
    psim = max(psim,-ten)
    psimz = max(psimz,-ten)
    psih = psim
    psihz = psimz
  ! mechanically driven turbulence
  else if ((rib < r0_2) .and. (rib > zero)) then
    psim = ( -five * rib) * gzsoz0 / (r1_1 - five*rib)  
    psimz = ( -five * rib) * gzzoz0 / (r1_1 - five*rib)  
    psim = max(psim,-ten)
    psimz = max(psimz,-ten)
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
    ust = von_karman * sqrt(V2) / (gzsoz0 - psim)
    ! heat flux factor
    mol = von_karman * (th1 - thg)/(gzsoz0 - psih)
    ! ratio of PBL height to Monin-Obukhov length
    if (ust < 0.01_kind_real) then
      hol = rib * gzsoz0
    else
      hol = von_karman * grav * phi1 * mol / (th1 * ust * ust)
    end if
    hol = min(hol,zero)
    hol = max(hol,-ten)
    holz = (obshgt / phi1) * hol 
    holz = min(holz,zero)
    holz = max(holz,-ten)

    xx = (one - r16 * hol) ** quarter 
    yy = log((one+xx*xx)/two) 
    psim = two * log((one+xx)/two) + yy - two * atan(xx) + cc
    psih = two * yy

    xx = (one - r16 * holz) ** quarter
    yy = log((one+xx*xx)/two) 
    psimz = two * log((one+xx)/two) + yy - two * atan(xx) + cc
    psihz = two * yy

    psim = min(psim,r0_9*gzsoz0)
    psimz = min(psimz, r0_9*gzzoz0)
    psih = min(psih,r0_9*gzsoz0)
    psihz = min(psihz,r0_9*gzzoz0)

  end if


end subroutine calc_psi_vars_gsi

!--------------------------------------------------------------------------

subroutine calc_conv_vel_gsi(u1, v1, thvg, thv1, V2)
  ! compute convective velocity for use in computing psi vars
  use kinds
  implicit none
  real(kind_real), intent(in) :: u1, v1, thvg, thv1
  real(kind_real), intent(out) :: V2
  real(kind_real) :: wspd2, Vc2

  wspd2 = u1*u1 + v1*v1
  if (thvg >= thv1) then
    Vc2 = 4.0_kind_real * (thvg - thv1)
  else
    Vc2 = 0.0_kind_real
  end if

  V2 = 1e-6_kind_real + wspd2 + Vc2

end subroutine calc_conv_vel_gsi

!--------------------------------------------------------------------------

end module atmsfc_mod
