module constants
!
! (C) Copyright 2019 UCAR
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 

!> Fortran module adopted from GSI to store constants used by GSI modules

  use kinds, only: r_kind => kind_float, i_kind => kind_int, &
                   r_single => kind_single, r_quad => kind_quad, &
                   i_long => kind_long
  implicit none

! set default as private
  private
! set passed variables to public
  public :: one,two,half,zero,deg2rad,pi,three,quarter,one_tenth
  public :: rad2deg,zero_quad,r3600,r1000,r60inv,five,four,rd_over_cp,grav
  public :: rd,rv,rozcon,rearth_equator,zero_single,tiny_r_kind,tiny_single,ten
  public :: omega,rcp,rearth,fv,h300,cp,cg_term,tpwcon,xb,ttp,psatk,xa,tmix
  public :: xai,xbi,psat,eps,omeps,wgtlim,one_quad,two_quad,epsq,climit,epsm1,hvap
  public :: hsub,cclimit,el2orc,elocp,h1000,cpr,pcpeff0,pcpeff2,delta,pcpeff1
  public :: factor1,c0,pcpeff3,factor2,dx_inv,dx_min,rhcbot,rhctop,hfus,ke2
  public :: rrow,cmr,cws,r60,huge_i_kind,huge_r_kind,t0c,rd_over_cp_mass
  public :: somigliana,grav_equator,grav_ratio,flattening,semi_major_axis
  public :: n_b,n_a,eccentricity,huge_single,constoz,g_over_rd,amsua_clw_d2
  public :: amsua_clw_d1,n_c,rd_over_g,zero_ilong
  public :: r10,r100,sqrt_tiny_r_kind,r2000,r4000
  public :: r0_01,r0_02,r0_03,r0_04,r0_05,r400,r2400
  public :: cpf_a0, cpf_a1, cpf_a2, cpf_b0, cpf_b1, cpf_c0, cpf_c1, cpf_d, cpf_e
  public :: psv_a, psv_b, psv_c, psv_d
  public :: ef_alpha, ef_beta, ef_gamma
  public :: max_varname_length
  public :: z_w_max,tfrozen
  public :: qmin,qcmin,tgmin
  public :: i_missing, r_missing

! Declare derived constants
  integer(i_kind):: huge_i_kind
  integer(i_kind), parameter :: max_varname_length=32
  real(r_single):: tiny_single, huge_single
  real(r_kind):: xai, xa, xbi, xb, dldt, rozcon,ozcon,fv, tpwcon,eps, rd_over_g
  real(r_kind):: el2orc, g_over_rd, rd_over_cp, cpr, omeps, epsm1, factor2
  real(r_kind):: factor1, huge_r_kind, tiny_r_kind, deg2rad, pi, rad2deg, cg_term
  real(r_kind):: eccentricity_linear, cv, rv, rd_over_cp_mass, cliq, rd, cp_mass
  real(r_kind):: eccentricity, grav, rearth, r60inv
  real(r_kind):: sqrt_tiny_r_kind
  real(r_kind):: n_a, n_b, n_c

! Define constants common to global and regional applications
  real(r_kind),parameter::  rearth_equator= 6.37813662e6_r_kind  ! equatorial earth radius          (m)
  real(r_kind),parameter::  omega  = 7.2921e-5_r_kind            !  angular velocity of earth       (1/s)
  real(r_kind),parameter::  cp     = 1.0046e+3_r_kind            !  specific heat of air @pressure  (J/kg/K)
  real(r_kind),parameter::  cvap   = 1.8460e+3_r_kind            !  specific heat of h2o vapor      (J/kg/K)
  real(r_kind),parameter::  csol   = 2.1060e+3_r_kind            !  specific heat of solid h2o (ice)(J/kg/K)
  real(r_kind),parameter::  hvap   = 2.5000e+6_r_kind            !  latent heat of h2o condensation (J/kg)
  real(r_kind),parameter::  hfus   = 3.3358e+5_r_kind            !  latent heat of h2o fusion       (J/kg)
  real(r_kind),parameter::  psat   = 6.1078e+2_r_kind            !  pressure at h2o triple point    (Pa)
  real(r_kind),parameter::  t0c    = 2.7315e+2_r_kind            !  temperature at zero celsius     (K)
  real(r_kind),parameter::  ttp    = 2.7316e+2_r_kind            !  temperature at h2o triple point (K)
  real(r_kind),parameter::  jcal   = 4.1855e+0_r_kind            !  joules per calorie              ()
! real(r_kind),parameter::  stndrd_atmos_ps = 1013.25e2_r_kind   ! 1976 US standard atmosphere ps   (Pa)

! Numeric constants

  integer(i_long),parameter::  zero_ilong = 0_i_long

  real(r_single),parameter::  zero_single= 0.0_r_single

  real(r_kind),parameter::  zero      = 0.0_r_kind
  real(r_kind),parameter::  r0_01     = 0.01_r_kind
  real(r_kind),parameter::  r0_02     = 0.02_r_kind
  real(r_kind),parameter::  r0_03     = 0.03_r_kind
  real(r_kind),parameter::  r0_04     = 0.04_r_kind
  real(r_kind),parameter::  r0_05     = 0.05_r_kind
  real(r_kind),parameter::  one_tenth = 0.10_r_kind
  real(r_kind),parameter::  quarter   = 0.25_r_kind
  real(r_kind),parameter::  one       = 1.0_r_kind
  real(r_kind),parameter::  two       = 2.0_r_kind
  real(r_kind),parameter::  three     = 3.0_r_kind
  real(r_kind),parameter::  four      = 4.0_r_kind
  real(r_kind),parameter::  five      = 5.0_r_kind
  real(r_kind),parameter::  ten       = 10.0_r_kind
  real(r_kind),parameter::  r10       = 10.0_r_kind
  real(r_kind),parameter::  r60       = 60._r_kind
  real(r_kind),parameter::  r100      = 100.0_r_kind
  real(r_kind),parameter::  r400      = 400.0_r_kind
  real(r_kind),parameter::  r1000     = 1000.0_r_kind
  real(r_kind),parameter::  r2000     = 2000.0_r_kind
  real(r_kind),parameter::  r2400     = 2400.0_r_kind
  real(r_kind),parameter::  r4000     = 4000.0_r_kind
  real(r_kind),parameter::  r3600     = 3600.0_r_kind

  real(r_kind),parameter:: z_w_max    = 30.0_r_kind     ! maximum diurnal thermocline thickness
  real(r_kind),parameter:: tfrozen    = 271.2_r_kind    ! sea water frozen point temperature

  real(r_quad),parameter::  zero_quad = 0.0_r_quad
  real(r_quad),parameter::  one_quad  = 1.0_r_quad
  real(r_quad),parameter::  two_quad  = 2.0_r_quad

! Constants for compressibility factor (Davis et al 1992)
  real(r_kind),parameter::  cpf_a0 =  1.58123e-6_r_kind ! K/Pa
  real(r_kind),parameter::  cpf_a1 = -2.9331e-8_r_kind  ! 1/Pa
  real(r_kind),parameter::  cpf_a2 =  1.1043e-10_r_kind ! 1/K 1/Pa
  real(r_kind),parameter::  cpf_b0 =  5.707e-6_r_kind   ! K/Pa
  real(r_kind),parameter::  cpf_b1 = -2.051e-8_r_kind   ! 1/Pa
  real(r_kind),parameter::  cpf_c0 =  1.9898e-4_r_kind  ! K/Pa
  real(r_kind),parameter::  cpf_c1 = -2.376e-6_r_kind   ! 1/Pa
  real(r_kind),parameter::  cpf_d  =  1.83e-11_r_kind   ! K2/Pa2
  real(r_kind),parameter::  cpf_e  = -0.765e-8_r_kind   ! K2/Pa2

! Constants for vapor pressure at saturation
  real(r_kind),parameter::  psv_a =  1.2378847e-5_r_kind       !  (1/K2)
  real(r_kind),parameter::  psv_b = -1.9121316e-2_r_kind       !  (1/K)
  real(r_kind),parameter::  psv_c = 33.93711047_r_kind         !
  real(r_kind),parameter::  psv_d = -6.3431645e+3_r_kind       !  (K)

! Constants for enhancement factor to calculating the mole fraction of water vapor
  real(r_kind),parameter::  ef_alpha = 1.00062_r_kind           !
  real(r_kind),parameter::  ef_beta  = 3.14e-8_r_kind           !  (1/Pa)
  real(r_kind),parameter::  ef_gamma = 5.6e-7_r_kind            !  (1/K2)

! Parameters below from WGS-84 model software inside GPS receivers.
  real(r_kind),parameter::  semi_major_axis = 6378.1370e3_r_kind     !                     (m)
  real(r_kind),parameter::  semi_minor_axis = 6356.7523142e3_r_kind  !                     (m)
  real(r_kind),parameter::  grav_polar      = 9.8321849378_r_kind    !                     (m/s2)
  real(r_kind),parameter::  grav_equator    = 9.7803253359_r_kind    !                     (m/s2) 
  real(r_kind),parameter::  earth_omega     = 7.292115e-5_r_kind     !                     (rad/s)
  real(r_kind),parameter::  grav_constant   = 3.986004418e14_r_kind  !                     (m3/s2)

! Derived geophysical constants
  real(r_kind),parameter::  flattening = (semi_major_axis-semi_minor_axis)/semi_major_axis
  real(r_kind),parameter::  somigliana = &
       (semi_minor_axis/semi_major_axis) * (grav_polar/grav_equator) - one
  real(r_kind),parameter::  grav_ratio = (earth_omega*earth_omega * &
       semi_major_axis*semi_major_axis * semi_minor_axis) / grav_constant 

! Derived thermodynamic constants
  real(r_kind),parameter::  dldti = cvap-csol
  real(r_kind),parameter::  hsub = hvap+hfus
  real(r_kind),parameter::  psatk = psat*0.001_r_kind
  real(r_kind),parameter::  tmix = ttp-20._r_kind
  real(r_kind),parameter::  elocp = hvap/cp
  real(r_kind),parameter::  rcp  = one/cp

! Constants used in GFS moist physics
  real(r_kind),parameter::  h300 = 300._r_kind
  real(r_kind),parameter::  half = 0.5_r_kind
  real(r_kind),parameter::  cclimit = 0.001_r_kind
  real(r_kind),parameter::  climit = 1.e-20_r_kind
  real(r_kind),parameter::  epsq = 2.e-12_r_kind
  real(r_kind),parameter::  h1000 = r1000
  real(r_kind),parameter::  rhcbot=0.85_r_kind
  real(r_kind),parameter::  rhctop=0.85_r_kind
  real(r_kind),parameter::  dx_max=-8.8818363_r_kind
  real(r_kind),parameter::  dx_min=-5.2574954_r_kind
  real(r_kind),parameter::  dx_inv=one/(dx_max-dx_min)
  real(r_kind),parameter::  c0=0.002_r_kind
  real(r_kind),parameter::  delta=0.6077338_r_kind
  real(r_kind),parameter::  pcpeff0=1.591_r_kind
  real(r_kind),parameter::  pcpeff1=-0.639_r_kind
  real(r_kind),parameter::  pcpeff2=0.0953_r_kind
  real(r_kind),parameter::  pcpeff3=-0.00496_r_kind
  real(r_kind),parameter::  cmr = one/0.0003_r_kind
  real(r_kind),parameter::  cws = 0.025_r_kind
  real(r_kind),parameter::  ke2 = 0.00002_r_kind
  real(r_kind),parameter::  row = r1000
  real(r_kind),parameter::  rrow = one/row
! real(r_kind),parameter::  qmin = 1.e-7_r_kind  !lower bound on ges_q

! Constant used to process ozone
  real(r_kind),parameter::  constoz = 604229.0_r_kind

! Constants used in cloud liquid water correction for AMSU-A
! brightness temperatures
  real(r_kind),parameter::  amsua_clw_d1 = 0.754_r_kind
  real(r_kind),parameter::  amsua_clw_d2 = -2.265_r_kind

! Constants used for variational qc
  real(r_kind),parameter::  wgtlim = quarter  ! Cutoff weight for concluding that obs has been
                                     ! rejected by nonlinear qc. This limit is arbitrary
                                     ! and DOES NOT affect nonlinear qc. It only affects
                                     ! the printout which "counts" the number of obs that
                                     ! "fail" nonlinear qc.  Observations counted as failing
                                     ! nonlinear qc are still assimilated.  Their weight
                                     ! relative to other observations is reduced. Changing
                                     ! wgtlim does not alter the analysis, only
                                     ! the nonlinear qc data "count"

! Minimum values for water vapor, cloud water mixing ratio, and trace gases
  real(r_kind),parameter:: qmin   = 1.e-07_r_kind   ! lower bound on ges_q
  real(r_kind),parameter:: qcmin  = 0.0_r_kind      ! lower bound on ges_cw
  real(r_kind),parameter:: tgmin  = 1.e-15_r_kind   ! lower bound on trace gases

! Constant used to detect missing input value
  integer(i_kind),parameter:: i_missing=-9999
  real(r_kind),parameter:: r_missing=-9999._r_kind

end module constants
