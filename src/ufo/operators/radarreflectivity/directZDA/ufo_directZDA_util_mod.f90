! (C) Copyright 2017-2018 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

!> Fortran module for directZDA observation operator util

module ufo_directZDA_util_mod

use kinds

 implicit none
 private 

 public :: coef4dbzfwrd
 public :: calc_coeffs_dry_snow_tm
 public :: Cr,     Pr
 public :: Cs_dry, Ps_dry, Cs_wet, Ps_wet
 public :: Cg_dry, Pg_dry, Cg_wet, Pg_wet

! coefficients and raise-to-power inder numebrs used in
! reflectivity obs forward operator equations (for single moment MP scheme)
 real(kind_real)  :: Cr        ! coefficient for rain drop
 real(kind_real)  :: Cs_dry    !             for dry snow
 real(kind_real)  :: Cs_wet    !             for wet snow (melting)
 real(kind_real)  :: Cg_dry    !             for drygraupel/hail
 real(kind_real)  :: Cg_wet    !             for drygraupel/hail
 real(kind_real)  :: Pr        ! power index for rain drop
 real(kind_real)  :: Ps_dry    !             for dry snow
 real(kind_real)  :: Ps_wet    !             for wet snow (melting)
 real(kind_real)  :: Pg_dry    !             for drygraupel/hail
 real(kind_real)  :: Pg_wet    !             for drygraupel/hail

contains

subroutine coef4dbzfwrd(mphyopt,iret)
!----------------------------------------------------------------------!
!
! purpose:
! this subroutine is based on the program to calcualte coefficients for
! reflecivity formula used in arps3dvar and arpsenkf (Mingjing Tong's version)
! ref. Smith et al. 1975
! It is modified to work with WSM6 considering the density of graupel, instead
! of the density of hail in Lin scheme.
!
! prgmmr:
!   2018-02-19  g.zhao  initialization of the code
!----------------------------------------------------------------------!

      implicit none

      integer :: mphyopt, iret
!-----------------------------------------------------------------------
! Declare local parameters.
!-----------------------------------------------------------------------

      real(kind_real), parameter :: ki2 = 0.176_kind_real ! Dielectric factor for ice if other
                                                          ! than melted drop diameters are used.
      real(kind_real), parameter :: kw2=0.93_kind_real    ! Dielectric factor for water.

      real(kind_real), parameter :: degKtoC=273.15_kind_real ! Conversion factor from degrees K to
                                                             !   degrees C

      real(kind_real), parameter :: m3todBZ=1.0E+18_kind_real ! Conversion factor from m**3 to
                                                              !   mm**6 m**-3.

      real(kind_real), parameter :: pi=3.1415926_kind_real ! Pi.

      real(kind_real), parameter :: pipowf=7.0_kind_real/4.0_kind_real ! Power to which pi is raised.

      real(kind_real), parameter :: N0r_0=8.0E+06_kind_real ! Intercept parameter in 1/(m^4) for rain.
      real(kind_real), parameter :: N0s_0=3.0E+06_kind_real ! Intercept parameter in 1/(m^4) for snow.
      real(kind_real), parameter :: N0g_0=4.0E+05_kind_real ! Intercept parameter in 1/(m^4) for graupel. (<=4.0E+6)
      real(kind_real), parameter :: N0h_0=4.0E+04_kind_real ! Intercept parameter in 1/(m^4) for hail.

      real(kind_real), parameter :: N0xpowf=3.0_kind_real/4.0_kind_real ! Power to which N0r,N0s & N0h are
                                                                        ! raised.

      real(kind_real), parameter :: approxpow_0=0.95_kind_real ! Approximation power for hail/graupel
                                                               ! integral.

      real(kind_real), parameter :: rqrpowf=7.0_kind_real/4.0_kind_real  ! Power to which product rho * qr
                                                                         ! is raised.
      real(kind_real), parameter :: rqsnpowf=7.0_kind_real/4.0_kind_real ! Power to which product rho * qs
                                                                         ! is raised (dry snow).
      real(kind_real), parameter :: rqsppowf=7.0_kind_real/4.0_kind_real ! Power to which product rho * qs
                                                                         ! is raised (wet snow).
      real(kind_real), parameter :: rqhnpowf=7.0_kind_real/4.0_kind_real ! Power to which product rho * qh
                                                                         ! is raised (dry hail)

      real(kind_real), parameter :: rhoi_0=917._kind_real  ! Density of ice (kg m**-3)
      real(kind_real), parameter :: rhor_0=1000._kind_real ! Density of rain (kg m**-3)
      real(kind_real), parameter :: rhos_0=100._kind_real  ! Density of snow (kg m**-3)
      real(kind_real), parameter :: rhoh_0=913._kind_real  ! Density of hail (kg m**-3)
      real(kind_real), parameter :: rhog_0=400._kind_real  ! Density of graupel (kg m**-3)

      real(kind_real), parameter :: rhoipowf=2.0_kind_real ! Power to which rhoi is raised.
      real(kind_real), parameter :: rhospowf=1.0_kind_real/4.0_kind_real ! Power to which rhos is raised.
      real(kind_real), parameter :: rhoxpowf=7.0_kind_real/4.0_kind_real ! Power to which rhoh is raised.

      real(kind_real), parameter :: Zefact=720.0_kind_real ! Multiplier for Ze components.

      real(kind_real), parameter :: lg10mul=10.0_kind_real ! Log10 multiplier

! ------------------------------------------------------------------------------

!-----------------------------------------------------------------------
! Declare local variables.
!-----------------------------------------------------------------------
      real(kind_real)            :: rqhppowf     ! Power to which product rho * qh is raised (wet hail).

      real(kind_real)            :: rcomp,scomp,hcomp,sumcomp
      integer(kind_int)          :: i,j,k
      real(kind_real)            :: Zerf,Zesnegf,Zesposf,Zegf, Zehf
      real(kind_real)            :: Zehnegf,Zehposf

      integer(kind_int)          :: ios
      real(kind_real)            :: N0r, N0s, N0g, N0h, N0i
      real(kind_real)            :: rhor, rhos, rhog, rhoh, rhoi
      real(kind_real)            :: approxpow

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
! Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!-----------------------------------------------------------------------
! Initialization of Namelist variables
      N0r=N0r_0
      N0s=N0s_0
      N0g=N0g_0
      N0h=N0h_0
      N0i=N0r_0                 ! N_not for ice is not used for now
      rhor=rhor_0
      rhos=rhos_0
      rhog=rhog_0
      rhoi=rhoi_0
      rhoh=rhoh_0
      approxpow=approxpow_0
!-----------------------------------------------------------------------
      rqhppowf=(7.0_kind_real/4.0_kind_real)*approxpow      ! Power to which product rho * qh is raised (wet hail).
!-----------------------------------------------------------------------
! First gather all the constants together.  (These are treated as
! variables because Fortran 90 does not allow real exponents when
! calculating parameters).
!-----------------------------------------------------------------------

      Zerf = (m3todBZ * Zefact) /  &
          ((pi ** pipowf) * (N0r ** N0xpowf) *  &
          (rhor ** rhoxpowf))
      Zesnegf = ((m3todBZ * Zefact   * Ki2 * (rhos ** rhospowf)) /  &
          ((pi ** pipowf) * Kw2 * (N0s ** N0xpowf) *  &
          (rhoi ** rhoipowf)))
      Zesposf = ((m3todBZ * Zefact) /  &
          ((pi ** pipowf) * (N0s ** N0xpowf) *  &
          (rhos ** rhoxpowf)))

      select case (mphyopt)
      case (2,3,4)
          ! hail 
          ! changed the density of hail to graupel to make it consistent with WRF Lin scheme
           Zehnegf = ((m3todBZ * Zefact) * Ki2  /  &
               ((pi ** pipowf) * Kw2 * (4.0E+06_kind_real ** N0xpowf) *  &
               (400._kind_real ** rhoxpowf)))
           Zehposf = (((m3todBZ * Zefact) /  &
               ((pi ** pipowf) * (4.0E+06_kind_real ** N0xpowf) *  &
               (400._kind_real ** rhoxpowf))) ** approxpow)
           iret = 0
      case (5,6,7)
          ! graupel
          Zehnegf = ((m3todBZ * Zefact) * Ki2  /  &
              ((pi ** pipowf) * Kw2 * (N0g ** N0xpowf) *  &
              (rhog ** rhoxpowf)))
          Zehposf = (((m3todBZ * Zefact) /  &
              ((pi ** pipowf) * (N0g ** N0xpowf) *  &
              (rhog ** rhoxpowf))) ** approxpow)
          iret = 0
      case default
          write(6,*) ' subroutine COEF4DBZFWRD: warning --> invalid mphyopt for single moment scheme'
          iret = -1
      end select

!     coefficients used in dbz formula
      Cr        = Zerf
      Cs_dry    = Zesnegf                            ! <= 273.15K
      Cs_wet    = Zesposf                            ! >  273.15K (melting)
      Cg_dry    = Zehnegf                            ! <= 273.15K
      Cg_wet    = Zehposf                            ! >  273.15K (melting)
!     raise-to-power index numbers used in dbz formula
      Pr        = rqrpowf
      Ps_dry    = rqsnpowf                           ! <= 273.15K
      Ps_wet    = rqsppowf                           ! >  273.15K (melting)
      Pg_dry    = rqhnpowf                           ! <= 273.15K
      Pg_wet    = rqhppowf                           ! >  273.15K (melting)

      WRITE(*,*)'*****************************************************************'
      write(*,*)'COEF4DBZFWRD: mphyopt==',mphyopt
      !rite(*,*)'COEF4DBZFWRD: rain:    ',' Cr    =',Cr,    '  Pr    =',Pr
      !rite(*,*)'COEF4DBZFWRD: snow:    ',' Cs_dry=',Cs_dry,'  Ps_dry=',Ps_dry
      !rite(*,*)'COEF4DBZFWRD:          ',' Cs_wet=',Cs_wet,'  Ps_wet=',Ps_wet
      !rite(*,*)'COEF4DBZFWRD: graupel: ',' Cg_dry=',Cg_dry,'  Pg_dry=',Pg_dry
      !rite(*,*)'COEF4DBZFWRD:          ',' Cg_wet=',Cg_wet,'  Pg_wet=',Pg_wet
      WRITE(*,*)'*****************************************************************'

      return

  end subroutine coef4dbzfwrd

  subroutine calc_coeffs_dry_snow_tm(temp,a_out,b_out)
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:    calc_coeffs_dry_snow_tm - calculate coefficients for dry snow
!                                          dependent on temperature 
! program history log:
!   2021-08-30  CAPS (J. Park) - Initial code; based on WRF v4.0
!                              - (from module_mp_thompson.F)
!                              - added kind attribute and simpified
!   input argument list:
!     temp     - air temperature (K)
!
!   output argument list:
!     a_out - coefficient for gamma distribution of dry snow
!     b_out - coefficient for gamma distribution of dry snow
!
    implicit none

    !..For snow moments conversions (from Field et al. 2005)
    real(kind_real), intent(in   ) :: temp
    real(kind_real), intent(  out) :: a_out, b_out
    ! LOCAL
    real(kind_real) :: a_, b_, loga_, tc0
    real(kind_real), dimension(10), parameter:: &
    sa = (/ 5.065339_kind_real, -0.062659_kind_real, -3.032362_kind_real, &
            0.029469_kind_real, -0.000285_kind_real, &
            0.31255_kind_real,   0.000204_kind_real,  0.003199_kind_real, &
            0.0_kind_real,      -0.015952_kind_real/)
    real(kind_real), dimension(10), parameter:: &
    sb = (/ 0.476221_kind_real, -0.015896_kind_real,  0.165977_kind_real, &
            0.007468_kind_real, -0.000141_kind_real, &
            0.060366_kind_real,  0.000079_kind_real,  0.000594_kind_real, &
            0.0_kind_real,      -0.003577_kind_real/)
    real(kind_real), parameter:: mu_s = 0.6357_kind_real
    real(kind_real), parameter:: bm_s = 2.0_kind_real
    real(kind_real), parameter:: bv_s = 0.55_kind_real

    real(kind_real), dimension(18):: cse

    cse(1) = bm_s + 1._kind_real
    cse(2) = bm_s + 2._kind_real
    cse(3) = bm_s*2._kind_real
    cse(4) = bm_s + bv_s + 1._kind_real
    cse(5) = bm_s*2._kind_real + bv_s + 1._kind_real
    cse(6) = bm_s*2._kind_real + 1._kind_real
    cse(7) = bm_s + mu_s + 1._kind_real
    cse(8) = bm_s + mu_s + 2._kind_real
    cse(9) = bm_s + mu_s + 3._kind_real
    cse(10) = bm_s + mu_s + bv_s + 1._kind_real
    cse(11) = bm_s*2._kind_real + mu_s + bv_s + 1._kind_real
    cse(12) = bm_s*2._kind_real + mu_s + 1._kind_real
    cse(13) = bv_s + 2._kind_real
    cse(14) = bm_s + bv_s
    cse(15) = mu_s + 1._kind_real
    cse(16) = 1.0_kind_real + (1.0_kind_real + bv_s)/2._kind_real
    cse(17) = cse(16) + mu_s + 1._kind_real
    cse(18) = bv_s + mu_s + 3._kind_real

!..Sum of two gamma distrib for snow (Field et al. 2005).
    tc0 = MIN(-0.1_kind_real, temp-273.15_kind_real)

    loga_ = sa(1) + sa(2)*tc0 + sa(3)*cse(3) &
     &     + sa(4)*tc0*cse(3) + sa(5)*tc0*tc0 &
     &     + sa(6)*cse(3)*cse(3) + sa(7)*tc0*tc0*cse(3) &
     &     + sa(8)*tc0*cse(3)*cse(3) + sa(9)*tc0*tc0*tc0 &
     &     + sa(10)*cse(3)*cse(3)*cse(3)
    a_ = 10.0_kind_real**loga_
    b_ = sb(1)+ sb(2)*tc0 + sb(3)*cse(3) + sb(4)*tc0*cse(3) &
     &    + sb(5)*tc0*tc0 + sb(6)*cse(3)*cse(3) &
     &    + sb(7)*tc0*tc0*cse(3) + sb(8)*tc0*cse(3)*cse(3) &
     &    + sb(9)*tc0*tc0*tc0 + sb(10)*cse(3)*cse(3)*cse(3)

    a_out = a_
    b_out = b_
  end subroutine calc_coeffs_dry_snow_tm


end module ufo_directZDA_util_mod
