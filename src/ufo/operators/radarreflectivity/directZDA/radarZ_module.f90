! Adopted from RADAREMUL library of ARPS system
! Used for direct reflectity DA capability
!
! !DESCRIPTION: module containing radaremul library
!               which is used for TM operator of EnKF application
!
! !REVISION HISTORY:
!   2021-05-10  J. Park(CAPS) - modified with the GSI coding standards
!                             - radaremul_dualpolepara.f90 and
!                               radaremul_convert2radar.f90 are merged.
!   2021-05-17  J. Park(CAPS) - radaremul renamed to radarz
!   2021-10-20  J. Park(CAPS) - update radarZ with the latest GSI code
!   2022-04     J. Park(CAPS) - update radar Z libary to support NSSL

MODULE RADARZ_MODULE
!########################################################################
!########################################################################
!#########                                                      #########
!#########                  module radarz_module                #########
!#########                                                      #########
!#########                     Developed by                     #########
!#########     Center for Analysis and Prediction of Storms     #########
!#########                University of Oklahoma                #########
!#########                                                      #########
!########################################################################
!########################################################################

  use kinds, only: kind_real,kind_single,kind_int
  use radarz_iface

!-----------------------------------------------------------------------
!
! PURPOSE:
!
! Declare some constants used for calculaion of dual polarization
! parameters such as Zhh, Zdr, and Kdp. (It can be expanded...)
!
!-----------------------------------------------------------------------
!
! AUTHOR:  Youngsun Jung, 12/3/2004
!
!-----------------------------------------------------------------------
! Declare parameters.
!-----------------------------------------------------------------------

  IMPLICIT NONE
  SAVE
  PRIVATE

  PUBLIC :: calcMDR
  PUBLIC :: calcMu
  PUBLIC :: cal_rho
  PUBLIC :: init_fox
  PUBLIC :: init_fox_no_grpl
  PUBLIC :: init_fox_no_hail
  PUBLIC :: model_dsd
  PUBLIC :: coeff_hail
  PUBLIC :: coeff_grpl
  PUBLIC :: assign_Refl
  PUBLIC :: init_Refl
  PUBLIC :: init_para_dsd
  PUBLIC :: assign_para_dsd_TM
  PUBLIC :: rainIceRefl
  PUBLIC :: calculate_obs
  PUBLIC :: snow_alpha_a
  PUBLIC :: snow_alpha_b
  PUBLIC :: hail_alpha_a
  PUBLIC :: hail_alpha_b
  PUBLIC :: grpl_alpha_a
  PUBLIC :: grpl_alpha_b
  PUBLIC :: partialRefRain
  PUBLIC :: partialRhoRain
  PUBLIC :: partialRefIce
  PUBLIC :: partialRhoIce
  PUBLIC :: fractionWater
  PUBLIC :: fractionWater_md
  PUBLIC :: fractionWater_temperature_snow
  PUBLIC :: fractionWater_temperature_hail
  PUBLIC :: power_mom
  PUBLIC :: calc_N0x_mp
  PUBLIC :: calc_N0x_melt
  PUBLIC :: gamma
  PUBLIC :: cal_N0
  PUBLIC :: calc_lamda_mp
  PUBLIC :: cal_Nt
  PUBLIC :: cal_lamda
  PUBLIC :: set_dsd_para
  PUBLIC :: rdr_obs


!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
! SUBROUTINES AND FUNCTIONS
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  CONTAINS


  SUBROUTINE calcMDR()
!-----------------------------------------------------------------------
!
! PURPOSE:
!
! Calculates mass-diameter relation based on MP scheme. 
!
!-----------------------------------------------------------------------
!
! Author:  Bryan Putnam, 4/16/2013
!
!-----------------------------------------------------------------------

  IMPLICIT NONE

  c_x(1) = (pi/6._kind_real)*rhor
  c_x(2) = (pi/6._kind_real)*rhor
  c_x(3) = 440.0_kind_real

  SELECT CASE (mp_option)
    CASE(1:12,14,106,109:110,116)
      c_x(4) = (pi/6._kind_real)*rhos
    CASE(108)
      c_x(4) = .069_kind_real 
  END SELECT 

  c_x(5) = (pi/6._kind_real)*rhog
  c_x(6) = (pi/6._kind_real)*rhoh

  END SUBROUTINE calcMDR


  SUBROUTINE calcMu()

!-----------------------------------------------------------------------
!
! PURPOSE:
!
! Set the values for Mu depending upon microphysics scheme
!
!-----------------------------------------------------------------------
!
! AUTHOR:  Jonathan Labriola, 10/03/2016 
!
!-----------------------------------------------------------------------
! Force explicit declarations.
!-----------------------------------------------------------------------

  IMPLICIT NONE

  T_mur = 1.0_kind_real/3.0_kind_real
  T_mus = 1.0_kind_real/3.0_kind_real
  T_mug = 1.0_kind_real/3.0_kind_real
  T_muh = 1.0_kind_real/3.0_kind_real

  IF (mp_option == 14 .AND. murain == 3)  T_mur = 1.0_kind_real
  IF (mp_option == 14 .AND. musnow == 3)  T_mus = 1.0_kind_real

  END SUBROUTINE calcMu

  SUBROUTINE init_fox()

!-----------------------------------------------------------------------
!
! PURPOSE:
!
! Setup default maximum fraction of water in the melting ice.
!
!-----------------------------------------------------------------------
!
! AUTHOR:  Youngsun Jung, 1/20/2006
!
!-----------------------------------------------------------------------
! Force explicit declarations.
!-----------------------------------------------------------------------

  IMPLICIT NONE

!-----------------------------------------------------------------------
! Variables can vary depend on whether graupel/hail exists. 
!-----------------------------------------------------------------------
  fos = 0.3_kind_real             ! Maximum fraction of rain-snow mixture
  foh = 0.2_kind_real             ! Maximum fraction of rain-hail mixture
  fog = 0.25_kind_real            ! Maximum fraction of rain-hail mixture

  END SUBROUTINE init_fox

  SUBROUTINE init_fox_no_grpl()

!-----------------------------------------------------------------------
!
! PURPOSE:
!
! Setup default maximum fraction of water in the melting ice 
! when graupel is suppressed.
!
!-----------------------------------------------------------------------
!
! AUTHOR:  Youngsun Jung, 1/20/2006
!
!-----------------------------------------------------------------------
! Force explicit declarations.
!-----------------------------------------------------------------------

  IMPLICIT NONE

!-----------------------------------------------------------------------
! Variables can be changed by parameter retrieval
!-----------------------------------------------------------------------
  fos = 0.5_kind_real         ! Maximum fraction of rain-snow mixture
  foh = 0.3_kind_real         ! Maximum fraction of rain-hail mixture
  fog = 0.0_kind_real         ! Maximum fraction of rain-hail mixture
      
  END SUBROUTINE init_fox_no_grpl


  SUBROUTINE init_fox_no_hail() 

!-----------------------------------------------------------------------
!
! PURPOSE:  
!
!  Setup default maximum fraction of water in the melting ice 
!  when hail is suprressed. 
!
!-----------------------------------------------------------------------
!
! AUTHOR: Bryan Putnam, 12/14/10
!
!-----------------------------------------------------------------------
! Force explicit declarations. 
!-----------------------------------------------------------------------

  IMPLICIT NONE 

!-----------------------------------------------------------------------
! Variables can be changed by parameter retrieval 
!-----------------------------------------------------------------------

  fos = 0.5_kind_real      ! Maximum fraction of rain-snow mixture
  foh = 0.0_kind_real      ! Maximum fraction of rain-hail mixture
  fog = 0.3_kind_real      ! Maximum fraction of rain-hail mixture

  END SUBROUTINE init_fox_no_hail

  SUBROUTINE model_dsd(n0rain,n0snow,n0hail,n0grpl,rhosnow,rhohail,rhogrpl, &
             alpharain,alphasnow,alphagrpl,alphahail)

!-----------------------------------------------------------------------
!
! PURPOSE:
!
! Set dsd values to those used in the arps forecasts
!
!-----------------------------------------------------------------------
!
! AUTHOR:  Youngsun Jung, 1/20/2006
!
!-----------------------------------------------------------------------
! Force explicit declarations.
!-----------------------------------------------------------------------

  IMPLICIT NONE

!-----------------------------------------------------------------------
! Define variables:
!-----------------------------------------------------------------------


  REAL(kind_real), INTENT(IN    ) :: n0rain,n0snow,n0hail,n0grpl
  REAL(kind_real), INTENT(IN    ) :: rhosnow,rhohail,rhogrpl
  REAL(kind_real), INTENT(IN    ) :: alpharain,alphasnow,alphagrpl,alphahail

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
! Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  N0r=n0rain
  N0s=n0snow
  N0h=n0hail
  N0g=n0grpl

  rhos=rhosnow
  rhoh=rhohail
  rhog=rhogrpl

  alphar = alpharain
  alphas = alphasnow
  alphag = alphagrpl
  alphah = alphahail

  END SUBROUTINE model_dsd

  SUBROUTINE coeff_hail(fw)

!-----------------------------------------------------------------------
!
! PURPOSE:
!
! Scattering matrix coefficient for hail
!
!-----------------------------------------------------------------------
!
! AUTHOR:  Youngsun Jung, 1/27/2007
!
!-----------------------------------------------------------------------
! Force explicit declarations.
!-----------------------------------------------------------------------

  IMPLICIT NONE

!-----------------------------------------------------------------------
! Define variables:
!-----------------------------------------------------------------------

  REAL(kind_real), INTENT(IN   ) :: fw
  
  sf = 0.8_kind_real

  sigmah=pi/3_kind_real*(1_kind_real-sf*fw)
  Ah=.125_kind_real*(3_kind_real+4_kind_real*exp(-2_kind_real*sigmah**2)+exp(-8_kind_real*sigmah**2))
  Bh=.125_kind_real*(3_kind_real-4_kind_real*exp(-2_kind_real*sigmah**2)+exp(-8_kind_real*sigmah**2))
  Ch=.125_kind_real*(1_kind_real-exp(-8_kind_real*sigmah**2))
  Dh=.125_kind_real*(3_kind_real+exp(-8_kind_real*sigmah**2))
  Ckh=exp(-2_kind_real*sigmah**2)

  END SUBROUTINE coeff_hail

SUBROUTINE coeff_grpl(fw)

!-----------------------------------------------------------------------
!
! PURPOSE:
!
! Scattering matrix coefficient for graupel
!
!-----------------------------------------------------------------------
!
! AUTHOR:  Youngsun Jung, 1/27/2007
!
! MODIFIED: Dan Dawson, 02/16/2012
!           Made separate version for graupel.
!-----------------------------------------------------------------------
! Force explicit declarations.
!-----------------------------------------------------------------------

  IMPLICIT NONE

!-----------------------------------------------------------------------
! Define variables:
!-----------------------------------------------------------------------

  REAL(kind_real), INTENT(IN   ) :: fw

  sf = 0.8_kind_real

  sigmag=pi/3_kind_real*(1_kind_real-sf*fw)
  Ag=.125_kind_real*(3_kind_real+4_kind_real*exp(-2_kind_real*sigmag**2)+exp(-8_kind_real*sigmag**2))
  Bg=.125_kind_real*(3_kind_real-4_kind_real*exp(-2_kind_real*sigmag**2)+exp(-8_kind_real*sigmag**2))
  Cg=.125_kind_real*(1_kind_real-exp(-8_kind_real*sigmag**2))
  Dg=.125_kind_real*(3_kind_real+exp(-8_kind_real*sigmag**2))
  Ckg=exp(-2_kind_real*sigmag**2)

  END SUBROUTINE coeff_grpl

SUBROUTINE cal_rho(rho,q,v,rhox,cx)
!-----------------------------------------------------------------------
!
! PURPOSE:
!
! Scattering matrix coefficient for graupel
!
!-----------------------------------------------------------------------
!
! AUTHOR:  Jon Labriola 10/29/2016
!
!-----------------------------------------------------------------------
! Force explicit declarations.
!-----------------------------------------------------------------------

  IMPLICIT NONE

!-----------------------------------------------------------------------
! Define variables:
!-----------------------------------------------------------------------
  REAL(kind_real),INTENT(IN) :: rho,q,v
  REAL(kind_real),INTENT(INOUT) :: rhox,cx

  rhox = (rho*q)/v
  cx = (pi/6._kind_real)*rhox

  END SUBROUTINE cal_rho



  TYPE(T_obs_dual) FUNCTION assign_Refl(var1,var2,var3,var4)
       REAL(kind_real), INTENT(IN   ) :: var1,var2,var3,var4
       assign_Refl%T_sum_ref_h = var1
       assign_Refl%T_sum_ref_v = var2
       assign_Refl%T_log_zdr = var3
       assign_Refl%T_log_ref = var4
  END FUNCTION assign_Refl

  TYPE(T_obs_dual) FUNCTION init_Refl()
       init_Refl%T_sum_ref_h = 0._kind_real
       init_Refl%T_sum_ref_v = 0._kind_real
       init_Refl%T_log_zdr = missing
       init_Refl%T_log_ref = 0._kind_real
       init_Refl%T_sum_ref_hv = 0._kind_real
       init_Refl%T_kdp = 0._kind_real
       init_Refl%T_Ahh = 0._kind_real
       init_Refl%T_Avv = 0._kind_real
       init_Refl%T_ref_r_h = 0._kind_real
       init_Refl%T_ref_s_h = 0._kind_real
       init_Refl%T_ref_h_h = 0._kind_real
       init_Refl%T_ref_g_h = 0._kind_real
       init_Refl%T_ref_rs_h = 0._kind_real
       init_Refl%T_ref_rh_h = 0._kind_real
       init_Refl%T_ref_rg_h = 0._kind_real
       init_Refl%T_ref_r_v = 0._kind_real
       init_Refl%T_ref_s_v = 0._kind_real
       init_Refl%T_ref_h_v = 0._kind_real
       init_Refl%T_ref_g_v = 0._kind_real
       init_Refl%T_ref_rs_v = 0._kind_real
       init_Refl%T_ref_rh_v = 0._kind_real
       init_Refl%T_ref_rg_v = 0._kind_real
  END FUNCTION init_Refl

  TYPE(T_para_dsd) FUNCTION init_para_dsd()
    init_para_dsd%T_qr = 0.0_kind_real
    init_para_dsd%T_qs = 0.0_kind_real
    init_para_dsd%T_qh = 0.0_kind_real
    init_para_dsd%T_qg = 0.0_kind_real
    init_para_dsd%T_Ntr = 0.0_kind_real
    init_para_dsd%T_Nts = 0.0_kind_real
    init_para_dsd%T_Nth = 0.0_kind_real
    init_para_dsd%T_Ntg = 0.0_kind_real
    init_para_dsd%T_alfr = 0.0_kind_real
    init_para_dsd%T_alfs = 0.0_kind_real
    init_para_dsd%T_alfh = 0.0_kind_real
    init_para_dsd%T_alfg = 0.0_kind_real
    init_para_dsd%T_vg= 0.0_kind_real
    init_para_dsd%T_vh= 0.0_kind_real
  END FUNCTION init_para_dsd

  TYPE(T_para_dsd) FUNCTION assign_para_dsd_TM(var1,var2,var3,var4, &
                            var5,var6,var7,var8,var9,var10,var11,var12,var13,var14)
    REAL(kind_real), INTENT(IN   ) :: var1,var2,var3,var4,var5,var6,var7,var8
    REAL(kind_real), INTENT(IN   ) :: var9,var10,var11,var12,var13,var14

    assign_para_dsd_TM%T_qr = var1
    assign_para_dsd_TM%T_qs = var2
    assign_para_dsd_TM%T_qh = var3
    assign_para_dsd_TM%T_qg = var4
    assign_para_dsd_TM%T_Ntr = var5
    assign_para_dsd_TM%T_Nts = var6
    assign_para_dsd_TM%T_Nth = var7
    assign_para_dsd_TM%T_Ntg = var8
    assign_para_dsd_TM%T_alfr = var9
    assign_para_dsd_TM%T_alfs = var10
    assign_para_dsd_TM%T_alfh = var11
    assign_para_dsd_TM%T_alfg = var12
    ! Volume JDl
    assign_para_dsd_TM%T_vg = var13
    assign_para_dsd_TM%T_vh = var14
  END FUNCTION assign_para_dsd_TM


TYPE(T_obs_dual) FUNCTION rainIceRefl(var_dsd,rho,flg)
!-----------------------------------------------------------------------
!
! PURPOSE:
!
! This subroutine calculates the partial reflectivity factor
! of melting(wet) snow/hail at horizontal polarization
! and compute total reflectivity as a sum of those.
! The same formula used in shfactor is used with different
! alpha and beta coefficients that contain the effect of the fraction
! of water in the melting snow to take the melting layer into account.
!
!-----------------------------------------------------------------------
!
! AUTHOR:  Youngsun Jung, 1/29/2007
!
!-----------------------------------------------------------------------
! Force explicit declarations.
!-----------------------------------------------------------------------

  IMPLICIT NONE

!-----------------------------------------------------------------------
! Declare local variables.
!-----------------------------------------------------------------------
 
  TYPE(T_para_dsd), INTENT(IN   ) :: var_dsd
  INTEGER(kind_int),  INTENT(IN   ) :: flg
  REAL(kind_real),     INTENT(IN   ) :: rho

  REAL(kind_real) :: qr,qs,qh,qg,ntr,nts,nth,ntg 
  REAL(kind_real) :: vg,vh
  REAL(kind_real) :: rainIceRefl_hh,rainIceRefl_vv,rainIceRefl_hv,zdr
  REAL(kind_real) :: fracqrs,fracqrh,fracqrg
  REAL(kind_real) :: fracqs,fracqh,fracqg
  REAL(kind_real) :: fms,fmh,fmg,fws,fwh,fwg,rhoms,rhomh,rhomg
  REAL(kind_real) :: qrf,qsf,qhf,qgf
  REAL(kind_real) :: alphaa_ws,alphab_ws,alphaa_wh,alphab_wh,alphaa_wg,alphab_wg
  REAL(kind_real) :: alphak_ws,alphak_wh,alphak_wg
  REAL(kind_real) :: rainReflH,ZdrysnowH,ZwetsnowH
  REAL(kind_real) :: rainReflV,ZdrysnowV,ZwetsnowV
  REAL(kind_real) :: ZdryhailH,ZwethailH,ZdrygrplH,ZwetgrplH
  REAL(kind_real) :: ZdryhailV,ZwethailV,ZdrygrplV,ZwetgrplV
  REAL(kind_real) :: rainReflHV,ZdrysnowHV,ZwetsnowHV
  REAL(kind_real) :: ZdryhailHV,ZwethailHV,ZdrygrplHV,ZwetgrplHV
  REAL(kind_real) :: log_ref
  REAL(kind_real) :: rho_0rs,rho_0rh,rho_0rg,temp
  REAL(kind_real) :: temH,temV,temHV

  REAL(kind_real) :: tair_C
  REAL(kind_real) :: z_snow_thom
  REAL(kind_real) :: ntms, ntmh, ntmg
  REAL(kind_real) :: gamma,exp_term,gam_term
  logical :: firstcall = .true.
  SAVE firstcall

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
! Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  IF(firstcall) THEN
    SELECT CASE (qgh_opt)
      CASE (1)
        fos = 0.5_kind_real; foh = 0.0_kind_real; fog = 0.0_kind_real
      CASE (2)
        CALL init_fox_no_grpl()
      CASE (3)
        CALL init_fox_no_hail()
      CASE (4)
        CALL init_fox()
    END SELECT

    firstcall = .false. 
  END IF

  qrf = 0._kind_real; qsf = 0._kind_real; qhf = 0._kind_real; qgf = 0._kind_real
  fracqs = 0._kind_real; fracqh = 0._kind_real; fracqg = 0._kind_real
  fracqrs = 0._kind_real; fracqrh = 0._kind_real; fracqrg = 0._kind_real

  fms = 0._kind_real; fmh = 0._kind_real; fmg = 0._kind_real
  fws = 0._kind_real; fwh = 0._kind_real; fwg = 0._kind_real
  rhoms = 100._kind_real; rhomh = 913._kind_real; rhomg = 400._kind_real

  rainReflH = 0._kind_real
  rainReflV = 0._kind_real
  rainReflHV = 0._kind_real
  ZdrysnowH = 0._kind_real
  ZdrysnowV = 0._kind_real
  ZdrysnowHV = 0._kind_real
  ZwetsnowH = 0._kind_real
  ZwetsnowV = 0._kind_real
  ZwetsnowHV = 0._kind_real
  ZdryhailH = 0._kind_real
  ZdryhailV = 0._kind_real
  ZdryhailHV = 0._kind_real
  ZwethailH = 0._kind_real
  ZwethailV = 0._kind_real
  ZwethailHV = 0._kind_real
  ZdrygrplH = 0._kind_real
  ZdrygrplV = 0._kind_real
  ZdrygrplHV = 0._kind_real
  ZwetgrplH = 0._kind_real
  ZwetgrplV = 0._kind_real
  ZwetgrplHV = 0._kind_real

  temH = 0._kind_real
  temV = 0._kind_real
  temHV = 0._kind_real 

  rainIceRefl_hh = 0._kind_real
  rainIceRefl_vv = 0._kind_real
  rainIceRefl_hv = 0._kind_real
  zdr = missing
  log_ref = 0._kind_real

  rho_0rs = rho_0rsf
  rho_0rh = rho_0rhf
  rho_0rg = rho_0rgf

  qr = var_dsd%T_qr
  qs = var_dsd%T_qs
  qh = var_dsd%T_qh
  qg = var_dsd%T_qg
  ntr = var_dsd%T_Ntr
  nts = var_dsd%T_Nts
  nth = var_dsd%T_Nth
  ntg = var_dsd%T_Ntg
  vg = var_dsd%T_vg
  vh = var_dsd%T_vh
  
  if(qr < 0.0_kind_real) qr =0.0_kind_real
  if(qs < 0.0_kind_real) qs =0.0_kind_real
  if(qh < 0.0_kind_real) qh =0.0_kind_real
  if(qg < 0.0_kind_real) qg =0.0_kind_real

 ntms = 0._kind_real; ntmh = 0._kind_real; ntmg = 0._kind_real


!-----------------------------------------------------------------------
! Calculate hydrometeor density if necessary
!-----------------------------------------------------------------------
  IF (mp_option == 14) THEN
    IF (qg > 0._kind_real .AND. vg > 0._kind_real .AND. ntg > 0) CALL cal_rho(rho,qg,vg,rhog,c_x(5))
    IF (qh > 0._kind_real .AND. vh > 0._kind_real .AND. nth > 0) CALL cal_rho(rho,qh,vh,rhoh,c_x(6))
    IF(qh >  500._kind_real .OR. vh > 20._kind_real .OR. qg > 500._kind_real .OR. vg > 20._kind_real) PRINT *,"JDL it looks like and ERROR-WARNING"
  END IF

!-----------------------------------------------------------------------
! Calculate the fraction of water and ice.
!   qrf  pure rain water mixing ratio
!   qsf  dry snow mixing ratio
!   qhf  dry hail mixing ratio
!   qgf  dry graupel mixing ratio
!   fms  wet snow mixing ratio
!   fmh  wet hail mixing ratio
!   fmg  wet graupel mixing ratio
!   rhoms  density of wet snow (kg/m-3)
!   rhomh  density of wet hail (kg/m-3)
!   rhomg  density of wet graupel (kg/m-3)
!-----------------------------------------------------------------------

  IF (MFflg == 0) THEN

    CALL fractionWater(qr,qs,fos,rhos,fracqrs,fracqs,fms,fws,rhoms)
    IF(hail_ON == 1)  &
      CALL fractionWater(qr,qh,foh,rhoh,fracqrh,fracqh,fmh,fwh,rhomh)
    IF(graupel_ON == 1) &
      CALL fractionWater(qr,qg,fog,rhog,fracqrg,fracqg,fmg,fwg,rhomg)

    qrf = qr - fracqrs - fracqrh - fracqrg
    if(qrf < 0.0_kind_real) qrf = 0.0_kind_real
    qsf = qs - fracqs
    if(qsf < 0.0_kind_real) qsf = 0.0_kind_real
    qhf = qh - fracqh
    if(qhf < 0.0_kind_real) qhf = 0.0_kind_real
    qgf = qg - fracqg
    if(qgf < 0.0_kind_real) qgf = 0.0_kind_real

    ntms = nts
    ntmh = nth
    ntmg = ntg

  ELSE IF (MFflg == 2) THEN

    qrf = qr
    qsf = qs; fms = 0.0_kind_real;
    IF(hail_ON == 1)   qhf = qh; fmh = 0.0_kind_real;
    IF(graupel_ON == 1) qgf = qg; fmg = 0.0_kind_real;

    ntms = nts
    ntmh = nth
    ntmg = ntg

  ELSE IF (MFflg == 3) THEN    ! Temperature-based melting.

    tair_C = ta - degKtoC

    CALL fractionWater_temperature_snow(qs,rhos,fms,fws,rhoms,tair_C)
    IF(hail_ON == 1)  &
    CALL fractionWater_temperature_hail(qh,rhoh,fmh,fwh,rhomh,tair_C)
    IF(graupel_ON == 1)  &
    CALL fractionWater_temperature_hail(qg,rhog,fmg,fwg,rhomg,tair_C)

    qrf = qr
    qsf = qs-fms
    qhf = qh-fmh
    qgf = qg-fmg

    ntms = nts
    ntmh = nth
    ntmg = ntg

  ! may need MFflg4 - !fractionWater with median diameter preservation
  ELSE IF (MFflg == 4) THEN    ! Temperature-based melting.
    CALL fractionWater(qr,qs,fos,rhos,fracqrs,fracqs,fms,fws,rhoms)
    IF(hail_ON == 1)  &
      CALL fractionWater_md(qr,qh,foh,rhoh,fracqrh,fracqh,fmh,fwh,rhomh,rho,nth,ntmh)
!    according to the FV3 code, hail is not considered for wet fraction.
    IF(graupel_ON == 1) &
      CALL fractionWater_md(qr,qg,fog,rhog,fracqrg,fracqg,fmg,fwg,rhomg,rho,ntg,ntmg)

     qrf = qr - fracqrs - fracqrh - fracqrg
     if(qrf < 0.0_kind_real) qrf = 0.0_kind_real
     qsf = qs - fracqs
     if(qsf < 0.0_kind_real) qsf = 0.0_kind_real
     qhf = qh - fracqh
     if(qhf < 0.0_kind_real) qhf = 0.0_kind_real
     qgf = qg - fracqg
     if(qgf < 0.0_kind_real) qgf = 0.0_kind_real

    ntms = nts

  END IF

!  qrf = qr - fracqrs - fracqrh - fracqrg
  if(qrf < 0.0_kind_real) qrf = 0.0_kind_real
!  qsf = qs - fracqs
  if(qsf < 0.0_kind_real) qsf = 0.0_kind_real
!  qhf = qh - fracqh
  if(qhf < 0.0_kind_real) qhf = 0.0_kind_real
!  qgf = qg - fracqg
  if(qgf < 0.0_kind_real) qgf = 0.0_kind_real

!-----------------------------------------------------------------------
! Calculate the matrix coefficient for hail (Ah,Bh,Ch,Ckh)
!-----------------------------------------------------------------------
! fmh and fmg is not calculated.argument is missing
  IF(hail_ON == 1)   CALL coeff_hail(fwh)
  IF(graupel_ON == 1) CALL coeff_grpl(fwg)

!-----------------------------------------------------------------------
! Calculate alpha values
!-----------------------------------------------------------------------
  IF(fms > 0._kind_real) THEN
    alphaa_ws = snow_alpha_a(fws)
    alphab_ws = snow_alpha_b(fws)
    alphak_ws = alphaa_ws - alphab_ws
  ENDIF

  IF(hail_ON == 1 .and. fmh > 0._kind_real) THEN
    alphaa_wh = hail_alpha_a(fwh)
    alphab_wh = hail_alpha_b(fwh)
    alphak_wh = alphaa_wh - alphab_wh
  ENDIF

  IF(graupel_ON == 1 .and. fmg > 0._kind_real) THEN
    alphaa_wg = grpl_alpha_a(fwg)
    alphab_wg = grpl_alpha_b(fwg)
    alphak_wg = alphaa_wg - alphab_wg
  ENDIF

!-----------------------------------------------------------------------
! Calculate rho_0rs, rho_0rh, and rho_0rg
!-----------------------------------------------------------------------
  IF(flg > 2 .and. fms > 0._kind_real) THEN
    temp=rho*fms*1.e3_kind_real
    if(temp > 1._kind_real) then
      rho_0rs = rho_0rsi
    else if (1.e-2_kind_real > temp .and. temp <= 1._kind_real) then
      rho_0rs = rho_0rsi - .5_kind_real*log10(temp)*(rho_0rsf-rho_0rsi)
    endif
  ENDIF

  IF(hail_ON == 1 .and. flg > 2 .and. fmh > 0._kind_real) THEN
    temp=rho*fmh*1.e3_kind_real
    if(temp > 1._kind_real) then
      rho_0rh = rho_0rhi
    else if (1.e-2_kind_real > temp .and. temp <= 1._kind_real) then
      rho_0rh = rho_0rhi - .5_kind_real*log10(temp)*(rho_0rhf-rho_0rhi)
    endif
  ENDIF

  IF(graupel_ON == 1 .and. flg > 2 .and. fmg > 0._kind_real) THEN
    temp=rho*fmg*1.e3_kind_real
    if(temp > 1._kind_real) then
      rho_0rg = rho_0rgi
    else if (1.e-2_kind_real > temp .and. temp <= 1._kind_real) then
      rho_0rg = rho_0rgi - .5_kind_real*log10(temp)*(rho_0rgf-rho_0rgi)
    endif
  ENDIF

!-----------------------------------------------------------------------
! Calculate reflectivity (Zhh and Zvv (and Zhv, if necessary))
!-----------------------------------------------------------------------

  CALL calc_N0x_mp(rho,rhoms,rhomh,rhomg,ntr,nts,ntms,nth,ntmh,ntg,ntmg, &
                   qrf,qsf,fms,qhf,fmh,qgf,fmg)
  CALL calc_lamda_mp(rho,rhoms,rhomh,rhomg,ntr,nts,ntms,nth,ntmh,ntg,ntmg,  &
                    qrf,qsf,fms,qhf,fmh,qgf,fmg)

  SELECT CASE (mp_option)
    CASE(1:12,14,106,109:110,116)
      IF(lamdas > 0._kind_real .and. N0s > 0._kind_real) THEN
        CALL partialRefIce(N0s,dble(alphas),As,Bs,Cs,alphaa_ds,       &
                           alphab_ds,beta_sa, beta_sb,dble(lamdas),   & 
                           ZdrysnowH,ZdrysnowV,dble(T_mus))
        IF(flg > 2) THEN
          CALL partialRhoIce(N0s,dble(alphas),Cs,Ds,alphaa_ds,        &
                           alphab_ds,beta_sa,beta_sb,rho_0s,    &
                           dble(lamdas),ZdrysnowHV,dble(T_mus))
        ENDIF
      ENDIF
      IF(lamdams > 0._kind_real .and. N0ms > 0._kind_real) THEN
        CALL partialRefIce(N0ms,dble(alphas),As,Bs,Cs,alphaa_ws,       &
                           alphab_ws,beta_sa,beta_sb,dble(lamdams),    &
                           ZwetsnowH,ZwetsnowV,dble(T_mus))
        IF(flg > 2) THEN
          CALL partialRhoIce(N0ms,dble(alphas),Cs,Ds,alphaa_ws,        &
                           alphab_ws,beta_sa,beta_sb,rho_0rs,    &
                           dble(lamdams),ZwetsnowHV,dble(T_mus))
        ENDIF
      ENDIF
    CASE(108)
      IF(lamdas > 0._kind_real .and. N0s > 0._kind_real .and. qsf > 0._kind_real) THEN
        CALL partialRefIce(N0s,dble(alphas),As,Bs,Cs,alphaa_tom_ds,     &
                          alphab_tom_ds,beta_tom_dsa,beta_tom_dsb,  &
                          dble(lamdas),ZdrysnowH,ZdrysnowV,dble(T_mus))

        CALL partialRefIce(N0s2,dble(alphas2),As,Bs,Cs,alphaa_tom_ds,     &
                          alphab_tom_ds,beta_tom_dsa,beta_tom_dsb,  &
                          dble(lamdas2),temH,temV,dble(T_mus)) 

        ZdrysnowH = ZdrysnowH + temH
        ZdrysnowV = ZdrysnowV + temV

        IF(flg > 2) THEN
          CALL partialRhoIce(N0s,dble(alphas),Cs,Ds,alphaa_ds,        &
                          alphab_ds,beta_tom_dsa,beta_tom_dsb,    &
                          rho_0s,dble(lamdas),ZdrysnowHV,dble(T_mus))
          CALL partialRhoIce(N0s2,dble(alphas2),Cs,Ds,alphaa_ds,        &
                          alphab_ds,beta_tom_dsa,beta_tom_dsb,    &
                          rho_0s,dble(lamdas2),temHV,dble(T_mus))

          ZdrysnowHV = ZdrysnowHV + temHV
        ENDIF
      END IF 
      IF(lamdams > 0._kind_real .and. N0ms > 0._kind_real .and. fms > 0._kind_real) THEN
         CALL partialRefIce(N0ms,dble(alphas),As,Bs,Cs,alphaa_ws,      &
                           alphab_ws,beta_sa,beta_sb,dble(lamdams),    &
                           ZwetsnowH,ZwetsnowV,dble(T_mus)) 
         IF(flg > 2) THEN
         CALL partialRhoIce(N0ms,dble(alphas),Cs,Ds,alphaa_ws,           &
                           alphab_ws,beta_sa,beta_sb,rho_0rs,      &
                           dble(lamdams),ZwetsnowHV,dble(T_mus))
         END IF
      ENDIF 
  END SELECT 


  IF(hail_ON == 1) THEN
    IF(lamdah > 0._kind_real .and. N0h > 0._kind_real .and. qhf > 0._kind_real)THEN
      CALL partialRefIce(N0h,dble(alphah),Ahd,Bhd,Chd,alphaa_dh,    &
                         alphab_dh,beta_ha,beta_hb,dble(lamdah),    &
                         ZdryhailH,ZdryhailV,dble(T_muh))
      IF(flg > 2) THEN
        CALL partialRhoIce(N0h,dble(alphah),Chd,Dhd,alphaa_dh,      &
                         alphab_dh,beta_ha,beta_hb,rho_0h,    &
                         dble(lamdah),ZdryhailHV,dble(T_muh))
      ENDIF
    ENDIF
    IF(lamdamh > 0._kind_real .and. N0mh > 0._kind_real .and. fmh > 0._kind_real) THEN
      CALL partialRefIce(N0mh,dble(alphah),Ah,Bh,Ch,alphaa_wh,       &
                         alphab_wh,beta_ha,beta_hb,dble(lamdamh),    &
                         ZwethailH,ZwethailV,dble(T_muh))
      IF(flg > 2) THEN
        CALL partialRhoIce(N0mh,dble(alphah),Ch,Dh,alphaa_wh,        &
                         alphab_wh,beta_ha,beta_hb,rho_0rh,    &
                         dble(lamdamh),ZwethailHV,dble(T_muh))
      ENDIF
    ENDIF
  ENDIF 

  IF(graupel_ON == 1) THEN
    IF(lamdag > 0._kind_real .and. N0g > 0._kind_real .and.  qgf > 0._kind_real)THEN
      CALL partialRefIce(N0g,dble(alphag),Agd,Bgd,Cgd,alphaa_dg,    &
                         alphab_dg,beta_ga, beta_gb,dble(lamdag),   &
                         ZdrygrplH,ZdrygrplV,dble(T_mug))
      IF(flg > 2) THEN
        CALL partialRhoIce(N0g,dble(alphag),Cgd,Dgd,alphaa_dg,      &
                         alphab_dg,beta_ga,beta_gb,rho_0g,    &
                         dble(lamdag),ZdrygrplHV,dble(T_mug))
      ENDIF
    ENDIF
     IF(lamdamg > 0._kind_real .and. N0mg > 0._kind_real .and. fmg > 0._kind_real) THEN 
      CALL partialRefIce(N0mg,dble(alphag),Ag,Bg,Cg,alphaa_wg,       &
                         alphab_wg,beta_ga,beta_gb,dble(lamdamg),    &
                         ZwetgrplH,ZwetgrplV,dble(T_mug))
      IF(flg > 2) THEN
        CALL partialRhoIce(N0mg,dble(alphag),Cg,Dg,alphaa_wg,        &
                         alphab_wg,beta_ga,beta_gb,rho_0rg,    &
                         dble(lamdamg),ZwetgrplHV,dble(T_mug))
      ENDIF
    ENDIF
  ENDIF

  IF(lamdar > 0._kind_real .and. N0r > 0._kind_real) THEN
    CALL partialRefRain(N0r,dble(alphar),alphaa,alphab,beta_ra,beta_rb,  &
                       dble(lamdar),rainReflH,rainReflV,dble(T_mur))
    rainReflV = MIN(rainReflV,rainReflH)
    IF(flg > 2) THEN

    CALL partialRhoRain(N0r,dble(alphar),alphaa,alphab,beta_ra,beta_rb,  &
                        dble(lamdar),rainReflHV,dble(T_mur))
    ENDIF
  ENDIF

  rainIceRefl_hh=rainReflH+ZdrysnowH+ZwetsnowH+ZdryhailH+ZwethailH &
                 +ZdrygrplH+ZwetgrplH

  log_ref = 10._kind_real*LOG10(MAX(1.0_kind_real,rainIceRefl_hh))

  IF(flg == 1) THEN
    rainIceRefl = assign_Refl(rainIceRefl_hh,rainIceRefl_vv,zdr,log_ref)

  ELSE IF(flg > 1) THEN
!-----------------------------------------------------------------------
! Calculate differential reflectivity (Zdr)
!-----------------------------------------------------------------------
    rainIceRefl_vv=rainReflV+ZdrysnowV+ZwetsnowV+ZdryhailV+ZwethailV &
                  +ZdrygrplV+ZwetgrplV

    if(rainIceRefl_vv > 0._kind_real) then
      zdr = 10._kind_real*LOG10(MAX(1.0_kind_real,rainIceRefl_hh/rainIceRefl_vv))
    endif

    rainIceRefl = assign_Refl(rainIceRefl_hh,rainIceRefl_vv,zdr,log_ref)

    IF(flg > 2) THEN

      rainIceRefl_hv=rainReflHV+ZdrysnowHV+ZwetsnowHV                  &
                   +ZdryhailHV+ZwethailHV+ZdrygrplHV+ZwetgrplHV

!-----------------------------------------------------------------------
! Safety block to ensure r_hv <= 1.
!-----------------------------------------------------------------------
      IF(rainIceRefl_hv > SQRT(rainIceRefl_hh*rainIceRefl_vv)) &
         rainIceRefl_hv = SQRT(rainIceRefl_hh*rainIceRefl_vv)
!-----------------------------------------------------------------------

      rainIceRefl%T_sum_ref_hv = rainIceRefl_hv

    ENDIF
  ENDIF

END FUNCTION rainIceRefl

!########################################################################
!########################################################################
!#########                                                      #########
!#########              FUNCTION calculate_obs                  #########
!#########                                                      #########
!#########                     Developed by                     #########
!#########     Center for Analysis and Prediction of Storms     #########
!#########                University of Oklahoma                #########
!#########                                                      #########
!########################################################################
!########################################################################

 TYPE(T_obs_dual) FUNCTION calculate_obs(rho,var_dsd,flg)

!-----------------------------------------------------------------------
!
! AUTHOR:  Youngsun Jung, 2/27/2007
!
! flg == (1: Zh, 2: Zdr, 3: rho_hv)
!-----------------------------------------------------------------------
! Force explicit declarations.
!-----------------------------------------------------------------------

  IMPLICIT NONE

!-----------------------------------------------------------------------
! Declare arguments.
!-----------------------------------------------------------------------

  REAL(kind_real),  INTENT(IN   ) :: rho ! Air density (kg m**-3)

  TYPE(T_para_dsd), INTENT(IN   ) :: var_dsd

  INTEGER(kind_int),  INTENT(IN   ) :: flg   ! flag for ref(1) and zdr(2)

  REAL(kind_real) :: qr
  REAL(kind_real) :: qs
  REAL(kind_real) :: qh
  REAL(kind_real) :: qg

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
! Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  calculate_obs = init_Refl()

!-----------------------------------------------------------------------
! Check for bad air density value.
!-----------------------------------------------------------------------
 
  qr = var_dsd%T_qr
  qs = var_dsd%T_qs
  qh = var_dsd%T_qh
  qg = var_dsd%T_qg


  IF (rho > 0.0_kind_real .and. &
      (qr > 0._kind_real .or. qs > 0._kind_real .or. qh > 0._kind_real .or. qg > 0._kind_real)) THEN 
     calculate_obs = rainIceRefl(var_dsd,rho,flg)
  END IF

END FUNCTION  calculate_obs

!--- from convert2radar
REAL(kind_real) FUNCTION snow_alpha_a(fw)

!-----------------------------------------------------------------------
!
! PURPOSE:
!
! This user defined function calculates alpha for Zhh
! for dry/melting snow.
!
!-----------------------------------------------------------------------
!
! AUTHOR:  Youngsun Jung, 5/30/2006
!
!-----------------------------------------------------------------------
! Force explicit declarations.
!-----------------------------------------------------------------------

  IMPLICIT NONE

!-----------------------------------------------------------------------
! Declare local variables.
!-----------------------------------------------------------------------

  REAL(kind_real),INTENT(IN   ) :: fw

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
! Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  snow_alpha_a = (0.194_kind_real + 7.094_kind_real*fw &
                  + 2.135_kind_real*fw**2._kind_real   &
                  - 5.225_kind_real*fw**3._kind_real)*10._kind_real**(-4)

END FUNCTION snow_alpha_a

REAL(kind_real) FUNCTION snow_alpha_b(fw)

!-----------------------------------------------------------------------
!
! PURPOSE:
!
! This user defined function calculates alpha for Zvv
! for dry/melting snow.
!
!-----------------------------------------------------------------------
!
! AUTHOR:  Youngsun Jung, 5/30/2006
!
!-----------------------------------------------------------------------
! Force explicit declarations.
!-----------------------------------------------------------------------

  IMPLICIT NONE

!-----------------------------------------------------------------------
! Declare local variables.
!-----------------------------------------------------------------------

  REAL(kind_real),INTENT(IN   ) :: fw

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
! Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  snow_alpha_b = (0.191_kind_real + 6.916_kind_real*fw &
                  - 2.841_kind_real*fw**2._kind_real &
                  - 1.160_kind_real*fw**3._kind_real)*10._kind_real**(-4)

END FUNCTION snow_alpha_b

!########################################################################
!########################################################################
!#########                                                      #########
!#########                  FUNCTIONS for hail                  #########
!#########                                                      #########
!#########                     Developed by                     #########
!#########     Center for Analysis and Prediction of Storms     #########
!#########                University of Oklahoma                #########
!#########                                                      #########
!########################################################################
!########################################################################

REAL(kind_real) FUNCTION hail_alpha_a(fw)

!-----------------------------------------------------------------------
!
! PURPOSE:
!
! This user defined function calculates alpha for Zhh
! for dry/melting hail.
!
!-----------------------------------------------------------------------
!
! AUTHOR:  Youngsun Jung, 1/22/2007
!
!-----------------------------------------------------------------------
! Force explicit declarations.
!-----------------------------------------------------------------------

  IMPLICIT NONE

!-----------------------------------------------------------------------
! Declare local variables.
!-----------------------------------------------------------------------

  REAL(kind_real),INTENT(IN   ) :: fw

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
! Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  hail_alpha_a = (0.191_kind_real + 2.39_kind_real*fw &
                  - 12.57_kind_real*fw**2._kind_real  &
                  + 38.71_kind_real*fw**3._kind_real  &
                  - 65.53_kind_real*fw**4._kind_real  &
                  + 56.16_kind_real*fw**5._kind_real  &
                  - 18.98_kind_real*fw**6._kind_real)*10._kind_real**(-3)

END FUNCTION hail_alpha_a

REAL(kind_real) FUNCTION hail_alpha_b(fw)

!-----------------------------------------------------------------------
!
! PURPOSE:
!
! This user defined function calculates alpha for Zvv
! for dry/melting hail.
!
!-----------------------------------------------------------------------
!
! AUTHOR:  Youngsun Jung, 1/22/2007
!
!-----------------------------------------------------------------------
! Force explicit declarations.
!-----------------------------------------------------------------------

  IMPLICIT NONE

!-----------------------------------------------------------------------
! Declare local variables.
!-----------------------------------------------------------------------

  REAL(kind_real),INTENT(IN   ) :: fw

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
! Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  hail_alpha_b = (0.165_kind_real + 1.72_kind_real*fw &
                  - 9.92_kind_real*fw**2._kind_real   &
                  + 32.15_kind_real*fw**3._kind_real  &
                  - 56.0_kind_real*fw**4._kind_real   &
                  + 48.83_kind_real*fw**5._kind_real  &
                  - 16.69_kind_real*fw**6._kind_real)*10._kind_real**(-3)

END FUNCTION hail_alpha_b
!
!########################################################################
!########################################################################
!#########                                                      #########
!#########                  FUNCTIONS for graupel               #########
!#########                                                      #########
!#########                     Developed by                     #########
!#########     Center for Analysis and Prediction of Storms     #########
!#########                University of Oklahoma                #########
!#########                                                      #########
!########################################################################
!########################################################################
!

REAL(kind_real) FUNCTION grpl_alpha_a(fw)

!-----------------------------------------------------------------------
!
! PURPOSE:
!
! This user defined function calculates alpha for Zhh
! for dry/melting graupel.
!
!-----------------------------------------------------------------------
!
! AUTHOR:  Youngsun Jung, 3/10/2010
!
!-----------------------------------------------------------------------
! Force explicit declarations.
!-----------------------------------------------------------------------

  IMPLICIT NONE

!-----------------------------------------------------------------------
! Declare local variables.
!-----------------------------------------------------------------------
  REAL(kind_real),INTENT(IN   ) :: fw

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
! Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  grpl_alpha_a = (0.081_kind_real + 2.04_kind_real*fw &
                  - 7.39_kind_real*fw**2._kind_real   &
                  + 18.14_kind_real*fw**3._kind_real  &
                  - 26.02_kind_real*fw**4._kind_real  &
                  + 19.37_kind_real*fw**5._kind_real  &
                  - 5.75_kind_real*fw**6._kind_real)*10._kind_real**(-3)

END FUNCTION grpl_alpha_a

REAL(kind_real) FUNCTION grpl_alpha_b(fw)

!-----------------------------------------------------------------------
!
! PURPOSE:
!
! This user defined function calculates alpha for Zvv
! for dry/melting graupel.
!
!-----------------------------------------------------------------------
!
! AUTHOR:  Youngsun Jung, 3/10/2010
!
!-----------------------------------------------------------------------
! Force explicit declarations.
!-----------------------------------------------------------------------

  IMPLICIT NONE

!-----------------------------------------------------------------------
! Declare local variables.
!-----------------------------------------------------------------------

  REAL(kind_real),INTENT(IN   ) :: fw

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
! Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  grpl_alpha_b = (0.076_kind_real + 1.74_kind_real*fw &
                  - 7.52_kind_real*fw**2._kind_real   &
                  + 20.22_kind_real*fw**3._kind_real  &
                  - 30.42_kind_real*fw**4._kind_real  &
                  + 23.31_kind_real*fw**5._kind_real  &
                  - 7.06_kind_real*fw**6._kind_real)*10._kind_real**(-3)

END FUNCTION grpl_alpha_b


! Adopted from RADREMUL library of ARPS system
! Used for direct reflectity DA capability
!########################################################################
!########################################################################
!#########                                                      #########
!#########                     convert2radar                    #########
!#########                                                      #########
!#########                     Developed by                     #########
!#########     Center for Analysis and Prediction of Storms     #########
!#########                University of Oklahoma                #########
!#########                                                      #########
!########################################################################
!########################################################################

SUBROUTINE partialRefRain(N0,alpha,alp_a,alp_b,beta_a,beta_b,lamda,    &
                           refRainHH,refRainVV,mu_x)

!-----------------------------------------------------------------------
!
! PURPOSE:
!
! This function calculates the partial reflectivity for rain species
!
!-----------------------------------------------------------------------
!
! AUTHOR:  Bryan Putnam, 4/16/2013
!
!-----------------------------------------------------------------------

  IMPLICIT NONE

!-----------------------------------------------------------------------
! Declare variables.
!-----------------------------------------------------------------------

  REAL(kind_real),INTENT(IN   ) :: N0,alpha,lamda,mu_x
  REAL(kind_real),INTENT(IN   ) :: alp_a,alp_b
  REAL(kind_real),INTENT(IN   ) :: beta_a,beta_b
  REAL(kind_real),INTENT(  OUT) :: refRainHH,refRainVV

  !local variables
  REAL(kind_real) :: gamma
  REAL(kind_real) :: N0_units,lamda_units_h,lamda_units_v
  REAL(kind_real) :: gamma_h,expon_h,gamma_v,expon_v

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
! Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  !N0_units = (1.e-3_kind_real)**(4.0_kind_real+alpha)

  lamda_units_h = (1.e-3_kind_real)**(-(alpha+2.0_kind_real*dble(beta_a)+1.0_kind_real))
  lamda_units_v = (1.e-3_kind_real)**(-(alpha+2.0_kind_real*dble(beta_b)+1.0_kind_real))
  N0_units = (1.e-3_kind_real)**(4.0_kind_real+alpha)

  gamma_h = gamma((alpha+2.0_kind_real*dble(beta_a)+1.0_kind_real)/(3.0_kind_real*mu_x))
  expon_h = -(alpha+2.0_kind_real*dble(beta_a)+1.0_kind_real)/(3.0_kind_real*mu_x)
  gamma_v = gamma((alpha+2.0_kind_real*dble(beta_b)+1.0_kind_real)/(3.0_kind_real*mu_x))
  expon_v = -(alpha+2.0_kind_real*dble(beta_b)+1.0_kind_real)/(3.0_kind_real*mu_x)

  refRainHH =sngl(dble(mm3todBZ*radar_const*alp_a**2.0_kind_real)*((N0/(3.0_kind_real*mu_x))*N0_units)*gamma_h*&
               (lamda**expon_h)*lamda_units_h)

  refRainVV =sngl(dble(mm3todBZ*radar_const*alp_b**2.0_kind_real)*((N0/(3.0_kind_real*mu_x))*N0_units)*gamma_v*&
               (lamda**expon_v)*lamda_units_v)


END SUBROUTINE partialRefRain


SUBROUTINE partialRhoRain(N0,alpha,alp_a,alp_b,beta_a,beta_b,         &
                          lamda,refRainHV,mu_x)

!-----------------------------------------------------------------------
!
! PURPOSE:
!
! This function calculates the cross components, Z_hv, for rain species
! for rho_hv calculation.
!
!-----------------------------------------------------------------------
!
! AUTHOR:  Bryan Putnam, 4/16/2013
!
!-----------------------------------------------------------------------

  IMPLICIT NONE

!-----------------------------------------------------------------------
! Declare variables.
!-----------------------------------------------------------------------

  REAL(kind_real),INTENT(IN   ) :: alp_a,alp_b,beta_a,beta_b
  REAL(kind_real),INTENT(IN   ) :: N0,alpha,lamda,mu_x

  REAL(kind_real),INTENT(  OUT) :: refRainHV

  !local variables
  REAL(kind_real) :: gamma
  REAL(kind_real) :: expon_hv,gamma_hv
  REAL(kind_real) :: N0_units,lamda_units_hv

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
! Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  lamda_units_hv = (1.e-3_kind_real)**(-(alpha+2.0_kind_real*dble(beta_a)+1.0_kind_real))
  N0_units = (1.e-3_kind_real)**(4.0_kind_real+alpha)

  gamma_hv = gamma((dble(beta_a)+dble(beta_b)+alpha+1.0_kind_real)/(3.0_kind_real*mu_x))
  expon_hv = -(alpha+dble(beta_a)+dble(beta_b)+1.0_kind_real)/(3.0_kind_real*mu_x)

  refRainHV = sngl(dble(mm3todBZ*radar_const*alp_a*alp_b)*((N0/(3.0_kind_real*mu_x))*N0_units)*&
                gamma_hv*(lamda**expon_hv)*lamda_units_hv)

END SUBROUTINE partialRhoRain


SUBROUTINE partialRefIce(N0,alpha,Ai,Bi,Ci,alp_a,alp_b,beta_a,beta_b,   &
                         lamda,refIceHH,refIceVV,mu_x)

!-----------------------------------------------------------------------
!
! PURPOSE:
!
! This function calculates the partial reflectivity for each species
!
!-----------------------------------------------------------------------
!
! AUTHOR:  Youngsun Jung, 1/22/2007
!
!-----------------------------------------------------------------------
! Force explicit declarations.
!-----------------------------------------------------------------------

  IMPLICIT NONE

!-----------------------------------------------------------------------
! Declare variables.
!-----------------------------------------------------------------------
  REAL(kind_real),INTENT(IN   ) :: Ai,Bi,Ci,alp_a,alp_b,beta_a,beta_b
  REAL(kind_real),INTENT(IN   ) :: N0, alpha, lamda,mu_x

  REAL(kind_real),INTENT(  OUT) :: refIceHH,refIceVV

  !local variables
  REAL(kind_real) :: gamma_h, gamma_v, expon_h, expon_v
  REAL(kind_real) :: N0_units, lamda_units_h, lamda_units_v
  REAL(kind_real) :: gamma

  gamma_h = gamma((alpha + 2.0_kind_real*dble(beta_a)+1.0_kind_real)/(3.0_kind_real*mu_x))
  expon_h = -(alpha+2*dble(beta_a)+1.0_kind_real)/(3.0_kind_real*mu_x)
  gamma_v = gamma((alpha + 2.0_kind_real*dble(beta_b)+1.0_kind_real)/(3.0_kind_real*mu_x))
  expon_v = -(alpha+2*dble(beta_b)+1.0_kind_real)/(3.0_kind_real*mu_x)

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
! Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  N0_units = (1.e-3_kind_real)**(4.0_kind_real+alpha)
  lamda_units_h = (1.e-3)**(expon_h*(3.d0*mu_x)) ! suggested by JDL
  lamda_units_v = (1.e-3)**(expon_v*(3.d0*mu_x)) ! suggested by JDL

  refIceHH = sngl(dble(mm3toDBZ*radar_Const)*gamma_h*(N0/(3.0_kind_real*mu_x))*N0_units*&
              dble(Ai*alp_a**2+Bi*alp_b**2+2*Ci*alp_a*alp_b)*                 &
             (lamda)**expon_h*lamda_units_h)

  refIceVV = sngl(dble(mm3toDBZ*radar_Const)*gamma_v*(N0/(3.0_kind_real*mu_x))*N0_units*&
             dble(Bi*alp_a**2+Ai*alp_b**2+2*Ci*alp_a*alp_b)*                  &
             (lamda)**expon_v*lamda_units_v)

END SUBROUTINE partialRefIce

SUBROUTINE partialRhoIce(N0,alpha,Ci,Di,alp_a,alp_b,beta_a,beta_b,rho_0,lamda,refIceHV,mu_x)

!-----------------------------------------------------------------------
!
! PURPOSE:
!
! This function calculates the cross components, Z_hv, for each species
! for rho_hv calculation.
!
!-----------------------------------------------------------------------
!
! AUTHOR:  Youngsun Jung, 12/16/2007
!
!-----------------------------------------------------------------------
! Force explicit declarations.
!-----------------------------------------------------------------------
  IMPLICIT NONE

!-----------------------------------------------------------------------
! Declare variables.
!-----------------------------------------------------------------------
  REAL(kind_real),INTENT(IN   ) :: rho_0,Ci,Di,alp_a,alp_b,beta_a,beta_b
  REAL(kind_real),INTENT(IN   ) :: N0,mu_x,lamda,alpha

  REAL(kind_real),INTENT(  OUT) :: refIceHV

  !local variables
  REAL(kind_real) :: gamma_hv, expon
  REAL(kind_real) :: N0_units,lamda_units
  REAL(kind_real) :: gamma

   gamma_hv = gamma((dble(beta_a)+dble(beta_b)+alpha + 1.0_kind_real)/(3.0_kind_real*mu_x))
   expon = -(alpha + dble(beta_a) + dble(beta_b) + 1.0_kind_real)/(3.0_kind_real*mu_x)

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
! Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

   lamda_units = lamda**(1.0_kind_real/(3.0_kind_real*mu_x)) *(1.e-3_kind_real)
   N0_units = N0**(1.0_kind_real/(1.0_kind_real+alpha)) * (1.e-3_kind_real)**4.0_kind_real

   refIceHV = sngl(dble(mm3todBZ*radar_Const)*gamma_hv*(N0_units/(3.0_kind_real*mu_x))*&
               dble(Ci*alp_a**2+Ci*alp_b**2+2*Di*alp_a*alp_b*rho_0)*    &
               (lamda_units)**expon)

END SUBROUTINE partialRhoIce

SUBROUTINE fractionWater(qr,qi,fo,density_ice,fracqr,fracqi,fm,fw,rhom)

!-----------------------------------------------------------------------
!
! PURPOSE:
!
! This subroutine calculates the fractions of water, dry ice (snow or
! hail), the mixture. It also calculate the density of mixture.
!
!-----------------------------------------------------------------------
!
! AUTHOR:  Youngsun Jung, 5/30/2006
!
!-----------------------------------------------------------------------
! Force explicit declarations.
!-----------------------------------------------------------------------

  IMPLICIT NONE

!-----------------------------------------------------------------------
! Declare variables.
!-----------------------------------------------------------------------

  REAL(kind_real),INTENT(IN   ) :: qr, qi, fo, density_ice
  REAL(kind_real),INTENT(  OUT) :: fracqr, fracqi, fm, fw, rhom

  REAL(kind_real) :: fr

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
! Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  fr = 0._kind_real
  fw = 0._kind_real
  fracqr = 0._kind_real
  fracqi = 0._kind_real
  fm = 0._kind_real
  rhom = 0._kind_real

!-----------------------------------------------------------------------
! Calculate the fraction of mleting ice (fr) based on the ratio between
! qr and qi. fo is the maximum allowable fraction of melting snow.
!-----------------------------------------------------------------------
  IF (qr > 0._kind_real .AND. qi > 0._kind_real) THEN
    fr = fo*(MIN(qi/qr,qr/qi))**.3_kind_real
  ENDIF

!-----------------------------------------------------------------------
! Calculate the faction of water and ice.
! fracqr : the mass of water in the melting ice
! fracqi : the mass of ice in the melting ice
! fm     : total mass of melting ice
! fw     : the fraction of water within melting ice
! rhom   : density of mixture
!-----------------------------------------------------------------------
  fracqr = fr * qr
  fracqi = fr * qi
  fm = fracqr + fracqi

  IF (fm == 0._kind_real .AND. qr > 0._kind_real) THEN
    fw = 1._kind_real
  ELSE IF (fm > 0._kind_real) THEN
    fw = fracqr/fm
  ENDIF

  rhom = 1000._kind_real*fw**2._kind_real + (1._kind_real-fw**2._kind_real)*density_ice

END SUBROUTINE fractionWater

SUBROUTINE fractionWater_temperature_snow (qi,density_ice,fm,fw,rhom,tair_C)

!-----------------------------------------------------------------------
!
! PURPOSE:
!
! This subroutine calculates the fractions of water, dry ice (snow or
! hail), the mixture based on the air temperature.
! It also calculate the density of mixture.
!
!-----------------------------------------------------------------------
!
! AUTHOR:  Youngsun Jung, 7/25/2014
!
!-----------------------------------------------------------------------
! Force explicit declarations.
!-----------------------------------------------------------------------

  IMPLICIT NONE

!-----------------------------------------------------------------------
! Declare variables.
!-----------------------------------------------------------------------

  REAL(kind_real),INTENT(IN   ) :: qi, density_ice, tair_C
  REAL(kind_real),INTENT(  OUT) :: fm, fw, rhom
  REAL(kind_real) :: layer_tmax, layer_tmin

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
! Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  fw = 0._kind_real
  fm = 0._kind_real
  rhom = 0._kind_real

!-----------------------------------------------------------------------
! Calculate the faction of water.
! fm     : total mass of melting ice
! fw     : the fraction of water within melting ice
! rhom   : density of mixture
!-----------------------------------------------------------------------

  fm = qi

! Compute the degree of wet in percentage based on air temperature
  layer_tmax = 2.5_kind_real
  layer_tmin = -2.5_kind_real
  if(tair_C >= layer_tmin .and. tair_C < layer_tmax) then
    fw = (tair_C - layer_tmin)/(layer_tmax-layer_tmin)
  else if(tair_C >= layer_tmax) then
    fw = 1._kind_real
  else
    fm = 0._kind_real
    fw = 0._kind_real
  endif

  rhom = 1000._kind_real*fw**2._kind_real + (1._kind_real-fw**2._kind_real)*density_ice

END SUBROUTINE fractionWater_temperature_snow

SUBROUTINE fractionWater_temperature_hail(qi,density_ice,fm,fw,rhom,tair_C)
  
!-----------------------------------------------------------------------
! 
! PURPOSE:
!
! This subroutine calculates the fractions of water, dry ice (snow or
! hail), the mixture based on the air temperature.
! It also calculate the density of mixture.
! 
!-----------------------------------------------------------------------
!
! AUTHOR:  Youngsun Jung, 7/25/2014
! 
!-----------------------------------------------------------------------
! Force explicit declarations.
!-----------------------------------------------------------------------

  IMPLICIT NONE

!-----------------------------------------------------------------------
! Declare variables.
!-----------------------------------------------------------------------

  REAL(kind_real),INTENT(IN   ) :: qi, density_ice, tair_C
  REAL(kind_real),INTENT(  OUT) :: fm, fw, rhom
  REAL(kind_real) :: layer_tmax, layer_tmin
  REAL(kind_real) :: maxfrac

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
! Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  fw = 0._kind_real
  fm = 0._kind_real
  rhom = 0._kind_real
  maxfrac = 0.6_kind_real

!-----------------------------------------------------------------------
! Calculate the faction of water.
! fm     : total mass of melting ice
! fw     : the fraction of water within melting ice
! rhom   : density of mixture
!-----------------------------------------------------------------------

  fm = qi

! Compute the degree of wet in percentage based on air temperature
  layer_tmax = 5.0_kind_real
  layer_tmin = 0.0_kind_real
  if(tair_C >= layer_tmin .and. tair_C < layer_tmax) then
    fw = (tair_C - layer_tmin)/(layer_tmax-layer_tmin) * maxfrac
  else if(tair_C >= layer_tmax) then
    fw = maxfrac
  else
    fm = 0._kind_real
    fw = 0._kind_real
  endif

  rhom = 1000._kind_real*fw**2._kind_real + (1._kind_real-fw**2._kind_real)*density_ice

END SUBROUTINE fractionWater_temperature_hail

SUBROUTINE fractionWater_md(qr,qi,fo,density_ice,fracqr,fracqi,fm,fw,rhom, &
                         rhoa,nti,ntm)

!-----------------------------------------------------------------------
!
! PURPOSE:
!
! This subroutine calculates the fractions of water, dry ice (snow or
! hail), in the mixture following fractionWater but accounts for preserving 
! median diamater by altering number concentration. 
!
!-----------------------------------------------------------------------
!
! AUTHOR:  Jon Labriola, 12/6/2017
!
!-----------------------------------------------------------------------
! Force explicit declarations.
!-----------------------------------------------------------------------

  IMPLICIT NONE

!-----------------------------------------------------------------------
! Declare variables.
!-----------------------------------------------------------------------

  REAL(kind_real),INTENT(IN ) :: qr, qi, fo, density_ice, rhoa
  REAL(kind_real),INTENT(OUT) :: fracqr, fracqi, fm, fw, rhom
  REAL(kind_real) :: fr

  REAL(kind_real) :: dm
  REAL(kind_real), INTENT(INOUT) :: ntm,nti
  REAL(kind_real), PARAMETER :: pisix = 0.523599_kind_real
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
! Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  fr = 0._kind_real
  fw = 0._kind_real
  fracqr = 0._kind_real
  fracqi = 0._kind_real
  fm = 0._kind_real
  rhom = 0._kind_real

!-----------------------------------------------------------------------
! Calculate the fraction of mleting ice (fr) based on the ratio between
! qr and qi. fo is the maximum allowable fraction of melting snow.
!-----------------------------------------------------------------------
  IF (qr > 0._kind_real .AND. qi > 0._kind_real) THEN
    fr = fo*(MIN(qi/qr,qr/qi))**.3_kind_real
  ENDIF

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
! Calculate the faction of water and ice.
! fracqr : the mass of water in the melting ice
! fracqi : the mass of ice in the melting ice
! fm     : total mass of melting ice
! fw     : the fraction of water within melting ice
! rhom   : density of mixture
!-----------------------------------------------------------------------
  fracqr = fr * qr
  fracqi = fr * qi
  fm = fracqr + fracqi

  IF (fm .EQ. 0._kind_real .AND. qr > 0._kind_real) THEN
    fw = 1._kind_real
  ELSE IF (fm > 0._kind_real) THEN
    fw = fracqr/fm
  ENDIF

  rhom = 1000._kind_real*fw**2._kind_real + (1._kind_real-fw**2._kind_real)*density_ice

! JDL FINAL - Modify wet number concentration to preserve mmdi
  IF (nti > 1.E-8_kind_real) THEN
    dm   = ((rhoa*qi)/(pisix*density_ice*nti))**(1._kind_real/3._kind_real)
    ntm  = (rhoa*fm)/((pisix*rhom*dm**3._kind_real))
  ELSE
    ntm = 0._kind_real
    fm = 0._kind_real
  ENDIF
END SUBROUTINE fractionWater_md

SUBROUTINE power_mom(power,cx,t,rhoa,q,moment)

!-----------------------------------------------------------------------
!
! PURPOSE:
!
! Calculates moments of the PSD based on the Field et al. 2005 power law
! relations. Used for Thompson scheme.
!
!
! AUTHOR:  Bryan Putnam, 4/16/2013
!
!-----------------------------------------------------------------------

  IMPLICIT NONE

!-----------------------------------------------------------------------
! Declare arguments.
!-----------------------------------------------------------------------

  REAL(kind_real),    INTENT(IN   ) :: rhoa
  INTEGER(kind_int), INTENT(IN   ) :: power
  REAL(kind_real),    INTENT(IN   ) :: t,q,cx
  REAL(kind_real),    INTENT(  OUT) :: moment

!-----------------------------------------------------------------------
! Declare local variables.
!-----------------------------------------------------------------------

  REAL(kind_real) :: a,b
  REAL(kind_real) :: rpower  
  REAL(kind_real) :: log_a
  REAL(kind_real) :: second_moment
  REAL(kind_real) :: T_c

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
! Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@


  T_c = t-273.16_kind_real

  SELECT CASE (mp_option)
    CASE(108)

      second_moment = rhoa * (q/cx)

      IF(power == 2) THEN
        moment = second_moment
      ELSE
 
        rpower = REAL(power)

        log_a = dble(5.065339_kind_real-.062659_kind_real*T_c                    &
                     - 3.032362_kind_real*rpower + 0.029469_kind_real*T_c*rpower &
                     - 0.000285_kind_real*(T_c**2._kind_real)                    &
                     + 0.312550_kind_real*(rpower**2._kind_real)                 &
                     + 0.000204_kind_real*(T_c**2._kind_real)*rpower             &
                     + 0.003199_kind_real*T_c*(rpower**2._kind_real)             &
                     + 0.000000_kind_real*(T_c**3._kind_real)                    &
                     - 0.015952_kind_real*(rpower**3._kind_real))

        a = sngl(10.0_kind_real**log_a)

        b = 0.476221_kind_real - 0.015896_kind_real*T_c + 0.165977_kind_real*rpower &
            + 0.007468_kind_real*T_c*rpower - 0.000141_kind_real*(T_c**2._kind_real)&
            + 0.060366_kind_real*(rpower**2._kind_real)                          &
            + 0.000079_kind_real*(T_c**2._kind_real)*rpower                      &
            + 0.000594_kind_real*T_c*(rpower**2._kind_real)                      &
            + 0.000000_kind_real*(T_c**3._kind_real)                             &
            - 0.003577_kind_real*(rpower**3._kind_real)

        moment = a*(second_moment)**b
      END IF

  END SELECT

END SUBROUTINE

SUBROUTINE calc_N0x_mp(rhoa,rhoms,rhomh,rhomg,ntr,nts,ntms,nth,ntmh,ntg, &
                       ntmg,qrf,qsf,fms,qhf,fmh,qgf,fmg)
!-----------------------------------------------------------------------
!
! PURPOSE:
!
! This subroutine calculates intercep parameter based on MP scheme.
!
!
!
! AUTHOR:  Bryan Putnam, 4/16/2013
!
!-----------------------------------------------------------------------

  IMPLICIT NONE

!-----------------------------------------------------------------------
! Declare arguments.
!-----------------------------------------------------------------------

  REAL(kind_real) :: rhoa,rhoms,rhomh,rhomg
  REAL(kind_real) :: ntr,nts,nth,ntg
  REAL(kind_real) :: qrf,qsf,qhf,qgf
  REAL(kind_real) :: fms,fmh,fmg
  REAL(kind_real) :: ntms,ntmg,ntmh  ! JDL ADD

!-----------------------------------------------------------------------
! Declare local variables.
!-----------------------------------------------------------------------
  REAL(kind_real) :: moma,momb
  REAL(kind_real) :: no_value = missing

  REAL(kind_real), PARAMETER :: D0r = 50.e-5_kind_real
  REAL(kind_real), PARAMETER :: R1 = 1.e-12_kind_real
  REAL(kind_real), PARAMETER :: R2 = 1.e-6_kind_real
  REAL(kind_real), PARAMETER :: gonv_min = 1.e4_kind_real
  REAL(kind_real), PARAMETER :: gonv_max = 3.e6_kind_real
  REAL(kind_real), PARAMETER :: bm_g = 3.0_kind_real

  LOGICAL :: L_qr
  REAL(kind_real) :: mvd_r
  REAL(kind_real) :: dble_alfr
  REAL(kind_real) :: lamr
  REAL(kind_real) :: gamma
  REAL(kind_real) :: xslwq,ygra1,zans1
  REAL(kind_real) :: N0_exp,N0_min
  REAL(kind_real) :: rg,am_g,oge1,cgg_1,cgg_2,cgg_3,ogg1,ogg2,ogmg,cge_1
  REAL(kind_real) :: lam_exp,lamg,ilamg

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
! Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

   SELECT CASE (mp_option)
   CASE(9:12,14,109)
     CALL calc_N0x_melt(rhoa,rhoms,rhomh,rhomg,ntr,nts,ntms,nth,ntmh, &
                         ntg,ntmg,qrf,qsf,fms,qhf,fmh,qgf,fmg)
   CASE(106)

    N0r = dble(8.0E06_kind_real)

    N0g = dble(4.0E06_kind_real)
    N0mg = N0g

    N0s = dble(2.0E06_kind_real*exp((.12_kind_real*(273.16_kind_real-ta))))
    N0ms = N0s
   CASE(108)

!     CALL calc_N0x_melt(rhoa,no_value,no_value,no_value,ntr,no_value, &
!                        no_value,no_value,qrf,no_value,no_value,      &
!                        no_value,no_value,no_value,no_value)    ! 18 arugment
      CALL calc_N0x_melt(rhoa,no_value,no_value,no_value,ntr,no_value, &
                         no_value,no_value,no_value,no_value,no_value, qrf, &
                         no_value,no_value,no_value,no_value,no_value,no_value)

     IF(qrf > R1) THEN
       L_qr = .true.
       dble_alfr = dble(alphar)
       lamr = 0.0_kind_real
       CALL cal_lamda(rhoa,qrf,ntr,rhor,dble_alfr,lamr,T_mur)
       mvd_r = (3.0_kind_real + alphar + 0.672_kind_real)/sngl(lamr)
       IF(mvd_r > 2.5e-3_kind_real) THEN
         mvd_r = 2.5e-3_kind_real
       ELSE IF(mvd_r < ((D0r)*(0.75_kind_real))) THEN
         mvd_r =  D0r*0.75_kind_real
       END IF
     ELSE
       L_qr = .false.
       qrf = 0.0_kind_real
     END IF

     IF(qgf > R1) THEN
       rg = qgf * rhoa
     ELSE
       rg = R1
     END IF

     IF((ta < 270.65_kind_real) .and. L_qr .and. (mvd_r > 100.0e-6_kind_real)) THEN
        xslwq = 4.01_kind_real + log10(mvd_r)
     ELSE
        xslwq = 0.01_kind_real
     END IF

     N0_min = gonv_max
     ygra1 = 4.31_kind_real + log10(max(5.e-5_kind_real,rg))
     zans1 = 3.1_kind_real + (100.0_kind_real/(300.0_kind_real*xslwq*ygra1/(10.0_kind_real/xslwq+1.0_kind_real+0.25_kind_real* &
             ygra1)+30.0_kind_real+10.0_kind_real*ygra1))
     N0_exp = 10.0_kind_real**zans1
     N0_exp = MAX(gonv_min,MIN(N0_exp,gonv_max))
     N0_min = MIN(N0_exp,N0_min)
     N0_exp = N0_min
     am_g = c_x(5)
     oge1 = 1._kind_real/(bm_g + 1._kind_real)
     cgg_1 = sngl(gamma(dble(bm_g) + 1.0_kind_real))
     cgg_2 = sngl(gamma(dble(alphag) + 1.0_kind_real))
     cgg_3 = sngl(gamma(dble(bm_g) + dble(alphag) + 1.0_kind_real))
     ogg1 = 1._kind_real/cgg_1
     ogg2 = 1._kind_real/cgg_2
     ogmg = 1._kind_real/bm_g
     cge_1 = alphag + 1.0_kind_real
     lam_exp = (N0_exp*am_g*cgg_1/rg)**oge1
     lamg = lam_exp*(cgg_3*ogg2*ogg1)**ogmg
     N0g = dble(N0_exp/(cgg_2*lam_exp)*lamg**cge_1)

     IF(fmg > R1) THEN
       rg = fmg * rhoa
     ELSE
       rg = R1
     END IF

     N0_min = gonv_max
     ygra1 = 4.31_kind_real + log10(max(5.e-5_kind_real,rg))
     zans1 = 3.1_kind_real + (100.0_kind_real/(300.0_kind_real*xslwq*ygra1/(10.0_kind_real/xslwq+1.0_kind_real+0.25_kind_real* &
             ygra1)+30.0_kind_real+10.0_kind_real*ygra1))
     N0_exp = 10.0_kind_real**zans1
     N0_exp = MAX(gonv_min,MIN(N0_exp,gonv_max))
     N0_min = MIN(N0_exp,N0_min)
     N0_exp = N0_min
     am_g = c_x(5)
     oge1 = 1._kind_real/(bm_g + 1._kind_real)
     cgg_1 = sngl(gamma(dble(bm_g) + 1.0_kind_real))
     cgg_2 = sngl(gamma(dble(alphag) + 1.0_kind_real))
     cgg_3 = sngl(gamma(dble(bm_g) + dble(alphag) + 1.0_kind_real))
     ogg1 = 1._kind_real/cgg_1
     ogg2 = 1._kind_real/cgg_2
     ogmg = 1._kind_real/bm_g
     cge_1 = alphag + 1.0_kind_real
     lam_exp = (N0_exp*am_g*cgg_1/rg)**oge1
     lamg = lam_exp*(cgg_3*ogg2*ogg1)**ogmg
     N0mg = dble(N0_exp/(cgg_2*lam_exp)*lamg**cge_1)

     IF(qsf >= 1.e-14_kind_real) THEN

       CALL  power_mom(2,c_x(4),ta,rhoa,qsf,moma)
       CALL  power_mom(3,c_x(4),ta,rhoa,qsf,momb)

       N0s = ((dble(moma)**4.0_kind_real)/(dble(momb)**3.0_kind_real))*dble(thom_k0)
       N0s2 =((dble(moma)**4.0_kind_real)/(dble(momb)**3.0_kind_real))*dble(thom_k1)*      &
                  ((dble(moma)/dble(momb))**dble(alphas2))

     ELSE
       N0s = dble(3.0E06_kind_real)
       N0s2 = dble(3.0E06_kind_real)
     END IF

   END SELECT

END SUBROUTINE calc_N0x_mp

SUBROUTINE calc_N0x_melt(rhoa,rhoms,rhomh,rhomg,ntr,nts,ntms,nth,ntmh, &
                         ntg,ntmg,qrf,qsf,fms,qhf,fmh,qgf,fmg)

!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Calculate intercept parameter including melting species
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Bryan Putnam
!    04/16/2013.
!
!-----------------------------------------------------------------------

  IMPLICIT NONE
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
  REAL(kind_real), INTENT(IN   )   :: rhoa,rhoms,rhomh,rhomg
  REAL(kind_real), INTENT(INOUT)  :: ntr,nts,nth,ntg
  REAL(kind_real), INTENT(INOUT)  :: qrf,qsf,qhf,qgf,fms,fmh,fmg
  REAL(kind_real), INTENT(INOUT)  :: ntms,ntmh,ntmg
  REAL(kind_real) :: db_alfr,db_alfs,db_alfh,db_alfg
  REAL(kind_real) :: db_mur,db_mus,db_muh,db_mug
  REAL(kind_real) :: pow1,pow2
  REAL(kind_real), PARAMETER :: epsQ  = 1.e-12_kind_real
  REAL(kind_real), PARAMETER :: epsN  = 1.e-3_kind_real
  REAL(kind_real), PARAMETER :: maxN0 = 4.e+37_kind_real


!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  db_alfr = dble(alphar); db_alfs = dble(alphas); db_alfh = dble(alphah);
  db_alfg = dble(alphag)

  db_mur = dble(T_mur); db_mus = dble(T_mus); db_muh = dble(T_muh);
  db_mug = dble(T_mug)

  IF(qrf >= epsQ .AND. ntr >= epsN) THEN
     CALL cal_N0(rhoa,qrf,ntr,rhor,db_alfr,N0r,db_mur)
     N0r = MIN(maxN0,sngl(N0r))
  ELSE
     qrf = 0.0_kind_real; ntr=0.0_kind_real
  ENDIF

  IF(qsf >= epsQ .AND. nts >= epsN) THEN
     CALL cal_N0(rhoa,qsf,nts,rhos,db_alfs,N0s,db_mus)
     N0s = MIN(maxN0,sngl(N0s))
  ELSE
     qsf = 0.0_kind_real; nts=0.0_kind_real
  ENDIF

  IF(fms >= epsQ .AND. ntms >= epsN) THEN
     CALL cal_N0(rhoa,fms,ntms,rhoms,db_alfs,N0ms,db_mus)
     N0ms = MIN(maxN0,sngl(N0ms))
  ELSE
     fms = 0.0_kind_real; ntms = 0.0_kind_real
  ENDIF

  IF(qhf >= epsQ .AND. nth >= epsN) THEN
     CALL cal_N0(rhoa,qhf,nth,rhoh,db_alfh,N0h,db_muh)
     N0h = MIN(maxN0,sngl(N0h))
  ELSE
     qhf = 0.0_kind_real; nth = 0.0_kind_real
  ENDIF

  IF(fmh >= epsQ .AND. ntmh >= epsN) THEN
     CALL cal_N0(rhoa,fmh,ntmh,rhomh,db_alfh,N0mh,db_muh)
     N0mh = MIN(maxN0,sngl(N0mh))
  ELSE
     fmh = 0.0_kind_real; ntmh = 0.0_kind_real
  ENDIF

  IF(qgf >= epsQ .AND. ntg >= epsN) THEN
     CALL cal_N0(rhoa,qgf,ntg,rhog,db_alfg,N0g,db_mug)
     N0g = MIN(maxN0,sngl(N0g))
  ELSE
     qgf = 0.0_kind_real; ntg = 0.0_kind_real
  ENDIF

  IF(fmg >= epsQ .AND. ntmg >= epsN) THEN
     CALL cal_N0(rhoa,fmg,ntmg,rhomg,db_alfg,N0mg,db_mug)
     N0mg = MIN(maxN0,sngl(N0mg))
  ELSE
     fmg = 0.0_kind_real; ntmg = 0.0_kind_real
  ENDIF

END SUBROUTINE calc_N0x_melt

FUNCTION gamma(xx)

!  Modified from "Numerical Recipes"

  IMPLICIT NONE

! PASSING PARAMETERS:
  REAL(kind_real), INTENT(IN   ) :: xx

! LOCAL PARAMETERS:
  REAL(kind_real) :: gamma
  REAL(kind_real) :: ser,stp,tmp,x,y,cof(6)
  INTEGER(kind_int)  :: j


  SAVE cof,stp
  DATA cof,stp/76.180091729471460_kind_real,-86.505320329416770_kind_real,               &
       24.014098240830910_kind_real,-1.2317395724501550_kind_real,.1208650973866179e-2_kind_real,  &
       -.5395239384953e-5_kind_real,2.50662827463100050_kind_real/
  x=xx
  y=x
  tmp=x+5.50_kind_real
  tmp=(x+0.50_kind_real)*log(tmp)-tmp
  ser=1.0000000001900150_kind_real
  do j=1,4
     y=y+1.0_kind_real
     ser=ser+cof(j)/y
  enddo
  gamma=tmp+log(stp*ser/x)
  gamma= exp(gamma)

END FUNCTION gamma

SUBROUTINE cal_N0(rhoa,q,Ntx,rhox,alpha,N0,mu_x)
!
!-----------------------------------------------------------------------
!  PURPOSE:  Calculates intercept parameter and "effective" intercept parameter
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Dan Dawson
!  (02/06/2008)
!
!  MODIFICATION HISTORY:
!
!  (03/26/2008)
!  Recast N0 as a double precision variable, and used double precision for
!  all intermediate calculations.  The calling subroutine should
!  also define it as double precision.  For situations with large alpha,
!  N0 can become very large, and loss of precision can result.
!  Also tweaked the calculation of N0 a bit to avoid overflow, in keeping
!  With Jason Milbrandt's calculation of N0 just before evaporation in
!  the multi-moment code.
!
!  (10/04/2016)
!  Added additional shape parameter mu, i.e. for a gamma distribution of
!  the form N(D) = N0*D^alpha*exp(-lambda*D^3mu).  Setting mu to 1/3
!  retrieves the original formulation for the standard gamma distribution.
!  This value assumes that the mass diameter relationship for the hydrometeor
!  distribution in question is of the form m(D) = c*D^3 and won't be
!  correct for other values of the exponent
!
!-----------------------------------------------------------------------
!  Variable Declarations:
!-----------------------------------------------------------------------
!
  REAL(kind_real),   INTENT(IN   ) :: rhoa,q,Ntx
  REAL(kind_real),   INTENT(IN   ) :: rhox
  REAL(kind_real), INTENT(IN   ) :: alpha,mu_x
  REAL(kind_real), INTENT(  OUT) :: N0

  REAL(kind_real), PARAMETER :: pi = 3.141592_kind_real   ! pi
  REAL(kind_real) :: gamma1, gamma4
  REAL(kind_real) :: gamma

  REAL(kind_real):: lamda

  gamma1 = gamma((1.0_kind_real+alpha)/(3.0_kind_real*mu_x))
  gamma4 = gamma((4.0_kind_real+alpha)/(3.0_kind_real*mu_x))

  IF(rhoa > 0.0_kind_real .and. q > 0.0_kind_real) THEN
    lamda = ((gamma4/gamma1)*dble(pi/6._kind_real*rhox)*dble(Ntx)/(dble(rhoa)*  &
        dble(q)))**mu_x

  ELSE
    lamda = 0.0_kind_real
  END IF


  N0 = 3.0_kind_real*mu_x*dble(Ntx)*lamda**(0.50_kind_real*((1.0_kind_real+alpha)/(3.0_kind_real*mu_x)))*&
              (1.0_kind_real/gamma1)*lamda**(0.50_kind_real*((1.0_kind_real+alpha))/(3.0_kind_real*mu_x))

END SUBROUTINE cal_N0

SUBROUTINE calc_lamda_mp(rhoa,rhoms,rhomh,rhomg,ntr,nts,ntms,nth,ntmh, &
                             ntg,ntmg,qrf,qsf,fms,qhf,fmh,qgf,fmg)

!-----------------------------------------------------------------------
!
! PURPOSE:
!
! Calculate slope parameter for PSD based on MP scheme.
!
!-----------------------------------------------------------------------
!
! AUTHOR:  Bryan Putnam, 4/16/2013
!
!-----------------------------------------------------------------------

  IMPLICIT NONE

!-----------------------------------------------------------------------
! Declare arguments.
!-----------------------------------------------------------------------

  REAL(kind_real), INTENT(IN   ) :: rhoa,rhoms,rhomh,rhomg
  REAL(kind_real), INTENT(IN   ) :: ntr,nts,nth,ntg
  REAL(kind_real), INTENT(IN   ) :: ntms,ntmh,ntmg ! JDL ADD
  REAL(kind_real), INTENT(IN   ) :: qrf,qsf,fms,qhf,fmh,qgf,fmg

!-----------------------------------------------------------------------
! Declare local variables.
!-----------------------------------------------------------------------


  REAL(kind_real) :: dble_alfr,dble_alfs,dble_alfg,dble_alfh
  REAL(kind_real) :: dble_mur,dble_mus,dble_muh,dble_mug
  REAL(kind_real) :: lamr,lams,lamrs,lamh,lamrh,lamg,lamrg
  REAL(kind_real) :: Ntw,Ntd

  REAL(kind_real) :: tem1,tem2

  REAL(kind_real), PARAMETER :: epsQ  = 1.e-12_kind_real ! coped from other

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
! Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  dble_alfr = dble(alphar)
  dble_alfs = dble(alphas)
  dble_alfg = dble(alphag)
  dble_alfh = dble(alphah)
  dble_mur = dble(T_mur)
  dble_mus = dble(T_mus)
  dble_muh = dble(T_muh)
  dble_mug = dble(T_mug)

  if(qrf > max(1e-20_kind_real, epsQ)) then

   Ntw = 0._kind_real
   if(ntr > 0.0_kind_real) then
    Ntw = ntr
    CALL cal_lamda(rhoa,qrf,Ntw,rhor,dble_alfr,lamr,dble_mur)
     lamdar = sngl(lamr)
   else
    !db_N0 = dble(N0r)
    CALL cal_Nt(rhoa,qrf,N0r,c_x(2),dble_alfr,Ntw,dble_mur)
    CALL cal_lamda(rhoa,qrf,Ntw,rhor,dble_alfr,lamr,dble_mur)
    lamdar = sngl(lamr)
   end if
  else
   lamdar = 0.0_kind_real
  end if

  SELECT CASE (mp_option)
  CASE(1:11,14,106,109,110,116)
   if(qsf > max(1e-20_kind_real, epsQ)) then
    Ntd = 0._kind_real
    if (nts > 0.0_kind_real) then
     Ntd = nts
     CALL cal_lamda(rhoa,qsf,Ntd,rhos,dble_alfs,lams,dble_mus)
     lamdas = sngl(lams)
    else
     CALL cal_Nt(rhoa,qsf,N0s,c_x(4),dble_alfs,Ntd,dble_mus)
     CALL cal_lamda(rhoa,qsf,Ntd,rhos,dble_alfs,lams,dble_mus)
     lamdas = sngl(lams)
    end if
   else
    lamdas = 0.0_kind_real
   end if

    if(fms > 0.0_kind_real) then
     Ntw = 0._kind_real
     if(ntms > 0.0_kind_real) then
      Ntw = ntms
      CALL cal_lamda(rhoa,fms,Ntw,rhoms,dble_alfs,lamrs,dble_mus)
      lamdams = sngl(lamrs)
    else
     CALL cal_Nt(rhoa,fms,N0s,c_x(4),dble_alfs,ntw,dble_mus)
     CALL cal_lamda(rhoa,fms,Ntw,rhoms,dble_alfs,lamrs,dble_mus)
     lamdams = sngl(lamrs)
    end if
   else
    lamdams = 0.0_kind_real
   end if

  CASE(108)
   if(qsf > max(1e-20_kind_real, epsQ)) then

    CALL power_mom(2,c_x(4),ta,rhoa,qsf,tem1)
    CALL power_mom(3,c_x(4),ta,rhoa,qsf,tem2)
    lamdas = (tem1/tem2)*thom_lam0
    lamdas2  = (tem1/tem2)*thom_lam1
   else
    lamdas = 0.0_kind_real
    lamdas2 = 0.0_kind_real
   end if

   if(fms > 0.0_kind_real) then
     Ntw = 0._kind_real
     if(ntms > 0.0_kind_real) then
      Ntw = ntms
      CALL cal_lamda(rhoa,fms,Ntw,rhoms,dble_alfs,lamrs,dble_mus)
      lamdams = sngl(lamrs)
    else
     CALL cal_Nt(rhoa,fms,N0ms,c_x(4),dble_alfs,ntw,dble_mus)
     CALL cal_lamda(rhoa,fms,Ntw,rhoms,dble_alfs,lamrs,dble_mus)
     lamdams = sngl(lamrs)
    end if
   end if

  END SELECT

 if(hail_ON == 1) then
   if(qhf > max(1e-20_kind_real, epsQ)) then
    Ntd = 0._kind_real
    if(nth > 0.0_kind_real) then
     Ntd = nth
     CALL cal_lamda(rhoa,qhf,Ntd,rhoh,dble_alfh,lamh,dble_muh)
     lamdah = sngl(lamh)
    else
     CALL cal_Nt(rhoa,qhf,N0h,c_x(6),dble_alfh,Ntd,dble_muh)
     CALL cal_lamda(rhoa,qhf,Ntd,rhoh,dble_alfh,lamh,dble_muh)
     lamdah = sngl(lamh)
    end if
   else
    lamdah = 0.0_kind_real
   end if

   if(fmh > 0._kind_real) then
    Ntw = 0._kind_real
    if(ntmh > 0.0_kind_real) then
     Ntw = ntmh
     CALL cal_lamda(rhoa,fmh,Ntw,rhomh,dble_alfh,lamrh,dble_muh)
     lamdamh = sngl(lamrh)
    else
     CALL cal_Nt(rhoa,fmh,N0mh,c_x(6),dble_alfh,Ntw,dble_muh)
     CALL cal_lamda(rhoa,fmh,Ntw,rhomh,dble_alfh,lamrh,dble_muh)
     lamdamh = sngl(lamrh)
    end if
   else
    lamdamh = 0.0_kind_real
   end if
 end if

 if(graupel_ON == 1) then

   if(qgf > max(1e-20_kind_real, epsQ)) then
    Ntd = 0._kind_real
    if(ntg > 0.0_kind_real) then
     Ntd = ntg
     CALL cal_lamda(rhoa,qgf,Ntd,rhog,dble_alfg,lamg,dble_mug)
     lamdag = sngl(lamg)
    else
     CALL cal_Nt(rhoa,qgf,N0g,c_x(5),dble_alfg,Ntd,dble_mug)
     CALL cal_lamda(rhoa,qgf,Ntd,rhog,dble_alfg,lamg,dble_mug)
     lamdag = sngl(lamg)
    end if
  else
   lamdag = 0.0_kind_real
  end if

   if(fmg > 0._kind_real) then
    Ntw = 0._kind_real
    if(ntmg > 0.0_kind_real) then
     Ntw = ntmg
     CALL cal_lamda(rhoa,fmg,Ntw,rhomg,dble_alfg,lamrg,dble_mug)
     lamdamg = sngl(lamrg)
    else
     CALL cal_Nt(rhoa,fmg,N0mg,c_x(5),dble_alfg,Ntw,dble_mug)
     CALL cal_lamda(rhoa,fmg,Ntw,rhomg,dble_alfg,lamrg,dble_mug)
     lamdamg = sngl(lamrg)
    end if
   else
    lamdamg = 0.0_kind_real
   end if
  end if

END SUBROUTINE calc_lamda_mp

SUBROUTINE cal_Nt(rhoa,q,N0,cx,alpha,Ntx,mu_x)
!
!-----------------------------------------------------------------------
!  PURPOSE:  Calculates number concentration at scalar points
!-----------------------------------------------------------------------
!
!  AUTHOR: Dan Dawson
!  (02/06/2008)
!
!  MODIFICATION HISTORY:
!
!  03/31/08 - converted intermediate calculations to double precision
!             as well as a few of the input arguments.
!
!  (10/04/2016)
!  Added additional shape parameter mu, i.e. for a gamma distribution of
!  the form N(D) = N0*D^alpha*exp(-lambda*D^3mu).  Setting mu to 1/3
!  retrieves the original formulation for the standard gamma distribution.
!  This value assumes that the mass diameter relationship for the hydrometeor
!  distribution in question is of the form m(D) = c*D^3 and won't be
!  correct for other values of the exponent.
!
!-----------------------------------------------------------------------
!  Variable Declarations:
!-----------------------------------------------------------------------
!
  REAL(kind_real),  INTENT(IN   ) :: rhoa,q
  REAL(kind_real),INTENT(IN   ) :: alpha,N0,mu_x
  REAL(kind_real),  INTENT(IN   ) :: cx
  REAL(kind_real),  INTENT(  OUT) :: Ntx

  REAL(kind_real) :: gamma1,gamma4
  REAL(kind_real) :: gamma

  gamma1 = gamma((1.0_kind_real+alpha)/(3.0_kind_real*mu_x))
  gamma4 = gamma((4.0_kind_real+alpha)/(3.0_kind_real*mu_x))

  Ntx = sngl((N0*gamma1/(3.0_kind_real*mu_x))**(3.0_kind_real/(4.0_kind_real+alpha))*   &
                   ((gamma1/gamma4)*dble(rhoa)* &
                   dble(q)/dble(cx))**((1.0_kind_real+alpha)/(4.0_kind_real+alpha)))

END SUBROUTINE 

SUBROUTINE cal_lamda(rhoa,q,Ntx,rhox,alpha,lamda,mu_x)
!
!-----------------------------------------------------------------------
!  PURPOSE:  Calculates slope parameter lamda
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Dan Dawson!  (02/06/2008)
!
!  MODIFICATION HISTORY:
!  (03/31/2008)
!  Converted intermediate calculations and arrays alpha and lamda to
!  double precision.
!
!  (10/04/2016)
!  Added additional shape parameter mu, i.e. for a gamma distribution of
!  the form N(D) = N0*D^alpha*exp(-lambda*D^3mu).  Setting mu to 1/3
!  retrieves the original formulation for the standard gamma distribution.
!  This value assumes that the mass diameter relationship for the hydrometeor
!  distribution in question is of the form m(D) = c*D^3 and won't be
!  correct for other values of the exponent.
!
!-----------------------------------------------------------------------
!  Variable Declarations:
!-----------------------------------------------------------------------
!

  REAL(kind_real),  INTENT(IN   ) :: rhoa,q
  REAL(kind_real),INTENT(IN   ) :: alpha,mu_x
  REAL(kind_real),  INTENT(IN   ) :: Ntx
  REAL(kind_real),  INTENT(IN   ) :: rhox
  REAL(kind_real),INTENT(  OUT) :: lamda
  REAL(kind_real), PARAMETER :: pi = 3.141592_kind_real   ! pi

  REAL(kind_real) :: gamma1, gamma4
  REAL(kind_real) :: gamma

  gamma1 = gamma((1.0_kind_real+alpha)/(3.0_kind_real*mu_x))
  gamma4 = gamma((4.0_kind_real+alpha)/(3.0_kind_real*mu_x))


  IF(rhoa > 0.0_kind_real .and. q > 0.0_kind_real) THEN
    lamda = ((gamma4/gamma1)*dble(pi/6._kind_real*rhox)*dble(Ntx)/(dble(rhoa)*  &
          dble(q)))**mu_x
  ELSE
    lamda = 0.0_kind_real
  END IF

END SUBROUTINE cal_lamda

!
!########################################################################
!########################################################################
!#########                                                      #########
!#########               SUBROUTINE RSET_DSD_PARA               #########
!#########                                                      #########
!#########                     Developed by                     #########
!#########     Center for Analysis and Prediction of Storms     #########
!#########                University of Oklahoma                #########
!#########                                                      #########
!########################################################################
!########################################################################

SUBROUTINE set_dsd_para()

!-----------------------------------------------------------------------
!
! PURPOSE:
!
! This subroutine sets intercept parameters for rain/snow/hail and
! densities for snow/hail based on values in history dump.
!
!-----------------------------------------------------------------------
!
! AUTHOR:  Youngsun Jung, Spring 2010
!
!-----------------------------------------------------------------------
! Include global variables only for dual-pol calculations
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
! Force explicit declarations.
!-----------------------------------------------------------------------

  IMPLICIT NONE

!-----------------------------------------------------------------------
! Include files.
!-----------------------------------------------------------------------

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
! Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  CALL model_dsd(n0rain,n0snow,n0hail,n0grpl,rhosnow,rhohail,rhogrpl,  &
       alpharain,alphasnow,alphagrpl,alphahail)

  IF (rhos <= 0.0_kind_real) THEN
    rhos = 100._kind_real
  END IF

  IF (rhoh <= 0.0_kind_real) THEN
    SELECT CASE (mp_option)
    CASE(1:12,106,108:110,116)
      rhoh = 913._kind_real
    CASE(13:14)
      rhoh = 900._kind_real
    END SELECT
  END IF

  IF (rhog <= 0.0_kind_real) THEN

    SELECT CASE (mp_option)
    CASE(1:4,8:12,108:110)
      rhog = 400._kind_real
    CASE(5:7,13:14,106,116)
      rhog = 500._kind_real
    END SELECT
  END IF

  IF (N0r <= 0.0_kind_real) THEN
    N0r = 8.0E+06_kind_real
    SELECT CASE (mp_option)
    CASE(13:14)
       N0r = 8.0E+05_kind_real
    END SELECT
  END IF

  IF (N0s <= 0.0) THEN
    N0s = 3.0E+06_kind_real
  ENDIF

  SELECT CASE (mp_option)
  CASE(1:14,106,109:110,116)
    N0s2 = 0.0_kind_real
  CASE(108)
    N0s2 = 3.0E+06_kind_real
  END SELECT

  IF (N0h <= 0.0_kind_real) THEN
    N0h = 4.0E+04_kind_real
  END IF

  IF (N0g <= 0.0_kind_real) THEN
    SELECT CASE (mp_option)
    CASE(1:4,8:14,108:110)
      N0g = 4.0E+05_kind_real
    CASE(5:7,106,116)
      N0g = 4.0E+06_kind_real
    END SELECT
  END IF

  N0ms = N0s
  N0ms2 = N0s2
  N0mh = N0h
  N0mg = N0g

  IF (alphar <= 0.0_kind_real) THEN
    SELECT CASE (mp_option)
      CASE(1:12,106,108:110)
        alphar = 0.0_kind_real
    END SELECT
  END IF

  IF (alphas <= 0.0_kind_real) THEN
    alphas = 0.0_kind_real
  END IF

  SELECT CASE (mp_option)
    CASE(1:12,106,109,110,116)
      alphas2 = 0.0_kind_real
    CASE(108)
      alphas2 = 0.6357_kind_real
  END SELECT

  IF (alphah <= 0.0_kind_real) THEN
    alphah = 0.0_kind_real
  END IF

  IF (alphag <= 0.0_kind_real) THEN
     alphag = 0.0_kind_real
  END IF

  lamdar = 0.0_kind_real
  lamdas = 0.0_kind_real
  lamdas2 = 0.0_kind_real
  lamdams = 0.0_kind_real
  lamdams2 = 0.0_kind_real
  lamdag = 0.0_kind_real
  lamdamg = 0.0_kind_real
  lamdah = 0.0_kind_real
  lamdamh = 0.0_kind_real

  RETURN
END SUBROUTINE set_dsd_para

!########################################################################
!########################################################################
!#########                                                      #########
!#########               SUBROUTINE rdr_obs                     #########
!#########                                                      #########
!#########                     Developed by                     #########
!#########     Center for Analysis and Prediction of Storms     #########
!#########                University of Oklahoma                #########
!#########                                                      #########
!########################################################################
!########################################################################

SUBROUTINE rdr_obs (rho,qscalar,obs_dual,var_dsd,    &
                       var_idx,dualpol)

!-----------------------------------------------------------------------
!
! PURPOSE:
!
! A shell subroutine to assign DSD parameters for the simulated
! radar parameters using parameterized formula based on Jung et al.(2008a).
!
!-----------------------------------------------------------------------
!
! AUTHOR:  Youngsun Jung, 12/14/2010
!
! MODIFICATION HISTORY:
!
!  Bryan Putnam 4/16/2013: Added in information for all radar parameters and
!  all operators, replaces rdr_obs_SM.
!
!-----------------------------------------------------------------------
! Include global variables only for dual-pol calculations
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
! Force explicit declarations.
!-----------------------------------------------------------------------

  IMPLICIT NONE

!-----------------------------------------------------------------------
! Declare arguments.
!-----------------------------------------------------------------------
  REAL(kind_real), INTENT(IN   ) :: qscalar(nscalar)
  REAL(kind_real), INTENT(IN   ) :: rho

  INTEGER(kind_int), INTENT(IN   ) :: var_idx,dualpol

  TYPE(T_obs_dual), INTENT(INOUT) :: obs_dual
  TYPE(T_para_dsd), INTENT(INOUT) :: var_dsd

 !local variables
  REAL(kind_real) :: no_value = missing

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
! Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  SELECT CASE (mp_option)
  CASE(2:8,106)  ! single moment schemes 
    SELECT CASE (qgh_opt)
      CASE (1)                       ! graupel off, hail off
        var_dsd = assign_para_dsd_TM(qscalar(P_qr),qscalar(P_qs),no_value, &
              no_value, no_value, no_value, no_value, no_value, alphar,    &
              alphas,no_value,no_value,no_value,no_value)
      CASE (2)                       ! graupel off, hail on
        var_dsd = assign_para_dsd_TM(qscalar(P_qr),qscalar(P_qs),        &
                   qscalar(P_qh),no_value,no_value,no_value,no_value,  &
                   no_value,alphar,alphas,alphah,no_value,no_value,no_value)
      CASE (3)                       ! graupel on, hail off
        var_dsd = assign_para_dsd_TM(qscalar(P_qr),qscalar(P_qs),no_value, &
                   qscalar(P_qg),no_value,no_value,no_value,no_value,    &
                   alphar,alphas,no_value,alphag,no_value,no_value)
      CASE (4)                       ! graupel on, hail on
        var_dsd = assign_para_dsd_TM(qscalar(P_qr),qscalar(P_qs),        &
                  qscalar(P_qh),qscalar(P_qg),no_value,no_value,        &
                  no_value,no_value,alphar,alphas,alphah,alphag,no_value,no_value)
    END SELECT
  CASE(108,116) ! double moment for rain only
    SELECT CASE (qgh_opt)
      CASE(1)
        var_dsd = assign_para_dsd_TM(qscalar(P_qr),qscalar(P_qs),        &
                   no_value,no_value,qscalar(P_nr),no_value,no_value,    &
                   no_value,alphar,alphas,no_value,no_value,no_value,no_value)
      CASE(2)
        var_dsd = assign_para_dsd_TM(qscalar(P_qr),qscalar(P_qs),        &
                    qscalar(P_qh),no_value,qscalar(P_nr),no_value,     &
                    no_value,no_value,alphar,alphas,         &
                    alphah,no_value,no_value,no_value)
      CASE(3)
        var_dsd = assign_para_dsd_TM(qscalar(P_qr),qscalar(P_qs),        &
                   no_value,qscalar(P_qg),qscalar(P_nr),no_value,      &
                   no_value,no_value,alphar,alphas,no_value,  &
                   alphag,no_value,no_value)
      CASE(4)
       var_dsd = assign_para_dsd_TM(qscalar(P_qr),qscalar(P_qs),         &
                   qscalar(P_qh),qscalar(P_qg),qscalar(P_nr),        &
                   no_value,no_value,no_value,alphar,alphas,  &
                   alphah,alphag,no_value,no_value)
    END SELECT

  CASE(14)  ! NSSL Scheme also accounts for Volume
    var_dsd =  assign_para_dsd_TM(qscalar(P_QR),qscalar(P_QS),        &
                   qscalar(P_QH),qscalar(P_QG),qscalar(P_NR),        &
                   qscalar(P_NS),qscalar(P_NH),qscalar(P_NG),        &
                   alphar,alphas,alphah,alphag,qscalar(P_VG),        &
                   qscalar(P_VH))

  END SELECT

  dualpol_opt = dualpol
  IF(dualpol == 1) THEN
     IF(var_idx <= 3) THEN
        obs_dual = calculate_obs(rho,var_dsd,var_idx)
     END IF
   END IF

  RETURN

END SUBROUTINE rdr_obs

END MODULE RADARZ_MODULE
