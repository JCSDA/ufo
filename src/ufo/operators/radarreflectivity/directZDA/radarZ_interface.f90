MODULE radarz_iface
!
! !DESCRIPTION: initialize variables required for radarz library
!               required for direct reflectivity DA capabilities
!
! !REVISION HISTORY:
!   2019-xx-xx  CAPS - initial commit
!   2021-05-17  J. Park(CAPS) - radaremul renamed to radarz
!                             - deleted unnecessary mphyopt options
!   2021-06-30  J. Park(CAPS) - refactor this for UFO
!   2021-10-20  J. Park(CAPS) - update radarZ with the latest GSI code
!

  use kinds, only: kind_real,kind_int

  IMPLICIT NONE
  SAVE
  PRIVATE

  PUBLIC :: init_mphyopt
  PUBLIC :: get_qgh_opt
  PUBLIC :: T_para_dsd
  PUBLIC :: T_obs_dual

  real(kind_real),PUBLIC,PARAMETER :: PI = 3.1415926536
  REAL(kind_real),PUBLIC,PARAMETER :: degKtoC = 273.15_kind_real ! Convert Kelvin to Centigrade
  REAL(kind_real),PUBLIC, PARAMETER :: missing = -9999.0_kind_real
  REAL(kind_real),PUBLIC :: ta = 273.16_kind_real

  integer(kind_int), public :: mp_option = 0
  integer(kind_int), public :: MFflg = 2
  integer(kind_int), public :: hail_ON = 0, graupel_ON = 0, qgh_opt = 0
  integer(kind_int), public :: nscalar = 0
  REAL(kind_real), public :: grpl_miss
  REAL(kind_real), public :: hl_miss 

  ! Place/index markers within arrays for each species mixing ratio and number, etc.
  integer(kind_int), public :: P_qc=0, P_qr=0, P_qi=0, P_qs=0, P_qh=0, P_qg=0,    &
     &                         P_nc=0, P_nr=0, P_ni=0, P_ns=0, P_nh=0, P_ng=0,    &
     &                                 P_zr=0, P_zi=0, P_zs=0, P_zh=0, P_zg=0,    &
     &                                                 P_vh=0, P_vg=0

  ! Default values of Y-intercept parameters, gamma shape parameters, and density
  real(kind_real), public :: n0rain = 8.0e6_kind_real
  real(kind_real), public :: n0snow = 3.0e6_kind_real
  real(kind_real), public :: n0hail = 4.0e4_kind_real
  real(kind_real), public :: n0grpl = 4.0e5_kind_real
  real(kind_real), public :: alpharain = 0.0_kind_real
  real(kind_real), public :: alphasnow = 0.0_kind_real
  real(kind_real), public :: alphahail = 0.0_kind_real
  real(kind_real), public :: alphagrpl = 0.0_kind_real
  real(kind_real), public :: rhosnow = 100.0_kind_real
  real(kind_real), public :: rhohail = 913.0_kind_real
  real(kind_real), public :: rhogrpl = 500.0_kind_real
  REAL(kind_real),PUBLIC,PARAMETER :: rhoi = 917._kind_real  ! Density of ice (kg m**-3)
 
  ! A few specific snow distribution values for Thompson MP
  REAL(kind_real),PUBLIC,PARAMETER :: thom_lam0 = 20.78_kind_real
  REAL(kind_real),PUBLIC,PARAMETER :: thom_lam1 = 3.29_kind_real
  REAL(kind_real),PUBLIC,PARAMETER :: thom_k0 = 490.6_kind_real
  REAL(kind_real),PUBLIC,PARAMETER :: thom_k1 = 17.46_kind_real

  ! A few specific distribution values for NSSL MP
  integer(kind_int), public :: murain = 1 ! NSSL default
  integer(kind_int), public :: musnow = 3 ! NSSL default

  integer(kind_int), public :: dsdparaopt = 0
  integer(kind_int), public :: nen = 1

  ! Precalculated complete gamma function values
  REAL(kind_real),PUBLIC,PARAMETER :: gamma7_08 = 836.7818_kind_real
  REAL(kind_real),PUBLIC,PARAMETER :: gamma6_81 = 505.8403_kind_real
  REAL(kind_real),PUBLIC,PARAMETER :: gamma6_54 = 309.3308_kind_real
  REAL(kind_real),PUBLIC,PARAMETER :: gamma5_63 = 64.6460_kind_real
  REAL(kind_real),PUBLIC,PARAMETER :: gamma4_16 = 7.3619_kind_real
  REAL(kind_real),PUBLIC,PARAMETER :: gamma3_97 = 5.7788_kind_real

  ! Variables to can be changed by parameter retrieval
  REAL(kind_real), PUBLIC :: N0r        ! Intercept parameter in 1/(m^4) for rain
  REAL(kind_real), PUBLIC :: N0s        ! Intercept parameter in 1/(m^4) for snow
  REAL(kind_real), PUBLIC :: N0h        ! Intercept parameter in 1/(m^4) for hail
  REAL(kind_real), PUBLIC :: N0g        ! Intercept parameter in 1/(m^4) for hail
  REAL(kind_real), PUBLIC :: N0s2       ! Second intercept parameter in 1/(m^4) for snow

  REAL(kind_real), PUBLIC :: N0ms       !Intercept parameter for melting species 
  REAL(kind_real), PUBLIC :: N0ms2 
  REAL(kind_real), PUBLIC :: N0mh
  REAL(kind_real), PUBLIC :: N0mg 

  REAL(kind_real), PUBLIC :: rhor = 1000._kind_real ! Density of rain (kg m**-3)
  REAL(kind_real), PUBLIC :: rhoh                ! Density of hail (kg m**-3)
  REAL(kind_real), PUBLIC :: rhos                ! Density of snow (kg m**-3)
  REAL(kind_real), PUBLIC :: rhog                ! Density of graupel (kg m**-3)

  REAL(kind_real), PUBLIC :: alphar     !Shape parameter for rain
  REAL(kind_real), PUBLIC :: alphas     !Shape parameter for snow
  REAL(kind_real), PUBLIC :: alphah     !SHape parameter for hail
  REAL(kind_real), PUBLIC :: alphag     !SHape parameter for graupel
  REAL(kind_real), PUBLIC :: alphas2

  REAL(kind_real), PUBLIC :: lamdar     !slope parameter for rain (1/m)
  REAL(kind_real), PUBLIC :: lamdas     
  REAL(kind_real), PUBLIC :: lamdas2  
  REAL(kind_real), PUBLIC :: lamdams
  REAL(kind_real), PUBLIC :: lamdams2  
  REAL(kind_real), PUBLIC :: lamdag
  REAL(kind_real), PUBLIC :: lamdamg
  REAL(kind_real), PUBLIC :: lamdah
  REAL(kind_real), PUBLIC :: lamdamh 

  ! Variables to can be changed for meling ice
  REAL(kind_real), PUBLIC :: fos        ! Maximum fraction of rain-snow mixture
  REAL(kind_real), PUBLIC :: foh        ! Maximum fraction of rain-hail mixture
  REAL(kind_real), PUBLIC :: fog        ! Maximum fraction of rain-hail mixture

  ! Additional shape parameter that accomodate gamma-in-volume 
  REAL(kind_real), PUBLIC :: T_mur,T_mus,T_mug,T_muh

  ! Constants related to radar parameters

  integer(kind_int), public, parameter :: rfopt = 1
  integer(kind_int), public :: attn_ON
  integer(kind_int), public :: dualpol_opt
  real(kind_real), public, parameter :: wavelen = 107.0_kind_real
  real(kind_real), public, parameter :: Kw2 = 0.93_kind_real  ! Dielectric factor for water.
  real(kind_real), public, parameter :: alphak = 3.88e-4_kind_real   ! differential forward scattering (rain)
  real(kind_real), public, parameter :: radar_const = (4._kind_real * wavelen**4._kind_real)/(PI**4 * Kw2)
  real(kind_real), public, parameter :: constKdpr = 180._kind_real * wavelen  * alphak * 1.0e6_kind_real / PI
  real(kind_real), public, parameter :: kdpCoefIce = (180 * wavelen * 1.e6_kind_real) / PI

  REAL(kind_real),PUBLIC,PARAMETER :: alphaa = 4.28e-4_kind_real   ! backscattering amplitude constant (rain)
  REAL(kind_real),PUBLIC,PARAMETER :: alphask = 8.53e-7_kind_real  ! differential forward scattering (snow)
  REAL(kind_real),PUBLIC,PARAMETER :: beta_ra = 3.04_kind_real
  REAL(kind_real),PUBLIC,PARAMETER :: alphab = 4.28e-4_kind_real   ! backscattering amplitude constant (rain)
  REAL(kind_real),PUBLIC,PARAMETER :: beta_rb = 2.77_kind_real
  REAL(kind_real),PUBLIC,PARAMETER :: alphaa_ds = 1.94e-5_kind_real ! for dry snow at horz plane
  REAL(kind_real),PUBLIC,PARAMETER :: alphab_ds = 1.91e-5_kind_real ! for dry snow at vert plane
  REAL(kind_real),PUBLIC,PARAMETER :: alphaa_tom_ds = 2.8e-5_kind_real !for dry snow at horz plane for Thomposon scheme
  REAL(kind_real),PUBLIC,PARAMETER :: alphab_tom_ds = 2.6e-5_kind_real !for dry snow at vert plane for the Thompson scheme 

  REAL(kind_real), PUBLIC, PARAMETER :: beta_sa = 3.0_kind_real
  REAL(kind_real), PUBLIC, PARAMETER :: beta_sb = 3.0_kind_real
  REAL(kind_real), PUBLIC, PARAMETER :: beta_tom_dsa = 1.95_kind_real  ! Special expon for Thompson scheme
  REAL(kind_real), PUBLIC, PARAMETER :: beta_tom_dsb = 1.965_kind_real ! Special expon for Thompson scheme

  REAL(kind_real),PUBLIC,PARAMETER :: alphaa_dh = 1.91e-4_kind_real ! for dry hail at horz plane
  REAL(kind_real),PUBLIC,PARAMETER :: alphab_dh = 1.65e-4_kind_real ! for dry hail at vert plane

  REAL(kind_real),PUBLIC,PARAMETER :: beta_ha = 3.0_kind_real
  REAL(kind_real),PUBLIC,PARAMETER :: beta_hb = 3.0_kind_real

  REAL(kind_real),PUBLIC,PARAMETER :: alphaa_dg = 0.81e-4_kind_real ! for dry graupel at horz plane
  REAL(kind_real),PUBLIC,PARAMETER :: alphab_dg = 0.76e-4_kind_real ! for dry graupel at vert plane

  REAL(kind_real),PUBLIC,PARAMETER :: beta_ga = 3.0_kind_real
  REAL(kind_real),PUBLIC,PARAMETER :: beta_gb = 3.0_kind_real

  REAL(kind_real),PUBLIC,PARAMETER :: alphak_ds = 0.03e-5_kind_real     ! alphaa_ds - alphab_ds
  REAL(kind_real),PUBLIC,PARAMETER :: alphak_tom_ds = 1.05e-6_kind_real !alphaa_ds - alphab_ds for Thompson scheme
  REAL(kind_real),PUBLIC,PARAMETER :: alphak_dh = 0.26e-4_kind_real ! alphaa_dh - alphab_dh
  REAL(kind_real),PUBLIC,PARAMETER :: alphak_dg = 0.05e-4_kind_real ! alphaa_dh - alphab_dh
  REAL(kind_real),PUBLIC,PARAMETER :: betak_s = 3.0_kind_real
  REAL(kind_real),PUBLIC,PARAMETER :: betak_tom_ds = 2.04_kind_real !For Thompson Scheme 
  REAL(kind_real),PUBLIC,PARAMETER :: betak_h = 3.0_kind_real
  REAL(kind_real),PUBLIC,PARAMETER :: betak_g = 3.0_kind_real

  REAL(kind_real),PUBLIC,PARAMETER :: rho_0r = 1.0_kind_real      ! rho_0 for rain
  REAL(kind_real),PUBLIC,PARAMETER :: rho_0s = 1.0_kind_real      ! rho_0 for snow
  REAL(kind_real),PUBLIC,PARAMETER :: rho_0h = 0.97_kind_real     ! rho_0 for hail
  REAL(kind_real),PUBLIC,PARAMETER :: rho_0g = 0.95_kind_real     ! rho_0 for hail
  REAL(kind_real),PUBLIC,PARAMETER :: rho_0rsi = 0.82_kind_real   ! lower limit of rho_0rs (rain-snow mixture)
  REAL(kind_real),PUBLIC,PARAMETER :: rho_0rsf = 0.95_kind_real   ! upper limit of rho_0rs (rain-snow mixture)
  REAL(kind_real),PUBLIC,PARAMETER :: rho_0rhi = 0.85_kind_real   ! lower limit of rho_0rh (rain-hail mixture)
  REAL(kind_real),PUBLIC,PARAMETER :: rho_0rhf = 0.95_kind_real   ! upper limit of rho_0rh (rain-hail mixture)
  REAL(kind_real),PUBLIC,PARAMETER :: rho_0rgi = 0.82_kind_real   ! lower limit of rho_0rg (rain-graupel mixture)
  REAL(kind_real),PUBLIC,PARAMETER :: rho_0rgf = 0.95_kind_real   ! upper limit of rho_0rg (rain-graupel mixture)

  REAL(kind_real),PUBLIC,PARAMETER :: mm3todBZ = 1.0E+9_kind_real ! Conversion factor from mm**3 to
                                                     !   mm**6 m**-3.

  REAL(kind_real),PUBLIC,PARAMETER :: unit_factor = 1.e-2_kind_real  ! Unit conversion factor not addressed
                                         ! in the T-matrix scattering amplitude (size D is in cm in T-matrix)

  REAL(kind_real), dimension(6), PUBLIC :: c_x   !(PI/6)*rho_qx

!-----------------------------------------------------------------------
! Scattering matrix coefficient for snow
!
! phi=0._kind_real       (Mean orientation)
! sigmas=pi/9
! As=1/8*(3+4*cos(2*phi)*exp(-2*sigmas**2)+cos(4*phi)*exp(-8*sigmas**2))
! Bs=1/8*(3-4*cos(2*phi)*exp(-2*sigmas**2)+cos(4*phi)*exp(-8*sigmas**2))
! Cs=1/8*(1-cos(4*phi)*exp(-8*sigmas**2))
! Ds=1/8*(3+cos(4*phi)*exp(-8*sigmas**2))
! Cks=cos(2*phi)*exp(-2*sigmas**2)
!-----------------------------------------------------------------------

  REAL(kind_real),PUBLIC,PARAMETER :: sigmas = 0.3491_kind_real
  REAL(kind_real),PUBLIC,PARAMETER :: As = 0.8140_kind_real
  REAL(kind_real),PUBLIC,PARAMETER :: Bs = 0.0303_kind_real
  REAL(kind_real),PUBLIC,PARAMETER :: Cs = 0.0778_kind_real
  REAL(kind_real),PUBLIC,PARAMETER :: Ds = 0.4221_kind_real
  REAL(kind_real),PUBLIC,PARAMETER :: Cks = 0.7837_kind_real

!-----------------------------------------------------------------------
! Scattering matrix coefficient for hail
!
! phi=0._kind_real     (Mean orientation)
! sigmah=pi/3*(1-sf*fw)
! Ah=1/8*(3+4*cos(2*phi)*exp(-2*sigmah**2)+cos(4*phi)*exp(-8*sigmah**2))
! Bh=1/8*(3-4*cos(2*phi)*exp(-2*sigmah**2)+cos(4*phi)*exp(-8*sigmah**2))
! Ch=1/8*(1-cos(4*phi)*exp(-8*sigmah**2))
! Dh=1/8*(3+cos(4*phi)*exp(-8*sigmah**2))
! Ckh=cos(2*phi)*exp(-2*sigmah**2)
!
! corresponding coefficient for dry hail: Ahd, Bhd, Chd, Dhd, Ckhd
!-----------------------------------------------------------------------

  REAL(kind_real),PUBLIC,PARAMETER :: sigmahd = 1.0472_kind_real
  REAL(kind_real),PUBLIC,PARAMETER :: Ahd = 0.4308_kind_real
  REAL(kind_real),PUBLIC,PARAMETER :: Bhd = 0.3192_kind_real
  REAL(kind_real),PUBLIC,PARAMETER :: Chd = 0.1250_kind_real
  REAL(kind_real),PUBLIC,PARAMETER :: Dhd = 0.3750_kind_real
  REAL(kind_real),PUBLIC,PARAMETER :: Ckhd = 0.1116_kind_real
  REAL(kind_real),PUBLIC :: sf
  REAL(kind_real),PUBLIC :: sigmah, Ah, Bh, Ch, Dh, Ckh

!-----------------------------------------------------------------------
! Scattering matrix coefficient for graupel
! 
! phi=0._kind_real     (Mean orientation)
! sigmag=pi/3*(1-sf*fw)
! Ag=1/8*(3+4*cos(2*phi)*exp(-2*sigmag**2)+cos(4*phi)*exp(-8*sigmag**2))
! Bg=1/8*(3-4*cos(2*phi)*exp(-2*sigmag**2)+cos(4*phi)*exp(-8*sigmag**2))
! Cg=1/8*(1-cos(4*phi)*exp(-8*sigmag**2))
! Dg=1/8*(3+cos(4*phi)*exp(-8*sigmag**2))
! Ckg=cos(2*phi)*exp(-2*sigmag**2)
! 
! corresponding coefficient for dry graupel: Agd, Bgd, Cgd, Dgd, Ckgd
!-----------------------------------------------------------------------
  
  REAL(kind_real),PUBLIC,PARAMETER :: sigmagd = 1.0472_kind_real
  REAL(kind_real),PUBLIC,PARAMETER :: Agd = 0.4308_kind_real
  REAL(kind_real),PUBLIC,PARAMETER :: Bgd = 0.3192_kind_real
  REAL(kind_real),PUBLIC,PARAMETER :: Cgd = 0.1250_kind_real
  REAL(kind_real),PUBLIC,PARAMETER :: Dgd = 0.3750_kind_real
  REAL(kind_real),PUBLIC,PARAMETER :: Ckgd = 0.1116_kind_real
  REAL(kind_real),PUBLIC :: sigmag, Ag, Bg, Cg, Dg, Ckg

!-----------------------------------------------------------------------
!  Declare new observation type
!-----------------------------------------------------------------------

  TYPE T_obs_dual
       REAL(kind_real) :: T_log_ref, T_sum_ref_h, T_sum_ref_v
       REAL(kind_real) :: T_log_zdr, T_sum_ref_hv, T_kdp
       REAL(kind_real) :: T_Ahh,     T_Avv
       REAL(kind_real) :: T_ref_r_h, T_ref_s_h, T_ref_h_h,T_ref_g_h
       REAL(kind_real) :: T_ref_rs_h,T_ref_rh_h,T_ref_rg_h
       REAL(kind_real) :: T_ref_r_v, T_ref_s_v, T_ref_h_v, T_ref_g_v
       REAL(kind_real) :: T_ref_rs_v, T_ref_rh_v, T_ref_rg_v
  END TYPE T_obs_dual

!-----------------------------------------------------------------------
!  Declare new DSD parameter data type
!-----------------------------------------------------------------------

  TYPE T_para_dsd
    REAL(kind_real) :: T_qr, T_qs, T_qh, T_qg, T_vg, T_vh
    REAL(kind_real) :: T_Ntr, T_Nts, T_Nth, T_Ntg
    REAL(kind_real) :: T_alfr,T_alfs,T_alfh,T_alfg
  END TYPE T_para_dsd


CONTAINS

!+---+-----------------------------------------------------------------+

integer(kind_int) function init_mphyopt(mp_option)

  IMPLICIT NONE

  integer, intent(in) :: mp_option
  integer:: iret = -1

  IF ( mp_option == 2 .OR. mp_option == 3 .OR. mp_option == 4 ) THEN
    nscalar = 5
    P_qc = 1; P_qr = 2; P_qi = 3; P_qs = 4; P_qh = 5

    graupel_ON = 0
    hail_ON = 1
  ELSE IF ( mp_option == 5 .OR. mp_option == 6 .OR. mp_option == 7 ) THEN
    nscalar = 5
    P_qc = 1; P_qr = 2; P_qi = 3; P_qs = 4; P_qg = 5

    graupel_ON = 1
    hail_ON = 0
  ELSE IF (mp_option == 102 .OR. mp_option == 106) THEN  ! linscheme, wsm6scheme
    P_qc = 1; P_qr = 2; P_qi = 3; P_qs = 4; P_qg = 5
    nscalar  = 5

    graupel_ON = 1
    hail_ON = 0
  ELSE IF (mp_option == 14) THEN                      ! NSSL
    P_qc = 1; P_qr =  2; P_qi =  3; P_qs =  4; P_qg =  5; P_qh =  6;
    P_nc = 7; P_nr =  8; P_ni =  9; P_ns = 10; P_ng = 11; P_nh = 12;
                                               P_vg = 13; P_vh = 14;
    nscalar = 14

    graupel_ON = 1
    hail_ON = 1
  ELSE IF (mp_option == 108 ) THEN  ! thompson
    P_qc = 1; P_qr = 2; P_qi = 3; P_qs = 4; P_qg = 5
              P_nr = 6; P_ni = 7
    nscalar = 7

    graupel_ON = 1
    hail_ON = 0
  END IF

  qgh_opt = get_qgh_opt(graupel_ON, hail_ON)

  SELECT CASE (mp_option)
  CASE(2,3,4)             ! LIN
      iret = 0
  CASE(5,6,7)             ! WSM6
      iret = 0
  CASE(14)                ! NSSL 2-mom
      MFflg = 2
      iret = 0
  CASE(108)               ! Thompson
      MFflg = 3
      iret = 0
  CASE DEFAULT            ! not ready for dbz operator
      iret = -1
  END SELECT

  init_mphyopt = iret

end function init_mphyopt

!+---+-----------------------------------------------------------------+

INTEGER(kind_int) FUNCTION get_qgh_opt(graupel_ON, hail_ON)

  INTEGER(kind_int), INTENT(IN) :: graupel_ON,hail_ON

  IF(graupel_ON == 0 .and. hail_ON == 0) THEN
    get_qgh_opt = 1
  ELSE IF(graupel_ON == 0 .and. hail_ON == 1) THEN
    get_qgh_opt = 2
  ELSE IF(graupel_ON == 1 .and. hail_ON == 0) THEN
    get_qgh_opt = 3
  ELSE IF(graupel_ON == 1 .and. hail_ON == 1) THEN
    get_qgh_opt = 4
  ENDIF

END FUNCTION get_qgh_opt

!+---+-----------------------------------------------------------------+

END MODULE
