!==========================================================================
module ufo_constants_mod
!==========================================================================

use kinds

implicit none
real(kind_real), parameter, public :: grav    = 9.80665e+0_kind_real  ! Gravity (m/s^2)
real(kind_real), parameter, public :: t0c     = 2.7315e+2_kind_real  ! temperature at zero celsius     (K)
real(kind_real), parameter, public :: rd     = 2.8705e2_kind_real  ! specific gas constant for dry air (J K^-1 kg^-1)
real(kind_real), parameter, public :: rv     = 4.6150e2_kind_real  ! specific gas constant for water vapor (J K^-1 kg^-1)
real(kind_real), parameter, public :: cp     = 1.0046e3_kind_real  ! heat capacity at constant pressure for dry air (J K^-1 kg^-1)
real(kind_real), parameter, public :: cv     = 7.1760e2_kind_real  ! heat capacity at constant volume for dry air (J K^-1 kg^-1)
real(kind_real), parameter, public :: avogadro   = 6.022e23_kind_real
real(kind_real), parameter, public :: gas_constant = 8.314_kind_real  ! R - universal gas constant
real(kind_real), parameter, public :: rd_over_rv = 0.62198_kind_real  ! ratio of molecular weights of water (18.01528 g/mol) to dry air (28.9645 g/mol)  (Unitless).
                                  ! Often called "epsilon" or the "gas constant ratio"
real(kind_real), parameter, public :: rd_over_cp = rd/cp
real(kind_real), parameter, public :: cv_over_cp = cv/cp
real(kind_real), parameter, public :: rv_over_rd = rv/rd
real(kind_real), parameter, public :: rd_over_g  = rd/grav
real(kind_real), parameter, public :: zero    = 0.0_kind_real
real(kind_real), parameter, public :: quarter = 0.25_kind_real
real(kind_real), parameter, public :: half    = 0.5_kind_real
real(kind_real), parameter, public :: one     = 1.0_kind_real
real(kind_real), parameter, public :: two     = 2.0_kind_real
real(kind_real), parameter, public :: three   = 3.0_kind_real
real(kind_real), parameter, public :: four    = 4.0_kind_real
real(kind_real), parameter, public :: five    = 5.0_kind_real
real(kind_real), parameter, public :: six     = 6.0_kind_real
real(kind_real), parameter, public :: ten     = 10.0_kind_real
real(kind_real), parameter, public :: k_t   = 0.65_kind_real       !> Thermal conductivity of water (W m^-1 K^-1)
real(kind_real), parameter, public :: L_e   = 2.26e+06_kind_real !> Latent heat of vaporization at 373.15K (J kg^-1)
real(kind_real), parameter, public :: eps   = 0.1_kind_real      !> Albedo of sea water
real(kind_real), parameter, public :: Stefan_Boltzmann_const = 5.67e-8_kind_real  !> Stefan-Boltzmann constant (W m^-2 K^-4)
real(kind_real), parameter, public :: alpha = 2.7e-4_kind_real !> Water thermal expansion coefficient (K^1)
real(kind_real), parameter, public :: cw    = 0.015_kind_real     !> Water specific heat (cal g^-1 degC^-1)
! 4184 J⋅kg−1⋅K−1 is given as specific heat of liquid water, need to fix this
real(kind_real), parameter, public :: v_w   = 0.8e-6_kind_real    !> Water kinematic viscosity (m^2/s)
real(kind_real), parameter, public :: c_virtual = 1./rd_over_rv-1. ! Related to gas-constant
real(kind_real), parameter, public :: S_B   = 0.026_kind_real   ! S is mean salinity, B is salinity expansion coefficient
                                                      ! SB is relatively constant value of 0.026 in ocean
real(kind_real), parameter, public :: Rou    = 1000.0_kind_real
real(kind_real), parameter, public :: DU    = 21.4e-6_kind_real !Dobson unit, kg O3/m**2
real(kind_real), parameter, public :: Lclr   = 0.0065_kind_real ! constant dry adiabatic lapse rate (K/m)
real(kind_real), parameter, public :: t2tv   = 0.608_kind_real ! constant moist adiabatic lapse rate (degC/km)
real(kind_real), parameter, public :: es_w_0 = 611.2_kind_real ! saturation vapor pressure of water at 0C
real(kind_real), parameter, public :: pi      = acos(-one)
real(kind_real), parameter, public :: deg2rad =  pi/180.0_kind_real
real(kind_real), parameter, public :: rad2deg = one/deg2rad
real(kind_real), parameter, public :: pref = 1.0e5_kind_real
real(kind_real), parameter, public :: hplanck = 6.62607015e-34_kind_real ! Planck's constant in J.s
real(kind_real), parameter, public :: speed_of_light = 2.99792458e8_kind_real ! Speed of light in m/s
real(kind_real), parameter, public :: kboltz = 1.380649e-23_kind_real ! Boltzmann's constant in J/K
! constants relating to WGS-84 ellipsoid and gravity above ellipsoid
real(kind_real), parameter, public :: ecc = 0.081819_kind_real            ! eccentricity
real(kind_real), parameter, public :: k_somig = 1.931853e-3_kind_real     ! Somigliana's constant
real(kind_real), parameter, public :: semi_major_axis = 6378.137e3_kind_real      ! semi-major axis of earth (m)
real(kind_real), parameter, public :: semi_minor_axis = 6356.7523142e3_kind_real  ! semi-minor axis of earth (m)
real(kind_real), parameter, public :: mean_earth_rad_m = ( (2 * semi_major_axis ) + semi_minor_axis)/ 3  ! Mean radius of the Earth (m)
real(kind_real), parameter, public :: grav_polar = 9.8321849378_kind_real     ! [m/s2]
real(kind_real), parameter, public :: grav_equator = 9.7803253359_kind_real   ! [m/s2] equatorial gravity
real(kind_real), parameter, public :: earth_omega = 7.292115e-5_kind_real     ! [rad/s]
real(kind_real), parameter, public :: grav_constant = 3.986004418e14_kind_real  ! [m3/s2]
real(kind_real), parameter ::  flattening  = (semi_major_axis-semi_minor_axis)/semi_major_axis
real(kind_real), parameter ::  somigliana  = (semi_minor_axis/semi_major_axis) * (grav_polar/grav_equator) - one
real(kind_real), parameter ::  grav_ratio  = (earth_omega*earth_omega * &
                                              semi_major_axis*semi_major_axis * semi_minor_axis) / grav_constant
real(kind_real), parameter ::  eccentricity = sqrt(semi_major_axis**2 - semi_minor_axis**2)/semi_major_axis

real(kind_real), parameter, public :: flatt = 0.003352811_kind_real       ! flattening of oblate ellipsoid
real(kind_real), parameter, public :: m_ratio= 0.003449787_kind_real      ! gravity ratio
   ! ratio of centrifugal to gravitational force on the equator

! constants used in CRTM interface
real(kind_real), parameter, public :: kg_to_g = 1000.0_kind_real
real(kind_real), parameter, public :: co2_rescale_to_ppmv = 1.e6_kind_real
real(kind_real), parameter, public :: co2_ppmv_value = 407.0_kind_real
real(kind_real), parameter, public :: midpoint_julday = -1.0_kind_real

! constants used in RTTOV interface
real(kind_real), parameter, public :: g_to_kg = 0.001_kind_real
real(kind_real), parameter, public :: m_to_km = 0.001_kind_real
real(kind_real), parameter, public :: Pa_to_hPa = 0.01_kind_real
real(kind_real), parameter, public :: min_q = 3.0e-6_kind_real
real(kind_real), parameter, public :: min_clw = zero
real(kind_real), parameter, public :: min_ciw = zero
real(kind_real), parameter, public :: RTTOV_ToA = 0.0001_kind_real ! hPa 

! constants used in Ops_QSat and Ops_QSatWat
real(kind_real), parameter, public :: ZeroDegC = 273.15_kind_real

! constants used in avgkernel
real(kind_real), parameter, public :: M_dryair = 0.0289645_kind_real ! molecular weight of dry air (kg.mol-1)
real(kind_real), parameter, public :: M_no2 = 0.0460055_kind_real ! molecular weight of nitrogen dioxide (kg.mol-1)
real(kind_real), parameter, public :: M_co = 0.0280101_kind_real ! molecular weight of carbon monoxide (kg.mol-1)
real(kind_real), parameter, public :: M_o3 = 0.047997_kind_real ! molecular weight of ozone (kg.mol-1)
real(kind_real), parameter, public :: M_hcho = 0.030031_kind_real ! molecular weight of formaldehyde (kg.mol-1)

! constants used in DirectZDA
real(kind_real), parameter, public :: T_melt=273.15_kind_real
real(kind_real), parameter, public :: rhor=1000._kind_real     ! Density of rain (kg m**-3)
real(kind_real), parameter, public :: rhoh=913._kind_real      ! Density of hail (kg m**-3)
real(kind_real), parameter, public :: rhos=100._kind_real      ! Density of snow (kg m**-3)
real(kind_real), parameter, public :: rhog=400._kind_real      ! Density of graupel (kg m**-3)
real(kind_real), parameter, public :: am_s = 0.069_kind_real   ! Empirical constant used in radar Z for dry snow;from Thompson MP

end module ufo_constants_mod

