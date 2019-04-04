MODULE fv3_utils

!parts taken from fv3 fv_sg.F90 and constants.F90

  USE platform_mod, ONLY: r8_kind
  IMPLICIT NONE
  PRIVATE

  REAL,               PUBLIC, PARAMETER :: RADIUS = 6.3712e+6_r8_kind           !< Radius of the Earth [m]
  REAL(kind=r8_kind), PUBLIC, PARAMETER :: PI_8   = 3.1415926535897931_r8_kind  !< Ratio of circle circumference to diameter [N/A]
  REAL,               PUBLIC, PARAMETER :: PI     = 3.1415926535897931_r8_kind  !< Ratio of circle circumference to diameter [N/A] (REAL(KIND=8))
  REAL,               PUBLIC, PARAMETER :: OMEGA  = 7.2921e-5_r8_kind   !< Rotation rate of the Earth [1/s]
  REAL,               PUBLIC, PARAMETER :: GRAV   = 9.80665_r8_kind     !< Acceleration due to gravity [m/s^2]
  REAL(kind=r8_kind), PUBLIC, PARAMETER :: GRAV_8 = 9.80665_r8_kind     !< Acceleration due to gravity [m/s^2] (REAL(KIND=8))
  REAL,               PUBLIC, PARAMETER :: RDGAS  = 287.05_r8_kind      !< Gas constant for dry air [J/kg/deg]
  REAL,               PUBLIC, PARAMETER :: RVGAS  = 461.50_r8_kind      !< Gas constant for water vapor [J/kg/deg]
! Extra:
  REAL,               PUBLIC, PARAMETER :: HLV      = 2.5e6_r8_kind     !< Latent heat of evaporation [J/kg]
  REAL,               PUBLIC, PARAMETER :: HLF      = 3.3358e5_r8_kind  !< Latent heat of fusion [J/kg]
  REAL,               PUBLIC, PARAMETER :: con_cliq = 4.1855e+3_r8_kind !< spec heat H2O liq [J/kg/K]
  REAL,               PUBLIC, PARAMETER :: con_csol = 2.1060e+3_r8_kind !< spec heat H2O ice [J/kg/K]
  REAL,               PUBLIC, PARAMETER :: CP_AIR = 1004.6_r8_kind      !< Specific heat capacity of dry air at constant pressure [J/kg/deg]
  REAL,               PUBLIC, PARAMETER :: KAPPA  = RDGAS/CP_AIR        !< RDGAS / CP_AIR [dimensionless]
  REAL,               PUBLIC, PARAMETER :: TFREEZE = 273.15_r8_kind     !< Freezing temperature of fresh water [K]
#else

  REAL(kind=8), PUBLIC, PARAMETER :: PI_8   = 3.14159265358979323846_r8_kind  !< Ratio of circle circumference to diameter [N/A]
  REAL,         PUBLIC, PARAMETER :: PI     = 3.14159265358979323846_r8_kind  !< Ratio of circle circumference to diameter [N/A]
  REAL,         PUBLIC, PARAMETER :: GRAV   = 9.80_r8_kind             !< Acceleration due to gravity [m/s^2]
  REAL,         PUBLIC, PARAMETER :: RDGAS  = 287.04_r8_kind           !< Gas constant for dry air [J/kg/deg]
  REAL,         PUBLIC, PARAMETER :: RVGAS  = 461.50_r8_kind           !< Gas constant for water vapor [J/kg/deg]
! Extra:
  REAL,         PUBLIC, PARAMETER :: HLV = 2.500e6_r8_kind             !< Latent heat of evaporation [J/kg]
  REAL,         PUBLIC, PARAMETER :: HLF = 3.34e5_r8_kind              !< Latent heat of fusion [J/kg]
  REAL,         PUBLIC, PARAMETER :: KAPPA  = 2.0_r8_kind/7.0_r8_kind  !< RDGAS / CP_AIR [dimensionless]
  REAL,         PUBLIC, PARAMETER :: CP_AIR = RDGAS/KAPPA              !< Specific heat capacity of dry air at constant pressure [J/kg/deg]
  REAL,         PUBLIC, PARAMETER :: TFREEZE = 273.16_r8_kind          !< Freezing temperature of fresh water [K]
#endif

  REAL, PUBLIC, PARAMETER :: STEFAN  = 5.6734e-8_r8_kind !< Stefan-Boltzmann constant [W/m^2/deg^4]

  REAL, PUBLIC, PARAMETER :: CP_VAPOR = 4.0_r8_kind*RVGAS      !< Specific heat capacity of water vapor at constant pressure [J/kg/deg]
  REAL, PUBLIC, PARAMETER :: CP_OCEAN = 3989.24495292815_r8_kind !< Specific heat capacity taken from McDougall (2002) 
!! "Potential Enthalpy ..." [J/kg/deg]
  REAL, PUBLIC, PARAMETER :: RHO0    = 1.035e3_r8_kind  !< Average density of sea water [kg/m^3]
  REAL, PUBLIC, PARAMETER :: RHO0R   = 1.0_r8_kind/RHO0 !< Reciprocal of average density of sea water [m^3/kg]
  REAL, PUBLIC, PARAMETER :: RHO_CP  = RHO0*CP_OCEAN    !< (kg/m^3)*(cal/kg/deg C)(joules/cal) = (joules/m^3/deg C) [J/m^3/deg]

  REAL, PUBLIC, PARAMETER :: ES0 = 1.0_r8_kind        !< Humidity factor. Controls the humidity content of the atmosphere through
!! the Saturation Vapour Pressure expression when using DO_SIMPLE. [dimensionless]
  REAL, PUBLIC, PARAMETER :: DENS_H2O = 1000._r8_kind !< Density of liquid water [kg/m^3]
  REAL, PUBLIC, PARAMETER :: HLS = HLV + HLF          !< Latent heat of sublimation [J/kg]

  REAL, PUBLIC, PARAMETER :: WTMAIR   = 2.896440E+01_r8_kind   !< Molecular weight of air [AMU]
  REAL, PUBLIC, PARAMETER :: WTMH2O   = WTMAIR*(RDGAS/RVGAS)   !< Molecular weight of water [AMU]
  REAL, PUBLIC, PARAMETER :: WTMOZONE =  47.99820_r8_kind      !< Molecular weight of ozone [AMU]
  REAL, PUBLIC, PARAMETER :: WTMC     =  12.00000_r8_kind      !< Molecular weight of carbon [AMU]
  REAL, PUBLIC, PARAMETER :: WTMCO2   =  44.00995_r8_kind      !< Molecular weight of carbon dioxide [AMU]
  REAL, PUBLIC, PARAMETER :: WTMCH4   =  16.0425_r8_kind       !< Molecular weight of methane [AMU]
  REAL, PUBLIC, PARAMETER :: WTMO2    =  31.9988_r8_kind       !< Molecular weight of molecular oxygen [AMU]
  REAL, PUBLIC, PARAMETER :: WTMCFC11 = 137.3681_r8_kind       !< Molecular weight of CFC-11 (CCl3F) [AMU]
  REAL, PUBLIC, PARAMETER :: WTMCFC12 = 120.9135_r8_kind       !< Molecular weight of CFC-21 (CCl2F2) [AMU]
  REAL, PUBLIC, PARAMETER :: WTMN     =  14.0067_r8_kind       !< Molecular weight of Nitrogen [AMU]
  REAL, PUBLIC, PARAMETER :: DIFFAC   = 1.660000E+00_r8_kind   !< Diffusivity factor [dimensionless]
  REAL, PUBLIC, PARAMETER :: AVOGNO   = 6.023000E+23_r8_kind   !< Avogadro's number [atoms/mole]
  REAL, PUBLIC, PARAMETER :: PSTD     = 1.013250E+06_r8_kind   !< Mean sea level pressure [dynes/cm^2]
  REAL, PUBLIC, PARAMETER :: PSTD_MKS = 101325.0_r8_kind       !< Mean sea level pressure [N/m^2]

  REAL, PUBLIC, PARAMETER :: SECONDS_PER_DAY    = 8.640000E+04_r8_kind !< Seconds in a day [s]
  REAL, PUBLIC, PARAMETER :: SECONDS_PER_HOUR   = 3600._r8_kind        !< Seconds in an hour [s]
  REAL, PUBLIC, PARAMETER :: SECONDS_PER_MINUTE = 60._r8_kind          !< Seconds in a minute [s]
  REAL, PUBLIC, PARAMETER :: RAD_TO_DEG         = 180._r8_kind/PI      !< Degrees per radian [deg/rad]
  REAL, PUBLIC, PARAMETER :: DEG_TO_RAD         = PI/180._r8_kind      !< Radians per degree [rad/deg]
  REAL, PUBLIC, PARAMETER :: RADIAN             = RAD_TO_DEG           !< Equal to RAD_TO_DEG for backward compatability. [rad/deg]
  REAL, PUBLIC, PARAMETER :: ALOGMIN            = -50.0_r8_kind        !< Minimum value allowed as argument to log function [N/A]
  REAL, PUBLIC, PARAMETER :: EPSLN              = 1.0e-40_r8_kind      !< A small number to prevent divide by zero exceptions [N/A]

  REAL, PUBLIC, PARAMETER :: RADCON = ((1.0E+02*GRAV)/(1.0E+04*CP_AIR))*SECONDS_PER_DAY !< Factor used to convert flux divergence to
!! heating rate in degrees per day [deg sec/(cm day)]
  REAL, PUBLIC, PARAMETER :: RADCON_MKS  = (GRAV/CP_AIR)*SECONDS_PER_DAY !< Factor used to convert flux divergence to
!! heating rate in degrees per day [deg sec/(m day)]
  REAL, PUBLIC, PARAMETER :: O2MIXRAT    = 2.0953E-01_r8_kind !< Mixing ratio of molecular oxygen in air [dimensionless]
  REAL, PUBLIC, PARAMETER :: RHOAIR      = 1.292269_r8_kind   !< Reference atmospheric density [kg/m^3]
  REAL, PUBLIC, PARAMETER :: VONKARM     = 0.40_r8_kind       !< Von Karman constant [dimensionless]
  REAL, PUBLIC, PARAMETER :: C2DBARS     = 1.e-4_r8_kind      !< Converts rho*g*z (in mks) to dbars: 1dbar = 10^4 (kg/m^3)(m/s^2)m [dbars]
  REAL, PUBLIC, PARAMETER :: KELVIN      = 273.15_r8_kind     !< Degrees Kelvin at zero Celsius [K]

  SUBROUTINE qsmith_init
    INTEGER, PARAMETER:: length=2621 
    INTEGER i

    IF( .NOT. ALLOCATED(table) ) THEN
!                            Generate es table (dT = 0.1 deg. C)

       ALLOCATE ( table(length) )
       ALLOCATE (  des (length) )

       CALL qs_table(length, table)

       DO i=1,length-1
          des(i) = table(i+1) - table(i)
       ENDDO
       des(length) = des(length-1)
    ENDIF

  END SUBROUTINE qsmith_init

  SUBROUTINE qsmith(im, km, k1, t, p, q, qs, dqdt)
! input T in deg K; p (Pa)
    INTEGER, INTENT(in):: im, km, k1
    REAL, INTENT(in),DIMENSION(im,km):: t, p, q
    REAL, INTENT(out),DIMENSION(im,km):: qs
    REAL, INTENT(out), OPTIONAL:: dqdt(im,km)
! Local:
    REAL es(im,km)
    REAL ap1, eps10
    REAL Tmin
    INTEGER i, k, it

    Tmin = tice-160.
    eps10  = 10.*esl

    IF( .NOT. ALLOCATED(table) ) CALL  qsmith_init

    DO k=k1,km
       DO i=1,im
          ap1 = 10.*DIM(t(i,k), Tmin) + 1.
          ap1 = MIN(2621., ap1)
          it = ap1
          es(i,k) = table(it) + (ap1-it)*des(it)
          qs(i,k) = esl*es(i,k)*(1.+zvir*q(i,k))/p(i,k)
       ENDDO
    ENDDO

    IF ( PRESENT(dqdt) ) THEN
       DO k=k1,km
          DO i=1,im
             ap1 = 10.*DIM(t(i,k), Tmin) + 1.
             ap1 = MIN(2621., ap1) - 0.5
             it  = ap1
             dqdt(i,k) = eps10*(des(it)+(ap1-it)*(des(it+1)-des(it)))*(1.+zvir*q(i,k))/p(i,k)
          ENDDO
       ENDDO
    ENDIF

  END SUBROUTINE qsmith

  SUBROUTINE qs_table(n,table)
    INTEGER, INTENT(in):: n
    REAL table (n)
    REAL:: dt=0.1
    REAL esbasw, tbasw, esbasi, tbasi, Tmin, tem, aa, b, c, d, e, esh20 
    REAL wice, wh2o
    INTEGER i
! Constants
    esbasw = 1013246.0
    tbasw =   373.16
    tbasi =   273.16
    Tmin = tbasi - 160.
!  Compute es over water
!  see smithsonian meteorological tables page 350.
    DO  i=1,n
       tem = Tmin+dt*REAL(i-1)
       aa  = -7.90298*(tbasw/tem-1)
       b   =  5.02808*alog10(tbasw/tem)
       c   = -1.3816e-07*(10**((1-tem/tbasw)*11.344)-1)
       d   =  8.1328e-03*(10**((tbasw/tem-1)*(-3.49149))-1)
       e   = alog10(esbasw)
       table(i)  = 0.1*10**(aa+b+c+d+e)
    ENDDO

  END SUBROUTINE qs_table

END  MODULE fv3_utils
