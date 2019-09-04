!==========================================================================
module ufo_constants_mod
!==========================================================================

use kinds
use iso_c_binding

implicit none
real(kind_real), parameter, public :: rad2deg = 57.29577954572      ! copy from fv3jedi_constants.f90
real(kind_real), parameter, public :: deg2rad =  0.01745329251
real(kind_real), parameter, public :: grav    = 9.80665e+0_kind_real
real(kind_real), parameter, public :: t0c     = 2.7315e+2_kind_real ! temperature at zero celsius     (K)
real(kind_real), parameter, public :: rd     = 2.8705e2_kind_real
real(kind_real), parameter, public :: rv     = 4.6150e2_kind_real
real(kind_real), parameter, public :: cp     = 1.0046e3_kind_real ! heat capacity at constant pressure for air
real(kind_real), parameter, public :: cv     = 7.1760e2_kind_real ! heat capacity at constant volume for air
real(kind_real), parameter, public :: rd_over_rv = rd/rv
real(kind_real), parameter, public :: rd_over_cp = rd/cp
real(kind_real), parameter, public :: cv_over_cp = cv/cp
real(kind_real), parameter, public :: rv_over_rd = rv/rd
real(kind_real), parameter, public :: rd_over_g  = rd/grav
real(kind_real), parameter, public :: mean_earth_rad = 6371.0
real(kind_real), parameter, public :: zero    = 0.0_kind_real
real(kind_real), parameter, public :: quarter = 0.25_kind_real
real(kind_real), parameter, public :: half    = 0.5_kind_real
real(kind_real), parameter, public :: one     = 1.0_kind_real
real(kind_real), parameter, public :: two     = 2.0_kind_real
real(kind_real), parameter, public :: four    = 4.0_kind_real
real(kind_real), parameter, public :: five    = 5.0_kind_real
real(kind_real), parameter, public :: ten     = 10.0_kind_real
real(kind_real), parameter, public :: k_t   = 0.65       !> Thermal conductivity of water
real(kind_real), parameter, public :: L_e   = 2.26e+06 !> Latent heat of vaporization
real(kind_real), parameter, public :: eps   = 0.1      !> Albedo of sea water
real(kind_real), parameter, public :: sig   = 5.67e-6  !> Stefan-Boltzmann constant
real(kind_real), parameter, public :: alpha = 2.7e-4 !> Water thermal expansion coefficient
real(kind_real), parameter, public :: cw    = 0.015     !> Water specific heat
real(kind_real), parameter, public :: v_w   = 0.8e-6     !> Water kinematic viscosity  
real(kind_real), parameter, public :: S_B   = 0.026
real(kind_real), parameter, public :: gr    = 9.81
real(kind_real), parameter, public :: Rou    = 1000.0
real(kind_real), parameter, public :: DU    = 21.4e-6 !Dobson unit, kg O3/m**2
real(kind_real), parameter, public :: Lclr   = 0.0065 ! constant lapse rate
real(kind_real), parameter, public :: t2tv   = 0.608 ! constant lapse rate
real(kind_real), parameter, public :: von_karman = 0.41_kind_real ! Von Karman Constant
real(kind_real), parameter, public :: es_w_0 = 611.2_kind_real ! saturation vapor pressure of water at 0C
end module ufo_constants_mod

