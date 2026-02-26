/*
 * (C) Copyright 2018-2019 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_UTILS_CONSTANTS_H_
#define UFO_UTILS_CONSTANTS_H_

#include <cmath>

//-------------------------------------------------------------------------------------------------

namespace ufo {

//-------------------------------------------------------------------------------------------------

/// Some useful constants
struct Constants {
  static constexpr double deg2rad        = M_PI / 180.;
  static constexpr double rad2deg        = 180. * M_1_PI;
  static constexpr double grav           = 9.80665e+0;    // Gravity m / s^2
  static constexpr double speed_of_light = 299792458;     // m / s
  static constexpr double t0c            = 2.7315e+2;     // temperature at zero celsius (K)
  static constexpr double ttp            = 2.7316e+2;     // temperature at h2o triple point (K)
  static constexpr double rd             = 2.8705e2;      // specific gas constant for dry air
                                                          // (J K^-1 kg^-1)
  static constexpr double rv             = 4.6150e2;      // specific gas constant for water vapor
                                                          // (J K^-1 kg^-1)
  static constexpr double cp             = 1.0046e3;      // heat capacity at constant pressure
                                                          // for dry air (J K^-1 kg^-1)
  static constexpr double cv             = 7.1760e2;      // heat capacity at constant volume
                                                          // for dry air (J K^-1 kg^-1)
  static constexpr double k_B            = 1.380649e-23;  // Boltzmann constant (J K-1)
  static constexpr double pref           = 1.0e5;         // Reference pressure for calculating
                                                          //      exner
  static constexpr double rd_over_rv     = 0.62198;
          // ratio of molecular weight of water (18.01528 g/mol) to dry air (28.9645 g/mol)
          // Often called "epsilon" or "gas constant ratio"
  static constexpr double rd_over_cp     = rd/cp;
  static constexpr double cv_over_cp     = cv/cp;
  static constexpr double rv_over_rd     = rv/rd;
  static constexpr double rd_over_g      = rd/grav;
  static constexpr double g_over_rd      = grav / rd;
  static constexpr double zero           = 0.0;
  static constexpr double quarter        = 0.25;
  static constexpr double half           = 0.5;
  static constexpr double one            = 1.0;
  static constexpr double two            = 2.0;
  static constexpr double four           = 4.0;
  static constexpr double five           = 5.0;
  static constexpr double ten            = 10.0;
  static constexpr double k_t            = 0.65;         // Thermal conductivity of water at 20degC
                                                         // (W m^-1 K^-1)
  static constexpr double L_e            = 2.26e+06;     // Latent heat of vaporization at 373.15K
                                                         // (J kg^-1)
  static constexpr double L_c            = 2.5e+06;      // Latent heat of condensation at 273.15K
                                                         // (J kg^-1)
  static constexpr double eps            = 0.1;          // Albedo of sea water
  static constexpr double Stefan_Boltzmann_const = 5.67e-8;
                                                         // Stefan-Boltzmann constant (W m^-2 K^-4)
  static constexpr double alpha          = 2.7e-4;       // Water thermal expansion coefficient
                                                         // (K^1)
  static constexpr double cw             = 0.015;        // Water specific heat (cal g^-1 degC^-1)
    // Water specific heat = 4184.0  (J kg^-1 K^-1) need to fix this
  static constexpr double v_w            = 0.8e-6;       // Water kinematic viscosity (m^2/s)
  static constexpr double Rou            = 1000.0;
  static constexpr double DU             = 21.4e-6;      // Dobson unit, kg O3/m**2
  static constexpr double Lclr           = 0.0065;       // constant dry adiabatic lapse rate (K/m)
  static constexpr double t2tv           = 0.608;        // constant moist adiabatic lapse rate
                                                         // (degC/km)
  static constexpr double es_w_0         = 611.2;        // saturation vapor pressure of water at
                                                         // 0degC (Pa)
  static constexpr double euzc_0         = 34.0;         // constant for estimating euphotic layer
  static constexpr double euzc_1         = -0.39;        // constant for estimating euphotic layer
  static constexpr double trop_lapse_rate   = 0.002;     // Lapse rate of tropopause (K/m)
                                                         // according to WMO definition

  // International Civil Aviation Organization (ICAO) atmosphere.
  // https://en.wikipedia.org/wiki/International_Standard_Atmosphere#ICAO_Standard_Atmosphere
  // with gpm = height in geopotential metres
  static constexpr double icao_lapse_rate_l = 6.5E-03;   // Lapse rate for levels up
                                                         //      to 11,000 gpm
                                                         //      (lower boundary of the
                                                         //      isothermal layer) [K/gpm]
  static constexpr double icao_lapse_rate_u = -1.0E-03;  // Lapse rate for levels above
                                                         //      20,000 gpm
                                                         //      (upper boundary of the
                                                         //      isothermal layer) [K/gpm]
  static constexpr double icao_height_l     = 11000.0;   // Height of bottom of isothermal
                                                         //      layer [gpm]
  static constexpr double icao_height_u     = 20000.0;   // Height of top of isothermal
                                                         //      layer [gpm]
  static constexpr double icao_temp_surface = 288.15;    // Surface temperature [K]
  static constexpr double icao_temp_isothermal_layer = 216.65;  // Temperature of isothermal
                                                                //      layer [K]
  static constexpr double icao_pressure_surface = 1013.25;  // Assumed surface pressure [hPa]
  static constexpr double icao_pressure_l   = 226.32;    // Assumed pressure at 11,000 gpm [hPa]
  static constexpr double icao_pressure_u   = 54.7487;   // Assumed pressure at 20,000 gpm [hPa]

  // Constants used for conversion between geopotential and geometric heights
  // using  MJ Mahoney's (2001)
  // Originally from src/ufo/operators/gnssro/utils/gnssro_mod_transform.F90
  static constexpr double semi_major_axis   = 6378.1370e3;     // [m]
  static constexpr double semi_minor_axis   = 6356.7523142e3;  // [m]
  static constexpr double mean_earth_rad_m  = ( (2 * semi_major_axis) + semi_minor_axis) / 3;
                                              // Mean radius of the Earth [m]
  static constexpr double grav_polar        = 9.8321849378;    // [m/s2]
  static constexpr double grav_equator      = 9.7803253359;    // [m/s2]
  static constexpr double earth_omega       = 7.292115e-5;     // [rad/s]
  static constexpr double grav_constant     = 3.986004418e14;  // [m3/s2]
  static constexpr double flattening   = (semi_major_axis - semi_minor_axis) / semi_major_axis;
  static constexpr double somigliana   = (semi_minor_axis / semi_major_axis) *
                                         (grav_polar / grav_equator) - 1.0;
  static constexpr double grav_ratio   = earth_omega * earth_omega *
                                         semi_major_axis * semi_major_axis * semi_minor_axis /
                                         grav_constant;
  static constexpr double eccentricity_sq = (semi_major_axis * semi_major_axis -
                                            semi_minor_axis * semi_minor_axis) /
                                            semi_major_axis / semi_major_axis;

  // Constants used for gnssro super refraction check
  static constexpr float superRefractionCritVal = 0.157;  // unit: N m^-1

  // Refractivity coefficients
  static constexpr double dry_refractivity_constant = 0.776;   // Dry refractivity constant a
  static constexpr double wet_refractivity_constant = 3.73e3;  // Wet refractivity constant b
};

//--------------------------------------------------------------------------------------------------
}  // namespace ufo

#endif  // UFO_UTILS_CONSTANTS_H_
