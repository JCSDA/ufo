/*
 * (C) Copyright 2018-2019 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_UTILS_CONSTANTS_H_
#define UFO_UTILS_CONSTANTS_H_

#include <cmath>

//------------------------------------------------------------------------------------------------------

namespace ufo {

//------------------------------------------------------------------------------------------------------

/// Some useful constants
struct Constants {
    static constexpr double deg2rad        = M_PI / 180.;
    static constexpr double rad2deg        = 180. * M_1_PI;
    static constexpr double grav           = 9.80665e+0;
    static constexpr double t0c            = 2.7315e+2;    // temperature at zero celsius (K)
    static constexpr double rd             = 2.8705e2;
    static constexpr double rv             = 4.6150e2;
    static constexpr double cp             = 1.0046e3;     // heat capacity at constant pressure
                                                           //      for air
    static constexpr double cv             = 7.1760e2;     // heat capacity at constant volume
                                                           //      for air
    static constexpr double rd_over_rv     = rd/rv;
    static constexpr double rd_over_cp     = rd/cp;
    static constexpr double cv_over_cp     = cv/cp;
    static constexpr double rv_over_rd     = rv/rd;
    static constexpr double rd_over_g      = rd/grav;
    static constexpr double mean_earth_rad = 6371.0;
    static constexpr double zero           = 0.0;
    static constexpr double quarter        = 0.25;
    static constexpr double half           = 0.5;
    static constexpr double one            = 1.0;
    static constexpr double two            = 2.0;
    static constexpr double four           = 4.0;
    static constexpr double five           = 5.0;
    static constexpr double ten            = 10.0;
    static constexpr double k_t            = 0.65;         // Thermal conductivity of water
    static constexpr double L_e            = 2.26e+06;     // Latent heat of vaporization
    static constexpr double eps            = 0.1;          // Albedo of sea water
    static constexpr double sig            = 5.67e-6;      // Stefan-Boltzmann constant
    static constexpr double alpha          = 2.7e-4;       // Water thermal expansion coefficient
    static constexpr double cw             = 0.015;        // Water specific heat
    static constexpr double v_w            = 0.8e-6;       // Water kinematic viscosity
    static constexpr double S_B            = 0.026;
    static constexpr double gr             = 9.81;
    static constexpr double Rou            = 1000.0;
    static constexpr double DU             = 21.4e-6;      // Dobson unit, kg O3/m**2
    static constexpr double Lclr           = 0.0065;       // constant lapse rate
    static constexpr double t2tv           = 0.608;        // constant lapse rate
    static constexpr double von_karman     = 0.41;         // Von Karman Constant
    static constexpr double es_w_0         = 611.2;        // saturation vapor pressure of water at
                                                       //    0degC
};
//------------------------------------------------------------------------------------------------------
}  // namespace ufo

#endif  // UFO_UTILS_CONSTANTS_H_
