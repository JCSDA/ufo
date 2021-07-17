/*
 * (C) Crown copyright 2020, Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_VARIABLETRANSFORMS_FORMULAS_H_
#define UFO_VARIABLETRANSFORMS_FORMULAS_H_

#include <algorithm>
#include <cmath>
#include <string>
#include <vector>

#include "oops/util/missingValues.h"
#include "ufo/utils/Constants.h"

namespace ufo {

namespace formulas {

/*! Various Methods and Formulations available */
enum MethodFormulation {
  // Methods: Met Center
  UKMO,   /*!< UKMO specific formulation */
  NCAR,   /*!< NCAR specific formulation */
  NOAA,   /*!< NOAA specific formulation */
  DEFAULT, /*!< DEFAULT formulation */

  // Formulations: Specific authors
  Murphy,
  Sonntag,
  LandoltBornstein,
  Walko,
  Rogers
};

MethodFormulation resolveMethods(const std::string& input);

MethodFormulation resolveFormulations(const std::string& input, const std::string& method);

// -------------------------------------------------------------------------------------
/*!
* \brief Calculates saturated vapour pressure from temperature
*
* \b Formulation \b available:
*      - UKMO:  as Sonntag 1997
*      - Sonntag:
*        Calculation is using the Eqn 7 Sonntag (1997)
*        Reference: "Sonntag, D., Advancements in the field of hygrometry,
*        Meteorol. Zeitschrift, N. F., 3, 51-66, 1994." .
*      - Walko:
*        Polynomial fit of Goff-Gratch (1946) formulation. (Walko, 1991)
*      - Murphy:
*        Alternative method to Walko 1991 (costs more CPU, more accurate)
*        Reference: "Murphy and Koop, Review of the vapour pressure of ice and
         supercooled water for atmospheric applications, Q. J. R.
         Meteorol. Soc (2005), 131, pp. 1539-1565."
*      - NCAR: as default
*      - NOAA: as default
*      - DEFAULT:
*        Classical formula from Rogers and Yau (1989; Eq2.17)
*
*
* \param temp_K
*     Temperature [k]
* \return saturated vapour pressure
*/
float SatVaporPres_fromTemp(const float temp_K,
                   const MethodFormulation formulation = formulas::MethodFormulation::DEFAULT);


// -------------------------------------------------------------------------------------
/*!
* \brief Calculates saturated vapour pressure from temperature
*
* \b Formulation \b available:
*      - NCAR: as Sonntag 1997
*      - NOAA: as Sonntag 1997
*      - UKMO:  as Sonntag 1997
*      - Sonntag:
*        Correct the the saturation vapour pressure of pure water vapour
*        using an enhancement factor needed for moist air.
*        Reference: eg eqns 20, 22 of Sonntag, but for consistency with QSAT the formula
         below is from eqn A4.6 of Adrian Gill's book.
*
*
* \param e_sub_s
*     saturation vapour pressure
* \param temp_K
*     temperature [k]
* \param pressure
*     air pressure [Pa]
* \return saturated vapour pressure
*/
float SatVaporPres_correction(float e_sub_s, float temp_K, float pressure,
                        const MethodFormulation formulation = formulas::MethodFormulation::DEFAULT);
// -------------------------------------------------------------------------------------
/*!
* \brief Calculates Saturated specific humidity or saturated vapour pressure using
* saturation vapour pressure.
*
*
* \b Formulation \b available:
*      - NCAR: as default
*      - NOAA: as default
*      - UKMO: as default
*      - DEFAULT:
*        Calculation is using the Sonntag (1994) formula. (With fix at low
*pressure)
*
* \param Psat
*      saturation vapour pressure of pure water vapour
* \param P
*      Pressure
* \return Saturated specific humidity or saturated vapour pressure
*/
float Qsat_From_Psat(float Psat, float P,
                     MethodFormulation formulation = formulas::MethodFormulation::DEFAULT);

// -------------------------------------------------------------------------------------
/*!
* \brief Derive Virtual Temperature from saturation vapour pressure, pressure and temperature
*
* \b Formulation \b available:
*      - NCAR: as default
*      - NOAA: as default
*      - UKMO: as default
*      - DEFAULT: \f$ Tv = T * ((P + Psat / \epsilon) / (P + Psat)) \f$
*
* \param Psat
*      saturation vapour pressure of pure water vapour
* \param T
*      Temperature
* \param P
*      Pressure
* \return Virtual Temperature
*/
float VirtualTemp_From_Psat_P_T(float Psat, float P, float T,
                          MethodFormulation formulation = formulas::MethodFormulation::DEFAULT);

// -------------------------------------------------------------------------------------
/*!
* \brief Derive Virtual Tempreture using Relative humidity, sat. vapour pressure, pressure
* and temperature
*
* \b Formulation \b available:
*      - NCAR: as default
*      - NOAA: as default
*      - UKMO: as default
*      - DEFAULT: \f$ \alpha = Psat * Rh * 0.01 \f$
*
* \param Rh
*     Relative humidity
* \param Psat
*     saturation vapour pressure of pure water vapour
* \param T
*     Temperature
* \param P
*     Pressure
* \return Virtual Temperature
*/
float VirtualTemp_From_Rh_Psat_P_T(float Rh, float Psat, float P, float T,
                            MethodFormulation formulation = formulas::MethodFormulation::DEFAULT);

// -------------------------------------------------------------------------------------
/*!
* \brief Converts height to pressure using the International Civil Aviation Organization
* (ICAO) atmosphere.
*
* \b Formulation \b available:
*      - NCAR: as default
*      - NOAA: as default
*      - UKMO: as default
*      - DEFAULT: using ICAO standard
*
* \param height
*     observation height in geopotential metres [gpm]
* \return pressure
*/
float Height_To_Pressure_ICAO_atmos(float Height,
                            MethodFormulation formulation = formulas::MethodFormulation::DEFAULT);

// -------------------------------------------------------------------------------------
/*!
* \brief Converts u and v wind component into wind direction.
* Wind direction is defined such that a northerly wind is 0°, an easterly wind is 90°,
* a southerly wind is 180°, and a westerly wind is 270°.
*
* \param u
*     eastward (u) wind component[m/s]
* \param v
*     northward (v) wind component[m/s]
* \return windDirection
*/
float GetWindDirection(float u, float v);

// -------------------------------------------------------------------------------------
/*!
* \brief Converts u and v wind component into wind speed.
*
* \param u
*     eastward (u) wind component[m/s]
* \param v
*     northward (v) wind component[m/s]
* \return windSpeed
*/
float GetWindSpeed(float u, float v);

// -------------------------------------------------------------------------------------
/*!
* \brief Get eastward (u) wind component from wind speed and direction.
*
* \param windSpeed
*     wind speed [m/s]
* \param windFromDirection
*     wind direction [degree]
* \return u
*/
float GetWind_U(float windSpeed, float windFromDirection);

// -------------------------------------------------------------------------------------
/*!
* \brief Get northward (v) wind component from wind speed and direction.
*
* \param windSpeed
*     wind speed [m/s]
* \param windFromDirection
*     wind direction [degree]
* \return v
*/
float GetWind_V(float windSpeed, float windFromDirection);

// -------------------------------------------------------------------------------------
/*!
* \brief Get renumbered scan position 1,2,3,... for satellite instrument
* which has been spatially resampled and for which scan position is 2,5,8,...
*
* \param scanpos
*     satellite instrument scan position
* \return newpos
*/
int RenumberScanPosition(int scanpos);

}  // namespace formulas
}  // namespace ufo

#endif  // UFO_VARIABLETRANSFORMS_FORMULAS_H_
