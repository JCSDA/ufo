/*
 * (C) Crown copyright 2020, Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_CALCULATE_FORMULAS_H_
#define UFO_CALCULATE_FORMULAS_H_

#include <algorithm>
#include <cmath>
#include <string>
#include <vector>

#include "oops/util/missingValues.h"
#include "ufo/utils/Constants.h"

namespace ufo {

namespace formulas {

/*! Various formulas available */
enum Method {
  UKMO,   /*!< UKMO specific formulation */
  NCAR,   /*!< NCAR specific formulation */
  NOAA,   /*!< NOAA specific formulation */
  DEFAULT /*!< DEFAULT specific formulation */
};

Method resolveMethods(const std::string& input);

// -------------------------------------------------------------------------------------
/*!
* \brief Calculates saturated vapour pressure from temperature
*
* \b Formulation \b available:
*      - NCAR: as default
*      - NOAA: as default
*      - UKMO: as default
*      - DEFAULT:
*        Calculation is using the Eqn 7 Sonntag (1997)
*        Reference: "Sonntag, D., Advancements in the field of hygrometry,
*        Meteorol. Zeitschrift, N. F., 3, 51-66, 1994." .
*
*
* \param T
*     Temperature
* \return saturated vapour pressure
*/
float Psat_From_T(const float T,
                   const Method method = formulas::Method::DEFAULT);

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
                     Method method = formulas::Method::DEFAULT);

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
                                Method method = formulas::Method::DEFAULT);

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
                                   Method method = formulas::Method::DEFAULT);
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
                                    Method method = formulas::Method::DEFAULT);

}  // namespace formulas
}  // namespace ufo

#endif  // UFO_CALCULATE_FORMULAS_H_
