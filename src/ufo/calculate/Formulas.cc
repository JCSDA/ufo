/*
 * (C) Crown copyright 2020, Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <algorithm>
#include <cmath>
#include <string>
#include <vector>
#include "ufo/calculate/Formulas.h"
#include "ufo/utils/Constants.h"

namespace ufo {

namespace formulas {

Method resolveMethods(const std::string& input) {
  if (input == "UKMO") return UKMO;
  if (input == "NCAR") return NCAR;
  if (input == "NOAA") return NOAA;

  return DEFAULT;
}

/* -------------------------------------------------------------------------------------*/
float Psat_From_T(float T, Method method) {
  float Psat = util::missingValue(1.0f);  // Saturation vapour pressure (Pa)

  switch (method) {
    case formulas::Method::NCAR:
    case formulas::Method::NOAA:
    case formulas::Method::UKMO:
    default: {
      /*! Note-1
       * Source: Eqn 7, Sonntag, D., Advancements in the field of hygrometry,
       *                Meteorol. Zeitschrift, N. F., 3, 51-66, 1994.
       * Most radiosonde manufacturers use Wexler, or Hyland and Wexler
       * or Sonntag formulations, which are all very similar (Holger Vomel,
       * pers. comm., 2011)
      */
      if (T > 0.0f) {
        Psat = std::exp(-6096.9385f / T + 21.2409642f - 2.711193E-2f * T +
                        1.673952E-5f * T * T + 2.433502f * std::log(T));
      } else {
        Psat = 0.0;
      }

      float FsubW;  // Enhancement factor - see Note-2
      float P = -1.0;
      /*! Note-2
       * Psat above is the saturation vapour pressure of pure water vapour
       * FsubW (~ 1.005 at 1000 hPa) is the enhancement factor needed for moist
       * air
       * (eg eqns 20, 22 of Sonntag, but for consistency with QSAT the formula
       * below is from eqn A4.6 of Adrian Gill's book)
      */
      FsubW =
          1.0f + 1.0E-8f * P * (4.5f + 6.0E-4f * (T - Constants::t0c) *
                                           (T - Constants::t0c));

      Psat = Psat * FsubW;
      break;
    }
  }
  return Psat;
}

/* -------------------------------------------------------------------------------------*/

float Qsat_From_Psat(float Psat, float P, Method method) {
  float QSat = util::missingValue(1.0f);  // Saturated specific humidity or
                                          // saturated vapour pressure (if P<0)

  switch (method) {
    case formulas::Method::NCAR:
    case formulas::Method::NOAA:
    case formulas::Method::UKMO:
    default: {
      // Calculation using the Sonntag (1994) formula. (With fix at low
      // pressure)
      QSat = (Constants::epsilon * Psat) /
             (std::max(P, Psat) - (1.0f - Constants::epsilon) * Psat);
      break;
    }
  }
  return QSat;
}

/* -------------------------------------------------------------------------------------*/

float VirtualTemp_From_Psat_P_T(float Psat, float P, float T, Method method) {
  float Tv = util::missingValue(1.0f);  // virtual temperature

  switch (method) {
    case formulas::Method::NCAR:
    case formulas::Method::NOAA:
    case formulas::Method::UKMO:
    default: {
      Tv = T * ((P + Psat / Constants::epsilon) / (P + Psat));
      break;
    }
  }
  return Tv;
}

/* -------------------------------------------------------------------------------------*/

float VirtualTemp_From_Rh_Psat_P_T(float Rh, float Psat, float P, float T,
                                   Method method) {
  float Tv = util::missingValue(1.0f);  // virtual temperature

  switch (method) {
    case formulas::Method::NCAR:
    case formulas::Method::NOAA:
    case formulas::Method::UKMO:
    default: {
      float PsatRh  = Psat * Rh * 0.01f;
      Tv = VirtualTemp_From_Psat_P_T(PsatRh , P, T, method);
      break;
    }
  }

  return Tv;
}

/* -------------------------------------------------------------------------------------*/

float Height_To_Pressure_ICAO_atmos(float height, Method method) {
  const float missingValueFloat = util::missingValue(1.0f);
  float Pressure = missingValueFloat;

  switch (method) {
    case formulas::Method::NCAR:
    case formulas::Method::NOAA:
    case formulas::Method::UKMO:
    default: {
      float RepT_Bot, RepT_Top, ZP1, ZP2;

      RepT_Bot = 1.0 / Constants::icao_temp_surface;
      RepT_Top = 1.0 / Constants::icao_temp_isothermal_layer;
      ZP1 = Constants::g_over_rd / Constants::icao_lapse_rate_l;
      ZP2 = Constants::g_over_rd / Constants::icao_lapse_rate_u;

      if (height <= missingValueFloat) {
        Pressure = missingValueFloat;
      } else if (height < -5000.0) {
        // TODO(david simonin): The original code has this test.
        // Not sure why! Are we expecting very negative height value??
        Pressure = missingValueFloat;
      } else if (height < Constants::icao_height_l) {
        // Heights up to 11,000 geopotential heigh in meter [gpm]
        Pressure = Constants::icao_lapse_rate_l * height * RepT_Bot;
        Pressure = std::pow((1.0 - Pressure), ZP1);
        Pressure = 100.0 * Pressure * Constants::icao_pressure_surface;
      } else if (height < Constants::icao_height_u) {
        // Heights between 11,000 and 20,000 geopotential heigh in meter [gpm]
        Pressure = Constants::g_over_rd * (height - Constants::icao_height_l) * RepT_Top;
        Pressure = std::log(Constants::icao_pressure_l) - Pressure;
        Pressure = 100.0 * std::exp(Pressure);
      } else {
        // Heights above 20,000 geopotential heigh in meter [gpm]
        Pressure = Constants::icao_lapse_rate_u * RepT_Top *
                   (height - Constants::icao_height_u);
        Pressure = 100.0 * Constants::icao_pressure_u *
                   std::pow((1.0 - Pressure), ZP2);
      }
      break;
    }
  }
  return Pressure;
}

}  // namespace formulas
}  // namespace ufo
