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
#include "eckit/exception/Exceptions.h"
#include "oops/util/Logger.h"
#include "ufo/utils/Constants.h"
#include "ufo/variabletransforms/Formulas.h"
#include "ufo/variabletransforms/LookupTable.h"

namespace ufo {

namespace formulas {

MethodFormulation resolveMethods(const std::string& input) {
  // Met Centers
  if (input == "UKMOmixingratio") return UKMOmixingratio;
  if (input == "UKMO") return UKMO;
  if (input == "NCAR") return NCAR;
  if (input == "NOAA") return NOAA;
  return DEFAULT;
}

MethodFormulation resolveFormulations(const std::string& input, const std::string& method) {
  if (input.empty()) {
    //  if we do not have any formulation set, then we assign
    // the method as formulation
    return resolveMethods(method);
  } else {
    // Formulation authors
    if (input == "Sonntag") return Sonntag;
    if (input == "Walko") return Walko;
    if (input == "Murphy") return Murphy;
    return DEFAULT;
  }
}

/* -------------------------------------------------------------------------------------*/
float SatVaporPres_fromTemp(float temp_K, MethodFormulation formulation) {
  const float missingValueFloat = util::missingValue<float>();
  const float t0c = static_cast<float>(ufo::Constants::t0c);
  float e_sub_s = missingValueFloat;  // Saturation vapour pressure (Pa)

  switch (formulation) {
    case formulas::MethodFormulation::UKMO:
    case formulas::MethodFormulation::Sonntag: {
      /* I. Source: Eqn 7, Sonntag, D., Advancements in the field of hygrometry,
       *     Meteorol. Zeitschrift, N. F., 3, 51-66, 1994.
       *     Most radiosonde manufacturers use Wexler, or Hyland and Wexler
       *     or Sonntag formulations, which are all very similar (Holger Vomel,
       *     pers. comm., 2011)
      */
      if (temp_K != missingValueFloat) {
        e_sub_s = std::exp(-6096.9385f / temp_K + 21.2409642f - 2.711193E-2f * temp_K +
                        1.673952E-5f * temp_K * temp_K + 2.433502f * std::log(temp_K));
      } else {
        e_sub_s = 0.0f;
      }
      break;
    }
    case formulas::MethodFormulation::UKMOmixingratio:
    case formulas::MethodFormulation::LandoltBornstein: {
      /* Returns a saturation mixing ratio given a temperature and pressure
         using saturation vapour pressures caluclated using the Goff-Gratch
         formulae, adopted by the WMO as taken from Landolt-Bornstein, 1987
         Numerical Data and Functional relationships in Science and
         Technology.  Group V/Vol 4B Meteorology.  Physical and Chemical
         properties of Air, P35.
      */
      if (temp_K != missingValueFloat) {
        float adj_Temp;
        float lookup_a;
        int lookup_i;

        const float Low_temp_thd = 183.15;   // Lowest temperature for which look-up table is valid
        const float High_temp_thd = 338.15;  // Highest temperature for which look-up table is valid
        const float Delta_Temp = 0.1;        // Temperature increment of look-up table

        //  Use the lookup table to find saturated vapour pressure.
        adj_Temp = std::max(Low_temp_thd, temp_K);
        adj_Temp = std::min(High_temp_thd, adj_Temp);

        lookup_a = (adj_Temp - Low_temp_thd + Delta_Temp) / Delta_Temp;
        lookup_i = static_cast<int>(lookup_a);
        lookup_a = lookup_a - lookup_i;
        e_sub_s = (1.0 - lookup_a) *
                  lookuptable::LandoltBornstein_lookuptable[lookup_i] +
                  lookup_a *
                  lookuptable::LandoltBornstein_lookuptable[lookup_i + 1];
      } else {
        e_sub_s = 0.0f;
      }
      break;
    }
    case formulas::MethodFormulation::Walko: {
          // Polynomial fit of Goff-Gratch (1946) formulation. (Walko, 1991)
          float x = std::max(-80.0f, temp_K-t0c);
          const std::vector<float> c{610.5851f, 44.40316f, 1.430341f, 0.2641412e-1f,
            0.2995057e-3f, 0.2031998e-5f, 0.6936113e-8f, 0.2564861e-11f, -0.3704404e-13f};
          e_sub_s = c[0]+x*(c[1]+x*(c[2]+x*(c[3]+x*(c[4]+x*(c[5]+x*(c[6]+x*(c[7]+x*c[8])))))));
          break;
    }
    case formulas::MethodFormulation::Murphy: {
      // ALTERNATIVE (costs more CPU, more accurate than Walko, 1991)
      // Source: Murphy and Koop, Review of the vapour pressure of ice and
      //       supercooled water for atmospheric applications, Q. J. R.
      //       Meteorol. Soc (2005), 131, pp. 1539-1565.
      e_sub_s = std::exp(54.842763f - 6763.22f / temp_K - 4.210f * std::log(temp_K)
                     + 0.000367f * temp_K + std::tanh(0.0415f * (temp_K - 218.8f))
                     * (53.878f - 1331.22f / temp_K - 9.44523f * std::log(temp_K)
                     + 0.014025f * temp_K));
      break;
    }
    case formulas::MethodFormulation::NCAR:
    case formulas::MethodFormulation::NOAA:
    case formulas::MethodFormulation::Rogers:
    default: {
      // Classical formula from Rogers and Yau (1989; Eq2.17)
      e_sub_s = 1000. * 0.6112 * std::exp(17.67f * (temp_K - t0c) / (temp_K - 29.65f));
      break;
    }
  }
  return e_sub_s;
}

/* -------------------------------------------------------------------------------------*/
float SatVaporPres_correction(float e_sub_s, float temp_K, float pressure,
                              MethodFormulation formulation) {
  const float t0c = static_cast<float>(ufo::Constants::t0c);

  switch (formulation) {
    case formulas::MethodFormulation::NCAR:
    case formulas::MethodFormulation::NOAA:
    case formulas::MethodFormulation::UKMOmixingratio:
    case formulas::MethodFormulation::UKMO:
    case formulas::MethodFormulation::Sonntag: {
      /* e_sub_s is the saturation vapour pressure of pure water vapour.FsubW (~ 1.005
         at 1000 hPa) is the enhancement factor needed for moist air.
         If P is set to -1: then eg eqns 20, 22 of Sonntag is used
         If P > 0 then eqn A4.6 of Adrian Gill's book is used to guarantee consistency with the
         saturated specific humidity.
      */
      float FsubW;  // Enhancement factor
      FsubW = 1.0f + 1.0E-8f * pressure * (4.5f + 6.0E-4f * (temp_K - t0c) *(temp_K - t0c));
      e_sub_s = e_sub_s * FsubW;

      break;
    }
    default: {
      std::string errString = "Aborting, no method matches enum formulas::MethodFormulation";
      oops::Log::error() << errString;
      throw eckit::BadValue(errString);
    }
  }
  return e_sub_s;
}
/* -------------------------------------------------------------------------------------*/

float Qsat_From_Psat(float Psat, float P, MethodFormulation formulation) {
  float QSat = util::missingValue<float>();  // Saturated specific humidity or
                                             // saturated vapour pressure (if P<0)

  switch (formulation) {
    case formulas::MethodFormulation::NCAR:
    case formulas::MethodFormulation::NOAA:
    case formulas::MethodFormulation::UKMOmixingratio:
    case formulas::MethodFormulation::UKMO:
    default: {
      // Calculation using the Sonntag (1994) formula. (With fix at low
      // pressure)
      //  Note that at very low pressures we apply a fix, to prevent a
      //     singularity (Qsat tends to 1.0 kg/kg).
      QSat = (Constants::epsilon * Psat) /
             (std::max(P, Psat) - (1.0f - Constants::epsilon) * Psat);
      break;
    }
  }
  return QSat;
}

/* -------------------------------------------------------------------------------------*/

// VirtualTemperature()
float VirtualTemp_From_Psat_P_T(float Psat, float P, float T, MethodFormulation formulation) {
  float Tv = util::missingValue<float>();  // virtual temperature

  switch (formulation) {
    case formulas::MethodFormulation::NCAR:
    case formulas::MethodFormulation::NOAA:
    case formulas::MethodFormulation::UKMO:
    default: {
      Tv = T * ((P + Psat / Constants::epsilon) / (P + Psat));
      break;
    }
  }
  return Tv;
}

/* -------------------------------------------------------------------------------------*/

float VirtualTemp_From_Rh_Psat_P_T(float Rh, float Psat, float P, float T,
                                   MethodFormulation formulation) {
  float Tv = util::missingValue<float>();  // virtual temperature

  switch (formulation) {
    case formulas::MethodFormulation::NCAR:
    case formulas::MethodFormulation::NOAA:
    case formulas::MethodFormulation::UKMO:
    default: {
      float PsatRh  = Psat * Rh * 0.01f;
      Tv = VirtualTemp_From_Psat_P_T(PsatRh , P, T, formulation);
      break;
    }
  }

  return Tv;
}

/* -------------------------------------------------------------------------------------*/

float Height_To_Pressure_ICAO_atmos(float height, MethodFormulation formulation) {
  const float missingValueFloat = util::missingValue<float>();
  float Pressure = missingValueFloat;

  switch (formulation) {
    case formulas::MethodFormulation::NCAR:
    case formulas::MethodFormulation::NOAA:
    case formulas::MethodFormulation::UKMO:
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

/* -------------------------------------------------------------------------------------*/

float Pressure_To_Height(float pressure, MethodFormulation method) {
  const float missingValueFloat = util::missingValue<float>();
  const float pressure_hPa = pressure * 0.01;
  float height = missingValueFloat;

  switch (method) {
    case formulas::MethodFormulation::NCAR:
      // The NCAR-RAL method: a fast approximation for pressures > 120 hPa.
      // Above 120hPa (~15km) use the ICAO atmosphere.
      if (pressure == missingValueFloat || pressure <= 0.0f) {
        height = missingValueFloat;
      } else if (pressure_hPa <= 120.0f &&
                 pressure_hPa > Constants::icao_pressure_u) {
        pressure = std::log(Constants::icao_pressure_l) - std::log(pressure_hPa);
        height = pressure * Constants::icao_temp_isothermal_layer / Constants::g_over_rd +
          Constants::icao_height_l;
      } else if (pressure_hPa <= Constants::icao_pressure_u) {
        pressure = 1.0 - std::pow(pressure_hPa / Constants::icao_pressure_u,
                                  Constants::icao_lapse_rate_u / Constants::g_over_rd);
        height = pressure * Constants::icao_temp_isothermal_layer / Constants::icao_lapse_rate_u +
          Constants::icao_height_u;
      } else {
        height = 44307.692 * (1.0 - std::pow(pressure / 101325.0, 0.190));
      }
      break;

    case formulas::MethodFormulation::NOAA:
    case formulas::MethodFormulation::UKMO:
    default: {
      if (pressure == missingValueFloat || pressure <= 0.0f) {
        height = missingValueFloat;
      } else if (pressure_hPa > Constants::icao_pressure_l) {
        pressure = 1.0 - std::pow(pressure_hPa / Constants::icao_pressure_surface,
                                  Constants::icao_lapse_rate_l / Constants::g_over_rd);
        height = pressure * Constants::icao_temp_surface / Constants::icao_lapse_rate_l;
      } else if (pressure_hPa <= Constants::icao_pressure_l &&
                 pressure_hPa > Constants::icao_pressure_u) {
        pressure = std::log(Constants::icao_pressure_l) - std::log(pressure_hPa);
        height = pressure * Constants::icao_temp_isothermal_layer / Constants::g_over_rd +
          Constants::icao_height_l;
      } else {
        pressure = 1.0 - std::pow(pressure_hPa / Constants::icao_pressure_u,
                                  Constants::icao_lapse_rate_u / Constants::g_over_rd);
        height = pressure * Constants::icao_temp_isothermal_layer / Constants::icao_lapse_rate_u +
          Constants::icao_height_u;
      }
      break;
    }
  }
  return height;
}

float GetWindDirection(float u, float v) {
  const float missing = util::missingValue<float>();
  float windDirection = missing;  // wind direction

  if (u != missing && v != missing) {
    if (u != 0.0f || v != 0.0f) {
        windDirection = fmod((270.0f - atan2(v, u) * Constants::rad2deg), 360.0f);
    } else {
        windDirection = 0.0f;
    }
  }
  return windDirection;
}

float GetWindSpeed(float u, float v) {
  const float missing = util::missingValue<float>();
  float windSpeed = missing;  // wind speed

  if (u != missing && v != missing) {
        windSpeed = hypot(u, v);
  }
  return windSpeed;
}

float GetWind_U(float windSpeed, float windFromDirection) {
  const float missing = util::missingValue<float>();
  float u = missing;  // wind speed

  if (windFromDirection != missing
      && windSpeed != missing && windSpeed >= 0) {
    u = -windSpeed * sin(windFromDirection * Constants::deg2rad);
  }
  return u;
}


float GetWind_V(float windSpeed, float windFromDirection) {
  const float missing = util::missingValue<float>();
  float v = missing;  // wind speed

  if (windFromDirection != missing
      && windSpeed != missing && windSpeed >= 0) {
    v = -windSpeed * cos(windFromDirection * Constants::deg2rad);
  }
  return v;
}

/* -------------------------------------------------------------------------------------
This formula takes a radiance (W / (m^2.sr.m^-1)) and a wavenumber (m^-1) and outputs
a brightness temperature. where:
h - planck constant (m2kgs-1)
c - speed of light (ms-1)
t_b - boltzman constant (m2kgs-2K-1)
planck1 = (2*h*c*c)    - (1.191042972e-16 W / (m^2.sr.m-4)) - forward declaration has the default
                         arguments to allow for rounding differences in porting.
planck2 = (h*c / T_b)  - (1.4387769e-2 m.K) - forward declaration has the default arguments
                         to allow for rounding differences in porting.
*/
double inversePlanck(const double radiance, const double wavenumber,
                     const double planck1, const double planck2) {
  const double p1 = planck1 * wavenumber * wavenumber * wavenumber;
  const double p2 = planck2 * wavenumber;
  double BT = p2 / std::log(1.0 + p1 / radiance);
  return BT;
}

/* -------------------------------------------------------------------------------------*/

int RenumberScanPosition(const int scanpos, const int numFOV, const bool floorRemap) {
  // std::ceil, std::floor have floats as input and output,
  // therefore static casts required
  int newpos;
  if (floorRemap) {
    const float scanpos_numFOV = static_cast<float>(scanpos + 1)/static_cast<float>(numFOV);
    newpos = static_cast<int>(std::floor(scanpos_numFOV));
  } else {
    const float scanpos_numFOV = static_cast<float>(scanpos)/static_cast<float>(numFOV);
    newpos = static_cast<int>(std::ceil(scanpos_numFOV));
  }
  return newpos;
}

/* -------------------------------------------------------------------------------------*/

void horizontalDrift
(const std::vector<size_t> & locs,
 const std::vector<bool> & apply,
 const std::vector<float> & lat_in,
 const std::vector<float> & lon_in,
 const std::vector<util::DateTime> & time_in,
 const std::vector<float> & height,
 const std::vector<float> & windspd,
 const std::vector<float> & winddir,
 std::vector<float> & lat_out,
 std::vector<float> & lon_out,
 std::vector<util::DateTime> & time_out,
 MethodFormulation formulation,
 const util::DateTime * const window_end) {
  const float missingValueFloat = util::missingValue<float>();

  switch (formulation) {
  case formulas::MethodFormulation::NCAR:
  case formulas::MethodFormulation::NOAA:
  case formulas::MethodFormulation::UKMO:
  default: {
    // Location of the first entry in the profile.
    const size_t loc0 = locs.front();

    // Values of latitude, longitude and datetime at the first entry of the profile.
    const double lat0 = lat_in[loc0];
    const double lon0 = lon_in[loc0];
    const util::DateTime time0 = time_in[loc0];

    // The drift computation is not performed for very high latitude sites.
    if (std::abs(lat0) >= 89.0) return;

    // Fill vector of valid locations.
    std::vector<size_t> locs_valid;
    for (size_t jloc : locs) {
      // If not selected by the where clause.
      if (!apply[jloc]) continue;
      // The location is classed as valid if the wind speed and height are not missing.
      if (windspd[jloc] != missingValueFloat && height[jloc] != missingValueFloat)
        locs_valid.push_back(jloc);
    }

    // If there are zero or one valid locations, exit the routine.
    if (locs_valid.size() < 2) return;

    // Average ascent speed (m/s).
    const double ascent_speed = 5.16;

    // Cumulative values of change in time.
    // This value is converted to a util::Duration object rather than performing the
    // same conversion to each individual change in time.
    // This avoids a loss in precision given util::Duration is accurate to the nearest second.
    double dt_cumul = 0.0;

    for (size_t k = 0; k < locs_valid.size() - 1; ++k) {
      // Locations of the current and next valid observations in the profile.
      const size_t loc_current = locs_valid[k];
      const size_t loc_next = locs_valid[k + 1];

      // Compute changes in latitude, longitude and time between adjacent valid levels.
      // Change in height.
      const double dh = height[loc_next] - height[loc_current];
      // Change in time.
      const double dt = dh / ascent_speed;
      // Average eastward and northward wind between the two levels.
      // 180 degrees is subtracted from the wind direction in order to account for the different
      // conventions used in the observations and in this calculation.
      const double avgu = 0.5 *
        (windspd[loc_current] * std::sin((winddir[loc_current] - 180.0) * Constants::deg2rad) +
         windspd[loc_next] * std::sin((winddir[loc_next] - 180.0) * Constants::deg2rad));
      const double avgv = 0.5 *
        (windspd[loc_current] * std::cos((winddir[loc_current] - 180.0) * Constants::deg2rad) +
         windspd[loc_next] * std::cos((winddir[loc_next] - 180.0) * Constants::deg2rad));
      // Total height of the observation above the centre of the Earth.
      const double totalheight = ufo::Constants::mean_earth_rad * 1000.0 + height[loc_current];
      // Change in latitude.
      const double dlat = ufo::Constants::rad2deg * avgv * dt / totalheight;
      // Change in longitude.
      const double dlon = ufo::Constants::rad2deg * avgu * dt /
        (totalheight * std::cos(lat_out[loc_current] * ufo::Constants::Constants::deg2rad));

      // Fill output values.
      lat_out[loc_next] = lat_out[loc_current] + dlat;
      lon_out[loc_next] = lon_out[loc_current] + dlon;
      // Convert the cumulative change in time to a util::Duration.
      dt_cumul += dt;
      // Calculate the level datetime, keeping it within the assimilation window if required.
      const util::DateTime t_cumul = time0 + util::Duration(static_cast<int64_t>(dt_cumul));
      if (window_end)
        time_out[loc_next] = t_cumul > *window_end ? *window_end : t_cumul;
      else
        time_out[loc_next] = t_cumul;
    }

    // Copy latitude, longitude and time at each valid location to all invalid
    // locations that lie between the current valid location and the next one above it.
    double lat = lat0;
    double lon = lon0;
    util::DateTime time = time0;
    for (size_t jloc : locs) {
      if (std::find(locs_valid.begin(), locs_valid.end(), jloc) != locs_valid.end()) {
        lat = lat_out[jloc];
        lon = lon_out[jloc];
        time = time_out[jloc];
      } else {
        lat_out[jloc] = lat;
        lon_out[jloc] = lon;
        time_out[jloc] = time;
      }
    }
    break;
  }
  }
}

/* -------------------------------------------------------------------------------------*/

float BackgroundPressure(float PSurfParamA, float  PSurfParamB, float height) {
  float BkP = util::missingValue<float>();
  double ToRaise =  (PSurfParamA - height)/PSurfParamB;
  if (ToRaise > 0.0) {
    BkP = pow(ToRaise, (Constants::grav/(Constants::Lclr*Constants::rd)));
  }
  return BkP;
}

/* -------------------------------------------------------------------------------------*/

float Geometric_to_Geopotential_Height(float latitude, float geomH) {
  float geopH = util::missingValue<float>();
  double sino = std::pow(std::sin(Constants::deg2rad * latitude), 2);
  double termg = Constants::grav_equator * ((1.0 + Constants::somigliana * sino) /
                 std::sqrt(1.0 - Constants::eccentricity_sq * sino));
  double termr = Constants::semi_major_axis / (1.0 + Constants::flattening +
                 Constants::grav_ratio - 2.0 * Constants::flattening * sino);
  geopH = (termg / Constants::grav) * ((termr * geomH) / (termr + geomH));  // [m]
  return geopH;
}

/* -------------------------------------------------------------------------------------*/

float Geopotential_to_Geometric_Height(float latitude, float geopH) {
  float geomH = util::missingValue<float>();
  double sino = std::pow(std::sin(Constants::deg2rad * latitude), 2);
  double termg = Constants::grav_equator * ((1.0 + Constants::somigliana * sino) /
                 std::sqrt(1.0 - Constants::eccentricity_sq * sino));
  double termr = Constants::semi_major_axis / (1.0 + Constants::flattening +
                 Constants::grav_ratio - 2.0 * Constants::flattening * sino);
  double termrg = termg / Constants::grav * termr;
  geomH = termr * geopH / (termrg - geopH);  // [m]
  return geomH;
}
}  // namespace formulas
}  // namespace ufo
