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

#include "oops/util/DateTime.h"
#include "oops/util/Duration.h"
#include "oops/util/missingValues.h"
#include "ufo/utils/Constants.h"

namespace ufo {

namespace formulas {

/*! Various Methods and Formulations available */
enum MethodFormulation {
  // Methods: Met Center
  UKMOmixingratio,  /*!< UKMO mixing ratio specfic formulation */
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
*      - UKMOmixingratio:  as Landolt-Bornstein 1987
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
*      - UKMOmixingratio: as default
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
*      - UKMOmixingratio: as default
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
* \brief Converts pressure to height.
*
* \b Formulation \b available:
*      - NCAR: uses a fast approximation for pressures > 120 hPa and the ICAO atmosphere otherwise.
*      - NOAA: as default
*      - UKMO: as default
*      - DEFAULT: uses the ICAO atmosphere for all pressures.
*
* \param pressure
*     observation pressure in Pa
* \return height in geopotential metres
*/
float Pressure_To_Height(float pressure,
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

// -------------------------------------------------------------------------------------------
/*!
* \brief Calculate the brightness temperature for an input radiance.  To minimize
*        differences this is done in double precision.
*
* \param radiance
*     The input satellite radiance in (W / (m^2 sr m^-1)).
* \param wavenumber
*     The input wavenumber in m^-1
* \param planck1
*     2*h*c*c - (1.191042972e-16 W / (m^2.sr.m-4)) - this has been made optional
*     to allow for rounding differences when porting.
* \param planck2
*     (h*c / T_b) - (1.4387769e-2 m.K) - this has been made optional
*     to allow for rounding differences when porting.
* \return
*     Brightness temperature in K.
*/
double inversePlanck(const double radiance, const double wavenumber,
                     double planck1 = 1.191042972e-16,  // (W / (m^2.sr.m-4))
                     double planck2 = 1.4387769e-2);    // (m.K)

// -------------------------------------------------------------------------------------
/*!
* \brief Get renumbered scan position for satellite instrument which has
* been spatially resampled. By default use the ceiling method of the number
* of fields of view:
*        numpos = std::ceil(scanpos/numFOV)
* where std::ceil calculates the maximum integer from a float calculation.
* Optionally, set floorRemap to true in order to use a variant floor method:
*        numpos = std::floor((scanpos+1)/numFOV)
*
* \param scanpos
*     satellite instrument scan position
* \param numFOV
*     satellite instrument number of fields of fov for an instrument.  For IASI
*     this is 4 as an example.
* \param floorRemap
*     (boolean) use floor method instead of default ceiling method.
* \return newpos
*/
int RenumberScanPosition(int scanpos, int numFOV, bool floorRemap);

// -------------------------------------------------------------------------------------
/*!
* \brief Compute horizontal drift latitude, longitude and time for an atmospheric profile.
* This formula accepts input and output vectors that correspond to the entire data sample
* as well as a vector of the locations of the current profile in the entire sample.
* The latter vector is used to select the relevant observations from the entire sample
* prior to running the horizontal drift algorithm.
* The outputs of the algorithm are placed in the relevant locations in the entire sample.
*
* \b Formulations \b available:
*      - DEFAULT:
*        Calculation is using Eq. 4. of Laroche and Sarrazin (2013).
*        Reference: "Laroche, S. and Sarrazin, R.,
*                   Impact of Radiosonde Balloon Drift on
*                   Numerical Weather Prediction and Verification,
*                   Weather and Forecasting, 28(3), 772-782, 2013."
*      - UKMO: as default
*      - NCAR: as default
*      - NOAA: as default
*
* \param locs
*     Vector of locations of the current profile in the entire data sample.
* \param apply
*     Vector specifying whether a location should be used or not (governed by the where clause).
* \param lat_in
*     Vector of input latitudes in the entire sample [degrees].
* \param lon_in
*     Vector of input longitudes in the entire sample [degrees].
* \param time_in
*     Vector of input datetimes in the entire sample [ISO 8601 format].
* \param height
*     Vector of input heights in the entire sample [m].
* \param windspd
*     Vector of input wind speeds in the entire sample [m/s].
* \param winddir
*     Vector of input wind directions in the entire sample [degrees].
*     Wind direction is defined such that a northerly wind is 0°, an easterly wind is 90°,
*     a southerly wind is 180°, and a westerly wind is 270°.
* \param [out] lat_out
*     Vector of output latitudes in the entire sample [degrees].
* \param [out] lon_out
*     Vector of output longitudes in the entire sample [degress].
* \param [out] time_out
*     Vector of output datetimes in the entire sample [ISO 8601 format].
* \param [in, optional] formulation
*     Method used to determine the horizontal drift positions.
* \param [in, optional] window_end
*     DateTime at the end of the observation window. If set, computed DateTimes that
*     are larger than this value are set to this value.
*/
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
 MethodFormulation formulation = formulas::MethodFormulation::DEFAULT,
 const util::DateTime * const window_end = nullptr);

// -------------------------------------------------------------------------------------
/*!
* \brief Get background pressure at specified height (station, standard, pmsl).
*
* \param PSurfParamA
*     surf_param_a GeoVaL.
* \param PSurfParamB
*     surf_param_b GeoVaL
* \param height
*     Height of the surface observation for which equivalent background pressure
*     is required.
* \return BkP
*/
float BackgroundPressure(float PSurfParamA, float  PSurfParamB, float height);

// -------------------------------------------------------------------------------------
/*!
* \brief Conversion from geometric heights to geopotential heights using MJ Mahoney's (2001).
*
* \parm latitude
*     Vector of input latitudes
* \parm geomH
*     Vector of input geometric height (m)
* \return geopotential height (m)
*/
float Geometric_to_Geopotential_Height(float latitude, float geomH);

// -------------------------------------------------------------------------------------
/*!
* \brief Conversion from geopotential heights to geometric heights using MJ Mahoney's (2001).
*
* \parm latitude
*     Vector of input latitudes
* \parm geopH
*     Vector of input geopotential height (m)
* \return geometric height (m)
*/
float Geopotential_to_Geometric_Height(float latitude, float geopH);
}  // namespace formulas
}  // namespace ufo

#endif  // UFO_VARIABLETRANSFORMS_FORMULAS_H_
