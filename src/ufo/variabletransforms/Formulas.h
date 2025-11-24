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

enum class Method {
  UKMOmixingratio,  ///< For humidity only: UKMO method where, for the relative
    ///< humidity transform, the relative humidity is approximated to be the
    ///< mixing ratio is divided by the saturation specific humidity.
  UKMO,   ///< UKMO specific methods
  NCAR,   ///< NCAR specific methods
  NOAA,   ///< NOAA specific methods
  UKMOQsatWater,  ///< For specific humidity only: UKMO method where the
    ///< saturation vapour pressure of pure water vapour over a plane water
    ///< surface calculated using the GoffGratchLandoltBornsteinWater
    ///< formulation, with other calculations the same as UKMOQSatIceWater (see
    ///< the method implementation for details).
  UKMOQsatIceWater,  ///< For specific humidity only: UKMO method where the
    ///< saturation vapour pressure of pure water vapour over a plane water
    ///< surface calculated using the GoffGratchLandoltBornsteinIceWater
    ///< formulation, with other calculations the same as UKMOQSatIceWater (see
    ///< the method implementation for details).
  Sonntag,  ///< For humidity only: use the default methods, calculating
     ///< saturation mixing ratio from temperature using Eqn 7, Sonntag, D.,
     ///< Advancements in the field of hygrometry, Meteorol. Zeitschrift, N. F., 3,
     ///< 51-66, 1994.
  Walko,  ///< For humidity only: use the default methods, calculating saturation
    ///< mixing ratio from temperature using polynomial fit of Goff-Gratch (1946)
    ///< formulation (Walko, 1991)
  Murphy,  ///< For humidity only: use the default methods, calculating
    ///< saturation mixing ratio from temperature using alternative method to Walko
    ///< 1991 (costs more CPU, more accurate)
    ///< Reference: Eq. 10 of "Murphy and Koop, Review of the vapour pressure of ice
    ///< and supercooled water for atmospheric applications, Q. J. R. Meteorol. Soc
    ///< (2005), 131, pp. 1539-1565."
  GoffGratchLandoltBornsteinIceWater,  ///< For humidity only: use the default
    ///< methods, calculating saturation mixing ratio from temperature using the
    ///< Goff-Gratch formulae over ice below and including 273.15 K (0 C) and over
    ///< water above 273.15 K, with calculations taken from "Landolt-Bornstein,
    ///< 1987, Numerical Data and Functional relationships in Science and
    ///< Technology. Group V/Vol 4B Meteorology. Physical and Chemical properties of
    ///< Air, P35."
  GoffGratchLandoltBornsteinWater,  ///< For humidity only: use the default
    ///< methods, calculating saturation mixing ratio from temperature using the
    ///< Goff-Gratch formulae over water, with calculations taken from
    ///< "Landolt-Bornstein, 1987, Numerical Data and Functional relationships
    ///< in Science and Technology. Group V/Vol 4B Meteorology. Physical and
    ///< Chemical properties of Air, P35."
  Rogers,  ///< For humidity only: use the default methods, calculating
    ///< saturation mixing ratio from temperature using classical formula from
    ///< Rogers and Yau (1989; Eq2.17)
  DEFAULT  ///< Default methods with default formulations (where applicable)
};

/*! Various Formulations available - Specified by author */
enum class Formulation {
  DEFAULT, /*!< DEFAULT formulation - see descriptions of formulae */

  // Humidity Formulations: Specific authors
  Murphy,
  Sonntag,
  GoffGratchLandoltBornsteinIceWater,
  GoffGratchLandoltBornsteinWater,
  Gill,
  GillUKMO,
  Walko,
  Rogers,
  Clark2008,
  UKMOUnififiedModelOPS,

  // Other Formulations
  NCARRAL,
  ICAO,
  LarocheSarrazin,
};

Method resolveMethods(const std::string& method);

// -------------------------------------------------------------------------------------
/*!
* \brief Calculates saturated vapour pressure from temperature
*
* \b Formulation \b available:
*      - Sonntag:
*        Calculation is using the Eqn 7 Sonntag (1994)
*        Reference: "Sonntag, D., Advancements in the field of hygrometry,
*        Meteorol. Zeitschrift, N. F., 3, 51-66, 1994." .
*      - Walko:
*        Polynomial fit of Goff-Gratch (1946) formulation. (Walko, 1991)
*      - Murphy:
*        Alternative method to Walko 1991 (costs more CPU, more accurate)
*        Reference: Eq. 10 of "Murphy and Koop, Review of the vapour pressure of ice and
*        supercooled water for atmospheric applications, Q. J. R.
*        Meteorol. Soc (2005), 131, pp. 1539-1565."
*      - GoffGratchLandoltBornsteinIceWater:
*        Using the Goff-Gratch formulae over ice below and including 273.15 K
*        (0 C) and over water above 273.15 K, with calculations taken from
*        "Landolt-Bornstein, 1987, Numerical Data and Functional relationships
*        in Science and Technology. Group V/Vol 4B Meteorology. Physical and
*        Chemical properties of Air, P35."
*      - GoffGratchLandoltBornsteinWater:
*        As above but over water for all temperatures.
*      - Rogers:
*        Classical formula from Rogers and Yau (1989; Eq2.17)
*      - DEFAULT:
*        Uses Rogers formulation.
*
*
* \param temp_K
*     Temperature [k]
* \return saturated vapour pressure
*/
float SatVaporPres_fromTemp(const float temp_K,
                   const Formulation formulation = Formulation::DEFAULT);


// -------------------------------------------------------------------------------------
/*!
* \brief Calculates a saturation vapour pressure for moist air from the
* saturation vapour pressure of pure water vapour at a give temperature and
* pressure.
*
* \b Formulation \b available:
*      - Gill:
*        Corrects the saturation vapour pressure of pure water vapour
*        using an enhancement factor needed for moist air.
*        Reference: Eq. A4.6 of Gill (1982) "Atmosphere-Ocean Dynamics",
*        Academic Press. This approximates table 89 of the Smithsonian
*        Meteorological Tables correct to 2 parts in 10^4.
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
                        const Formulation formulation = Formulation::Gill);
// -------------------------------------------------------------------------------------
/*!
* \brief Calculates Saturated specific humidity or saturated vapour pressure using
* saturation vapour pressure.
*
*
* \b Formulation \b available:
*      - GillUKMO:
*        Qsat is given by rearranging equation A4.3 in Gill (1982) with the vapour
*        pressure `e'` replaced by `Psat` and specific humidity `q` by `Qsat`:
*           `Qsat = epsilon Psat / (P - (1 - epsilon) Psat)`
*        where `epsilon = 0.62198` (the ratio of molecular weights of water and
*        dry air). To avoid asymptotic behaviour for `P < Psat`, the denominator
*        is modified to ensure `QSat = 1 kg/kg` for `P < Psat` (this appears to be
*        UKMO specific hence the formulation name).
*
* \param Psat
*      saturation vapour pressure of pure water vapour
* \param P
*      Pressure
* \return Saturated specific humidity or saturated vapour pressure
*/
float Qsat_From_Psat(float Psat, float P,
                     Formulation formulation = Formulation::GillUKMO);

// -------------------------------------------------------------------------------------
/*!
* \brief Derive Virtual Temperature from saturation vapour pressure, pressure and temperature
*
* \b Formulation \b available:
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
                          Formulation formulation = Formulation::DEFAULT);

// -------------------------------------------------------------------------------------
/*!
* \brief Derive Virtual Tempreture using Relative humidity, sat. vapour pressure, pressure
* and temperature
*
* \b Formulation \b available:
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
                            Formulation formulation = Formulation::DEFAULT);

// -------------------------------------------------------------------------------------
/*!
* \brief Derive Virtual Temperature using Specific Humidity, and air temperature.
* Formula used here is from variabletransforms, cal_VirtualTemperature method
* in the UFO code base.
*
* \param Sh
*     Specific humidity
* \param AT
*     Air Temperature
*
* \return Virtual Temperature
*/
float VirtualTemp_From_Sh_AT(float Sh, float At);

// -------------------------------------------------------------------------------------
/*!
* \brief Converts height to pressure using the International Civil Aviation Organization
* (ICAO) atmosphere.
*
* \b Formulation \b available:
*      - ICAO: using ICAO standard
*
* \param height
*     observation height in geopotential metres [gpm]
* \return pressure
*/
float Height_To_Pressure_ICAO_atmos(float Height,
                            Formulation formulation = Formulation::ICAO);

// -------------------------------------------------------------------------------------
/*!
* \brief Converts pressure to height.
*
* \b Formulation \b available:
*      - NCARRAL: NCAR-RAL is a fast approximation for pressures > 120 hPa.
*        Below 120hPa (~15km) use the ICAO atmosphere.
*      - ICAO: uses the ICAO atmosphere for all pressures (default)
*
* \param pressure
*     observation pressure in Pa
* \return height in geopotential metres
*/
float Pressure_To_Height(float pressure,
                         Formulation formulation = Formulation::ICAO);

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
*      - LarocheSarrazin:
*        Calculation is using Eq. 4. of Laroche and Sarrazin (2013).
*        Reference: "Laroche, S. and Sarrazin, R.,
*                   Impact of Radiosonde Balloon Drift on
*                   Numerical Weather Prediction and Verification,
*                   Weather and Forecasting, 28(3), 772-782, 2013."
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
 Formulation formulation = Formulation::LarocheSarrazin,
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

/**
 * @brief Calculates the standard dry aerosol mass mixing ratio (kg kg^-1).
 *
 * \b Formulation \b available:
 *    - Clark2008:
 *      Eq. 1 in Clark et al. (2008).
 *
 * @param rho_aerosol The aerosol density (kg m^-3).
 * @param rho_air The air density (kg m^-3).
 * @param N_0 The standard aerosol number density (m^-3).
 * @param r_0 The standard aerosol radius (m).
 * @param formulation The formulation to use.
 * @return The calculated standard dry aerosol mass mixing ratio m_0 (kg kg^-1).
 */
float dryAerosolMassMixingRatio(
    const float rho_aerosol, const float rho_air, const float N_0,
    const float r_0, const Formulation formulation = Formulation::Clark2008);

/**
 * @brief Calculates the mean dry aerosol radius (m) which is assumed to
 * vary as a power of the aerosol mass concentration.
 *
 * \b Formulation \b available:
 *    - Clark2008:
 *      Eq. 2 in Clark et al. (2008).
 *
 * @param r_0 The standard aerosol radius (m).
 * @param m The aerosol mass concentration (kg kg^-1).
 * @param m_0 The standard dry aerosol mass mixing ratio (kg kg^-1).
 * @param p The aerosol size distribution power law exponent.
 * @param formulation The formulation to use.
 * @return The calculated mean dry aerosol radius r_md (m).
 */
float dryMeanParticleVolumeRadius(
    const float r_0, const float m, const float m_0, const float p,
    const Formulation formulation = Formulation::Clark2008);

/**
 * @brief Calculates the aerosol number density (m^-3).
 *
 * \b Formulation \b available:
 *    - Clark2008:
 *      Eq. 3 in Clark et al. (2008).
 *
 * @param N_0 The standard aerosol number density (m^-3).
 * @param m The aerosol mass concentration (kg kg^-1).
 * @param m_0 The standard dry aerosol mass mixing ratio (kg kg^-1).
 * @param p The aerosol size distribution power law exponent.
 * @param formulation The formulation to use.
 * @return The calculated aerosol number density N (m^-3).
 */
float aerosolNumberDensity(
    const float N_0, const float m, const float m_0, const float p,
    const Formulation formulation = Formulation::Clark2008);

/**
 * @brief Calculates the approximate activation droplet radius (m), which is
 * the radius at which the peak of the Koehler curve occurs.
 *
 * \b Formulation \b available:
 *    - Clark2008:
 *      Eq. 11 in Clark et al. (2008).
 *
 * @param A The Koehler curve constant A (m).
 * @param B The Koehler curve constant B (unitless).
 * @param r_md The dry aerosol radius (m).
 * @param formulation The formulation to use.
 * @return The calculated approximate activation droplet radius r_act (m).
 */
float approximateActivationDropletRadius(
    const float A, const float B, const float r_md,
    const Formulation formulation = Formulation::Clark2008);

/**
 * @brief Calculates the approximate aerosol extinction coefficient factor
 * (unitless) as a function of the extinction efficiency and particle size.
 *
 * \b Formulation \b available:
 *    - Clark2008:
 *      Eq. 18 in Clark et al. (2008).
 *
 * @param Q The extinction efficiency factor (unitless).
 * @param eta The particle size weighting factor (unitless).
 * @param formulation The formulation to use.
 * @return The calculated approximate aerosol extinction coefficient factor
 * beta_0 (unitless).
 */
float approximateAerosolExtinctionCoefficientFactor(
    const float Q, const float eta,
    const Formulation formulation = Formulation::Clark2008);

/**
 * @brief Calculates the cloud cover fraction (unitless) as a function of the
 * relative humidity.
 *
 * \b Formulation \b available:
 *    - UKMOUnififiedModelOPS:
 *      Corresponds to the RH_TO_CC diagnostic taken from the UK Met Office
 *      Unified Model, as used in the Met Office Observation Processing System
 *      (OPS).
 *
 * @param rh The relative humidity (fraction).
 * @param rh_crit The critical relative humidity (fraction).
 * @param ccp1 The cloud cover parameter 1 (unitless) - corresponds to PC1 in
 * RH_TO_CC.
 * @param ccp2 The cloud cover parameter 2 (unitless) - corresponds to PC2 in
 * RH_TO_CC.
 * @param ccp3 The cloud cover parameter 3 (unitless) - corresponds to PC3 in
 * RH_TO_CC.
 * @param formulation The formulation to use.
 * @return The calculated cloud cover fraction (unitless).
 */
float cloudCoverFraction(
    float rh, const float rh_crit, const float ccp1, const float ccp2,
    const float ccp3,
    const Formulation formulation = Formulation::UKMOUnififiedModelOPS);

/**
 * @brief Calculates a total water relative humidity (unitless) as a
 * function of the cloud cover fraction.
 *
 * \b Formulation \b available:
 *    - UKMOUnififiedModelOPS:
 *      Corresponds to the CC_TO_RHTOT diagnostic taken from the Met Office
 *      Unified Model, as used in the Met Office Observation Processing System
 *      (OPS) where total water relative humidity is associated with cloud
 *      cover via a piecewise function, derived from a triangular distribution.
 *      This is consistent with the Smith (1990) cloud scheme
 *      (https://doi.org/10.1002/qj.49711649210).
 *
 * @param cc The cloud cover fraction (unitless).
 * @param rh_crit The critical relative humidity: that at which liquid water
 * droplets are considered to form (fraction).
 * @param rh_tot_p1 The total water relative humidity parameter 1 - The
 * cloud cover fraction at which the piecewise function transitions.
 * @param rh_tot_p2 The scaling to apply to the cloud cover fraction within
 * the first part of the piecewise function.
 * @param rh_tot_p3 The scaling to apply to (1 - cc), where cc is the cloud
 * cover fraction, within the second part of the piecewise function.
 * @param formulation The formulation to use.
 * @return The calculated total water relative humidity (unitless).
 */
float totalWaterRelativeHumidity(
    float cc, const float rh_crit, const float rh_tot_p1, const float rh_tot_p2,
    const float rh_tot_p3,
    const Formulation formulation = Formulation::UKMOUnififiedModelOPS);

/**
 * @brief Calculates the equilibrium relative humidity (fraction) as a
 * function of the normalized droplet radius, expressed as a Koehler curve.
 *
 * \b Formulation \b available:
 *    - Clark2008:
 *      Corresponds to Eq. 10 in Clark et al. (2008).
 *
 * @param g The normalized droplet radius (r_m / r_md) (m).
 * @param r_md The dry aerosol radius (m).
 * @param A The Koehler curve constant A (m).
 * @param B The Koehler curve constant B (unitless).
 * @param formulation The formulation to use.
 * @return The calculated equilibrium relative humidity (fraction).
 */
float equilibriumRelativeHumidity(
    const float g, const float r_md, const float A, const float B,
    const Formulation formulation = Formulation::Clark2008);

/**
 * @brief Calculates the derivative, with respect to the normalized droplet
 * radius, of the equilibrium relative humidity expressed as a Koehler curve.
 *
 * @see equilibriumRelativeHumidity
 *
 * \b Formulation \b available:
 *    - Clark2008:
 *      Corresponds to the derivative Eq. 10 in Clark et al. (2008).
 *
 * @param g The normalized droplet radius (r_m / r_md) (m m^-1).
 * @param r_md The dry aerosol radius (m).
 * @param A The Koehler curve constant A (m).
 * @param B The Koehler curve constant B (unitless).
 * @param formulation The formulation to use.
 * @return The calculated derivative of the relative humidity with respect to g
 * (m m^-1).
 */
float equilibriumRelativeHumidityDerivative(
    const float g, const float r_md, const float A, const float B,
    const Formulation formulation = Formulation::Clark2008);

/**
 * @brief Calculates the specific cloud liquid water content (kg kg^-1) as a
 * function of the normalized droplet radius.
 *
 * \b Formulation \b available:
 *    - Clark2008:
 *      Corresponds to Eq. 13 in Clark et al. (2008).
 *
 * @param g The normalized droplet radius (r_m / r_md) (m m^-1).
 * @param r_md The dry aerosol radius (m).
 * @param N The number concentration of the aerosol particles (m^-3).
 * @param rho_wat The density of water (kg m^-3).
 * @param formulation The formulation to use.
 * @return The calculated liquid water content q_l (kg kg^-1).
 */
float specificCloudWaterContent(
    const float g, const float r_md, const float N, const float rho_wat,
    const Formulation formulation = Formulation::Clark2008);

/**
 * @brief Calculates the derivative, with respect to the normalized droplet
 * radius, of the specific cloud water content.
 *
 * @see specificCloudWaterContent
 *
 * \b Formulation \b available:
 *    - Clark2008:
 *      Corresponds to the derivative Eq. 13 in Clark et al. (2008).
 *
 * @param g The normalized droplet radius (r_m / r_md) (m m^-1).
 * @param r_md The dry aerosol radius (m).
 * @param N The number concentration of the aerosol particles (m^-3).
 * @param rho_wat The density of water (kg m^-3).
 * @param formulation The formulation to use.
 * @return The calculated derivative of the liquid water content with respect to
 * g (kg kg^-1).
 */
float specificCloudWaterContentDerivative(
    const float g, const float r_md, const float N, const float rho_wat,
    const Formulation formulation = Formulation::Clark2008);

/**
 * @brief Calculates an approximate aerosol extinction coefficient (m^-1) as
 * a function of the number density and radius of (potentially wet) aerosol
 * particles.
 *
 * \b Formulation \b available:
 *    - Clark2008:
 *      Corresponds to Eq. 17 in Clark et al. (2008).
 *
 * @param beta_0 An aerosol extinction coefficient factor (unitless).
 * @param N The number density of the aerosol particles (m^-3).
 * @param r_m The aerosol particle radius (m).
 * @param formulation The formulation to use.
 * @return The calculated aerosol extinction coefficient beta (m^-1).
 */
float aerosolExtinctionCoefficient(
    const float beta_0, const float N, const float r_m,
    const Formulation formulation = Formulation::Clark2008);

/**
 * @brief Calculates a visibility (m) from a liminal contrast and aerosol
 * extinction coefficient, ignoring the effects of clean air.
 *
 * \b Formulation \b available:
 *    - Clark2008:
 *      Corresponds to Eq. 20 in Clark et al. (2008).
 *
 * @param epsilon The liminal contrast factor (unitless).
 * @param beta The aerosol extinction coefficient (m^-1).
 * @param formulation The formulation to use.
 * @return The calculated aerosol visibility (m).
 */
float aerosolVisibility(const float epsilon, const float beta,
                        const Formulation formulation = Formulation::Clark2008);

}  // namespace formulas
}  // namespace ufo

#endif  // UFO_VARIABLETRANSFORMS_FORMULAS_H_
