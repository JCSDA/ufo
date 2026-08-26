/*
 * (C) Crown copyright 2020, Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ioda/ObsSpace.h"
#include "ufo/utils/Constants.h"
#include "ufo/variabletransforms/Cal_Humidity.h"

namespace ufo {
/**************************************************************************************************/
//  Cal_RelativeHumidity
/**************************************************************************************************/

static TransformMaker<Cal_RelativeHumidity>
    makerCal_RelativeHumidity_("RelativeHumidity");

Cal_RelativeHumidity::Cal_RelativeHumidity(
    const Parameters_ &options,
    const ObsFilterData &data,
    ioda::ObsDataVector<int> &flags,
    ioda::ObsDataVector<float> &obserr)
  : TransformBase(options, data, flags, obserr),
      relativeHumidityUnits_(options.RelativeHumidityUnits),
      allowSuperSaturation_(options.AllowSuperSaturation),
      specifichumidityvariable_(options.SpecificHumidityVariable),
      pressurevariable_(options.PressureVariable),
      pressureat2mvariable_(options.PressureAt2MVariable),
      pressuregroupvariable_(options.PressureGroupVariable),
      temperaturevariable_(options.TemperatureVariable),
      temperatureat2mvariable_(options.TemperatureAt2MVariable),
      relativehumidityvariable_(options.RelativeHumidityVariable),
      relativehumidityat2mvariable_(options.RelativeHumidityAt2MVariable),
      watervapormixingratiovariable_(options.WaterVaporMixingRatioVariable),
      dewpointtemperaturevariable_(options.DewPointTemperatureVariable),
      dewpointtemperatureat2mvariable_(options.DewPointTemperature2MVariable)
      {}

/**************************************************************************************************/

void Cal_RelativeHumidity::runTransform(const std::vector<bool> &apply) {
  oops::Log::trace() << " --> Retrieve Relative humidity"
            << std::endl;
  oops::Log::trace() << "      --> method: " << options_.Method.value() << std::endl;
  oops::Log::trace() << "      --> obsName: " << obsName() << std::endl;

  // Get the right method
  switch (method()) {
    case formulas::Method::UKMO:
      methodUKMO(apply);
      break;
    case formulas::Method::UKMOmixingratio:
      methodUKMOmixingratio(apply);
      break;
    case formulas::Method::NCAR:
      methodDEFAULT(apply, formulas::Formulation::Rogers);
      break;
    case formulas::Method::NOAA:
      methodDEFAULT(apply, formulas::Formulation::Rogers);
      break;
    case formulas::Method::Sonntag:
      methodDEFAULT(apply, formulas::Formulation::Sonntag);
      break;
    case formulas::Method::Walko:
      methodDEFAULT(apply, formulas::Formulation::Walko);
      break;
    case formulas::Method::Murphy:
      methodDEFAULT(apply, formulas::Formulation::Murphy);
      break;
    case formulas::Method::GoffGratchLandoltBornsteinIceWater:
      methodDEFAULT(apply,
                    formulas::Formulation::GoffGratchLandoltBornsteinIceWater);
      break;
    case formulas::Method::GoffGratchLandoltBornsteinWater:
      methodDEFAULT(apply,
                    formulas::Formulation::GoffGratchLandoltBornsteinWater);
      break;
    case formulas::Method::Rogers:
      methodDEFAULT(apply, formulas::Formulation::Rogers);
      break;
    default:
      methodDEFAULT(apply, formulas::Formulation::Rogers);
      break;
  }
}

/**************************************************************************************************/
/*
Calculates relative humidity (RH_ice) from dew point temperature,
or converts RH_water to RH_ice

Method: -
  Saturated specific humidity at the dew point (ie w.r.t water), and saturated
  specific humidity (over ice below 0C) at the air temperature are calculated.
  Relative humidity is then calculated using :

         RH = (QSAT(DEW POINT)/QSAT(DRY BULB))*relativeHumidityAtSaturation()

   For some temperatures (e.g. when dew point = temperature), supersaturation
   w.r.t ice may occur.
   The option AllowSuperSaturation (false by default) controls whether upper air relative humidity
   is capped at 100%.
   If the pressure, dew point or temperature take extreme
   values or are missing, the relative humidity is set to missing data.
*/
void Cal_RelativeHumidity::methodUKMO(const std::vector<bool> &apply) {
  const size_t nlocs = obsdb_.nlocs();

  float Q_sub_s_w_ice, Q_sub_s_w;
  float pressure, temperature, dewPoint;
  std::vector<float> airTemperature;
  std::vector<float> dewPointTemperature;
  std::vector<float> relativeHumidity;
  std::vector<float> airPressure;
  bool surfaceData = true;

  // Here we can only use data that have not been rejected by quality control
  // so making sure UseValidDataOnly_ is set to True
  SetUseValidDataOnly(true);

  // 0. Innitialise the ouput array
  // -----------------------------------------------------------------------------------------------

  // 1. get the right variables
  // -----------------------------------------------------------------------------------------------
  // Compulsory surface observation
  //     First looking for surface observation
  //     Then looking for upperair data
  if (obsdb_.has("ObsValue", pressureat2mvariable_) &&
      obsdb_.has("ObsValue", temperatureat2mvariable_) &&
      obsdb_.has("ObsValue", dewpointtemperatureat2mvariable_)) {
    getObservation("ObsValue", pressureat2mvariable_,
                   airPressure, true);
    getObservation("ObsValue", temperatureat2mvariable_,
                   airTemperature, true);
    getObservation("ObsValue", dewpointtemperatureat2mvariable_,
                   dewPointTemperature, true);
    getObservation("ObsValue", relativehumidityat2mvariable_,
                   relativeHumidity);
  } else {
    getObservation(pressuregroupvariable_, pressurevariable_,
                   airPressure, true);
    getObservation("ObsValue", temperaturevariable_,
                   airTemperature, true);
    getObservation("ObsValue", dewpointtemperaturevariable_,
                   dewPointTemperature, true);
    getObservation("ObsValue", relativehumidityvariable_,
                   relativeHumidity);
    surfaceData = false;
  }
  if (relativeHumidity.empty()) {
    relativeHumidity.assign(nlocs, missingValueFloat);
  }

  // 2. making sure we have what we need is here
  // -----------------------------------------------------------------------------------------------
  if (!oops::allVectorsSameSize(airPressure, airTemperature, dewPointTemperature)) {
    oops::Log::warning() << "Vector sizes: "
                         << oops::listOfVectorSizes(airPressure, airTemperature,
                                                    dewPointTemperature)
                         << std::endl;
    throw eckit::BadValue("At least one variable vector is the wrong size or empty out of "
                          "pressure, air temperature and dew point temperature", Here());
  }

  // Lambda function to calcualate saturation specific humidity: the output
  // variables are placed into Q_sub_s_w (calculated over water for all
  // temperatures) and Q_sub_s_w_ice (calculated over water above 0C and over
  // ice below and including 0C).
  // temp_1: airTemprature or dewPointTemperature
  // temp_2: airTemperature
  // -----------------------------------------------------------------------------------------------
  auto evaluateSatSpecHumidity = [&](float temp_1, float temp_2) {
    float e_sub_s_w, e_sub_s_w_ice;
    // sat. vapor pressure from Dewpoint temperature - over water.
    e_sub_s_w = formulas::SatVaporPres_fromTemp(
        temp_1, formulas::Formulation::Sonntag);
    e_sub_s_w = formulas::SatVaporPres_correction(e_sub_s_w, temp_1, pressure,
                                                  formulas::Formulation::Gill);
    // Convert sat. vapor pressure to sat. specific humidity (both over water)
    Q_sub_s_w = formulas::Qsat_From_Psat(e_sub_s_w, pressure, formulas::Formulation::GillUKMO);

    // sat. vapor pressure from Drybulb temperature - over ice if temp2 <= 0C
    e_sub_s_w_ice = formulas::SatVaporPres_fromTemp(
        temp_2, formulas::Formulation::GoffGratchLandoltBornsteinIceWater);
    e_sub_s_w_ice = formulas::SatVaporPres_correction(
        e_sub_s_w_ice, temp_2, pressure, formulas::Formulation::Gill);
    // Convert sat. vapor pressure to saturated specific humidity (both over
    // ice if temp2 < 0C)
    Q_sub_s_w_ice = formulas::Qsat_From_Psat(
        e_sub_s_w_ice, pressure, formulas::Formulation::GillUKMO);
  };

  // 3. Loop over each record
  // -----------------------------------------------------------------------------------------------
  for (ioda::ObsSpace::RecIdxIter irec = obsdb_.recidx_begin();
       irec != obsdb_.recidx_end(); ++irec) {
    const std::vector<std::size_t> &rSort = obsdb_.recidx_vector(irec);

    // 3.1 Loop over each record
    for (size_t iloc : rSort) {
      // if the data have been excluded by the where statement
      if (!apply[iloc]) continue;

      // store some variables
      pressure = airPressure[iloc];
      temperature = airTemperature[iloc];
      dewPoint = dewPointTemperature[iloc];

      // There is very little sensitivity of calculated RH to P
      // (less than 0.1% to change from 1000 to 800 hPa)
      // --> so for surface observation that do not have any airPressure
      //     we set it to 1000 hPa.
      if (pressure == missingValueFloat && surfaceData) {
        pressure = 100000.0;  // default pressure in Pascal
      }

      // if dewpoint temperature is reported (most stations)
      if (dewPoint > 1.0 && temperature != missingValueFloat && pressure > 1.0) {
        // calculate saturation specific humidity from dewpoint (always over
        // water) and from temperature (over water if temperature > 0C and over
        // ice if temperature <= 0C).
        evaluateSatSpecHumidity(dewPoint, temperature);

        // if saturated specific humidities are positive calculate relative
        // humidity
        if (Q_sub_s_w > 0 && Q_sub_s_w_ice > 0) {
          relativeHumidity[iloc] =
              (Q_sub_s_w / Q_sub_s_w_ice) * relativeHumidityAtSaturation();
          if (!allowSuperSaturation_)
            relativeHumidity[iloc] = std::min(relativeHumidityAtSaturation(),
                                              relativeHumidity[iloc]);
        }
      // if relative humidity (Rh) is reported (small minority of stations)
      // update from Rh over water (reported) to Rh over ice if necessary
      } else if (relativeHumidity[iloc] != missingValueFloat) {
        if (temperature == missingValueFloat) {
          relativeHumidity[iloc] = missingValueFloat;
        } else if (temperature < ufo::Constants::t0c) {
          // calculate saturation specific humidity over water and ice given
          // this temperature to allow adjustment
          evaluateSatSpecHumidity(temperature, temperature);
          relativeHumidity[iloc] *= (Q_sub_s_w / Q_sub_s_w_ice);
          if (!allowSuperSaturation_)
            relativeHumidity[iloc] = std::min(relativeHumidityAtSaturation(),
                                              relativeHumidity[iloc]);
        }
      }
    }
  }

  // assign the derived relative humidity as DerivedObsValue
  if (surfaceData) {
    putObservation(relativehumidityat2mvariable_, relativeHumidity);
  } else {
    putObservation(relativehumidityvariable_, relativeHumidity);
  }
}

/**************************************************************************************************/
/*
Approximate relative humidity as the mixing ratio divided by the saturation specific humidity

Method: -
  Saturated specific humidity at the air temperature is obtained
  above water for temperatures greater than 0 and above ice below zero.
  Relative humidity is then calculated using :

    RH = (Mixing_Ratio/QSAT(DRY BULB))*relativeHumidityAtSaturation()

  The option AllowSuperSaturation (false by default) controls whether upper air relative humidity
  is capped at relativeHumidityAtSaturation() (1 or 100 depending on chosen units).

  Example:
  \code{.unparsed}
  obs filter:
  - filter: Variable Transforms
    Transform: RelativeHumidity
    observation relative humidity units: percentage
    Method: UKMOmixingratio
    AllowSuperSaturation: false
  \endcode
*/
void Cal_RelativeHumidity::methodUKMOmixingratio(const std::vector<bool> &apply) {
  const size_t nlocs = obsdb_.nlocs();

  std::vector<float> airTemperature;
  std::vector<float> mixingRatio;
  std::vector<float> relativeHumidity;
  std::vector<float> airPressure;

  // Here we can only use data that has not been rejected by quality control
  // so making sure UseValidDataOnly_ is set to True
  SetUseValidDataOnly(true);

  // Get variables
  getObservation(pressuregroupvariable_, pressurevariable_,
                 airPressure, true);
  getObservation("ObsValue", temperaturevariable_,
                 airTemperature, true);
  getObservation("ObsValue", watervapormixingratiovariable_,
                 mixingRatio, true);
  getObservation("ObsValue", relativehumidityvariable_,
                 relativeHumidity);

  if (relativeHumidity.empty()) {
    relativeHumidity.assign(nlocs, missingValueFloat);
  }

  if (!oops::allVectorsSameSize(airPressure, airTemperature, mixingRatio)) {
    oops::Log::warning() << "Vector sizes: "
                         << oops::listOfVectorSizes(airPressure, airTemperature,
                                                    mixingRatio)
                         << std::endl;
    throw eckit::BadValue("At least one variable vector is the wrong size or empty out of "
                          "pressure, air temperature and mixing ratio", Here());
  }

  for (ioda::ObsSpace::RecIdxIter irec = obsdb_.recidx_begin();
       irec != obsdb_.recidx_end(); ++irec) {
    const std::vector<std::size_t> &rSort = obsdb_.recidx_vector(irec);

    // Loop over each record
    for (size_t iloc : rSort) {
      if (!apply[iloc]) continue;

      // Store some variables
      const float pressure = airPressure[iloc];
      const float temperature = airTemperature[iloc];
      const float mixRatio = mixingRatio[iloc];

      if (pressure > 0  &&
          pressure != missingValueFloat &&
          temperature != missingValueFloat &&
          mixRatio != missingValueFloat) {
        // Sat. vapor pressure from drybulb temperature and pressure
        float e_sub_s = formulas::SatVaporPres_fromTemp(
            temperature, formulas::Formulation::GoffGratchLandoltBornsteinIceWater);
        e_sub_s = formulas::SatVaporPres_correction(
            e_sub_s, temperature, pressure, formulas::Formulation::Gill);
        // Convert sat. vapor pressure to sat. specific humidity
        const float Q_sub_s =
            formulas::Qsat_From_Psat(e_sub_s, pressure, formulas::Formulation::GillUKMO);

        // Calculate RH
        if (mixRatio >= 0 && Q_sub_s > 0) {
          relativeHumidity[iloc] =
              (mixRatio / Q_sub_s) * relativeHumidityAtSaturation();
        } else {
          relativeHumidity[iloc] = missingValueFloat;
        }
        if (relativeHumidity[iloc] != missingValueFloat &&
            !allowSuperSaturation_) {
          relativeHumidity[iloc] =
              std::min(relativeHumidityAtSaturation(), relativeHumidity[iloc]);
        }
      }
    }
  }
  // Assign the derived relative humidity as DerivedObsValue
  putObservation(relativehumidityvariable_, relativeHumidity);
}

/**************************************************************************************************/

void Cal_RelativeHumidity::methodDEFAULT(
    const std::vector<bool> &apply,
    formulas::Formulation SatVaporPres_fromTemp_form) {
  const size_t nlocs = obsdb_.nlocs();

  float esat, qvs, qv, satVaporPres;

  std::vector<float> specificHumidity;
  std::vector<float> airTemperature;
  std::vector<float> pressure;
  std::vector<float> relativeHumidity(nlocs);

  getObservation("ObsValue", specifichumidityvariable_,
                 specificHumidity, true);
  getObservation("ObsValue", temperaturevariable_,
                 airTemperature, true);
  getObservation(pressuregroupvariable_, pressurevariable_,
                 pressure, false);
  if (pressure.empty()) {
    getObservation("ObsValue", pressureat2mvariable_,
                   pressure, true);
  }

  if (!oops::allVectorsSameSize(specificHumidity, airTemperature, pressure)) {
    oops::Log::warning() << "Vector sizes: "
                         << oops::listOfVectorSizes(specificHumidity, airTemperature,
                                                    pressure)
                         << std::endl;
    throw eckit::BadValue("At least one variable vector is the wrong size or empty out of "
                          "specific humidity, air temperature and pressure", Here());
  }

  // Initialise this vector with missing value
  relativeHumidity.assign(nlocs, missingValueFloat);

  // Loop over all obs
  for (size_t jobs = 0; jobs < nlocs; ++jobs) {
    // if the data have been excluded by the where statement
    if (!apply[jobs]) continue;

    if (specificHumidity[jobs] != missingValueFloat &&
        airTemperature[jobs] != missingValueFloat && pressure[jobs] != missingValueFloat) {
      // Calculate saturation vapor pressure from temperature according to requested formulation
      // Double-check result is always lower than 15% of incoming pressure.
      satVaporPres = formulas::SatVaporPres_fromTemp(
          airTemperature[jobs], SatVaporPres_fromTemp_form);
      esat = std::min(pressure[jobs]*0.15f, satVaporPres);

      // Convert sat. vapor pressure to sat water vapor mixing ratio
      qvs = Constants::rd_over_rv * esat/(pressure[jobs]-esat);

      // Convert specific humidity to water vapor mixing ratio
      qv = std::max(1.0e-12f, specificHumidity[jobs]/(1.0f-specificHumidity[jobs]));

      // Final RH (which can be greater than fully saturated) is q/qsat, but set
      // sensible lowest limit
      relativeHumidity[jobs] = std::max(1.0e-6f, qv/qvs) * relativeHumidityAtSaturation();
    }
  }
  putObservation(relativehumidityvariable_, relativeHumidity);
}

/************************************************************************************/
//  Cal_SpecificHumidity
/************************************************************************************/
static TransformMaker<Cal_SpecificHumidity>
    makerCal_SpecificHumidity_("SpecificHumidity");

Cal_SpecificHumidity::Cal_SpecificHumidity(
    const Parameters_ &options,
    const ObsFilterData &data,
    ioda::ObsDataVector<int> &flags,
    ioda::ObsDataVector<float> &obserr)
    : TransformBase(options, data, flags, obserr),
      relativeHumidityUnits_(options.RelativeHumidityUnits),
      specifichumidityvariable_(options.SpecificHumidityVariable),
      pressurevariable_(options.PressureVariable),
      pressureat2mvariable_(options.PressureAt2MVariable),
      pressuregroupvariable_(options.PressureGroupVariable),
      temperaturevariable_(options.TemperatureVariable),
      relativehumidityvariable_(options.RelativeHumidityVariable),
      dewpointtemperaturevariable_(options.DewPointTemperatureVariable)
{}

/************************************************************************************/

void Cal_SpecificHumidity::runTransform(const std::vector<bool> &apply) {
  oops::Log::trace() << " Retrieve Specific Humidity" << std::endl;
  oops::Log::trace() << "      --> method: " << options_.Method.value() << std::endl;
  oops::Log::trace() << "      --> obsName: " << obsName() << std::endl;

  switch (method()) {
    case formulas::Method::NCAR:
      methodDEFAULT(apply, formulas::Formulation::Rogers);
      break;
    case formulas::Method::NOAA:
      methodDEFAULT(apply, formulas::Formulation::Rogers);
      break;
    case formulas::Method::UKMO:
      methodDEFAULT(apply, formulas::Formulation::Sonntag);
      break;
    case formulas::Method::UKMOmixingratio:
      methodDEFAULT(apply,
                    formulas::Formulation::GoffGratchLandoltBornsteinIceWater);
      break;
    case formulas::Method::UKMOQsatWater:
      methodQSat(apply, formulas::Formulation::GoffGratchLandoltBornsteinWater);
      break;
    case formulas::Method::UKMOQsatIceWater:
      methodQSat(apply,
                 formulas::Formulation::GoffGratchLandoltBornsteinIceWater);
      break;
    case formulas::Method::Sonntag:
      methodDEFAULT(apply, formulas::Formulation::Sonntag);
      break;
    case formulas::Method::Walko:
      methodDEFAULT(apply, formulas::Formulation::Walko);
      break;
    case formulas::Method::Murphy:
      methodDEFAULT(apply, formulas::Formulation::Murphy);
      break;
    case formulas::Method::GoffGratchLandoltBornsteinIceWater:
      methodDEFAULT(apply,
                    formulas::Formulation::GoffGratchLandoltBornsteinIceWater);
    case formulas::Method::GoffGratchLandoltBornsteinWater:
      methodDEFAULT(apply,
                    formulas::Formulation::GoffGratchLandoltBornsteinWater);
      break;
    case formulas::Method::Rogers:
      methodDEFAULT(apply, formulas::Formulation::Rogers);
      break;
    default:
      methodDEFAULT(apply, formulas::Formulation::DEFAULT);
      break;
  }
}

/************************************************************************************/

void Cal_SpecificHumidity::methodQSat(
    const std::vector<bool> &apply,
    formulas::Formulation SatVaporPres_fromTemp_form) {
  const size_t nlocs = obsdb_.nlocs();
  std::vector<float> relativeHumidity;
  std::vector<float> airTemperature;
  std::vector<float> pressure;
  std::vector<float> specificHumidity(nlocs);
  /// FUTURE: Allow dewpoint temperature to be used as an alternative to
  /// relative humidity

  getObservation("ObsValue", relativehumidityvariable_, relativeHumidity, true);
  getObservation("ObsValue", temperaturevariable_, airTemperature, true);
  if (obsdb_.has("ObsValue", pressureat2mvariable_)) {
    getObservation("ObsValue", pressureat2mvariable_, pressure, true);
  } else {
    getObservation(pressuregroupvariable_, pressurevariable_, pressure, false);
  }
  if (!oops::allVectorsSameSize(relativeHumidity, airTemperature, pressure)) {
    oops::Log::warning() << "Vector sizes: "
                         << oops::listOfVectorSizes(relativeHumidity,
                                                    airTemperature, pressure)
                         << std::endl;
    throw eckit::BadValue(
        "At least one variable vector is the wrong size or empty out of "
        "relative humidity, air temperature and pressure",
        Here());
  }
  specificHumidity.assign(nlocs, missingValueFloat);

  for (size_t jobs = 0; jobs < nlocs; ++jobs) {
    if (!apply[jobs]) {
      continue;  // data have been excluded by the where statement
    }
    if (pressure[jobs] != missingValueFloat &&
        airTemperature[jobs] != missingValueFloat &&
        relativeHumidity[jobs] != missingValueFloat) {
      float esat_pure_water = formulas::SatVaporPres_fromTemp(
          airTemperature[jobs], SatVaporPres_fromTemp_form);
      float esat_air = formulas::SatVaporPres_correction(
          esat_pure_water, airTemperature[jobs], pressure[jobs],
          formulas::Formulation::Gill);
      float qsat = formulas::Qsat_From_Psat(esat_air, pressure[jobs],
                                            formulas::Formulation::GillUKMO);
      specificHumidity[jobs] = qsat * relativeHumidity[jobs] / relativeHumidityAtSaturation();
    }
  }
  putObservation(specifichumidityvariable_, specificHumidity);
}

/************************************************************************************/

void Cal_SpecificHumidity::methodDEFAULT(
    const std::vector<bool> &apply,
    formulas::Formulation SatVaporPres_fromTemp_form) {

  const size_t nlocs = obsdb_.nlocs();
  float esat, qvs, qv, satVaporPres;
  std::vector<float> relativeHumidity;
  std::vector<float> dewPointTemperature;
  std::vector<float> airTemperature;
  std::vector<float> pressure;
  std::vector<float> specificHumidity(nlocs);
  bool have_dewpoint = false;

  if (obsdb_.has("ObsValue", relativehumidityvariable_)) {
    getObservation("ObsValue", relativehumidityvariable_,
                   relativeHumidity, true);
  } else {
    oops::Log::debug() << "Looks like relative humidity is not present, "
                       << "so looking for dewpoint instead: "
                       << dewpointtemperaturevariable_ << std::endl;
    getObservation("ObsValue", dewpointtemperaturevariable_,
                   dewPointTemperature, true);
    if (dewPointTemperature.empty()) {
      oops::Log::warning() << "Neither relative humidity or dewpoint temperature exists. "
                           << "Must have one or the other." << std::endl;
      throw eckit::BadValue("Must have either relative humidity or dewpoint temperature"
                            " to calculate specific humidity", Here());
    }
    have_dewpoint = true;
  }
  getObservation("ObsValue", temperaturevariable_,
                 airTemperature, true);
  getObservation(pressuregroupvariable_, pressurevariable_,
                 pressure, false);
  if (pressure.empty()) {
    getObservation("ObsValue", pressureat2mvariable_,
                   pressure, true);
  }

  if (have_dewpoint) {
    if (!oops::allVectorsSameSize(dewPointTemperature, pressure)) {
      oops::Log::warning() << "Vector sizes: "
                           << oops::listOfVectorSizes(dewPointTemperature, pressure) << std::endl;
      throw eckit::BadValue("At least one variable vector is the wrong size or empty out of "
                            "dew point temperature and pressure", Here());
    }
  } else {
    if (!oops::allVectorsSameSize(relativeHumidity, airTemperature, pressure)) {
      oops::Log::warning() << "Vector sizes: "
                           << oops::listOfVectorSizes(relativeHumidity, airTemperature,
                                                      pressure)
                           << std::endl;
      throw eckit::BadValue("At least one variable vector is the wrong size or empty out of "
                            "relative humidity, air temperature and pressure", Here());
    }
  }

  // Initialise this vector with missing value
  specificHumidity.assign(nlocs, missingValueFloat);

  // From dewpoint and pressure, we get vapor pressure that easily converts to mixing ratio
  // then to final answer of specific humidity.  Otherwise, with RH, we must have temp
  // to compute saturated vapor pressure, convert this to mixing ratio using relative
  // humidity, then end up with final conversion of mixing ratio to specific humdity.

  if (have_dewpoint) {
    for (size_t jobs = 0; jobs < nlocs; ++jobs) {
      // if the data have been excluded by the where statement
      if (!apply[jobs]) continue;
      if (pressure[jobs] != missingValueFloat && dewPointTemperature[jobs] != missingValueFloat) {
        satVaporPres = formulas::SatVaporPres_fromTemp(
            dewPointTemperature[jobs], SatVaporPres_fromTemp_form);
        esat = std::min(pressure[jobs]*0.15f, satVaporPres);
        qv = Constants::rd_over_rv * esat/(pressure[jobs]-esat);
        specificHumidity[jobs] = std::max(1.0e-12f, qv/(1.0f+qv));
      }
    }
  } else {
    for (size_t jobs = 0; jobs < nlocs; ++jobs) {
      // if the data have been excluded by the where statement
      if (!apply[jobs]) continue;
      if (pressure[jobs] != missingValueFloat && airTemperature[jobs] != missingValueFloat &&
                relativeHumidity[jobs] != missingValueFloat) {
        satVaporPres = formulas::SatVaporPres_fromTemp(
            airTemperature[jobs], SatVaporPres_fromTemp_form);
        esat = std::min(pressure[jobs]*0.15f, satVaporPres);
        qvs = Constants::rd_over_rv * esat/(pressure[jobs]-esat);
        qv = std::max(1.0e-12f, relativeHumidity[jobs] * qvs /
                                    relativeHumidityAtSaturation());
        specificHumidity[jobs] = std::max(1.0e-12f, qv/(1.0f+qv));
      }
    }
  }
  putObservation(specifichumidityvariable_, specificHumidity);
}

/**************************************************************************************************/
//  Cal_VirtualTemperature
/**************************************************************************************************/

static TransformMaker<Cal_VirtualTemperature>
    makerCal_VirtualTemperature_("VirtualTemperature");

Cal_VirtualTemperature::Cal_VirtualTemperature(
    const Parameters_ &options,
    const ObsFilterData &data,
    ioda::ObsDataVector<int> &flags,
    ioda::ObsDataVector<float> &obserr)
  : TransformBase(options, data, flags, obserr),
      specifichumidityvariable_(options.SpecificHumidityVariable),
      temperaturevariable_(options.TemperatureVariable),
      virtualtempvariable_(options.VirtualTempVariable)
{}

/**************************************************************************************************/

void Cal_VirtualTemperature::runTransform(const std::vector<bool> &apply) {
  oops::Log::trace() << " --> Retrieve Virtual Temperature"
            << std::endl;
  oops::Log::trace() << "      --> method: " << options_.Method.value() << std::endl;
  oops::Log::trace() << "      --> obsName: " << obsName() << std::endl;

  // Get the right method
  switch (method()) {
    default: {
      methodDEFAULT(apply);
      break;
    }
  }
}

/**************************************************************************************************/

void Cal_VirtualTemperature::methodDEFAULT(const std::vector<bool> &apply) {
  const size_t nlocs = obsdb_.nlocs();

  float qv;

  std::vector<float> specificHumidity;
  std::vector<float> airTemperature;
  std::vector<float> virtualTemperature(nlocs);

  getObservation("ObsValue", specifichumidityvariable_, specificHumidity, true);
  getObservation("ObsValue", temperaturevariable_, airTemperature, true);

  if (!oops::allVectorsSameSize(specificHumidity, airTemperature)) {
    oops::Log::warning() << "Vector sizes: "
                         << oops::listOfVectorSizes(specificHumidity, airTemperature)
                         << std::endl;
    throw eckit::BadValue("At least one variable vector is the wrong size or empty out of "
                          "specific humidity and air temperature", Here());
  }

  // Initialise this vector with missing value
  virtualTemperature.assign(nlocs, missingValueFloat);

  // Loop over all obs
  for (size_t jobs = 0; jobs < nlocs; ++jobs) {
    // if the data have been excluded by the where statement
    if (!apply[jobs]) continue;

    if (specificHumidity[jobs] != missingValueFloat && airTemperature[jobs] != missingValueFloat) {
      // Convert specific humidity to water vapor mixing ratio
      qv = std::max(1.0e-12f, specificHumidity[jobs]/(1.0f-specificHumidity[jobs]));
      virtualTemperature[jobs] = airTemperature[jobs]*(1.0f + 0.61f*qv);
    }
  }
  putObservation(virtualtempvariable_, virtualTemperature);
}

/**************************************************************************************************/

}  // namespace ufo
