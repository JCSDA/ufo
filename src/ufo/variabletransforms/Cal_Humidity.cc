/*
 * (C) Crown copyright 2020, Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/utils/Constants.h"
#include "ufo/variabletransforms/Cal_Humidity.h"

namespace ufo {
/**************************************************************************************************/
//  Cal_RelativeHumidity
/**************************************************************************************************/

static TransformMaker<Cal_RelativeHumidity>
    makerCal_RelativeHumidity_("RelativeHumidity");

Cal_RelativeHumidity::Cal_RelativeHumidity(
    const VariableTransformsParameters &options, ioda::ObsSpace &os,
    const std::shared_ptr<ioda::ObsDataVector<int>> &flags,
    const std::vector<bool> &apply)
    : TransformBase(options, os, flags, apply) {}

/**************************************************************************************************/

void Cal_RelativeHumidity::runTransform() {
  oops::Log::trace() << " --> Retrieve Relative humidity"
            << std::endl;
  oops::Log::trace() << "      --> method: " << method() << std::endl;
  oops::Log::trace() << "      --> formulation: " << formulation() << std::endl;
  oops::Log::trace() << "      --> obsName: " << obsName() << std::endl;

  // Get the right method
  switch (method()) {
    case formulas::MethodFormulation::UKMO: {
      methodUKMO();
      break;
    }
    case formulas::MethodFormulation::NCAR:
    case formulas::MethodFormulation::NOAA:
    default: {
      methodDEFAULT();
      break;
    }
  }
}

/**************************************************************************************************/
/*
Calculates relative humidity (RH_ice) from dew point temperature,
or converts RH_water to RH_ice

Method: -
  Saturated specific humidity at the dew point (ie w.r.t water), and saturated
  specific humidity (w.r.t ice below 0C) at the air temperature are calculated.
  Relative humidity is then calculated using :
  
         RH = (QSAT(DEW POINT)/QSAT(DRY BULB))*100
  
   For some temperatures (e.g. when dew point = temperature), supersaturation
   w.r.t ice may occur.
   The option AllowSuperSaturation (false by default) controls whether upper air relative humidity
   is capped at 100%.
   If the pressure, dew point or temperature take extreme
   values or are missing, the relative humidity is set to missing data.
*/
void Cal_RelativeHumidity::methodUKMO() {
  const size_t nlocs_ = obsdb_.nlocs();
  // return if no data
  if (obsdb_.nlocs() == 0) {
    return;
  }

  float Q_sub_s_ice, Q_sub_s_w;
  float pressure, temperature, dewPoint;
  std::vector<float> airTemperature;
  std::vector<float> dewPointTemperature;
  std::vector<float> relativeHumidity;
  std::vector<float> airPressure;
  bool surfaceData = true;
  bool hasBeenUpdated = false;

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
  if (obsdb_.has("ObsValue", "pressure_surface") &&
      obsdb_.has("ObsValue", "air_temperature_surface") &&
      obsdb_.has("ObsValue", "dew_point_temperature_surface")) {
    getObservation("ObsValue", "pressure_surface",
                   airPressure, true);
    getObservation("ObsValue", "air_temperature_surface",
                   airTemperature, true);
    getObservation("ObsValue", "dew_point_temperature_surface",
                   dewPointTemperature, true);
    getObservation("ObsValue", "relative_humidity_surface",
                   relativeHumidity);
  } else {
    getObservation("ObsValue", "air_pressure",
                   airPressure, true);
    getObservation("ObsValue", "air_temperature",
                   airTemperature, true);
    getObservation("ObsValue", "dew_point_temperature",
                   dewPointTemperature, true);
    getObservation("ObsValue", "relative_humidity",
                   relativeHumidity);
    surfaceData = false;
  }
  if (relativeHumidity.empty()) {
    relativeHumidity.assign(nlocs_, missingValueFloat);
  }

  // 2. making sure we have what we need is here
  // -----------------------------------------------------------------------------------------------
  if (!oops::allVectorsSameNonZeroSize(airPressure, airTemperature, dewPointTemperature)) {
    oops::Log::warning() << "Vector sizes: "
                         << oops::listOfVectorSizes(airPressure, airTemperature,
                                                    dewPointTemperature)
                         << std::endl;
    throw eckit::BadValue("At least one vector is the wrong size or empty out of "
                          "P, T and Td", Here());
  }

  // Lambda function to evaluate Saturated specific humidity
  // temp_1: airTemprature or dewPointTemperature
  // temp_2: airTemperature
  // -----------------------------------------------------------------------------------------------
  auto evaluateSatSpecHumidity = [&](float temp_1, float temp_2){
    float e_sub_s_w, e_sub_s_ice;
    // sat. vapor pressure from Dewpoint  temperature - wrt water
    e_sub_s_w = formulas::SatVaporPres_fromTemp(temp_1,
                                                formulation());
    e_sub_s_w = formulas::SatVaporPres_correction(e_sub_s_w,
                                                  temp_1,
                                                  pressure,
                                                  formulation());
    // Convert sat. vapor pressure (wrt water) to saturated specific humidity (water)
    Q_sub_s_w = formulas::Qsat_From_Psat(e_sub_s_w, pressure);

    // sat. vapor pressure from Drybulb temperature  - wrt ice
    e_sub_s_ice = formulas::SatVaporPres_fromTemp(temp_2,
                                                 formulas::MethodFormulation::LandoltBornstein);
    e_sub_s_ice = formulas::SatVaporPres_correction(e_sub_s_ice,
                                                    temp_2,
                                                    pressure,
                                                    formulation());
    // Convert sat. vapor pressure (wrt ice) to saturated specific humidity (ice)
    Q_sub_s_ice = formulas::Qsat_From_Psat(e_sub_s_ice, pressure);
  };

  // 3. Loop over each record
  // -----------------------------------------------------------------------------------------------
  for (ioda::ObsSpace::RecIdxIter irec = obsdb_.recidx_begin();
       irec != obsdb_.recidx_end(); ++irec) {
    const std::vector<std::size_t> &rSort = obsdb_.recidx_vector(irec);

    // 3.1 Loop over each record
    for (size_t iloc : rSort) {
      // if the data have been excluded by the where statement
      if (!apply_[iloc]) continue;

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

      // Cycle if observatoin are valid
      if (temperature == missingValueFloat ||
          pressure <= 1.0) continue;

      // if dewpoint temperature is reported (most stations)
      if (dewPoint > 1.0) {
        // calculate saturated specific humidity wrt water and ice
        evaluateSatSpecHumidity(dewPoint, temperature);

        // if saturated specific humidity wrt water and ice are positive
        // calculate relative humidity
        if (Q_sub_s_w > 0 && Q_sub_s_ice > 0) {
          relativeHumidity[iloc] = (Q_sub_s_w / Q_sub_s_ice) * 100.0;
          if (!AllowSuperSaturation())
            relativeHumidity[iloc] = std::min(100.0f, relativeHumidity[iloc]);
          hasBeenUpdated = true;
        }
      // if relative humidity (Rh) is reported (small minority of stations)
      // update from Rh wrt water to Rh wrt ice for temperatures below freezing
      } else if (relativeHumidity[iloc] != missingValueFloat &&
                temperature < ufo::Constants::t0c) {
        // calculate saturated specific humidity wrt water and ice
        evaluateSatSpecHumidity(temperature, temperature);
        // update from Rh wrt water to Rh wrt ice for temperatures below freezing
        relativeHumidity[iloc] *= (Q_sub_s_w / Q_sub_s_ice);
        if (!AllowSuperSaturation())
          relativeHumidity[iloc] = std::min(100.0f, relativeHumidity[iloc]);

        hasBeenUpdated = true;
      }
    }
  }

  // assign the derived relative humidity as DerivedValue
  if (hasBeenUpdated) {
    if (surfaceData) {
      obsdb_.put_db(outputTag, "relative_humidity_surface", relativeHumidity);
    } else {
      obsdb_.put_db(outputTag, "relative_humidity", relativeHumidity);
    }
  }
}

/**************************************************************************************************/

void Cal_RelativeHumidity::methodDEFAULT() {
  const size_t nlocs = obsdb_.nlocs();

  float esat, qvs, qv, satVaporPres;

  std::vector<float> specificHumidity;
  std::vector<float> airTemperature;
  std::vector<float> pressure;
  std::vector<float> relativeHumidity(nlocs);

  getObservation("ObsValue", "specific_humidity",
                 specificHumidity, true);
  getObservation("ObsValue", "air_temperature",
                 airTemperature, true);
  getObservation("MetaData", "air_pressure",
                 pressure, false);
  if (pressure.empty()) {
    getObservation("ObsValue", "surface_pressure",
                   pressure, true);
  }

  if (!oops::allVectorsSameNonZeroSize(specificHumidity, airTemperature, pressure)) {
    oops::Log::warning() << "Vector sizes: "
                         << oops::listOfVectorSizes(specificHumidity, airTemperature,
                                                    pressure)
                         << std::endl;
    throw eckit::BadValue("At least one vector is the wrong size or empty out of "
                          "specific_humidity, air_temperature and pressure", Here());
  }

  // Initialise this vector with missing value
  relativeHumidity.assign(nlocs, missingValueFloat);

  // Loop over all obs
  for (size_t jobs = 0; jobs < nlocs; ++jobs) {
    // if the data have been excluded by the where statement
    if (!apply_[jobs]) continue;

    if (specificHumidity[jobs] != missingValueFloat &&
        airTemperature[jobs] != missingValueFloat && pressure[jobs] != missingValueFloat) {
      // Calculate saturation vapor pressure from temperature according to requested formulation
      // Double-check result is always lower than 15% of incoming pressure.
      satVaporPres = formulas::SatVaporPres_fromTemp(airTemperature[jobs], formulation());
      esat = std::min(pressure[jobs]*0.15f, satVaporPres);

      // Convert sat. vapor pressure to sat water vapor mixing ratio
      qvs = 0.622 * esat/(pressure[jobs]-esat);

      // Convert specific humidity to water vapor mixing ratio
      qv = std::max(1.0e-12f, specificHumidity[jobs]/(1.0f-specificHumidity[jobs]));

      // Final RH (which can be greater than 100%) is q/qsat, but set sensible lowest limit
      relativeHumidity[jobs] = std::max(1.0e-6f, qv/qvs);
    }
  }

  obsdb_.put_db("DerivedValue", "relative_humidity", relativeHumidity);
}

/************************************************************************************/
//  Cal_SpecificHumidity
/************************************************************************************/
static TransformMaker<Cal_SpecificHumidity>
    makerCal_SpecificHumidity_("SpecificHumidity");

Cal_SpecificHumidity::Cal_SpecificHumidity(
    const VariableTransformsParameters &options, ioda::ObsSpace &os,
    const std::shared_ptr<ioda::ObsDataVector<int>> &flags,
    const std::vector<bool> &apply)
    : TransformBase(options, os, flags, apply) {}

/************************************************************************************/

void Cal_SpecificHumidity::runTransform() {
  oops::Log::trace() << " Retrieve Specific Humidity" << std::endl;
  oops::Log::trace() << "      --> method: " << method() << std::endl;
  oops::Log::trace() << "      --> formulation: " << formulation() << std::endl;
  oops::Log::trace() << "      --> obsName: " << obsName() << std::endl;

  // Get the right method
  switch (method()) {
    case formulas::MethodFormulation::NCAR:
    case formulas::MethodFormulation::NOAA:
    case formulas::MethodFormulation::UKMO:
    default: {
      methodDEFAULT();
      break;
    }
  }
}
/************************************************************************************/

void Cal_SpecificHumidity::methodDEFAULT() {
  const size_t nlocs = obsdb_.nlocs();
  float esat, qvs, qv, satVaporPres;
  std::vector<float> relativeHumidity;
  std::vector<float> airTemperature;
  std::vector<float> pressure;
  std::vector<float> specificHumidity(nlocs);

  getObservation("ObsValue", "relative_humidity",
                 relativeHumidity, true);
  getObservation("ObsValue", "air_temperature",
                 airTemperature, true);
  getObservation("MetaData", "air_pressure",
                 pressure, false);
  if (pressure.empty()) {
    getObservation("ObsValue", "surface_pressure",
                   pressure, true);
  }

  if (!oops::allVectorsSameNonZeroSize(relativeHumidity, airTemperature, pressure)) {
    oops::Log::warning() << "Vector sizes: "
                         << oops::listOfVectorSizes(relativeHumidity, airTemperature,
                                                    pressure)
                         << std::endl;
    throw eckit::BadValue("At least one vector is the wrong size or empty out of "
                          "relative_humidity, air_temperature and pressure", Here());
  }

  // Initialise this vector with missing value
  specificHumidity.assign(nlocs, missingValueFloat);

  // Loop over all obs
  for (size_t jobs = 0; jobs < nlocs; ++jobs) {
    if (relativeHumidity[jobs] != missingValueFloat &&
        airTemperature[jobs] != missingValueFloat && pressure[jobs] != missingValueFloat) {
      // Calculate saturation vapor pressure from temperature according to requested formulation
      // Double-check result is always lower than 15% of incoming pressure.
      satVaporPres = formulas::SatVaporPres_fromTemp(airTemperature[jobs], formulation());
      esat = std::min(pressure[jobs]*0.15f, satVaporPres);

      // Convert sat. vapor pressure to sat water vapor mixing ratio
      qvs = 0.622 * esat/(pressure[jobs]-esat);

      // Using RH, calculate water vapor mixing ratio
      qv = std::max(1.0e-12f, relativeHumidity[jobs]*qvs);

      // Final RH (which can be greater than 100%) is q/qsat, but set sensible lowest limit
      specificHumidity[jobs] = std::max(1.0e-12f, qv/(1.0f+qv));
    }
  }
  obsdb_.put_db("DerivedValue", "specific_humidity", specificHumidity);
}



}  // namespace ufo

