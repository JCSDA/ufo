/*
 * (C) Crown copyright 2020, Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/filters/ProfileConsistencyCheckParameters.h"
#include "ufo/utils/Constants.h"
#include "ufo/variabletransforms/Cal_Humidity.h"

namespace ufo {
/************************************************************************************/
//  Cal_RelativeHumidity
/************************************************************************************/

static TransformMaker<Cal_RelativeHumidity>
    makerCal_RelativeHumidity_("RelativeHumidity");

Cal_RelativeHumidity::Cal_RelativeHumidity(
    const VariableTransformsParameters &options, ioda::ObsSpace &os,
    const std::shared_ptr<ioda::ObsDataVector<int>> &flags)
    : TransformBase(options, os, flags) {}

/************************************************************************************/

void Cal_RelativeHumidity::runTransform() {
  oops::Log::trace() << " --> Retrieve Relative humidity"
            << std::endl;
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
    const std::shared_ptr<ioda::ObsDataVector<int>> &flags)
    : TransformBase(options, os, flags) {}

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

