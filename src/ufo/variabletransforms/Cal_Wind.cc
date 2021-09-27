/*
 * (C) Crown copyright 2020, Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/variabletransforms/Cal_Wind.h"
#include "ufo/utils/Constants.h"

#include "ufo/filters/VariableTransformParametersBase.h"


namespace ufo {

/************************************************************************************/
//  Cal_WindSpeedAndDirection
/************************************************************************************/

static TransformMaker<Cal_WindSpeedAndDirection>
    makerCal_WindSpeedAndDirection_("WindSpeedAndDirection");

Cal_WindSpeedAndDirection::Cal_WindSpeedAndDirection(
    const GenericVariableTransformParameters &options, const ObsFilterData &data,
    const std::shared_ptr<ioda::ObsDataVector<int>> &flags)
    : TransformBase(options, data, flags) {
}

/************************************************************************************/

void Cal_WindSpeedAndDirection::runTransform(const std::vector<bool> &apply) {
  oops::Log::trace() << " --> Retrieve wind speed and direction"
            << std::endl;
  oops::Log::trace() << "      --> method: " << method() << std::endl;
  oops::Log::trace() << "      --> obsName: " << obsName() << std::endl;

  const size_t nlocs = obsdb_.nlocs();

  std::vector<float> u, v;
  getObservation("ObsValue", "eastward_wind",
                 u, true);
  getObservation("ObsValue", "northward_wind",
                 v, true);

  if (!oops::allVectorsSameNonZeroSize(u, v)) {
    oops::Log::warning() << "Vector sizes: "
                         << oops::listOfVectorSizes(u, v)
                         << std::endl;
    throw eckit::BadValue("At least one vector is the wrong size or empty out of "
                          "U, and V", Here());
  }

  std::vector<float> windSpeed(nlocs), windFromDirection(nlocs);
  windSpeed.assign(nlocs, missingValueFloat);
  windFromDirection.assign(nlocs, missingValueFloat);

  // Loop over all obs
  for (size_t jobs = 0; jobs < nlocs; ++jobs) {
    // if the data have been excluded by the where statement
    if (!apply[jobs]) continue;

    // Calculate wind vector
    if (u[jobs] != missingValueFloat && v[jobs] != missingValueFloat) {
       windFromDirection[jobs] = formulas::GetWindDirection(u[jobs], v[jobs]);
       windSpeed[jobs] = formulas::GetWindSpeed(u[jobs], v[jobs]);
    }
  }
  // put new variable at existing locations
  putObservation("wind_speed", windSpeed);
  putObservation("wind_from_direction", windFromDirection);
}

/************************************************************************************/
//  Cal_WindComponents
/************************************************************************************/
static TransformMaker<Cal_WindComponents>
    makerCal_WindComponents_("WindComponents");

Cal_WindComponents::Cal_WindComponents(
    const GenericVariableTransformParameters &options,
    const ObsFilterData &data,
    const std::shared_ptr<ioda::ObsDataVector<int>> &flags)
    : TransformBase(options, data, flags) {}

/************************************************************************************/

void Cal_WindComponents::runTransform(const std::vector<bool> &apply) {
  oops::Log::trace() << " --> Retrieve wind component"
            << std::endl;
  oops::Log::trace() << "      --> method: " << method() << std::endl;
  oops::Log::trace() << "      --> obsName: " << obsName() << std::endl;

  const size_t nlocs = obsdb_.nlocs();

  std::vector<float> windSpeed, windFromDirection;
  getObservation("ObsValue", "wind_speed",
                 windSpeed, true);
  getObservation("ObsValue", "wind_from_direction",
                 windFromDirection, true);

  if (!oops::allVectorsSameNonZeroSize(windSpeed, windFromDirection)) {
    oops::Log::warning() << "Vector sizes: "
                         << oops::listOfVectorSizes(windSpeed, windFromDirection)
                         << std::endl;
    throw eckit::BadValue("At least one vector is the wrong size or empty out of "
                          "wind speed, and direction", Here());
  }

  std::vector<float> u(nlocs), v(nlocs);
  u.assign(nlocs, missingValueFloat);
  v.assign(nlocs, missingValueFloat);

  // Loop over all obs
  for (size_t jobs = 0; jobs < nlocs; ++jobs) {
    // if the data have been excluded by the where statement
    if (!apply[jobs]) continue;

    // Calculate wind vector
    if (windFromDirection[jobs] != missingValueFloat &&
          windSpeed[jobs] != missingValueFloat && windSpeed[jobs] >= 0) {
      u[jobs] = formulas::GetWind_U(windSpeed[jobs], windFromDirection[jobs]);
      v[jobs] = formulas::GetWind_V(windSpeed[jobs], windFromDirection[jobs]);
    }
  }

  // put new variable at existing locations
  putObservation("eastward_wind", u);
  putObservation("northward_wind", v);
}
}  // namespace ufo

