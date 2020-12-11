/*
 * (C) Copyright 2020 Met Office UK
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/filters/variabletransforms/WindSpeedAndDirection.h"

#include <cmath>
#include <vector>

#include "eckit/config/Configuration.h"

#include "ioda/ObsDataVector.h"
#include "ioda/ObsSpace.h"

#include "oops/util/Logger.h"

#include "ufo/filters/Variable.h"
#include "ufo/utils/Constants.h"

namespace ufo {

// -----------------------------------------------------------------------------

WindSpeedAndDirection::WindSpeedAndDirection(ioda::ObsSpace & obsdb,
                                             const eckit::Configuration & config,
                                             std::shared_ptr<ioda::ObsDataVector<int> > flags,
                                             std::shared_ptr<ioda::ObsDataVector<float> > obserr)
  : FilterBase(obsdb, config, flags, obserr)

{
  oops::Log::trace() << "WindSpeedAndDirection constructor starting" << std::endl;
}

// -----------------------------------------------------------------------------

WindSpeedAndDirection::~WindSpeedAndDirection() {
  oops::Log::trace() << "WindSpeedAndDirection destructed" << std::endl;
}

// -----------------------------------------------------------------------------

void WindSpeedAndDirection::applyFilter(const std::vector<bool> & apply,
                                        const Variables & filtervars,
                                        std::vector<std::vector<bool>> & flagged) const {
  oops::Log::trace() << "WindSpeedAndDirection applyFilter" << std::endl;

  const float missing = util::missingValue(missing);
  const float deg = static_cast<float>(Constants::rad2deg);
  const size_t nlocs = obsdb_.nlocs();

  std::vector<float> windSpeed(nlocs), windFromDirection(nlocs);
  std::vector<float> u(nlocs), v(nlocs);

  obsdb_.get_db("ObsValue", "eastward_wind", u);
  obsdb_.get_db("ObsValue", "northward_wind", v);
  // wind vector is missing unless a valid value calculated below
  windSpeed.assign(nlocs, missing);
  windFromDirection.assign(nlocs, missing);

// Loop over all obs
  for (size_t jobs = 0; jobs < nlocs; ++jobs) {
    if (apply[jobs]) {
      // Calculate wind vector
      if (u[jobs] != missing && v[jobs] != missing) {
        windSpeed[jobs] = hypot(u[jobs], v[jobs]);
        if (u[jobs] == 0 && v[jobs] ==0) {
            windFromDirection[jobs] = 0;
        } else {
            windFromDirection[jobs] = fmod((270.0 - atan2(v[jobs], u[jobs]) * deg), 360.0);
        }
        oops::Log::debug() << "eastward_wind, northward_wind:"
                           << u[jobs] << ", " << v[jobs]
                           << " wind_speed=" << windSpeed[jobs]
                           << " wind_from_direction=" << windFromDirection[jobs]
                           << std::endl;
      }
    }
  }

  // put new variable at existing locations
  obsdb_.put_db("ObsValue", "wind_speed", windSpeed);
  obsdb_.put_db("ObsValue", "wind_from_direction", windFromDirection);
}

// -----------------------------------------------------------------------------

void WindSpeedAndDirection::print(std::ostream & os) const {
  os << "WindSpeedAndDirection";
}

// -----------------------------------------------------------------------------

}  // namespace ufo
