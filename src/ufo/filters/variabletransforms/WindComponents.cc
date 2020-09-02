/*
 * (C) Copyright 2020 Met Office UK
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/filters/variabletransforms/WindComponents.h"

#include <cmath>
#include <vector>

#include "eckit/config/Configuration.h"

#include "ioda/ObsDataVector.h"
#include "ioda/ObsSpace.h"

#include "oops/util/Logger.h"

namespace ufo {

// -----------------------------------------------------------------------------

WindComponents::WindComponents(ioda::ObsSpace & obsdb, const eckit::Configuration & config,
                                 boost::shared_ptr<ioda::ObsDataVector<int> > flags,
                                 boost::shared_ptr<ioda::ObsDataVector<float> > obserr)
  : FilterBase(obsdb, config, flags, obserr)

{
  oops::Log::trace() << "WindComponents contructor starting" << std::endl;
}

// -----------------------------------------------------------------------------

WindComponents::~WindComponents() {
  oops::Log::trace() << "WindComponents destructed" << std::endl;
}

// -----------------------------------------------------------------------------

// currently none of the filter arguments are used
void WindComponents::applyFilter(const std::vector<bool> & apply,
                                  const Variables & filtervars,
                                  std::vector<std::vector<bool>> & flagged) const {
  oops::Log::trace() << "WindComponents priorFilter" << std::endl;

  const float missing = util::missingValue(missing);
  const size_t nlocs = obsdb_.nlocs();

  std::vector<float> u(nlocs), v(nlocs);
  std::vector<float> Zfff(nlocs), Zddd(nlocs);

  obsdb_.get_db("ObsValue", "wind_speed", Zfff);
  obsdb_.get_db("ObsValue", "wind_from_direction", Zddd);
  // wind components are missing unless a valid value calculated below
  u.assign(nlocs, missing);
  v.assign(nlocs, missing);

// Loop over all obs
  for (size_t jobs = 0; jobs < nlocs; ++jobs) {
      // Have removed 'apply' here: 'where' not applicable for a unit transformation
      // Calculate wind components
      // OPS checks for extremevalueflag (speed < 0, direction > 360 or < 0).
      if (Zddd[jobs] != missing && Zfff[jobs] != missing) {
        u[jobs] = -Zfff[jobs] * sin(Zddd[jobs] * static_cast<float>(Constants::deg2rad));
        v[jobs] = -Zfff[jobs] * cos(Zddd[jobs] * static_cast<float>(Constants::deg2rad));
        oops::Log::debug() << "wind_speed, wind_from_direction: " << Zfff[jobs] << ", "
                           << Zddd[jobs] << ", eastward_wind=" << u[jobs] << ", northward_wind=" << v[jobs]
                           << std::endl;
      }
  }
  // put new variable at existing locations
  obsdb_.put_db("ObsValue", "eastward_wind", u);
  obsdb_.put_db("ObsValue", "northward_wind", v);
  }

// -----------------------------------------------------------------------------

void WindComponents::print(std::ostream & os) const {
  os << "WindComponents::print not yet implemented ";
}

// -----------------------------------------------------------------------------

}  // namespace ufo
