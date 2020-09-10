/*
 * (C) Copyright 2020 Met Office UK
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/filters/variabletransforms/WindVector.h"

#include <cmath>
#include <vector>

#include "eckit/config/Configuration.h"

#include "ioda/ObsDataVector.h"
#include "ioda/ObsSpace.h"

#include "oops/util/Logger.h"

namespace ufo {

// -----------------------------------------------------------------------------

WindVector::WindVector(ioda::ObsSpace & obsdb, const eckit::Configuration & config,
                                 boost::shared_ptr<ioda::ObsDataVector<int> > flags,
                                 boost::shared_ptr<ioda::ObsDataVector<float> > obserr)
  : FilterBase(obsdb, config, flags, obserr)

{
  oops::Log::trace() << "WindVector contructor starting" << std::endl;
}

// -----------------------------------------------------------------------------

WindVector::~WindVector() {
  oops::Log::trace() << "WindVector destructed" << std::endl;
}

// -----------------------------------------------------------------------------

// currently none of the filter arguments are used
void WindVector::applyFilter(const std::vector<bool> & apply,
                                  const Variables & filtervars,
                                  std::vector<std::vector<bool>> & flagged) const {
  oops::Log::trace() << "WindVector priorFilter" << std::endl;

  const float missing = util::missingValue(missing);
  const float deg = static_cast<float>(Constants::rad2deg);
  const size_t nlocs = obsdb_.nlocs();

  std::vector<float> Zfff(nlocs), Zddd(nlocs);
  std::vector<float> u(nlocs), v(nlocs);

  obsdb_.get_db("ObsValue", "eastward_wind", u);
  obsdb_.get_db("ObsValue", "northward_wind", v);
  // wind vector is missing unless a valid value calculated below
  Zfff.assign(nlocs, missing);
  Zddd.assign(nlocs, missing);

// Loop over all obs
  for (size_t jobs = 0; jobs < nlocs; ++jobs) {
    if (apply[jobs]) {
      oops::Log::trace() << "WindVector priorFilter -jobs = " << jobs << std::endl;
      // Calculate wind vector
      if (u[jobs] != missing && v[jobs] != missing) {
        Zfff[jobs] = sqrt(pow(u[jobs], 2) + pow(v[jobs], 2));
        if (u[jobs] == 0 && v[jobs] ==0) {
            Zddd[jobs] = 0;
        }
        else {
            Zddd[jobs] = fmod((270.0 - atan2 (v[jobs], u[jobs]) * deg), 360.0);
        }
        oops::Log::debug() << "eastward_wind, northward_wind:" << u[jobs] << ", " << v[jobs]
                           << " wind_speed=" << Zfff[jobs] << " wind_from_direction=" << Zddd[jobs]
                           << std::endl;
      }
    }
  }

  // put new variable at existing locations
  obsdb_.put_db("ObsValue", "wind_speed", Zfff);
  obsdb_.put_db("ObsValue", "wind_from_direction", Zddd);
  }

// -----------------------------------------------------------------------------

void WindVector::print(std::ostream & os) const {
  os << "WindVector::print not yet implemented ";
}

// -----------------------------------------------------------------------------

}  // namespace ufo
