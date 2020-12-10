/*
 * (C) Copyright 2017-2018 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/filters/RatioCheck.h"

#include <cmath>
#include <vector>

#include "eckit/config/Configuration.h"

#include "ioda/ObsDataVector.h"
#include "ioda/ObsSpace.h"

#include "oops/util/Logger.h"

namespace ufo {

// -----------------------------------------------------------------------------

RatioCheck::RatioCheck(ioda::ObsSpace & obsdb, const eckit::Configuration & config,
                                 std::shared_ptr<ioda::ObsDataVector<int> > flags,
                                 std::shared_ptr<ioda::ObsDataVector<float> > obserr)
  : FilterBase(obsdb, config, flags, obserr),
    ref_(config_.getString("reference")), val_(config_.getString("value"))
{
  oops::Log::trace() << "RatioCheck contructor starting" << std::endl;
  allvars_ += ref_;
  allvars_ += val_;
}

// -----------------------------------------------------------------------------

RatioCheck::~RatioCheck() {
  oops::Log::trace() << "RatioCheck destructed" << std::endl;
}

// -----------------------------------------------------------------------------

void RatioCheck::applyFilter(const std::vector<bool> & apply,
                                  const Variables & filtervars,
                                  std::vector<std::vector<bool>> & flagged) const {
  oops::Log::trace() << "RatioCheck priorFilter" << std::endl;

  const float missing = util::missingValue(missing);
  const size_t nlocs = obsdb_.nlocs();

// min/max value setup
  float vmin = config_.getFloat("minvalue", missing);
  float vmax = config_.getFloat("maxvalue", missing);

// check if threshold should be absolute or not
  const bool absval = config_.getBool("absolute", false);

// Get reference values and values to compare (as floats)
  std::vector<float> ref, val;
  data_.get(ref_, ref);
  data_.get(val_, val);
  ASSERT(ref.size() == val.size());

// Loop over all obs
  for (size_t jobs = 0; jobs < nlocs; ++jobs) {
    if (apply[jobs]) {
      // check to see if one of the reference or value is missing
      if (val[jobs] == missing || ref[jobs] == missing || ref[jobs] == 0) {
        for (size_t jv = 0; jv < filtervars.nvars(); ++jv) {
          flagged[jv][jobs] = true;
        }
      } else {
// Check if the ratio of val/ref is within min/max range and set flag
        float ratio = val[jobs] / ref[jobs];
        if (absval) {
          ratio = fabs(ratio);
        }
        for (size_t jv = 0; jv < filtervars.nvars(); ++jv) {
          if (vmin != missing && ratio < vmin) flagged[jv][jobs] = true;
          if (vmax != missing && ratio > vmax) flagged[jv][jobs] = true;
        }
      }
    }
  }
}

// -----------------------------------------------------------------------------

void RatioCheck::print(std::ostream & os) const {
  os << "RatioCheck::print not yet implemented ";
}

// -----------------------------------------------------------------------------

}  // namespace ufo
