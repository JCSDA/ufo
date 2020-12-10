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

// ref is denominator, val is numerator ---------------------------------

RatioCheck::RatioCheck(ioda::ObsSpace & obsdb, const eckit::Configuration & config,
                                 std::shared_ptr<ioda::ObsDataVector<int> > flags,
                                 std::shared_ptr<ioda::ObsDataVector<float> > obserr)
  : FilterBase(obsdb, config, flags, obserr),
    den_(config_.getString("denominator")), num_(config_.getString("numerator"))
{
  oops::Log::trace() << "RatioCheck contructor starting" << std::endl;
  allvars_ += den_;
  allvars_ += num_;
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


// Get denominator values and numerators to compare (as floats)
  std::vector<float> den, num;
  data_.get(den_, den);
  data_.get(num_, num);
  ASSERT(den.size() == num.size());

// Loop over all obs
  for (size_t jobs = 0; jobs < nlocs; ++jobs) {
    if (apply[jobs]) {
      // check to see if one of the denominator or numerator is missing
      if (num[jobs] == missing || den[jobs] == missing || den[jobs] == 0.0) {
        for (size_t jv = 0; jv < filtervars.nvars(); ++jv) {
          flagged[jv][jobs] = true;
        }
      } else {
// Check if ratio is within min/max value range and set
        float ratio = num[jobs] / den[jobs];
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
