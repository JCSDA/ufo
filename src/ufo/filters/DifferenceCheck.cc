/*
 * (C) Copyright 2017-2018 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/filters/DifferenceCheck.h"

#include <cmath>
#include <vector>

#include "eckit/config/Configuration.h"

#include "ioda/ObsDataVector.h"
#include "ioda/ObsSpace.h"

#include "oops/util/Logger.h"

namespace ufo {

// -----------------------------------------------------------------------------

DifferenceCheck::DifferenceCheck(ioda::ObsSpace & obsdb, const eckit::Configuration & config,
                                 std::shared_ptr<ioda::ObsDataVector<int> > flags,
                                 std::shared_ptr<ioda::ObsDataVector<float> > obserr)
  : FilterBase(obsdb, config, flags, obserr),
    ref_(config_.getString("reference")), val_(config_.getString("value"))
{
  oops::Log::trace() << "DifferenceCheck contructor starting" << std::endl;
  allvars_ += ref_;
  allvars_ += val_;
}

// -----------------------------------------------------------------------------

DifferenceCheck::~DifferenceCheck() {
  oops::Log::trace() << "DifferenceCheck destructed" << std::endl;
}

// -----------------------------------------------------------------------------

void DifferenceCheck::applyFilter(const std::vector<bool> & apply,
                                  const Variables & filtervars,
                                  std::vector<std::vector<bool>> & flagged) const {
  oops::Log::trace() << "DifferenceCheck priorFilter" << std::endl;

  const float missing = util::missingValue(missing);
  const size_t nlocs = obsdb_.nlocs();

// min/max value setup
  float vmin = config_.getFloat("minvalue", missing);
  float vmax = config_.getFloat("maxvalue", missing);

// check for threshold and if exists, set vmin and vmax appropriately
  const float thresh = config_.getFloat("threshold", missing);
  if (thresh != missing) {
    vmin = -thresh;
    vmax = thresh;
  }

// Get reference values and values to compare (as floats)
  std::vector<float> ref, val;
  data_.get(ref_, ref);
  data_.get(val_, val);
  ASSERT(ref.size() == val.size());

// Loop over all obs
  for (size_t jobs = 0; jobs < nlocs; ++jobs) {
    if (apply[jobs]) {
      // check to see if one of the reference or value is missing
      if (val[jobs] == missing || ref[jobs] == missing) {
        for (size_t jv = 0; jv < filtervars.nvars(); ++jv) {
          flagged[jv][jobs] = true;
        }
      } else {
// Check if difference is within min/max value range and set flag
        float diff = val[jobs] - ref[jobs];
        for (size_t jv = 0; jv < filtervars.nvars(); ++jv) {
          if (vmin != missing && diff < vmin) flagged[jv][jobs] = true;
          if (vmax != missing && diff > vmax) flagged[jv][jobs] = true;
        }
      }
    }
  }
}

// -----------------------------------------------------------------------------

void DifferenceCheck::print(std::ostream & os) const {
  os << "DifferenceCheck::print not yet implemented ";
}

// -----------------------------------------------------------------------------

}  // namespace ufo
