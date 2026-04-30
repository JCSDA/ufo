/*
 * (C) Crown copyright 2026, Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/filters/obsfunctions/CircularDifference.h"

#include <cmath>

#include "ioda/ObsDataVector.h"
#include "oops/util/Logger.h"
#include "oops/util/missingValues.h"
#include "ufo/filters/ObsFilterData.h"

namespace ufo {

// Register this ObsFunction with the factory
static ObsFunctionMaker<CircularDifference>
    makerCircularDifference_("CircularDifference");

// -----------------------------------------------------------------------------

CircularDifference::CircularDifference(const eckit::LocalConfiguration &conf)
    : invars_() {
  oops::Log::trace() << "CircularDifference constructor" << std::endl;
  // Validate and deserialize parameters
  options_.validateAndDeserialize(conf);

  // Add required variables
  invars_ += options_.variableStart.value();
  invars_ += options_.variableEnd.value();

  oops::Log::debug() << "CircularDifference: configured with period="
                     << options_.circularPeriod.value()
                     << ", signed=" << options_.computeSignedDifference.value()
                     << std::endl;
}

// -----------------------------------------------------------------------------

CircularDifference::~CircularDifference() {
  oops::Log::trace() << "CircularDifference destructor" << std::endl;
}

// -----------------------------------------------------------------------------

void CircularDifference::compute(const ObsFilterData &in,
                                 ioda::ObsDataVector<float> &out) const {
  oops::Log::trace() << "CircularDifference compute start" << std::endl;
  const size_t nlocs = in.nlocs();
  const float missing = util::missingValue<float>();
  const float period = options_.circularPeriod.value();
  const bool signedDiff = options_.computeSignedDifference.value();

  // Validate period
  if (period <= 0.0f) {
    throw eckit::BadValue("Circular period must be positive", Here());
  }

  // Get input variable data
  const Variable &varStart = options_.variableStart.value();
  const Variable &varEnd = options_.variableEnd.value();

  ioda::ObsDataVector<float> startData(in.obsspace(),
                                       varStart.toOopsObsVariables());
  in.get(varStart, startData);

  ioda::ObsDataVector<float> endData(in.obsspace(),
                                     varEnd.toOopsObsVariables());
  in.get(varEnd, endData);

  // Compute circular difference for each variable and location
  for (size_t ivar = 0; ivar < out.nvars(); ++ivar) {
    for (size_t iloc = 0; iloc < nlocs; ++iloc) {
      const float start = startData[ivar][iloc];
      const float end = endData[ivar][iloc];

      // Check for missing values
      if (start == missing || end == missing) {
        out[ivar][iloc] = missing;
        continue;
      }

      // Compute raw difference
      float diff = end - start;

      // Normalize difference to [-period/2, period/2)
      // This gives the signed circular difference
      // Uses the formula: ((b - a + half) % period) - half
      const float half = period / 2.0f;
      diff = std::fmod(diff + half, period);
      if (diff < 0.0f) {
        diff += period;  // Ensure non-negative result from fmod (consistent with Python's behavior)
      }
      diff -= half;

      // Apply signed/unsigned option
      if (signedDiff) {
        out[ivar][iloc] = diff;
      } else {
        out[ivar][iloc] = std::abs(diff);
      }
    }
  }
  oops::Log::trace() << "CircularDifference compute complete" << std::endl;
}

// -----------------------------------------------------------------------------

const ufo::Variables &CircularDifference::requiredVariables() const {
  return invars_;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
