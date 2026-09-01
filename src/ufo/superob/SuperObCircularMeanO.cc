/*
 * (C) Crown Copyright 2026 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/superob/SuperObCircularMeanO.h"

#include <optional>

#include <cassert>
#include <cmath>
#include <vector>

#include "oops/util/missingValues.h"

namespace ufo {

static SuperObMaker<SuperObCircularMeanO> makerCircularMeanO_(
    "circular mean obs");

SuperObCircularMeanO::SuperObCircularMeanO(
    const SuperObCircularMeanOParameters& params, const ObsFilterData& data,
    const std::vector<bool>& apply, const Variables& filtervars,
    const ioda::ObsDataVector<int>& flags,
    std::vector<std::vector<bool>>& flagged)
    : SuperObBase(params, data, apply, filtervars, flags, flagged),
      params_(params) {
  oops::Log::trace() << "SuperObCircularMeanO constructor" << std::endl;
}

bool SuperObCircularMeanO::computeSuperOb(const std::vector<std::size_t>& locs,
                                          const std::vector<float>& obs,
                                          const std::vector<float>& /*hofx*/,
                                          const ioda::ObsDataRow<int>& flags,
                                          std::vector<float>& superobs,
                                          std::vector<bool>& flagged) const {
  const float missing = util::missingValue<float>();

  // Get bounds for circular mean
  const float lowerBound = params_.lowerBound.value();
  const float exclusiveUpperBound = params_.exclusiveUpperBound.value();
  const float range = exclusiveUpperBound - lowerBound;
  if (range <= 0.0f) {
    throw eckit::BadParameter(
        "\"exclusive upper bound\" must be greater than \"lower bound\" for circular mean obs",
        Here());
  }
  const bool assignToAllValues = params_.assignToAllValues.value();

  // Deduplicate locations if grouping values are provided
  const std::vector<std::size_t>& locsForComputation =
      params_.groupingVariable.value()
          ? getUniqueLocations(locs, *params_.groupingVariable.value(), obs,
                               flags)
          : locs;

  bool superObComputed = false;
  if (locsForComputation.empty()) return superObComputed;

  std::vector<float> validValues;
  // Set on first valid contributing location; optional avoids implying a
  // default location.
  std::optional<size_t> superobloc;

  for (std::size_t jloc : locsForComputation) {
    const float obsValue = obs[jloc];
    // Only consider locations which have valid observation values
    // and are either passing QC or have been marked as passive (i.e. H(x) is
    // computed but the observation is not assimilated).
    if ((flags[jloc] != QCflags::pass && flags[jloc] != QCflags::passive) ||
        obsValue == missing) {
      continue;
    }
    if (!superObComputed) {
      superobloc = jloc;
      superObComputed = true;
    }
    validValues.push_back(obsValue);
  }

  if (validValues.empty()) {
    assert(!superObComputed);
    return false;
  }

  // Calculate circular mean
  double sumSin = 0.0;
  double sumCos = 0.0;

  // Check if normalization is needed (skip if bounds are 0 and 2π)
  const bool needsNormalization =
      !(std::abs(lowerBound) < 1e-5f &&
        std::abs(exclusiveUpperBound - static_cast<float>(2.0 * M_PI)) < 1e-5f);

  if (needsNormalization) {
    // Convert values to [0, 2π) range relative to [lowerBound,
    // exclusiveUpperBound)
    for (float value : validValues) {
      // Normalize value to [0, 1) range
      const double valueDouble = static_cast<double>(value);
      const double normalizedValue = (valueDouble - lowerBound) / range;
      const double radians = normalizedValue * 2.0 * M_PI;
      // Get unit circle components
      sumSin += std::sin(radians);
      sumCos += std::cos(radians);
    }
  } else {
    // Values are already in [0, 2π) range
    for (float value : validValues) {
      // Get unit circle components
      const double valueDouble = static_cast<double>(value);
      sumSin += std::sin(valueDouble);
      sumCos += std::cos(valueDouble);
    }
  }

  // Angle is in range [-π, π], needs conversion back to original range
  double meanAngleRadians = std::atan2(sumSin, sumCos);

  float circularMean;
  if (needsNormalization) {
    // Convert back to [lowerBound, exclusiveUpperBound) range - first normalize
    // to [-0.5, 0.5]
    double normalizedMean = meanAngleRadians / (2.0 * M_PI);
    // Shift negative values to [0, 1) range where 0 corresponds to lowerBound
    // and 1 to exclusiveUpperBound
    if (normalizedMean < 0.0) {
      normalizedMean += 1.0;
    }
    // Scale to output
    circularMean = static_cast<float>(lowerBound + normalizedMean * range);
  } else {
    // Shift negative values to [0, 2π) range
    if (meanAngleRadians < 0.0) {
      meanAngleRadians += 2.0 * M_PI;
    }
    circularMean = static_cast<float>(meanAngleRadians);
  }

  // Assign superob value to locations
  assert(superobloc.has_value());
  assignSuperObToLocations(locs, circularMean, *superobloc, assignToAllValues,
                           superobs, flagged);

  return true;
}

}  // namespace ufo
