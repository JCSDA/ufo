/*
 * (C) Crown copyright 2025, Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/superob/SuperObCountO.h"
#include <optional>
#include <cassert>

namespace ufo {

static SuperObMaker<SuperObCount> makerSuperObCount_("count obs");

SuperObCount::SuperObCount(const Parameters_& params,
                           const ObsFilterData& obsdb,
                           const std::vector<bool>& apply,
                           const Variables& filtervars,
                           const ioda::ObsDataVector<int>& flags,
                           std::vector<std::vector<bool>>& flagged)
    : SuperObBase(params, obsdb, apply, filtervars, flags, flagged),
      params_(params) {}

bool SuperObCount::computeSuperOb(const std::vector<std::size_t>& locs,
                                  const std::vector<float>& obs,
                                  const std::vector<float>& /*hofx*/,
                                  const ioda::ObsDataRow<int>& flags,
                                  std::vector<float>& superobs,
                                  std::vector<bool>& flagged) const {
  const float missing = util::missingValue<float>();

  const bool assignToAllValues = params_.assignToAllValues.value();

  // Deduplicate locations if grouping values are provided
  const std::vector<std::size_t>& locsForComputation =
      params_.groupingVariable.value()
          ? getUniqueLocations(locs, *params_.groupingVariable.value(), obs,
                               flags)
          : locs;

  // The only way a superob can fail to be computed is if there are no valid
  // locations to use in the computation
  if (locsForComputation.empty()) return false;

  size_t count = 0;
  // Set on first valid contributing location; optional avoids implying a
  // default location.
  std::optional<size_t> superobloc;
  bool haveSuperObLoc = false;

  for (std::size_t jloc : locsForComputation) {
    const float obsValue = obs[jloc];
    // Only consider locations which have valid observation values
    // and are either passing QC or have been marked as passive (i.e. H(x) is
    // computed but the observation is not assimilated).
    if ((flags[jloc] != QCflags::pass && flags[jloc] != QCflags::passive) ||
        obsValue == missing) {
      continue;
    }
    count++;
    if (!haveSuperObLoc) {
      haveSuperObLoc = true;
      superobloc = jloc;
    }
  }
  if (!haveSuperObLoc) {
    assert(count == 0);
    // We always need to return a value for the superob, so set the superobloc
    // to the first location in locsForComputation, even though it doesn't have
    // a valid observation value.
    superobloc = locsForComputation.front();
  }
  assert(superobloc.has_value());
  assignSuperObToLocations(locs, static_cast<float>(count), *superobloc,
                           assignToAllValues, superobs, flagged);
  return true;
}

}  // namespace ufo
