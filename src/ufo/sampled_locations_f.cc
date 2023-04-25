/*
 * (C) Copyright 2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/sampled_locations_f.h"

#include <algorithm>
#include <vector>

#include "eckit/exception/Exceptions.h"

namespace ufo {

// -----------------------------------------------------------------------------
std::size_t sampled_locations_get_npaths_f(const SampledLocations & locs) {
  return locs.size();
}
// -----------------------------------------------------------------------------
void sampled_locations_get_lons_f(const SampledLocations & locs,
                                  const std::size_t & npaths, double * lons) {
  ASSERT(npaths == locs.size());
  const std::vector<float> & data = locs.lons();
  std::copy(data.begin(), data.end(), lons);
}
// -----------------------------------------------------------------------------
void sampled_locations_get_lats_f(const SampledLocations & locs,
                                  const std::size_t & npaths, double * lats) {
  ASSERT(npaths == locs.size());
  const std::vector<float> & data = locs.lats();
  std::copy(data.begin(), data.end(), lats);
}
// -----------------------------------------------------------------------------
void sampled_locations_get_timemask_f(const SampledLocations & locs,
                                      const util::DateTime & t1,
                                      const util::DateTime & t2,
                                      const std::size_t & npaths, bool * mask) {
  ASSERT(npaths == locs.size());
  std::vector<bool> data = locs.isInTimeWindow(t1, t2);
  std::copy(data.begin(), data.end(), mask);
}
// -----------------------------------------------------------------------------
}  // namespace ufo
