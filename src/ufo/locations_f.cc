/*
 * (C) Copyright 2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/locations_f.h"

#include <algorithm>
#include <vector>

#include "eckit/exception/Exceptions.h"

namespace ufo {

// -----------------------------------------------------------------------------
std::size_t locations_get_nlocs_f(const Locations & locs) {
  return locs.size();
}
// -----------------------------------------------------------------------------
void locations_get_lons_f(const Locations & locs,
                          const std::size_t & nlocs, double * lons) {
  ASSERT(nlocs == locs.size());
  const std::vector<float> & data = locs.lons();
  std::copy(data.begin(), data.end(), lons);
}
// -----------------------------------------------------------------------------
void locations_get_lats_f(const Locations & locs,
                          const std::size_t & nlocs, double * lats) {
  ASSERT(nlocs == locs.size());
  const std::vector<float> & data = locs.lats();
  std::copy(data.begin(), data.end(), lats);
}
// -----------------------------------------------------------------------------
void locations_get_timemask_f(const Locations & locs,
                              const util::DateTime & t1, const util::DateTime & t2,
                              const std::size_t & nlocs, bool * mask) {
  ASSERT(nlocs == locs.size());
  std::vector<bool> data = locs.isInTimeWindow(t1, t2);
  std::copy(data.begin(), data.end(), mask);
}
// -----------------------------------------------------------------------------
}  // namespace ufo
