/*
 * (C) Copyright 2020 Met Office UK
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_FILTERS_METOFFICEBUDDYCOLLECTORV1_H_
#define UFO_FILTERS_METOFFICEBUDDYCOLLECTORV1_H_

#include "MetOfficeBuddyCollector.h"

#include <vector>

namespace ufo {

/// \brief An implementation of the MetOfficeBuddyCollector interface intended to produce the
/// same results as Met Office's OPS system.
///
/// A drawback of this implementation is that counters of buddies of certain types (e.g. those from
/// the same zonal band) are sometimes incremented prematurely -- before a candidate buddy has been
/// fully vetted -- and as a result the collector may sometimes fail to collect as many valid
/// as it should.
class MetOfficeBuddyCollectorV1 : public MetOfficeBuddyCollector {
 public:
  MetOfficeBuddyCollectorV1(const MetOfficeBuddyCheckParameters &options,
                            const std::vector<float> &latitudes,
                            const std::vector<float> &longitudes,
                            const std::vector<int> &stationIds);

  void examinePotentialBuddy(int obsIdB) override;

  void appendBuddyPairsTo(std::vector<MetOfficeBuddyPair> &buddyPairs) const override;

  void reset(int obsIdA) override;

 private:
  std::vector<int> potentialBuddies_;
};

}  // namespace ufo

#endif  // UFO_FILTERS_METOFFICEBUDDYCOLLECTORV1_H_
