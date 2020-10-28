/*
 * (C) Copyright 2020 Met Office UK
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_FILTERS_METOFFICEBUDDYCOLLECTORV2_H_
#define UFO_FILTERS_METOFFICEBUDDYCOLLECTORV2_H_

#include "MetOfficeBuddyCollector.h"

#include <vector>

namespace ufo {

/// \brief An implementation of the MetOfficeBuddyCollector interface correcting the deficiencies
/// of MetOfficeBuddyCollectorV1.
class MetOfficeBuddyCollectorV2 : public MetOfficeBuddyCollector {
 public:
  MetOfficeBuddyCollectorV2(const MetOfficeBuddyCheckParameters &options,
                            const std::vector<float> &latitudes,
                            const std::vector<float> &longitudes,
                            const std::vector<int> &stationIds);

  void examinePotentialBuddy(int obsIdB) override;

  void appendBuddyPairsTo(std::vector<MetOfficeBuddyPair> &buddyPairs) const override;

  void reset(int obsIdA) override;

 private:
  std::vector<MetOfficeBuddyPair> buddyPairs_;
};

}  // namespace ufo

#endif  // UFO_FILTERS_METOFFICEBUDDYCOLLECTORV2_H_
