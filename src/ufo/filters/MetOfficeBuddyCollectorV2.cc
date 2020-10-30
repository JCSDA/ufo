/*
 * (C) Copyright 2020 Met Office UK
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/filters/MetOfficeBuddyCollectorV2.h"

#include "ufo/filters/MetOfficeBuddyPair.h"

namespace ufo {

MetOfficeBuddyCollectorV2::MetOfficeBuddyCollectorV2(const MetOfficeBuddyCheckParameters &options,
                                                     const std::vector<float> &latitudes,
                                                     const std::vector<float> &longitudes,
                                                     const std::vector<int> &stationIds)
  : MetOfficeBuddyCollector(options, latitudes, longitudes, stationIds)
{}

void MetOfficeBuddyCollectorV2::examinePotentialBuddy(int obsIdB) {
  assert(numBuddiesWithSameStationId_ <= options_.maxNumBuddiesWithSameStationId);

  const bool sameStationId = (stationIds_[obsIdA_] == stationIds_[obsIdB]);
  if (sameStationId && numBuddiesWithSameStationId_ == options_.maxNumBuddiesWithSameStationId)
    return;

  if (std::abs(latitudes_[obsIdA_] - latitudes_[obsIdB]) <= maxLatDifferenceBetweenBuddiesInDeg_) {
    double deltaLatInRad, deltaLonInRad, distanceInKm;
    calcDeltaLatLonAndDistanceTo(obsIdB, deltaLatInRad, deltaLonInRad, distanceInKm);

    if (distanceInKm <= options_.searchRadius) {
      buddyPairs_.push_back(createBuddyPair(obsIdB, deltaLatInRad, deltaLonInRad, distanceInKm));
      ++totalNumBuddies_;
      ++numBuddiesInCurrentBand_;
      if (sameStationId)
        ++numBuddiesWithSameStationId_;
    }
  }
}

void MetOfficeBuddyCollectorV2::appendBuddyPairsTo(
    std::vector<MetOfficeBuddyPair> &buddyPairs) const {
  buddyPairs.insert(buddyPairs.end(), buddyPairs_.begin(), buddyPairs_.end());
}

void MetOfficeBuddyCollectorV2::reset(int obsIdA) {
  MetOfficeBuddyCollector::reset(obsIdA);
  buddyPairs_.clear();
}


}  // namespace ufo
