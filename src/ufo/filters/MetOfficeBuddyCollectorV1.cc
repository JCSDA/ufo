/*
 * (C) Copyright 2020 Met Office UK
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/filters/MetOfficeBuddyCollectorV1.h"

#include "ufo/filters/MetOfficeBuddyPair.h"

namespace ufo {

MetOfficeBuddyCollectorV1::MetOfficeBuddyCollectorV1(const MetOfficeBuddyCheckParameters &options,
                                                     const std::vector<float> &latitudes,
                                                     const std::vector<float> &longitudes,
                                                     const std::vector<int> &stationIds)
  : MetOfficeBuddyCollector(options, latitudes, longitudes, stationIds)
{}

void MetOfficeBuddyCollectorV1::examinePotentialBuddy(int obsIdB) {
  if (stationIds_[obsIdA_] == stationIds_[obsIdB]) {
    if (numBuddiesWithSameStationId_ >= options_.maxNumBuddiesWithSameStationId)
      return;
    // wsmigaj: This counter is incremented here in the original Fortran implementation,
    // but it would be more logical to increment it only if observation B turns out
    // to be a buddy of observation A.
    ++numBuddiesWithSameStationId_;
  }

  if (std::abs(latitudes_[obsIdA_] - latitudes_[obsIdB]) <= maxLatDifferenceBetweenBuddiesInDeg_) {
    potentialBuddies_.push_back(obsIdB);
    ++numBuddiesInCurrentBand_;
    ++totalNumBuddies_;
  }
}

void MetOfficeBuddyCollectorV1::appendBuddyPairsTo(
    std::vector<MetOfficeBuddyPair> &buddyPairs) const {
  // Calculate horizontal distance between the two obs for each pair.
  // Discard if > options_.searchRadius, otherwise calculate bearing of A from B
  // and reciprocal bearing of B from A. Store the results in 'buddyPairs'.

  for (int obsIdB : potentialBuddies_) {
    double deltaLatInRad, deltaLonInRad, distanceInKm;
    calcDeltaLatLonAndDistanceTo(obsIdB, deltaLatInRad, deltaLonInRad, distanceInKm);

    // wsmigaj: This check is done here in the original Fortran implementation, but I think it
    // would be better to do it already in examinePotentialBuddy(). That would ensure that the
    // search for buddies is terminated only when enough *genuine* (rather than potential) buddies
    // are found.
    if (distanceInKm <= options_.searchRadius)
      buddyPairs.push_back(createBuddyPair(obsIdB, deltaLatInRad, deltaLonInRad, distanceInKm));
  }
}

void MetOfficeBuddyCollectorV1::reset(int obsIdA) {
  MetOfficeBuddyCollector::reset(obsIdA);
  potentialBuddies_.clear();
}

}  // namespace ufo
