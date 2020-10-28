/*
 * (C) Copyright 2020 Met Office UK
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/filters/MetOfficeBuddyCollector.h"

#include <algorithm>
#include <cmath>

#include "oops/util/sqr.h"
#include "ufo/filters/MetOfficeBuddyPair.h"
#include "ufo/utils/Constants.h"

namespace ufo {

MetOfficeBuddyCollector::MetOfficeBuddyCollector(const MetOfficeBuddyCheckParameters &options,
                                                 const std::vector<float> &latitudes,
                                                 const std::vector<float> &longitudes,
                                                 const std::vector<int> &stationIds)
  : options_(options), latitudes_(latitudes), longitudes_(longitudes), stationIds_(stationIds)
{
  // eqn 3.1
  maxLatDifferenceBetweenBuddiesInDeg_ =
      Constants::rad2deg * options.searchRadius / Constants::mean_earth_rad;
}

void MetOfficeBuddyCollector::calcDeltaLatLonAndDistanceTo(int obsIdB,
                                                           double &deltaLatInRad,
                                                           double &deltaLonInRad,
                                                           double &distanceInKm) const {
  deltaLatInRad = (latitudes_[obsIdB] - latitudes_[obsIdA_]) * Constants::deg2rad;

  deltaLonInRad = (longitudes_[obsIdB] - longitudes_[obsIdA_]) * Constants::deg2rad;
  if (deltaLonInRad > M_PI)
    deltaLonInRad -= 2 * M_PI;
  else if (deltaLonInRad < -M_PI)
    deltaLonInRad += 2 * M_PI;

  distanceInKm = Constants::mean_earth_rad *
      std::sqrt(util::sqr(deltaLatInRad) +
                4.0 * util::sqr(std::sin(0.5 * deltaLonInRad)) *
                std::cos(latitudes_[obsIdA_] * Constants::deg2rad) *
                std::cos(latitudes_[obsIdB] * Constants::deg2rad));    // eqn 3.3
}

MetOfficeBuddyPair MetOfficeBuddyCollector::createBuddyPair(int obsIdB,
                                                            double deltaLatInRad,
                                                            double deltaLonInRad,
                                                            double distanceInKm) const {
  double rotA, rotB;
  if (distanceInKm < 10.0) {
    rotA = 0.0;  // the transformation is undefined
    rotB = 0.0;  // for distanceInKm = 0.0, use u,v components
  } else {
    // eqn 3.5
    double alpha = 0.5 * std::sin(latitudes_[obsIdA_] * Constants::deg2rad) * deltaLonInRad;
    // eqn 3.6
    double sinBeta = Constants::mean_earth_rad * deltaLatInRad * std::cos(alpha) / distanceInKm;
    sinBeta = std::min(1.0, std::max(-1.0, sinBeta));
    double beta = std::asin(sinBeta);
    rotA = alpha + beta;  // eqn 3.7
    rotB = beta - alpha;  // eqn 3.8
  }

  return MetOfficeBuddyPair(obsIdA_, obsIdB, distanceInKm, rotA, rotB);
}

}  // namespace ufo
