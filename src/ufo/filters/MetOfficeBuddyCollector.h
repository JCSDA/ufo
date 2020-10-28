/*
 * (C) Copyright 2020 Met Office UK
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_FILTERS_METOFFICEBUDDYCOLLECTOR_H_
#define UFO_FILTERS_METOFFICEBUDDYCOLLECTOR_H_

#include <cassert>
#include <vector>

#include "MetOfficeBuddyCheckParameters.h"

namespace ufo {

struct MetOfficeBuddyPair;

/// \brief Interface of classes used by MetOfficeBuddyPairFinder to select observations used as
/// buddies of other observations during the buddy check.
class MetOfficeBuddyCollector {
 public:
  /// \brief Constructor.
  ///
  /// \param options
  ///   Buddy check parameters.
  /// \param latitudes
  ///   Observation latitudes.
  /// \param longitudes
  ///   Observation longitudes.
  /// \param stationIds
  ///   IDs of stations that have collected the observations.
  MetOfficeBuddyCollector(const MetOfficeBuddyCheckParameters &options,
                          const std::vector<float> &latitudes,
                          const std::vector<float> &longitudes,
                          const std::vector<int> &stationIds);

  virtual ~MetOfficeBuddyCollector() = default;
  // This class is expected to be used polymorphically -- via a pointer or reference --
  // so there doesn't seem to be much need for a copy or move constructor. Delete them for now.
  MetOfficeBuddyCollector(const MetOfficeBuddyCollector &) = delete;
  MetOfficeBuddyCollector(MetOfficeBuddyCollector &&) = delete;
  MetOfficeBuddyCollector & operator=(const MetOfficeBuddyCollector &) = delete;
  MetOfficeBuddyCollector & operator=(MetOfficeBuddyCollector &&) = delete;

  /// \brief Prepare the object for examination of potential buddies of the observation with ID
  /// \p obsIdA.
  ///
  /// \note This function must be called before any calls to other member functions.
  virtual void reset(int obsIdA) = 0;

  /// \brief Check if the observation with ID \p obsIdB can be selected as a buddy of obsIdB. If
  /// so, record its ID internally.
  virtual void examinePotentialBuddy(int obsIdB) = 0;

  /// \brief Called to indicate that observations passed to subsequent calls to
  /// examinePotentialBuddy() will belong to a new zonal band.
  void startProcessingNextBand();

  /// \brief Returns true if the number of observations selected as buddies of the observation
  /// passed to reset() since the last call to startProcessingNextBand() has reached the
  /// limit set by the \c max_num_buddies_from_single_band parameter.
  bool foundEnoughBuddiesInCurrentBand() const;

  /// \brief Returns true if the number of observations selected as buddies of the observation
  /// passed to reset() since the last call to startProcessingNextBand() has reached the
  /// limit set by the \c max_total_num_buddies parameter.
  bool foundEnoughBuddies() const;

  /// \brief Extend \p buddyPairs with MetOfficeBuddyPair objects storing the properies of all buddy
  /// pairs found since the last call to reset().
  virtual void appendBuddyPairsTo(std::vector<MetOfficeBuddyPair> &buddyPairs) const = 0;

 protected:
  void calcDeltaLatLonAndDistanceTo(int obsIdB,
                                      double &deltaLatInRad, double &deltaLonInRad,
                                      double &distanceInKm) const;

  MetOfficeBuddyPair createBuddyPair(int obsIdB,
                                     double deltaLatInRad, double deltaLonInRad,
                                     double distanceInKm) const;

 protected:
  const MetOfficeBuddyCheckParameters &options_;
  const std::vector<float> &latitudes_;
  const std::vector<float> &longitudes_;
  const std::vector<int> &stationIds_;
  double maxLatDifferenceBetweenBuddiesInDeg_;
  int obsIdA_ = 0;
  int numBuddiesInCurrentBand_ = 0;
  int numBuddiesWithSameStationId_ = 0;
  int totalNumBuddies_ = 0;
};


inline void MetOfficeBuddyCollector::startProcessingNextBand() {
  numBuddiesInCurrentBand_ = 0;
}

inline bool MetOfficeBuddyCollector::foundEnoughBuddiesInCurrentBand() const {
  assert(numBuddiesInCurrentBand_ <= options_.maxNumBuddiesFromSingleBand);
  return numBuddiesInCurrentBand_ == options_.maxNumBuddiesFromSingleBand;
}

inline bool MetOfficeBuddyCollector::foundEnoughBuddies() const {
  assert(totalNumBuddies_ <= options_.maxTotalNumBuddies);
  return totalNumBuddies_ == options_.maxTotalNumBuddies;
}

inline void MetOfficeBuddyCollector::reset(int obsIdA) {
  obsIdA_ = obsIdA;
  numBuddiesInCurrentBand_ = 0;
  numBuddiesWithSameStationId_ = 0;
  totalNumBuddies_ = 0;
}

}  // namespace ufo

#endif  // UFO_FILTERS_METOFFICEBUDDYCOLLECTOR_H_
