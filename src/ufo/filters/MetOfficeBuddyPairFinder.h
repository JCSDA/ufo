/*
 * (C) Copyright 2020 Met Office UK
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_FILTERS_METOFFICEBUDDYPAIRFINDER_H_
#define UFO_FILTERS_METOFFICEBUDDYPAIRFINDER_H_

#include <memory>
#include <vector>

#include "ufo/filters/MetOfficeBuddyPair.h"

namespace util {
class DateTime;
}

namespace ufo {

class MetOfficeBuddyCheckParameters;
class MetOfficeBuddyCollector;

/// \brief Finds pairs of close observations ("buddies") to check against each other.
class MetOfficeBuddyPairFinder {
 public:
  /// \brief Constructor.
  ///
  /// \param pressures Optional -- may be null.
  MetOfficeBuddyPairFinder(const MetOfficeBuddyCheckParameters &options,
                           const std::vector<float> &latitudes,
                           const std::vector<float> &longitudes,
                           const std::vector<util::DateTime> &datetimes,
                           const std::vector<float> *pressures,
                           const std::vector<int> &stationIds);

  /// \brief Returns a list of MetOfficeBuddyPair objects representing pairs of "buddy" observations
  /// that should be checked against each other.
  std::vector<MetOfficeBuddyPair> findBuddyPairs(const std::vector<size_t> &validObsIds);

 private:
  /// \brief Sorts observations in an order facilitating rapid search for buddies.
  ///
  /// See the OPS Scientific Documentation Paper 2, section 3.3.
  ///
  /// \param[in] validObsIds
  ///   IDs of valid observations.
  /// \param[out] validObsIdsInSortOrder
  ///   IDs of valid observations sorted by zonal band index, longitude, -latitude, air pressure
  ///   (if available) and time.
  /// \param[out] bandLbounds
  ///   On output, a vector of length (options_.numZonalBands + 1) such that
  ///   [\p bandLbounds[i], \p bandLbounds[i + 1]) is the half-open range of indices of elements of
  ///   \p validObsIdsInSortOrder representing the IDs of observations from ith zonal band.
  void sortObservations(const std::vector<size_t> &validObsIds,
                        std::vector<int> &validObsIdsInSortOrder,
                        std::vector<int> &bandLbounds);

  /// \brief Finds pairs of observations to be considered as buddies. Calculates the distance and
  /// mutual orientation of each pair of buddies.
  ///
  /// See the OPS Scientific Documentation Paper 2, sections 3.4 and 3.5.
  ///
  /// \param validObsIdsInSortOrder, bandLbounds
  ///   Outputs produced by sortObservations().
  ///
  /// \returns Vector of pairs of observations to be considered as buddies.
  std::vector<MetOfficeBuddyPair> pairObservations(const std::vector<int> &validObsIdsInSortOrder,
                                                   const std::vector<int> &bandLbounds);

  std::unique_ptr<MetOfficeBuddyCollector> makeBuddyCollector() const;

  float getLongitudeSearchRangeHalfWidth(int bandIndex, float bandWidth) const;

 private:
  const MetOfficeBuddyCheckParameters &options_;
  const std::vector<float> &latitudes_;
  const std::vector<float> &longitudes_;
  const std::vector<util::DateTime> &datetimes_;
  const std::vector<float> *pressures_;
  const std::vector<int> &stationIds_;
};

}  // namespace ufo

#endif  // UFO_FILTERS_METOFFICEBUDDYPAIRFINDER_H_
