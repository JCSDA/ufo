/*
 * (C) Copyright 2020 Met Office UK
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/filters/MetOfficeBuddyPairFinder.h"

#include <algorithm>

#include "oops/util/Logger.h"
#include "ufo/filters/MetOfficeBuddyCheckParameters.h"
#include "ufo/filters/MetOfficeBuddyCollectorV1.h"
#include "ufo/filters/MetOfficeBuddyCollectorV2.h"
#include "ufo/utils/RecursiveSplitter.h"

#include <boost/make_unique.hpp>

namespace ufo {

namespace {

float zonalBandWidth(int numBands) {
  return 180.0f / numBands;
}

}  // namespace

MetOfficeBuddyPairFinder::MetOfficeBuddyPairFinder(const MetOfficeBuddyCheckParameters &options,
                           const std::vector<float> &latitudes,
                           const std::vector<float> &longitudes,
                           const std::vector<util::DateTime> &datetimes,
                           const std::vector<float> *pressures,
                           const std::vector<int> &stationIds)
  : options_(options), latitudes_(latitudes), longitudes_(longitudes),
    datetimes_(datetimes), pressures_(pressures), stationIds_(stationIds)
{}

std::vector<MetOfficeBuddyPair> MetOfficeBuddyPairFinder::findBuddyPairs(
    const std::vector<size_t> & validObsIds) {
  std::vector<int> validObsIdsInSortOrder;
  std::vector<int> bandLbounds;
  sortObservations(validObsIds, validObsIdsInSortOrder, bandLbounds);
  return pairObservations(validObsIdsInSortOrder, bandLbounds);
}

void MetOfficeBuddyPairFinder::sortObservations(const std::vector<size_t> & validObsIds,
                                                std::vector<int> &validObsIdsInSortOrder,
                                                std::vector<int> &bandLbounds)
{
  // Initialize output parameters

  validObsIdsInSortOrder.clear();
  validObsIdsInSortOrder.reserve(validObsIds.size());
  bandLbounds.assign(options_.numZonalBands + 1, 0);

  // Identify the band containing each valid observation

  const float bandWidth = zonalBandWidth(options_.numZonalBands);
  std::vector<int> bandIndices(validObsIds.size());
  for (size_t validObsIndex = 0; validObsIndex < validObsIds.size(); ++validObsIndex) {
    const size_t obsId = validObsIds[validObsIndex];
    // Use 89.999 to round up values, e.g. ranges are 70.0 to 74.999
    // (There may be a more elegant way to do this, but use this one for compatibility with OPS.)
    bandIndices[validObsIndex] =
        std::max(0, static_cast<int>((89.999f - latitudes_[obsId]) / bandWidth));
  }

  // Sort observations

  RecursiveSplitter splitter(validObsIds.size());
  splitter.groupBy(bandIndices);
  splitter.sortGroupsBy(
        [this, &validObsIds](size_t obsIndex)
        {
          size_t obsId = validObsIds[obsIndex];
          return std::make_tuple(longitudes_[obsId], -latitudes_[obsId],
                                 pressures_ ? (*pressures_)[obsId] : 0.0f, datetimes_[obsId]);
        });

  // Fill the validObsIdsInSortOrder and bandLbounds vectors

  int currentBandIndex = -1;
  for (auto group : splitter.groups()) {
    for (size_t validObsIndex : group) {
      const size_t obsId = validObsIds[validObsIndex];
      for (int band = currentBandIndex + 1; band <= bandIndices[validObsIndex]; ++band) {
        bandLbounds[band] = validObsIdsInSortOrder.size();
      }
      currentBandIndex = bandIndices[validObsIndex];
      validObsIdsInSortOrder.push_back(obsId);
    }
  }
  for (int band = currentBandIndex + 1; band < bandLbounds.size(); ++band) {
    bandLbounds[band] = validObsIds.size();
  }

  oops::Log::trace() << "Buddy check: " << validObsIds.size() << " input observations" << std::endl;
}

std::vector<MetOfficeBuddyPair> MetOfficeBuddyPairFinder::pairObservations(
    const std::vector<int> &validObsIdsInSortOrder,
    const std::vector<int> &bandLbounds) {

  std::vector<MetOfficeBuddyPair> pairs;

  // Initialise variables
  const float bandWidth = zonalBandWidth(options_.numZonalBands);
  // eqn 3.1
  const float searchDLat = Constants::rad2deg * options_.searchRadius / Constants::mean_earth_rad;
  const float searchDLatB = searchDLat + 0.5f * bandWidth;
  const int numSearchBands = static_cast<int>(searchDLat / bandWidth) + 1;

  typedef std::vector<int>::const_iterator ObsIdIt;
  typedef std::vector<int>::const_reverse_iterator ObsIdRevIt;

  std::vector<ObsIdIt> bandBegins(options_.numZonalBands);
  std::vector<ObsIdIt> bandEnds(options_.numZonalBands);

  for (int bandIndex = 0; bandIndex < options_.numZonalBands; ++bandIndex) {
    bandBegins[bandIndex] = validObsIdsInSortOrder.begin() + bandLbounds[bandIndex];
    bandEnds[bandIndex] = validObsIdsInSortOrder.begin() + bandLbounds[bandIndex + 1];
  }

  std::vector<ObsIdIt> firstObsToCheckInBands(options_.numZonalBands);

  // Collects buddies of a single observation. When we're done with that observation, the collected
  // list of buddies is extracted into 'pairs' and the collector is reset.
  std::unique_ptr<MetOfficeBuddyCollector> buddyCollector = makeBuddyCollector();

  // Iterate over all bands
  for (int jBandA = 0; jBandA < options_.numZonalBands; ++jBandA) {
    const float lonSearchRangeHalfWidth = getLongitudeSearchRangeHalfWidth(jBandA, bandWidth);

    const int firstBandToSearch = jBandA;
    const int lastBandToSearch = std::min(options_.numZonalBands.value() - 1,
                                          jBandA + numSearchBands);

    firstObsToCheckInBands = bandBegins;

    // Iterate over observations in (jBandA)th band
    for (ObsIdIt obsIdItA = bandBegins[jBandA]; obsIdItA != bandEnds[jBandA]; ++obsIdItA) {
      const int obsIdA = *obsIdItA;

      const float minLonToCheck = longitudes_[obsIdA] - lonSearchRangeHalfWidth;
      const float maxLonToCheck = longitudes_[obsIdA] + lonSearchRangeHalfWidth;

      buddyCollector->reset(obsIdA);

      // Iterate over bands that may contain buddies
      for (int jBandB = firstBandToSearch; jBandB <= lastBandToSearch; ++jBandB) {
        float midBandLatB = 90.0f - bandWidth * (jBandB + 0.5f);
        if (std::abs(latitudes_[obsIdA] - midBandLatB) > searchDLatB)
          continue;

        ObsIdIt firstObsIdItB = firstObsToCheckInBands[jBandB];
        if (jBandA == jBandB)
          firstObsIdItB = obsIdItA + 1;  // Iterate only over observations following observation A.

        // First loop: look for buddies at longitudes [minLonToCheck, maxLonToCheck]
        bool firstObsInSearchRangeFound = false;
        for (ObsIdIt obsIdItB = firstObsIdItB; obsIdItB != bandEnds[jBandB]; ++obsIdItB) {
          const int obsIdB = *obsIdItB;
          if (longitudes_[obsIdB] < minLonToCheck)
            continue;
          if (longitudes_[obsIdB] > maxLonToCheck)
            break;
          if (!firstObsInSearchRangeFound) {
            firstObsToCheckInBands[jBandB] = obsIdItB;
            firstObsInSearchRangeFound = true;
          }

          buddyCollector->examinePotentialBuddy(obsIdB);
          if (buddyCollector->foundEnoughBuddies())
            goto FinishProcessingObsA;
          if (buddyCollector->foundEnoughBuddiesInCurrentBand())
            goto FinishProcessingBandB;
        }

        // Optional second loop: look for buddies at longitudes [-180, maxLonToCheck - 360]
        if (maxLonToCheck > 180 && (jBandA != jBandB || lonSearchRangeHalfWidth < 180)) {
          // Observation A is near band end (+180); wrap around and check the band start too.
          float wrappedMaxLonToCheck = maxLonToCheck - 360;
          for (ObsIdIt obsIdItB = bandBegins[jBandB]; obsIdItB != bandEnds[jBandB]; ++obsIdItB) {
            const int obsIdB = *obsIdItB;
            if (longitudes_[obsIdB] > wrappedMaxLonToCheck ||
                longitudes_[obsIdB] >= minLonToCheck /* visited already in the first loop */)
              break;

            buddyCollector->examinePotentialBuddy(obsIdB);
            if (buddyCollector->foundEnoughBuddies())
              goto FinishProcessingObsA;
            if (buddyCollector->foundEnoughBuddiesInCurrentBand())
              goto FinishProcessingBandB;
          }
        }

        // Optional third loop: look for buddies at longitudes [minLonToCheck + 360, 180]
        if (minLonToCheck < -180 && jBandA != jBandB) {
          // Observation A is near band start (-180); wrap around and check the band end too.
          float wrappedMinLonToCheck = minLonToCheck + 360;
          for (ObsIdRevIt obsIdItB(bandEnds[jBandB]), reverseBandEnd(bandBegins[jBandB]);
               obsIdItB != reverseBandEnd; ++obsIdItB) {
            const int obsIdB = *obsIdItB;
            if (longitudes_[obsIdB] < wrappedMinLonToCheck ||
                longitudes_[obsIdB] <= maxLonToCheck /* visited already in the first loop */)
              break;

            buddyCollector->examinePotentialBuddy(obsIdB);
            if (buddyCollector->foundEnoughBuddies())
              goto FinishProcessingObsA;
            if (buddyCollector->foundEnoughBuddiesInCurrentBand())
              goto FinishProcessingBandB;
          }
        }

FinishProcessingBandB:
        buddyCollector->startProcessingNextBand();
      }  // end of secondary loop over bands (jBandB)
FinishProcessingObsA:
      buddyCollector->appendBuddyPairsTo(pairs);
    }  // end of main loop over observations (obsIdItA)
  }  // end of main loop over bands (jBandA)

  oops::Log::trace() << "Found " << pairs.size() << " buddy pairs.\n";

  return pairs;
}

std::unique_ptr<MetOfficeBuddyCollector> MetOfficeBuddyPairFinder::makeBuddyCollector() const {
  if (options_.useLegacyBuddyCollector)
    return boost::make_unique<MetOfficeBuddyCollectorV1>(options_, latitudes_,
                                                         longitudes_, stationIds_);
  else
    return boost::make_unique<MetOfficeBuddyCollectorV2>(options_, latitudes_,
                                                         longitudes_, stationIds_);
}

float MetOfficeBuddyPairFinder::getLongitudeSearchRangeHalfWidth(int bandIndex,
                                                                 float bandWidth) const {
  const float earthRadius = Constants::mean_earth_rad;
  const float deg2rad = static_cast<float>(Constants::deg2rad);
  const float rad2deg = static_cast<float>(Constants::rad2deg);

  const float midBandLatA = 90.0f - bandWidth * (bandIndex + 0.5f);
  // eqn 3.2a
  const float rad = earthRadius * std::cos((std::abs(midBandLatA) + bandWidth * 0.5f) * deg2rad);
  if (rad <= 10.0f)
    return 360.0f;  // Adjacent to pole
  else
    return rad2deg * options_.searchRadius / rad;  // eqn 3.2b
}

}  // namespace ufo

