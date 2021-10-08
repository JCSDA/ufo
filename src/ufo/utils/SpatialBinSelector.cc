/*
 * (C) Copyright 2019 Met Office UK
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <algorithm>
#include <cmath>
#include <vector>

#include "ufo/utils/Constants.h"
#include "ufo/utils/SpatialBinSelector.h"
#include "ufo/utils/TruncatingEquispacedBinSelector.h"

namespace ufo {

SpatialBinSelector::SpatialBinSelector(IndexType numLatitudeBins,
                                       SpatialBinCountRoundingMode roundingMode,
                                       float horizontalMesh,
                                       bool metOfficeOpsCompatibilityMode,
                                       bool partitionLongitudeBinsUsingMesh)
  : metOfficeOpsCompatibilityMode_(metOfficeOpsCompatibilityMode),
    latitudeBinSelector_(latitudeLowerBound_, latitudeUpperBound_, numLatitudeBins) {
  longitudeBinSelectors_.reserve(numLatitudeBins);
  for (IndexType latBin = 0; latBin < numLatitudeBins; ++latBin) {
    ValueType latBinCenter = latitudeBinCenter(latBin);

    const int equatorToMeridianLengthRatio = 2;
    float tentativeNumLongitudeBins;
    if (partitionLongitudeBinsUsingMesh) {
      const float meridianLength = M_PI * Constants::mean_earth_rad;
      tentativeNumLongitudeBins =
            equatorToMeridianLengthRatio * (meridianLength/horizontalMesh) *
            std::cos(latBinCenter * static_cast<float>(Constants::deg2rad));
    } else {
      tentativeNumLongitudeBins =
            equatorToMeridianLengthRatio * numLatitudeBins *
            std::cos(latBinCenter * static_cast<float>(Constants::deg2rad));
    }
    const IndexType numLonBins = roundNumBins(tentativeNumLongitudeBins, roundingMode);

    if (metOfficeOpsCompatibilityMode)
      longitudeBinSelectors_.emplace_back(
            static_cast<ValueType>(opsCompatibilityModeLongitudeLowerBound_),
            static_cast<ValueType>(opsCompatibilityModeLongitudeUpperBound_),
            opsCompatibilityModeRelativeLongitudeRange_ * numLonBins);
    else
      longitudeBinSelectors_.emplace_back(
            static_cast<ValueType>(longitudeLowerBound_),
            static_cast<ValueType>(longitudeUpperBound_), numLonBins);
  }
}

SpatialBinSelector::SpatialBinSelector(IndexType numLatitudeBins, IndexType numLongitudeBins,
                                       bool metOfficeOpsCompatibilityMode)
  : metOfficeOpsCompatibilityMode_(metOfficeOpsCompatibilityMode),
    latitudeBinSelector_(latitudeLowerBound_, latitudeUpperBound_, numLatitudeBins),
    longitudeBinSelectors_(numLatitudeBins,
                           metOfficeOpsCompatibilityMode ?
                             TruncatingEquispacedBinSelector(
                               opsCompatibilityModeLongitudeLowerBound_,
                               opsCompatibilityModeLongitudeUpperBound_,
                               opsCompatibilityModeRelativeLongitudeRange_ * numLongitudeBins) :
                             TruncatingEquispacedBinSelector(
                               longitudeLowerBound_,
                               longitudeUpperBound_,
                               numLongitudeBins))
{}

SpatialBinSelector::IndexType SpatialBinSelector::totalNumBins() const {
  size_t n = 0;
  for (const TruncatingEquispacedBinSelector & selector : longitudeBinSelectors_)
    n += *selector.numBins();
  if (metOfficeOpsCompatibilityMode_)
    n /= opsCompatibilityModeRelativeLongitudeRange_;
  return n;
}

SpatialBinSelector::IndexType SpatialBinSelector::roundNumBins(
    float idealNumBins, SpatialBinCountRoundingMode roundingMode) {
  IndexType numBins = static_cast<IndexType>(
        roundingMode == SpatialBinCountRoundingMode::DOWN ?
          idealNumBins : std::round(idealNumBins));
  return std::max(1, numBins);
}

}  // namespace ufo
