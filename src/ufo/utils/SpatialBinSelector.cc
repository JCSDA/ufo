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
#include "ufo/utils/EquispacedBinSelector.h"
#include "ufo/utils/SpatialBinSelector.h"

namespace ufo {

SpatialBinSelector::SpatialBinSelector(IndexType numLatitudeBins,
                                       SpatialBinCountRoundingMode roundingMode)
  : latitudeBinSelector_(latitudeLowerBound_, latitudeUpperBound_, numLatitudeBins) {
  longitudeBinSelectors_.reserve(numLatitudeBins);
  for (IndexType latBin = 0; latBin < numLatitudeBins; ++latBin) {
    ValueType latBinCenter = latitudeBinCenter(latBin);
    // NOTE: original code was rounding *down*.

    const int equatorToMeridianLengthRatio = 2;
    const float tentativeNumLongitudeBins =
        equatorToMeridianLengthRatio * numLatitudeBins *
        std::cos(latBinCenter * static_cast<float>(Constants::deg2rad));
    IndexType numLonBins = static_cast<IndexType>(
          roundingMode == SpatialBinCountRoundingMode::DOWN ?
            tentativeNumLongitudeBins :
            std::round(tentativeNumLongitudeBins));
    numLonBins = std::max(1, numLonBins);

    longitudeBinSelectors_.emplace_back(
          static_cast<ValueType>(longitudeLowerBound_),
          static_cast<ValueType>(longitudeUpperBound_), numLonBins);
  }
}

SpatialBinSelector::SpatialBinSelector(IndexType numLatitudeBins, IndexType numLongitudeBins)
  : latitudeBinSelector_(latitudeLowerBound_, latitudeUpperBound_, numLatitudeBins),
    longitudeBinSelectors_(numLatitudeBins,
                           EquispacedBinSelector(longitudeLowerBound_, longitudeUpperBound_,
                                                 numLongitudeBins))
{}

SpatialBinSelector::IndexType SpatialBinSelector::totalNumBins() const {
  size_t n = 0;
  for (const EquispacedBinSelector & selector : longitudeBinSelectors_)
    n += selector.numBins();
  return n;
}

}  // namespace ufo
