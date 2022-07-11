/*
 * (C) Crown copyright 2021, Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_UTILS_ROUNDINGEQUISPACEDBINSELECTOR_H_
#define UFO_UTILS_ROUNDINGEQUISPACEDBINSELECTOR_H_

#include <algorithm>
#include <cmath>

#include "eckit/exception/Exceptions.h"
#include "ufo/utils/EquispacedBinSelectorBase.h"

namespace ufo
{

/// \brief Represents an infinite set of consecutive intervals (_bins_) of the same width.
///
/// Bin 0 is open; bins 1, 2 etc. are right-open; bins -1, -2 etc. are left-open.
///
/// Call the bin() function to find the bin containing a particular value.
class RoundingEquispacedBinSelector : public EquispacedBinSelectorBase {
 public:
  // If necessary, these could be made template parameters.
  typedef double ValueType;
  typedef int IndexType;

  /// \brief Partition the real axis into bins of width \p binWidth, with the center of bin 0
  /// located at \p bin0Center.
  explicit RoundingEquispacedBinSelector(ValueType binWidth, ValueType bin0Center = 0)
    : binWidth_(binWidth),
      inverseBinWidth_(1 / binWidth_),
      bin0Center_(bin0Center)
  {
    ASSERT_MSG(binWidth > 0, "Bin width must be positive");
  }

  IndexType bin(ValueType value) const override {
    IndexType binIndex = std::lround((value - bin0Center_) * inverseBinWidth_);
    return binIndex;
  }

  boost::optional<IndexType> numBins() const override {
    return boost::none;
  }

  ValueType binWidth() const override {
    return binWidth_;
  }

  ValueType inverseBinWidth() const override {
    return inverseBinWidth_;
  }

  ValueType binCenter(IndexType bin) const override {
    return bin0Center_ + binWidth_ * bin;
  }

 private:
  ValueType binWidth_;
  ValueType inverseBinWidth_;
  ValueType bin0Center_;
};

}  // namespace ufo

#endif  // UFO_UTILS_ROUNDINGEQUISPACEDBINSELECTOR_H_
