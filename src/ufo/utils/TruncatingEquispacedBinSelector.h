/*
 * (C) Copyright 2019 Met Office UK
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_UTILS_TRUNCATINGEQUISPACEDBINSELECTOR_H_
#define UFO_UTILS_TRUNCATINGEQUISPACEDBINSELECTOR_H_

#include <algorithm>

#include "eckit/exception/Exceptions.h"
#include "ufo/utils/EquispacedBinSelectorBase.h"

namespace ufo
{

/// \brief Represents a finite set of consecutive intervals (_bins_) of the same width, each closed
/// from the left and open from the right.
///
/// Call the bin() function to find the bin containing a particular value.
class TruncatingEquispacedBinSelector : public EquispacedBinSelectorBase {
 public:
  // If necessary, these could be made template parameters.
  typedef double ValueType;
  typedef int IndexType;

  /// \brief Partition the interval [\p lowerBound, \p upperBound) into \p numBins
  /// bins of the same width.
  TruncatingEquispacedBinSelector(ValueType lowerBound, ValueType upperBound, IndexType numBins)
    : lowerBound_(lowerBound),
      binWidth_((upperBound - lowerBound) / numBins),
      inverseBinWidth_(1 / binWidth_),
      numBins_(numBins)
  {
    ASSERT_MSG(upperBound > lowerBound, "Upper bound must be larger than lower bound");
    ASSERT_MSG(numBins > 0, "Number of bins must be positive");
  }

  IndexType bin(ValueType value) const override {
    IndexType binIndex = static_cast<IndexType>((value - lowerBound_) * inverseBinWidth_);
    binIndex = std::max(IndexType(0), binIndex);
    binIndex = std::min(numBins_ - 1, binIndex);
    return binIndex;
  }

  boost::optional<IndexType> numBins() const override {
    return numBins_;
  }

  ValueType binWidth() const override {
    return binWidth_;
  }

  ValueType inverseBinWidth() const override {
    return inverseBinWidth_;
  }

  ValueType binCenter(IndexType bin) const override {
    return lowerBound_ + binWidth_ * (bin + ValueType(0.5));
  }

 private:
  ValueType lowerBound_;
  ValueType binWidth_;
  ValueType inverseBinWidth_;
  IndexType numBins_;
};

}  // namespace ufo

#endif  // UFO_UTILS_TRUNCATINGEQUISPACEDBINSELECTOR_H_
