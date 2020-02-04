/*
 * (C) Copyright 2019 Met Office UK
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_UTILS_EQUISPACEDBINSELECTOR_H_
#define UFO_UTILS_EQUISPACEDBINSELECTOR_H_

#include <algorithm>

#include "eckit/exception/Exceptions.h"

namespace ufo
{

/// \brief Represents a set of consecutive intervals (_bins_) of the same width.
///
/// Call the bin() function to find the bin containing a particular value.
class EquispacedBinSelector
{
 public:
  // If necessary, these could be made template parameters.
  typedef float ValueType;
  typedef int IndexType;

  /// \brief Partition the interval [\p lowerBound, \p upperBound] into \p numBins
  /// bins of the same width.
  EquispacedBinSelector(ValueType lowerBound, ValueType upperBound, IndexType numBins)
    : lowerBound_(lowerBound),
      binWidth_((upperBound - lowerBound) / numBins),
      inverseBinWidth_(1 / binWidth_),
      numBins_(numBins)
  {
    ASSERT_MSG(upperBound > lowerBound, "Upper bound must be larger than lower bound");
    ASSERT_MSG(numBins > 0, "Number of bins must be positive");
  }

  /// \brief Return the (0-based) index of the bin containing \p value, or the nearest bin
  /// if \p value lies outside all bins.
  IndexType bin(ValueType value) const {
    IndexType binIndex = static_cast<IndexType>((value - lowerBound_) * inverseBinWidth_);
    binIndex = std::max(IndexType(0), binIndex);
    binIndex = std::min(numBins_ - 1, binIndex);
    return binIndex;
  }

  /// \brief Return the number of bins.
  IndexType numBins() const {
    return numBins_;
  }

  /// \brief Return the width of each bin.
  ValueType binWidth() const {
    return binWidth_;
  }

  /// \brief Return the inverse of the width of each bin.
  ValueType inverseBinWidth() const {
    return inverseBinWidth_;
  }

  /// \brief Return the value lying at the center of the bin with index \p bin.
  ValueType binCenter(IndexType bin) const {
    return lowerBound_ + binWidth_ * (bin + ValueType(0.5));
  }

 private:
  ValueType lowerBound_;
  ValueType binWidth_;
  ValueType inverseBinWidth_;
  IndexType numBins_;
};

}  // namespace ufo

#endif  // UFO_UTILS_EQUISPACEDBINSELECTOR_H_
