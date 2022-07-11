/*
 * (C) Crown copyright 2021, Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_UTILS_EQUISPACEDBINSELECTORBASE_H_
#define UFO_UTILS_EQUISPACEDBINSELECTORBASE_H_

#include <boost/optional.hpp>

namespace ufo
{

/// \brief A finite or infinite collection of non-overlapping intervals (_bins_) of the same width.
///
/// Call the bin() function to find the bin containing a particular value.
class EquispacedBinSelectorBase {
 public:
  // If necessary, these could be made template parameters.
  typedef double ValueType;
  typedef int IndexType;

  virtual ~EquispacedBinSelectorBase() {}

  /// \brief Return the index of the bin containing \p value, or the nearest bin
  /// if \p value lies outside all bins.
  virtual IndexType bin(ValueType value) const = 0;

  /// \brief Return the number of bins or boost::none if the bin collection is infinite.
  virtual boost::optional<IndexType> numBins() const = 0;

  /// \brief Return the width of each bin.
  virtual ValueType binWidth() const = 0;

  /// \brief Return the inverse of the width of each bin.
  virtual ValueType inverseBinWidth() const = 0;

  /// \brief Return the value lying at the center of the bin with index \p bin.
  virtual ValueType binCenter(IndexType bin) const = 0;
};

}  // namespace ufo

#endif  // UFO_UTILS_EQUISPACEDBINSELECTORBASE_H_
