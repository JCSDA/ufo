/*
 * (C) Copyright 2019 Met Office UK
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_UTILS_SPATIALBINSELECTOR_H_
#define UFO_UTILS_SPATIALBINSELECTOR_H_

#include <cmath>
#include <vector>

#include "ufo/utils/Constants.h"
#include "ufo/utils/TruncatingEquispacedBinSelector.h"

namespace ufo {

enum class SpatialBinCountRoundingMode {
  DOWN,    //< Round down
  NEAREST  //< Round to the nearest integer
};

/// \brief Represents a partition of a sphere into a number of subsets (_bins_).
class SpatialBinSelector {
 public:
  // If necessary, these could be made template parameters.
  typedef float ValueType;
  typedef int IndexType;

 private:
  static constexpr ValueType latitudeLowerBound_ = -90;
  static constexpr ValueType latitudeUpperBound_ = 90;
  static constexpr ValueType longitudeLowerBound_ = 0;
  static constexpr ValueType longitudeUpperBound_ = 360;
  static constexpr ValueType opsCompatibilityModeLongitudeLowerBound_ = -360;
  static constexpr ValueType opsCompatibilityModeLongitudeUpperBound_ =  720;
  static constexpr int opsCompatibilityModeRelativeLongitudeRange_ = 3;

 public:
  /// \brief Partition a sphere into bins whose centers lie on a reduced Gaussian grid.
  ///
  /// \param numLatitudeBins
  ///   The number of zonal bands of bins into which the sphere is split.
  /// \param roundingMode
  ///   - If set to NEAREST, the number of bins in each zonal band is chosen so that the bin width
  ///     in the zonal direction is as close as possible to that in the meridional direction.
  ///   - If set to DOWN, the number of bins is chosen so that the bin width in the zonal direction
  ///     is as small as possible, but no smaller than in the meridional direction.
  /// \param horizontalMesh
  ///   Cell size (in km) along the meridians.
  /// \param metOfficeOpsCompatibilityMode
  ///   If true, the compatibility mode with the Met Office OPS system will be activated.
  ///   This will stop longitudeBin() from clamping longitudes to the interval [0, 360] deg.
  ///
  ///   Activation of this mode in new code is discouraged; even if input longitudes are e.g. in
  ///   the range [-180, 180] deg, it is normally better to wrap them to the [0, 360] deg range
  ///   before passing them to longitudeBin(). Otherwise there's a risk that points lying exactly
  ///   at -180 or 180 deg will be put into a bin of their own.
  /// \param partitionLongitudeBinsUsingMesh
  ///   - If true, the partitioning of longitude bins per latitude bin index is calculated
  ///     using the horizontal mesh cell size.
  ///   - If false, the integer number of latitude bins is used to calculate how longitude bins are
  ///     partitioned.
  SpatialBinSelector(IndexType numLatitudeBins, SpatialBinCountRoundingMode roundingMode,
                     float horizontalMesh, bool metOfficeOpsCompatibilityMode = false,
                     bool partitionLongitudeBinsUsingMesh = false);

  /// \brief Partition a sphere into bins whose centers lie on a regular Gaussian grid.
  ///
  /// \param numLatitudeBins
  ///   The number of zonal bands of bins into which the sphere is split.
  /// \param numLongitudeBins
  ///   The number of meridional bands of bins into which the sphere is split.
  SpatialBinSelector(IndexType numLatitudeBins, IndexType numLongitudeBins,
                     bool metOfficeOpsCompatibilityMode = false);

  /// \brief Return the index of the zonal band of bins containing points with a given latitude
  /// (in degrees, assumed to lie in the interval [-90, 90]).
  IndexType latitudeBin(ValueType latitude) const {
    return latitudeBinSelector_.bin(latitude);
  }

  /// \brief Return the index of the bin within the zonal band of index \p latitudeBin
  /// containing points with a given longitude (in degrees).
  ///
  /// The longitude is assumed to lie in the interval [0, 360] unless the compatibility mode with
  /// the Met Office OPS system is in effect.
  IndexType longitudeBin(IndexType latitudeBin, ValueType longitude) const {
    return longitudeBinSelectors_[latitudeBin].bin(longitude);
  }

  /// \brief Return the latitude at the center of the zonal band of bins with index \p latitudeBin.
  ValueType latitudeBinCenter(IndexType latitudeBin) const {
    return latitudeBinSelector_.binCenter(latitudeBin);
  }

  /// \brief Return the longitude at the center of the bin with index \p latitudeBin lying in the
  /// zonal band of bins with index \p latitudeBin.
  ValueType longitudeBinCenter(IndexType latitudeBin, IndexType longitudeBin) const {
    return longitudeBinSelectors_[latitudeBin].binCenter(longitudeBin);
  }

  /// \brief Return the number of bins into which the sphere is split.
  IndexType totalNumBins() const;

  /// \brief Return the width of each zonal band of bins.
  ValueType latitudeBinWidth() const {
    return latitudeBinSelector_.binWidth();
  }

  /// \brief Return the zonal width of each bin in the band of bins with index \p latitudeBin.
  ValueType longitudeBinWidth(IndexType latitudeBin) const {
    return longitudeBinSelectors_[latitudeBin].binWidth();
  }

  /// \brief Return the inverse of the width of each zonal band of bins.
  ValueType inverseLatitudeBinWidth() const {
    return latitudeBinSelector_.inverseBinWidth();
  }

  /// \brief Return the inverse of the zonal width of each bin in the band of bins with index
  /// \p latitudeBin.
  ValueType inverseLongitudeBinWidth(IndexType latitudeBin) const {
    return longitudeBinSelectors_[latitudeBin].inverseBinWidth();
  }

  /// \brief Return \p idealNumBins rounded to a positive integer according to the rounding
  /// strategy \p roundingMode.
  static IndexType roundNumBins(float idealNumBins, SpatialBinCountRoundingMode roundingMode);

 private:
  bool metOfficeOpsCompatibilityMode_;
  TruncatingEquispacedBinSelector latitudeBinSelector_;
  std::vector<TruncatingEquispacedBinSelector> longitudeBinSelectors_;
};

}  // namespace ufo

#endif  // UFO_UTILS_SPATIALBINSELECTOR_H_
