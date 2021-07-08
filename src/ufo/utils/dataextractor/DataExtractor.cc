/*
 * (C) Copyright 2021 Met Office UK
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <algorithm>           // sort
#include <functional>          // greater
#include <limits>              // std::numeric_limits
#include <list>                // list
#include <sstream>             // stringstream
#include <utility>             // pair

#include <boost/make_unique.hpp>

#include "eckit/exception/Exceptions.h"
#include "eckit/utils/StringTools.h"

#include "ioda/Misc/StringFuncs.h"

#include "oops/util/Logger.h"
#include "oops/util/missingValues.h"

#include "ufo/utils/dataextractor/DataExtractor.h"
#include "ufo/utils/dataextractor/DataExtractorCSVBackend.h"
#include "ufo/utils/dataextractor/DataExtractorInput.h"
#include "ufo/utils/dataextractor/DataExtractorNetCDFBackend.h"

namespace ufo {

namespace {

/// \brief Boost visitor which allows us to sort a vector.
class SortUpdateVisitor : public boost::static_visitor<void> {
 public:
  explicit SortUpdateVisitor(ufo::RecursiveSplitter &splitter) : splitter(splitter) {}

  template <typename T>
  void operator()(const std::vector<T> &coord) {
    splitter.groupBy(coord);
  }

  void operator()(const std::vector<float> &coord) {
    splitter.sortGroupsBy(
      [&coord](int indexA, int indexB) {
        return coord[static_cast<size_t>(indexA)] < coord[static_cast<size_t>(indexB)];
      });
  }

  ufo::RecursiveSplitter &splitter;
};


/// \brief Boost visitor which allows us to sort a vector.
class SortVisitor : public boost::static_visitor<void> {
 public:
  explicit SortVisitor(const ufo::RecursiveSplitter &splitter) : splitter(splitter) {}

  template <typename T>
  void operator()(std::vector<T> &coord) {
    std::vector<T> newCoord;
    newCoord.reserve(coord.size());
    for (const auto &group : splitter.groups()) {
      for (const auto &index : group) {
        newCoord.push_back(coord[index]);
      }
    }
    // Replace the coordinate with the sorted one.
    coord = std::move(newCoord);
  }

  const ufo::RecursiveSplitter &splitter;
};


/// \brief Update our extract constraint based on an exact match against the specified coordinate
/// indexing a dimension of the payload array.
///
/// \param[in] varName
///   Name of the coordinate to match against.
/// \param[in] varValues
///   Vector of values of that coordinate.
/// \param[in] obVal
///   Value to match.
/// \param[inout] range
///   On input, the range of slices of the payload array along the dimension indexed by
///   `varValues` that matches all constraints considered so far. On output, the subrange of
///   slices matching also the current constraint.
template<typename T>
void exactMatch(const std::string &varName,
                const std::vector<T> &varValues,
                const T &obVal,
                ConstrainedRange &range) {
  // Find the first and last matching index
  auto bounds = std::equal_range(varValues.begin() + range.begin(),
                                 varValues.begin() + range.end(),
                                 obVal);
  if (bounds.first == bounds.second) {
    // No matching coordinate found. If the coordinate contains a 'missing value' entry,
    // use it as a fallback. (If it doesn't, the 'bounds' range will stay empty, so an error will
    // be reported).
    bounds = std::equal_range(varValues.begin() + range.begin(),
                              varValues.begin() + range.end(),
                              util::missingValue(obVal));
  }

  range.constrain(static_cast<int>(bounds.first - varValues.begin()),
                  static_cast<int>(bounds.second - varValues.begin()));

  if (range.begin() == range.end()) {
    std::stringstream msg;
    msg << "No match found for exact match extraction of value '" << obVal
        << "' of the variable '" << varName << "'";
    throw eckit::Exception(msg.str(), Here());
  }
  oops::Log::debug() << "Exact match; name: " << varName << " range: " <<
    range.begin() << "," << range.end() << std::endl;
}


/// \brief Update our extract constraint based on a nearest match against the specified
/// coordinate indexing a dimension of the payload array.
/// \details
///
/// Method:
/// - Find **first** discovered nearest value in our loop.
/// - Determine which indices match this nearest value.
///   (more than one index could have this one value).
///
///   [1, 1, 2, 3, 4, 5]
///
/// Nearest neighbour extraction of “1”, has more than one neighbour.
/// That is, more than one index with the same value have the same distance:
///
///   [1, 1]  i.e. range=(0, 2)
///
/// - Note that an alternative implementation could consider equidistant
///   values, though it was decided this was not desirable behaviour:
///
///   [1, 1, 2, 3, 4, 5]
///
/// Nearest neighbour extraction of “1.5” could be then considered to have 3
/// equidistant neighbours (1, 1, 2).  That is, two different values with the
/// same distance.
///
/// [1, 1, 2] i.e. range=(0, 3)
///
/// \param[in] varName
///   Name of the coordinate to match against.
/// \param[in] varValues
///   Vector of values of that coordinate.
/// \param[in] obVal
///   Value to match.
/// \param[inout] range
///   On input, the range of slices of the payload array along the dimension indexed by
///   `varValues` that matches all constraints considered so far. On output, the subrange of
///   slices matching also the current constraint.
template<typename T>
void nearestMatch(const std::string &varName,
                  const std::vector<T> &varValues,
                  const T &obVal,
                  ConstrainedRange &range) {
  // Find first index of varValues >= obVal
  int nnIndex = std::lower_bound(varValues.begin() + range.begin(),
                                 varValues.begin() + range.end(),
                                 obVal) - varValues.begin();
  if (nnIndex >= range.end()) {
    nnIndex = range.end() - 1;
  }

  // Now fetch the nearest neighbour index (lower index prioritised for different values with
  // same distance)
  T dist = std::abs(varValues[nnIndex] - obVal);
  if ((varValues[nnIndex] > obVal) && (nnIndex > range.begin()) &&
      (std::abs(varValues[nnIndex - 1] - obVal) <= dist))
    nnIndex--;

  // Now find **same value** equidistant neighbours
  auto bounds = std::equal_range(varValues.begin() + range.begin(),
                                 varValues.begin() + range.end(),
                                 varValues[nnIndex]);
  range.constrain(static_cast<int>(bounds.first - varValues.begin()),
                  static_cast<int>(bounds.second - varValues.begin()));
  oops::Log::debug() << "Nearest match; name: " << varName << " range: " <<
    range.begin() << "," << range.end() << std::endl;
}


void nearestMatch(const std::string &varName,
                  const std::vector<std::string> &varValues,
                  const std::string &obVal,
                  ConstrainedRange &range) {
  throw eckit::UserError("The 'nearest' method cannot be used for string variables.", Here());
}


/// \brief Update our extract constraint based on a least-upper-bound match against the specified
/// coordinate indexing a dimension of the payload array.
///
/// \param[in] varName
///   Name of the coordinate to match against.
/// \param[in] varValues
///   Vector of values of that coordinate.
/// \param[in] obVal
///   Value to match.
/// \param[inout] range
///   On input, the range of slices of the payload array along the dimension indexed by
///   `varValues` that matches all constraints considered so far. On output, the subrange of
///   slices matching also the current constraint.
template<typename T>
void leastUpperBoundMatch(const std::string &varName,
                          const std::vector<T> &varValues,
                          const T &obVal,
                          ConstrainedRange &range) {
  // Find index of the first varValues >= obVal
  typedef typename std::vector<T>::const_iterator It;
  const It rangeBegin(varValues.begin() + range.begin());
  const It rangeEnd(varValues.begin() + range.end());

  const It leastUpperBoundIt = std::lower_bound(rangeBegin, rangeEnd, obVal);
  if (leastUpperBoundIt == rangeEnd) {
    std::stringstream msg;
    msg << "No match found for 'least upper bound' extraction of value '" << obVal
        << "' of the variable '" << varName << "'";
    throw eckit::Exception(msg.str(), Here());
  }

  // Find the range of items with the same value of this coordinate
  const auto bounds = std::equal_range(rangeBegin, rangeEnd, *leastUpperBoundIt);
  range.constrain(static_cast<int>(bounds.first - varValues.begin()),
                  static_cast<int>(bounds.second - varValues.begin()));
  oops::Log::debug() << "Least upper bound match; name: " << varName << " range: "
                     << range.begin() << "," << range.end() << std::endl;
}

void leastUpperBoundMatch(const std::string &varName,
                          const std::vector<std::string> &varValues,
                          const std::string &obVal,
                          ConstrainedRange &range) {
  throw eckit::UserError("The 'least upper bound' method cannot be used for string variables.",
                         Here());
}

/// \brief Update our extract constraint based on a greatest-lower-bound match against the
/// specified coordinate indexing a dimension of the payload array.
///
/// \param[in] varName
///   Name of the coordinate to match against.
/// \param[in] varValues
///   Vector of values of that coordinate.
/// \param[in] obVal
///   Value to match.
/// \param[inout] range
///   On input, the range of slices of the payload array along the dimension indexed by
///   `varValues` that matches all constraints considered so far. On output, the subrange of
///   slices matching also the current constraint.
template<typename T>
void greatestLowerBoundMatch(const std::string &varName,
                             const std::vector<T> &varValues,
                             const T &obVal,
                             ConstrainedRange &range) {
  // Find index of the last varValues <= obVal
  typedef typename std::vector<T>::const_reverse_iterator ReverseIt;
  typedef std::greater<T> Compare;
  const ReverseIt reverseRangeBegin(varValues.begin() + range.end());
  const ReverseIt reverseRangeEnd(varValues.begin() + range.begin());

  const ReverseIt greatestLowerBoundIt =
      std::lower_bound(reverseRangeBegin, reverseRangeEnd, obVal, Compare());
  if (greatestLowerBoundIt == reverseRangeEnd) {
    std::stringstream msg;
    msg << "No match found for 'greatest lower bound' extraction of value '" << obVal
        << "' of the variable '" << varName << "'";
    throw eckit::Exception(msg.str(), Here());
  }

  // Find the range of items with the same value of this coordinate
  const auto bounds = std::equal_range(varValues.begin() + range.begin(),
                                       varValues.begin() + range.end(),
                                       *greatestLowerBoundIt);
  range.constrain(static_cast<int>(bounds.first - varValues.begin()),
                  static_cast<int>(bounds.second - varValues.begin()));
  oops::Log::debug() << "Greatest lower bound match; name: " << varName << " range: "
                     << range.begin() << "," << range.end() << std::endl;
}

void greatestLowerBoundMatch(const std::string &varName,
                             const std::vector<std::string> &varValues,
                             const std::string &obVal,
                             ConstrainedRange &range) {
  throw eckit::UserError("The 'greatest lower bound' method cannot be used for string variables.",
                         Here());
}


/// \brief Restrict `range` to the subrange of `varValues` matching `obVal` according to the
/// criteria of `method`.
///
/// \param[in] method
///   Matching method.
/// \param[in] varName
///   Name of the coordinate to match against.
/// \param[in] varValues
///   Vector of values of that coordinate.
/// \param[in] obVal
///   Value to match.
/// \param[inout] range
///   On input, the range of slices of the payload array along the dimension indexed by
///   `varValues` that matches all constraints considered so far. On output, the subrange of
///   slices matching also the current constraint.
template <typename T>
void match(InterpMethod method,
           const std::string &varName,
           const std::vector<T> &varValues,
           const T &obVal,
           ConstrainedRange &range) {
  switch (method) {
    case InterpMethod::EXACT:
      exactMatch(varName, varValues, obVal, range);
      break;
    case InterpMethod::NEAREST:
      nearestMatch(varName, varValues, obVal, range);
      break;
    case InterpMethod::LEAST_UPPER_BOUND:
      leastUpperBoundMatch(varName, varValues, obVal, range);
      break;
    case InterpMethod::GREATEST_LOWER_BOUND:
      greatestLowerBoundMatch(varName, varValues, obVal, range);
      break;
    default:
      throw eckit::BadParameter("Unrecognized interpolation method", Here());
  }
}


/// \brief Perform piecewise linear interpolation along axis `dimIndex` at the 'location' `obVal`.
///
/// \details It is assumed that previous calls to extract() have extracted a single 1D slice of
/// the interpolated array parallel to the axis `varName`. This function returns the value
/// produced by piecewise linear interpolation of this slice at the point `obVal`.
///
/// \param[in] dimIndex
///   Axis of `interpolatedArray` along which to interpolate.
/// \param[in] ranges
///   Ranges of slices of `interpolatedArray` matching all constraints considered so far.
/// \param[in] varName
///   Name of the coordinate along which to interpolate.
/// \param[in] varValues
///   Vector of values of that coordinate.
/// \param[in] obVal
///   Interpolation location.
/// \param[in] interpolatedArray
///   Interpolated array.
template <typename CoordinateValue>
float linearInterpolation(size_t dimIndex,
                          const std::array<ConstrainedRange, 2> &ranges,
                          const std::string &varName,
                          const std::vector<CoordinateValue> &varValues,
                          const CoordinateValue &obVal,
                          const boost::multi_array<float, 2> &interpolatedArray) {
  // Sanity check constraint
  int sizeDim0 = ranges[0].size();
  int sizeDim1 = ranges[1].size();
  if ((dimIndex == 1 && !(sizeDim1 > 1 && sizeDim0 == 1)) ||
      (dimIndex == 0 && !(sizeDim0 > 1 && sizeDim1 == 1))) {
    throw eckit::Exception("Linear interpolation failed - data must be 1D.", Here());
  }

  // Constrain our index range in the relevant dimension.
  const ConstrainedRange &range = ranges[static_cast<size_t>(dimIndex)];

  if ((obVal > varValues[range.end() - 1]) || (obVal < varValues[range.begin()])) {
    throw eckit::Exception("Linear interpolation failed, value is beyond grid extent."
                           "No extrapolation supported.",
                           Here());
  }
  // Find first index of varValues >= obVal
  int nnIndex = std::lower_bound(varValues.begin() + range.begin(),
                                 varValues.begin() + range.end(),
                                 obVal) - varValues.begin();

  // Determine upper or lower indices from this
  if (varValues[nnIndex] == obVal) {
    // No interpolation required (is equal)
    float res;
    if (dimIndex == 1) {
      res = interpolatedArray[ranges[0].begin()][nnIndex];
    } else {
      res = interpolatedArray[nnIndex][ranges[1].begin()];
    }
    return res;
  }
  // Linearly interpolate between these two indices.
  float zUpper;
  float zLower;
  if (dimIndex == 1) {
    zUpper = interpolatedArray[ranges[0].begin()][nnIndex];
    zLower = interpolatedArray[ranges[0].begin()][nnIndex-1];
  } else {
    zUpper = interpolatedArray[nnIndex][ranges[1].begin()];
    zLower = interpolatedArray[nnIndex-1][ranges[1].begin()];
  }
  float res = ((static_cast<float>(obVal - varValues[nnIndex-1]) /
                static_cast<float>(varValues[nnIndex] - varValues[nnIndex-1])) *
               (zUpper - zLower)) + zLower;
  return res;
}


float linearInterpolation(size_t dimIndex,
                          const std::array<ConstrainedRange, 2> &ranges,
                          const std::string &varName,
                          const std::vector<std::string> &varValues,
                          const std::string &obVal,
                          const boost::multi_array<float, 2> &interpolatedArray) {
  throw eckit::UserError("Linear interpolation cannot be performed along coordinate axes indexed "
                         "by string variables such as " + varName + ".", Here());
}

}  // namespace


template <typename ExtractedValue>
DataExtractor<ExtractedValue>::DataExtractor(const std::string &filepath,
                                             const std::string &group) {
  // Read the data from the file
  load(filepath, group);
  // Start by constraining to the full range of our data
  resetExtract();
  // Initialise splitter for both dimensions
  splitter_.emplace_back(ufo::RecursiveSplitter(interpolatedArray2D_.shape()[0]));
  splitter_.emplace_back(ufo::RecursiveSplitter(interpolatedArray2D_.shape()[1]));
}


template <typename ExtractedValue>
void DataExtractor<ExtractedValue>::load(const std::string &filepath,
                                         const std::string &interpolatedArrayGroup) {
  std::unique_ptr<DataExtractorBackend<ExtractedValue>> backend = createBackendFor(filepath);
  DataExtractorInput<ExtractedValue> input = backend->loadData(interpolatedArrayGroup);
  coord2DimMapping_ = std::move(input.coord2DimMapping);
  dim2CoordMapping_ = std::move(input.dim2CoordMapping);
  coordsVals_ = std::move(input.coordsVals);
  interpolatedArray2D_.resize(boost::extents[input.payloadArray.shape()[0]]
                                            [input.payloadArray.shape()[1]]);
  interpolatedArray2D_ = std::move(input.payloadArray);
  // Set the unconstrained size of matching ranges along both axes of the payload array.
  for (size_t i = 0; i < constrainedRanges_.size(); ++i)
    constrainedRanges_[i] = ConstrainedRange(input.payloadArray.shape()[i]);
}


template <typename ExtractedValue>
std::unique_ptr<DataExtractorBackend<ExtractedValue>>
DataExtractor<ExtractedValue>::createBackendFor(const std::string &filepath) {
  const std::string lowercasePath = eckit::StringTools::lower(filepath);
  if (eckit::StringTools::endsWith(lowercasePath, ".nc") ||
      eckit::StringTools::endsWith(lowercasePath, ".nc4"))
    return boost::make_unique<DataExtractorNetCDFBackend<ExtractedValue>>(filepath);
  else if (eckit::StringTools::endsWith(lowercasePath, ".csv"))
    return boost::make_unique<DataExtractorCSVBackend<ExtractedValue>>(filepath);
  else
    throw eckit::BadValue("File '" + filepath + "' has an unrecognized extension. "
                          "The supported extensions are .nc, .nc4 and .csv", Here());
}


template <typename ExtractedValue>
void DataExtractor<ExtractedValue>::sort() {
  boost::multi_array<ExtractedValue, 2> sortedArray = interpolatedArray2D_;
  nextCoordToExtractBy_ = coordsToExtractBy_.begin();

  for (size_t dim = 0; dim < dim2CoordMapping_.size(); ++dim) {
    // Reorder coordinates
    for (auto &coord : dim2CoordMapping_[dim]) {
      auto &coordVal = coordsVals_[coord];
      SortVisitor visitor(splitter_[dim]);
      boost::apply_visitor(visitor, coordVal);
    }
    // Reorder the array to be interpolated
    if (dim == 0) {
      int ind = -1;
      for (const auto &group : splitter_[dim].groups()) {
        for (const auto &index : group) {
          ind++;
          oops::Log::debug() << "Sort index dim0; index-from: " << ind << " index-to: " <<
            index << std::endl;
          for (size_t j = 0; j < interpolatedArray2D_.shape()[1]; j++) {
            sortedArray[ind][j] = interpolatedArray2D_[index][j];
          }
        }
      }
      // Replace the unsorted array with the sorted one.
      interpolatedArray2D_ = sortedArray;
    } else if (dim == 1) {
      int ind = -1;
      for (const auto &group : splitter_[dim].groups()) {
        for (const auto &index : group) {
          ind++;
          oops::Log::debug() << "Sort index dim1; index-from: " << ind << " index-to: " <<
            index << std::endl;
          for (size_t i = 0; i < interpolatedArray2D_.shape()[0]; i++) {
            sortedArray[i][ind] = interpolatedArray2D_[i][index];
          }
        }
      }
      // Replace the unsorted array with the sorted one.
      interpolatedArray2D_ = sortedArray;
    } else {
        throw eckit::Exception("Unable to reorder the array to be interpolated: "
                               "it has more than 2 dimensions.", Here());
    }
  }
}


template <typename ExtractedValue>
void DataExtractor<ExtractedValue>::scheduleSort(const std::string &varName,
                                                 const InterpMethod &method) {
  if (method == InterpMethod::LINEAR && !std::is_floating_point<ExtractedValue>::value)
    throw eckit::BadParameter("Linear interpolation can be used when extracting floating-point "
                              "values, but not integers or strings.", Here());

  // Map any names of the form var@Group to Group/var
  const std::string canonicalVarName = ioda::convertV1PathToV2Path(varName);

  const CoordinateValues &coordVal = coordsVals_.at(canonicalVarName);
  const int dimIndex = coord2DimMapping_.at(canonicalVarName);

  SortUpdateVisitor visitor(splitter_[static_cast<size_t>(dimIndex)]);
  boost::apply_visitor(visitor, coordVal);

  // Update our map between coordinate (variable) and interpolation/extract method
  coordsToExtractBy_.emplace_back(Coordinate{varName, coordVal, method, dimIndex});
}


template <typename ExtractedValue>
void DataExtractor<ExtractedValue>::extract(float obVal) {
  extractImpl(obVal);
}


template <typename ExtractedValue>
void DataExtractor<ExtractedValue>::extract(int obVal) {
  extractImpl(obVal);
}


template <typename ExtractedValue>
void DataExtractor<ExtractedValue>::extract(const std::string &obVal) {
  extractImpl(obVal);
}


template <typename ExtractedValue>
template <typename T>
void DataExtractor<ExtractedValue>::extractImpl(const T &obVal) {
  if (nextCoordToExtractBy_ == coordsToExtractBy_.cend()) {
    throw eckit::UserError("Too many extract() calls made for the expected number of variables.",
                           Here());
  }

  // Perform the extraction using the selected method
  if (nextCoordToExtractBy_->method == InterpMethod::LINEAR)
    maybeExtractByLinearInterpolation(obVal);
  else
    match(nextCoordToExtractBy_->method,
          nextCoordToExtractBy_->name,
          boost::get<std::vector<T>>(nextCoordToExtractBy_->values),
          obVal,
          constrainedRanges_[nextCoordToExtractBy_->payloadDim]);

  ++nextCoordToExtractBy_;
}


// Primary template, used for all ExtractedValue types except float.
template <typename ExtractedValue>
template <typename T>
void DataExtractor<ExtractedValue>::maybeExtractByLinearInterpolation(const T &obVal) {
  // Should never be called -- this error should be detected earlier.
  throw eckit::BadParameter("Linear interpolation can be used when extracting floating-point "
                            "values, but not integers or strings.", Here());
}


// Specialization for ExtractedValue = float.
template <>
template <typename T>
void DataExtractor<float>::maybeExtractByLinearInterpolation(const T &obVal) {
  result_ = linearInterpolation(nextCoordToExtractBy_->payloadDim,
                                constrainedRanges_,
                                nextCoordToExtractBy_->name,
                                boost::get<std::vector<T>>(nextCoordToExtractBy_->values),
                                obVal,
                                interpolatedArray2D_);
  resultSet_ = true;
}


// Primary template, used for all ExtractedValue types except float.
template <typename ExtractedValue>
ExtractedValue DataExtractor<ExtractedValue>::getResult() {
  // Fetch the result
  ExtractedValue res = getUniqueMatch();
  resetExtract();
  return res;
}


// Specialization adding support for linear interpolation.
template <>
float DataExtractor<float>::getResult() {
  // Fetch the result
  if (resultSet_) {
    // This was derived from linear interpolation so return it.
    resetExtract();
    return result_;
  }

  float res = getUniqueMatch();
  resetExtract();
  return res;
}


template <typename ExtractedValue>
ExtractedValue DataExtractor<ExtractedValue>::getUniqueMatch() const {
  // This function should be called only if linear interpolation is not used within the
  // extraction process.
  ASSERT(!resultSet_);

  int sizeDim0 = constrainedRanges_[0].size();
  int sizeDim1 = constrainedRanges_[1].size();
  if (sizeDim0 != 1 || sizeDim1 != 1) {
    throw eckit::Exception("Previous calls to extract() have failed to identify "
                           "a single value to return.", Here());
  }
  return interpolatedArray2D_[constrainedRanges_[0].begin()]
                             [constrainedRanges_[1].begin()];
}


template <typename ExtractedValue>
void DataExtractor<ExtractedValue>::resetExtract() {
  for (ConstrainedRange &range : constrainedRanges_)
    range.reset();
  resultSet_ = false;
  nextCoordToExtractBy_ = coordsToExtractBy_.begin();
}


// Explicit instantiations
template class DataExtractor<float>;
template class DataExtractor<int>;
template class DataExtractor<std::string>;

}  // namespace ufo
