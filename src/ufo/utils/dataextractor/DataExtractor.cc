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

#include "eckit/utils/StringTools.h"

#include "ioda/Misc/StringFuncs.h"

#include "oops/util/Logger.h"

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
    splitter.sortGroupsBy([&coord](int index) { return coord[static_cast<size_t>(index)]; });
  }

  ufo::RecursiveSplitter &splitter;
};

/// \brief Boost visitor which allows us to apply a coordinate transformation.
class CoordinateTransformationVisitor : public boost::static_visitor<void> {
 public:
  explicit CoordinateTransformationVisitor(const CoordinateTransformation &coordTrans,
                                           const std::string &varName)
    : coordTrans_(coordTrans), varName_(varName) {}

  // Do not transform any non-float coordinates
  template <typename T>
  void operator()(std::vector<T> &coord) {
  }

  void operator()(std::vector<float> &coord) {
    if (coordTrans_ == CoordinateTransformation::LOGLINEAR)
      applyLogLinearTransform(varName_, coord);
  }

  const CoordinateTransformation &coordTrans_;
  const std::string &varName_;
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
                              util::missingValue<T>());
  }

  range.constrain(static_cast<int>(bounds.first - varValues.begin()),
                  static_cast<int>(bounds.second - varValues.begin()));

  if (range.begin() == range.end()) {
    std::stringstream msg;
    msg << "No match found for exact match extraction of value '" << obVal
        << "' of the variable '" << varName << "'";
    throw eckit::Exception(msg.str(), Here());
  }
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
/// \ param[in] equidistantChoice
///   Choice of index (first or last) to use if there is a value which is equidistant between
///   two lookup values.
/// \param[inout] range
///   On input, the range of slices of the payload array along the dimension indexed by
///   `varValues` that matches all constraints considered so far. On output, the subrange of
///   slices matching also the current constraint.
template<typename T>
void nearestMatch(const std::string &varName,
                  const std::vector<T> &varValues,
                  const T &obVal,
                  const EquidistantChoice &equidistantChoice,
                  ConstrainedRange &range) {
  if (isOutOfBounds(obVal, varValues, range)) {
    std::stringstream msg;
    msg << "No match found for 'nearest' extraction of value '" << obVal << "' of the variable '"
        << varName << "'.  Value is out of bounds.  Consider using extrapolation.";
    throw eckit::Exception(msg.str(), Here());
  }

  // Find first index of varValues >= obVal
  int nnIndex = std::lower_bound(varValues.begin() + range.begin(),
                                 varValues.begin() + range.end(),
                                 obVal) - varValues.begin();
  if (nnIndex >= range.end())
    nnIndex = range.end() - 1;

  // Now fetch the nearest neighbour index.  In the case of an equidistance result then
  // either the first or last index is chosen depending on the input.
  T dist = std::abs(varValues[nnIndex] - obVal);
  if (equidistantChoice == EquidistantChoice::FIRST) {
      if ((varValues[nnIndex] > obVal) && (nnIndex > range.begin()) &&
          (std::abs(varValues[nnIndex - 1] - obVal) <= dist))
      nnIndex--;
  } else if (equidistantChoice == EquidistantChoice::LAST) {
      if ((varValues[nnIndex] > obVal) && (nnIndex > range.begin()) &&
          (std::abs(varValues[nnIndex - 1] - obVal) < dist))
      nnIndex--;
  }

  // Now find **same value** equidistant neighbours
  auto bounds = std::equal_range(varValues.begin() + range.begin(),
                                 varValues.begin() + range.end(),
                                 varValues[nnIndex]);
  range.constrain(static_cast<int>(bounds.first - varValues.begin()),
                  static_cast<int>(bounds.second - varValues.begin()));
}


void nearestMatch(const std::string &varName,
                  const std::vector<std::string> &varValues,
                  const std::string &obVal,
                  const EquidistantChoice &equidstantChoice,
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
        << "' of the variable '" << varName << "'.  Value is out of bounds.  Consider using "
           "extrapolation.";
    throw eckit::Exception(msg.str(), Here());
  }

  // Find the range of items with the same value of this coordinate
  const auto bounds = std::equal_range(rangeBegin, rangeEnd, *leastUpperBoundIt);
  range.constrain(static_cast<int>(bounds.first - varValues.begin()),
                  static_cast<int>(bounds.second - varValues.begin()));
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
        << "' of the variable '" << varName << "'.  Value is out of bounds.  Consider using "
           "extrapolation.";
    throw eckit::Exception(msg.str(), Here());
  }

  // Find the range of items with the same value of this coordinate
  const auto bounds = std::equal_range(varValues.begin() + range.begin(),
                                       varValues.begin() + range.end(),
                                       *greatestLowerBoundIt);
  range.constrain(static_cast<int>(bounds.first - varValues.begin()),
                  static_cast<int>(bounds.second - varValues.begin()));
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
/// \ param[in] equidistantChoice
///   Value to select whether first or last index is used when a value is equidistant to
///   two lookup indices with the nearest method.
/// \param[inout] range
///   On input, the range of slices of the payload array along the dimension indexed by
///   `varValues` that matches all constraints considered so far. On output, the subrange of
///   slices matching also the current constraint.
template <typename T>
void match(const InterpMethod method,
           const std::string &varName,
           const std::vector<T> &varValues,
           const T &obVal,
           const EquidistantChoice equidistantChoice,
           const CoordinateTransformation coordTransformation,
           ConstrainedRange &range) {
  switch (method) {
    case InterpMethod::EXACT:
      exactMatch(varName, varValues, obVal, range);
      break;
    case InterpMethod::NEAREST:
      nearestMatch(varName, varValues, obVal, equidistantChoice, range);
      break;
    case InterpMethod::LEAST_UPPER_BOUND:
      leastUpperBoundMatch(varName, varValues, obVal, range);
      break;
    case InterpMethod::GREATEST_LOWER_BOUND:
      greatestLowerBoundMatch(varName, varValues, obVal, range);
      break;
    default:
      throw eckit::BadParameter("Unrecognized interpolation method for '" + varName + "'", Here());
  }
}


/// \brief Perform piecewise linear interpolation of the provided array `varValues` at 'location'
/// `obVal`.
///
/// \details It is assumed that the provided 1D array is described by coordinate `varName`.
/// This function returns the value produced by piecewise linear interpolation of this array at
/// the point `obVal`.
///
/// \param[in] range
///   Defines how to constrain (slice) `varValues` along with `interpolatedArray`.
/// \param[in] varName
///   Name of the coordinate along which to interpolate.
/// \param[in] varValues
///   Vector of values of that coordinate.
/// \param[in] obVal
///   Interpolation location.
/// \param[in] interpolatedArray
///   Interpolated array.
template <typename CoordinateValue>
float linearInterpolation(
    const std::string &varName,
    const std::vector<CoordinateValue> &varValues,
    const CoordinateValue &obVal,
    const ConstrainedRange &range,
    const DataExtractorPayload<float>::const_array_view<1>::type &interpolatedArray) {

  if (obVal == util::missingValue<CoordinateValue>()) {
    const float missing = util::missingValue<float>();
    return missing;
  }

  if (isOutOfBounds(obVal, varValues, range)) {
      std::stringstream msg;
      msg << "No match found for 'linear' interpolation of value '" << obVal
          << "' of the variable '" << varName << "'.  Value is out of bounds.  Consider using "
          << "extrapolation.";
      throw eckit::Exception(msg.str(), Here());
  }
  // Find first index of varValues >= obVal
  int nnIndex = std::lower_bound(varValues.begin() + range.begin(),
                                 varValues.begin() + range.end(),
                                 obVal) - varValues.begin();

  // No interpolation required (is equal)
  if (varValues[nnIndex] == obVal)
    return interpolatedArray[nnIndex];

  // Linearly interpolate between these two indices.
  const float zUpper = interpolatedArray[nnIndex];
  const float zLower = interpolatedArray[nnIndex-1];
  float res = ((static_cast<float>(obVal - varValues[nnIndex-1]) /
                static_cast<float>(varValues[nnIndex] - varValues[nnIndex-1])) *
               (zUpper - zLower)) + zLower;
  return res;
}


float linearInterpolation(
    const std::string &varName,
    const std::vector<std::string> &varValues,
    const std::string &obVal,
    const ConstrainedRange &range,
    const DataExtractorPayload<float>::const_array_view<1>::type &interpolatedArray) {
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
  // Initialise splitter for each dimension
  splitter_.emplace_back(ufo::RecursiveSplitter(interpolatedArray_.shape()[0]));
  splitter_.emplace_back(ufo::RecursiveSplitter(interpolatedArray_.shape()[1]));
  splitter_.emplace_back(ufo::RecursiveSplitter(interpolatedArray_.shape()[2]));
}


template <typename ExtractedValue>
void DataExtractor<ExtractedValue>::load(const std::string &filepath,
                                         const std::string &interpolatedArrayGroup) {
  std::unique_ptr<DataExtractorBackend<ExtractedValue>> backend = createBackendFor(filepath);
  DataExtractorInput<ExtractedValue> input = backend->loadData(interpolatedArrayGroup);
  coord2DimMapping_ = std::move(input.coord2DimMapping);
  dim2CoordMapping_ = std::move(input.dim2CoordMapping);
  coordNDims_ = std::move(input.coordNDims);
  coordsVals_ = std::move(input.coordsVals);
  interpolatedArray_.resize(boost::extents[input.payloadArray.shape()[0]]
                                            [input.payloadArray.shape()[1]]
                                            [input.payloadArray.shape()[2]]);
  interpolatedArray_ = std::move(input.payloadArray);
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
  DataExtractorPayload<ExtractedValue> sortedArray = interpolatedArray_;
  nextCoordToExtractBy_ = coordsToExtractBy_.begin();

  for (size_t dim = 0; dim < dim2CoordMapping_.size(); ++dim) {
    if (interpolatedArray_.shape()[dim] == 1)  // Avoid sorting scalar coordinates
      continue;

    // Reorder coordinates
    for (auto &coord : dim2CoordMapping_[dim]) {
      auto &coordVal = coordsVals_[coord];
      SortVisitor visitor(splitter_[dim]);
      boost::apply_visitor(visitor, coordVal);
    }

    // Reorder the array to be interpolated
    int ind = 0;
    std::array<size_t, 2> otherDims;
    for (size_t odim = 0; odim < interpolatedArray_.dimensionality; ++odim) {
      if (odim != dim) {
        otherDims[ind] = odim;
        ind++;
      }
    }

    ind = 0;
    for (const auto &group : splitter_[dim].groups()) {
      for (const auto &index : group) {
        for (size_t j = 0; j < interpolatedArray_.shape()[otherDims[0]]; j++) {
          for (size_t k = 0; k < interpolatedArray_.shape()[otherDims[1]]; k++) {
            if (dim == 0) {
              sortedArray[ind][j][k] = interpolatedArray_[index][j][k];
            } else if (dim == 1) {
              sortedArray[j][ind][k] = interpolatedArray_[j][index][k];
            } else if (dim == 2) {
              sortedArray[j][k][ind] = interpolatedArray_[j][k][index];
            } else {
              // We shouldn't ever end up here (exception should be thrown eariler).
              throw eckit::Exception("Unable to reorder the array to be interpolated: "
                                     "it has more than 3 dimensions.", Here());
            }
          }
        }
        ind++;
      }
    }
    // Replace the unsorted array with the sorted one.
    interpolatedArray_ = sortedArray;
  }
}


template <typename ExtractedValue>
void DataExtractor<ExtractedValue>::scheduleSort
    (const std::string &varName,
     const InterpMethod &method,
     const ExtrapolationMode &extrapMode,
     const EquidistantChoice &equidistantChoice,
     const CoordinateTransformation &coordinateTransformation) {
  if (!std::is_floating_point<ExtractedValue>::value) {
      std::string msg = "interpolation can be used when extracting floating-point values, but not "
                        "integers or strings.";
      if (method == InterpMethod::LINEAR) {
        throw eckit::BadParameter("Linear " + msg, Here());
      } else if (method == InterpMethod::BILINEAR) {
        throw eckit::BadParameter("Bilinear " + msg, Here());
      } else if (method == InterpMethod::TRILINEAR) {
        throw eckit::BadParameter("Trilinear " + msg, Here());
      }
  }

  // Map any names of the form var@Group to Group/var
  const std::string canonicalVarName = ioda::convertV1PathToV2Path(varName);

  CoordinateValues &coordVal = coordsVals_.at(canonicalVarName);
  const std::vector<int> &dimIndices = coord2DimMapping_.at(canonicalVarName);
  const size_t coordDim = coordNDims_.at(canonicalVarName);

  if (coordDim != dimIndices.size())
    throw eckit::Exception("Variable: '" + varName + "' has one or more dimension mappings not "
                           "shared by the payload variable.", Here());
  if (coordDim > 1) {
    std::stringstream msg;
    msg << "Variable: '" + varName + "' is a '" << coordDim << + "' dimensional coordinate."
        << "  Only 1D coordinates currently supported.";
    throw eckit::Exception(msg.str(), Here());
  }

  const int dimIndex = dimIndices[0];

  SortUpdateVisitor visitor(splitter_[static_cast<size_t>(dimIndex)]);
  boost::apply_visitor(visitor, coordVal);

  // Apply coordinate transformations.
  CoordinateTransformationVisitor coordTransVisitor(coordinateTransformation, varName);
  boost::apply_visitor(coordTransVisitor, coordVal);

  // Update our map between coordinate (variable) and interpolation/extract method
  coordsToExtractBy_.emplace_back(Coordinate{varName, coordVal, method, extrapMode,
        equidistantChoice, coordinateTransformation, dimIndex});
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
  if (nextCoordToExtractBy_ == coordsToExtractBy_.cend())
    throw eckit::UserError("Too many extract() calls made for the expected number of variables.",
                           Here());
  T obValN = applyExtrapolation(obVal);
  if (resultSet_)
    return;

  // Perform the extraction using the selected method
  if (nextCoordToExtractBy_->method == InterpMethod::LINEAR)
    maybeExtractByLinearInterpolation(obValN);
  else
    match(nextCoordToExtractBy_->method, nextCoordToExtractBy_->name,
          boost::get<std::vector<T>>(nextCoordToExtractBy_->values), obValN,
          nextCoordToExtractBy_->equidistantChoice,
          nextCoordToExtractBy_->coordinateTransformation,
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
  int dimIndex = nextCoordToExtractBy_->payloadDim;
  const auto &interpolatedArray = get1DSlice(interpolatedArray_,
                                             dimIndex,
                                             constrainedRanges_);
  result_ = linearInterpolation(nextCoordToExtractBy_->name,
                                boost::get<std::vector<T>>(nextCoordToExtractBy_->values),
                                obVal, constrainedRanges_[dimIndex], interpolatedArray);
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
    // This was derived from linear/bilinear/trilinear interpolation so return it.
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

  for (size_t dim=0; dim < constrainedRanges_.size(); dim++) {
    if (constrainedRanges_[dim].size() == 0)
      throw eckit::Exception("No match found in the interpolation array.", Here());
    else if (constrainedRanges_[dim].size() > 1)
      throw eckit::Exception("Extraction criteria were not sufficient to identify a unique match "
                             "in the interpolation array.", Here());
  }
  return interpolatedArray_[constrainedRanges_[0].begin()]
                           [constrainedRanges_[1].begin()]
                           [constrainedRanges_[2].begin()];
}


template <typename ExtractedValue>
void DataExtractor<ExtractedValue>::resetExtract() {
  for (ConstrainedRange &range : constrainedRanges_)
    range.reset();
  resultSet_ = false;
  nextCoordToExtractBy_ = coordsToExtractBy_.begin();
}


void applyLogLinearTransform(const std::string &varName,
                             std::vector<float> &varValues) {
  const float missing = util::missingValue<float>();
  if (std::any_of(varValues.cbegin(), varValues.cend(),
                  [missing](float c){return c != missing && c < 0.0f;})) {
    std::stringstream msg;
    msg << "Cannot perform log-linear interpolation for variable '" << varName << "'. "
        << "At least one value is negative.";
    throw eckit::Exception(msg.str(), Here());
  }

  // Transform data.
  std::transform(varValues.cbegin(), varValues.cend(), varValues.begin(),
                 [missing](float c){return c != missing ? std::log(c) : missing;});
}

void applyLogLinearTransform(const std::string &varName,
                             float &varValue) {
  std::vector<float> vecTemp{varValue};
  applyLogLinearTransform(varName, vecTemp);
  varValue = vecTemp[0];
}

float trilinearInterpolation(
    const std::string &varName0,
    const std::vector<float> &varValues0,
    float obVal0,
    const ConstrainedRange &range0,
    const CoordinateTransformation &coordTrans0,
    const std::string &varName1,
    const std::vector<float> &varValues1,
    float obVal1,
    const ConstrainedRange &range1,
    const CoordinateTransformation &coordTrans1,
    const std::string &varName2,
    const std::vector<float> &varValues2,
    float obVal2,
    const ConstrainedRange &range2,
    const CoordinateTransformation &coordTrans2,
    const DataExtractorPayload<float>::const_array_view<3>::type &interpolatedArray)
{
  const float missing = util::missingValue<float>();

  if (obVal0 == missing || obVal1 == missing || obVal2 == missing) {
    return missing;
  }

  if (isOutOfBounds(obVal0, varValues0, range0)) {
      std::stringstream msg;
      msg << "No match found for 'trilinear' interpolation of value '" << obVal0
          << "' of the variable '" << varName0 << "'.  Value is out of bounds.  Consider using "
          << "extrapolation.";
      throw eckit::Exception(msg.str(), Here());
  }
  if (isOutOfBounds(obVal1, varValues1, range1)) {
      std::stringstream msg;
      msg << "No match found for 'trilinear' interpolation of value '" << obVal1
          << "' of the variable '" << varName1 << "'.  Value is out of bounds.  Consider using "
          << "extrapolation.";
      throw eckit::Exception(msg.str(), Here());
  }
  if (isOutOfBounds(obVal2, varValues2, range2)) {
      std::stringstream msg;
      msg << "No match found for 'trilinear' interpolation of value '" << obVal2
          << "' of the variable '" << varName2 << "'.  Value is out of bounds.  Consider using "
          << "extrapolation.";
      throw eckit::Exception(msg.str(), Here());
  }

  const int nnIndex0 =
    std::lower_bound(varValues0.begin() + range0.begin(),
                     varValues0.begin() + range0.end(), obVal0) - varValues0.begin();
  const int nnIndex1 =
    std::lower_bound(varValues1.begin() + range1.begin(),
                     varValues1.begin() + range1.end(), obVal1) - varValues1.begin();
  const int nnIndex2 =
    std::lower_bound(varValues2.begin() + range2.begin(),
                     varValues2.begin() + range2.end(), obVal2) - varValues2.begin();

  if ((varValues0[nnIndex0] == obVal0) &&
      (varValues1[nnIndex1] == obVal1) &&
      (varValues2[nnIndex2] == obVal2)) {
    // No interpolation required.
    return interpolatedArray[nnIndex0][nnIndex1][nnIndex2];
  }

  // Setup points
  // ------------
  // Indices
  const int ix1 = nnIndex0-1 >= 0 ? nnIndex0-1 : nnIndex0;
  const int ix2 = ix1 + 1;
  const int iy1 = nnIndex1-1 >= 0 ? nnIndex1-1 : nnIndex1;
  const int iy2 = iy1 + 1;
  const int iz1 = nnIndex2-1 >= 0 ? nnIndex2-1 : nnIndex2;
  const int iz2 = iz1 + 1;

  // Coord locations (surrounding the coordinates of the observation)
  const float x1 = varValues0[ix1], x2 = varValues0[ix2];
  const float y1 = varValues1[iy1], y2 = varValues1[iy2];
  const float z1 = varValues2[iz1], z2 = varValues2[iz2];

  // Coordinate values at these locations
  // Variable names are of the form qxyz
  const float q111 = interpolatedArray[ix1][iy1][iz1],
    q112 = interpolatedArray[ix1][iy1][iz2],
    q121 = interpolatedArray[ix1][iy2][iz1],
    q122 = interpolatedArray[ix1][iy2][iz2],
    q211 = interpolatedArray[ix2][iy1][iz1],
    q212 = interpolatedArray[ix2][iy1][iz2],
    q221 = interpolatedArray[ix2][iy2][iz1],
    q222 = interpolatedArray[ix2][iy2][iz2];

  // Missing data handling
  // - Pick nearest non-missing neighbour of the moore neighbourhood (no diagonals) if any of these
  //   8 neighbours are missing.
  // - Return missing if all 8 of the moore neighbourhood (no diagonals) are missing.
  const int nmissing = (q111 == missing) + (q112 == missing) + (q121 == missing) +
    (q122 == missing) + (q211 == missing) + (q212 == missing) + (q221 == missing) +
    (q222 == missing);
  if (nmissing > 0) {
    if (nmissing == 8)
      return missing;
    const std::array<float, 8> xd = {x1, x1, x1, x1, x2, x2, x2, x2};
    const std::array<float, 8> yd = {y1, y1, y2, y2, y1, y1, y2, y2};
    const std::array<float, 8> zd = {z1, z2, z1, z2, z1, z2, z1, z2};
    const std::array<float, 8> qd = {q111, q112, q121, q122, q211, q212, q221, q222};

    float minq, mindist = std::numeric_limits<float>::max();
    for (size_t ind=0; ind < qd.size(); ind++) {
      // Euclidean distance between two points.
      // todo: std::hypot accepts three arguments in C++17
      const float currdist = std::sqrt((obVal0-xd[ind]) * (obVal0-xd[ind]) +
                                       (obVal1-yd[ind]) * (obVal1-yd[ind]) +
                                       (obVal2-zd[ind]) * (obVal2-zd[ind]));
      if (currdist < mindist && qd[ind] != missing) {
        mindist = currdist;
        minq = qd[ind];
      }
    }
    return minq;
  }

  // Interpolate at obVal0, obVal1 and obVal2
  // Scaled coordinate values
  const float xscaled = (obVal0 - x1) / (x2 - x1);
  const float yscaled = (obVal1 - y1) / (y2 - y1);
  const float zscaled = (obVal2 - z1) / (z2 - z1);

  // Interpolate along x (variable names are of the form qyz)
  const float q11 = q111 * (1.0 - xscaled) + q211 * xscaled;
  const float q12 = q112 * (1.0 - xscaled) + q212 * xscaled;
  const float q21 = q121 * (1.0 - xscaled) + q221 * xscaled;
  const float q22 = q122 * (1.0 - xscaled) + q222 * xscaled;

  // Interpolate along y (variable names are of the form qz)
  const float q1 = q11 * (1.0 - yscaled) + q21 * yscaled;
  const float q2 = q12 * (1.0 - yscaled) + q22 * yscaled;

  // Interpolate along z
  return q1 * (1.0 - zscaled) + q2 * zscaled;
}

template<>
void DataExtractor<float>::maybeExtractByTriLinearInterpolation(
  const float &obValDim0, const float &obValDim1, const float &obValDim2) {
  auto &ranges = constrainedRanges_;

  // Perform requested transformation and extrapolation to the value along axis 0.
  float obValDim0N = obValDim0;
  const CoordinateTransformation &coordTrans0 = nextCoordToExtractBy_->coordinateTransformation;
  if (coordTrans0 == CoordinateTransformation::LOGLINEAR)
    applyLogLinearTransform(nextCoordToExtractBy_->name, obValDim0N);
  obValDim0N = applyExtrapolation(obValDim0N);

  if (resultSet_)
    return;
  const size_t dimIndex0 = nextCoordToExtractBy_->payloadDim;
  const std::string &varName0 = nextCoordToExtractBy_->name;
  const std::vector<float> &varValues0 =
    boost::get<std::vector<float>>(nextCoordToExtractBy_->values);
  ++nextCoordToExtractBy_;  // Consume variable

  if (nextCoordToExtractBy_->method != InterpMethod::TRILINEAR)
    throw eckit::BadParameter("Second parameter provided to the Trilinear interpolator is not of "
                              "method 'trilinear'.", Here());

  // Perform requested transformation and extrapolation to the value along axis 1.
  float obValDim1N = obValDim1;
  const CoordinateTransformation &coordTrans1 = nextCoordToExtractBy_->coordinateTransformation;
  if (coordTrans1 == CoordinateTransformation::LOGLINEAR)
    applyLogLinearTransform(nextCoordToExtractBy_->name, obValDim1N);
  obValDim1N = applyExtrapolation(obValDim1N);
  if (resultSet_)
    return;
  const size_t dimIndex1 = nextCoordToExtractBy_->payloadDim;
  const std::string &varName1 = nextCoordToExtractBy_->name;
  const std::vector<float> &varValues1 =
    boost::get<std::vector<float>>(nextCoordToExtractBy_->values);
  ++nextCoordToExtractBy_;  // Consume variable

  if (nextCoordToExtractBy_->method != InterpMethod::TRILINEAR)
    throw eckit::BadParameter("Third parameter provided to the Trilinear interpolator is not of "
                              "method 'trilinear'.", Here());

  // Perform requested transformation and extrapolation to the value along axis 2.
  float obValDim2N = obValDim2;
  const CoordinateTransformation &coordTrans2 = nextCoordToExtractBy_->coordinateTransformation;
  if (coordTrans2 == CoordinateTransformation::LOGLINEAR)
    applyLogLinearTransform(nextCoordToExtractBy_->name, obValDim2N);
  obValDim2N = applyExtrapolation(obValDim2N);
  if (resultSet_)
    return;
  const size_t dimIndex2 = nextCoordToExtractBy_->payloadDim;
  const std::string &varName2 = nextCoordToExtractBy_->name;
  const std::vector<float> &varValues2 =
    boost::get<std::vector<float>>(nextCoordToExtractBy_->values);

  // convert interpolatedArray_ to the required type
  typename DataExtractorPayload<float>::index_gen indices;
  typedef boost::multi_array_types::index_range range_t;
  typedef boost::multi_array_ref<float, 3>::array_view<3>::type interparray;
  const interparray interpolatedArray =
    interpolatedArray_[indices[range_t()][range_t()][range_t()]];

  // Check the order of the variables, and make sure that they match the order in the file
  // that we are drawing from.  Hard to avoid long-and-complicated if-block.
  if (dimIndex0 == 0 && dimIndex1 == 1) {
    result_ =
      trilinearInterpolation(varName0, varValues0, obValDim0N, ranges[dimIndex0], coordTrans0,
                             varName1, varValues1, obValDim1N, ranges[dimIndex1], coordTrans1,
                             varName2, varValues2, obValDim2N, ranges[dimIndex2], coordTrans2,
                             interpolatedArray);
  } else if (dimIndex0 == 0 && dimIndex1 == 2) {
    result_ =
      trilinearInterpolation(varName0, varValues0, obValDim0N, ranges[dimIndex0], coordTrans0,
                             varName2, varValues2, obValDim2N, ranges[dimIndex2], coordTrans2,
                             varName1, varValues1, obValDim1N, ranges[dimIndex1], coordTrans1,
                             interpolatedArray);
  } else if (dimIndex0 == 1 && dimIndex1 == 0) {
    result_ =
      trilinearInterpolation(varName1, varValues1, obValDim1N, ranges[dimIndex1], coordTrans1,
                             varName0, varValues0, obValDim0N, ranges[dimIndex0], coordTrans0,
                             varName2, varValues2, obValDim2N, ranges[dimIndex2], coordTrans2,
                             interpolatedArray);
  } else if (dimIndex0 == 2 && dimIndex1 == 0) {
    result_ =
      trilinearInterpolation(varName1, varValues1, obValDim1N, ranges[dimIndex1], coordTrans1,
                             varName2, varValues2, obValDim2N, ranges[dimIndex2], coordTrans2,
                             varName0, varValues0, obValDim0N, ranges[dimIndex0], coordTrans0,
                             interpolatedArray);
  } else if (dimIndex0 == 1 && dimIndex1 == 2) {
    result_ =
      trilinearInterpolation(varName2, varValues2, obValDim2N, ranges[dimIndex2], coordTrans2,
                             varName0, varValues0, obValDim0N, ranges[dimIndex0], coordTrans0,
                             varName1, varValues1, obValDim1N, ranges[dimIndex1], coordTrans1,
                             interpolatedArray);
  } else if (dimIndex0 == 2 && dimIndex1 == 1) {
    result_ =
      trilinearInterpolation(varName2, varValues2, obValDim2N, ranges[dimIndex2], coordTrans2,
                             varName1, varValues1, obValDim1N, ranges[dimIndex1], coordTrans1,
                             varName0, varValues0, obValDim0N, ranges[dimIndex0], coordTrans0,
                             interpolatedArray);
  } else {
    // This should never be executed
    throw eckit::BadParameter("Logic for variable order has failed", Here());
  }

  resultSet_ = true;
}



// Explicit instantiations
template class DataExtractor<float>;
template class DataExtractor<int>;
template class DataExtractor<std::string>;

}  // namespace ufo
