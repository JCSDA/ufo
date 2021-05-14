/*
 * (C) Copyright 2021 Met Office UK
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_UTILS_DATAEXTRACTOR_DATAEXTRACTOR_H_
#define UFO_UTILS_DATAEXTRACTOR_DATAEXTRACTOR_H_

#include <algorithm>           // sort
#include <functional>          // greater
#include <limits>              // std::numeric_limits
#include <list>                // list
#include <memory>              // unique_ptr
#include <sstream>             // stringstream
#include <string>
#include <unordered_map>
#include <utility>             // pair
#include <vector>

#include "boost/variant.hpp"
#include "Eigen/Dense"         // Eigen Arrays and Matrices

#include "eckit/exception/Exceptions.h"
#include "oops/util/Logger.h"
#include "oops/util/missingValues.h"
#include "ufo/filters/Variables.h"
#include "ufo/utils/RecursiveSplitter.h"

namespace ioda {
class Variable;
}


namespace ufo {

class DataExtractorBackend;

/// \brief Method used by the DataExtractor to map the value of an ObsSpace variable to
/// a range of slices of the interpolated array along the dimension indexed by that variable.
enum class InterpMethod {
  /// \brief Select slices where the indexing coordinate matches exactly the value of the
  /// corresponding ObsSpace variable.
  ///
  /// If no match is found, an exception is thrown unless there are slices where the indexing
  /// coordinate is set to the missing value placeholder; in this case these slices are selected
  /// instead. This can be used to define a fallback value (used if there is no exact match).
  ///
  /// This is the only method that can be used for variables of type 'string'.
  EXACT,

  /// \brief Select slices where the indexing coordinate is closest to the value of the
  /// corresponding ObsSpace variable.
  ///
  /// In case of a tie (e.g. if the value of the ObsSpace variable is 3 and the coordinate contains
  /// values 2 and 4, but not 3), the smaller of the candidate coordinate values is used (in this
  /// example, 2).
  NEAREST,

  /// \brief Select slices corresponding to the least value of the indexing coordinate
  /// greater than or equal to the value of the corresponding ObsSpace variable.
  LEAST_UPPER_BOUND,

  /// \brief Select slices corresponding to the greatest value of the indexing coordinate
  /// less than or equal to the value of the corresponding ObsSpace variable.
  GREATEST_LOWER_BOUND,

  /// \brief Perform a piecewise linear interpolation along the dimension indexed by the ObsSpace
  /// variable.
  ///
  /// This method can only be used for the last indexing variable.
  LINEAR
};


/// \brief This class makes it possible to extract and interpolate data loaded from a file.
///
/// Currently the following file formats are supported: NetCDF and CSV. See
/// DataExtractorNetCDFBackend and DataExtractorCSVBackend for more information about these
/// formats.
///
/// In both cases, the files will contain a _payload array_ (the array from which data will be
/// extracted) and one or more _coordinate arrays_ indexing the payload array. This class makes it
/// possible to rapidly extract a value from the payload array corresponding to particular values
/// of the coordinates, or to interpolate multiple values from this array. Coordinate matching can
/// be exact or approximate (looking for the nearest match). It is also possible to perform a
/// piecewise linear interpolation of the data along one coordinate axis. (Multidimensional linear
/// interpolation is not supported.)
///
/// Here's how to use this class:
///
/// * Call the constructor to load data from an input file.
/// * Call scheduleSort() one or more times, each time passing to it the name of a coordinate
///   array and the matching method to be used for this coordinate. This will determine the order
///   in which coordinates will be matched.
/// * Call sort(). This will prepare internal data structures required for a rapid search for
///   matching coordinates.
/// * To extract a value from the payload array for a particular data point, pass the values of
///   successive coordinates of that point to calls to extract() (in the order matching the order of
///   the preceding calls to scheduleSort()). Then call getResult() to retrieve the extracted value.
///
/// Here is a summary of particulars to the extraction/interpolation algorithms available:
/// - Nearest neighbour 'interpolation' chooses the **first** nearest value to be found in the case
///   of equidistant neighbours of different values. The algorithm will then return the one or more
///   locations corresponding to this nearest value.  Let's illustrate by extracting the nearest
///   neighbours to 1.5:
///
///   [1, 1, 2, 3, 4, 5]
///
///   Both 1 and 2 are equidistant, but we take the first found equidistant neighbour (1), then
///   return all indices matching this neighbour (indices 0 and 1 in this case).
class DataExtractor
{
 public:
  /// \brief Create an object that can be used to extract data loaded from a file.
  /// \details This object is capable of sorting the data from this file extracting the relevant
  /// values for a given observation as well as performing linear interpolation to derive the final
  /// value.
  /// \param[in] filepath Path to the input file.
  /// \param[in] group Group containing the payload variable.
  explicit DataExtractor(const std::string &filepath, const std::string &group);

  /// \brief Update the instruction on how to sort the data for the provided variable name.
  /// \details This works iteratively by further splitting the RecursiveSplitter sub-groups
  /// according to the variable name provided.  By this, it is possible to sort the data in
  /// such a way to ensure that extraction always results in 1 contiguous chunk.
  /// No data or coordinates are actually physically sorted yet.  Special treatment for float
  /// type variables, where this is used to sort each of the sub-groups.
  /// \param[in] varName is the name of the coordinate axis to sort.
  /// \param[in] method is the interpolation/extraction method to use for this coordinate.
  /// \internal This member function call corresponds to a RecursiveSplitter.groupBy call, useful
  /// to sort according to nearesr/exact match variables.  In the special case of float type,
  /// RecursiveSplitter.sortGroupsBy is used.
  void scheduleSort(const std::string &varName, const InterpMethod &method);

  /// \brief Finalise the sort, sorting each of the coordinates indexing the axes of the array to
  /// be interpolated, as well as that array itself.
  /// \details Utilising the instructions provided by the user calling the scheduleSort() member
  /// function, we now physically sort the array itself along with all coordinates which
  /// describe it.
  /// \internal Applies the RecursiveSplitter object and necessarily creates copies to achieve
  /// this sort.
  void sort();

  /// \brief Perform extract, given an observation value for the coordinate associated with this
  /// extract iteration.
  /// \details Calls the relevant extract methood (linear, nearest or exact), corresponding to the
  /// coordinate associated with this extract iteration (along with the associated interpolation
  /// method).  The extract then utilises the observation value ('obVal') to perform this
  /// extraction.
  /// \param[in] obVal is the observation value used for the extract operation.
  template<typename T>
  void extract(const T &obVal) {
    if (nextCoordToExtractBy_ == coordsToExtractBy_.cend()) {
      throw eckit::UserError("Too many extract() calls made for the expected number of variables.",
                             Here());
    }

    // Fetch coordinate values
    const std::vector<T> &coordValues = boost::get<std::vector<T>>(nextCoordToExtractBy_->values);
    // Perform the extraction using the selected method
    switch (nextCoordToExtractBy_->method) {
      case InterpMethod::EXACT:
        exactMatch(nextCoordToExtractBy_->name, coordValues,
                   nextCoordToExtractBy_->payloadDim, obVal);
        break;
      case InterpMethod::LINEAR:
        result_ = getResult(nextCoordToExtractBy_->name, coordValues,
                            nextCoordToExtractBy_->payloadDim, obVal);
        resultSet_ = true;
        break;
      case InterpMethod::NEAREST:
        nearestMatch(nextCoordToExtractBy_->name, coordValues,
                     nextCoordToExtractBy_->payloadDim, obVal);
        break;
      case InterpMethod::LEAST_UPPER_BOUND:
        leastUpperBoundMatch(nextCoordToExtractBy_->name, coordValues,
                             nextCoordToExtractBy_->payloadDim, obVal);
        break;
      case InterpMethod::GREATEST_LOWER_BOUND:
        greatestLowerBoundMatch(nextCoordToExtractBy_->name, coordValues,
                                nextCoordToExtractBy_->payloadDim, obVal);
        break;
      default:
        throw eckit::UserError("Only 'linear', 'nearest', 'exact', 'least upper bound' and "
                               "'greatest lower bound' interpolation methods supported.",
                               Here());
    }
    ++nextCoordToExtractBy_;
  }

  /// \brief Fetch the final interpolated value.
  /// \details This will only be succesful if previous calls to extract() have produced a single
  /// value to return.
  float getResult();

 private:
  /// \brief Fetch the final interpolated value at the 'location' `obVal`.
  ///
  /// \details It is assumed that previous calls to extract() have extracted a single 1D slice of
  /// the interpolated array parallel to the axis `varName`. This function returns the value
  /// produced by piecewise linear interpolation of this slice at the point `obVal`.
  ///
  /// \param[in] varName is the name of the coordinate along which to interpolate.
  /// \param[in] varValues is the vector of values of that coordinate.
  /// \param[in] dimIndex is the dimension of the payload array indexed by the coordinate.
  /// \param[in] obVal is the interpolation location.
  template<typename T>
  float getResult(const std::string &varName, const std::vector<T> &varValues,
                  int dimIndex, const T &obVal) {
    // Sanity check constraint
    int sizeDim0 = constrainedRanges_[0].end - constrainedRanges_[0].begin;
    int sizeDim1 = constrainedRanges_[1].end - constrainedRanges_[1].begin;
    if ((dimIndex == 1 && !(sizeDim1 > 1 && sizeDim0 == 1)) ||
        (dimIndex == 0 && !(sizeDim0 > 1 && sizeDim1 == 1))) {
      throw eckit::Exception("Linear interpolation failed - data must be 1D.", Here());
    }

    // Constrain our index range in the relevant dimension.
    const Range &range = constrainedRanges_[static_cast<size_t>(dimIndex)];

    if ((obVal > varValues[range.end - 1]) || (obVal < varValues[range.begin])) {
      throw eckit::Exception("Linear interpolation failed, value is beyond grid extent."
                             "No extrapolation supported.",
                             Here());
    }
    // Find first index of varValues >= obVal
    int nnIndex = std::lower_bound(varValues.begin() + range.begin,
                                     varValues.begin() + range.end,
                                     obVal) - varValues.begin();

    // Determine upper or lower indices from this
    if (varValues[nnIndex] == obVal) {
      // No interpolation required (is equal)
      float res;
      if (dimIndex == 1) {
        res = static_cast<float>(interpolatedArray2D_(constrainedRanges_[0].begin, nnIndex));
      } else {
        res = static_cast<float>(interpolatedArray2D_(nnIndex, constrainedRanges_[1].begin));
      }
      return res;
    }
    // Linearly interpolate between these two indices.
    auto zUpper = *(interpolatedArray2D_.data());
    auto zLower = *(interpolatedArray2D_.data());
    if (dimIndex == 1) {
      zUpper = interpolatedArray2D_(constrainedRanges_[0].begin, nnIndex);
      zLower = interpolatedArray2D_(constrainedRanges_[0].begin, nnIndex-1);
    } else {
      zUpper = interpolatedArray2D_(nnIndex, constrainedRanges_[1].begin);
      zLower = interpolatedArray2D_(nnIndex-1, constrainedRanges_[1].begin);
    }
    float res = ((static_cast<float>(obVal - varValues[nnIndex-1]) /
                  static_cast<float>(varValues[nnIndex] - varValues[nnIndex-1])) *
                 (zUpper - zLower)) + zLower;
    return res;
  }

  float getResult(const std::string &varName, const std::vector<std::string> &varValues,
                  int dimIndex, const std::string &obVal) {
    throw eckit::UserError("VarName: " + varName +
                           " - linear interpolation not compatible with string type.", Here());
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
  /// \param[in] varName is the name of the coordinate to match against.
  /// \param[in] varValues is the vector of values of that coordinate.
  /// \param[in] dimIndex is the dimension of the payload array indexed by the coordinate.
  /// \param[in] obVal is the value to match.
  template<typename T>
  void nearestMatch(const std::string &varName, const std::vector<T> &varValues,
                    int dimIndex, const T &obVal) {
    // Constrain our index range in the relevant dimension.
    Range &range = constrainedRanges_[static_cast<size_t>(dimIndex)];

    // Find first index of varValues >= obVal
    int nnIndex = std::lower_bound(varValues.begin() + range.begin,
                                   varValues.begin() + range.end,
                                   obVal) - varValues.begin();
    if (nnIndex >= range.end) {
      nnIndex = range.end - 1;
    }

    // Now fetch the nearest neighbour index (lower index prioritised for different values with
    // same distance)
    T dist = std::abs(varValues[nnIndex] - obVal);
    if ((varValues[nnIndex] > obVal) && (nnIndex > range.begin) &&
        (std::abs(varValues[nnIndex - 1] - obVal) <= dist))
      nnIndex--;

    // Now find **same value** equidistant neighbours
    auto bounds = std::equal_range(varValues.begin() + range.begin,
                                   varValues.begin() + range.end,
                                   varValues[nnIndex]);
    range = {static_cast<int>(bounds.first - varValues.begin()),
             static_cast<int>(bounds.second - varValues.begin())};
    oops::Log::debug() << "Nearest match; name: " << varName << " range: " <<
      range.begin << "," << range.end << std::endl;
  }

  void nearestMatch(const std::string &varName, const std::vector<std::string> &varValues,
                    int dimIndex, const std::string &obVal) {
    throw eckit::UserError("Nearest match not compatible with string type.", Here());
  }

  /// \brief Update our extract constraint based on an exact match against the specified coordinate
  /// indexing a dimension of the payload array.
  ///
  /// \param[in] varName is the name of the coordinate to match against.
  /// \param[in] varValues is the vector of values of that coordinate.
  /// \param[in] dimIndex is the dimension of the payload array indexed by the coordinate.
  /// \param[in] obVal is the value to match.
  template<typename T>
  void exactMatch(const std::string &varName, const std::vector<T> &varValues,
                  int dimIndex, const T &obVal) {
    // Constrain our index range in the relevant dimension.
    Range &range = constrainedRanges_[static_cast<size_t>(dimIndex)];

    // Find the first and last matching index
    auto bounds = std::equal_range(varValues.begin() + range.begin,
                                   varValues.begin() + range.end,
                                   obVal);
    if (bounds.first == bounds.second) {
      // No matching coordinate found. If the coordinate contains a 'missing value' entry,
      // use it as a fallback. (If it doesn't, the 'bounds' range will stay empty, so an error will
      // be reported).
      bounds = std::equal_range(varValues.begin() + range.begin,
                                varValues.begin() + range.end,
                                util::missingValue(obVal));
    }

    range = {static_cast<int>(bounds.first - varValues.begin()),
             static_cast<int>(bounds.second - varValues.begin())};

    if (range.begin == range.end) {
      std::stringstream msg;
      msg << "No match found for exact match extraction of value '" << obVal
          << "' of the variable '" << varName << "'";
      throw eckit::Exception(msg.str(), Here());
    }
    oops::Log::debug() << "Exact match; name: " << varName << " range: " <<
      range.begin << "," << range.end << std::endl;
  }

  /// \brief Update our extract constraint based on a least-upper-bound match against the specified
  /// coordinate indexing a dimension of the payload array.
  ///
  /// \param[in] varName is the name of the coordinate to match against.
  /// \param[in] varValues is the vector of values of that coordinate.
  /// \param[in] dimIndex is the dimension of the payload array indexed by the coordinate.
  /// \param[in] obVal is the value to match.
  template<typename T>
  void leastUpperBoundMatch(const std::string &varName, const std::vector<T> &varValues,
                            int dimIndex, const T &obVal) {
    // Constrain our index range in the relevant dimension.
    Range &range = constrainedRanges_[static_cast<size_t>(dimIndex)];

    // Find index of the first varValues >= obVal
    typedef typename std::vector<T>::const_iterator It;
    const It rangeBegin(varValues.begin() + range.begin);
    const It rangeEnd(varValues.begin() + range.end);

    const It leastUpperBoundIt = std::lower_bound(rangeBegin, rangeEnd, obVal);
    if (leastUpperBoundIt == rangeEnd) {
      std::stringstream msg;
      msg << "No match found for 'least upper bound' extraction of value '" << obVal
          << "' of the variable '" << varName << "'";
      throw eckit::Exception(msg.str(), Here());
    }

    // Find the range of items with the same value of this coordinate
    const auto bounds = std::equal_range(rangeBegin, rangeEnd, *leastUpperBoundIt);
    range = {static_cast<int>(bounds.first - varValues.begin()),
             static_cast<int>(bounds.second - varValues.begin())};
    oops::Log::debug() << "Least upper bound match; name: " << varName << " range: "
                       << range.begin << "," << range.end << std::endl;
  }

  void leastUpperBoundMatch(const std::string &varName, const std::vector<std::string> &varValues,
                            int dimIndex, const std::string &obVal) {
    throw eckit::UserError("The 'least upper bound' method cannot be used for string variables.",
                           Here());
  }

  /// \brief Update our extract constraint based on a greatest-lower-bound match against the
  /// specified coordinate indexing a dimension of the payload array.
  ///
  /// \param[in] varName is the name of the coordinate to match against.
  /// \param[in] varValues is the vector of values of that coordinate.
  /// \param[in] dimIndex is the dimension of the payload array indexed by the coordinate.
  /// \param[in] obVal is the value to match, against the NetCDF coordinate of name 'varName'.
  template<typename T>
  void greatestLowerBoundMatch(const std::string &varName, const std::vector<T> &varValues,
                               int dimIndex, const T &obVal) {
    // Constrain our index range in the relevant dimension.
    Range &range = constrainedRanges_[static_cast<size_t>(dimIndex)];

    // Find index of the last varValues <= obVal

    typedef typename std::vector<T>::const_reverse_iterator ReverseIt;
    typedef std::greater<T> Compare;
    const ReverseIt reverseRangeBegin(varValues.begin() + range.end);
    const ReverseIt reverseRangeEnd(varValues.begin() + range.begin);

    const ReverseIt greatestLowerBoundIt =
        std::lower_bound(reverseRangeBegin, reverseRangeEnd, obVal, Compare());
    if (greatestLowerBoundIt == reverseRangeEnd) {
      std::stringstream msg;
      msg << "No match found for 'greatest lower bound' extraction of value '" << obVal
          << "' of the variable '" << varName << "'";
      throw eckit::Exception(msg.str(), Here());
    }

    // Find the range of items with the same value of this coordinate
    const auto bounds = std::equal_range(varValues.begin() + range.begin,
                                         varValues.begin() + range.end,
                                         *greatestLowerBoundIt);
    range = {static_cast<int>(bounds.first - varValues.begin()),
             static_cast<int>(bounds.second - varValues.begin())};
    oops::Log::debug() << "Greatest lower bound match; name: " << varName << " range: "
                       << range.begin << "," << range.end << std::endl;
  }

  void greatestLowerBoundMatch(const std::string &varName,
                               const std::vector<std::string> &varValues,
                               int dimIndex, const std::string &obVal) {
    throw eckit::UserError("The 'greatest lower bound' method cannot be used for string variables.",
                           Here());
  }

  /// \brief Reset the extraction range for this object.
  /// \details Each time an exactMatch, nearestMatch, leastUpperBoundMatch or
  /// greatestLowerBoundMatch call is made for one or more variable,
  /// the extraction range is further constrained to match our updated match conditions.  After
  /// the final 'extract' is made (i.e. an interpolated value is derived) it is desirable to reset
  /// the extraction range by calling this method.
  /// \internal This is called by the getObsErrorValue member functions just before returning the
  /// interpolated value.
  void resetExtract();

  /// \brief Fetch the coordinate data of specified name.
  /// \details We cache this data so as not to require performing multiple loads from disk.
  template <typename T>
  T& get(const std::string &key) {
    try {
      return boost::get<T> (coordsVals_[key]);
    } catch (boost::bad_get) {
      throw eckit::BadParameter("Unable to find coordinate with this type", Here());
    }
  }

  /// \brief Load all data from the input file.
  void load(const std::string &filepath, const std::string &interpolatedArrayGroup);

  /// \brief Create a backend able to read file \p filepath.
  static std::unique_ptr<DataExtractorBackend> createBackendFor(const std::string &filepath);

  // Object represent the extraction range in both dimensions.
  struct Range {int begin, end;};
  std::array<Range, 2> constrainedRanges_;

  // Container holding coordinate arrays (of all supported types) loaded from the input file.
  typedef boost::variant<std::vector<int>,
                         std::vector<float>,
                         std::vector<std::string>> CoordinateValues;
  std::unordered_map<std::string, CoordinateValues> coordsVals_;
  // The array to be interpolated (the payload array).
  Eigen::ArrayXXf interpolatedArray2D_;
  float result_;
  bool resultSet_;
  // Container for re-ordering our data
  std::vector<ufo::RecursiveSplitter> splitter_;

  /// Maps coordinate names to dimensions (0 or 1) of the payload array
  std::unordered_map<std::string, int> coord2DimMapping_;
  /// Maps dimensions of the payload array (0 or 1) to coordinate names
  std::vector<std::vector<std::string>> dim2CoordMapping_;

  /// Coordinate used for data extraction from the payload array.
  struct Coordinate {
    /// Coordinate name
    std::string name;
    /// Coordinate values
    const CoordinateValues &values;
    /// Extraction method to use
    InterpMethod method;
    /// Axis of the payload array indexed by the coordinate (0 or 1)
    int payloadDim;
  };

  /// Coordinates to use in successive calls to extract().
  std::vector<Coordinate> coordsToExtractBy_;
  std::vector<Coordinate>::const_iterator nextCoordToExtractBy_;
};

}  // namespace ufo

#endif  // UFO_UTILS_DATAEXTRACTOR_DATAEXTRACTOR_H_
