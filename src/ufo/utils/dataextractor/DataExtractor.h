/*
 * (C) Copyright 2021 Met Office UK
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_UTILS_DATAEXTRACTOR_DATAEXTRACTOR_H_
#define UFO_UTILS_DATAEXTRACTOR_DATAEXTRACTOR_H_

#include <array>
#include <limits>              // std::numeric_limits
#include <memory>              // unique_ptr
#include <string>
#include <unordered_map>
#include <vector>

#include <boost/multi_array.hpp>
#include <boost/optional.hpp>
#include <boost/variant.hpp>
#include <boost/variant/multivisitors.hpp>

#include "eckit/exception/Exceptions.h"

#include "oops/util/missingValues.h"

#include "ufo/utils/dataextractor/ConstrainedRange.h"
#include "ufo/utils/RecursiveSplitter.h"


namespace ioda {
class Variable;
struct Named_Variable;
}


namespace ufo {

template <typename T>
using DataExtractorPayload = boost::multi_array<T, 3>;

/// \brief Coordinate transformation used by the DataExtractor prior to interpolation.
enum class CoordinateTransformation {
  /// \brief No transformation of coordinates.
  NONE,

  /// \brief Log-linear transformation of coordinates, ignoring missing values.
  LOGLINEAR
};


/// \brief Determine whether the value provided is out-of-bounds for the given vector and
/// associated range (constraint).
///
/// \param[in] obVal
///   Interpolation location.
/// \param[in] varValues
///   Vector of values (coordinate).
/// \param[in] range
///   Defines how to constrain (slice) `varValues`.
template <typename T>
bool isOutOfBounds(const T obVal, const std::vector<T> &varValues, const ConstrainedRange &range) {
  if ((obVal > varValues[range.end()-1]) || (obVal < varValues[range.begin()]))
    return true;
  return false;
}


/// \brief Fetch a 1D sliced view of a boost multi_array object.
///
/// \details Given a set of constraints (ranges), return a 1D sliced view of the array.
/// This will be along the dimension specified.  The constraints must constrain all but the
/// dimension that the 1D varies.  Note that the 1D slice returned is not itself constrained or
/// reordered.

/// \param[in] array
///   Some boost array for us to return its sliced view.
/// \param[in] dimIndex
///   Dimension index corresponding to this 1D slice.  Though this argument might be considered
///   superfluous with also passing 'ranges', it does serve to ensure that the 1D slice returned
///   corresponds to the dimension index you were expecting.
/// \param[in] ranges
///   Constraints for the array provided, used for extracting a 1D slice.
///
/// Example:
/// \code
///   // Define some 3D boost array
///   boost::multi_array<float, 3> array (boost::extents[2][3][2]);
///   for (int i=0; i<2; i++)
///     for (int j=0; j<3; j++)
///       for (int k=0; k<2; k++)
///       array[i][j][k] = ...
///
///   // Let's fetch a 1D slice, constraining our first and last dimensions.
///   std::array<ConstrainedRange, 3> ranges {ConstrainedRange(2), ConstrainedRange(3),
///                                           ConstrainedRange(2)};
///   ranges[0].constrain(0, 1);  // Let's pick the very first index...
///   ranges[2].constrain(0, 1);  // Let's pick the very first index...
///
///   auto arraySliceView = get1DSlice(array, 1, ranges);
/// \endcode
template <typename T>
const typename DataExtractorPayload<T>::template const_array_view<1>::type get1DSlice(
    const DataExtractorPayload<T> &array, const size_t dimIndex,
    const std::array<ConstrainedRange, 3> &ranges) {

  // Sanity check constraints
  for (size_t dim=0; dim < ranges.size(); dim++)
    if ((dim != dimIndex) && (ranges[dim].size() > 1))
      throw eckit::Exception(
        "Unable to fetch a 1D array slice with the provided constraints.", Here());

  typedef boost::multi_array_types::index_range range_t;
  typedef typename DataExtractorPayload<T>::template const_array_view<1>::type view1D;
  typename DataExtractorPayload<T>::index_gen indices;

  if (dimIndex == 0) {
    return array[indices[range_t()][ranges[1].begin()][ranges[2].begin()]];
  } else if (dimIndex == 1) {
    return array[indices[ranges[0].begin()][range_t()][ranges[2].begin()]];
  } else if (dimIndex == 2) {
    return array[indices[ranges[0].begin()][ranges[1].begin()][range_t()]];
  } else {
    // We shouldn't ever end up here (exception should be thrown earlier).
    throw eckit::UserError("Invalid dimension index, expecting value mappings corresponding to 1 "
                           "of 3 axes.", Here());
  }
}


/// \brief Fetch a 2D sliced view of a boost multi_array object.
///
/// \details Given a set of constraints (ranges), return a 2D sliced view of the array.
/// This will be along the two dimensions specified.  The constraints must constrain all
/// but the two dimensions specified.  Note that the 2D slice returned is not itself constrained
/// or reordered.
///
/// \param[in] array
///   Some boost array for us to return its sliced view.
/// \param[in] dimIndex0
///   First of two 'array' dimension indices corresponding to our 2D slice.
///   Though this argument might be considered superfluous with also passing 'ranges', it does
///   serve to ensure that the 2D slice returned corresponds to the dimension index you were
///   expecting.
/// \param[in] dimIndex1
///   Second of two 'array' dimension indices corresponding to our 2D slice.
///   Though this argument might be considered superfluous with also passing 'ranges', it does
///   serve to ensure that the 2D slice returned corresponds to the dimension index you were
///   expecting.
/// \param[in] ranges
///   Constraints for the array provided, used for extracting a 2D slice.
///
/// Example:
/// \code
///   // Define some 3D boost array
///   boost::multi_array<float, 3> array (boost::extents[2][3][2]);
///   for (int i=0; i<2; i++)
///     for (int j=0; j<3; j++)
///       for (int k=0; k<2; j++)
///       array[i][j][k] = ...
///
///   // Let's fetch a 2D slice, constraining our first dimension to a specific index.
///   std::array<ConstrainedRange, 3> ranges {ConstrainedRange(2), ConstrainedRange(3),
///                                           ConstrainedRange(2)};
///   ranges[0].constrain(0, 1);  // Let's pick the very first index...
///
///   auto arraySliceView = get2DSlice(array, 1, 2, ranges);
/// \endcode
template <typename T>
const typename DataExtractorPayload<T>::template const_array_view<2>::type get2DSlice(
    const DataExtractorPayload<T> &array, const size_t dimIndex0, const size_t dimIndex1,
    const std::array<ConstrainedRange, 3> &ranges) {

  // Sanity check constraints
  for (size_t dim=0; dim < ranges.size(); dim++)
    if ((dim != dimIndex0) && (dim != dimIndex1) && (ranges[dim].size() > 1))
      throw eckit::Exception(
        "Unable to fetch a 2D array slice with the provided constraints.", Here());

  typedef boost::multi_array_types::index_range range_t;
  typedef typename DataExtractorPayload<T>::template const_array_view<2>::type view2D;
  typename DataExtractorPayload<T>::index_gen indices;

  size_t sumIndex = dimIndex0 + dimIndex1;
  if (sumIndex == 1) {
    return array[indices[range_t()][range_t()][ranges[2].begin()]];
  } else if (sumIndex == 2) {
    return array[indices[range_t()][ranges[1].begin()][range_t()]];
  } else if (sumIndex == 3) {
    return array[indices[ranges[0].begin()][range_t()][range_t()]];
  } else {
    // We shouldn't ever end up here (exception should be thrown earlier).
    throw eckit::UserError("Invalid dimension index, expecting value mappings corresponding to two "
                           "of 3 axes.", Here());
  }
}


/// \brief Perform bilinear interpolation at 'location' at location `obVal0` and `obVal1`,
/// corresponding to the first and second axes.
///
/// \details This function returns the value produced by a bilinear interpolation at point
/// `obVal0, obVal1`.  Where any of the neighbouring points used in the calculation are missing,
/// return the nearest of the non-missing neighbouring points.  Where all 4 neighbours are missing
/// or the coordinate value is out of bounds, return a missing value.
/// This must be a final call in the sequence calls to extract.  This template handles all types
/// except string (see below overloads).
///
/// \param[in] varName0
///   Name of the coordinate describing the first dimension of 'interpolatedArray' over which we
///   are to interpolate.
/// \param[in] varValues0
///   Vector of values of the `varName0` coordinate.
/// \param[in] obVal0
///   Interpolation location along the axis corresponding to `varName0`.
/// \param[in] range0
///   Defines how to constrain (slice) `varValues0` and the corresponding first dimension of
///   'interpolatedArray'.
/// \param[in] varName1
///   Name of the coordinate describing the second dimension of 'interpolatedArray' over which we
///   are to interpolate.
/// \param[in] varValues1
///   Vector of values of the `varName1` coordinate.
/// \param[in] obVal1
///   Interpolation location along the axis corresponding to `varName1`.
/// \param[in] range1
///   Defines how to constrain (slice) `varValues1` and the corresponding second dimension of
/// 'interpolatedArray'.
/// \param[in] interpolatedArray
///   2D interpolated array view.
template <typename T, typename R>
float bilinearInterpolation(
    const std::string &varName0,
    const std::vector<T> &varValues0,
    const T &obVal0,
    const ConstrainedRange &range0,
    const std::string &varName1,
    const std::vector<R> &varValues1,
    const R &obVal1,
    const ConstrainedRange &range1,
    const DataExtractorPayload<float>::const_array_view<2>::type &interpolatedArray) {

  const float missing = util::missingValue<float>();

  if (obVal0 == util::missingValue<T>() || obVal1 == util::missingValue<R>()) {
    return missing;
  }

  if (isOutOfBounds(obVal0, varValues0, range0)) {
      std::stringstream msg;
      msg << "No match found for 'bilinear' interpolation of value '" << obVal0
          << "' of the variable '" << varName0 << "'.  Value is out of bounds.  Consider using "
          << "extrapolation.";
      throw eckit::Exception(msg.str(), Here());
  }
  if (isOutOfBounds(obVal1, varValues1, range1)) {
      std::stringstream msg;
      msg << "No match found for 'bilinear' interpolation of value '" << obVal1
          << "' of the variable '" << varName1 << "'.  Value is out of bounds.  Consider using "
          << "extrapolation.";
      throw eckit::Exception(msg.str(), Here());
  }

  const int nnIndex0 = std::lower_bound(varValues0.begin() + range0.begin(),
                                        varValues0.begin() + range0.end(), obVal0) -
    varValues0.begin();
  const int nnIndex1 = std::lower_bound(varValues1.begin() + range1.begin(),
                                        varValues1.begin() + range1.end(), obVal1) -
    varValues1.begin();

  if ((varValues0[nnIndex0] == obVal0) && (varValues1[nnIndex1] == obVal1)) {
    // No interpolation required.
    return interpolatedArray[nnIndex0][nnIndex1];
  }

  // Setup points
  // ------------
  // Indices
  const int ix1 = nnIndex0-1 >= 0 ? nnIndex0-1 : nnIndex0;
  const int ix2 = ix1 + 1;
  const int iy1 = nnIndex1-1 >= 0 ? nnIndex1-1 : nnIndex1;
  const int iy2 = iy1 + 1;

  // Coord locations
  const T x1 = varValues0[ix1], x2 = varValues0[ix2];
  const R y1 = varValues1[iy1], y2 = varValues1[iy2];

  // Z values at these locations
  const float q11 = interpolatedArray[ix1][iy1], q12 = interpolatedArray[ix1][iy2],
      q22 = interpolatedArray[ix2][iy2], q21 = interpolatedArray[ix2][iy1];

  // Missing data handling
  // - Pick nearest non-missing neighbour of the moore neighbourhood (no diagonals) if any of these
  //   4 neighbours are missing.
  // - Return missing if all 4 of the moore neighbourhood (no diagonals) are missing.
  const int nmissing = (q11 == missing) + (q12 == missing) + (q22 == missing) + (q21 == missing);
  if (nmissing > 0) {
    if (nmissing == 4)
      return missing;
    const std::array<T, 4> xd = {x1, x1, x2, x2};
    const std::array<R, 4> yd = {y1, y2, y2, y1};
    const std::array<float, 4> qd = {q11, q12, q22, q21};

    float minq, mindist = std::numeric_limits<float>::max();
    for (size_t ind=0; ind < qd.size(); ind++) {
      // Euclidean distance between two points.
      const float currdist = std::hypot(obVal0-xd[ind], obVal1-yd[ind]);
      if (currdist < mindist && qd[ind] != missing) {
        mindist = currdist;
        minq = qd[ind];
      }
    }
    return minq;
  }

  // Interpolate at obVal0 and obVal1
  const float denom = static_cast<float>((x2 - x1)*(y2 - y1));
  float res = (((x2 - obVal0)*(y2 - obVal1))) * q11;
  res +=      (((obVal0 - x1)*(y2 - obVal1))) * q21;
  res +=      (((x2 - obVal0)*(obVal1 - y1))) * q12;
  res +=      (((obVal0 - x1)*(obVal1 - y1))) * q22;
  return res/denom;
}


/// \brief Bilinear interpolation template, templated coord1, string coord2.
template <typename T>
float bilinearInterpolation(
    const std::string &varName0,
    const std::vector<T> &varValues0,
    const T &obVal0,
    const ConstrainedRange &range0,
    const std::string &varName1,
    const std::vector<std::string> &varValues1,
    const std::string &obVal1,
    const ConstrainedRange &range1,
    const DataExtractorPayload<float>::const_array_view<2>::type &interpolatedArray) {
  throw eckit::UserError("Bilinear interpolation cannot be performed along coordinate axes indexed "
                         "by string variables such as " + varName1 + ".", Here());
}


/// \brief Bilinear interpolation template, string coord1, templated coord2.
template <typename T>
float bilinearInterpolation(
    const std::string &varName0,
    const std::vector<std::string> &varValues0,
    const std::string &obVal0,
    const ConstrainedRange &range0,
    const std::string &varName1,
    const std::vector<T> &varValues1,
    const T &obVal1,
    const ConstrainedRange &range1,
    const DataExtractorPayload<float>::const_array_view<2>::type &interpolatedArray) {
  throw eckit::UserError("Bilinear interpolation cannot be performed along coordinate axes indexed "
                         "by string variables such as " + varName0 + ".", Here());
}


/// Bilinear interpolation template, string coord1, string coord2.
inline float bilinearInterpolation(
    const std::string &varName0,
    const std::vector<std::string> &varValues0,
    const std::string &obVal0,
    const ConstrainedRange &range0,
    const std::string &varName1,
    const std::vector<std::string> &varValues1,
    const std::string &obVal1,
    const ConstrainedRange &range1,
    const DataExtractorPayload<float>::const_array_view<2>::type &interpolatedArray) {
  throw eckit::UserError("Bilinear interpolation cannot be performed along coordinate axes indexed "
                         "by string variables such as " + varName0 + " or " +
                         varName1 + ".", Here());
}


template <typename ExtractedValue>
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
  /// This method can only be used for the last indexing variable.  It is supported only when
  /// DataExtractor produces values of type `float`, but not `int` or `std::string`.
  LINEAR,

  /// \brief Perform a bilinear interpolation along two dimensions indexed by the ObsSpace
  /// variables.
  ///
  /// This method can only be used for the final two indexing variables (see
  /// `bilinearInterpolation`).  It is supported only when DataExtractor produces values of type
  /// `float`, but not `int` or `std::string`.
  BILINEAR,

  /// \brief Perform a trilinear interpolation along three dimensions indexed by the ObsSpace
  /// variables.
  ///
  /// This method can only be used for the final three indexing variables (see
  /// `trilinearInterpolation`).  It is supported only when DataExtractor produces values of type
  /// `float`, but not `int` or `std::string`.
  TRILINEAR
};


/// \brief Extrapolation method used by the DataExtractor interpolation method for a given
/// variable-interpmethod pair.
///
/// Extrapolation is where the value to be extracted/interpolated is out-of-bounds.
enum class ExtrapolationMode {
  /// \brief Pick nearest index when out-of-bounds.
  NEAREST,

  /// \brief Throw an exception when out-of-bounds.
  ERROR,

  /// \brief Return 'missing', meaning that any subsequent extraction stages are then ignored.
  MISSING
};

/// \brief Equidistance choice used by the DataExtractor to determine whether to use the first
/// or last index in the case where a value is equidistant between two lookup values.
///
enum class EquidistantChoice {
  /// \brief Select the first matching index in the case of a tiebreak.
  FIRST,

  /// \brief Select the last matching index in the case of a tiebreak.
  LAST
};

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
    const DataExtractorPayload<float>::const_array_view<3>::type &interpolatedArray);

/// \brief Apply log-linear transformation to a vector of values.
void applyLogLinearTransform(const std::string &varName,
                             std::vector<float> &varValues);
/// \brief Apply log-linear transformation to a single value.
void applyLogLinearTransform(const std::string &varName,
                             float &varValues);

/// \brief This class makes it possible to extract and interpolate data loaded from a file.
///
/// \tparam ExtractedValue
///   Type of the values to extract. Must be `float`, `int` or `std::string`.
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
/// piecewise linear interpolation of the data along one coordinate axis or bilinear interpolation
/// along two axes.
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
template <typename ExtractedValue>
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
  /// \param[in] extrapMode is the extrapolation mode which applies when the value to be
  /// extracted is out-of-bounds.
  /// \ param[in] equidistantChoice is the choice of index when a value is equidistant between
  /// two indices.
  /// \ param[in] coordinateTransformation is the transformation applied to the coordinates.
  /// \internal This member function call corresponds to a RecursiveSplitter.groupBy call, useful
  /// to sort according to nearesr/exact match variables.  In the special case of float type,
  /// RecursiveSplitter.sortGroupsBy is used.
  void scheduleSort(const std::string &varName, const InterpMethod &method,
                    const ExtrapolationMode &extrapMode,
                    const EquidistantChoice &equidistantChoice,
                    const CoordinateTransformation &coordinateTransformation);

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
  /// \details Calls the relevant extract method (linear, nearest or exact), corresponding to the
  /// coordinate associated with this extract iteration (along with the associated interpolation
  /// method).  The extract then utilises the observation value ('obVal') to perform this
  /// extraction.
  /// \param[in] obVal is the observation value used for the extract operation.
  void extract(float obVal);
  void extract(int obVal);
  void extract(const std::string &obVal);

  /// \brief Perform extract, given the observation values associated with the current extract
  /// iteration and the next.
  /// \details Calls the relevant extract method (bilinear), corresponding to the coordinate
  /// associated with this extract iteration and the next (along with the associated interpolation
  /// method).  This method actually functions as two iterations, passing the current iteration
  /// coordinate and the next iteration coordinate.  These are passed to the underlying binary
  /// operation (for example bilinear interpolation).
  /// \param[in] obValDim0 is the observation value used for the extract operation corresponding
  /// to the first coordinate utilised by the underlying method.
  /// \param[in] obValDim1 is the observation value used for the extract operation corresponding
  /// to the second coordinate utilised by the underlying method.
  template <typename T, typename R>
  void extract(T obValDim0, R obValDim1) {
    if (nextCoordToExtractBy_ == coordsToExtractBy_.cend())
      throw eckit::UserError("Too many extract() calls made for the expected number of variables.",
                             Here());

    // Perform the extraction using the selected method
    if (nextCoordToExtractBy_->method == InterpMethod::BILINEAR)
      maybeExtractByBiLinearInterpolation(obValDim0, obValDim1);
    else
        throw eckit::UserError("Only bilinear method supports two variables as arguments.", Here());
    ++nextCoordToExtractBy_;
  }

  /// \brief Perform extract, given the observation values associated with the current extract
  /// iteration and the next.
  /// \details Calls the relevant extract method (trilinear), corresponding to the coordinates
  /// associated with this extract iteration (along with the associated interpolation
  /// method).  This method actually functions as three iterations, passing the current iteration
  /// coordinate and the next two iteration coordinates.
  /// These are passed to the underlying operation.
  /// \param[in] obValDim0 is the observation value used for the extract operation corresponding
  /// to the first coordinate utilised by the underlying method.
  /// \param[in] obValDim1 is the observation value used for the extract operation corresponding
  /// to the second coordinate utilised by the underlying method.
  /// \param[in] obValDim2 is the observation value used for the extract operation corresponding
  /// to the third coordinate utilised by the underlying method.
  void extract(float obValDim0, float obValDim1, float obValDim2) {
    if (nextCoordToExtractBy_ == coordsToExtractBy_.cend())
      throw eckit::UserError("Too many extract() calls made for the expected number of variables.",
                             Here());

    // Perform the extraction using the selected method
    if (nextCoordToExtractBy_->method == InterpMethod::TRILINEAR)
      maybeExtractByTriLinearInterpolation(obValDim0, obValDim1, obValDim2);
    else
      throw eckit::UserError("Only trilinear method supports three variables as arguments.",
                             Here());
    ++nextCoordToExtractBy_;
  }

  /// \brief Fetch the final interpolated value.
  /// \details This will only be successful if previous calls to extract() have produced a single
  /// value to return.
  ExtractedValue getResult();

 private:
  /// \brief Common implementation of the overloaded public function extract().
  template <typename T>
  void extractImpl(const T &obVal);

  /// \brief Apply extrapolation stage.
  ///
  /// Return a new obVal after applying the relevant extrapolation mode for the given value.
  /// It can also populate `result_` (as indicated by `resultSet_`), the final interpolation
  /// result where applicable (where out-of-bounds should return missing).  Where `result_`
  /// is populated, all subsequent extraction stages are then ignored.
  template <typename T>
  T applyExtrapolation(const T &obVal) {
    if (resultSet_)
      return obVal;

    const std::vector<T> varValues = boost::get<std::vector<T>>(nextCoordToExtractBy_->values);
    ConstrainedRange &range = constrainedRanges_[nextCoordToExtractBy_->payloadDim];
    const std::string &varName = nextCoordToExtractBy_->name;

    // Handle extrapolation mode
    T obValN = obVal;
    if ((nextCoordToExtractBy_->method == ufo::InterpMethod::EXACT) &&
        (nextCoordToExtractBy_->extrapMode != ufo::ExtrapolationMode::ERROR))
      throw eckit::BadParameter("Only 'error' extrapolation mode supported for 'exact' method "
                                "extract.", Here());
    switch (nextCoordToExtractBy_->extrapMode) {
      case ufo::ExtrapolationMode::ERROR:
        // Error (no extrapolation) is the default behaviour of all methods so no action required.
        break;
      case ufo::ExtrapolationMode::NEAREST:
        if (obVal > varValues[range.end()-1]) {
          obValN = varValues[range.end()-1];
        } else if (obVal < varValues[range.begin()]) {
          obValN = varValues[range.begin()];
        }
        break;
      case ufo::ExtrapolationMode::MISSING:
        resultSet_ = true;
        result_ = util::missingValue<ExtractedValue>();
        return obValN;
      default:
        throw eckit::Exception("Unrecognised extrapolation mode for '" + varName + "', please "
                               "choose either 'missing', 'nearest' or 'error'.", Here());
    }
    return obValN;
  }

  /// \brief Perform extraction using piecewise linear interpolation, if it's compatible with the
  /// ExtractedValue type in use; otherwise throw an exception.
  template <typename T>
  void maybeExtractByLinearInterpolation(const T &obVal);

  template <typename T, typename R>
  void maybeExtractByBiLinearInterpolation(const T &obValDim0, const R &obValDim1) {
    // Should never be called -- this error should be detected earlier (scheduleSort).
    throw eckit::BadParameter("Bilinear interpolation can be used when extracting floating-point "
                              "values, but not integers or strings.", Here());
  }

  void maybeExtractByTriLinearInterpolation(const float &obValDim0,
                                            const float &obValDim1,
                                            const float &obValDim2) {
    // Should never be called -- this error should be detected earlier (scheduleSort).
    throw eckit::BadParameter("Trilinear interpolation can be used when extracting floating-point "
                              "values, but not integers or strings.", Here());
  }

  /// \brief Fetch the result produced by previous calls to extract(), none of which may have
  /// used linear interpolation.
  ///
  /// An exception is thrown if these calls haven't produced a unique match of the extraction
  /// criteria.
  ExtractedValue getUniqueMatch() const;

  /// \brief Reset the extraction range for this object.
  /// \details Each time an exactMatch, nearestMatch, leastUpperBoundMatch or
  /// greatestLowerBoundMatch call is made for one or more variable,
  /// the extraction range is further constrained to match our updated match conditions.  After
  /// the final 'extract' is made (i.e. an interpolated value is derived) it is desirable to reset
  /// the extraction range by calling this method.
  /// \internal This is called by the getResult member function just before returning the
  /// interpolated value.
  void resetExtract();

  /// \brief Load all data from the input file.
  void load(const std::string &filepath, const std::string &interpolatedArrayGroup);

  /// \brief Create a backend able to read file \p filepath.
  static std::unique_ptr<DataExtractorBackend<ExtractedValue>> createBackendFor(
      const std::string &filepath);

  // Object represent the extraction range in both dimensions.
  std::array<ConstrainedRange, 3> constrainedRanges_;

  // Container holding coordinate arrays (of all supported types) loaded from the input file.
  typedef boost::variant<std::vector<int>,
                         std::vector<float>,
                         std::vector<std::string>> CoordinateValues;
  std::unordered_map<std::string, CoordinateValues> coordsVals_;
  // The array to be interpolated (the payload array).
  DataExtractorPayload<ExtractedValue> interpolatedArray_;
  // Interpolation result.  Used when utilising linear/bilinear interpolation and also with
  // extrapolation (where applicable).
  ExtractedValue result_;
  // Set to true if result_ is a valid value.
  bool resultSet_;
  // Container for re-ordering our data
  std::vector<ufo::RecursiveSplitter> splitter_;

  /// Maps coordinate names to dimensions (0 or 1) of the payload array
  std::unordered_map<std::string, std::vector<int>> coord2DimMapping_;
  /// Maps dimensions of the payload array (0 or 1) to coordinate names
  std::vector<std::vector<std::string>> dim2CoordMapping_;
  /// Maps coordinate names to their dimensionality.
  std::unordered_map<std::string, size_t> coordNDims_;

  /// Coordinate used for data extraction from the payload array.
  struct Coordinate {
    /// Coordinate name
    std::string name;
    /// Coordinate values
    const CoordinateValues &values;
    /// Extraction method to use
    InterpMethod method;
    /// Extrapolation mode to use
    ExtrapolationMode extrapMode;
    /// Equidistant choice to use
    EquidistantChoice equidistantChoice;
    /// Coordinate transformation to use
    CoordinateTransformation coordinateTransformation;
    /// Axis of the payload array indexed by the coordinate (0 or 1)
    int payloadDim;
  };

  /// Coordinates to use in successive calls to extract().
  std::vector<Coordinate> coordsToExtractBy_;
  typename std::vector<Coordinate>::const_iterator nextCoordToExtractBy_;
};


// Specialization for ExtractedValue = float.
template <>
template <typename T, typename R>
void DataExtractor<float>::maybeExtractByBiLinearInterpolation(
    const T &obValDim0, const R &obValDim1) {
  auto &ranges = constrainedRanges_;

  T obValDim0N = applyExtrapolation(obValDim0);
  if (resultSet_)
    return;
  const size_t dimIndex0 = nextCoordToExtractBy_->payloadDim;
  const std::string &varName0 = nextCoordToExtractBy_->name;
  const std::vector<T> &varValues0 = boost::get<std::vector<T>>(nextCoordToExtractBy_->values);
  ++nextCoordToExtractBy_;  // Consume variable

  if (nextCoordToExtractBy_->method != InterpMethod::BILINEAR)
    throw eckit::BadParameter("Second parameter provided to the Bilinear interpolator is not of "
                              "method 'bilinear'.", Here());
  R obValDim1N = applyExtrapolation(obValDim1);
  if (resultSet_)
    return;
  const size_t dimIndex1 = nextCoordToExtractBy_->payloadDim;
  const std::string &varName1 = nextCoordToExtractBy_->name;
  const std::vector<R> &varValues1 = boost::get<std::vector<R>>(nextCoordToExtractBy_->values);

  auto interpolatedArray = get2DSlice(interpolatedArray_, dimIndex0, dimIndex1,
                                      ranges);
  if (dimIndex1 > dimIndex0) {
    result_ = bilinearInterpolation(varName0, varValues0, obValDim0N, ranges[dimIndex0],
                                    varName1, varValues1, obValDim1N, ranges[dimIndex1],
                                    interpolatedArray);
  } else {
    result_ = bilinearInterpolation(varName1, varValues1, obValDim1N, ranges[dimIndex1],
                                    varName0, varValues0, obValDim0N, ranges[dimIndex0],
                                    interpolatedArray);
  }
  resultSet_ = true;
}



// Specialization of trilinear interpolation for ExtractedValue = float.
template <>
void DataExtractor<float>::maybeExtractByTriLinearInterpolation
  (const float &obValDim0, const float &obValDim1, const float &obValDim2);

}  // namespace ufo

#endif  // UFO_UTILS_DATAEXTRACTOR_DATAEXTRACTOR_H_
