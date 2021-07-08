/*
 * (C) Copyright 2021 Met Office UK
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_UTILS_DATAEXTRACTOR_DATAEXTRACTOR_H_
#define UFO_UTILS_DATAEXTRACTOR_DATAEXTRACTOR_H_

#include <memory>              // unique_ptr
#include <string>
#include <unordered_map>
#include <vector>

#include <boost/multi_array.hpp>
#include <boost/variant.hpp>

#include "ufo/utils/dataextractor/ConstrainedRange.h"
#include "ufo/utils/RecursiveSplitter.h"

namespace ioda {
class Variable;
struct Named_Variable;
}


namespace ufo {

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
  /// This method can only be used for the last indexing variable. It is supported only when the
  /// DataExtractor produces values of type `float`, but not `int` or `std::string`.
  LINEAR
};


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
  void extract(float obVal);
  void extract(int obVal);
  void extract(const std::string &obVal);

  /// \brief Fetch the final interpolated value.
  /// \details This will only be succesful if previous calls to extract() have produced a single
  /// value to return.
  ExtractedValue getResult();

 private:
  /// \brief Common implementation of the overloaded public function extract().
  template <typename T>
  void extractImpl(const T &obVal);

  /// \brief Perform extraction using piecewise linear interpolation, if it's compatible with the
  /// ExtractedValue type in use; otherwise throw an exception.
  template <typename T>
  void maybeExtractByLinearInterpolation(const T &obVal);

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
  std::array<ConstrainedRange, 2> constrainedRanges_;

  // Container holding coordinate arrays (of all supported types) loaded from the input file.
  typedef boost::variant<std::vector<int>,
                         std::vector<float>,
                         std::vector<std::string>> CoordinateValues;
  std::unordered_map<std::string, CoordinateValues> coordsVals_;
  // The array to be interpolated (the payload array).
  boost::multi_array<ExtractedValue, 2> interpolatedArray2D_;
  // Linear interpolation result. Only used when ExtractedValue = float. interpolation.
  float result_;
  // Set to true if result_ is a valid value.
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
  typename std::vector<Coordinate>::const_iterator nextCoordToExtractBy_;
};

}  // namespace ufo

#endif  // UFO_UTILS_DATAEXTRACTOR_DATAEXTRACTOR_H_
