/*
 * (C) Crown copyright 2021, Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_UTILS_DATAEXTRACTOR_DATAEXTRACTORINPUT_H_
#define UFO_UTILS_DATAEXTRACTOR_DATAEXTRACTORINPUT_H_

#include <string>
#include <unordered_map>
#include <vector>

#include <boost/multi_array.hpp>
#include <boost/variant.hpp>

namespace ufo
{

/// \brief Parts of the input data for the DataExtractor that don't depend on the type of the
/// extracted values.
///
/// Note: the names of all coordinates are expected to be of the form `Group/var`.
struct DataExtractorInputBase {
  /// \brief A coordinate indexing a dimension of the payload array, i.e. the array from which
  /// a DataExtractor will extract data.
  typedef boost::variant<std::vector<int>,
                         std::vector<float>,
                         std::vector<std::string>
                        > Coordinate;
  /// \brief A collection of named coordinate vectors.
  typedef std::unordered_map<std::string, Coordinate> Coordinates;

  /// Coordinates indexing the payload array
  Coordinates coordsVals;

  /// Maps coordinate names to dimensions (0 or 1) of the payload array.
  std::unordered_map<std::string, std::vector<int>> coord2DimMapping;
  /// Maps coordinate names to their dimensionality.
  std::unordered_map<std::string, size_t> coordNDims;
  /// Maps dimensions of the payload array (0 or 1) to coordinate names
  std::vector<std::vector<std::string>> dim2CoordMapping;
};

/// \brief Input data for the DataExtractor.
///
/// \tparam ExtractedValue
///   Type of the values to be extracted. Must be `float`, `int` or `std::string`.
///
/// Note: the names of all coordinates are expected to be of the form `Group/var`.
template <typename ExtractedValue>
struct DataExtractorInput : public DataExtractorInputBase {
  /// Array from which values will be extracted
  boost::multi_array<ExtractedValue, 3> payloadArray;
};

}  // namespace ufo

#endif  // UFO_UTILS_DATAEXTRACTOR_DATAEXTRACTORINPUT_H_
