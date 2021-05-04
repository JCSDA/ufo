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

#include <boost/variant.hpp>
// cpplint misclassifies this file as a c system include
#include <Eigen/Core>  // NOLINT(build/include_order)

namespace ufo
{

/// \brief Input data for the DataExtractor.
///
/// Note: the names of all coordinates are expected to be of the form `Group/var` (ioda-v2 style)
/// rather than `var@Group` (ioda-v1 style).
struct DataExtractorInput {
  /// Array from which values will be extracted
  Eigen::ArrayXXf payloadArray;

  typedef boost::variant<std::vector<int>,
                         std::vector<float>,
                         std::vector<std::string>
                        > Coordinate;
  typedef std::unordered_map<std::string, Coordinate> Coordinates;
  /// Coordinates indexing payloadArray
  Coordinates coordsVals;

  /// Maps coordinate names to dimensions (0 or 1) of the payload array
  std::unordered_map<std::string, int> coord2DimMapping;
  /// Maps dimensions of the payload array (0 or 1) to coordinate names
  std::vector<std::vector<std::string>> dim2CoordMapping;
};

}  // namespace ufo

#endif  // UFO_UTILS_DATAEXTRACTOR_DATAEXTRACTORINPUT_H_
