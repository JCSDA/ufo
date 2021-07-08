/*
 * (C) Crown copyright 2021, Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_UTILS_DATAEXTRACTOR_DATAEXTRACTORBACKEND_H_
#define UFO_UTILS_DATAEXTRACTOR_DATAEXTRACTORBACKEND_H_

#include <string>

namespace ufo
{

template <typename ExtractedValue>
struct DataExtractorInput;

/// \brief Provides data to the DataExtractor.
///
/// \tparam ExtractedValue
///   Type of values extracted by the DataExtractor. Must be `float`, `int` or `std::string`.
///
template <typename ExtractedValue>
class DataExtractorBackend {
 public:
  virtual ~DataExtractorBackend() = default;

  /// \brief Load data for subsequent extraction.
  ///
  /// \param payloadGroup
  ///   Group (e.g. ObsBias or ErrorVariance) containing the payload variable,
  ///   i.e. the variable that will be interpolated. The data source
  ///   must contain exactly one variable from this group.
  ///
  /// \returns An object encapsulating the payload variable, all coordinates indexing it
  /// and the mapping between dimensions of the payload array and coordinates.
  virtual DataExtractorInput<ExtractedValue> loadData(const std::string &payloadGroup) const = 0;
};

}  // namespace ufo

#endif  // UFO_UTILS_DATAEXTRACTOR_DATAEXTRACTORBACKEND_H_
