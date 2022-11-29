/*
 * (C) Crown copyright 2021, Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_UTILS_DATAEXTRACTOR_DATAEXTRACTORNETCDFBACKEND_H_
#define UFO_UTILS_DATAEXTRACTOR_DATAEXTRACTORNETCDFBACKEND_H_

#include <string>

#include "ufo/utils/dataextractor/DataExtractorBackend.h"

namespace ufo
{

/// \brief Produces input for a DataExtractor by loading data from a NetCDF file.
///
/// The input NetCDF file should have the following structure:
///   - The file should contain exactly one variable in a group whose name
///     is the value of the `payloadGroup` parameter passed to the loadData() method.
///     This is the variable containing the values to be extracted (also known as the _payload_).
///     Its type needs to match the template parameter `ExtractedValue`, which must be
///     set to either `float`, `int` or `std::string`.
///   - For most types of data, the above array should be 1D, 2D or 3D. As a special case,
///     if this class is used to extract variances, the file may contain a full 2D covariance
///     matrix or a stack of such matrices stored as a 3D array, with variances located on the
///     diagonal of each matrix. In that case, the array should be equipped with a `full` attribute
///     of type `string` set to `"true"`. Only the diagonal elements (`A[i, i]` for a 2D array and
///     `A[i, i, k]` for a 3D array) will be used during the data extraction process.
///   - 1D coordinates must be defined for each dimension of the payload array.
///     Coordinates can be of type `float`, `int` or `string`. Datetimes should be represented
///     as ISO 8601 strings. Auxiliary coordinates are supported, i.e. there can be more than one
///     coordinate per dimension. Coordinate names should correspond to names of ObsSpace variables.
///     Use the name `MetaData/sensorChannelNumber` for channel numbers (for which there's no
///     dedicated ObsSpace variable).
///   - Coordinate variables should be placed in appropriate groups, e.g. `MetaData`. Because
///     of the limitations of the NetCDF file format, these variables can only be used as auxiliary
///     coordinates of the payload variable (listed in its `coordinates` attribute).
///
/// Example 1: the following NetCDF metadata describe a file encoding the dependence
/// of the diagonal elements of a covariance matrix (with rows and columns corresponding to certain
/// air pressures) on the observation type:
///
/// \code
/// dimensions:
///   rows = 10 ;
///   index = 5 ;
/// variables:
///   int rows(rows) ;
///   int index(index) ;
///
/// group: MetaData {
///   variables:
///     int pressure(rows) ;
///     int observationSubtypeNum(index) ;
/// }
///
/// group: ErrorVariance {
///   variables:
///     float airTemperature(rows, index) ;
///       airTemperature:coordinates = "/MetaData/pressure /MetaData/observationSubTypeNum" ;
/// }
/// \endcode

/// Example 2: the following NetCDF metadata describe a ioda-v1-style file encoding the dependence
/// of a full covariance matrix (with rows and columns corresponding to channel numbers) on
/// multiple variables (`MetaData/latitudeBand`, `MetaData/processingCenter` and
/// `MetaData/satelliteIdentifier`):
///
/// \code
/// dimensions:
///   sensorChannelNumber = 10 ;
///   MetaData/sensorChannelNumber = 10 ;
///   index = 8 ;
/// variables:
///   int index(index) ;
///   int sensorChannelNumber(sensorChannelNumber) ;
///
/// group: MetaData {
///   variables:
///     int sensorChannelNumber(MetaData/sensorChannelNumber) ;
///     int latitudeBand(index) ;
///     int preocessingCenter(index) ;
///     int satelliteIdentifier(index) ;
/// }
/// group: ErrorVariance {
///   variables:
///     float airTemperature(sensorChannelNumber, MetaData/sensorChannelNumber, index) ;
///       airTemperature:coordinates = "MetaData/latitudeBand \
///       MetaData/processingCenter MetaData/satelliteIdentifier" ;
/// }
/// \endcode
/// Notice how a sensorChannelNumber describes the outermost two dimensions.  The very outermost
/// dimension is collapsed and discarded after extracting the diagonal elements of the covariance
/// matrices, so its name can be arbitrary. (This lets us stay conforming to the CF convention,
/// which forbids multiple axes of an array to be indexed by the same coordinate.)
///
template <typename ExtractedValue>
class DataExtractorNetCDFBackend : public DataExtractorBackend<ExtractedValue> {
 public:
  /// \brief Create a new instance.
  ///
  /// \param filepath Path to the NetCDF file that will be read by loadData().
  explicit DataExtractorNetCDFBackend(const std::string &filepath);

  DataExtractorInput<ExtractedValue> loadData(const std::string &payloadGroup) const override;

 private:
  std::string filepath_;
};

}  // namespace ufo

#endif  // UFO_UTILS_DATAEXTRACTOR_DATAEXTRACTORNETCDFBACKEND_H_
