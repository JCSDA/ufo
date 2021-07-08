/*
 * (C) Copyright 2021 Met Office UK
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_FILTERS_OBSFUNCTIONS_DRAWOBSERRORFROMFILE_H_
#define UFO_FILTERS_OBSFUNCTIONS_DRAWOBSERRORFROMFILE_H_

#include <string>

#include "ufo/filters/obsfunctions/DrawValueFromFile.h"
#include "ufo/filters/obsfunctions/ObsFunctionBase.h"

namespace ufo {

/// \brief Derive observation error from a file representing the variance or the covariance
/// matrix.
/// \details See DataExtractor for details on the format of this file.
///
/// ### example configurations: ###
///
/// \code{.yaml}
///     - Filter: Perform Action
///       filter variables:
///       - name: air_temperature
///         channels: &all_channels 1-3
///       action:
///         name: assign error
///         error function:
///           name: DrawObsErrorFromFile@ObsFunction
///           channels: *all_channels
///           options:
///             file: <filepath>
///             channels: *all_channels
///             interpolation:
///             - name: satellite_id@MetaData
///               method: exact
///             - name: processing_center@MetaData
///               method: exact
///             - name: air_pressure@MetaData
///               method: linear
/// \endcode
///
/// Note that channel number extraction is implicit, using the channels selected and performed as
/// an exact match before any user defined interpolation takes place.
class DrawObsErrorFromFile : public ObsFunctionBase<float> {
 public:
  explicit DrawObsErrorFromFile(const eckit::LocalConfiguration &);

  void compute(const ObsFilterData &, ioda::ObsDataVector<float> &) const;
  const ufo::Variables & requiredVariables() const;

 private:
  DrawValueFromFile<float> drawValueFromFile_;
};

}  // namespace ufo

#endif  // UFO_FILTERS_OBSFUNCTIONS_DRAWOBSERRORFROMFILE_H_
