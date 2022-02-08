/*
 * (C) Copyright 2021 Met Office UK
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_FILTERS_OBSFUNCTIONS_DRAWOBSERRORFROMFILE_H_
#define UFO_FILTERS_OBSFUNCTIONS_DRAWOBSERRORFROMFILE_H_

#include <memory>
#include <string>

#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/Parameters.h"

#include "ufo/filters/obsfunctions/DrawValueFromFile.h"
#include "ufo/filters/obsfunctions/ObsFunctionBase.h"

namespace ufo {

enum class DispersionMeasure {
  STDEV, VARIANCE
};

struct DispersionMeasureParameterTraitsHelper {
  typedef DispersionMeasure EnumType;
  static constexpr char enumTypeName[] = "DispersionMeasure";
  static constexpr util::NamedEnumerator<DispersionMeasure> namedValues[] = {
    { DispersionMeasure::STDEV, "standard deviation" },
    { DispersionMeasure::VARIANCE, "variance" }
  };
};
}  // namespace ufo


namespace oops {

template <>
struct ParameterTraits<ufo::DispersionMeasure> :
    public EnumParameterTraits<ufo::DispersionMeasureParameterTraitsHelper>
{};

}  // namespace oops


namespace ufo {

/// \brief Options controlling the DrawObsErrorFromFile ObsFunction
class DrawObsErrorFromFileParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(DrawObsErrorFromFileParameters, Parameters)

 public:
  /// Group
  oops::Parameter<std::string> group{"group", "ErrorVariance", this};
  /// Measure of dispersion (standard deviation or variance)
  oops::Parameter<DispersionMeasure> dispersionMeasure{"dispersion measure",
      DispersionMeasure::VARIANCE, this};
};


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
  std::unique_ptr<DrawValueFromFile<float>> drawValueFromFile_;
  std::unique_ptr<DrawObsErrorFromFileParameters> options_;
};

}  // namespace ufo

#endif  // UFO_FILTERS_OBSFUNCTIONS_DRAWOBSERRORFROMFILE_H_
