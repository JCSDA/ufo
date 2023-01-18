/*
 * (C) Copyright 2021 Met Office UK
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_FILTERS_OBSFUNCTIONS_DRAWOBSERRORFROMFILE_H_
#define UFO_FILTERS_OBSFUNCTIONS_DRAWOBSERRORFROMFILE_H_

#include <memory>
#include <set>
#include <string>
#include <vector>

#include "oops/util/parameters/OptionalParameter.h"
#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/Parameters.h"

#include "ufo/filters/obsfunctions/DrawValueFromFile.h"
#include "ufo/filters/obsfunctions/ObsFunctionBase.h"
#include "ufo/filters/Variables.h"

namespace ufo {

enum class DispersionMeasure {
  STDEV, VARIANCE, NORMALIZED, FRACTIONAL
};

struct DispersionMeasureParameterTraitsHelper {
  typedef DispersionMeasure EnumType;
  static constexpr char enumTypeName[] = "DispersionMeasure";
  static constexpr util::NamedEnumerator<DispersionMeasure> namedValues[] = {
    { DispersionMeasure::STDEV, "standard deviation" },
    { DispersionMeasure::VARIANCE, "variance" },
    { DispersionMeasure::NORMALIZED, "normalized standard deviation" },
    { DispersionMeasure::FRACTIONAL, "fractional standard deviation" },
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
  /// Which variable to multiply, when Dispersion Measure is normalized standard deviation
  oops::OptionalParameter<Variable> normvariable{"normalization variable", this};
  /// Set a minimum value for the observation uncertainty (default to zero)
  oops::Parameter<float> minValue{"minimum value", 0, this};
  /// List of channel numbers (then deriving an observation error per channel)
  /// If this option is provided, the channel number is implicitly prepended to the list of
  /// interpolation variables and matched exactly.
  oops::OptionalParameter<std::set<int>> chlist{"channels", this};
};


/// \brief Derive the observation uncertainties from a file.
/// \details The file can contain either the standard deviation or variance of each uncertainty,
/// or the full observation-error covariance matrix.
/// Variances are converted to standard deviations in the ObsFunction.
/// See DataExtractor for details on the format of the file.
/// If the optional "normalization variable" exists and "disperson measure" is
/// "normalized standard deviation", then the final output obserr is
/// calculated as a result of normalized standard deviation multiplied by the observation
/// value of the normzliation variable.
///
/// ### example configurations: ###
///
/// \code{.yaml}
///     - Filter: Perform Action
///       filter variables:
///       - name: airTemperature
///         channels: &all_channels 1-3
///       action:
///         name: assign error
///         error function:
///           name: ObsFunction/DrawObsErrorFromFile
///           channels: *all_channels
///           options:
///             file: <filepath>
///             channels: *all_channels
///             interpolation:
///             - name: MetaData/satelliteIdentifier
///               method: exact
///             - name: MetaData/processingCenter
///               method: exact
///             - name: MetaData/pressure
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
  ufo::Variables invars_;
  std::unique_ptr<DrawValueFromFile<float>> drawValueFromFile_;
  DrawObsErrorFromFileParameters options_;
  bool multiplicative_ = false;
  std::vector<int> channels_;
};

}  // namespace ufo

#endif  // UFO_FILTERS_OBSFUNCTIONS_DRAWOBSERRORFROMFILE_H_
