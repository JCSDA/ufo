/*
 * (C) Crown copyright 2021, Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_FILTERS_OBSFUNCTIONS_OBSFUNCTIONEXPONENTIAL_H_
#define UFO_FILTERS_OBSFUNCTIONS_OBSFUNCTIONEXPONENTIAL_H_

#include <vector>

#include "oops/util/missingValues.h"
#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/RequiredParameter.h"

#include "ufo/filters/obsfunctions/ObsFunctionBase.h"
#include "ufo/filters/Variable.h"
#include "ufo/filters/Variables.h"
#include "ufo/utils/parameters/ParameterTraitsVariable.h"

namespace ufo {

class ObsFilterData;

/// \brief Options controlling ObsFunctionExponential ObsFunction
class ExponentialParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(ExponentialParameters, Parameters)

 public:
  /// Input variable of the exponential function
  oops::RequiredParameter<std::vector<Variable>> variables{"variables", this};
  /// coefficients associated with the above variable (x), such that
  /// y(x) = A*exp(B*x)+C, or y(x) = D if x missing
  oops::Parameter<float> coeffA{"coeffA", 1.0, this};
  oops::Parameter<float> coeffB{"coeffB", 1.0, this};
  oops::Parameter<float> coeffC{"coeffC", 0.0, this};
  oops::Parameter<float> coeffD{"coeffD", util::missingValue<float>(), this};
};

// -----------------------------------------------------------------------------

/// \brief Outputs an exponential function y(x) = A*exp(B*x)+C, or y(x) = D if x missing,
///   where default values are A = 1.0, B = 1.0; C = 0.0; D = missing.
/// \details Example 1
///
///  obs function:
///    name: ObsFunction/Exponential
///    options:
///      variables: [ObsValue/windSpeed]
///      coeffA: 0.1
///      coeffB: -0.3
///      coeffC: 0.2
///      coeffD: 0.12
///
/// will return 0.1 * exp(-0.3 * ObsValue/windSpeed) + 0.2
/// or if ObsValue/windSpeed is missing, 0.12 in that location.
///
/// Example 2 - multi-channel (not sure who might need this, but...)
///
///  obs function:
///    name: ObsFunction/Exponential
///    channels: &select_chans 1-3
///    options:
///      variables:
///      - name: ObsValue/brightnessTemperature
///        channels: *select_chans
///      coeffB: 0.5
///      coeffC: -1.0
///
/// will return 1.0 * exp(0.5 * ObsValue/brightnessTemperature[channel]) -1.0,
/// or missing at locations where ObsValue/brightnessTemperature[channel] is missing,
/// since default values are coeffA: 1.0, coeffB = 1.0; coeffC = 0.0; coeffD = missing.

class Exponential : public ObsFunctionBase<float> {
 public:
  explicit Exponential(const eckit::LocalConfiguration &);
  ~Exponential();

  void compute(const ObsFilterData &,
               ioda::ObsDataVector<float> &) const;
  const ufo::Variables & requiredVariables() const;
 private:
  ExponentialParameters options_;
  ufo::Variables invars_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo

#endif  // UFO_FILTERS_OBSFUNCTIONS_OBSFUNCTIONEXPONENTIAL_H_
