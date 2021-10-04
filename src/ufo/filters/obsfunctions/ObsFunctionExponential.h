/*
 * (C) Crown copyright 2021, Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_FILTERS_OBSFUNCTIONS_OBSFUNCTIONEXPONENTIAL_H_
#define UFO_FILTERS_OBSFUNCTIONS_OBSFUNCTIONEXPONENTIAL_H_

#include <vector>

#include "oops/util/parameters/RequiredParameter.h"

#include "ufo/filters/ObsFilterData.h"
#include "ufo/filters/obsfunctions/ObsFunctionBase.h"
#include "ufo/filters/Variable.h"
#include "ufo/filters/Variables.h"
#include "ufo/utils/parameters/ParameterTraitsVariable.h"

namespace ufo {

/// \brief Options controlling ObsFunctionExponential ObsFunction
class ExponentialParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(ExponentialParameters, Parameters)

 public:
  /// Input variable of the exponential function
  oops::RequiredParameter<std::vector<Variable>> variables{"variables", this};
  /// coefficients associated with the above variable
  oops::RequiredParameter<std::vector<float>> coefs{"coefs", this};
};

// -----------------------------------------------------------------------------

/// \brief Outputs an exponential function y(x) = A*exp(B*x)+C, or y(x) = D if x missing
/// \details Example 1
///
///  obs function:
///    name: Exponential@ObsFunction
///    options:
///      variables: [wind_speed@ObsValue]
///      coefs: [0.1,-0.3,0.2,0.12]
///
/// will return 0.1 * exp(-0.3 * wind_speed@ObsValue) + 0.2
/// or if wind_speed@ObsValue is missing, 0.12 in that location.
/// Whereas with only 3 coefficients,
///      coefs: [0.1,-0.3,0.2]
/// a missing value indicator is returned at locations where wind_speed@ObsValue is missing.
///
/// Example 2 - multi-channel (not sure who might need this, but...)
///
///  obs function:
///    name: Exponential@ObsFunction
///    channels: &select_chans 1-3
///    options:
///      variables:
///      - name: brightness_temperature@ObsValue
///        channels: *select_chans
///      coefs: [1.0,0.5,-1.0]
///
/// will return 1.0 * exp(0.5 * brightness_temperature_<channel>@ObsValue) -1.0,
/// or missing at locations where brightness_temperature_<channel>@ObsValue is missing.

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
