/*
 * (C) Copyright 2021 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_FILTERS_OBSFUNCTIONS_OBSFUNCTIONLINEARCOMBINATION_H_
#define UFO_FILTERS_OBSFUNCTIONS_OBSFUNCTIONLINEARCOMBINATION_H_

#include <vector>

#include "oops/util/parameters/OptionalParameter.h"
#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/RequiredParameter.h"

#include "ufo/filters/obsfunctions/ObsFunctionBase.h"
#include "ufo/filters/Variable.h"
#include "ufo/filters/Variables.h"
#include "ufo/utils/parameters/ParameterTraitsVariable.h"

namespace ufo {

class ObsFilterData;

/// \brief Options controlling ObsFunctionLinearCombination ObsFunction
template <typename FunctionValue>
class LinearCombinationParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(LinearCombinationParameters, Parameters)

 public:
  /// Input variables of the linear combination
  oops::RequiredParameter<std::vector<Variable>> variables{"variables", this};
  /// coefficient associated with the above variables
  oops::RequiredParameter<std::vector<FunctionValue>> coefs{"coefs", this};
  /// Adds the option to add an intercept or initial value
  oops::OptionalParameter<FunctionValue> intercept{"intercept", this};
  /// Use channel number in the calculation not the value from the variable
  oops::Parameter<bool> useChannelNumber{"use channel numbers", false, this};
};

// -----------------------------------------------------------------------------

/// \brief Outputs a linear combination of variables
///
/// Example 1
///
///  obs function:
///    name: LinearCombination@ObsFunction
///    options:
///      variables: [representation_error@GeoVaLs,
///                  sea_water_temperature@ObsError]
///      coefs: [0.1, 1.0]
///
/// will return 0.1 * representation_error@GeoVaLs +
///             1.0 * sea_water_temperature@ObsError
///
/// Example 2 - multi-channel
///
///  obs function:
///    name: LinearCombination@ObsFunction
///    channels: &select_chans 6-15, 18-22 # this line may be needed depending on the filter used
///    options:
///      variables:
///      - name: brightness_temperature@ObsValue
///        channels: *select_chans
///      - name: brightness_temperature@ObsError
///        channels: *select_chans
///      coefs: [1.0, 0.5]
///
/// will return 1.0 * brightness_temperature_<channel>@ObsValue +
///             0.5 * brightness_temperature_<channel>@ObsError
///
/// Example 3 - multi-channel with intercept and using channel numbers
///
///  obs function:
///    name: LinearCombination@ObsFunction
///    channels: &select_chans 6-15, 18-22 # this line may be needed depending on the filter used
///    options:
///      variables:
///      - name: brightness_temperature@ObsValue
///        channels: *select_chans
///      coefs: [0.5]
///      intercept: 3.6
///      use channel numbers: true
///
/// will return 3.6 +
///             0.5 * channels
///

template <typename FunctionValue>
class LinearCombination : public ObsFunctionBase<FunctionValue> {
 public:
  explicit LinearCombination(const eckit::LocalConfiguration &);

  void compute(const ObsFilterData &,
               ioda::ObsDataVector<FunctionValue> &) const;
  const ufo::Variables & requiredVariables() const;
 private:
  LinearCombinationParameters<FunctionValue> options_;
  ufo::Variables invars_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo

#endif  // UFO_FILTERS_OBSFUNCTIONS_OBSFUNCTIONLINEARCOMBINATION_H_
