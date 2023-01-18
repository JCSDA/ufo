/*
 * (C) Copyright 2021 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_FILTERS_OBSFUNCTIONS_OBSFUNCTIONARITHMETIC_H_
#define UFO_FILTERS_OBSFUNCTIONS_OBSFUNCTIONARITHMETIC_H_

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

/// \brief Options controlling ObsFunctionArithmetic ObsFunction
template <typename FunctionValue>
class ArithmeticParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(ArithmeticParameters, Parameters)

 public:
  /// Input variables of the linear combination
  oops::RequiredParameter<std::vector<Variable>> variables{"variables", this};
  /// coefficient associated with the above variables
  oops::OptionalParameter<std::vector<FunctionValue>> coefs{"coefs", this};
  /// exponent associated with the above variables
  oops::OptionalParameter<std::vector<FunctionValue>> exponents{"exponents", this};
  /// total exponent
  oops::OptionalParameter<FunctionValue> total_exponent{"total exponent", this};
  /// total multiplicative coefficeint
  oops::OptionalParameter<FunctionValue> total_coeff{"total coefficient", this};
  /// Adds the option to add an intercept or initial value
  oops::OptionalParameter<FunctionValue> intercept{"intercept", this};
  /// Use channel number in the calculation not the value from the variable
  oops::Parameter<bool> useChannelNumber{"use channel numbers", false, this};
};

// -----------------------------------------------------------------------------

/// \brief Outputs an arithmetic combination of variables
///
/// Example
///
///  obs function:
///    name: ObsFunction/Arithmetic
///    options:
///      variables: [ObsValue/variable1,
///                  ObsValue/variable2,
///                  ObsValue/variable3]
///      coefficients: [0.1, 0.2, 0.3]
///      exponents: [1, 2, 3]
///      total coefficient: 4
///      total exponent: 5
///      additive constant: 6
///
/// will return 4 * (0.1 * (ObsValue/variable1)^1 +
///                  0.2 * (ObsValue/variable2)^2 +
///                  0.3 * (ObsValue/variable3)^3)^5 + 6
///
/// Can be also be used with the name LinearCombination
/// to output a linear combination of variables
///
/// Example 1
///
///  obs function:
///    name: ObsFunction/LinearCombination
///    options:
///      variables: [GeoVaLs/representation_error,
///                  ObsError/waterTemperature]
///      coefs: [0.1, 1.0]
///
/// will return 0.1 * GeoVaLs/representation_error +
///             1.0 * ObsError/waterTemperature
///
/// Example 2 - multi-channel
///
///  obs function:
///    name: ObsFunction/LinearCombination
///    channels: &select_chans 6-15, 18-22 # this line may be needed depending on the filter used
///    options:
///      variables:
///      - name: ObsValue/brightnessTemperature
///        channels: *select_chans
///      - name: ObsError/brightnessTemperature
///        channels: *select_chans
///      coefs: [1.0, 0.5]
///
/// will return 1.0 * ObsValue/brightnessTemperature[channel] +
///             0.5 * ObsError/brightnessTemperature[channel]
///
/// Example 3 - multi-channel with intercept and using channel numbers
///
///  obs function:
///    name: ObsFunction/LinearCombination
///    channels: &select_chans 6-15, 18-22 # this line may be needed depending on the filter used
///    options:
///      variables:
///      - name: ObsValue/brightnessTemperature
///        channels: *select_chans
///      coefs: [0.5]
///      intercept: 3.6
///      use channel numbers: true
///
/// will return 3.6 +
///             0.5 * channels
///
///
template <typename FunctionValue>
class Arithmetic : public ObsFunctionBase<FunctionValue> {
 public:
  explicit Arithmetic(const eckit::LocalConfiguration &);

  void compute(const ObsFilterData &,
               ioda::ObsDataVector<FunctionValue> &) const;
  FunctionValue power(FunctionValue, FunctionValue) const;
  const ufo::Variables & requiredVariables() const;
 private:
  ArithmeticParameters<FunctionValue> options_;
  ufo::Variables invars_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo

#endif  // UFO_FILTERS_OBSFUNCTIONS_OBSFUNCTIONARITHMETIC_H_
