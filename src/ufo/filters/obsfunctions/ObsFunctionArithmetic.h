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
///    name: Arithmetic@ObsFunction
///    options:
///      variables: [variable1@ObsValue,
///                  variable2@ObsValue,
///                  variable3@ObsValue]
///      coefficients: [0.1, 0.2, 0.3]
///      exponents: [1, 2, 3]
///      total coefficient: 4
///      total exponent: 5
///      additive constant: 6
///
/// will return 4 * (0.1 * (variable1@ObsValue)^1 +
///                  0.2 * (variable2@ObsValue)^2 +
///                  0.3 * (variable3@ObsValue)^3)^5 + 6
///
/// Can be also be used with the name LinearCombination
/// to output a linear combination of variables
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
