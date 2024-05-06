/*
 * (C) Copyright 2024 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_FILTERS_OBSFUNCTIONS_OBSFUNCTIONELEMENTMULTIPLY_H_
#define UFO_FILTERS_OBSFUNCTIONS_OBSFUNCTIONELEMENTMULTIPLY_H_

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

/// \brief Options controlling ObsFunctionElementMultiply ObsFunction
template <typename FunctionValue>
class ElementMultiplyParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(ElementMultiplyParameters, Parameters)

 public:
  /// Input variables of the linear combination
  oops::RequiredParameter<std::vector<Variable>> variables{"variables", this};
  /// exponent associated with the above variables. Defaults to 1 if no value is entered.
  oops::OptionalParameter<std::vector<FunctionValue>> exponents{"exponents", this};
  /// If \c true (the default), will raise exceptions when trying to perform
  /// invalid or unsupported operations. Otherwise will set output values as
  /// missing.
  oops::OptionalParameter<bool> abortIfInvalid{"abort if invalid operation",
                                               this};
};

// -----------------------------------------------------------------------------

/// \brief Outputs an elementwise multiplication of variables which can be raised
///  to powers allowing implementation of division.
///
/// Example
///
///  obs function:
///    name: ObsFunction/ElementMultiply
///    options:
///      variables: [ObsValue/variable1,
///                  ObsValue/variable2,
///                  ObsValue/variable3]
///      exponents: [1, -1, 3]
///
/// will return   (ObsValue/variable1)^1
///             * (ObsValue/variable2)^-1
///             * (ObsValue/variable3)^3
///
///
///
template <typename FunctionValue>
class ElementMultiply : public ObsFunctionBase<FunctionValue> {
 public:
  explicit ElementMultiply(const eckit::LocalConfiguration &);

  void compute(const ObsFilterData &,
               ioda::ObsDataVector<FunctionValue> &) const;
  FunctionValue power(FunctionValue, FunctionValue) const;
  const ufo::Variables & requiredVariables() const;
 private:
  ElementMultiplyParameters<FunctionValue> options_;
  ufo::Variables invars_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo

#endif  // UFO_FILTERS_OBSFUNCTIONS_OBSFUNCTIONELEMENTMULTIPLY_H_
