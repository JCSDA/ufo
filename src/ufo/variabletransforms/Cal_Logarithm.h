/*
 * (C) Crown copyright 2024, Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <cmath>
#include <memory>
#include <ostream>
#include <string>
#include <vector>

#include "ufo/variabletransforms/TransformBase.h"

namespace ufo {

/// Configuration parameters for the logarithm conversion.
class Cal_LogarithmParameters : public VariableTransformParametersBase {
  OOPS_CONCRETE_PARAMETERS(Cal_LogarithmParameters,
                           VariableTransformParametersBase);

 public:
  /// Variable to take the logarithm of.
  oops::RequiredParameter<std::string> variable{"variable", this};

  /// Variable group.
  oops::RequiredParameter<std::string> group{"group", this};

  /// Logarithm base (default e, i.e. natural logarithm).
  oops::OptionalParameter<float> base{"base", this};

  /// Optional output variable name.
  oops::OptionalParameter<std::string> outputVariable{"output variable", this};

  /// Optional output group name.
  oops::OptionalParameter<std::string> outputGroup{"output group", this};
};

/*!
* \brief Take the logarithm using the chosen base.
*/
class Cal_Logarithm : public TransformBase {
 public:
  typedef Cal_LogarithmParameters Parameters_;

  Cal_Logarithm(const Parameters_ &options, const ObsFilterData &data,
                const std::shared_ptr<ioda::ObsDataVector<int>> &flags,
                const std::shared_ptr<ioda::ObsDataVector<float>> &obserr);
  // Run check
  void runTransform(const std::vector<bool> &apply) override;

 private:
  /// Variable to take the logarithm of.
  std::string variable_;

  /// Group within which the variable is found.
  std::string group_;

  /// The base of the logarithm.
  float base_ = 0.0;  // 0.0 used to indicate natural logarithm

  /// Optional output variable name.
  std::string outputVariable_;

  /// Optional output group name.
  std::string outputGroup_;
};
}  // namespace ufo
