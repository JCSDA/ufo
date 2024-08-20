/*
 * (C) Copyright 2023- UCAR
 * (C) Crown Copyright 2024 - UK Met Office.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <string>
#include <vector>

#include "oops/util/parameters/OptionalParameter.h"
#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/RequiredParameter.h"
#include "ufo/filters/Variable.h"
#include "ufo/ObsOperatorParametersBase.h"
#include "ufo/utils/parameters/ParameterTraitsVariable.h"

namespace ufo {

class ObsLogarithmParameters : public ObsOperatorParametersBase {
  OOPS_CONCRETE_PARAMETERS(ObsLogarithmParameters, ObsOperatorParametersBase)

 public:
  /// An optional `variables` parameter, which controls which ObsSpace
  /// variables will be simulated. This option should only be set if this
  /// operator is used as a component of the Composite operator. If `variables`
  /// is not set, the operator will simulate all ObsSpace variables. Please see
  /// the documentation of the Composite operator for further details.
  oops::OptionalParameter<std::vector<ufo::Variable>> variables{
      "variables", "List of variables to be simulated", this};

  /// The logarithm base to use. If not set, the natural logarithm (base e) is
  /// used.
  oops::OptionalParameter<float> logBase{"base", this};
};

}  // namespace ufo
