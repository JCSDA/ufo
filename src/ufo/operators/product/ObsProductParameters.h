/*
 * (C) Copyright 2023- UCAR
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

class ObsProductParameters : public ObsOperatorParametersBase {
  OOPS_CONCRETE_PARAMETERS(ObsProductParameters, ObsOperatorParametersBase)

 public:
  /// An optional `variables` parameter, which controls which ObsSpace
  /// variables will be simulated. This option should only be set if this operator is used as a
  /// component of the Composite operator. If `variables` is not set, the operator will simulate
  /// all ObsSpace variables. Please see the documentation of the Composite operator for further
  /// details.
  oops::OptionalParameter<std::vector<ufo::Variable>> variables{
     "variables",
     "List of variables to be simulated",
     this};

  /// h(x) = x(lowest level) * variableToScaleHofxBy
  oops::RequiredParameter<std::string> variableNameToScaleHofxBy{"geovals to scale hofx by", this};

  // Optionally specify the group for the scaling variable (default is GeoVaLs)
  oops::Parameter<std::string> variableGroupToScaleHofxBy{"group of geovals to scale hofx by",
                                                    "GeoVaLs", this};

  /// Optional parameter to raise the scaling variable to a power, h(x) = x(lowest level) *
  /// (variable)^a
  oops::OptionalParameter<float> scalingVariableExponent{"geovals exponent", this};

  /// Optional parameter to specify name of geoval H(x) is to act on, if different from
  /// the simulated variable
  oops::OptionalParameter<std::string> geovalVariable{"geovals to act on", this};
};

}  // namespace ufo

